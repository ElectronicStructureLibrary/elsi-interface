!==================================================================================================!
! example3 : k points, spin unpolarized                                                            !
!                                                                                                  !
! This example demonstrates a typical calculation with an outer loop of two MD iterations, and an  !
! inner loop of five SCF iterations per MD step. For convenience, the Hamiltonian and overlap      !
! matrices at each point have been pre-generated and are read in from file; in a real code,        !
! however, they should be calculated at each step based on the output density matrix given by      !
! libOMM.                                                                                          !
!                                                                                                  !
! This example is for a system with 988 basis orbitals and 152 occupied states at each k point     !
! (generated from a 76-atom C nanotube). There are four k points, weighted equally. At the end of  !
! each SCF iteration, the Kohn-Sham energy 0.25*(e_min(1)+e_min(2)+e_min(3)+e_min(4)) is printed   !
! out, together with the first element of the density matrix at each k point as an extra check. At !
! the end of each MD iteration, the energy-weighted density matrix at each k point is also         !
! calculated, and its first element is printed out.                                                !
!                                                                                                  !
! Things to note:                                                                                  !
!   1. The eigenspectrum shift parameter eta in this case can be set to 0 (see example1/2).        !
!   2. There are no optional arguments in the call to libOMM. Therefore, matrices which are not    !
!      needed should simpy be passed without having been allocated by MatrixSwitch (m_allocate     !
!      routine). Other variables which are not needed will be ignored by libOMM.                   !
!   3. Try enabling Cholesky factorization (precon=1) or preconditioning (precon=3) in the call to !
!      libOMM to see how the convergence speed is affected. Preconditioning is even more effective !
!      if a T matrix is provided (scale_T should be set around 10 Ry). Take care that, for         !
!      Cholesky factorization, S and H will be overwritten by U and U^(-T)*H*U^(-1).               !
!   4. The dealloc variable should only be .true. for the very last call. This is because libOMM   !
!      stores and reuses internal information from one call to the next.                           !
!                                                                                                  !
! Sample output can be found in example3.out and example3_libOMM.log                               !
!==================================================================================================!
program example3
  use MatrixSwitch

  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!

  character(5) :: m_storage
  character(3) :: m_operation
  character(21) :: file_name

  logical :: new_S, dealloc

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: m, n, nk, imd, iscf, i, j, k, l, iostat

  real(dp) :: eta
  real(dp), allocatable :: he(:), se(:), e_min(:)

  complex(dp) :: cmplx_he, cmplx_se, el

  type(matrix) :: T
  type(matrix), allocatable :: H(:), S(:), D_min(:), ED_min(:), C_min(:)

  !**********************************************!

#ifdef MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  call ms_scalapack_setup(mpi_size,1,'c',32)

  m_storage='pzdbc'
  m_operation='lap'
#else
  mpi_rank=0

  m_storage='szden'
  m_operation='ref'
#endif

  m=988
  n=152
  nk=4

  allocate(e_min(nk))
  allocate(he(2*nk))
  allocate(se(2*nk))
  allocate(H(nk))
  allocate(S(nk))
  allocate(D_min(nk))
  allocate(ED_min(nk))
  allocate(C_min(nk))

  do i=1,nk
    call m_allocate(H(i),m,m,m_storage)
    call m_allocate(S(i),m,m,m_storage)
    call m_allocate(D_min(i),m,m,m_storage)
    call m_allocate(ED_min(i),m,m,m_storage)
    call m_allocate(C_min(i),n,m,m_storage)
  end do

  eta=0.0_dp

  do imd=1,2

    do iscf=1,5

      if (mpi_rank==0) then
        print('(a)'), '//////////////////////////'
        print('(a,1x,i1,a,1x,i1,a)'), '/ MD STEP', imd, ' - SCF STEP', iscf, ' /'
        print('(a)'), '//////////////////////////'
        print('()')
      end if

      write(file_name,'(a,i1,a,i1,a)') 'data/example3_', imd, '_', iscf, '.dat'
      open(10,file=file_name)
      do i=1,nk
        call m_set(H(i),'a',0.0_dp,0.0_dp)
        if (iscf==1) call m_set(S(i),'a',0.0_dp,0.0_dp)
      end do
      do i=1,m*m
        read(10,'(2(1x,i5),16(1x,f21.15))',iostat=iostat) k, l, he(1:2), se(1:2), &
                                                                he(3:4), se(3:4), &
                                                                he(5:6), se(5:6), &
                                                                he(7:8), se(7:8)
        if (iostat/=0) exit
        do j=1,nk
          cmplx_se=cmplx(se(2*j-1),se(2*j),dp)
          cmplx_he=cmplx(he(2*j-1),he(2*j),dp)-eta*cmplx_se
          call m_set_element(H(j),k,l,cmplx_he)
          if (iscf==1) call m_set_element(S(j),k,l,cmplx_se)
        end do
      end do
      close(10)

      if (iscf==1) then
        new_S=.true.
      else
        new_S=.false.
      end if
      do i=1,nk
        if (mpi_rank==0) print('(a,1x,i1,a,1x,i1,a)'), 'k point', i, ' of', nk, '...'
        if ((imd>0) .and. (iscf==1) .and. (i>1)) then
          call m_add(C_min(i-1),'n',C_min(i),1.0_dp,0.0_dp,m_operation)
          call omm(m,n,H(i),S(i),new_S,e_min(i),D_min(i),.false.,eta,C_min(i),.true.,T,0.0_dp,1,nk,i,-1.0_dp,.true.,.false.,&
                   m_storage,m_operation,mpi_rank)
        else
          call omm(m,n,H(i),S(i),new_S,e_min(i),D_min(i),.false.,eta,C_min(i),.false.,T,0.0_dp,1,nk,i,-1.0_dp,.true.,.false.,&
                   m_storage,m_operation,mpi_rank)
        end if
      end do

      if (mpi_rank==0) then
        print('()')
        print('(a,f21.15)'), 'e_min : ', 0.25_dp*sum(e_min(:))
      end if

      do i=1,nk
        call m_get_element(D_min(i),1,1,el)

        if (mpi_rank==0) print('(2(a,f21.15),a,1x,i1,a)'), 'D_11  : ', real(el,dp), ' , ', aimag(el), ' (k point', i, ')'
      end do
      if (mpi_rank==0) print('()')

    end do

    do i=1,nk
      if ((imd==2) .and. (i==nk)) then
        dealloc=.true.
      else
        dealloc=.false.
      end if
      call omm(m,n,H(i),S(i),.false.,e_min(i),ED_min(i),.true.,eta,C_min(i),.false.,T,0.0_dp,1,nk,i,-1.0_dp,.true.,dealloc,&
               m_storage,m_operation,mpi_rank)

      call m_get_element(ED_min(i),1,1,el)
    
      if (mpi_rank==0) print('(2(a,f21.15),a,1x,i1,a)'), 'ED_11  : ', real(el,dp), ' , ', aimag(el), ' (k point', i, ')'
    end do
    if (mpi_rank==0) print('()')

  end do

  do i=1,nk
    call m_deallocate(C_min(i))
    call m_deallocate(ED_min(i))
    call m_deallocate(D_min(i))
    call m_deallocate(S(i))
    call m_deallocate(H(i))
  end do

  deallocate(C_min)
  deallocate(ED_min)
  deallocate(D_min)
  deallocate(S)
  deallocate(H)
  deallocate(se)
  deallocate(he)
  deallocate(e_min)

#ifdef MPI
  call mpi_finalize(mpi_err)
#endif

end program example3
