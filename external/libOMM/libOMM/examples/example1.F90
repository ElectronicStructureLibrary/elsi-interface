!==================================================================================================!
! example1 : Gamma-point only, spin unpolarized                                                    !
!                                                                                                  !
! This example demonstrates a typical calculation with an outer loop of two MD iterations, and an  !
! inner loop of five SCF iterations per MD step. For convenience, the Hamiltonian and overlap      !
! matrices at each point have been pre-generated and are read in from file; in a real code,        !
! however, they should be calculated at each step based on the output density matrix given by      !
! libOMM.                                                                                          !
!                                                                                                  !
! This example is for a system with 832 basis orbitals and 128 occupied states (generated from a   !
! 64-atom supercell of bulk Si). At the end of each SCF iteration, the Kohn-Sham energy 2*e_min is !
! printed out, together with the first element of the density matrix as an extra check. At the end !
! of each MD iteration, the energy-weighted density matrix is also calculated, and its first       !
! element is printed out.                                                                          !
!                                                                                                  !
! Things to note:                                                                                  !
!   1. The eigenspectrum shift parameter eta has to be set larger than 0 for convergence, since    !
!      the occupied spectrum extends beyond 0 (this is typically a sign of a poor basis). Try      !
!      changing eta to see how the convergence speed is affected. However, also take into account: !
!   2. The Hamiltonian has to be provided to libOMM *already shifted by eta* (H -> H-eta*S)        !
!   3. There are no optional arguments in the call to libOMM. Therefore, matrices which are not    !
!      needed should simpy be passed without having been allocated by MatrixSwitch (m_allocate     !
!      routine). Other variables which are not needed will be ignored by libOMM.                   !
!   4. Try enabling Cholesky factorization (precon=1) or preconditioning (precon=3) in the call to !
!      libOMM to see how the convergence speed is affected. Preconditioning is even more effective !
!      if a T matrix is provided (scale_T should be set around 10 Ry). Take care that, for         !
!      Cholesky factorization, S and H will be overwritten by U and U^(-T)*H*U^(-1).               !
!   5. The dealloc variable should only be .true. for the very last call. This is because libOMM   !
!      stores and reuses internal information from one call to the next.                           !
!                                                                                                  !
! Sample output can be found in example1.out and example1_libOMM.log                               !
!==================================================================================================!
program example1
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
  integer :: m, n, imd, iscf, i, j, k, l, iostat

  real(dp) :: he, se, el, e_min, eta

  type(matrix) :: H, S, D_min, ED_min, C_min, T

  !**********************************************!

#ifdef MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  call ms_scalapack_setup(mpi_rank,mpi_size,1,'c',32)

  m_storage='pddbc'
  m_operation='lap'
#else
  mpi_rank=0

  m_storage='sdden'
  m_operation='ref'
#endif

  m=832
  n=128

  call m_allocate(H,m,m,m_storage)
  call m_allocate(S,m,m,m_storage)
  call m_allocate(D_min,m,m,m_storage)
  call m_allocate(ED_min,m,m,m_storage)
  call m_allocate(C_min,n,m,m_storage)

  eta=0.36748994224532205_dp

  do imd=1,2

    do iscf=1,5

      if (mpi_rank==0) then
        print('(a)'), '//////////////////////////'
        print('(a,1x,i1,a,1x,i1,a)'), '/ MD STEP', imd, ' - SCF STEP', iscf, ' /'
        print('(a)'), '//////////////////////////'
        print('()')
      end if

      write(file_name,'(a,i1,a,i1,a)') 'data/example1_', imd, '_', iscf, '.dat'
      open(10,file=file_name)
      call m_set(H,'a',0.0_dp,0.0_dp)
      if (iscf==1) call m_set(S,'a',0.0_dp,0.0_dp)
      do i=1,m*m
        read(10,'(2(1x,i5),2(1x,f21.15))',iostat=iostat) k, l, he, se
        if (iostat/=0) exit
        call m_set_element(H,k,l,he-eta*se)
        if (iscf==1) call m_set_element(S,k,l,se)
      end do
      close(10)

      if (iscf==1) then
        new_S=.true.
      else
        new_S=.false.
      end if
      call omm(m,n,H,S,new_S,e_min,D_min,.false.,eta,C_min,.false.,T,0.0_dp,0,1,1,-1.0_dp,.true.,.false.,m_storage,m_operation,&
               mpi_rank)

      if (mpi_rank==0) then
        print('(a,f21.15)'), 'e_min : ', 2.0_dp*e_min
      end if

      call m_get_element(D_min,1,1,el)

      if (mpi_rank==0) then
        print('(a,f21.15)'), 'D_11  : ', el
        print('()')
      end if

    end do

    if (imd==2) then
      dealloc=.true.
    else
      dealloc=.false.
    end if
    call omm(m,n,H,S,.false.,e_min,ED_min,.true.,eta,C_min,.false.,T,0.0_dp,0,1,1,-1.0_dp,.true.,dealloc,m_storage,m_operation,&
             mpi_rank)

    call m_get_element(ED_min,1,1,el)

    if (mpi_rank==0) then
      print('(a,f21.15)'), 'ED_11 : ', el
      print('()')
    end if

  end do

  call m_deallocate(C_min)
  call m_deallocate(ED_min)
  call m_deallocate(D_min)
  call m_deallocate(S)
  call m_deallocate(H)

#ifdef MPI
  call mpi_finalize(mpi_err)
#endif

end program example1
