!==================================================================================================!
! example9 : k points, spin unpolarized                                                            !
!                                                                                                  !
! This example is the same as example3, but demonstrates the use of the wrappers and sparse        !
! matrices.                                                                                        !
!                                                                                                  !
! Sample output can be found in example3.out and example3_libOMM.log                               !
!                                                                                                  !
! flavour of the OMM functional:                                                                   !
! 0 for basic                                                                                      !
! 1 for Cholesky factorization, S provided,                                                        !
!   not available in the current version of sparse libomm,                                         !
!   TODO: add SuperLU and keep sparse factors instead of computing their products                  !
! 2 for Cholesky factorization, U provided,                                                        !
!   not available in the current version of sparse libomm,                                         !
!   TODO: keep sparse factors instead of computing their products                                  !
! 3 for preconditioning, S provided (T optional),                                                  !
!   not available in the current version of sparse libomm,                                         !
!   TODO: level1 pspBLAS                                                                           !
!==================================================================================================!
program example9

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

  character(1) :: order
  character(21) :: file_name

  logical :: new_S, dealloc

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: nprow, npcol, bs_def, icontxt
  integer :: H_dim(2), S_dim(2), D_min_dim(2), ED_min_dim(2), C_min_dim(2), T_dim(2), info
  integer :: desc_H(9), desc_S(9), desc_D_min(9), desc_ED_min(9), desc_C_min(9), desc_T(9)
  integer :: m, n, nk, imd, iscf, i, j, k, l, iostat, flavour

  real(dp) :: eta
  real(dp), allocatable :: he(:), se(:), e_min(:)

  complex(dp) :: cmplx_he, cmplx_se

  complex(dp), allocatable :: H(:,:,:), S(:,:,:), D_min(:,:,:), ED_min(:,:,:), C_min(:,:,:), T(:,:)

  !**********************************************!

  integer, external :: numroc

  !**********************************************!

#ifdef MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  nprow=INT(sqrt(DBLE(mpi_size)))
  npcol=nprow
  order='r' ! important TODO: check how to adapt to different orders
  bs_def=16
  call blacs_get(-1,0,icontxt)
  call blacs_gridinit(icontxt,order,nprow,npcol)
#else
  mpi_rank=0
#endif

  m=988
  n=152
  nk=4
  flavour=0

  allocate(e_min(nk))
  allocate(he(2*nk))
  allocate(se(2*nk))
#ifdef MPI
  call blacs_gridinfo(icontxt,i,j,k,l)
  H_dim(1)=numroc(m,bs_def,k,0,nprow)
  H_dim(2)=numroc(m,bs_def,l,0,npcol)
  call descinit(desc_H,m,m,bs_def,bs_def,0,0,icontxt,H_dim(1),info)
  allocate(H(H_dim(1),H_dim(2),nk))
  S_dim=H_dim
  desc_S=desc_H
  allocate(S(S_dim(1),S_dim(2),nk))
  D_min_dim=H_dim
  desc_D_min=desc_H
  allocate(D_min(D_min_dim(1),D_min_dim(2),nk))
  ED_min_dim=H_dim
  desc_ED_min=desc_H
  allocate(ED_min(ED_min_dim(1),ED_min_dim(2),nk))
  C_min_dim(1)=numroc(n,bs_def,k,0,nprow)
  C_min_dim(2)=numroc(m,bs_def,l,0,npcol)
  call descinit(desc_C_min,n,m,bs_def,bs_def,0,0,icontxt,C_min_dim(1),info)
  allocate(C_min(C_min_dim(1),C_min_dim(2),nk))
  T_dim=1
  allocate(T(T_dim(1),T_dim(2)))
#else
  allocate(H(m,m,nk))
  allocate(S(m,m,nk))
  allocate(D_min(m,m,nk))
  allocate(ED_min(m,m,nk))
  allocate(C_min(n,m,nk))
  allocate(T(m,m))
#endif

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
      H=cmplx_0
      if (iscf==1) S=cmplx_0
      do i=1,m*m
        read(10,'(2(1x,i5),16(1x,f21.15))',iostat=iostat) k, l, he(1:2), se(1:2), &
                                                                he(3:4), se(3:4), &
                                                                he(5:6), se(5:6), &
                                                                he(7:8), se(7:8)
        if (iostat/=0) exit
        do j=1,nk
          cmplx_se=cmplx(se(2*j-1),se(2*j),dp)
          cmplx_he=cmplx(he(2*j-1),he(2*j),dp)-eta*cmplx_se
#ifdef MPI
        call pzelset(H(1,1,j),k,l,desc_H,cmplx_he)
        if (iscf==1) call pzelset(S(1,1,j),k,l,desc_S,cmplx_se)
#else
          H(k,l,j)=cmplx_he
          if (iscf==1) S(k,l,j)=cmplx_se
#endif
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
#ifdef MPI
          C_min(1:C_min_dim(1),1:C_min_dim(2),i)=C_min(1:C_min_dim(1),1:C_min_dim(2),i-1)
          call omm_pzdbc2pzcsc_psp(m,n,H_dim,H(1,1,i),desc_H,.true.,S_dim,S(1,1,i),desc_S,new_S,e_min(i),D_min_dim,D_min(1,1,i),&
                   desc_D_min,.false.,eta,C_min_dim,C_min(1,1,i),desc_C_min,.true.,.false.,T_dim,T,desc_T,0.0_dp,flavour,nk,&
                   i,-1.0_dp,.true.,.false.,mpi_rank,mpi_size,nprow,order,bs_def,icontxt)
#else
          C_min(1:n,1:m,i)=C_min(1:n,1:m,i-1)
          call omm_szden_ref(m,n,H(1,1,i),.true.,S(1,1,i),new_S,e_min(i),D_min(1,1,i),.false.,eta,C_min(1,1,i),.true.,.false.,T,&
                   0.0_dp,flavour,nk,i,-1.0_dp,.true.,.false.,mpi_rank)
#endif
        else
#ifdef MPI
          call omm_pzdbc2pzcsc_psp(m,n,H_dim,H(1,1,i),desc_H,.true.,S_dim,S(1,1,i),desc_S,new_S,e_min(i),D_min_dim,D_min(1,1,i),&
                   desc_D_min,.false.,eta,C_min_dim,C_min(1,1,i),desc_C_min,.false.,.false.,T_dim,T,desc_T,0.0_dp,flavour,nk,&
                   i,-1.0_dp,.true.,.false.,mpi_rank,mpi_size,nprow,order,bs_def,icontxt)
#else
          call omm_szden_ref(m,n,H(1,1,i),.true.,S(1,1,i),new_S,e_min(i),D_min(1,1,i),.false.,eta,C_min(1,1,i),.false.,.false.,T,&
                   0.0_dp,flavour,nk,i,-1.0_dp,.true.,.false.,mpi_rank)
#endif
        end if
      end do

      if (mpi_rank==0) then
        print('()')
        print('(a,f21.15)'), 'e_min : ', 0.25_dp*sum(e_min(:))
      end if

      do i=1,nk
        if (mpi_rank==0) print('(2(a,f21.15),a,1x,i1,a)'), 'D_11  : ', real(D_min(1,1,i),dp), ' , ', aimag(D_min(1,1,i)), &
                             ' (k point', i, ')'
      end do
      if (mpi_rank==0) print('()')

    end do

    do i=1,nk
      if ((imd==2) .and. (i==nk)) then
        dealloc=.true.
      else
        dealloc=.false.
      end if
#ifdef MPI
      call omm_pzdbc2pzcsc_psp(m,n,H_dim,H(1,1,i),desc_H,.true.,S_dim,S(1,1,i),desc_S,.false.,e_min(i),ED_min_dim,ED_min(1,1,i),&
               desc_ED_min,.true.,eta,C_min_dim,C_min(1,1,i),desc_C_min,.false.,.false.,T_dim,T,desc_T,0.0_dp,flavour,nk,i,&
               -1.0_dp,.true.,dealloc,mpi_rank,mpi_size,nprow,order,bs_def,icontxt)
#else
      call omm_szden_ref(m,n,H(1,1,i),.true.,S(1,1,i),.false.,e_min(i),ED_min(1,1,i),.true.,eta,C_min(1,1,i),.false.,.false.,T,&
               0.0_dp,flavour,nk,i,-1.0_dp,.true.,dealloc,mpi_rank)
#endif
    
      if (mpi_rank==0) print('(2(a,f21.15),a,1x,i1,a)'), 'ED_11  : ', real(ED_min(1,1,i),dp), ' , ', aimag(ED_min(1,1,i)), &
                           ' (k point', i, ')'
    end do
    if (mpi_rank==0) print('()')

  end do

  deallocate(T)
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

end program example9
