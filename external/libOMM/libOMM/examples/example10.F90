!==================================================================================================!
! example10 : Gamma-point only, spin unpolarized                                                   !
!                                                                                                  !
! This example is the same as example7, but demonstrates the use of the wrappers and sparse        !
! matrices in a CSC format                                                                         !
!                                                                                                  !
! Sample output can be found in example1.out and example1_libOMM.log                               !
!==================================================================================================!
program example10
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

  character(1) :: order
  character(21) :: file_name

  logical :: new_S, dealloc

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: nprow, npcol, bs_def, icontxt
  integer :: H_dim(2), S_dim(2), D_min_dim(2), ED_min_dim(2), C_min_dim(2), T_dim(2), info
  integer :: desc_H(9), desc_S(9), desc_D_min(9), desc_ED_min(9), desc_C_min(9), desc_T(9)
  integer :: m, n, imd, iscf, i, j, k, l, iostat

  real(dp) :: he, se, e_min, eta
  real(dp), allocatable :: H(:,:), S(:,:), D_min(:,:), ED_min(:,:), C_min(:,:), T(:,:)
  type(matrix) :: Hsp, Ssp

  !**********************************************!

  integer, external :: numroc

  !**********************************************!

#ifdef MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  nprow=4!INT(sqrt(DBLE(mpi_size)))
  npcol=mpi_size/nprow
  order='r' ! important TODO: check how to adapt to different orders
  bs_def=16
  call blacs_get(-1,0,icontxt)
  call blacs_gridinit(icontxt,order,nprow,npcol)
  call ms_scalapack_setup(mpi_rank,mpi_size,nprow,order,bs_def,icontxt=icontxt)
#else
  mpi_rank=0
#endif

  m=832
  n=128

#ifdef MPI
  call blacs_gridinfo(icontxt,i,j,k,l)
  H_dim(1)=numroc(m,bs_def,k,0,nprow)
  H_dim(2)=numroc(m,bs_def,l,0,npcol)
  call descinit(desc_H,m,m,bs_def,bs_def,0,0,icontxt,H_dim(1),info)
  allocate(H(H_dim(1),H_dim(2)))
  S_dim=H_dim
  desc_S=desc_H
  allocate(S(S_dim(1),S_dim(2)))
  D_min_dim=H_dim
  desc_D_min=desc_H
  allocate(D_min(D_min_dim(1),D_min_dim(2)))
  ED_min_dim=H_dim
  desc_ED_min=desc_H
  allocate(ED_min(ED_min_dim(1),ED_min_dim(2)))
  C_min_dim(1)=numroc(n,bs_def,k,0,nprow)
  C_min_dim(2)=numroc(m,bs_def,l,0,npcol)
  call descinit(desc_C_min,n,m,bs_def,bs_def,0,0,icontxt,C_min_dim(1),info)
  allocate(C_min(C_min_dim(1),C_min_dim(2)))
  T_dim=1
  allocate(T(T_dim(1),T_dim(2)))
#else
  allocate(H(m,m))
  allocate(S(m,m))
  allocate(D_min(m,m))
  allocate(ED_min(m,m))
  allocate(C_min(n,m))
  allocate(T(m,m))
#endif

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
        H=0.0_dp
        if (iscf==1) S=0.0_dp
        do i=1,m*m
           read(10,'(2(1x,i5),2(1x,f21.15))',iostat=iostat) k, l, he, se
           if (iostat/=0) exit
#ifdef MPI
           call pdelset(H,k,l,desc_H,he-eta*se)
           if (iscf==1) call pdelset(S,k,l,desc_S,se)
#else
           H(k,l)=he-eta*se
           S(k,l)=se
#endif
        end do
        close(10)

        if (iscf==1) then
           new_S=.true.
        else
           new_S=.false.
        end if

        ! If the sparse matrix is not stored locally in a CSC format, then convert it into a CSC format and input it in omm_pddbc_spm
        call m_register_psp_thre(Hsp,H,desc_H,'csc',0.0_dp)
        call m_register_psp_thre(Ssp,S,desc_S,'csc',0.0_dp)


#ifdef MPI
        call omm_pddbc_spm(m,n,H_dim,Hsp,desc_H,.true.,S_dim,Ssp,desc_S,new_S,e_min,D_min_dim,D_min,desc_D_min,.false.,&
             eta,C_min_dim,C_min,desc_C_min,.false.,.false.,T_dim,T,desc_T,0.0_dp,0,1,1,-1.0_dp,.true.,.false.,mpi_rank,mpi_size,&
             nprow,order,bs_def,icontxt)
#else
        call omm_sdden_ref(m,n,H,.true.,S,new_S,e_min,D_min,.false.,eta,C_min,.false.,.false.,T,0.0_dp,0,1,1,-1.0_dp,.true.,.false.,&
             mpi_rank)
#endif

        if (mpi_rank==0) then
           print('(a,f21.15)'), 'e_min : ', 2.0_dp*e_min
        end if

        if (mpi_rank==0) then
           print('(a,f21.15)'), 'D_11  : ', D_min(1,1)
           print('()')
        end if

     end do

     if (imd==2) then
        dealloc=.true.
     else
        dealloc=.false.
     end if
#ifdef MPI
     call omm_pddbc_spm(m,n,H_dim,Hsp,desc_H,.true.,S_dim,Ssp,desc_S,.false.,e_min,ED_min_dim,ED_min,desc_ED_min,.true.,eta,&
          C_min_dim,C_min,desc_C_min,.false.,.false.,T_dim,T,desc_T,0.0_dp,0,1,1,-1.0_dp,.true.,dealloc,mpi_rank,mpi_size,&
          nprow,order,bs_def,icontxt)
#else
     call omm_sdden_ref(m,n,H,.true.,S,.false.,e_min,ED_min,.true.,eta,C_min,.false.,.false.,T,0.0_dp,0,1,1,-1.0_dp,.true.,dealloc,&
          mpi_rank)
#endif

     if (mpi_rank==0) then
        print('(a,f21.15)'), 'ED_11 : ', ED_min(1,1)
        print('()')
     end if

  end do

  deallocate(T)
  deallocate(C_min)
  deallocate(ED_min)
  deallocate(D_min)
  deallocate(S)
  deallocate(H)
  call m_deallocate(Hsp)
  call m_deallocate(Ssp)

#ifdef MPI
  call mpi_finalize(mpi_err)
#endif

end program example10
