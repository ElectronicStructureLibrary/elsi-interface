!==================================================================================================!
! example5 : Gamma-point only, spin polarized                                                      !
!                                                                                                  !
! This example is the same as example2, but demonstrates the use of the wrappers.                  !
!                                                                                                  !
! Sample output can be found in example2.out and example2_libOMM.log                               !
!==================================================================================================!
program example5

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
  integer :: H_dim(2), S_dim(2), D_min_dim(2), ED_min_dim(2), C_min_up_dim(2), C_min_down_dim(2), T_dim(2), info
  integer :: desc_H(9), desc_S(9), desc_D_min(9), desc_ED_min(9), desc_C_min_up(9), desc_C_min_down(9), desc_T(9)
  integer :: m, n_up, n_down, imd, iscf, i, j, k, l, iostat

  real(dp) :: hue, hde, se, e_min_up, e_min_down, eta

  real(dp), allocatable :: H_up(:,:), H_down(:,:), S(:,:), D_min_up(:,:), D_min_down(:,:), ED_min_up(:,:), ED_min_down(:,:)
  real(dp), allocatable :: C_min_up(:,:), C_min_down(:,:), T(:,:)

  !**********************************************!

  integer, external :: numroc

  !**********************************************!

#ifdef MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  nprow=1
  npcol=mpi_size/nprow
  order='c'
  bs_def=32
  call blacs_get(-1,0,icontxt)
  call blacs_gridinit(icontxt,order,nprow,npcol)
#else
  mpi_rank=0
#endif

  m=832
  n_up=131
  n_down=125

#ifdef MPI
  call blacs_gridinfo(icontxt,i,j,k,l)
  H_dim(1)=numroc(m,bs_def,k,0,nprow)
  H_dim(2)=numroc(m,bs_def,l,0,npcol)
  call descinit(desc_H,m,m,bs_def,bs_def,0,0,icontxt,H_dim(1),info)
  allocate(H_up(H_dim(1),H_dim(2)))
  allocate(H_down(H_dim(1),H_dim(2)))
  S_dim=H_dim
  desc_S=desc_H
  allocate(S(S_dim(1),S_dim(2)))
  D_min_dim=H_dim
  desc_D_min=desc_H
  allocate(D_min_up(D_min_dim(1),D_min_dim(2)))
  allocate(D_min_down(D_min_dim(1),D_min_dim(2)))
  ED_min_dim=H_dim
  desc_ED_min=desc_H
  allocate(ED_min_up(ED_min_dim(1),ED_min_dim(2)))
  allocate(ED_min_down(ED_min_dim(1),ED_min_dim(2)))
  C_min_up_dim(1)=numroc(n_up,bs_def,k,0,nprow)
  C_min_up_dim(2)=numroc(m,bs_def,l,0,npcol)
  call descinit(desc_C_min_up,n_up,m,bs_def,bs_def,0,0,icontxt,C_min_up_dim(1),info)
  allocate(C_min_up(C_min_up_dim(1),C_min_up_dim(2)))
  C_min_down_dim(1)=numroc(n_down,bs_def,k,0,nprow)
  C_min_down_dim(2)=numroc(m,bs_def,l,0,npcol)
  call descinit(desc_C_min_down,n_down,m,bs_def,bs_def,0,0,icontxt,C_min_down_dim(1),info)
  allocate(C_min_down(C_min_down_dim(1),C_min_down_dim(2)))
  T_dim=1
  allocate(T(T_dim(1),T_dim(2)))
#else
  allocate(H_up(m,m))
  allocate(H_down(m,m))
  allocate(S(m,m))
  allocate(D_min_up(m,m))
  allocate(D_min_down(m,m))
  allocate(ED_min_up(m,m))
  allocate(ED_min_down(m,m))
  allocate(C_min_up(n_up,m))
  allocate(C_min_down(n_down,m))
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

      write(file_name,'(a,i1,a,i1,a)') 'data/example2_', imd, '_', iscf, '.dat'
      open(10,file=file_name)
      H_up=0.0_dp
      H_down=0.0_dp
      if (iscf==1) S=0.0_dp
      do i=1,m*m
        read(10,'(2(1x,i5),3(1x,f21.15))',iostat=iostat) k, l, hue, hde, se
        if (iostat/=0) exit
#ifdef MPI
        call pdelset(H_up,k,l,desc_H,hue-eta*se)
        call pdelset(H_down,k,l,desc_H,hde-eta*se)
        if (iscf==1) call pdelset(S,k,l,desc_S,se)
#else
        H_up(k,l)=hue-eta*se
        H_down(k,l)=hde-eta*se
        if (iscf==1) S(k,l)=se
#endif
      end do
      close(10)

      if (iscf==1) then
        new_S=.true.
      else
        new_S=.false.
      end if
      if (mpi_rank==0) print('(a)'), 'up spin...'
#ifdef MPI
      call omm_pddbc_lap(m,n_up,H_dim,H_up,desc_H,.true.,S_dim,S,desc_S,new_S,e_min_up,D_min_dim,D_min_up,desc_D_min,.false.,eta,&
               C_min_up_dim,C_min_up,desc_C_min_up,.false.,.false.,T_dim,T,desc_T,0.0_dp,0,2,1,-1.0_dp,.true.,.false.,mpi_rank,&
               mpi_size,nprow,order,bs_def,icontxt)
#else
      call omm_sdden_ref(m,n_up,H_up,.true.,S,new_S,e_min_up,D_min_up,.false.,eta,C_min_up,.false.,.false.,T,0.0_dp,0,2,1,-1.0_dp,&
               .true.,.false.,mpi_rank)
#endif
      if (mpi_rank==0) print('(a)'), 'down spin...'
#ifdef MPI
      call omm_pddbc_lap(m,n_down,H_dim,H_down,desc_H,.true.,S_dim,S,desc_S,new_S,e_min_down,D_min_dim,D_min_down,desc_D_min,&
               .false.,eta,C_min_down_dim,C_min_down,desc_C_min_down,.false.,.false.,T_dim,T,desc_T,1.0_dp,0,2,2,-1.0_dp,.true.,&
               .false.,mpi_rank,mpi_size,nprow,order,bs_def,icontxt)
#else
      call omm_sdden_ref(m,n_down,H_down,.true.,S,new_S,e_min_down,D_min_down,.false.,eta,C_min_down,.false.,.false.,T,1.0_dp,0,2,&
               2,-1.0_dp,.true.,.false.,mpi_rank)
#endif

      if (mpi_rank==0) then
        print('()')
        print('(a,f21.15)'), 'e_min : ', e_min_up+e_min_down
      end if

      if (mpi_rank==0) then
        print('(a,f21.15,a,f21.15,a)'), 'D_11  : ', D_min_up(1,1), ' (up), ', D_min_down(1,1), ' (down)'
        print('()')
      end if

    end do

#ifdef MPI
    call omm_pddbc_lap(m,n_up,H_dim,H_up,desc_H,.true.,S_dim,S,desc_S,.false.,e_min_up,ED_min_dim,ED_min_up,desc_ED_min,.true.,eta,&
             C_min_up_dim,C_min_up,desc_C_min_up,.false.,.false.,T_dim,T,desc_T,0.0_dp,0,2,1,-1.0_dp,.true.,.false.,mpi_rank,&
             mpi_size,nprow,order,bs_def,icontxt)
#else
    call omm_sdden_ref(m,n_up,H_up,.true.,S,.false.,e_min_up,ED_min_up,.true.,eta,C_min_up,.false.,.false.,T,0.0_dp,0,2,1,-1.0_dp,&
             .true.,.false.,mpi_rank)
#endif
    if (imd==2) then
      dealloc=.true.
    else
      dealloc=.false.
    end if
#ifdef MPI
    call omm_pddbc_lap(m,n_down,H_dim,H_down,desc_H,.true.,S_dim,S,desc_S,.false.,e_min_down,ED_min_dim,ED_min_down,desc_ED_min,&
             .true.,eta,C_min_down_dim,C_min_down,desc_C_min_down,.false.,.false.,T_dim,T,desc_T,0.0_dp,0,2,2,-1.0_dp,.true.,&
             dealloc,mpi_rank,mpi_size,nprow,order,bs_def,icontxt)
#else
    call omm_sdden_ref(m,n_down,H_down,.true.,S,.false.,e_min_down,ED_min_down,.true.,eta,C_min_down,.false.,.false.,T,0.0_dp,0,2,&
             2,-1.0_dp,.true.,dealloc,mpi_rank)
#endif

    if (mpi_rank==0) then
      print('(a,f21.15,a,f21.15,a)'), 'ED_11 : ', ED_min_up(1,1), ' (up), ', ED_min_down(1,1), ' (down)'
      print('()')
    end if

  end do

  deallocate(T)
  deallocate(C_min_down)
  deallocate(C_min_up)
  deallocate(ED_min_down)
  deallocate(ED_min_up)
  deallocate(D_min_down)
  deallocate(D_min_up)
  deallocate(S)
  deallocate(H_down)
  deallocate(H_up)

#ifdef MPI
  call mpi_finalize(mpi_err)
#endif

end program example5
