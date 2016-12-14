! This code test the pzgemm in scalapack

program pzgemmScaling
  use pspBLAS


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

  logical :: dealloc

  integer :: mpi_err, mpi_size, mpi_rank, m, n, iostat, counti, countf, count_rate
  integer :: nprow, npcol, bs_def_row, bs_def_col, icontxt, iprow, ipcol, niter, i, idxr, idxc
  integer :: H_dim(2), S_dim(2), D_dim(2), info
  integer :: desc_H(9), desc_S(9), desc_D(9)

  complex(dp) :: he, se
  real(dp) :: her,ser, alpha, beta
  real(dp) :: rpt, cpt
  real(dp), allocatable :: H(:,:), S(:,:), D(:,:), Dtrue(:,:)
  complex(dp), allocatable :: zH(:,:), zS(:,:), zD(:,:), zDtrue(:,:)
  real(dp) :: dtime, t0, t1, mm_err, err
  integer*4 timeArray(3)    ! Holds the hour, minute, and second
  complex(dp) :: zalpha, zbeta

  type(psp_matrix_spm) :: Hsp, Ssp, Dsp, Dtruesp ! sparse matrices in pspBLAS
  type(psp_matrix_spm) :: zHsp, zSsp, zDsp, zDtruesp ! sparse matrices in pspBLAS
  character(3) :: fmtH, fmtS, fmtD ! storage type of the sparse matrix, 'coo' or 'csc'
  real(dp) :: thre ! threshold parameter for converting a dense matrix to a sparse amtrix
  integer, allocatable :: idx1(:), idx2(:)! vectors for sparse matrices
  real(dp), allocatable :: val(:)
  complex(dp), allocatable :: zval(:)
  integer :: desc_Hsp(9), desc_Ssp(9), desc_Dsp(9)
  !**********************************************!

  integer, external :: numroc ! it is a function to compute local size

  !**********************************************!
  ! initialization first

  ! initialize information in MPI
  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  ! set up parameters for parallel computing test
  niter=1
  nprow=INT(sqrt(DBLE(mpi_size)))
  npcol=mpi_size/nprow
  order='r' ! important TODO: check how to adapt to different orders
  bs_def_row=10
  bs_def_col=10

  ! initialize information in scalapack
  call blacs_get(-1,0,icontxt)
  call blacs_gridinit(icontxt,order,nprow,npcol)
  ! obtain grid information, working on the (k,l) processor in a grid i by j
  call blacs_gridinfo(icontxt,nprow,npcol,iprow,ipcol)

  ! initialized grid information in pspBLAS
  call psp_gridinit_2D(mpi_comm_world,mpi_size,nprow,order,bs_def_row,bs_def_col,icontxt)

  !*************************************************************************!
  ! generate test matrices

  !***********************************
  if (MPI_rank==0) print *, 'initialize and assign distributed dense matrices in Scalapack'

  if (.true.) then
     ! random matrices
     m=2000 ! global matrix size
     n=1000
     H_dim(1)=numroc(m,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
     H_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
     call descinit(desc_H,m,n,bs_def_row,bs_def_col,0,0,icontxt,H_dim(1),info) ! initialize the descriptor of the global matrix H
     allocate(H(H_dim(1),H_dim(2)))! allocate matrix H
     H=0.0_dp ! assign zero values
     allocate(zH(H_dim(1),H_dim(2)))! allocate matrix H
     zH=cmplx_0 ! assign zero values

     ! initialize and allocate S, D, Dtrue similarly
     S_dim(1)=numroc(m,bs_def_row,iprow,0,nprow)
     S_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol)
     call descinit(desc_S,m,n,bs_def_row,bs_def_col,0,0,icontxt,S_dim(1),info) !
     allocate(S(S_dim(1),S_dim(2)))
     S=0.0_dp
     allocate(zS(S_dim(1),S_dim(2)))
     zS=cmplx_0

     D_dim(1)=numroc(m,bs_def_row,iprow,0,nprow)
     D_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol)
     call descinit(desc_D,m,n,bs_def_row,bs_def_col,0,0,icontxt,D_dim(1),info) !
     allocate(D(D_dim(1),D_dim(2)))
     D=0.0_dp
     allocate(zD(D_dim(1),D_dim(2)))
     zD=cmplx_0

     allocate(Dtrue(D_dim(1),D_dim(2)))
     Dtrue=0.0_dp
     allocate(zDtrue(D_dim(1),D_dim(2)))
     zDtrue=cmplx_0

     ! generate random matrices
     if (.true.) then
        call init_random_seed()
        call RANDOM_NUMBER(H)
        call RANDOM_NUMBER(S)
        do idxc=1,n
           do idxr=1,m
              call RANDOM_NUMBER(rpt)
              call RANDOM_NUMBER(cpt)
              he=CMPLX(rpt,cpt)
              call pzelset(zH,idxr,idxc,desc_H,he)
           end do
        end do
        do idxc=1,n
           do idxr=1,m
              call RANDOM_NUMBER(rpt)
              call RANDOM_NUMBER(cpt)
              se=CMPLX(rpt,cpt)
              call pzelset(zS,idxr,idxc,desc_S,se)
           end do
        end do
     end if
  end if

  !***********************************
  if (MPI_rank==0) print *, 'initialize and assign distributed sparse matrices in pspBLAS'

  if (MPI_rank==0) print *, 'test sparse summation in coo format'
  fmtH='coo' ! specify storage format, 'coo' or 'csc'
  fmtS='coo'
  fmtD='coo'
  alpha=1.5_dp
  beta=1.2_dp
  zalpha=(1.2_dp,1.0_dp)
  zbeta=(1.0_dp,2.5_dp)

  ! first method to generate a sparse matrix: thresholding a dense matrix in MatrixSwitch
  thre = 0.5_dp
  call psp_den2sp_m(H,desc_H,Hsp,fmtH,thre)   ! need to initialize all matrices
  call psp_den2sp_m(S,desc_S,Ssp,fmtS,thre)   ! need to initialize all matrices
  call psp_den2sp_m(zH,desc_H,zHsp,fmtH,thre) ! need to initialize all matrices
  call psp_den2sp_m(zS,desc_S,zSsp,fmtS,thre) ! need to initialize all matrices
  call psp_spm_zeros(Dsp,m,n,fmtD,.true.)     ! need to initialize all matrices
  call psp_spm_zeros(zDsp,m,n,fmtD,.false.)   ! need to initialize all matrices

  do idxc=1,H_dim(2)
     do idxr=1,H_dim(1)
        if (abs(H(idxr,idxc))<thre) H(idxr,idxc) = 0.0_dp
        if (abs(S(idxr,idxc))<thre) S(idxr,idxc) = 0.0_dp
        if (abs(zH(idxr,idxc))<thre) zH(idxr,idxc) = cmplx_0
        if (abs(zS(idxr,idxc))<thre) zS(idxr,idxc) = cmplx_0

        Dtrue(idxr,idxc)=alpha*H(idxr,idxc)+beta*S(idxr,idxc)
        zDtrue(idxr,idxc)=zalpha*zH(idxr,idxc)+zbeta*zS(idxr,idxc)
     end do
  end do

  if (MPI_rank==0) print *, 'begin spm+m'
  call psp_sum_spmm(Hsp,S,D,alpha,beta)
  call psp_sum_spmm(zHsp,zS,zD,zalpha,zbeta)

  if (MPI_rank==0) print *, 'check results'
  mm_err=MAXVAL(abs(Dtrue-D))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (MPI_rank==0) print *, 'real case              ', err
  mm_err=MAXVAL(abs(zDtrue-zD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (MPI_rank==0) print *, 'complex case           ', err


  if (MPI_rank==0) print *, 'begin spm+spm'
  call psp_sum_spmspm(Hsp,Ssp,Dsp,alpha,beta)
  call psp_sum_spmspm(zHsp,zSsp,zDsp,zalpha,zbeta)

  ! sparse to dense
  call psp_sp2den_m(Dsp,D,desc_D)
  call psp_sp2den_m(zDsp,zD,desc_D)

  if (MPI_rank==0) print *, 'check results'
  mm_err=MAXVAL(abs(Dtrue-D))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (MPI_rank==0) print *, 'real case              ', err
  mm_err=MAXVAL(abs(zDtrue-zD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (MPI_rank==0) print *, 'complex case           ', err

  if (MPI_rank==0) print *, 'test sparse summation in csc format'

  fmtH='csc' ! specify storage format, 'coo' or 'csc'
  fmtS='csc'
  fmtD='csc'

  ! first method to generate a sparse matrix: thresholding a dense matrix in MatrixSwitch
  thre = 0.9_dp
  call psp_den2sp_m(H,desc_H,Hsp,fmtH,thre)   ! need to initialize all matrices
  call psp_den2sp_m(S,desc_S,Ssp,fmtS,thre)   ! need to initialize all matrices
  call psp_den2sp_m(zH,desc_H,zHsp,fmtH,thre) ! need to initialize all matrices
  call psp_den2sp_m(zS,desc_S,zSsp,fmtS,thre) ! need to initialize all matrices
  call psp_spm_zeros(Dsp,m,n,fmtD,.true.)     ! need to initialize all matrices
  call psp_spm_zeros(zDsp,m,n,fmtD,.false.)   ! need to initialize all matrices

  !if (MPI_rank==0) print *, 'thresholding'

  do idxc=1,H_dim(2)
     do idxr=1,H_dim(1)
        if (abs(H(idxr,idxc))<thre) H(idxr,idxc) = 0.0_dp
        if (abs(S(idxr,idxc))<thre) S(idxr,idxc) = 0.0_dp
        if (abs(zH(idxr,idxc))<thre) zH(idxr,idxc) = cmplx_0
        if (abs(zS(idxr,idxc))<thre) zS(idxr,idxc) = cmplx_0

        Dtrue(idxr,idxc)=alpha*H(idxr,idxc)+beta*S(idxr,idxc)
        zDtrue(idxr,idxc)=zalpha*zH(idxr,idxc)+zbeta*zS(idxr,idxc)
     end do
  end do


  if (.true.) then
     if (MPI_rank==0) print *, 'begin spm+m'
     call psp_sum_spmm(Hsp,S,D,alpha,beta)
     call psp_sum_spmm(zHsp,zS,zD,zalpha,zbeta)

     if (MPI_rank==0) print *, 'check results'
     mm_err=MAXVAL(abs(Dtrue-D))
     err=0.0_dp
     call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
     if (MPI_rank==0) print *, 'real case              ', err
     mm_err=MAXVAL(abs(zDtrue-zD))
     err=0.0_dp
     call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
     if (MPI_rank==0) print *, 'complex case           ', err
  end if


  if (MPI_rank==0) print *, 'begin spm+spm'
  call psp_sum_spmspm(Hsp,Ssp,Dsp,alpha,beta)
  call psp_sum_spmspm(zHsp,zSsp,zDsp,zalpha,zbeta)
if (.true.) then
  ! sparse to dense
  call psp_sp2den_m(Dsp,D,desc_D)
  call psp_sp2den_m(zDsp,zD,desc_D)

  if (MPI_rank==0) print *, 'check results'
  mm_err=MAXVAL(abs(Dtrue-D))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (MPI_rank==0) print *, 'real case              ', err
  mm_err=MAXVAL(abs(zDtrue-zD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (MPI_rank==0) print *, 'complex case           ', err
else
call psp_sp2den_m(Dsp,D,desc_D)

if (MPI_rank==0) print *, 'check results'
mm_err=MAXVAL(abs(Dtrue-D))
err=0.0_dp
call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
if (MPI_rank==0) print *, 'real case              ', err


call psp_csc2coo(Dsp)
call psp_sp2den_m(Dsp,D,desc_D)

if (MPI_rank==0) print *, 'check results'
mm_err=MAXVAL(abs(Dtrue-D))
err=0.0_dp
call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
if (MPI_rank==0) print *, 'real case              ', err

call psp_coo2csc(Dsp)

! sparse to dense
call psp_sp2den_m(Dsp,D,desc_D)
call psp_sp2den_m(zDsp,zD,desc_D)

if (MPI_rank==0) print *, 'check results'
mm_err=MAXVAL(abs(Dtrue-D))
err=0.0_dp
call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
if (MPI_rank==0) print *, 'real case              ', err
mm_err=MAXVAL(abs(zDtrue-zD))
err=0.0_dp
call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
if (MPI_rank==0) print *, 'complex case           ', err
end if



  deallocate(H)
  deallocate(S)
  deallocate(D)
  deallocate(Dtrue)
  deallocate(zH)
  deallocate(zS)
  deallocate(zD)
  deallocate(zDtrue)

  call psp_deallocate_spm(Hsp)
  call psp_deallocate_spm(Ssp)
  call psp_deallocate_spm(Dsp)
  call psp_deallocate_spm(zHsp)
  call psp_deallocate_spm(zSsp)
  call psp_deallocate_spm(zDsp)

  call mpi_finalize(mpi_err)

end program pzgemmScaling
