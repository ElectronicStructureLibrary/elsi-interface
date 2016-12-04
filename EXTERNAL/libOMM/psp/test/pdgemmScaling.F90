! This code test the pdgemm in pspBLAS

program pdgemmScaling
  use pspBLAS


  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  !**** VARIABLES *******************************!

  character(1) :: order
  character(21) :: file_name

  logical :: dealloc

  integer :: mpi_err, mpi_size, mpi_rank, m, n, k, iostat, counti, countf, count_rate
  integer :: nprow, npcol, bs_def_row, bs_def_col, icontxt, iprow, ipcol, niter, i, idxr, idxc
  integer :: H_dim(2), S_dim(2), D_dim(2), Ht_dim(2), St_dim(2), Dt3_dim(2), info
  integer :: desc_H(9), desc_S(9), desc_D(9), desc_Ht(9), desc_St(9)

  real(dp) :: he, se, t0, t1, mm_err, err
  real(dp), allocatable :: H(:,:), S(:,:), D(:,:), DD(:,:), Ht(:,:), St(:,:)
  real(dp)    :: dtime, alpha, beta
  integer*4 timeArray(3)    ! Holds the hour, minute, and second

  type(psp_matrix_spm) :: Hsp, Ssp, Dsp, Htsp, Stsp ! sparse matrices in pspBLAS
  character(3) :: spm_storage ! storage type of the sparse matrix, 'coo' or 'csc'
  real(dp) :: thre ! threshold parameter for converting a dense matrix to a sparse amtrix
  integer, allocatable :: idx1(:), idx2(:)! vectors for sparse matrices
  real(dp), allocatable :: val(:)
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
  npcol=nprow
  order='r' ! important TODO: check how to adapt to different orders
  bs_def_row=100
  bs_def_col=100

  ! initialize information in scalapack
  call blacs_get(-1,0,icontxt)
  call blacs_gridinit(icontxt,order,nprow,npcol)
  ! obtain grid information, working on the (k,l) processor in a grid i by j
  call blacs_gridinfo(icontxt,nprow,npcol,iprow,ipcol)

  ! initialized grid information in pspBLAS
  call psp_gridinit_2D(mpi_size,nprow,order,bs_def_row,bs_def_col,icontxt)
  !**********************************************!
  ! open files for read and write

  open(11,file='dlap.txt',action='write',position='append') ! storing running time and errors
  open(13,file='dgemm.txt',action='write',position='append')
  open(14,file='dspmm.txt',action='write',position='append')
  open(15,file='dmspm.txt',action='write',position='append')
  open(16,file='dspmspm.txt',action='write',position='append')
5 format (a,i4,a,f16.10,a,f16.10)

  !*************************************************************************!
  ! generate test matrices

  !***********************************
  ! initialize and assign distributed dense matrices in Scalapack

  if (.true.) then
     ! random matrices
     m=1000 ! global matrix size
     n=1020
     k=833

     alpha = 1.5_dp
     beta = 0.5_dp

     H_dim(1)=numroc(m,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
     H_dim(2)=numroc(k,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
     call descinit(desc_H,m,k,bs_def_row,bs_def_col,0,0,icontxt,H_dim(1),info) ! initialize the descriptor of the global matrix H
     allocate(H(H_dim(1),H_dim(2)))! allocate matrix H
     H=0.0_dp ! assign zero values

     Ht_dim(1)=numroc(k,bs_def_row,iprow,0,nprow)
     Ht_dim(2)=numroc(m,bs_def_col,ipcol,0,npcol)
     call descinit(desc_Ht,k,m,bs_def_row,bs_def_col,0,0,icontxt,Ht_dim(1),info) ! Ht=transpose(H)
     allocate(Ht(Ht_dim(1),Ht_dim(2)))
     Ht=0.0_dp

     ! initialize and allocate S and D similarly
     S_dim(1)=numroc(k,bs_def_row,iprow,0,nprow)
     S_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol)
     call descinit(desc_S,k,n,bs_def_row,bs_def_col,0,0,icontxt,S_dim(1),info) !
     allocate(S(S_dim(1),S_dim(2)))
     S=0.0_dp

     St_dim(1)=numroc(n,bs_def_row,iprow,0,nprow)
     St_dim(2)=numroc(k,bs_def_col,ipcol,0,npcol)
     call descinit(desc_St,n,k,bs_def_row,bs_def_col,0,0,icontxt,St_dim(1),info) ! St=transpose(S)
     allocate(St(St_dim(1),St_dim(2)))
     St=0.0_dp

     D_dim(1)=numroc(m,bs_def_row,iprow,0,nprow)
     D_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol)
     call descinit(desc_D,m,n,bs_def_row,bs_def_col,0,0,icontxt,D_dim(1),info) !
     allocate(D(D_dim(1),D_dim(2)))
     D=0.0_dp
     allocate(DD(D_dim(1),D_dim(2)))
     DD=0.0_dp

     ! generate random matrices
     call init_random_seed()
     !call random_seed()
     call RANDOM_NUMBER(H)
     call RANDOM_NUMBER(Ht)

     call RANDOM_NUMBER(S)
     call RANDOM_NUMBER(St)

  end if

  !***********************************
  ! initialize and assign distributed sparse matrices in pspBLAS

  if (.false.) then
     ! test 'coo' format
     spm_storage='coo' ! specify storage format, 'coo' or 'csc'

     ! first method to generate a sparse matrix: thresholding a dense matrix in MatrixSwitch
     thre = 0.99_dp
     call psp_den2sp_m(H,desc_H,Hsp,spm_storage,thre)
     call psp_den2sp_m(S,desc_S,Ssp,spm_storage,thre)
     call psp_den2sp_m(Ht,desc_Ht,Htsp,spm_storage,thre)
     call psp_den2sp_m(St,desc_St,Stsp,spm_storage,thre)

     ! convert the sparse matrix back to dense and use them as comparison
     call psp_sp2den_m(Hsp,H,desc_H)
     call psp_sp2den_m(Ssp,S,desc_S)
     call psp_sp2den_m(Htsp,Ht,desc_Ht)
     call psp_sp2den_m(Stsp,St,desc_St)

     ! second method to generate a sparse matrix: use the information of the Sparse Triplet (ST) format (in this example, it's coo)
     allocate(idx1(Hsp%nnz))
     allocate(idx2(Hsp%nnz))
     allocate(val(Hsp%nnz))
     idx1=Hsp%row_ind ! row indices
     idx2=Hsp%col_ind ! column indices if 'coo', or pointers to position in row indices if 'csc'
     val=Hsp%dval     ! nonzero entries
     desc_Hsp=Hsp%desc ! use the same descriptor as the one in Scalapack or MatrixSwitch
     call psp_register_spm(Hsp,idx1,idx2,val,desc_Hsp,spm_storage,H_dim,nprow,npcol) ! assign a sparse matrix
  else
     ! test 'csc' format
     spm_storage='csc' ! specify storage format, 'coo' or 'csc'

     ! first method to generate a sparse matrix: thresholding a dense matrix in MatrixSwitch
     thre = 0.8_dp
     call psp_den2sp_m(H,desc_H,Hsp,spm_storage,thre)
     call psp_den2sp_m(S,desc_S,Ssp,spm_storage,thre)
     call psp_den2sp_m(Ht,desc_Ht,Htsp,spm_storage,thre)
     call psp_den2sp_m(St,desc_St,Stsp,spm_storage,thre)

     ! convert the sparse matrix back to dense and use them as comparison
     call psp_sp2den_m(Hsp,H,desc_H)
     call psp_sp2den_m(Ssp,S,desc_S)
     call psp_sp2den_m(Htsp,Ht,desc_Ht)
     call psp_sp2den_m(Stsp,St,desc_St)

     ! second method to generate a sparse matrix: use the information of the Sparse Triplet (ST) format (in this example, it's csc)
     allocate(idx1(Hsp%nnz))
     allocate(idx2(Hsp%loc_dim2+1))
     allocate(val(Hsp%nnz))
     idx1=Hsp%row_ind ! row indices
     idx2=Hsp%col_ptr ! column indices if 'coo', or pointers to position in row indices if 'csc'
     val=Hsp%dval     ! nonzero entries
     desc_Hsp=Hsp%desc ! use the same descriptor as the one in Scalapack or MatrixSwitch
     call psp_register_spm(Hsp,idx1,idx2,val,desc_Hsp,spm_storage,H_dim,nprow,npcol) ! assign a sparse matrix
  end if

  if (MPI_rank==0) print *, 'process grid', nprow, npcol
  if (MPI_rank==0) print *, 'block size', bs_def_row, bs_def_col

  !************************************************************************!
  if (mpi_rank==0) print *,  'Begin n n'
  t0 = MPI_Wtime()
  do i=1,niter
     call pdgemm('n','n',m,n,k,alpha,H,1,1,desc_H,S,1,1,desc_S,beta,D,1,1,desc_D)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of Scalapack = n n ', dtime
  if (mpi_rank==0)  write(11,5) 'n  ', mpi_size, ' type   nn   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemm(m,n,k,H,'n',S,'n',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gemm = n n ', dtime, 'error of psp_gemm: n n ', err
  if (mpi_rank==0)  write(13,5) 'n  ', mpi_size, ' type   nn   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gespmm(m,n,k,Hsp,'n',S,'n',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gespmm = n n ', dtime, 'error of psp_gespmm: n n ', err
  if (mpi_rank==0)  write(14,5) 'n  ', mpi_size, ' type   nn   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemspm(m,n,k,H,'n',Ssp,'n',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time psp_gemspm = n n ', dtime, 'error of psp_gemspm: n n ', err
  if (mpi_rank==0)  write(15,5) 'n  ', mpi_size, ' type   nn   time  ', dtime, '  err  ', err

  ! initialize the sparse solution
  call psp_spm_zeros(Dsp,m,n,spm_storage,.true.) ! need to initialize a sparse matrix before we use it
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gespmspm(m,n,k,Hsp,'n',Ssp,'n',Dsp,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  ! sparse to dense
  call psp_sp2den_m(Dsp,DD,desc_D)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gespmspm = n n ', dtime, 'error of psp_gespmspm: n n ', err
  if (mpi_rank==0)  write(16,5) 'n  ', mpi_size, ' type   nn   time  ', dtime, '  err  ', err

  !************************************************************************!
  if (mpi_rank==0) print *,  'Begin n t'
  D=0.0_dp
  t0=MPI_Wtime()
  do i=1,niter
     call pdgemm('n','t',m,n,k,alpha,H,1,1,desc_H,St,1,1,desc_St,beta,D,1,1,desc_D)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of Scalapack = n t ', dtime
  if (mpi_rank==0)  write(11,5) 'n  ', mpi_size, ' type   nt   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemm(m,n,k,H,'n',St,'t',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gemm = n t ', dtime, 'error of psp_gemm: n t ', err
  if (mpi_rank==0)  write(13,5) 'n  ', mpi_size, ' type   nt   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gespmm(m,n,k,Hsp,'n',St,'t',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gespmm = n t ', dtime, 'error of psp_gespmm: n t ', err
  if (mpi_rank==0)  write(14,5) 'n  ', mpi_size, ' type   nt   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemspm(m,n,k,H,'n',Stsp,'t',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time psp_gemspm = n t ', dtime, 'error of psp_gemspm: n t ', err
  if (mpi_rank==0)  write(15,5) 'n  ', mpi_size, ' type   nt   time  ', dtime, '  err  ', err

  ! initialize the sparse solution
  call psp_spm_zeros(Dsp,m,n,spm_storage,.true.) ! need to initialize a sparse matrix before we use it
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gespmspm(m,n,k,Hsp,'n',Stsp,'t',Dsp,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  ! sparse to dense
  call psp_sp2den_m(Dsp,DD,desc_D)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gespmspm = n t ', dtime, 'error of psp_gespmspm: n t ', err
  if (mpi_rank==0)  write(16,5) 'n  ', mpi_size, ' type   nt   time  ', dtime, '  err  ', err


  !************************************************************************!
  if (mpi_rank==0) print *,  'Begin t n'
  D=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call pdgemm('t','n',m,n,k,alpha,Ht,1,1,desc_Ht,S,1,1,desc_S,beta,D,1,1,desc_D)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of Scalapack = t n ', dtime
  if (mpi_rank==0)  write(11,5) 'n  ', mpi_size, ' type   tn   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemm(m,n,k,Ht,'t',S,'n',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gemm = t n ', dtime, 'error of psp_gemm: t n ', err
  if (mpi_rank==0)  write(13,5) 'n  ', mpi_size, ' type   tn   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gespmm(m,n,k,Htsp,'t',S,'n',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gespmm = t n ', dtime, 'error of psp_gespmm: t n ', err
  if (mpi_rank==0)  write(14,5) 'n  ', mpi_size, ' type   tn   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemspm(m,n,k,Ht,'t',Ssp,'n',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time psp_gemspm = t n ', dtime, 'error of psp_gemspm: t n ', err
  if (mpi_rank==0)  write(15,5) 'n  ', mpi_size, ' type   tn   time  ', dtime, '  err  ', err

  ! initialize the sparse solution
  call psp_spm_zeros(Dsp,m,n,spm_storage,.true.) ! need to initialize a sparse matrix before we use it
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gespmspm(m,n,k,Htsp,'t',Ssp,'n',Dsp,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  ! sparse to dense
  call psp_sp2den_m(Dsp,DD,desc_D)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gespmspm = t n ', dtime, 'error of psp_gespmspm: t n ', err
  if (mpi_rank==0)  write(16,5) 'n  ', mpi_size, ' type   tn   time  ', dtime, '  err  ', err


  !************************************************************************!
  if (mpi_rank==0) print *,  'Begin t t'
  D=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call pdgemm('t','t',m,n,k,alpha,Ht,1,1,desc_Ht,St,1,1,desc_St,beta,D,1,1,desc_D)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  if (mpi_rank==0) print *, 'time of Scalapack = t t ', dtime
  if (mpi_rank==0)  write(11,5) 'n  ', mpi_size, ' type   tt   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemm(m,n,k,Ht,'t',St,'t',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gemm = t t ', dtime, 'error of psp_gemm: t t ', err
  if (mpi_rank==0)  write(13,5) 'n  ', mpi_size, ' type   tt   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gespmm(m,n,k,Htsp,'t',St,'t',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time of psp_gespmm = t t ', dtime, 'error of psp_gespmm: t t ', err
  if (mpi_rank==0)  write(14,5) 'n  ', mpi_size, ' type   tt   time  ', dtime, '  err  ', err

  DD=0.0_dp
  t0 = MPI_Wtime()
  do i=1,niter
     call psp_gemspm(m,n,k,Ht,'t',Stsp,'t',DD,alpha,beta)
  enddo
  t1 = MPI_Wtime()
  dtime=(t1-t0)/DBLE(niter)
  mm_err=MAXVAL(abs(D-DD))
  err=0.0_dp
  call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE, MPI_SUM, 0, mpi_comm_world, mpi_err)
  if (mpi_rank==0) print *, 'time psp_gemspm = t t ', dtime, 'error of psp_gemspm: t t ', err
  if (mpi_rank==0)  write(15,5) 'n  ', mpi_size, ' type   tt   time  ', dtime, '  err  ', err

  close(11)
  close(13)
  close(14)
  close(15)
  close(16)

  deallocate(D)
  deallocate(DD)
  deallocate(S)
  deallocate(St)
  deallocate(H)
  deallocate(Ht)
  call psp_deallocate_spm(Hsp)
  call psp_deallocate_spm(Ssp)
  call psp_deallocate_spm(Htsp)
  call psp_deallocate_spm(Stsp)
  call psp_deallocate_spm(Dsp)

  call mpi_finalize(mpi_err)

end program pdgemmScaling
