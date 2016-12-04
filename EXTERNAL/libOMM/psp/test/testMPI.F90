program testMPI
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

  logical :: dealloc, isAlloc

  integer :: mpi_err, mpi_size, mpi_rank, m, n, iostat, counti, countf, count_rate
  integer :: nprow, npcol, bs_def_row, bs_def_col, icontxt, iprow, ipcol, niter, i, idxr, idxc
  integer :: H_dim(2), info, sizeidx1,sizeidx2,sizeval,sizeint, sizebuf, position, idx_loc, iprowNew
  integer, allocatable :: buf(:)
  integer :: desc_H(9)
  integer :: tempIP, numStep, cnt, cnt1, len, dest, src, tag=0, requestInt, requestSpm
  integer :: content, contentLoc, lenLoc, lenGlb, STATUS(MPI_STATUS_SIZE)
  complex(dp) :: he, se
  real(dp), allocatable :: H(:,:), Hrd(:,:)
  real(dp) :: dtime, t0, t1, mm_err, err
  type(psp_MPI_dspm) :: C_loc_MPI, CC_loc_MPI
  integer, allocatable :: C_loc_idx1(:), C_loc_idx2(:)
  real(dp), allocatable :: C_loc_val(:)
  integer :: psp_MPI_datatype_spm
  type(dList), pointer :: list, lastElem
  integer, dimension(6) :: tp=(/ MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_DOUBLE_PRECISION /)
  integer, dimension(6) :: addPos, disp, blockLen

  type(psp_matrix_spm) :: Hsp, Hrdsp ! sparse matrices in pspBLAS
  character(3) :: fmtH ! storage type of the sparse matrix, 'coo' or 'csc'
  real(dp) :: thre, temp ! threshold parameter for converting a dense matrix to a sparse amtrix
  integer, allocatable :: idx1(:), idx2(:)! vectors for sparse matrices
  real(dp), allocatable :: val(:)
  !**********************************************!

  integer, external :: numroc ! it is a function to compute local size

  !**** GLOBAL **********************************!
#ifdef MPI
  character(1) :: psp_proc_order

  integer :: psp_mpi_comm_world
  integer :: psp_mpi_size
  integer :: psp_nprow
  integer :: psp_npcol
  integer :: psp_bs_def_row
  integer :: psp_bs_def_col
  integer :: psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
  integer :: psp_bs_num
  integer :: psp_icontxt ! BLACS context handle used by psp
  integer :: psp_mpi_comm_cart, psp_mpi_comm_row, psp_mpi_comm_col


  common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
  common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
  common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
  common /psp_grid2D/ psp_bs_num
  common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
  common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row,  psp_mpi_comm_col
#endif

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
  call psp_gridinit_2D(mpi_size,nprow,order,bs_def_row,bs_def_col,icontxt)

  !*************************************************************************!
  ! generate test matrices

  !***********************************
  if (MPI_rank==0) print *, 'Construct H and Hsp'

  ! random matrices
  m=nprow*40  ! global matrix size
  n=npcol*40 ! guarantee that data to be passed have the same size
  if (MPI_rank==0) print *, 'm,n ',m, n
  H_dim(1)=numroc(m,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
  H_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
  call descinit(desc_H,m,n,bs_def_row,bs_def_col,0,0,icontxt,H_dim(1),info) ! initialize the descriptor of the global matrix H
  allocate(H(H_dim(1),H_dim(2)))! allocate matrix H
  H=0.0_dp ! assign zero values
  allocate(Hrd(H_dim(1),H_dim(2)))! allocate matrix H
  Hrd=0.0_dp ! assign zero values

  ! if (MPI_rank==0) print *, 'H_dim', H_dim, shape(H), 1

  ! generate random matrices
  call init_random_seed()
  call RANDOM_NUMBER(H)
  ! H=0.0_dp
  ! H(1,1)=1.0_dp

  fmtH='csc'
  thre = 0.5_dp
  call psp_den2sp_m(H,desc_H,Hsp,fmtH,thre)
  call psp_spm_zeros(Hrdsp,m,n,'csc',.true.)

  do idxc=1,H_dim(2)
     do idxr=1,H_dim(1)
        if (abs(H(idxr,idxc))<thre) H(idxr,idxc) = 0.0_dp
     end do
  end do


  !***********************************
  ! MPI_REDUCE to compute summation of scalars
  if (.false.) then
     if (MPI_rank==0) print *, 'H_dim', H_dim, shape(H), 2
     content=iprow
     tempIP=iprow
     numStep=0
     temp=real(nprow)
     do while (temp>1)
        numStep=numStep+1
        temp=temp/2
     end do
     do cnt=1,numStep
        if (MODULO(iprow,2**(cnt-1))==0) then
           ! get uniform data size

           ! send and receive A_loc
           if (MODULO(tempIP,2)==1) then
              ! set up A_loc and send
              contentLoc=content
              dest=iprow-2**(cnt-1)
              !if (ipcol==0) print *, 'dest', tempIP,dest,iprow,cnt,numStep
              call MPI_SEND(contentLoc,1,MPI_INT,dest,tag,psp_mpi_comm_col,mpi_err)
           else
              ! receive A_loc and and to A
              src=iprow+2**(cnt-1)
              if (src<nprow) then
                 !if (ipcol==0) print *, 'src', tempIP,src,iprow,cnt,numStep
                 call MPI_RECV(contentLoc,1,MPI_INT, src,tag,psp_mpi_comm_col,mpi_err)
                 content=content+contentLoc
              end if
           end if
           tempIP=tempIP/2
        end if
     end do
     if (iprow==0) print *, (0+nprow-1)*nprow/2,content,iprow
     ! TODO: Why H_dim is changed?

     if (MPI_rank==0) print *, 'H_dim', H_dim, shape(H), 3
  end if
  !***********************************
  ! MPI_REDUCE to compute summation of matrices
  idx_loc=0

  call MPI_REDUCE(H, Hrd, H_dim(1)*H_dim(2), MPI_DOUBLE_PRECISiON, MPI_SUM, 0, psp_mpi_comm_col, mpi_err)
  !if (MPI_rank==0) print *, H(1:5,1:5)
  if (MPI_rank==0) print *, 'true'
  if (MPI_rank==0) print *, Hrd(1:5,1:5)

  Hrd=0.0_dp
  if (.true.) then
     tempIP=iprow
     numStep=0
     temp=real(nprow)
     do while (temp>1)
        numStep=numStep+1
        temp=temp/2
     end do
     do cnt=1,numStep
        if (MODULO(iprow,2**(cnt-1))==0) then
           ! get uniform data size

           ! send and receive A_loc
           if (MODULO(tempIP,2)==1) then
              ! set up A_loc and send
              Hrd=H
              dest=iprow-2**(cnt-1)
              H_dim(1)=numroc(m,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
              H_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
              !if (ipcol==0) print *, 'H_dim', H_dim, shape(H), 'dest', tempIP,dest,iprow,cnt,numStep
              !if (ipcol==0) print *, 'dest', tempIP,dest,iprow,cnt,numStep
              call MPI_SEND(Hrd,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISiON,dest,tag,psp_mpi_comm_col,mpi_err)
           else
              ! receive A_loc and and to A
              src=iprow+2**(cnt-1)
              if (src<nprow) then
                 H_dim(1)=numroc(m,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
                 H_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
                 !if (ipcol==0) print *, 'src', tempIP,src,iprow,cnt,numStep
                 call MPI_RECV(Hrd,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISiON, src,tag,psp_mpi_comm_col,mpi_err)
                 H=H+Hrd
              end if
           end if
           tempIP=tempIP/2
        end if
     end do
     H_dim(1)=numroc(m,bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
     H_dim(2)=numroc(n,bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
     if (iprow==0) then
        call MPI_ISEND(H,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISION,idx_loc,tag,psp_mpi_comm_col,requestSpm,mpi_err)
        call MPI_WAIT(requestSpm,STATUS,mpi_err)
     end if
     if (iprow==idx_loc) then
        call MPI_IRECV(H,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISION,0,tag,psp_mpi_comm_col,requestSpm,mpi_err)
        call MPI_WAIT(requestSpm,STATUS,mpi_err)
     end if
  else
     call psp_MPI_REDUCE_den(H,m,n)
  end if
  if (MPI_rank==0) print *, 'finish dense'
  call MPI_Barrier(psp_mpi_comm_world,mpi_err)
  !***********************************
  ! MPI_REDUCE to compute summation of sparse matrices using user defined MPI structure

  if (.false.) then
     tempIP=iprow
     numStep=0
     temp=real(nprow)
     do while (temp>1)
        numStep=numStep+1
        temp=temp/2
     end do
     do cnt=1,numStep
        if (MODULO(iprow,2**(cnt-1))==0) then
           ! get uniform data size
           if (MODULO(tempIP,2)==1) then
              ! set up lenLoc and send
              lenLoc=Hsp%nnz
              dest=iprow-2**(cnt-1)
              src=iprow
              call MPI_ISEND(lenLoc,1,MPI_INT,dest,tag,psp_mpi_comm_col,requestInt,mpi_err)
              call MPI_WAIT(requestInt,STATUS,mpi_err)
           else
              ! receive lenLoc and update
              dest=iprow
              src=iprow+2**(cnt-1)
              if (src<nprow) then
                 call MPI_IRECV(lenLoc,1,MPI_INT, src,tag,psp_mpi_comm_col,requestInt,mpi_err)!need to use different tags to distinguish messages
                 call MPI_WAIT(requestInt,STATUS,mpi_err)
              end if
           end if

           if (src<nprow) then
              C_loc_MPI%nnz=Hsp%nnz
              C_loc_MPI%loc_dim1=Hsp%loc_dim1
              C_loc_MPI%loc_dim2=Hsp%loc_dim2
              if (allocated(C_loc_idx1)) deallocate(C_loc_idx1)
              if (allocated(C_loc_idx2)) deallocate(C_loc_idx2)
              if (allocated(C_loc_val)) deallocate(C_loc_val)
              allocate(C_loc_idx1(lenLoc))
              allocate(C_loc_idx2(C_loc_MPI%loc_dim2+1))
              allocate(C_loc_val(lenLoc))
              if (MODULO(tempIP,2)==1) then
                 if (Hsp%nnz>0) then
                    do cnt1=1,lenLoc
                       C_loc_idx1(cnt1)=Hsp%row_ind(cnt1)
                       C_loc_val(cnt1)=Hsp%dval(cnt1)
                    end do
                 end if
                 do cnt1=1,C_loc_MPI%loc_dim2+1
                    C_loc_idx2(cnt1)=Hsp%col_ptr(cnt1)
                 end do
              end if

              ! define new MPI datatype
              call MPI_GET_ADDRESS(C_loc_MPI%nnz, disp(1), mpi_err)
              call MPI_GET_ADDRESS(C_loc_MPI%loc_dim1, disp(2), mpi_err)
              call MPI_GET_ADDRESS(C_loc_MPI%loc_dim2, disp(3), mpi_err)
              call MPI_GET_ADDRESS(C_loc_idx1, disp(4), mpi_err)
              call MPI_GET_ADDRESS(C_loc_idx2, disp(5), mpi_err)
              call MPI_GET_ADDRESS(C_loc_val, disp(6), mpi_err)
              disp(2)=disp(2)-disp(1)
              disp(3)=disp(3)-disp(1)
              disp(4)=disp(4)-disp(1)
              disp(5)=disp(5)-disp(1)
              disp(6)=disp(6)-disp(1)
              disp(1)=0
              blockLen=1
              blockLen(4)=lenLoc
              blockLen(6)=lenLoc
              blockLen(5)=C_loc_MPI%loc_dim2+1
              call MPI_Type_struct(6, blockLen, disp, tp, psp_MPI_datatype_spm, mpi_err)
              call MPI_Type_commit(psp_MPI_datatype_spm, mpi_err)
           end if

           ! send and receive A_loc
           if (MODULO(tempIP,2)==1) then
              ! set up A_loc and send
              dest=iprow-2**(cnt-1)
              call MPI_ISEND(C_loc_MPI,1,psp_MPI_datatype_spm,dest,tag,psp_mpi_comm_col,requestSpm,mpi_err)
              call MPI_WAIT(requestSpm,STATUS,mpi_err)
           else
              ! receive A_loc and and to A
              src=iprow+2**(cnt-1)
              if (src<nprow) then
                 call MPI_IRECV(C_loc_MPI,1,psp_MPI_datatype_spm, src,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                 call MPI_WAIT(requestSpm,STATUS,mpi_err)
                 ! Hsp=C_loc_MPI+Hsp
                 if (C_loc_MPI%nnz>0) then
                    call psp_sst_sum_spmspm(Hsp%loc_dim1,Hsp%loc_dim2,1.0_dp,&
                         C_loc_idx1,C_loc_idx2,C_loc_val,'csc',1.0_dp,&
                         Hsp%row_ind,Hsp%col_ptr,Hsp%dval,'csc',Hsp%row_ind,Hsp%col_ptr,Hsp%dval,&
                         'csc',C_loc_MPI%nnz,Hsp%nnz)
                 end if
                 Hsp%nnz=size(Hsp%dval)
              end if
           end if
           tempIP=tempIP/2
        end if
        !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
     end do
  else
  end if
  if (.false.) then
     !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
     call psp_MPI_REDUCE_dspm_struct(Hsp, 0, .false.)
  end if
  !call MPI_Barrier(psp_mpi_comm_world,mpi_err)


  !***********************************
  ! MPI_REDUCE to compute summation of sparse matrices using user defined MPI_pack and MPI_unpack
  ! This method should be faster than user defined structure because each structure is used only once.
  ! If each structure is used repeatedly, then user defined structure is faster.


  if (.false.) then
     iprowNew = MODULO(iprow-idx_loc,nprow)
     tempIP=iprowNew
     numStep=0
     temp=real(nprow)
     do while (temp>1)
        numStep=numStep+1
        temp=temp/2
     end do
     do cnt=1,numStep
        if (MODULO(iprowNew,2**(cnt-1))==0) then
           ! get uniform data size
           if (MODULO(tempIP,2)==1) then
              ! set up lenLoc and send
              lenLoc=Hsp%nnz
              dest=MODULO(iprowNew-2**(cnt-1)+idx_loc,nprow)
              src=MODULO(iprowNew+idx_loc,nprow)
              isAlloc=.true.
              call MPI_ISEND(lenLoc,1,MPI_INT,dest,tag,psp_mpi_comm_col,requestInt,mpi_err)
              call MPI_WAIT(requestInt,STATUS,mpi_err)
           else
              ! receive lenLoc and update
              dest=MODULO(iprowNew+idx_loc,nprow)
              src=iprowNew+2**(cnt-1)
              if (src<nprow) then
                 src=MODULO(src+idx_loc,nprow)
                 isAlloc=.true.
                 call MPI_IRECV(lenLoc,1,MPI_INT, src,tag,psp_mpi_comm_col,requestInt,mpi_err)!need to use different tags to distinguish messages
                 call MPI_WAIT(requestInt,STATUS,mpi_err)
              else
                 isAlloc=.false.
              end if
           end if
           if (isAlloc) then
              C_loc_MPI%nnz=Hsp%nnz
              C_loc_MPI%loc_dim1=Hsp%loc_dim1
              C_loc_MPI%loc_dim2=Hsp%loc_dim2
              if (allocated(C_loc_idx1)) deallocate(C_loc_idx1)
              if (allocated(C_loc_idx2)) deallocate(C_loc_idx2)
              if (allocated(C_loc_val)) deallocate(C_loc_val)
              allocate(C_loc_idx1(lenLoc))
              allocate(C_loc_idx2(C_loc_MPI%loc_dim2+1))
              allocate(C_loc_val(lenLoc))
              if (MODULO(tempIP,2)==1) then
                 if (Hsp%nnz>0) then
                    do cnt1=1,lenLoc
                       C_loc_idx1(cnt1)=Hsp%row_ind(cnt1)
                       C_loc_val(cnt1)=Hsp%dval(cnt1)
                    end do
                 end if
                 do cnt1=1,C_loc_MPI%loc_dim2+1
                    C_loc_idx2(cnt1)=Hsp%col_ptr(cnt1)
                 end do
              end if

              ! decide buffer size
              call MPI_PACK_SIZE(1,MPI_INTEGER,psp_mpi_comm_col,sizeint,mpi_err)
              call MPI_PACK_SIZE(lenLoc,MPI_INTEGER,psp_mpi_comm_col,sizeidx1,mpi_err)
              call MPI_PACK_SIZE(lenLoc,MPI_DOUBLE_PRECISiON,psp_mpi_comm_col,sizeval,mpi_err)
              call MPI_PACK_SIZE(C_loc_MPI%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,sizeidx2,mpi_err)
              sizebuf=sizeint+sizeidx1+sizeidx2+sizeval
              if (allocated(buf)) deallocate(buf)
              allocate(buf(sizebuf))
           end if
           ! send and receive A_loc
           if (MODULO(tempIP,2)==1) then
              ! pack buffer
              position=0
              call MPI_PACK(C_loc_MPI%nnz,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
              !call MPI_PACK(C_loc_MPI%loc_dim1,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
              !call MPI_PACK(C_loc_MPI%loc_dim2,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
              call MPI_PACK(C_loc_idx1,lenLoc,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
              call MPI_PACK(C_loc_idx2,C_loc_MPI%loc_dim2+1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
              call MPI_PACK(C_loc_val,lenLoc,MPI_DOUBLE_PRECISiON,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)

              ! set up A_loc and send
              dest=MODULO(iprowNew-2**(cnt-1)+idx_loc,nprow)
              call MPI_ISEND(buf,sizebuf,MPI_PACKED,dest,tag,psp_mpi_comm_col,requestSpm,mpi_err)
              call MPI_WAIT(requestSpm,STATUS,mpi_err)
           else
              ! receive A_loc and add to A
              src=iprowNew+2**(cnt-1)
              if (src<nprow) then
                 src=MODULO(src+idx_loc,nprow)
                 call MPI_IRECV(buf,sizebuf,MPI_PACKED, src,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                 call MPI_WAIT(requestSpm,STATUS,mpi_err)

                 ! unpack buffer
                 position=0
                 call MPI_UNPACK(buf,sizebuf,position,C_loc_MPI%nnz,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                 !call MPI_UNPACK(buf,sizebuf,position,C_loc_MPI%loc_dim1,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                 !call MPI_UNPACK(buf,sizebuf,position,C_loc_MPI%loc_dim2,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                 call MPI_UNPACK(buf,sizebuf,position,C_loc_idx1,lenLoc,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                 call MPI_UNPACK(buf,sizebuf,position,C_loc_idx2,C_loc_MPI%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                 call MPI_UNPACK(buf,sizebuf,position,C_loc_val,lenLoc,MPI_DOUBLE_PRECISiON,psp_mpi_comm_col,mpi_err)
                 ! Hsp=C_loc_MPI+Hsp
                 if (C_loc_MPI%nnz>0) then
                    call psp_sst_sum_spmspm(Hsp%loc_dim1,Hsp%loc_dim2,1.0_dp,&
                         C_loc_idx1,C_loc_idx2,C_loc_val,'csc',1.0_dp,&
                         Hsp%row_ind,Hsp%col_ptr,Hsp%dval,'csc',Hsp%row_ind,Hsp%col_ptr,Hsp%dval,&
                         'csc',C_loc_MPI%nnz,Hsp%nnz)
                 end if
                 Hsp%nnz=size(Hsp%dval)
              end if
           end if
           tempIP=tempIP/2
        end if
        !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
     end do
  else
! copy psp_matrix_spm to psp_MPI_spm
     C_loc_MPI%nnz=Hsp%nnz
     C_loc_MPI%loc_dim1=Hsp%loc_dim1
     C_loc_MPI%loc_dim2=Hsp%loc_dim2
     if (allocated(C_loc_MPI%idx1)) deallocate(C_loc_MPI%idx1)
     if (allocated(C_loc_MPI%idx2)) deallocate(C_loc_MPI%idx2)
     if (allocated(C_loc_MPI%val)) deallocate(C_loc_MPI%val)
     allocate(C_loc_MPI%idx1(Hsp%nnz))
     allocate(C_loc_MPI%idx2(C_loc_MPI%loc_dim2+1))
     allocate(C_loc_MPI%val(Hsp%nnz))
     if (Hsp%nnz>0) then
        do cnt1=1,Hsp%nnz
           C_loc_MPI%idx1(cnt1)=Hsp%row_ind(cnt1)
           C_loc_MPI%val(cnt1)=Hsp%dval(cnt1)
        end do
     end if
     do cnt1=1,C_loc_MPI%loc_dim2+1
        C_loc_MPI%idx2(cnt1)=Hsp%col_ptr(cnt1)
     end do

! call psp_MPI_REDUCE
     call psp_MPI_REDUCE_spm_packed(C_loc_MPI, idx_loc, .false.)

! copy psp_MPI_spm to  psp_matrix_spm 
     Hsp%nnz=C_loc_MPI%nnz
     if (allocated(Hsp%row_ind)) deallocate(Hsp%row_ind)
     if (allocated(Hsp%col_ptr)) deallocate(Hsp%col_ptr)
     if (allocated(Hsp%dval)) deallocate(Hsp%dval)
     allocate(Hsp%row_ind(C_loc_MPI%nnz))
     allocate(Hsp%col_ptr(C_loc_MPI%loc_dim2+1))
     allocate(Hsp%dval(C_loc_MPI%nnz))
     if (C_loc_MPI%nnz>0) then
        do cnt1=1,Hsp%nnz
           Hsp%row_ind(cnt1)=C_loc_MPI%idx1(cnt1)
           Hsp%dval(cnt1)=C_loc_MPI%val(cnt1)
        end do
     end if
     do cnt1=1,C_loc_MPI%loc_dim2+1
        Hsp%col_ptr(cnt1)=C_loc_MPI%idx2(cnt1)
     end do

  end if

  ! sparse to dense
  call psp_sp2den_m(Hsp,Hrd,desc_H)
  mm_err=MAXVAL(abs(Hrd-H))
  !err=1.0_dp
  !call MPI_REDUCE(mm_err, err, 1, MPI_DOUBLE_PRECISiON, MPI_SUM, 0, psp_mpi_comm_col, mpi_err)
  if (iprow==idx_loc) print *, 'MPI_reduce_spm error', mm_err

  if (allocated(H)) deallocate(H)
  if (allocated(Hrd)) deallocate(Hrd)

  call psp_deallocate_spm(Hsp)
  call psp_deallocate_spm(Hrdsp)

  call mpi_finalize(mpi_err)



end program testMPI


