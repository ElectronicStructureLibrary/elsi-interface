MODULE pspMPI
  use pspUtility

#ifdef MPI
  use mpi
#endif

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**********************************************!

  !**** TYPES *************************************!

  type psp_MPI_dspm
     integer :: nnz      ! local number of nonzero entries
     integer :: loc_dim1 ! (local) row dimension size of the local matrix
     integer :: loc_dim2 ! (local) column dimension size of the local matrix
     integer, allocatable  :: idx1(:)
     integer, allocatable  :: idx2(:)
     real(dp), allocatable  :: val(:)  ! matrix elements for a real matrix
  end type psp_MPI_dspm

  type psp_MPI_zspm
     integer :: nnz      ! local number of nonzero entries
     integer :: loc_dim1 ! (local) row dimension size of the local matrix
     integer :: loc_dim2 ! (local) column dimension size of the local matrix
     integer, allocatable  :: idx1(:)
     integer, allocatable  :: idx2(:)
     complex(dp), allocatable  :: val(:)
  end type psp_MPI_zspm

  !**** INTERFACES ********************************!

  interface psp_gridinit_2D
     module procedure psp_gridinit_2D
  end interface psp_gridinit_2D

  interface psp_gridinit_3D
     module procedure psp_gridinit_3D
  end interface psp_gridinit_3D

  interface psp_MPI_REDUCE_dspm_struct ! TODO: bug inside
     module procedure psp_MPI_REDUCE_dspm_struct
  end interface psp_MPI_REDUCE_dspm_struct

  interface psp_MPI_REDUCE_zspm_struct ! TODO: bug inside
     module procedure psp_MPI_REDUCE_zspm_struct
  end interface psp_MPI_REDUCE_zspm_struct

  interface psp_MPI_REDUCE_spm_packed
     module procedure psp_MPI_REDUCE_dspm_packed
     module procedure psp_MPI_REDUCE_zspm_packed
  end interface psp_MPI_REDUCE_spm_packed

  interface psp_MPI_REDUCE_den
     module procedure psp_MPI_REDUCE_dden
     module procedure psp_MPI_REDUCE_zden
  end interface psp_MPI_REDUCE_den

  public :: psp_gridinit_2D
  public :: psp_gridinit_3D
  public :: psp_MPI_REDUCE_dspm_struct
  public :: psp_MPI_REDUCE_zspm_struct
  public :: psp_MPI_REDUCE_spm_packed
  public :: psp_MPI_dspm
  public :: psp_MPI_zspm
  public :: psp_MPI_REDUCE_den

contains

  !================================================!
  !       MPI_REDUCE for sparse matrices           !
  !================================================!
  subroutine psp_MPI_REDUCE_dspm_packed(A, idx_proc, isRow)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: idx_proc
    logical, intent(in) :: isRow
    ! isRow=.true., then reduce in processor rows
    ! otherwise, reduce in processor columns

    !**** INOUT ***********************************!

    type(psp_MPI_dspm), intent(inout ) :: A

    !**** INTERNAL ********************************!

    integer :: tempIP, numStep, cnt, cnt1, dest, src, requestInt, requestSpm, iprowNew, ipcolNew
    integer :: lenLoc, mpi_err, STATUS(MPI_STATUS_SIZE), nprow, npcol, iprow, ipcol, tag=0
    type(psp_MPI_dspm) :: A_loc
    real(dp) :: temp
    integer :: sizeidx1,sizeidx2,sizeval,sizeint, sizebuf, position
    integer, allocatable :: buf(:)
    logical :: isAlloc

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
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    if (isRow) then
       ipcolNew = MODULO(ipcol-idx_proc,npcol)
       tempIP=ipcolNew
       numStep=0
       temp=real(psp_npcol)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(ipcolNew,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=MODULO(ipcolNew-2**(cnt-1)+idx_proc,psp_npcol)
                src=MODULO(ipcolNew+idx_proc,psp_npcol)
                isAlloc=.true.
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_row,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=MODULO(ipcolNew+idx_proc,psp_npcol)
                src=ipcolNew+2**(cnt-1)
                if (src<psp_npcol) then
                   src=MODULO(src+idx_proc,psp_npcol)
                   isAlloc=.true.
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_row,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                else
                   isAlloc=.false.
                end if
             end if

             if (isAlloc) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%idx1(cnt1)
                         A_loc%val(cnt1)=A%val(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%idx2(cnt1)
                   end do
                end if

                ! decide buffer size
                call MPI_PACK_SIZE(1,MPI_INTEGER,psp_mpi_comm_row,sizeint,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_INTEGER,psp_mpi_comm_row,sizeidx1,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_DOUBLE_PRECISION,psp_mpi_comm_row,sizeval,mpi_err)
                call MPI_PACK_SIZE(A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_row,sizeidx2,mpi_err)
                sizebuf=sizeint+sizeidx1+sizeidx2+sizeval
                if (allocated(buf)) deallocate(buf)
                allocate(buf(sizebuf))
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! pack buffer
                position=0
                call MPI_PACK(A_loc%nnz,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim1,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim2,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx1,lenLoc,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%val,lenLoc,MPI_DOUBLE_PRECISION,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)

                ! set up A_loc and send
                dest=MODULO(ipcolNew-2**(cnt-1)+idx_proc,psp_npcol)
                call MPI_ISEND(buf,sizebuf,MPI_PACKED,dest,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=ipcolNew+2**(cnt-1)
                if (src<psp_npcol) then
                   src=MODULO(src+idx_proc,psp_npcol)
                   call MPI_IRECV(buf,sizebuf,MPI_PACKED, src,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)

                   ! unpack buffer
                   position=0
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%nnz,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim1,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim2,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx1,lenLoc,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%val,lenLoc,MPI_DOUBLE_PRECISION,psp_mpi_comm_col,mpi_err)

                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,1.0_dp,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',1.0_dp,&
                           A%idx1,A%idx2,A%val,'csc',A%idx1,A%idx2,A%val,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%val)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    else
       iprowNew=MODULO(iprow-idx_proc,psp_nprow)
       tempIP=iprowNew
       numStep=0
       temp=real(psp_nprow)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(iprowNew,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=MODULO(iprowNew-2**(cnt-1)+idx_proc,psp_nprow)
                src=MODULO(iprowNew+idx_proc,psp_nprow)
                isAlloc=.true.
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_col,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=MODULO(iprowNew+idx_proc,psp_nprow)
                src=iprowNew+2**(cnt-1)
                if (src<psp_nprow) then
                   src=MODULO(src+idx_proc,psp_nprow)
                   isAlloc=.true.
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_col,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                else
                   isAlloc=.false.
                end if
             end if

             if (isAlloc) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%idx1(cnt1)
                         A_loc%val(cnt1)=A%val(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%idx2(cnt1)
                   end do
                end if

                ! decide buffer size
                call MPI_PACK_SIZE(1,MPI_INTEGER,psp_mpi_comm_col,sizeint,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_INTEGER,psp_mpi_comm_col,sizeidx1,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_DOUBLE_PRECISION,psp_mpi_comm_col,sizeval,mpi_err)
                call MPI_PACK_SIZE(A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,sizeidx2,mpi_err)
                sizebuf=sizeint+sizeidx1+sizeidx2+sizeval
                if (allocated(buf)) deallocate(buf)
                allocate(buf(sizebuf))
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! pack buffer
                position=0
                call MPI_PACK(A_loc%nnz,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim1,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim2,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx1,lenLoc,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%val,lenLoc,MPI_DOUBLE_PRECISION,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)

                ! set up A_loc and send
                dest=MODULO(iprowNew-2**(cnt-1)+idx_proc,psp_nprow)
                call MPI_ISEND(buf,sizebuf,MPI_PACKED,dest,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=iprowNew+2**(cnt-1)
                if (src<psp_nprow) then
                   src=MODULO(src+idx_proc,psp_nprow)
                   call MPI_IRECV(buf,sizebuf,MPI_PACKED, src,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)

                   ! unpack buffer
                   position=0
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%nnz,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim1,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim2,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx1,lenLoc,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%val,lenLoc,MPI_DOUBLE_PRECISION,psp_mpi_comm_col,mpi_err)

                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,1.0_dp,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',1.0_dp,&
                           A%idx1,A%idx2,A%val,'csc',A%idx1,A%idx2,A%val,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%val)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    end if

    if (allocated(A_loc%val)) deallocate(A_loc%val)
    if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
    if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)

  end subroutine psp_MPI_REDUCE_dspm_packed

  subroutine psp_MPI_REDUCE_zspm_packed(A, idx_proc, isRow)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: idx_proc
    logical, intent(in) :: isRow
    ! isRow=.true., then reduce in processor rows
    ! otherwise, reduce in processor columns

    !**** INOUT ***********************************!

    type(psp_MPI_zspm), intent(inout ) :: A

    !**** INTERNAL ********************************!

    integer :: tempIP, numStep, cnt, cnt1, dest, src, requestInt, requestSpm, iprowNew, ipcolNew
    integer :: lenLoc, mpi_err, STATUS(MPI_STATUS_SIZE), nprow, npcol, iprow, ipcol, tag=0
    type(psp_MPI_zspm) :: A_loc
    real(dp) :: temp
    integer :: sizeidx1,sizeidx2,sizeval,sizeint, sizebuf, position
    integer, allocatable :: buf(:)
    logical :: isAlloc

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
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    if (isRow) then
       ipcolNew = MODULO(ipcol-idx_proc,npcol)
       tempIP=ipcolNew
       numStep=0
       temp=real(psp_npcol)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(ipcolNew,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=MODULO(ipcolNew-2**(cnt-1)+idx_proc,psp_npcol)
                src=MODULO(ipcolNew+idx_proc,psp_npcol)
                isAlloc=.true.
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_row,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=MODULO(ipcolNew+idx_proc,psp_npcol)
                src=ipcolNew+2**(cnt-1)
                if (src<psp_npcol) then
                   src=MODULO(src+idx_proc,psp_npcol)
                   isAlloc=.true.
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_row,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                else
                   isAlloc=.false.
                end if
             end if

             if (isAlloc) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%idx1(cnt1)
                         A_loc%val(cnt1)=A%val(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%idx2(cnt1)
                   end do
                end if

                ! decide buffer size
                call MPI_PACK_SIZE(1,MPI_INTEGER,psp_mpi_comm_row,sizeint,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_INTEGER,psp_mpi_comm_row,sizeidx1,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_DOUBLE_COMPLEX,psp_mpi_comm_row,sizeval,mpi_err)
                call MPI_PACK_SIZE(A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_row,sizeidx2,mpi_err)
                sizebuf=sizeint+sizeidx1+sizeidx2+sizeval
                if (allocated(buf)) deallocate(buf)
                allocate(buf(sizebuf))
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! pack buffer
                position=0
                call MPI_PACK(A_loc%nnz,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim1,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim2,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx1,lenLoc,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%val,lenLoc,MPI_DOUBLE_COMPLEX,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)

                ! set up A_loc and send
                dest=MODULO(ipcolNew-2**(cnt-1)+idx_proc,psp_npcol)
                call MPI_ISEND(buf,sizebuf,MPI_PACKED,dest,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=ipcolNew+2**(cnt-1)
                if (src<psp_npcol) then
                   src=MODULO(src+idx_proc,psp_npcol)
                   call MPI_IRECV(buf,sizebuf,MPI_PACKED, src,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)

                   ! unpack buffer
                   position=0
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%nnz,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim1,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim2,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx1,lenLoc,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%val,lenLoc,MPI_DOUBLE_COMPLEX,psp_mpi_comm_col,mpi_err)

                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,cmplx_1,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',cmplx_1,&
                           A%idx1,A%idx2,A%val,'csc',A%idx1,A%idx2,A%val,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%val)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    else
       iprowNew=MODULO(iprow-idx_proc,psp_nprow)
       tempIP=iprowNew
       numStep=0
       temp=real(psp_nprow)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(iprowNew,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=MODULO(iprowNew-2**(cnt-1)+idx_proc,psp_nprow)
                src=MODULO(iprowNew+idx_proc,psp_nprow)
                isAlloc=.true.
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_col,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=MODULO(iprowNew+idx_proc,psp_nprow)
                src=iprowNew+2**(cnt-1)
                if (src<psp_nprow) then
                   src=MODULO(src+idx_proc,psp_nprow)
                   isAlloc=.true.
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_col,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                else
                   isAlloc=.false.
                end if
             end if

             if (isAlloc) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%idx1(cnt1)
                         A_loc%val(cnt1)=A%val(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%idx2(cnt1)
                   end do
                end if

                ! decide buffer size
                call MPI_PACK_SIZE(1,MPI_INTEGER,psp_mpi_comm_col,sizeint,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_INTEGER,psp_mpi_comm_col,sizeidx1,mpi_err)
                call MPI_PACK_SIZE(lenLoc,MPI_DOUBLE_COMPLEX,psp_mpi_comm_col,sizeval,mpi_err)
                call MPI_PACK_SIZE(A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,sizeidx2,mpi_err)
                sizebuf=sizeint+sizeidx1+sizeidx2+sizeval
                if (allocated(buf)) deallocate(buf)
                allocate(buf(sizebuf))
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! pack buffer
                position=0
                call MPI_PACK(A_loc%nnz,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim1,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                !call MPI_PACK(A_loc%loc_dim2,1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx1,lenLoc,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)
                call MPI_PACK(A_loc%val,lenLoc,MPI_DOUBLE_COMPLEX,buf,sizebuf,position,psp_mpi_comm_col,mpi_err)

                ! set up A_loc and send
                dest=MODULO(iprowNew-2**(cnt-1)+idx_proc,psp_nprow)
                call MPI_ISEND(buf,sizebuf,MPI_PACKED,dest,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=iprowNew+2**(cnt-1)
                if (src<psp_nprow) then
                   src=MODULO(src+idx_proc,psp_nprow)
                   call MPI_IRECV(buf,sizebuf,MPI_PACKED, src,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)

                   ! unpack buffer
                   position=0
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%nnz,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim1,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   !call MPI_UNPACK(buf,sizebuf,position,A_loc%loc_dim2,1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx1,lenLoc,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%idx2,A_loc%loc_dim2+1,MPI_INTEGER,psp_mpi_comm_col,mpi_err)
                   call MPI_UNPACK(buf,sizebuf,position,A_loc%val,lenLoc,MPI_DOUBLE_COMPLEX,psp_mpi_comm_col,mpi_err)

                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,cmplx_1,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',cmplx_1,&
                           A%idx1,A%idx2,A%val,'csc',A%idx1,A%idx2,A%val,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%val)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    end if

    if (allocated(A_loc%val)) deallocate(A_loc%val)
    if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
    if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)

  end subroutine psp_MPI_REDUCE_zspm_packed

  !================================================!
  !       MPI_REDUCE for sparse matrices           !
  !================================================!
  subroutine psp_MPI_REDUCE_dspm_struct(A, idx_proc, isRow)
    implicit none
    ! TODO: there is a bug in MPI_Get_address

    !**** INPUT ***********************************!

    integer, intent(in) :: idx_proc
    logical, intent(in) :: isRow
    ! isRow=.true., then reduce in processor rows
    ! otherwise, reduce in processor columns

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout ) :: A

    !**** INTERNAL ********************************!

    integer :: tempIP, numStep, cnt, cnt1, dest, src, requestInt, requestSpm
    integer :: lenLoc, mpi_err, STATUS(MPI_STATUS_SIZE), nprow, npcol, iprow, ipcol, tag=0
    integer, dimension(6) :: disp, blockLen
    integer, dimension(6) :: tp=(/ MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_DOUBLE /)
    type(psp_MPI_dspm) :: A_loc
    real(dp) :: temp
    integer :: psp_MPI_datatype_spm

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
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)

    if (isRow) then
       tempIP=ipcol
       numStep=0
       temp=real(psp_npcol)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(ipcol,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=ipcol-2**(cnt-1)
                src=ipcol
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_row,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=ipcol
                src=ipcol+2**(cnt-1)
                if (src<psp_npcol) then
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_row,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                end if
             end if

             if (src<psp_npcol) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%row_ind(cnt1)
                         A_loc%val(cnt1)=A%dval(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%col_ptr(cnt1)
                   end do
                end if

                ! define new MPI datatype
                call MPI_ADDRESS(A_loc%nnz, disp(1), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim1, disp(2), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim2, disp(3), mpi_err)
                call MPI_ADDRESS(A_loc%idx1, disp(4), mpi_err)
                call MPI_ADDRESS(A_loc%idx2, disp(5), mpi_err)
                call MPI_ADDRESS(A_loc%val, disp(6), mpi_err)
                disp(2)=disp(2)-disp(1)
                disp(3)=disp(3)-disp(1)
                disp(4)=disp(4)-disp(1)
                disp(5)=disp(5)-disp(1)
                disp(6)=disp(6)-disp(1)
                disp(1)=0
                blockLen=1
                blockLen(4)=lenLoc
                blockLen(6)=lenLoc
                blockLen(5)=A_loc%loc_dim2+1
                call MPI_Type_struct(6, blockLen, disp, tp, psp_MPI_datatype_spm, mpi_err)
                call MPI_Type_commit(psp_MPI_datatype_spm, mpi_err)
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! set up A_loc and send
                dest=ipcol-2**(cnt-1)
                call MPI_ISEND(A_loc,1,psp_MPI_datatype_spm,dest,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=ipcol+2**(cnt-1)
                if (src<psp_npcol) then
                   call MPI_IRECV(A_loc,1,psp_MPI_datatype_spm, src,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)
                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,1.0_dp,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',1.0_dp,&
                           A%row_ind,A%col_ptr,A%dval,'csc',A%row_ind,A%col_ptr,A%dval,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%dval)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    else
       tempIP=iprow
       numStep=0
       temp=real(psp_nprow)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(iprow,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=iprow-2**(cnt-1)
                src=iprow
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_col,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=iprow
                src=iprow+2**(cnt-1)
                if (src<psp_nprow) then
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_col,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                end if
             end if

             if (src<psp_nprow) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%row_ind(cnt1)
                         A_loc%val(cnt1)=A%dval(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%col_ptr(cnt1)
                   end do
                end if

                ! define new MPI datatype
                call MPI_ADDRESS(A_loc%nnz, disp(1), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim1, disp(2), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim2, disp(3), mpi_err)
                call MPI_ADDRESS(A_loc%idx1, disp(4), mpi_err)
                call MPI_ADDRESS(A_loc%idx2, disp(5), mpi_err)
                call MPI_ADDRESS(A_loc%val, disp(6), mpi_err)
                disp(2)=disp(2)-disp(1)
                disp(3)=disp(3)-disp(1)
                disp(4)=disp(4)-disp(1)
                disp(5)=disp(5)-disp(1)
                disp(6)=disp(6)-disp(1)
                disp(1)=0
                blockLen=1
                blockLen(4)=lenLoc
                blockLen(6)=lenLoc
                blockLen(5)=A_loc%loc_dim2+1
                call MPI_Type_struct(6, blockLen, disp, tp, psp_MPI_datatype_spm, mpi_err)
                call MPI_Type_commit(psp_MPI_datatype_spm, mpi_err)
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! set up A_loc and send
                src=iprow
                dest=iprow-2**(cnt-1)
                call MPI_ISEND(A_loc,1,psp_MPI_datatype_spm,dest,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=iprow+2**(cnt-1)
                if (src<psp_nprow) then
                   call MPI_IRECV(A_loc,1,psp_MPI_datatype_spm, src,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)
                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,1.0_dp,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',1.0_dp,&
                           A%row_ind,A%col_ptr,A%dval,'csc',A%row_ind,A%col_ptr,A%dval,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%dval)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    end if

    if (allocated(A_loc%val)) deallocate(A_loc%val)
    if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
    if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)

  end subroutine psp_MPI_REDUCE_dspm_struct

  subroutine psp_MPI_REDUCE_zspm_struct(A, idx_proc, isRow)
    implicit none
    ! TODO: there is a bug in MPI_Get_address

    !**** INPUT ***********************************!

    integer, intent(in) :: idx_proc
    logical, intent(in) :: isRow
    ! isRow=.true., then reduce in processor rows
    ! otherwise, reduce in processor columns

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(inout ) :: A

    !**** INTERNAL ********************************!

    integer :: tempIP, numStep, cnt, cnt1, dest, src, requestInt, requestSpm
    integer :: lenLoc, mpi_err, STATUS(MPI_STATUS_SIZE), nprow, npcol, iprow, ipcol, tag=0
    integer, dimension(6) :: disp, blockLen
    integer, dimension(6) :: tp=(/ MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,&
         MPI_INTEGER,MPI_DOUBLE_COMPLEX /)
    type(psp_MPI_zspm) :: A_loc
    real(dp) :: temp
    integer :: psp_MPI_datatype_spm

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
    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    if (isRow) then
       tempIP=ipcol
       numStep=0
       temp=real(psp_npcol)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(ipcol,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=ipcol-2**(cnt-1)
                src=ipcol
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_row,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=ipcol
                src=ipcol+2**(cnt-1)
                if (src<psp_npcol) then
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_row,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                end if
             end if

             if (src<psp_npcol) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%row_ind(cnt1)
                         A_loc%val(cnt1)=A%zval(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%col_ptr(cnt1)
                   end do
                end if

                ! define new MPI datatype
                call MPI_ADDRESS(A_loc%nnz, disp(1), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim1, disp(2), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim2, disp(3), mpi_err)
                call MPI_ADDRESS(A_loc%idx1, disp(4), mpi_err)
                call MPI_ADDRESS(A_loc%idx2, disp(5), mpi_err)
                call MPI_ADDRESS(A_loc%val, disp(6), mpi_err)
                disp(2)=disp(2)-disp(1)
                disp(3)=disp(3)-disp(1)
                disp(4)=disp(4)-disp(1)
                disp(5)=disp(5)-disp(1)
                disp(6)=disp(6)-disp(1)
                disp(1)=0
                blockLen=1
                blockLen(4)=lenLoc
                blockLen(6)=lenLoc
                blockLen(5)=A_loc%loc_dim2+1
                call MPI_Type_struct(6, blockLen, disp, tp, psp_MPI_datatype_spm, mpi_err)
                call MPI_Type_commit(psp_MPI_datatype_spm, mpi_err)
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! set up A_loc and send
                dest=ipcol-2**(cnt-1)
                call MPI_ISEND(A_loc,1,psp_MPI_datatype_spm,dest,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=ipcol+2**(cnt-1)
                if (src<psp_npcol) then
                   call MPI_IRECV(A_loc,1,psp_MPI_datatype_spm, src,tag,psp_mpi_comm_row,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)
                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,cmplx_1,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',cmplx_1,&
                           A%row_ind,A%col_ptr,A%zval,'csc',A%row_ind,A%col_ptr,A%zval,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%zval)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    else
       tempIP=iprow
       numStep=0
       temp=real(psp_nprow)
       do while (temp>1)
          numStep=numStep+1
          temp=temp/2
       end do
       do cnt=1,numStep
          if (MODULO(iprow,2**(cnt-1))==0) then
             ! get uniform data size
             if (MODULO(tempIP,2)==1) then
                ! set up lenLoc and send
                lenLoc=A%nnz
                dest=iprow-2**(cnt-1)
                src=iprow
                call MPI_ISEND(lenLoc,1,MPI_INTEGER,dest,tag,psp_mpi_comm_col,requestInt,mpi_err)
                call MPI_WAIT(requestInt,STATUS,mpi_err)
             else
                ! receive lenLoc and update
                dest=iprow
                src=iprow+2**(cnt-1)
                if (src<psp_nprow) then
                   call MPI_IRECV(lenLoc,1,MPI_INTEGER, src,tag,psp_mpi_comm_col,requestInt,mpi_err)
                   call MPI_WAIT(requestInt,STATUS,mpi_err)
                end if
             end if

             if (src<psp_nprow) then
                A_loc%nnz=A%nnz
                A_loc%loc_dim1=A%loc_dim1
                A_loc%loc_dim2=A%loc_dim2
                if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
                if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)
                if (allocated(A_loc%val)) deallocate(A_loc%val)
                allocate(A_loc%idx1(lenLoc))
                allocate(A_loc%idx2(A_loc%loc_dim2+1))
                allocate(A_loc%val(lenLoc))
                if (MODULO(tempIP,2)==1) then
                   if (A%nnz>0) then
                      do cnt1=1,lenLoc
                         A_loc%idx1(cnt1)=A%row_ind(cnt1)
                         A_loc%val(cnt1)=A%zval(cnt1)
                      end do
                   end if
                   do cnt1=1,A_loc%loc_dim2+1
                      A_loc%idx2(cnt1)=A%col_ptr(cnt1)
                   end do
                end if

                ! define new MPI datatype
                call MPI_ADDRESS(A_loc%nnz, disp(1), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim1, disp(2), mpi_err)
                call MPI_ADDRESS(A_loc%loc_dim2, disp(3), mpi_err)
                call MPI_ADDRESS(A_loc%idx1, disp(4), mpi_err)
                call MPI_ADDRESS(A_loc%idx2, disp(5), mpi_err)
                call MPI_ADDRESS(A_loc%val, disp(6), mpi_err)
                disp(2)=disp(2)-disp(1)
                disp(3)=disp(3)-disp(1)
                disp(4)=disp(4)-disp(1)
                disp(5)=disp(5)-disp(1)
                disp(6)=disp(6)-disp(1)
                disp(1)=0
                blockLen=1
                blockLen(4)=lenLoc
                blockLen(6)=lenLoc
                blockLen(5)=A_loc%loc_dim2+1
                call MPI_Type_struct(6, blockLen, disp, tp, psp_MPI_datatype_spm, mpi_err)
                call MPI_Type_commit(psp_MPI_datatype_spm, mpi_err)
             end if

             ! send and receive A_loc
             if (MODULO(tempIP,2)==1) then
                ! set up A_loc and send
                dest=iprow-2**(cnt-1)
                call MPI_ISEND(A_loc,1,psp_MPI_datatype_spm,dest,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and add to A
                src=iprow+2**(cnt-1)
                if (src<psp_nprow) then
                   call MPI_IRECV(A_loc,1,psp_MPI_datatype_spm, src,tag,psp_mpi_comm_col,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)
                   ! A=A_loc+A
                   if (A_loc%nnz>0) then
                      call psp_sst_sum_spmspm(A%loc_dim1,A%loc_dim2,cmplx_1,&
                           A_loc%idx1,A_loc%idx2,A_loc%val,'csc',cmplx_1,&
                           A%row_ind,A%col_ptr,A%zval,'csc',A%row_ind,A%col_ptr,A%zval,&
                           'csc',A_loc%nnz,A%nnz)
                   end if
                   A%nnz=size(A%zval)
                end if
             end if
             tempIP=tempIP/2
          end if
          !call MPI_Barrier(psp_mpi_comm_world,mpi_err)
       end do
    end if

    if (allocated(A_loc%val)) deallocate(A_loc%val)
    if (allocated(A_loc%idx1)) deallocate(A_loc%idx1)
    if (allocated(A_loc%idx2)) deallocate(A_loc%idx2)

  end subroutine psp_MPI_REDUCE_zspm_struct


  !================================================!
  ! initialize 2D grid topology for psp            !
  !================================================!
  subroutine psp_gridinit_2D(mpi_comm_in,mpi_size,nprow,order,bs_def_row,bs_def_col,icontxt)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: order ! ordering of processor grid: 'r/R' or other for row-major, 'c/C' for column-major
    integer, intent(in) :: mpi_comm_in
    integer, intent(in) :: mpi_size ! total number of MPI processes for the processor grid
    integer, intent(in) :: nprow ! number of rows in the processor grid
    integer, intent(in) :: bs_def_row ! default block size in row
    integer, intent(in) :: bs_def_col ! default block size in column
    integer, intent(in), optional :: icontxt ! existing BLACS context handle in case ScaLAPACK is already initialized

    !**** INTERNAL ********************************!

    integer :: i, mpi_err
    logical :: remain_dims(2), periods(2), reorder
    integer :: mpi_rank, dims(2), coords(2)

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

    psp_mpi_comm_world=mpi_comm_in
    psp_mpi_size=mpi_size
    psp_nprow=nprow
    psp_npcol=mpi_size/nprow
    psp_proc_order=order
    psp_bs_def_row=bs_def_row
    psp_bs_def_col=bs_def_col
    if (bs_def_row.eq.bs_def_col) then
       psp_update_rank = bs_def_row
    else
       psp_update_rank=1
    end if

    if (present(icontxt)) then
       psp_icontxt=icontxt
    else
       psp_icontxt=psp_mpi_comm_world
       call blacs_gridinit(psp_icontxt,psp_proc_order,psp_nprow,psp_npcol)
    end if

    ! generate sub-topology for row and column boardcasting
    dims(1)=psp_nprow
    dims(2)=psp_npcol
    periods(1)=.true.
    periods(2)=.true.
    reorder=.false.
    call MPI_Cart_create(psp_mpi_comm_world, 2, dims, periods, reorder, psp_mpi_comm_cart,mpi_err)
    ! generate sub-topology for row and column boardcasting
    remain_dims(1) = .false.
    remain_dims(2) = .true.
    call MPI_Cart_Sub(psp_mpi_comm_cart,remain_dims,psp_mpi_comm_row,mpi_err)
    remain_dims(1) = .true.
    remain_dims(2) = .false.
    call MPI_Cart_Sub(psp_mpi_comm_cart,remain_dims,psp_mpi_comm_col,mpi_err)

  end subroutine psp_gridinit_2D

  !================================================!
  ! initialize grid topology for psp               !
  !================================================!
  subroutine psp_gridinit_3D(mpi_comm_in,mpi_size,nprow,npcol,nplay,order,bs_def_row,bs_def_col,icontxt)
    implicit none

    !**** INPUT ***********************************!

    character(1), intent(in) :: order ! ordering of processor grid: 'r/R' or other for row-major, 'c/C' for column-major
    integer, intent(in) :: mpi_comm_in
    integer, intent(in) :: mpi_size ! total number of MPI processes for the processor grid
    integer, intent(in) :: nprow ! number of rows in the processor grid
    integer, intent(in) :: npcol ! number of columns in the processor grid
    integer, intent(in) :: nplay ! number of procersor layers
    integer, intent(in) :: bs_def_row ! default block size in row
    integer, intent(in) :: bs_def_col ! default block size in column
    integer, intent(in), optional :: icontxt ! existing BLACS context handle in case ScaLAPACK is already initialized

    !**** INTERNAL ********************************!

    integer :: i, mpi_err
    logical :: remain_dims(3), periods(3), reorder
    integer :: mpi_rank, dims(3), coords(3)

    !**** GLOBAL **********************************!
#ifdef MPI
    character(1) :: psp_proc_order

    integer :: psp_mpi_comm_world
    integer :: psp_mpi_size
    integer :: psp_nprow
    integer :: psp_npcol
    integer :: psp_nplay
    integer :: psp_bs_def_row
    integer :: psp_bs_def_col
    integer :: psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    integer :: psp_bs_num
    integer :: psp_icontxt ! BLACS context handle used by psp
    integer :: psp_mpi_comm_cart
    integer :: psp_mpi_comm_row, psp_mpi_comm_col, psp_mpi_comm_lay


    common /psp_grid2D/ psp_mpi_comm_world, psp_mpi_size, psp_nprow, psp_npcol
    common /psp_grid3D/ psp_nplay
    common /psp_grid2D/ psp_bs_def_row, psp_bs_def_col
    common /psp_grid2D/ psp_update_rank ! In SUMMA for C=A*B, C is computed by rank psp_update_rank=1 local update
    common /psp_grid2D/ psp_bs_num
    common /psp_grid2D/ psp_icontxt ! BLACS context handle used by psp
    common /psp_grid2D/ psp_mpi_comm_cart, psp_mpi_comm_row,  psp_mpi_comm_col
    common /psp_grid3D/ psp_mpi_comm_lay
#endif

    !**********************************************!

    psp_mpi_comm_world=mpi_comm_in
    psp_mpi_size=mpi_size
    psp_nprow=nprow
    psp_npcol=npcol
    psp_nplay=nplay
    psp_proc_order=order
    psp_bs_def_row=bs_def_row
    psp_bs_def_col=bs_def_col
    if (bs_def_row.eq.bs_def_col) then
       psp_update_rank = bs_def_row
    else
       psp_update_rank=1
    end if

    if (present(icontxt)) then
       psp_icontxt=icontxt
    else
       psp_icontxt=psp_mpi_comm_world
    end if

    ! generate sub-topology for row and column boardcasting
    dims(1)=psp_nprow
    dims(2)=psp_npcol
    dims(3)=psp_nplay
    periods(1)=.true.
    periods(2)=.true.
    periods(3)=.true.
    reorder=.false.
    call MPI_Cart_create(psp_mpi_comm_world, 3, dims, periods, reorder, psp_mpi_comm_cart,mpi_err)
    ! generate sub-topology for row and column boardcasting in each layer
    remain_dims(1) = .false.
    remain_dims(2) = .true.
    remain_dims(3) = .false.
    call MPI_Cart_Sub(psp_mpi_comm_cart,remain_dims,psp_mpi_comm_row,mpi_err)
    remain_dims(1) = .true.
    remain_dims(2) = .false.
    remain_dims(3) = .false.
    call MPI_Cart_Sub(psp_mpi_comm_cart,remain_dims,psp_mpi_comm_col,mpi_err)
    remain_dims(1) = .false.
    remain_dims(2) = .false.
    remain_dims(3) = .true.
    call MPI_Cart_Sub(psp_mpi_comm_cart,remain_dims,psp_mpi_comm_lay,mpi_err)

  end subroutine psp_gridinit_3D

  subroutine psp_MPI_REDUCE_dden(H,m,n)

    !**** INOUT ***********************************!
    integer, intent(in) :: m, n
    real(dp), intent(inout ) :: H(:,:)

    !**** INTERNAL ********************************!

    integer :: tempIP, numStep, cnt, cnt1, dest, src, requestInt, requestSpm
    integer :: lenLoc, mpi_err, STATUS(MPI_STATUS_SIZE), nprow, npcol, iprow, ipcol
    integer, dimension(2) :: H_dim
    integer, dimension(6) :: disp, blockLen
    integer, dimension(6) :: tp=(/ MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_DOUBLE_COMPLEX /)
    real(dp), allocatable :: Hrd(:,:)
    real(dp) :: temp
    integer :: psp_MPI_datatype_spm

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

    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    H_dim=shape(H)
    allocate(Hrd(H_dim(1),H_dim(2)))
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
                H_dim(1)=numroc(m,psp_bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
                H_dim(2)=numroc(n,psp_bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
                call MPI_ISEND(Hrd,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISION,dest,0,psp_mpi_comm_col,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and and to A
                src=iprow+2**(cnt-1)
                dest=iprow
                if (src<nprow) then
                   H_dim(1)=numroc(m,psp_bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
                   H_dim(2)=numroc(n,psp_bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
                   call MPI_IRECV(Hrd,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISION, src,0,psp_mpi_comm_col,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)
                   H=H+Hrd
                end if
             end if
             tempIP=tempIP/2
          end if
       end do
    end if

    deallocate(Hrd)

  end subroutine psp_MPI_REDUCE_dden

  subroutine psp_MPI_REDUCE_zden(H,m,n)

    !**** INOUT ***********************************!
    integer, intent(in) :: m, n
    complex(dp), intent(inout ) :: H(:,:)

    !**** INTERNAL ********************************!

    integer :: tempIP, numStep, cnt, cnt1, dest, src, requestInt, requestSpm
    integer :: lenLoc, mpi_err, STATUS(MPI_STATUS_SIZE), nprow, npcol, iprow, ipcol
    integer, dimension(2) :: H_dim
    integer, dimension(6) :: disp, blockLen
    integer, dimension(6) :: tp=(/ MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_DOUBLE_COMPLEX /)
    complex(dp), allocatable :: Hrd(:,:)
    real(dp) :: temp
    integer :: psp_MPI_datatype_spm

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

    call blacs_gridinfo(psp_icontxt,nprow,npcol,iprow,ipcol)
    H_dim=shape(H)
    allocate(Hrd(H_dim(1),H_dim(2)))
    Hrd=cmplx_0
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
                H_dim(1)=numroc(m,psp_bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
                H_dim(2)=numroc(n,psp_bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
                call MPI_ISEND(Hrd,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISION,dest,0,psp_mpi_comm_col,requestSpm,mpi_err)
                call MPI_WAIT(requestSpm,STATUS,mpi_err)
             else
                ! receive A_loc and and to A
                src=iprow+2**(cnt-1)
                dest=iprow
                if (src<nprow) then
                   H_dim(1)=numroc(m,psp_bs_def_row,iprow,0,nprow) ! use 'numroc' to compute the size of a local matrix, row
                   H_dim(2)=numroc(n,psp_bs_def_col,ipcol,0,npcol) ! use 'numroc' to compute the size of a local matrix, column
                   call MPI_IRECV(Hrd,H_dim(1)*H_dim(2),MPI_DOUBLE_PRECISION, src,0,psp_mpi_comm_col,requestSpm,mpi_err)
                   call MPI_WAIT(requestSpm,STATUS,mpi_err)
                   H=H+Hrd
                end if
             end if
             tempIP=tempIP/2
          end if
       end do
    end if

    deallocate(Hrd)

  end subroutine psp_MPI_REDUCE_zden


END MODULE pspMPI
