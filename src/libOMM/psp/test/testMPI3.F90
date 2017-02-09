!
! Template for a Monte Carlo code.
!
! The root processor comes up with a list of seeds
!    which it ships to all processors.
!
! In a real application, each processor would compute
!    something and send it back.  Here they compute
!    a vector of random numbers and send it back.
!
! This version uses a single BCAST to distribute
!   both integer and double precision data.  It packs
!   integer vectors as well as double precision data
!   into an MPI_PACKED buffer.
!   This illustrates the use of MPI_PACK and MPI_UNPACK commands.
!
program monte3
use pspBLAS
  !
  include 'mpif.h'
  !
  integer my_rank
  integer p
  integer source
  integer dest
  integer tag
  integer iseed, initseed, initvec(200)
  integer status(MPI_STATUS_SIZE)
  integer ierr
  integer i  , j , root, sizeidx1,sizeidx2,sizeval,sizeint, sizebuf
  integer isize, position, kquad, nj
  integer itemp(100)
integer, allocatable :: BUF(:)
  real*8 ans(10), ans2(10) , temp(100)
  real*8 startim,entim
  real*8 rpt

  call MPI_Init(ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)
  startim = MPI_Wtime()
call init_random_seed()

  do j = 1,100

call MPI_PACK_SIZE(2,MPI_INTEGER,MPI_COMM_WORLD, sizeInt,IERR)
call MPI_PACK_SIZE(10,MPI_INTEGER,MPI_COMM_WORLD,sizeidx1, IERR)
call MPI_PACK_SIZE(100,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,sizeval, IERR)
sizebuf=sizeint+sizeidx1+sizeval
if (allocated(BUF)) deallocate(BUF)
allocate(BUF(sizebuf))

     if (my_rank.eq.0) then
        iseed = 2
        do i=1,p
            call RANDOM_NUMBER(rpt)
           initvec(i) = int(rpt*1000)
        end do
        ISIZE =  ( 26 + NJ + 2*KQUAD )
        POSITION = 0
        nj = 3
        kquad = 4 ! actually these might have been read from a file
        itemp(1) = isize
        itemp(2) = kquad
        do i=1,100
           temp(i) = 1.  ! more realistically we would read from a file
        end do

        CALL MPI_PACK ( ITEMP, 2, MPI_INTEGER, BUF, sizebuf,POSITION, MPI_COMM_WORLD, IERR)
        ! This call packs the integer vector itemp of length 2 into buf
        ! starting at position 0. It increments position.
        !--------------------------------------------------------------
        !  See below to see how to unpack.
        !  MPI_PACK and MPI_UNPACK are good for reducing the number of
        !  total calls.  These calls will allow us to pass multiple
        !  messages for the latency of one.  They are flexible in
        !  that the length of the unpack can be part of the message.
        !  The time of the calls to pack and unpack is not significant
        !  compared to the time to pass a message.
        !
        !  Disadvantage: On the Compaq AlphaServerS!the packing turned
        !        out to be pretty loose, i.e., there were empty bytes.  So combining
        !         short messages would speed thing , but long messages might get
        !          enough longer that they would take more time.
        !
        !--------------------------------------------------------------
        CALL MPI_PACK ( initvec, 10, MPI_INTEGER, BUF, sizebuf,POSITION, MPI_COMM_WORLD, IERR)
        CALL MPI_PACK ( TEMP, 100, MPI_DOUBLE_PRECISION, BUF, sizebuf,POSITION, MPI_COMM_WORLD, IERR)

     endif

     ! pack up data into one message


     root = 0
     call MPI_BCAST(BUF,sizebuf,MPI_PACKED,0,MPI_COMM_WORLD,IERR)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     IF ( my_rank.NE.0 ) THEN
        POSITION = 0
        !for the unpack, reset position to 0.  The unpack order is
        !   the same as the pack order.  But the order of arguments
        !   is changed.
        !The call below unpacks integer vector itemp of length 2 from
        !    the BUF buffer.
        !unpack itemp
        CALL MPI_UNPACK(BUF, sizebuf, POSITION,  ITEMP, 2, MPI_INTEGER,MPI_COMM_WORLD,IERR)
        isize  = ITEMP(1)
        kquad  = ITEMP(2)
        !unpack initvec
        CALL MPI_UNPACK (BUF, sizebuf, POSITION, initvec, 10,MPI_INTEGER, MPI_COMM_WORLD, IERR)
        myseed = initvec(my_rank)
        !unpack temp

        CALL MPI_UNPACK ( BUF, sizebuf, POSITION,  TEMP, 100,MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,IERR)
     ENDIF

     !-----------------------------------------------------------------------
     ! Left out -- a body of code that does a bunch of particle tracking
     !          stuff to produce the double precision vector ans
     !-----------------------------------------------------------------------
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     ans1 = rand(myseed)
     do i=1,10
        ans(i) = rand()  ! at least we initialize stuff to send back.
        ! But this call is something I had to change to get the code to
        ! run here.
     end do

     call MPI_REDUCE (ans,ans2, 10, MPI_DOUBLE_PRECISION,MPI_SUM, root, MPI_COMM_WORLD, ierr)
     !-------------------------------------------------------
     ! Get the (sum of) data back
     !      ans         -- arg1 -- message sent from each processor
     !      ans2        -- arg2 -- result deposited on root -- out
     !      10          -- arg3 -- length of message
     ! MPI_DOUBLE_PRECISION --arg4 - data type
     !  MPI_SUM         -- arg5 -- operation performed by reduce
     !      root        -- arg6 -- reduce deposits answer on root
     !                                same on all processors
     ! MPI_COMM_WORLD   -- arg7 -- all procs must have same communicator
     !      ierr        -- arg8 -- integer error--out (only in Fortran)
     !------------------------------------------------------
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     if(my_rank.eq.0) then
        !            do some stuff to process and output ans2
     endif
  end do
  entim =  MPI_Wtime() - startim

  call MPI_Finalize(ierr)
  print*,' elapsed time =',entim, '      my_rank=',my_rank
end program monte3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
