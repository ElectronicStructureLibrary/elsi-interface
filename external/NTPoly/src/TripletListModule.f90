

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Storing Lists of Triplets.
MODULE TripletListModule
  USE DataTypesModule, ONLY: NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX, &
       & MPINTINTEGER
  USE TripletModule, ONLY : Triplet_r, Triplet_c, CompareTriplets, &
       & ConvertTripletType
  USE MatrixMarketModule, ONLY : MM_SYMMETRIC, MM_SKEW_SYMMETRIC, MM_HERMITIAN
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  TYPE :: TripletList_r
     !> Internal representation of the data.
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: DATA
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
  END TYPE TripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A data type for a list of triplets.
  TYPE :: TripletList_c
     !> Internal representation of the data.
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: DATA
     !> Current number of elements in the triplet list
     INTEGER :: CurrentSize
  END TYPE TripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: TripletList_r
  PUBLIC :: TripletList_c
  PUBLIC :: ConstructTripletList
  PUBLIC :: DestructTripletList
  PUBLIC :: ResizeTripletList
  PUBLIC :: AppendToTripletList
  PUBLIC :: SetTripletAt
  PUBLIC :: GetTripletAt
  PUBLIC :: SortTripletList
  PUBLIC :: SymmetrizeTripletList
  PUBLIC :: GetTripletListSize
  PUBLIC :: RedistributeTripletLists
  PUBLIC :: ShiftTripletList
  PUBLIC :: ConvertTripletListType
  PUBLIC :: ConstructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE TripletList_r
     MODULE PROCEDURE ConstructTripletList_r
  END INTERFACE
  INTERFACE TripletList_c
     MODULE PROCEDURE ConstructTripletList_c
  END INTERFACE
  INTERFACE ConstructTripletList
     MODULE PROCEDURE ConstructTripletListSup_r
     MODULE PROCEDURE ConstructTripletListSup_c
  END INTERFACE
  INTERFACE DestructTripletList
     MODULE PROCEDURE DestructTripletList_r
     MODULE PROCEDURE DestructTripletList_c
  END INTERFACE
  INTERFACE ResizeTripletList
     MODULE PROCEDURE ResizeTripletList_r
     MODULE PROCEDURE ResizeTripletList_c
  END INTERFACE
  INTERFACE AppendToTripletList
     MODULE PROCEDURE AppendToTripletList_r
     MODULE PROCEDURE AppendToTripletList_c
  END INTERFACE
  INTERFACE SetTripletAt
     MODULE PROCEDURE SetTripletAt_r
     MODULE PROCEDURE SetTripletAt_c
  END INTERFACE
  INTERFACE GetTripletAt
     MODULE PROCEDURE GetTripletAt_r
     MODULE PROCEDURE GetTripletAt_c
  END INTERFACE
  INTERFACE SortTripletList
     MODULE PROCEDURE SortTripletList_r
     MODULE PROCEDURE SortTripletList_c
  END INTERFACE
  INTERFACE SortDenseTripletList
     MODULE PROCEDURE SortDenseTripletList_r
     MODULE PROCEDURE SortDenseTripletList_c
  END INTERFACE
  INTERFACE SymmetrizeTripletList
     MODULE PROCEDURE SymmetrizeTripletList_r
     MODULE PROCEDURE SymmetrizeTripletList_c
  END INTERFACE
  INTERFACE GetTripletListSize
     MODULE PROCEDURE GetTripletListSize_r
     MODULE PROCEDURE GetTripletListSize_c
  END INTERFACE
  INTERFACE RedistributeTripletLists
     MODULE PROCEDURE RedistributeTripletLists_r
     MODULE PROCEDURE RedistributeTripletLists_c
  END INTERFACE
  INTERFACE ShiftTripletList
     MODULE PROCEDURE ShiftTripletList_r
     MODULE PROCEDURE ShiftTripletList_c
  END INTERFACE
  INTERFACE ConvertTripletListType
     MODULE PROCEDURE ConvertTripletListToReal
     MODULE PROCEDURE ConvertTripletListToComplex
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for constructing a triplet list.
  PURE SUBROUTINE ConstructTripletListSup_r(this, size_in)
    !> The triplet list to construct.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> The length of the triplet list (default=0).
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    IF (PRESENT(size_in)) THEN
       this = ConstructTripletList_r(size_in)
    ELSE
       this = ConstructTripletList_r()
    END IF
  END SUBROUTINE ConstructTripletListSup_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for constructing a triplet list.
  PURE SUBROUTINE ConstructTripletListSup_c(this, size_in)
    !> The triplet list to construct.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> The length of the triplet list (default=0).
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

    IF (PRESENT(size_in)) THEN
       this = ConstructTripletList_c(size_in)
    ELSE
       this = ConstructTripletList_c()
    END IF
  END SUBROUTINE ConstructTripletListSup_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  PURE FUNCTION ConstructTripletList_r(size_in) RESULT(this)
    !> The triplet list to construct.
    TYPE(TripletList_r) :: this
    !> The length of the triplet list (default=0).
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

  !! Local data
  INTEGER :: size

  IF (PRESENT(size_in)) THEN
     size = size_in
  ELSE
     size = 0
  END IF

  IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
  this%CurrentSize = size

  ALLOCATE(this%data(size))

  END FUNCTION ConstructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a triplet list.
  PURE FUNCTION ConstructTripletList_c(size_in) RESULT(this)
    !> The triplet list to construct.
    TYPE(TripletList_c) :: this
    !> The length of the triplet list (default=0).
    INTEGER(kind=c_int), INTENT(IN), OPTIONAL :: size_in

  !! Local data
  INTEGER :: size

  IF (PRESENT(size_in)) THEN
     size = size_in
  ELSE
     size = 0
  END IF

  IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
  this%CurrentSize = size

  ALLOCATE(this%data(size))

  END FUNCTION ConstructTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  PURE SUBROUTINE DestructTripletList_r(this)
    !> The triplet list to destruct.
    TYPE(TripletList_r), INTENT(INOUT) :: this

  IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
  this%CurrentSize = 0

  END SUBROUTINE DestructTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructs a triplet list.
  PURE SUBROUTINE DestructTripletList_c(this)
    !> The triplet list to destruct.
    TYPE(TripletList_c), INTENT(INOUT) :: this

  IF (ALLOCATED(this%data)) DEALLOCATE(this%data)
  this%CurrentSize = 0

  END SUBROUTINE DestructTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  PURE SUBROUTINE ResizeTripletList_r(this, size)
    !> The triplet list to resize.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> Size to resize to.
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE :: temporary_data

  !! Temporary copy
  ALLOCATE(temporary_data(this%CurrentSize))
  temporary_data = this%DATA(:this%CurrentSize)

  !! Create new memory
  IF (ALLOCATED(this%DATA)) DEALLOCATE(this%DATA)
  ALLOCATE(this%DATA(size))

  !! Copy back
  this%DATA(:this%CurrentSize) = temporary_data

  !! Cleanup
  DEALLOCATE(temporary_data)

  END SUBROUTINE ResizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Increase the size of a triplet list.
  PURE SUBROUTINE ResizeTripletList_c(this, size)
    !> The triplet list to resize.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> Size to resize to.
    INTEGER(KIND=c_int), INTENT(IN) :: size
    !! Local Data
    TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE :: temporary_data

  !! Temporary copy
  ALLOCATE(temporary_data(this%CurrentSize))
  temporary_data = this%DATA(:this%CurrentSize)

  !! Create new memory
  IF (ALLOCATED(this%DATA)) DEALLOCATE(this%DATA)
  ALLOCATE(this%DATA(size))

  !! Copy back
  this%DATA(:this%CurrentSize) = temporary_data

  !! Cleanup
  DEALLOCATE(temporary_data)

  END SUBROUTINE ResizeTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  PURE SUBROUTINE AppendToTripletList_r(this, triplet_value)
    !> This the triplet list to append to.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> The value to append.
    TYPE(Triplet_r), INTENT(IN)        :: triplet_value

  !! Local data
  INTEGER :: new_size

  !! First, check if we need to allocate more memory
  IF (this%CurrentSize+1 .GT. SIZE(this%DATA)) THEN
     IF (SIZE(this%data) .EQ. 0) THEN
        new_size = 1
     ELSE IF (SIZE(this%data) .EQ. 1) THEN
        new_size = 2
     ELSE
        new_size = INT(SIZE(this%data)*1.5)
     END IF
     CALL ResizeTripletList(this,new_size)
  END IF

  !! Append
  this%CurrentSize = this%CurrentSize+1
  this%DATA(this%CurrentSize) = triplet_value

  END SUBROUTINE AppendToTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add a value to the end of the triplet list.
  PURE SUBROUTINE AppendToTripletList_c(this, triplet_value)
    !> This the triplet list to append to.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> The value to append.
    TYPE(Triplet_c), INTENT(IN)        :: triplet_value

  !! Local data
  INTEGER :: new_size

  !! First, check if we need to allocate more memory
  IF (this%CurrentSize+1 .GT. SIZE(this%DATA)) THEN
     IF (SIZE(this%data) .EQ. 0) THEN
        new_size = 1
     ELSE IF (SIZE(this%data) .EQ. 1) THEN
        new_size = 2
     ELSE
        new_size = INT(SIZE(this%data)*1.5)
     END IF
     CALL ResizeTripletList(this,new_size)
  END IF

  !! Append
  this%CurrentSize = this%CurrentSize+1
  this%DATA(this%CurrentSize) = triplet_value

  END SUBROUTINE AppendToTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  PURE SUBROUTINE SetTripletAt_r(this,index,triplet_value)
    !> The triplet list to set.
    TYPE(TripletList_r), INTENT(INOUT) :: this
    !> The index at which to set the triplet.
    INTEGER(KIND=c_int), INTENT(IN)    :: index
    !> The value of the triplet to set.
    TYPE(Triplet_r), INTENT(IN)        :: triplet_value

  this%data(index) = triplet_value
  END SUBROUTINE SetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Set the value of a triplet at a particular index.
  PURE SUBROUTINE SetTripletAt_c(this,index,triplet_value)
    !> The triplet list to set.
    TYPE(TripletList_c), INTENT(INOUT) :: this
    !> The index at which to set the triplet.
    INTEGER(KIND=c_int), INTENT(IN)    :: index
    !> The value of the triplet to set.
    TYPE(Triplet_c), INTENT(IN)        :: triplet_value

  this%data(index) = triplet_value
  END SUBROUTINE SetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  PURE SUBROUTINE GetTripletAt_r(this,index,triplet_value)
    !> The triplet list to get the value from.
    TYPE(TripletList_r), INTENT(IN) :: this
    !> The index from which to get the triplet.
    INTEGER(kind=c_int), INTENT(IN) :: index
    !> The extracted triplet value.
    TYPE(Triplet_r), INTENT(OUT)    :: triplet_value

  triplet_value = this%data(index)
  END SUBROUTINE GetTripletAt_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the value of a triplet at a particular index.
  PURE SUBROUTINE GetTripletAt_c(this,index,triplet_value)
    !> The triplet list to get the value from.
    TYPE(TripletList_c), INTENT(IN) :: this
    !> The index from which to get the triplet.
    INTEGER(kind=c_int), INTENT(IN) :: index
    !> The extracted triplet value.
    TYPE(Triplet_c), INTENT(OUT)    :: triplet_value

  triplet_value = this%data(index)
  END SUBROUTINE GetTripletAt_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !> Implementation is based on bucket sort. This is why it needs the number of
  !> matrix columns. Bubble sort is used within a bucket.
  PURE SUBROUTINE SortTripletList_r(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !> List to be sorted.
    TYPE(TripletList_r), INTENT(IN)  :: input_list
    !> This is the highest column value in the list.
    INTEGER, INTENT(IN) :: matrix_columns
    !> This is the highest row value in the list.
    INTEGER, INTENT(IN) :: matrix_rows
    !> A now sorted version of the list. This routine will allocate it.
    TYPE(TripletList_r), INTENT(OUT) :: sorted_list
    !> False if you do not need the final bubble sort.
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_r) :: temporary

  !! Local Data
  LOGICAL :: bubble
  LOGICAL :: swap_occured
  INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
  INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
  !! Counters and temporary variables
  INTEGER :: counter
  INTEGER :: temp_index
  INTEGER :: alloc_stat
  INTEGER :: list_length

  IF (PRESENT(bubble_in)) THEN
     bubble = bubble_in
  ELSE
     bubble = .TRUE.
  END IF

  list_length = input_list%CurrentSize

  IF (bubble .AND. list_length .GT. matrix_rows*matrix_columns*0.1) THEN
     CALL SortDenseTripletList(input_list, matrix_columns, matrix_rows, &
          & sorted_list)
  ELSE
     !! Data Allocation
     CALL ConstructTripletList(sorted_list, list_length)
     ALLOCATE(values_per_row(matrix_columns), stat=alloc_stat)
     ALLOCATE(offset_array(matrix_columns), stat=alloc_stat)
     ALLOCATE(inserted_per_row(matrix_columns), stat=alloc_stat)

     !! Initial one dimensional sort
     values_per_row = 0
     inserted_per_row = 0

     !! Do a first pass bucket sort
     DO counter = 1, input_list%CurrentSize
        values_per_row(input_list%data(counter)%index_column) = &
             & values_per_row(input_list%data(counter)%index_column) + 1
     END DO
     offset_array(1) = 1
     DO counter = 2, UBOUND(offset_array,dim=1)
        offset_array(counter) = offset_array(counter-1) + &
             & values_per_row(counter-1)
     END DO
     DO counter = 1, input_list%CurrentSize
        temp_index = input_list%data(counter)%index_column
        sorted_list%data(offset_array(temp_index)+inserted_per_row(temp_index))=&
             & input_list%data(counter)
        inserted_per_row(temp_index) = inserted_per_row(temp_index) + 1
     END DO

     !! Finish with bubble sort
     !! Not necessary for transposing or unpacking.
     swap_occured = .TRUE.
     IF (bubble) THEN
        DO WHILE (swap_occured .EQV. .TRUE.)
           swap_occured = .FALSE.
           DO counter = 2, sorted_list%CurrentSize
              IF (CompareTriplets(sorted_list%data(counter-1), &
                   & sorted_list%data(counter))) THEN
                 temporary = sorted_list%data(counter)
                 sorted_list%data(counter) = sorted_list%data(counter-1)
                 sorted_list%data(counter-1) = temporary
                 swap_occured = .TRUE.
              END IF
           END DO
        END DO
     END IF

     !! Cleanup
     DEALLOCATE(values_per_row)
     DEALLOCATE(offset_array)
     DEALLOCATE(inserted_per_row)
  END IF

  END SUBROUTINE SortTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sorts a triplet list by index values.
  !> Implementation is based on bucket sort. This is why it needs the number of
  !> matrix columns. Bubble sort is used within a bucket.
  PURE SUBROUTINE SortTripletList_c(input_list, matrix_columns, matrix_rows, &
       & sorted_list, bubble_in)
    !> List to be sorted.
    TYPE(TripletList_c), INTENT(IN)  :: input_list
    !> This is the highest column value in the list.
    INTEGER, INTENT(IN) :: matrix_columns
    !> This is the highest row value in the list.
    INTEGER, INTENT(IN) :: matrix_rows
    !> A now sorted version of the list. This routine will allocate it.
    TYPE(TripletList_c), INTENT(OUT) :: sorted_list
    !> False if you do not need the final bubble sort.
    LOGICAL, OPTIONAL, INTENT(IN) :: bubble_in
    !! Local Data
    TYPE(Triplet_c) :: temporary

  !! Local Data
  LOGICAL :: bubble
  LOGICAL :: swap_occured
  INTEGER, DIMENSION(:), ALLOCATABLE :: values_per_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: offset_array
  INTEGER, DIMENSION(:), ALLOCATABLE :: inserted_per_row
  !! Counters and temporary variables
  INTEGER :: counter
  INTEGER :: temp_index
  INTEGER :: alloc_stat
  INTEGER :: list_length

  IF (PRESENT(bubble_in)) THEN
     bubble = bubble_in
  ELSE
     bubble = .TRUE.
  END IF

  list_length = input_list%CurrentSize

  IF (bubble .AND. list_length .GT. matrix_rows*matrix_columns*0.1) THEN
     CALL SortDenseTripletList(input_list, matrix_columns, matrix_rows, &
          & sorted_list)
  ELSE
     !! Data Allocation
     CALL ConstructTripletList(sorted_list, list_length)
     ALLOCATE(values_per_row(matrix_columns), stat=alloc_stat)
     ALLOCATE(offset_array(matrix_columns), stat=alloc_stat)
     ALLOCATE(inserted_per_row(matrix_columns), stat=alloc_stat)

     !! Initial one dimensional sort
     values_per_row = 0
     inserted_per_row = 0

     !! Do a first pass bucket sort
     DO counter = 1, input_list%CurrentSize
        values_per_row(input_list%data(counter)%index_column) = &
             & values_per_row(input_list%data(counter)%index_column) + 1
     END DO
     offset_array(1) = 1
     DO counter = 2, UBOUND(offset_array,dim=1)
        offset_array(counter) = offset_array(counter-1) + &
             & values_per_row(counter-1)
     END DO
     DO counter = 1, input_list%CurrentSize
        temp_index = input_list%data(counter)%index_column
        sorted_list%data(offset_array(temp_index)+inserted_per_row(temp_index))=&
             & input_list%data(counter)
        inserted_per_row(temp_index) = inserted_per_row(temp_index) + 1
     END DO

     !! Finish with bubble sort
     !! Not necessary for transposing or unpacking.
     swap_occured = .TRUE.
     IF (bubble) THEN
        DO WHILE (swap_occured .EQV. .TRUE.)
           swap_occured = .FALSE.
           DO counter = 2, sorted_list%CurrentSize
              IF (CompareTriplets(sorted_list%data(counter-1), &
                   & sorted_list%data(counter))) THEN
                 temporary = sorted_list%data(counter)
                 sorted_list%data(counter) = sorted_list%data(counter-1)
                 sorted_list%data(counter-1) = temporary
                 swap_occured = .TRUE.
              END IF
           END DO
        END DO
     END IF

     !! Cleanup
     DEALLOCATE(values_per_row)
     DEALLOCATE(offset_array)
     DEALLOCATE(inserted_per_row)
  END IF

  END SUBROUTINE SortTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  PURE FUNCTION GetTripletListSize_r(triplet_list) RESULT(list_size)
    !> List to get the size of.
    TYPE(TripletList_r), INTENT(IN)  :: triplet_list
    !> The number of entries in the triplet list.
    INTEGER :: list_size

  list_size = triplet_list%CurrentSize

  END FUNCTION GetTripletListSize_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Get the number of entries in a triplet list.
  PURE FUNCTION GetTripletListSize_c(triplet_list) RESULT(list_size)
    !> List to get the size of.
    TYPE(TripletList_c), INTENT(IN)  :: triplet_list
    !> The number of entries in the triplet list.
    INTEGER :: list_size

  list_size = triplet_list%CurrentSize

  END FUNCTION GetTripletListSize_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !> Takes in a list of triplet lists, one list for each processor. Then the
  !> all to all redistribution is performed along the given communicator.
  SUBROUTINE RedistributeTripletLists_r(triplet_lists, comm, local_data_out)
    !> A list of triplet lists, one for each process.
    TYPE(TripletList_r), DIMENSION(:), INTENT(IN) :: triplet_lists
    !> The mpi communicator to redistribute along.
    INTEGER, INTENT(INOUT) :: comm
    !> The resulting local triplet list.
    TYPE(TripletList_r), INTENT(INOUT) :: local_data_out
    !! Local data (type specific)
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_r) :: temp_triplet

  !! Local Data - Offsets
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_process
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_offsets
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_process
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_offsets
  !! Local Data - Send/Recv Buffers
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
  !! ETC
  INTEGER :: num_processes
  INTEGER :: counter, inner_counter, insert_pt
  INTEGER :: mpi_error

  !! Allocate Size Buffers
  CALL MPI_COMM_SIZE(comm, num_processes, mpi_error)
  ALLOCATE(send_per_process(num_processes))
  ALLOCATE(send_offsets(num_processes))
  ALLOCATE(recv_per_process(num_processes))
  ALLOCATE(recv_offsets(num_processes))

  !! Figure Out How Much Data Gets Sent
  DO counter = 1, num_processes
     send_per_process(counter) = triplet_lists(counter)%CurrentSize
  END DO
  send_offsets(1) = 0
  DO counter = 2, num_processes
     send_offsets(counter) = send_offsets(counter-1) + &
          & send_per_process(counter-1)
  END DO

  !! Figure Out How Much Data Gets Received
  CALL MPI_ALLTOALL(send_per_process, 1, MPINTINTEGER, recv_per_process, 1, &
       & MPINTINTEGER, comm, mpi_error)
  recv_offsets(1) = 0
  DO counter = 2, num_processes
     recv_offsets(counter) = recv_offsets(counter-1) + &
          & recv_per_process(counter-1)
  END DO

  !! Allocate And Fill Send Buffers
  ALLOCATE(send_buffer_row(SUM(send_per_process)))
  ALLOCATE(send_buffer_col(SUM(send_per_process)))
  ALLOCATE(send_buffer_val(SUM(send_per_process)))
  ALLOCATE(recv_buffer_row(SUM(recv_per_process)))
  ALLOCATE(recv_buffer_col(SUM(recv_per_process)))
  ALLOCATE(recv_buffer_val(SUM(recv_per_process)))

  !! Fill Send Buffer
  insert_pt = 1
  DO counter = 1, num_processes
     DO inner_counter = 1, triplet_lists(counter)%CurrentSize
        CALL GetTripletAt(triplet_lists(counter), inner_counter, temp_triplet)
        send_buffer_row(insert_pt) = temp_triplet%index_row
        send_buffer_col(insert_pt) = temp_triplet%index_column
        send_buffer_val(insert_pt) = temp_triplet%point_value
        insert_pt = insert_pt + 1
     END DO
  END DO

  !! Do Actual Send
  CALL MPI_Alltoallv(send_buffer_col, send_per_process, send_offsets, &
       & MPINTINTEGER, recv_buffer_col, recv_per_process, recv_offsets, &
       & MPINTINTEGER, comm, mpi_error)
  CALL MPI_Alltoallv(send_buffer_row, send_per_process, send_offsets, &
       & MPINTINTEGER, recv_buffer_row, recv_per_process, recv_offsets, &
       & MPINTINTEGER, comm, mpi_error)
  CALL MPI_Alltoallv(send_buffer_val, send_per_process, send_offsets, &
       & MPINTREAL, recv_buffer_val, recv_per_process, recv_offsets, &
       & MPINTREAL, comm, mpi_error)

  !! Unpack Into The Output Triplet List
  CALL ConstructTripletList(local_data_out, size_in=SUM(recv_per_process))
  DO counter = 1, SUM(recv_per_process)
     local_data_out%data(counter)%index_column = recv_buffer_col(counter)
     local_data_out%data(counter)%index_row = recv_buffer_row(counter)
     local_data_out%data(counter)%point_value = recv_buffer_val(counter)
  END DO

  !! Cleanup
  DEALLOCATE(send_per_process)
  DEALLOCATE(send_offsets)
  DEALLOCATE(recv_per_process)
  DEALLOCATE(recv_offsets)
  DEALLOCATE(send_buffer_row)
  DEALLOCATE(send_buffer_col)
  DEALLOCATE(send_buffer_val)
  DEALLOCATE(recv_buffer_row)
  DEALLOCATE(recv_buffer_col)
  DEALLOCATE(recv_buffer_val)

  END SUBROUTINE RedistributeTripletLists_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Redistribute some triplet lists amongst a set of processors.
  !> Takes in a list of triplet lists, one list for each processor. Then the
  !> all to all redistribution is performed along the given communicator.
  SUBROUTINE RedistributeTripletLists_c(triplet_lists, comm, local_data_out)
    !> A list of triplet lists, one for each process.
    TYPE(TripletList_c), DIMENSION(:), INTENT(IN) :: triplet_lists
    !> The mpi communicator to redistribute along.
    INTEGER, INTENT(INOUT) :: comm
    !> The resulting local triplet list.
    TYPE(TripletList_c), INTENT(INOUT) :: local_data_out
    !! Local data (type specific)
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: send_buffer_val
    COMPLEX(NTCOMPLEX), DIMENSION(:), ALLOCATABLE :: recv_buffer_val
    TYPE(Triplet_c) :: temp_triplet

  !! Local Data - Offsets
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_per_process
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_offsets
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_per_process
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_offsets
  !! Local Data - Send/Recv Buffers
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: send_buffer_col
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_row
  INTEGER, DIMENSION(:), ALLOCATABLE :: recv_buffer_col
  !! ETC
  INTEGER :: num_processes
  INTEGER :: counter, inner_counter, insert_pt
  INTEGER :: mpi_error

  !! Allocate Size Buffers
  CALL MPI_COMM_SIZE(comm, num_processes, mpi_error)
  ALLOCATE(send_per_process(num_processes))
  ALLOCATE(send_offsets(num_processes))
  ALLOCATE(recv_per_process(num_processes))
  ALLOCATE(recv_offsets(num_processes))

  !! Figure Out How Much Data Gets Sent
  DO counter = 1, num_processes
     send_per_process(counter) = triplet_lists(counter)%CurrentSize
  END DO
  send_offsets(1) = 0
  DO counter = 2, num_processes
     send_offsets(counter) = send_offsets(counter-1) + &
          & send_per_process(counter-1)
  END DO

  !! Figure Out How Much Data Gets Received
  CALL MPI_ALLTOALL(send_per_process, 1, MPINTINTEGER, recv_per_process, 1, &
       & MPINTINTEGER, comm, mpi_error)
  recv_offsets(1) = 0
  DO counter = 2, num_processes
     recv_offsets(counter) = recv_offsets(counter-1) + &
          & recv_per_process(counter-1)
  END DO

  !! Allocate And Fill Send Buffers
  ALLOCATE(send_buffer_row(SUM(send_per_process)))
  ALLOCATE(send_buffer_col(SUM(send_per_process)))
  ALLOCATE(send_buffer_val(SUM(send_per_process)))
  ALLOCATE(recv_buffer_row(SUM(recv_per_process)))
  ALLOCATE(recv_buffer_col(SUM(recv_per_process)))
  ALLOCATE(recv_buffer_val(SUM(recv_per_process)))

  !! Fill Send Buffer
  insert_pt = 1
  DO counter = 1, num_processes
     DO inner_counter = 1, triplet_lists(counter)%CurrentSize
        CALL GetTripletAt(triplet_lists(counter), inner_counter, temp_triplet)
        send_buffer_row(insert_pt) = temp_triplet%index_row
        send_buffer_col(insert_pt) = temp_triplet%index_column
        send_buffer_val(insert_pt) = temp_triplet%point_value
        insert_pt = insert_pt + 1
     END DO
  END DO

  !! Do Actual Send
  CALL MPI_Alltoallv(send_buffer_col, send_per_process, send_offsets, &
       & MPINTINTEGER, recv_buffer_col, recv_per_process, recv_offsets, &
       & MPINTINTEGER, comm, mpi_error)
  CALL MPI_Alltoallv(send_buffer_row, send_per_process, send_offsets, &
       & MPINTINTEGER, recv_buffer_row, recv_per_process, recv_offsets, &
       & MPINTINTEGER, comm, mpi_error)
  CALL MPI_Alltoallv(send_buffer_val, send_per_process, send_offsets, &
       & MPINTCOMPLEX, recv_buffer_val, recv_per_process, recv_offsets, &
       & MPINTCOMPLEX, comm, mpi_error)

  !! Unpack Into The Output Triplet List
  CALL ConstructTripletList(local_data_out, size_in=SUM(recv_per_process))
  DO counter = 1, SUM(recv_per_process)
     local_data_out%data(counter)%index_column = recv_buffer_col(counter)
     local_data_out%data(counter)%index_row = recv_buffer_row(counter)
     local_data_out%data(counter)%point_value = recv_buffer_val(counter)
  END DO

  !! Cleanup
  DEALLOCATE(send_per_process)
  DEALLOCATE(send_offsets)
  DEALLOCATE(recv_per_process)
  DEALLOCATE(recv_offsets)
  DEALLOCATE(send_buffer_row)
  DEALLOCATE(send_buffer_col)
  DEALLOCATE(send_buffer_val)
  DEALLOCATE(recv_buffer_row)
  DEALLOCATE(recv_buffer_col)
  DEALLOCATE(recv_buffer_val)

  END SUBROUTINE RedistributeTripletLists_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !> Frequently, we have a triplet list that comes from the global matrix which
  !> we would like to shift into a local matrix. In that case, just pass
  !> the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_r(triplet_list, row_shift, column_shift)
    !> The triplet list to shift.
    TYPE(TripletList_r), INTENT(INOUT) :: triplet_list
    !> The row offset to shift by.
    INTEGER, INTENT(IN) :: row_shift
    !> The column offset to shift by.
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: counter

  !! Loop
  DO counter = 1, triplet_list%CurrentSize
     triplet_list%data(counter)%index_row = &
          triplet_list%data(counter)%index_row + row_shift
     triplet_list%data(counter)%index_column = &
          triplet_list%data(counter)%index_column + column_shift
  END DO

  END SUBROUTINE ShiftTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Shift the rows and columns of a triplet list by set values.
  !> Frequently, we have a triplet list that comes from the global matrix which
  !> we would like to shift into a local matrix. In that case, just pass
  !> the negative of the starting row and column (plus 1) to this routine.
  PURE SUBROUTINE ShiftTripletList_c(triplet_list, row_shift, column_shift)
    !> The triplet list to shift.
    TYPE(TripletList_c), INTENT(INOUT) :: triplet_list
    !> The row offset to shift by.
    INTEGER, INTENT(IN) :: row_shift
    !> The column offset to shift by.
    INTEGER, INTENT(IN) :: column_shift
    !! Local Variables
    INTEGER :: counter

  !! Loop
  DO counter = 1, triplet_list%CurrentSize
     triplet_list%data(counter)%index_row = &
          triplet_list%data(counter)%index_row + row_shift
     triplet_list%data(counter)%index_column = &
          triplet_list%data(counter)%index_column + column_shift
  END DO

  END SUBROUTINE ShiftTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !> dense.
  PURE SUBROUTINE SortDenseTripletList_r(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !> The list to sort.
    TYPE(TripletList_r), INTENT(IN)  :: input_list
    !> Number of columns for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_columns
    !> Number of rows for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_rows
    !> Sorted and ready to use for building matrices.
    TYPE(TripletList_r), INTENT(OUT) :: sorted_list
    !! Local Variables
    REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: value_buffer

  !! Local Data
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: dirty_buffer
  INTEGER :: list_length
  INTEGER :: row, col, ind
  INTEGER :: II, JJ

  !! Setup Memory
  ALLOCATE(value_buffer(matrix_rows,matrix_columns))
  ALLOCATE(dirty_buffer(matrix_rows,matrix_columns))
  value_buffer = 0
  dirty_buffer = 0
  list_length = input_list%CurrentSize
  CALL ConstructTripletList(sorted_list, list_length)

  !! Unpack
  DO II = 1, list_length
     row = input_list%data(II)%index_row
     col = input_list%data(II)%index_column
     value_buffer(row,col) = input_list%data(II)%point_value
     dirty_buffer(row,col) = 1
  END DO

  !! Repack
  ind = 1
  DO JJ = 1, matrix_columns
     DO II = 1, matrix_rows
        IF (dirty_buffer(II,JJ) .EQ. 1) THEN
           sorted_list%data(ind)%index_row = II
           sorted_list%data(ind)%index_column = JJ
           sorted_list%data(ind)%point_value = value_buffer(II,JJ)
           ind = ind + 1
        END IF
     END DO
  END DO

  !! Cleanup
  DEALLOCATE(value_buffer)
  DEALLOCATE(dirty_buffer)

  END SUBROUTINE SortDenseTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sort a triplet list assuming that the matrix it corresponds to is nearly
  !> dense.
  PURE SUBROUTINE SortDenseTripletList_c(input_list, matrix_columns, &
       & matrix_rows, sorted_list)
    !> The list to sort.
    TYPE(TripletList_c), INTENT(IN)  :: input_list
    !> Number of columns for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_columns
    !> Number of rows for the corresponding matrix.
    INTEGER, INTENT(IN) :: matrix_rows
    !> Sorted and ready to use for building matrices.
    TYPE(TripletList_c), INTENT(OUT) :: sorted_list
    !! Local Variables
    COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: value_buffer

  !! Local Data
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: dirty_buffer
  INTEGER :: list_length
  INTEGER :: row, col, ind
  INTEGER :: II, JJ

  !! Setup Memory
  ALLOCATE(value_buffer(matrix_rows,matrix_columns))
  ALLOCATE(dirty_buffer(matrix_rows,matrix_columns))
  value_buffer = 0
  dirty_buffer = 0
  list_length = input_list%CurrentSize
  CALL ConstructTripletList(sorted_list, list_length)

  !! Unpack
  DO II = 1, list_length
     row = input_list%data(II)%index_row
     col = input_list%data(II)%index_column
     value_buffer(row,col) = input_list%data(II)%point_value
     dirty_buffer(row,col) = 1
  END DO

  !! Repack
  ind = 1
  DO JJ = 1, matrix_columns
     DO II = 1, matrix_rows
        IF (dirty_buffer(II,JJ) .EQ. 1) THEN
           sorted_list%data(ind)%index_row = II
           sorted_list%data(ind)%index_column = JJ
           sorted_list%data(ind)%point_value = value_buffer(II,JJ)
           ind = ind + 1
        END IF
     END DO
  END DO

  !! Cleanup
  DEALLOCATE(value_buffer)
  DEALLOCATE(dirty_buffer)

  END SUBROUTINE SortDenseTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !> symmetry type.
  SUBROUTINE SymmetrizeTripletList_r(triplet_list, pattern_type)
    !> List to be symmetrized.
    TYPE(TripletList_r), INTENT(INOUT)  :: triplet_list
    !> Type of symmetry.
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(Triplet_r) :: temporary, temporary_transpose
    INTEGER :: counter
    INTEGER :: initial_size

    initial_size = triplet_list%CurrentSize
    SELECT CASE(pattern_type)
    CASE(MM_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    CASE(MM_SKEW_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = -1.0*temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    END SELECT
  END SUBROUTINE SymmetrizeTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Symmetrizes an unsymmetric triplet list according to the specified
  !> symmetry type.
  SUBROUTINE SymmetrizeTripletList_c(triplet_list, pattern_type)
    !> List to be symmetrized.
    TYPE(TripletList_c), INTENT(INOUT)  :: triplet_list
    !> Type of symmetry.
    INTEGER, INTENT(IN) :: pattern_type
    !! Local variables
    TYPE(Triplet_c) :: temporary, temporary_transpose
    INTEGER :: counter
    INTEGER :: initial_size

    initial_size = triplet_list%CurrentSize
    SELECT CASE(pattern_type)
    CASE(MM_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    CASE(MM_HERMITIAN)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = CONJG(temporary%point_value)
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    CASE(MM_SKEW_SYMMETRIC)
       DO counter = 1, initial_size
          CALL GetTripletAt(triplet_list,counter,temporary)
          IF (temporary%index_column .NE. temporary%index_row) THEN
             temporary_transpose%index_row = temporary%index_column
             temporary_transpose%index_column = temporary%index_row
             temporary_transpose%point_value = -1.0*temporary%point_value
             CALL AppendToTripletList(triplet_list,temporary_transpose)
          END IF
       END DO
    END SELECT
  END SUBROUTINE SymmetrizeTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a complex triplet list to a real triplet list.
  SUBROUTINE ConvertTripletListToReal(cin_triplet, rout_triplet)
    !> The starting triplet list.
    TYPE(TripletList_c), INTENT(IN)    :: cin_triplet
    !> Real valued triplet list.
    TYPE(TripletList_r), INTENT(INOUT) :: rout_triplet
    !! Local Variables
    INTEGER :: II

    CALL ConstructTripletList(rout_triplet, cin_triplet%CurrentSize)
    DO II = 1, cin_triplet%CurrentSize
       CALL ConvertTripletType(cin_triplet%data(II), rout_triplet%data(II))
    END DO

  END SUBROUTINE ConvertTripletListToReal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a real triplet to a complex triplet list.
  SUBROUTINE ConvertTripletListToComplex(rin_triplet, cout_triplet)
    !> The starting triplet list.
    TYPE(TripletList_r), INTENT(IN)    :: rin_triplet
    !> Complex valued triplet list.
    TYPE(TripletList_c), INTENT(INOUT) :: cout_triplet
    !! Local Variables
    INTEGER :: II

    CALL ConstructTripletList(cout_triplet, rin_triplet%CurrentSize)
    DO II = 1, rin_triplet%CurrentSize
       CALL ConvertTripletType(rin_triplet%data(II), cout_triplet%data(II))
    END DO

  END SUBROUTINE ConvertTripletListToComplex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE TripletListModule
