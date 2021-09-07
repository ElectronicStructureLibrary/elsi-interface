!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Helper Routines for Computing The Cholesky Decomposition
MODULE CholeskyModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, MPINTINTEGER
  USE DMatrixModule, ONLY : Matrix_ldr
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceAndComposeMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, FillMatrixFromTripletList
  USE ProcessGridModule, ONLY : ProcessGrid_t
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, TransposeMatrix, &
       & DestructMatrix
  USE SVectorModule, ONLY : DotSparseVectors
  USE TripletListModule, ONLY : TripletList_r, AppendToTripletList, &
       & DestructTripletList
  USE TripletModule, ONLY : Triplet_r
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: AppendToVector
  PUBLIC :: BroadcastVector
  PUBLIC :: ConstructDiag
  PUBLIC :: ConstructRankLookup
  PUBLIC :: DotAllHelper
  PUBLIC :: DotAllPivoted
  PUBLIC :: GatherMatrixColumn
  PUBLIC :: GetPivot
  PUBLIC :: UnpackCholesky
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE AppendToVector
     MODULE PROCEDURE AppendToVector_r
  END INTERFACE AppendToVector
  INTERFACE BroadcastVector
     MODULE PROCEDURE BroadcastVector_r
  END INTERFACE BroadcastVector
  INTERFACE ConstructDiag
     MODULE PROCEDURE ConstructDiag_r
  END INTERFACE ConstructDiag
  INTERFACE DotAllHelper
     MODULE PROCEDURE DotAllHelper_r
  END INTERFACE DotAllHelper
  INTERFACE DotAllPivoted
     MODULE PROCEDURE DotAllPivoted_r
  END INTERFACE DotAllPivoted
  INTERFACE GatherMatrixColumn
     MODULE PROCEDURE GatherMatrixColumn_r
  END INTERFACE GatherMatrixColumn
  INTERFACE GetPivot
     MODULE PROCEDURE GetPivot_r
  END INTERFACE GetPivot
  INTERFACE UnpackCholesky
     MODULE PROCEDURE UnpackCholesky_r
  END INTERFACE UnpackCholesky
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to insert a value into a sparse vector.
  PURE SUBROUTINE AppendToVector_r(values_per, indices, values, insert_row, &
       & insert_value)
    !> Values per row.
    INTEGER, INTENT(INOUT) :: values_per
    !> Indices associated with each value.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    !> Values.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    !> Row to insert into.
    INTEGER, INTENT(IN) :: insert_row
    !> Value to insert.
    REAL(NTREAL), INTENT(IN) :: insert_value


    values_per = values_per + 1
    indices(values_per) = insert_row
    values(values_per) = insert_value
  END SUBROUTINE AppendToVector_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine to broadcast a sparse vector
  SUBROUTINE BroadcastVector_r(num_values, indices, values, root, comm)
    !> Number of values we are broadcasting.
    INTEGER, INTENT(INOUT) :: num_values
    !> Indices to broadcast.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
    !> Values to broadcast.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: values
    !> Root from which we broadcast.
    INTEGER, INTENT(IN) :: root
    !> Communicator to broadcast along.
    INTEGER, INTENT(INOUT) :: comm
    !! Local
    INTEGER :: err


    CALL MPI_Bcast(num_values, 1, MPINTINTEGER, root, comm, err)
    CALL MPI_Bcast(indices(:num_values), num_values, MPINTINTEGER, root, &
         & comm, err)
    CALL MPI_Bcast(values(:num_values), num_values, MPINTREAL, root, comm, err)

  END SUBROUTINE BroadcastVector_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct the vector holding the accumulated diagonal values
  SUBROUTINE ConstructDiag_r(AMat, process_grid, dense_a, diag)
    !> AMat the matrix we're working on (for meta data).
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> The process grid we're operating on.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !> A dense representation of the values.
    TYPE(Matrix_ldr), INTENT(IN) :: dense_a
    !> Diagonal values computed.
    REAL(NTREAL), DIMENSION(:), INTENT(INOUT) :: diag


    !! Local Variables
    INTEGER :: fill_counter
    INTEGER :: II, JJ, local_JJ, local_II
    INTEGER, DIMENSION(process_grid%num_process_rows) :: diags_per_proc
    INTEGER, DIMENSION(process_grid%num_process_rows) :: diag_displ
    INTEGER :: ierr

    diag = 0
    fill_counter = 0
    DO JJ = AMat%start_column, AMat%end_column - 1
       IF (JJ .GE. AMat%start_row .AND. JJ .LT. AMat%end_row) THEN
          local_JJ = JJ - AMat%start_column + 1
          local_II = JJ - AMat%start_row + 1
          diag(local_JJ) = dense_a%DATA(local_II, local_JJ)
          fill_counter = fill_counter + 1
       END IF
    END DO
    diags_per_proc(process_grid%my_row+1) = fill_counter

    !! Duplicate the diagonal entries along the process column (across rows)
    CALL MPI_Allgather(MPI_IN_PLACE, 1, MPINTINTEGER, diags_per_proc, 1, &
         & MPINTINTEGER, process_grid%column_comm, ierr)
    diag_displ(1) = 0
    DO II = 2, process_grid%num_process_rows
       diag_displ(II) = diag_displ(II-1) + diags_per_proc(II-1)
    END DO
    CALL MPI_Allgatherv(MPI_IN_PLACE, diags_per_proc(process_grid%my_row+1), &
         & MPINTREAL, diag, diags_per_proc, diag_displ, MPINTREAL, &
         & process_grid%column_comm, ierr)
  END SUBROUTINE ConstructDiag_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct a lookup for columns
  SUBROUTINE ConstructRankLookup(AMat, process_grid, col_root_lookup)
    !> Matrix we're computing.
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> Grid we're computing along.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !> The lookup we're computing.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: col_root_lookup
    !! Local Variables
    INTEGER, DIMENSION(process_grid%num_process_columns) :: cols_per_proc
    INTEGER, DIMENSION(process_grid%num_process_columns) :: d_cols_per_proc
    INTEGER :: II
    INTEGER :: ierr

    cols_per_proc(process_grid%my_column+1) = AMat%local_columns
    CALL MPI_Allgather(MPI_IN_PLACE, 1, MPI_INTEGER, cols_per_proc, 1, &
         & MPI_INTEGER, process_grid%row_comm, ierr)
    d_cols_per_proc(1) = 0
    DO II = 2, process_grid%num_process_columns
       d_cols_per_proc(II) = d_cols_per_proc(II-1) + cols_per_proc(II-1)
    END DO
    col_root_lookup(AMat%start_column:AMat%end_column - 1) =  &
         & process_grid%row_rank
    CALL MPI_Allgatherv(MPI_IN_PLACE, AMat%local_columns, MPI_INTEGER, &
         & col_root_lookup, cols_per_proc, d_cols_per_proc, MPI_INTEGER, &
         & process_grid%row_comm, ierr)

  END SUBROUTINE ConstructRankLookup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine which computes sparse dot products across processors.
  !> Computes the dot product of one vector with several others.
  SUBROUTINE DotAllHelper_r(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, out_values, comm)
    !> The length of vector i.
    INTEGER, INTENT(IN) :: num_values_i
    !> Tn array with the length of vectors j.
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    !> The index value of the sparse vector i.
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    !> The indices of the vectors j.
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    !> The values of the sparse vector i.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    !> The values of the vectors j.
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    !> The dot product values for each vector j.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    !> The communicator to reduce along.
    INTEGER, INTENT(INOUT) :: comm


    !! Local Variables
    INTEGER :: err
    INTEGER :: counter
    INTEGER :: inner_len_j

    !! Local Dot
    !$omp parallel private(inner_len_j)
    !$omp do
    DO counter = 1, SIZE(num_values_j)
       inner_len_j = num_values_j(counter)
       out_values(counter) = DotSparseVectors(indices_i(:num_values_i), &
            & values_i(:num_values_i), indices_j(:inner_len_j, counter), &
            & values_j(:inner_len_j, counter))
    END DO
    !$omp end do
    !$omp end parallel

    !! Reduce Over Processes
    CALL MPI_Allreduce(MPI_IN_PLACE, out_values, SIZE(num_values_j), &
         & MPINTREAL, MPI_SUM, comm, err)

  END SUBROUTINE DotAllHelper_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Helper routine which computes sparse dot products across processors.
  !> Computes the dot product of one vector with several others.
  !> The pivoted version has the number of local pivots to work on as a
  !> parameter.
  SUBROUTINE DotAllPivoted_r(num_values_i, indices_i, values_i, num_values_j, &
       & indices_j, values_j, pivot_vector, num_local_pivots, out_values, comm)
    !> The length of vector i.
    INTEGER, INTENT(IN) :: num_values_i
    !> Tn array with the length of vectors j.
    INTEGER, DIMENSION(:), INTENT(IN) :: num_values_j
    !> The index value of the sparse vector i.
    INTEGER, DIMENSION(:), INTENT(IN) :: indices_i
    !> The indices of the vectors j.
    INTEGER, DIMENSION(:,:), INTENT(IN) :: indices_j
    !> The values of the sparse vector i.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_i
    !> The values of the vectors j.
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values_j
    !> Vector storing the pivot values.
    INTEGER, DIMENSION(:), INTENT(IN) :: pivot_vector
    !> Number of pivots.
    INTEGER, INTENT(IN) :: num_local_pivots
    !> The dot product values for each vector j.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: out_values
    !> The communicator to reduce along.
    INTEGER, INTENT(INOUT) :: comm


    !! Local Variables
    INTEGER :: err
    INTEGER :: counter
    INTEGER :: inner_len_j
    INTEGER :: local_pi_i

    !! Local Dot
    !$omp parallel private(inner_len_j, local_pi_i)
    !$omp do
    DO counter = 1, num_local_pivots
       local_pi_i = pivot_vector(counter)
       inner_len_j = num_values_j(local_pi_i)
       out_values(counter) = DotSparseVectors(indices_i(:num_values_i), &
            & values_i(:num_values_i), indices_j(:inner_len_j, local_pi_i), &
            & values_j(:inner_len_j, local_pi_i))
    END DO
    !$omp end do
    !$omp end parallel

    !! Reduce Over Processes
    CALL MPI_Allreduce(MPI_IN_PLACE, out_values, SIZE(num_values_j), &
         & MPINTREAL, MPI_SUM, comm, err)

  END SUBROUTINE DotAllPivoted_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A helper routine that gathers the matrices in the same column into one.
  SUBROUTINE GatherMatrixColumn_r(local_matrix, column_matrix, process_grid)
    !> The local matrix on each process.
    TYPE(Matrix_lsr), INTENT(IN) :: local_matrix
    !> The final result.
    TYPE(Matrix_lsr), INTENT(INOUT) :: column_matrix
    !> The process grid to operate on.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !! Local Variables
    TYPE(Matrix_lsr) :: local_matrixT


    CALL TransposeMatrix(local_matrix, local_matrixT)
    CALL ReduceAndComposeMatrix(local_matrixT, column_matrix, &
         & process_grid%column_comm)

    CALL DestructMatrix(local_matrixT)
  END SUBROUTINE GatherMatrixColumn_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the pivot vector.
  SUBROUTINE GetPivot_r(AMat, process_grid, start_index, pivot_vector, diag, &
       & index, VALUE, local_pivots, num_local_pivots)
    !> The matrix we're working on.
    TYPE(Matrix_ps), INTENT(IN) :: AMat
    !> The process grid to compute on.
    TYPE(ProcessGrid_t), INTENT(INOUT) :: process_grid
    !> The current pivot vector.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: pivot_vector
    !> The diagonal values.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: diag
    !> The start index to look
    INTEGER, INTENT(IN) :: start_index
    !> The pivot index selected.
    INTEGER, INTENT(OUT) :: index
    !> The pivot value.
    REAL(NTREAL), INTENT(OUT) :: VALUE
    !> The local pivot values to modify.
    INTEGER, DIMENSION(:), INTENT(INOUT) :: local_pivots
    !> Number of pivots stored locally.
    INTEGER, INTENT(OUT) :: num_local_pivots
    !! Local Variables
    REAL(NTREAL) :: temp_diag
    DOUBLE PRECISION, DIMENSION(2) :: max_diag


    !! Local Variables
    INTEGER :: pind
    INTEGER :: II
    INTEGER :: swap
    INTEGER :: ierr

    !! Search for the maximum pivot
    max_diag = [0, 0]
    DO II = start_index, AMat%actual_matrix_dimension
       pind = pivot_vector(II)
       IF (pind .GE. AMat%start_column .AND. pind .LT. AMat%end_column) THEN
          temp_diag = diag(pind - AMat%start_column + 1)
          IF (temp_diag .GT. max_diag(1)) THEN
             max_diag(1) = temp_diag
             max_diag(2) = II
          END IF
       END IF
    END DO
    CALL MPI_Allreduce(MPI_IN_PLACE, max_diag, 1, MPI_2DOUBLE_PRECISION, &
         & MPI_MAXLOC, process_grid%row_comm, ierr)
    index = INT(max_diag(2))
    VALUE = max_diag(1)

    swap = pivot_vector(index)
    pivot_vector(index) = pivot_vector(start_index)
    pivot_vector(start_index) = swap

    index = swap

    !! Determine local pivots
    num_local_pivots = 0
    DO II = start_index + 1, AMat%actual_matrix_dimension
       pind = pivot_vector(II)
       IF (pind .GE. AMat%start_column .AND. pind .LT. AMat%end_column) THEN
          num_local_pivots = num_local_pivots + 1
          local_pivots(num_local_pivots) = pind - AMat%start_column + 1
       END IF
    END DO

  END SUBROUTINE GetPivot_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Unpack to a global matrix.
  SUBROUTINE UnpackCholesky_r(values_per_column, index, values, LMat)
    !> The number of values in a column.
    INTEGER, DIMENSION(:), INTENT(IN) :: values_per_column
    !> Index values.
    INTEGER, DIMENSION(:,:), INTENT(IN) :: index
    !> Actual values.
    REAL(NTREAL), DIMENSION(:,:), INTENT(IN) :: values
    !> Matrix to unpack into.
    TYPE(Matrix_ps), INTENT(INOUT) :: LMat


    !! Local Variables
    INTEGER :: local_columns
    TYPE(TripletList_r) :: local_triplets
    TYPE(Triplet_r) :: temp
    INTEGER :: II, JJ

    local_columns = LMat%local_columns

    local_triplets = TripletList_r()
    IF (LMat%process_grid%my_slice .EQ. 0) THEN
       DO JJ = 1, local_columns
          !! note transpose
          temp%index_row = JJ + LMat%start_column - 1
          DO II = 1, values_per_column(JJ)
             !! note transpose
             temp%index_column = INDEX(II,JJ) + LMat%start_row - 1
             temp%point_value = values(II,JJ)
             CALL AppendToTripletList(local_triplets, temp)
          END DO
       END DO
    END IF
    CALL FillMatrixFromTripletList(LMat, local_triplets)

    !! Cleanup
    CALL DestructTripletList(local_triplets)
  END SUBROUTINE UnpackCholesky_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE CholeskyModule
