

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module allows one to convert a sparse matrix to a dense matrix. It also
!! supports dense the dense versions of core matrix routines. This module is
!! used in situations where matrices become too dense for good sparse matrix
!! performance.
MODULE DMatrixModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc
  USE TripletListModule
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a dense matrix.
  TYPE, PUBLIC :: Matrix_ldr
     REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE :: DATA !< values of the matrix.
     INTEGER :: rows !< Matrix dimension: rows.
     INTEGER :: columns !< Matrix dimension: columns.
  END TYPE Matrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A datatype for storing a dense matrix.
  TYPE, PUBLIC :: Matrix_ldc
     !> values of the matrix.
     COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE :: DATA
     INTEGER :: rows !< Matrix dimension: rows.
     INTEGER :: columns !< Matrix dimension: columns.
  END TYPE Matrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructEmptyMatrix
  PUBLIC :: ConstructMatrixDFromS
  PUBLIC :: ConstructMatrixSFromD
  PUBLIC :: CopyMatrix
  PUBLIC :: DestructMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SplitMatrix
  PUBLIC :: ComposeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MatrixNorm
  PUBLIC :: IncrementMatrix
  PUBLIC :: MultiplyMatrix
  PUBLIC :: TransposeMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE Matrix_ldr
     MODULE PROCEDURE ConstructEmptyMatrix_ldr
  END INTERFACE
  INTERFACE Matrix_ldc
     MODULE PROCEDURE ConstructEmptyMatrix_ldc
  END INTERFACE
  INTERFACE ConstructEmptyMatrix
     MODULE PROCEDURE ConstructEmptyMatrixSup_ldr
     MODULE PROCEDURE ConstructEmptyMatrixSup_ldc
  END INTERFACE
  INTERFACE ConstructMatrixDFromS
     MODULE PROCEDURE ConstructMatrixDFromS_ldr
     MODULE PROCEDURE ConstructMatrixDFromS_ldc
  END INTERFACE
  INTERFACE ConstructMatrixSFromD
     MODULE PROCEDURE ConstructMatrixSFromD_ldr
     MODULE PROCEDURE ConstructMatrixSFromD_ldc
  END INTERFACE
  INTERFACE CopyMatrix
     MODULE PROCEDURE CopyMatrix_ldr
     MODULE PROCEDURE CopyMatrix_ldc
  END INTERFACE
  INTERFACE DestructMatrix
     MODULE PROCEDURE DestructMatrix_ldr
     MODULE PROCEDURE DestructMatrix_ldc
  END INTERFACE
  INTERFACE SplitMatrix
     MODULE PROCEDURE SplitMatrix_ldr
     MODULE PROCEDURE SplitMatrix_ldc
  END INTERFACE
  INTERFACE ComposeMatrix
     MODULE PROCEDURE ComposeMatrix_ldr
     MODULE PROCEDURE ComposeMatrix_ldc
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_ldr
     MODULE PROCEDURE MatrixNorm_ldc
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_ldr
     MODULE PROCEDURE IncrementMatrix_ldc
  END INTERFACE
  INTERFACE MultiplyMatrix
     MODULE PROCEDURE MultiplyMatrix_ldr
     MODULE PROCEDURE MultiplyMatrix_ldc
  END INTERFACE
  INTERFACE TransposeMatrix
     MODULE PROCEDURE TransposeMatrix_ldr
     MODULE PROCEDURE TransposeMatrix_ldc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine wrapper for the empty constructor.
  PURE SUBROUTINE ConstructEmptyMatrixSup_ldr(this, rows, columns)
    !> The matrix to construct
    TYPE(Matrix_ldr), INTENT(INOUT) :: this
    !> Rows of the matrix
    INTEGER, INTENT(IN) :: rows
    !> Columns of the matrix
    INTEGER, INTENT(IN) :: columns

    this = ConstructEmptyMatrix_ldr(rows, columns)
  END SUBROUTINE ConstructEmptyMatrixSup_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty dense matrix with a set number of rows and columns
  PURE FUNCTION ConstructEmptyMatrix_ldr(rows, columns) RESULT(this)
    !> The matrix to construct.
    TYPE(Matrix_ldr) :: this
    !> Rows of the matrix
    INTEGER, INTENT(IN) :: rows
    !> Columns of the matrix.
    INTEGER, INTENT(IN) :: columns

  this%rows = rows
  this%columns = columns

  ALLOCATE(this%data(rows,columns))

  END FUNCTION ConstructEmptyMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  PURE SUBROUTINE ConstructMatrixDFromS_ldr(sparse_matrix, dense_matrix)
    !> The sparse matrix to convert.
    TYPE(Matrix_lsr), INTENT(IN) :: sparse_matrix
    !> Output. Must be preallocated.
    TYPE(Matrix_ldr), INTENT(INOUT) :: dense_matrix
    !! Helper Variables
    TYPE(Triplet_r) :: temporary

  !! Helper Variables
  INTEGER :: inner_counter, outer_counter
  INTEGER :: elements_per_inner
  INTEGER :: total_counter

  CALL ConstructEmptyMatrix(dense_matrix, sparse_matrix%rows, &
       & sparse_matrix%columns)

  !! Loop over elements.
  dense_matrix%data = 0
  total_counter = 1
  DO outer_counter = 1, sparse_matrix%columns
     elements_per_inner = sparse_matrix%outer_index(outer_counter+1) - &
          & sparse_matrix%outer_index(outer_counter)
     temporary%index_column = outer_counter
     DO inner_counter = 1, elements_per_inner
        temporary%index_row = sparse_matrix%inner_index(total_counter)
        temporary%point_value = sparse_matrix%values(total_counter)
        dense_matrix%data(temporary%index_row, temporary%index_column) = &
             & temporary%point_value
        total_counter = total_counter + 1
     END DO
  END DO

  END SUBROUTINE ConstructMatrixDFromS_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  PURE SUBROUTINE ConstructMatrixSFromD_ldr(dense_matrix, sparse_matrix, &
       & threshold_in)
    !> Matrix to convert.
    TYPE(Matrix_ldr), INTENT(IN) :: dense_matrix
    !> Output matrix.
    TYPE(Matrix_lsr), INTENT(INOUT) :: sparse_matrix
    !> Value for pruning values to zero.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local Variables
    TYPE(Triplet_r) :: temporary
    TYPE(TripletList_r) :: temporary_list

  !! Local Variables
  INTEGER :: inner_counter, outer_counter
  INTEGER :: columns, rows

  columns = dense_matrix%columns
  rows = dense_matrix%rows

  IF (PRESENT(threshold_in)) THEN
     CALL ConstructTripletList(temporary_list)
     DO outer_counter = 1, columns
        temporary%index_column = outer_counter
        DO inner_counter = 1, rows
           temporary%point_value = &
                & dense_matrix%data(inner_counter,outer_counter)
           IF (ABS(temporary%point_value) .GT. threshold_in) THEN
              temporary%index_row = inner_counter
              CALL AppendToTripletList(temporary_list,temporary)
           END IF
        END DO
     END DO
  ELSE
     CALL ConstructTripletList(temporary_list, rows*columns)
     DO outer_counter = 1, columns
        temporary%index_column = outer_counter
        DO inner_counter = 1, rows
           temporary%point_value = &
                & dense_matrix%data(inner_counter,outer_counter)
           temporary%index_row = inner_counter
           temporary_list%data(inner_counter+rows*(outer_counter-1)) = &
                & temporary
        END DO
     END DO
  END IF

  sparse_matrix = Matrix_lsr(temporary_list, rows, columns)

  END SUBROUTINE ConstructMatrixSFromD_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy the matrix A into the B.
  PURE SUBROUTINE CopyMatrix_ldr(matA, matB)
    !> The matrix to copy.
    TYPE(Matrix_ldr), INTENT(IN) :: matA
    !> matB = matA
    TYPE(Matrix_ldr), INTENT(INOUT) :: matB

  CALL ConstructEmptyMatrix(matB, matA%rows, matA%columns)
  matB%data = matA%data

  END SUBROUTINE CopyMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate the memory associated with this matrix.
  PURE SUBROUTINE DestructMatrix_ldr(this)
    !> The matrix to delete.
    TYPE(Matrix_ldr), INTENT(INOUT) :: this

  IF (ALLOCATED(this%data)) THEN
     DEALLOCATE(this%data)
  END IF

  END SUBROUTINE DestructMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> AXPY for dense matrices. B = B + alpha*A
  PURE SUBROUTINE IncrementMatrix_ldr(MatA,MatB,alpha_in)
    !> MatA is added
    TYPE(Matrix_ldr), INTENT(IN) :: MatA
    !> MatB is incremented.
    TYPE(Matrix_ldr), INTENT(INOUT) :: MatB
    !> A scaling parameter.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !! Temporary
    REAL(NTREAL) :: alpha

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0d+0
  ELSE
     alpha = alpha_in
  END IF

  MatB%data = MatB%data + alpha*MatA%data

  END SUBROUTINE IncrementMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a dense matrix.
  !> Computes the Frobenius norm.
  FUNCTION MatrixNorm_ldr(this) RESULT(norm)
    !! Parameters
    !> The matrix to compute the norm of.
    TYPE(Matrix_ldr), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    INTEGER :: II, JJ

    norm = 0
    DO II =1, this%rows
       DO JJ = 1,  this%columns
          norm = norm + this%data(II,JJ)**2
       END DO
    END DO
  END FUNCTION MatrixNorm_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a dense matrix.
  PURE SUBROUTINE TransposeMatrix_ldr(matA, matAT)
    !> matA the matrix to transpose.
    TYPE(Matrix_ldr), INTENT(IN) :: matA
    !> matAT = matA^T.
    TYPE(Matrix_ldr), INTENT(INOUT) :: matAT

  CALL ConstructEmptyMatrix(matAT, matA%columns, matA%rows)
  matAT%data = TRANSPOSE(matA%data)

  END SUBROUTINE TransposeMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !> to another.
  PURE SUBROUTINE ComposeMatrix_ldr(mat_array, block_rows, block_columns, &
       & out_matrix)
    !> The number of rows of the array of blocks.
    INTEGER, INTENT(IN) :: block_rows
    !> The number of columns of the array of blocks.
    INTEGER, INTENT(IN) :: block_columns
    !> 2d array of matrices to compose.
    TYPE(Matrix_ldr), DIMENSION(block_rows, block_columns), INTENT(IN) :: &
         & mat_array
    !> The composed matrix.
    TYPE(Matrix_ldr), INTENT(INOUT) :: out_matrix

  !! Local Data
  INTEGER, DIMENSION(block_rows+1) :: row_offsets
  INTEGER, DIMENSION(block_columns+1) :: column_offsets
  INTEGER :: out_rows, out_columns
  INTEGER :: II, JJ

  !! Determine the size of the big matrix
  out_columns = 0
  column_offsets(1) = 1
  out_columns = 0
  DO JJ = 1, block_columns
     column_offsets(JJ+1) = column_offsets(JJ) + mat_array(1,JJ)%columns
     out_columns = out_columns + mat_array(1,JJ)%columns
  END DO
  row_offsets(1) = 1
  out_rows = 0
  DO II = 1, block_rows
     row_offsets(II+1) = row_offsets(II) + mat_array(II,1)%rows
     out_rows = out_rows + mat_array(II,1)%rows
  END DO

  !! Allocate Memory
  CALL ConstructEmptyMatrix(out_matrix, out_columns, out_rows)

  DO JJ = 1, block_columns
     DO II = 1, block_rows
        out_matrix%data(row_offsets(II):row_offsets(II+1)-1, &
             & column_offsets(JJ):column_offsets(JJ+1)-1) = &
             & mat_array(II,JJ)%data
     END DO
  END DO

  END SUBROUTINE ComposeMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  PURE SUBROUTINE SplitMatrix_ldr(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !> The matrix to split.
    TYPE(Matrix_ldr), INTENT(IN) :: this
    !> Number of rows to split the matrix into.
    INTEGER, INTENT(IN) :: block_rows
    !> Number of columns to split the matrix into.
    INTEGER, INTENT(IN) :: block_columns
    !> A block_columns x block_rows array for the output to go into.
    TYPE(Matrix_ldr), DIMENSION(:,:), INTENT(INOUT) :: split_array
    !> Specifies the size of the  rows.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    !> Specifies the size of the columns.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in

  !! Local Data
  INTEGER, DIMENSION(block_rows) :: block_size_row
  INTEGER, DIMENSION(block_columns) :: block_size_column
  INTEGER, DIMENSION(block_rows+1) :: row_offsets
  INTEGER, DIMENSION(block_columns+1) :: column_offsets
  !! Temporary Variables
  INTEGER :: divisor_row, divisor_column
  INTEGER :: II, JJ

  !! Calculate the split sizes
  IF (PRESENT(block_size_row_in)) THEN
     block_size_row = block_size_row_in
  ELSE
     divisor_row = this%rows/block_rows
     block_size_row = divisor_row
     block_size_row(block_rows) = this%rows - divisor_row*(block_rows-1)
  END IF
  IF (PRESENT(block_size_column_in)) THEN
     block_size_column = block_size_column_in
  ELSE
     divisor_column = this%columns/block_columns
     block_size_column = divisor_column
     block_size_column(block_columns) = this%columns - &
          & divisor_column*(block_columns-1)
  END IF

  !! Copy the block offsets
  row_offsets(1) = 1
  DO II = 1, block_rows
     row_offsets(II+1) = row_offsets(II) + block_size_row(II)
  END DO
  column_offsets(1) = 1
  DO JJ = 1, block_columns
     column_offsets(JJ+1) = column_offsets(JJ) + block_size_column(JJ)
  END DO

  !! Copy
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        CALL ConstructEmptyMatrix(split_array(II,JJ), block_size_column(JJ), &
             & block_size_row(II))
        split_array(II,JJ)%data = &
             & this%data(row_offsets(II):row_offsets(II+1)-1, &
             & column_offsets(JJ):column_offsets(JJ+1)-1)
     END DO
  END DO

  END SUBROUTINE SplitMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  SUBROUTINE MultiplyMatrix_ldr(MatA,MatB,MatC)
    !> The first matrix.
    TYPE(Matrix_ldr), INTENT(IN) :: MatA
    !> The second matrix.
    TYPE(Matrix_ldr), INTENT(IN) :: MatB
    !> MatC = MatA*MatB.
    TYPE(Matrix_ldr), INTENT(INOUT) :: MatC
    !! Local variables
    CHARACTER, PARAMETER :: TRANSA = 'N'
    CHARACTER, PARAMETER :: TRANSB = 'N'
    INTEGER :: M
    INTEGER :: N
    INTEGER :: K
    DOUBLE PRECISION, PARAMETER :: ALPHA = 1.0
    INTEGER :: LDA
    INTEGER :: LDB
    DOUBLE PRECISION, PARAMETER :: BETA = 0.0
    INTEGER :: LDC

    MatC = Matrix_ldr(MatA%rows,MatB%columns)

    !! Setup Lapack
    M = MatA%rows
    N = MatB%columns
    K = MatA%columns
    LDA = M
    LDB = K
    LDC = M

    CALL DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%data, LDA, MatB%data, &
         & LDB, BETA, MatC%data, LDC)

  END SUBROUTINE MultiplyMatrix_ldr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A subroutine style wrapper for the constructor.
  PURE SUBROUTINE ConstructEmptyMatrixSup_ldc(this, rows, columns)
    !> The matrix to construct.
    TYPE(Matrix_ldc), INTENT(INOUT) :: this
    !> The number of rows of the matrix.
    INTEGER, INTENT(IN) :: rows
    !> The number of columns o the matrix.
    INTEGER, INTENT(IN) :: columns

    this = ConstructEmptyMatrix_ldc(rows, columns)
  END SUBROUTINE ConstructEmptyMatrixSup_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct an empty dense matrix with a set number of rows and columns
  PURE FUNCTION ConstructEmptyMatrix_ldc(rows, columns) RESULT(this)
    !> The matrix to construct.
    TYPE(Matrix_ldc) :: this
    !> Rows of the matrix
    INTEGER, INTENT(IN) :: rows
    !> Columns of the matrix.
    INTEGER, INTENT(IN) :: columns

  this%rows = rows
  this%columns = columns

  ALLOCATE(this%data(rows,columns))

  END FUNCTION ConstructEmptyMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a sparse matrix to a dense matrix.
  PURE SUBROUTINE ConstructMatrixDFromS_ldc(sparse_matrix, dense_matrix)
    !> The sparse matrix to convert.
    TYPE(Matrix_lsc), INTENT(IN) :: sparse_matrix
    !> Dense matrix output. Must be preallocated.
    TYPE(Matrix_ldc), INTENT(INOUT) :: dense_matrix
    !! Helper Variables
    TYPE(Triplet_c) :: temporary

  !! Helper Variables
  INTEGER :: inner_counter, outer_counter
  INTEGER :: elements_per_inner
  INTEGER :: total_counter

  CALL ConstructEmptyMatrix(dense_matrix, sparse_matrix%rows, &
       & sparse_matrix%columns)

  !! Loop over elements.
  dense_matrix%data = 0
  total_counter = 1
  DO outer_counter = 1, sparse_matrix%columns
     elements_per_inner = sparse_matrix%outer_index(outer_counter+1) - &
          & sparse_matrix%outer_index(outer_counter)
     temporary%index_column = outer_counter
     DO inner_counter = 1, elements_per_inner
        temporary%index_row = sparse_matrix%inner_index(total_counter)
        temporary%point_value = sparse_matrix%values(total_counter)
        dense_matrix%data(temporary%index_row, temporary%index_column) = &
             & temporary%point_value
        total_counter = total_counter + 1
     END DO
  END DO

  END SUBROUTINE ConstructMatrixDFromS_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A function that converts a dense matrix to a sparse matrix.
  PURE SUBROUTINE ConstructMatrixSFromD_ldc(dense_matrix, sparse_matrix, &
       & threshold_in)
    !> The matrix to convert.
    TYPE(Matrix_ldc), INTENT(IN) :: dense_matrix
    !> The sparse output matrix.
    TYPE(Matrix_lsc), INTENT(INOUT) :: sparse_matrix
    !> Value for pruning values to zero.
    REAL(NTREAL), INTENT(IN), OPTIONAL :: threshold_in
    !! Local Variables
    TYPE(Triplet_c) :: temporary
    TYPE(TripletList_c) :: temporary_list

  !! Local Variables
  INTEGER :: inner_counter, outer_counter
  INTEGER :: columns, rows

  columns = dense_matrix%columns
  rows = dense_matrix%rows

  IF (PRESENT(threshold_in)) THEN
     CALL ConstructTripletList(temporary_list)
     DO outer_counter = 1, columns
        temporary%index_column = outer_counter
        DO inner_counter = 1, rows
           temporary%point_value = &
                & dense_matrix%data(inner_counter,outer_counter)
           IF (ABS(temporary%point_value) .GT. threshold_in) THEN
              temporary%index_row = inner_counter
              CALL AppendToTripletList(temporary_list,temporary)
           END IF
        END DO
     END DO
  ELSE
     CALL ConstructTripletList(temporary_list, rows*columns)
     DO outer_counter = 1, columns
        temporary%index_column = outer_counter
        DO inner_counter = 1, rows
           temporary%point_value = &
                & dense_matrix%data(inner_counter,outer_counter)
           temporary%index_row = inner_counter
           temporary_list%data(inner_counter+rows*(outer_counter-1)) = &
                & temporary
        END DO
     END DO
  END IF

  sparse_matrix = Matrix_lsc(temporary_list, rows, columns)

  END SUBROUTINE ConstructMatrixSFromD_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Copy the matrix A into the B.
  PURE SUBROUTINE CopyMatrix_ldc(matA, matB)
    !> The matrix to copy.
    TYPE(Matrix_ldc), INTENT(IN) :: matA
    !> matB = matA
    TYPE(Matrix_ldc), INTENT(INOUT) :: matB

  CALL ConstructEmptyMatrix(matB, matA%rows, matA%columns)
  matB%data = matA%data

  END SUBROUTINE CopyMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Deallocate the memory associated with this matrix.
  PURE SUBROUTINE DestructMatrix_ldc(this)
    !> This the matrix to delete.
    TYPE(Matrix_ldc), INTENT(INOUT) :: this

  IF (ALLOCATED(this%data)) THEN
     DEALLOCATE(this%data)
  END IF

  END SUBROUTINE DestructMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> AXPY for dense matrices. B = B + alpha*A
  PURE SUBROUTINE IncrementMatrix_ldc(MatA,MatB,alpha_in)
    !> MatA is added
    TYPE(Matrix_ldc), INTENT(IN) :: MatA
    !> MatB is incremented.
    TYPE(Matrix_ldc), INTENT(INOUT) :: MatB
    !> A scaling parameter.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !! Temporary
    REAL(NTREAL) :: alpha

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0d+0
  ELSE
     alpha = alpha_in
  END IF

  MatB%data = MatB%data + alpha*MatA%data

  END SUBROUTINE IncrementMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a dense matrix.
  !> Computes the Frobenius norm.
  FUNCTION MatrixNorm_ldc(this) RESULT(norm)
    !> The matrix to compute the norm of.
    TYPE(Matrix_ldc), INTENT(IN) :: this
    !> The norm of the matrix.
    REAL(NTREAL) :: norm
    !! Local Variables
    INTEGER :: II, JJ

    norm = 0
    DO II =1, this%rows
       DO JJ = 1,  this%columns
          norm = norm + &
               & REAL(this%data(II,JJ)*CONJG(this%data(II,JJ)),KIND=NTREAL)
       END DO
    END DO
  END FUNCTION MatrixNorm_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Transpose a dense matrix.
  PURE SUBROUTINE TransposeMatrix_ldc(matA, matAT)
    !> The matrix to transpose.
    TYPE(Matrix_ldc), INTENT(IN) :: matA
    !> matAT = matA^T.
    TYPE(Matrix_ldc), INTENT(INOUT) :: matAT

  CALL ConstructEmptyMatrix(matAT, matA%columns, matA%rows)
  matAT%data = TRANSPOSE(matA%data)

  END SUBROUTINE TransposeMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Create a big matrix from an array of matrices by putting them one next
  !> to another.
  PURE SUBROUTINE ComposeMatrix_ldc(mat_array, block_rows, block_columns, &
       & out_matrix)
    !> The number of rows of the array of blocks.
    INTEGER, INTENT(IN) :: block_rows
    !> The number of columns of the array of blocks.
    INTEGER, INTENT(IN) :: block_columns
    !> 2d array of matrices to compose.
    TYPE(Matrix_ldc), DIMENSION(block_rows, block_columns), INTENT(IN) :: &
         & mat_array
    !> The composed matrix.
    TYPE(Matrix_ldc), INTENT(INOUT) :: out_matrix

  !! Local Data
  INTEGER, DIMENSION(block_rows+1) :: row_offsets
  INTEGER, DIMENSION(block_columns+1) :: column_offsets
  INTEGER :: out_rows, out_columns
  INTEGER :: II, JJ

  !! Determine the size of the big matrix
  out_columns = 0
  column_offsets(1) = 1
  out_columns = 0
  DO JJ = 1, block_columns
     column_offsets(JJ+1) = column_offsets(JJ) + mat_array(1,JJ)%columns
     out_columns = out_columns + mat_array(1,JJ)%columns
  END DO
  row_offsets(1) = 1
  out_rows = 0
  DO II = 1, block_rows
     row_offsets(II+1) = row_offsets(II) + mat_array(II,1)%rows
     out_rows = out_rows + mat_array(II,1)%rows
  END DO

  !! Allocate Memory
  CALL ConstructEmptyMatrix(out_matrix, out_columns, out_rows)

  DO JJ = 1, block_columns
     DO II = 1, block_rows
        out_matrix%data(row_offsets(II):row_offsets(II+1)-1, &
             & column_offsets(JJ):column_offsets(JJ+1)-1) = &
             & mat_array(II,JJ)%data
     END DO
  END DO

  END SUBROUTINE ComposeMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Split a sparse matrix into an array of sparse matrices.
  PURE SUBROUTINE SplitMatrix_ldc(this, block_rows, block_columns, &
       & split_array, block_size_row_in, block_size_column_in)
    !> The matrix to split.
    TYPE(Matrix_ldc), INTENT(IN) :: this
    !> Number of rows to split the matrix into.
    INTEGER, INTENT(IN) :: block_rows
    !> Number of columns to split the matrix into.
    INTEGER, INTENT(IN) :: block_columns
    !> A COLUMNxROW array for the output to go into.
    TYPE(Matrix_ldc), DIMENSION(:,:), INTENT(INOUT) :: split_array
    !> Specifies the size of the  rows.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_row_in
    !> Specifies the size of the columns.
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: block_size_column_in

  !! Local Data
  INTEGER, DIMENSION(block_rows) :: block_size_row
  INTEGER, DIMENSION(block_columns) :: block_size_column
  INTEGER, DIMENSION(block_rows+1) :: row_offsets
  INTEGER, DIMENSION(block_columns+1) :: column_offsets
  !! Temporary Variables
  INTEGER :: divisor_row, divisor_column
  INTEGER :: II, JJ

  !! Calculate the split sizes
  IF (PRESENT(block_size_row_in)) THEN
     block_size_row = block_size_row_in
  ELSE
     divisor_row = this%rows/block_rows
     block_size_row = divisor_row
     block_size_row(block_rows) = this%rows - divisor_row*(block_rows-1)
  END IF
  IF (PRESENT(block_size_column_in)) THEN
     block_size_column = block_size_column_in
  ELSE
     divisor_column = this%columns/block_columns
     block_size_column = divisor_column
     block_size_column(block_columns) = this%columns - &
          & divisor_column*(block_columns-1)
  END IF

  !! Copy the block offsets
  row_offsets(1) = 1
  DO II = 1, block_rows
     row_offsets(II+1) = row_offsets(II) + block_size_row(II)
  END DO
  column_offsets(1) = 1
  DO JJ = 1, block_columns
     column_offsets(JJ+1) = column_offsets(JJ) + block_size_column(JJ)
  END DO

  !! Copy
  DO JJ = 1, block_columns
     DO II = 1, block_rows
        CALL ConstructEmptyMatrix(split_array(II,JJ), block_size_column(JJ), &
             & block_size_row(II))
        split_array(II,JJ)%data = &
             & this%data(row_offsets(II):row_offsets(II+1)-1, &
             & column_offsets(JJ):column_offsets(JJ+1)-1)
     END DO
  END DO

  END SUBROUTINE SplitMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A wrapper for multiplying two dense matrices.
  SUBROUTINE MultiplyMatrix_ldc(MatA,MatB,MatC)
    !> The first matrix.
    TYPE(Matrix_ldc), INTENT(IN) :: MatA
    !> The second matrix.
    TYPE(Matrix_ldc), INTENT(IN) :: MatB
    !> MatC = MatA*MatB.
    TYPE(Matrix_ldc), INTENT(INOUT) :: MatC
    !! Local variables
    CHARACTER, PARAMETER :: TRANSA = 'N'
    CHARACTER, PARAMETER :: TRANSB = 'N'
    INTEGER :: M
    INTEGER :: N
    INTEGER :: K
    COMPLEX*16, PARAMETER :: ALPHA = 1.0
    INTEGER :: LDA
    INTEGER :: LDB
    COMPLEX*16, PARAMETER :: BETA = 0.0
    INTEGER :: LDC

    MatC = Matrix_ldc(MatA%rows,MatB%columns)

    !! Setup Lapack
    M = MatA%rows
    N = MatB%columns
    K = MatA%columns
    LDA = M
    LDB = K
    LDC = M

    CALL ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, MatA%data, LDA, MatB%data, &
         & LDB, BETA, MatC%data, LDC)

  END SUBROUTINE MultiplyMatrix_ldc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE DMatrixModule
