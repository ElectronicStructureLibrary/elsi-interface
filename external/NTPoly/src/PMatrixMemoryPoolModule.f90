

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for distributed matrix multiplication.
MODULE PMatrixMemoryPoolModule
  USE PSMatrixModule, ONLY : Matrix_ps
  USE MatrixMemoryPoolModule, ONLY : MatrixMemoryPool_lr, MatrixMemoryPool_lc, &
       & DestructMatrixMemoryPool
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !> this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: MatrixMemoryPool_p
     !> Grid of local pools.
     TYPE(MatrixMemoryPool_lr), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid_r
     !> Grid of local pools (complex).
     TYPE(MatrixMemoryPool_lc), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: grid_c
  END TYPE MatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: DestructMatrixMemoryPool
  PUBLIC :: CheckMemoryPoolValidity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixMemoryPool_p
     MODULE PROCEDURE ConstructMatrixMemoryPool_p
  END INTERFACE
  INTERFACE DestructMatrixMemoryPool
     MODULE PROCEDURE DestructMatrixMemoryPool_p
  END INTERFACE
  INTERFACE CheckMemoryPoolValidity
     MODULE PROCEDURE CheckMemoryPoolValidity_p
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Distributed Matrix Memory Pool object.
  PURE FUNCTION ConstructMatrixMemoryPool_p(matrix) RESULT(this)
    !> A constructed Matrix Memory Pool object.
    TYPE(MatrixMemoryPool_p) :: this
    !> The associated distributed sparse matrix.
    TYPE(Matrix_ps), INTENT(IN) :: matrix

    !! Allocate
    IF (matrix%is_complex) THEN
       ALLOCATE(this%grid_c(matrix%process_grid%number_of_blocks_rows, &
            & matrix%process_grid%number_of_blocks_columns))
    ELSE
       ALLOCATE(this%grid_r(matrix%process_grid%number_of_blocks_rows, &
            & matrix%process_grid%number_of_blocks_columns))
    END IF
  END FUNCTION ConstructMatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destruct a Distributed Matrix Memory Pool object.
  PURE SUBROUTINE DestructMatrixMemoryPool_p(this)
    !> Distributed Matrix Memory Pool object to destroy.
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: this
    !! Local Data
    INTEGER :: row_counter, column_counter

  !! Allocate
  IF (ALLOCATED(this%grid_r)) THEN
     DO column_counter = LBOUND(this%grid_r,2), UBOUND(this%grid_r,2)
        DO row_counter = LBOUND(this%grid_r,1), UBOUND(this%grid_r,1)
           CALL DestructMatrixMemoryPool( &
                & this%grid_r(row_counter, column_counter))
        END DO
     END DO
     DEALLOCATE(this%grid_r)
  END IF

  !! Allocate
  IF (ALLOCATED(this%grid_c)) THEN
     DO column_counter = LBOUND(this%grid_c,2), UBOUND(this%grid_c,2)
        DO row_counter = LBOUND(this%grid_c,1), UBOUND(this%grid_c,1)
           CALL DestructMatrixMemoryPool( &
                & this%grid_c(row_counter, column_counter))
        END DO
     END DO
     DEALLOCATE(this%grid_c)
  END IF

  END SUBROUTINE DestructMatrixMemoryPool_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given distributed memory pool has been validly allocated to
  !> handle the given parameters.
  PURE FUNCTION CheckMemoryPoolValidity_p(this, matrix) RESULT(isvalid)
    !> The memory pool to check.
    TYPE(MatrixMemoryPool_p), INTENT(IN) :: this
    !> The associated matrix to check against.
    TYPE(Matrix_ps), INTENT(IN) :: matrix
    !> True if the memory pool is valid.
    LOGICAL :: isvalid

    isvalid = .TRUE.
    !! Check allocation
    IF (matrix%is_complex) THEN
       IF (.NOT. ALLOCATED(this%grid_c)) isvalid = .FALSE.
    ELSE
       IF (.NOT. ALLOCATED(this%grid_r)) isvalid = .FALSE.
    END IF

  END FUNCTION CheckMemoryPoolValidity_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PMatrixMemoryPoolModule
