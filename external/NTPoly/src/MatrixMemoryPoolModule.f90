

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling scratch memory for matrix multiplication.
!> The purpose of this module is to avoid having to allocate memory on the
!> heap during a matrix multiply, and to manage the underlying hash table.
MODULE MatrixMemoryPoolModule
  USE DataTypesModule, ONLY: NTREAL, NTCOMPLEX
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : c_int
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !> this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: MatrixMemoryPool_lr
     PRIVATE
     !> Shape of matrix: columns
     INTEGER, PUBLIC :: columns
     !> Shape of matrix: rows
     INTEGER, PUBLIC :: rows
     !> storage for actual values added to the matrix.
     TYPE(Triplet_r), DIMENSION(:), ALLOCATABLE, PUBLIC :: pruned_list
     !> storage for potential values added to the matrix.
     REAL(NTREAL), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: value_array
     !> true if an element has been pushed to this part of the matrix.
     LOGICAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dirty_array
     !> Storage space for indices, hashed.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: hash_index
     !> Internal storage space for amount of items added to a bucket.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: inserted_per_bucket
     !> Size of the buckets.
     INTEGER, PUBLIC :: hash_size
  END TYPE MatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A memory pool datatype that can be reused for matrix matrix multiplication.
  !> this is to prevent excessive alloc/dealloc.
  TYPE, PUBLIC :: MatrixMemoryPool_lc
     PRIVATE
     !> Shape of matrix: columns
     INTEGER, PUBLIC :: columns
     !> Shape of matrix: rows
     INTEGER, PUBLIC :: rows
     !> storage for actual values added to the matrix.
     TYPE(Triplet_c), DIMENSION(:), ALLOCATABLE, PUBLIC :: pruned_list
     !> storage for potential values added to the matrix.
     COMPLEX(NTCOMPLEX), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: value_array
     !> true if an element has been pushed to this part of the matrix.
     LOGICAL, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: dirty_array
     !> Storage space for indices, hashed.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: hash_index
     !> Internal storage space for amount of items added to a bucket.
     INTEGER, DIMENSION(:,:), ALLOCATABLE, PUBLIC :: inserted_per_bucket
     !> Size of the buckets.
     INTEGER, PUBLIC :: hash_size
  END TYPE MatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ConstructMatrixMemoryPool
  PUBLIC :: DestructMatrixMemoryPool
  PUBLIC :: CheckMemoryPoolValidity
  PUBLIC :: SetPoolSparsity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixMemoryPool_lr
     MODULE PROCEDURE ConstructMatrixMemoryPool_lr
  END INTERFACE
  INTERFACE MatrixMemoryPool_lc
     MODULE PROCEDURE ConstructMatrixMemoryPool_lc
  END INTERFACE
  INTERFACE ConstructMatrixMemoryPool
     MODULE PROCEDURE ConstructMatrixMemoryPoolSub_lr
     MODULE PROCEDURE ConstructMatrixMemoryPoolSub_lc
  END INTERFACE
  INTERFACE DestructMatrixMemoryPool
     MODULE PROCEDURE DestructMatrixMemoryPool_lr
     MODULE PROCEDURE DestructMatrixMemoryPool_lc
  END INTERFACE
  INTERFACE CheckMemoryPoolValidity
     MODULE PROCEDURE CheckMemoryPoolValidity_lr
     MODULE PROCEDURE CheckMemoryPoolValidity_lc
  END INTERFACE
  INTERFACE SetPoolSparsity
     MODULE PROCEDURE SetPoolSparsity_lr
     MODULE PROCEDURE SetPoolSparsity_lc
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for the constructor.
  SUBROUTINE ConstructMatrixMemoryPoolSub_lr(this, columns, rows, sparsity_in)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lr), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

    IF (PRESENT(sparsity_in)) THEN
       this = MatrixMemoryPool_lr(columns, rows, sparsity_in)
    ELSE
       this = MatrixMemoryPool_lr(columns, rows)
    END IF
  END SUBROUTINE ConstructMatrixMemoryPoolSub_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Subroutine wrapper for the constructor.
  SUBROUTINE ConstructMatrixMemoryPoolSub_lc(this, columns, rows, sparsity_in)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lc), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

    IF (PRESENT(sparsity_in)) THEN
       this = MatrixMemoryPool_lc(columns, rows, sparsity_in)
    ELSE
       this = MatrixMemoryPool_lc(columns, rows)
    END IF
  END SUBROUTINE ConstructMatrixMemoryPoolSub_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  FUNCTION ConstructMatrixMemoryPool_lr(columns, rows, sparsity_in) RESULT(this)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lr), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

  !! Temporary variables
  INTEGER :: alloc_stat
  INTEGER :: num_buckets

  this%columns = columns
  this%rows = rows

  IF (.NOT. PRESENT(sparsity_in)) THEN
     this%hash_size = 1
  ELSE
     this%hash_size = INT(1.0/sparsity_in)
     IF (this%hash_size > columns) this%hash_size = columns
  END IF

  num_buckets = columns/this%hash_size + 1

  !! Allocate
  ALLOCATE(this%pruned_list(columns*rows), stat=alloc_stat)
  ALLOCATE(this%value_array(columns,rows), stat=alloc_stat)
  ALLOCATE(this%dirty_array(columns,rows), stat=alloc_stat)

  ALLOCATE(this%hash_index(columns,rows))
  ALLOCATE(this%inserted_per_bucket(columns,rows))

  this%value_array = 0
  this%hash_index = 0
  this%inserted_per_bucket = 0
  this%dirty_array = .FALSE.

  END FUNCTION ConstructMatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Construct Matrix Memory Pool object.
  FUNCTION ConstructMatrixMemoryPool_lc(columns, rows, sparsity_in) RESULT(this)
    !> The matrix to construct.
    TYPE(MatrixMemoryPool_lc), TARGET :: this
    !> Number of columns in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER(kind=c_int), INTENT(IN) :: rows
    !> Estimated sparsity (optional).
    REAL(NTREAL), INTENT(IN), OPTIONAL :: sparsity_in

  !! Temporary variables
  INTEGER :: alloc_stat
  INTEGER :: num_buckets

  this%columns = columns
  this%rows = rows

  IF (.NOT. PRESENT(sparsity_in)) THEN
     this%hash_size = 1
  ELSE
     this%hash_size = INT(1.0/sparsity_in)
     IF (this%hash_size > columns) this%hash_size = columns
  END IF

  num_buckets = columns/this%hash_size + 1

  !! Allocate
  ALLOCATE(this%pruned_list(columns*rows), stat=alloc_stat)
  ALLOCATE(this%value_array(columns,rows), stat=alloc_stat)
  ALLOCATE(this%dirty_array(columns,rows), stat=alloc_stat)

  ALLOCATE(this%hash_index(columns,rows))
  ALLOCATE(this%inserted_per_bucket(columns,rows))

  this%value_array = 0
  this%hash_index = 0
  this%inserted_per_bucket = 0
  this%dirty_array = .FALSE.

  END FUNCTION ConstructMatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  PURE SUBROUTINE DestructMatrixMemoryPool_lr(this)
    !> The matrix being destructed.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT) :: this

  !! Perform deallocations.
  IF (ALLOCATED(this%pruned_list)) DEALLOCATE(this%pruned_list)
  IF (ALLOCATED(this%value_array)) DEALLOCATE(this%value_array)
  IF (ALLOCATED(this%dirty_array)) DEALLOCATE(this%dirty_array)
  IF (ALLOCATED(this%hash_index)) DEALLOCATE(this%hash_index)
  IF (ALLOCATED(this%inserted_per_bucket)) DEALLOCATE(this%inserted_per_bucket)

  END SUBROUTINE DestructMatrixMemoryPool_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> A destructor for a matrix memory pool
  PURE SUBROUTINE DestructMatrixMemoryPool_lc(this)
    !> The matrix being destructed.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT) :: this

  !! Perform deallocations.
  IF (ALLOCATED(this%pruned_list)) DEALLOCATE(this%pruned_list)
  IF (ALLOCATED(this%value_array)) DEALLOCATE(this%value_array)
  IF (ALLOCATED(this%dirty_array)) DEALLOCATE(this%dirty_array)
  IF (ALLOCATED(this%hash_index)) DEALLOCATE(this%hash_index)
  IF (ALLOCATED(this%inserted_per_bucket)) DEALLOCATE(this%inserted_per_bucket)

  END SUBROUTINE DestructMatrixMemoryPool_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !> the given parameters.
  PURE FUNCTION CheckMemoryPoolValidity_lr(this, columns, rows) RESULT(isvalid)
    !> The memory pool to check.
    TYPE(MatrixMemoryPool_lr), INTENT(in) :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> true if the memory pool is valid.
    LOGICAL :: isvalid

  isvalid = .TRUE.
  !! Check allocation
  IF (.NOT. ALLOCATED(this%pruned_list)) isvalid = .FALSE.
  IF (.NOT. ALLOCATED(this%value_array)) isvalid = .FALSE.

  !! Check allocation size
  IF (.NOT. SIZE(this%value_array,dim=2) .EQ. rows) THEN
     isvalid = .FALSE.
  END IF
  IF (.NOT. SIZE(this%value_array,dim=1) .EQ. columns) THEN
     isvalid = .FALSE.
  END IF

  END FUNCTION CheckMemoryPoolValidity_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Checks if a given memory pool has been validly allocated to handle
  !> Checks if a given memory pool has been validly allocated to handle
  !> the given parameters.
  PURE FUNCTION CheckMemoryPoolValidity_lc(this, columns, rows) RESULT(isvalid)
    !> The memory pool to check.
    TYPE(MatrixMemoryPool_lc), INTENT(in) :: this
    !> Number of columns in the matrix.
    INTEGER, INTENT(IN) :: columns
    !> Number of rows in the matrix.
    INTEGER, INTENT(IN) :: rows
    !> true if the memory pool is valid.
    LOGICAL :: isvalid

  isvalid = .TRUE.
  !! Check allocation
  IF (.NOT. ALLOCATED(this%pruned_list)) isvalid = .FALSE.
  IF (.NOT. ALLOCATED(this%value_array)) isvalid = .FALSE.

  !! Check allocation size
  IF (.NOT. SIZE(this%value_array,dim=2) .EQ. rows) THEN
     isvalid = .FALSE.
  END IF
  IF (.NOT. SIZE(this%value_array,dim=1) .EQ. columns) THEN
     isvalid = .FALSE.
  END IF

  END FUNCTION CheckMemoryPoolValidity_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  SUBROUTINE SetPoolSparsity_lr(this,sparsity)
    !> The memory pool to set the sparsity of.
    TYPE(MatrixMemoryPool_lr), INTENT(INOUT), TARGET :: this
    !> The sparsity value.
    REAL(NTREAL), INTENT(IN) :: sparsity

  !! Local Variables
  INTEGER :: num_buckets

  this%hash_size = INT(1.0/sparsity)
  IF (this%hash_size > this%columns) this%hash_size = this%columns
  num_buckets = this%columns/this%hash_size + 1

  END SUBROUTINE SetPoolSparsity_lr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets the expected sparsity of the matrix, which helps with hashing.
  SUBROUTINE SetPoolSparsity_lc(this,sparsity)
    !> The memory pool to set the sparsity of.
    TYPE(MatrixMemoryPool_lc), INTENT(INOUT), TARGET :: this
    !> The sparsity value.
    REAL(NTREAL), INTENT(IN) :: sparsity

  !! Local Variables
  INTEGER :: num_buckets

  this%hash_size = INT(1.0/sparsity)
  IF (this%hash_size > this%columns) this%hash_size = this%columns
  num_buckets = this%columns/this%hash_size + 1

  END SUBROUTINE SetPoolSparsity_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMemoryPoolModule
