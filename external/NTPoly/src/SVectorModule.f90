

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for handling compressed vectors.
!> Compressed vectors are stored in two lists. The first is a list of indices,
!> the second a list of values.
MODULE SVectorModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: AddSparseVectors
  PUBLIC :: DotSparseVectors
  PUBLIC :: PairwiseMultiplyVectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE AddSparseVectors
     MODULE PROCEDURE AddSparseVectors_r
     MODULE PROCEDURE AddSparseVectors_c
  END INTERFACE
  INTERFACE DotSparseVectors
     MODULE PROCEDURE DotSparseVectors_r
     MODULE PROCEDURE DotSparseVectors_c
  END INTERFACE
  INTERFACE PairwiseMultiplyVectors
     MODULE PROCEDURE PairwiseMultiplyVectors_r
     MODULE PROCEDURE PairwiseMultiplyVectors_c
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add together two sparse vectors. C = A + alpha*B
  !> The values that are returned for C are only valid in the range
  !> (1:total_values_c). We do not do an automatic shrinking of the array
  !> to keep this routine low in overhead.
  PURE SUBROUTINE AddSparseVectors_r(inner_index_a,values_a,inner_index_b, &
       & values_b,inner_index_c,values_c,total_values_c, alpha_in, threshold_in)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_b
    !> List of values computed for C.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: values_c
    !> Value to scale VecB by. Optional, default is 1.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> for flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> The total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Local Variables
    REAL(NTREAL) :: working_value_a, working_value_b

  !! Local Data
  REAL(NTREAL) :: alpha
  REAL(NTREAL) :: threshold
  !! Temporary Variables
  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b, counter_c

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0d+0
  ELSE
     alpha = alpha_in
  END IF
  IF (.NOT. PRESENT(threshold_in)) THEN
     threshold = 0.0d+0
  ELSE
     threshold = threshold_in
  END IF

  counter_a = 1
  counter_b = 1
  counter_c = 1
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = alpha*values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        IF (ABS(working_value_a + working_value_b) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_a
           values_c(counter_c) = working_value_a + working_value_b
           counter_c = counter_c + 1
        END IF
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        IF (ABS(working_value_b) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_b
           values_c(counter_c) = working_value_b
           counter_c = counter_c + 1
        END IF
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_b
        IF (ABS(working_value_a) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_a
           values_c(counter_c) = working_value_a
           counter_c = counter_c + 1
        END IF
        counter_a = counter_a + 1
     END IF
  END DO sum_loop

  !! Handle case where one was blank
  cleanup_a: DO WHILE (counter_a .LE. SIZE(inner_index_a))
     inner_index_c(counter_c) = inner_index_a(counter_a)
     values_c(counter_c) = values_a(counter_a)*alpha
     counter_a = counter_a + 1
     counter_c = counter_c + 1
  END DO cleanup_a
  cleanup_b: DO WHILE (counter_b .LE. SIZE(inner_index_b))
     inner_index_c(counter_c) = inner_index_b(counter_b)
     values_c(counter_c) = values_b(counter_b)
     counter_b = counter_b + 1
     counter_c = counter_c + 1
  END DO cleanup_b

  total_values_c = counter_c - 1

  END SUBROUTINE AddSparseVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Add together two sparse vectors. C = A + alpha*B
  !> The values that are returned for C are only valid in the range
  !> (1:total_values_c). We do not do an automatic shrinking of the array
  !> to keep this routine low in overhead.
  PURE SUBROUTINE AddSparseVectors_c(inner_index_a,values_a,inner_index_b, &
       & values_b,inner_index_c,values_c,total_values_c, alpha_in, threshold_in)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_b
    !> List of values computed for C.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(OUT) :: values_c
    !> Value to scale VecB by. Optional, default is 1.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> for flushing values to zero. Default value is 0.0.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> The total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Local Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

  !! Local Data
  REAL(NTREAL) :: alpha
  REAL(NTREAL) :: threshold
  !! Temporary Variables
  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b, counter_c

  !! Process Optional Parameters
  IF (.NOT. PRESENT(alpha_in)) THEN
     alpha = 1.0d+0
  ELSE
     alpha = alpha_in
  END IF
  IF (.NOT. PRESENT(threshold_in)) THEN
     threshold = 0.0d+0
  ELSE
     threshold = threshold_in
  END IF

  counter_a = 1
  counter_b = 1
  counter_c = 1
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = alpha*values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        IF (ABS(working_value_a + working_value_b) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_a
           values_c(counter_c) = working_value_a + working_value_b
           counter_c = counter_c + 1
        END IF
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        IF (ABS(working_value_b) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_b
           values_c(counter_c) = working_value_b
           counter_c = counter_c + 1
        END IF
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_b
        IF (ABS(working_value_a) .GT. threshold) THEN
           inner_index_c(counter_c) = working_index_a
           values_c(counter_c) = working_value_a
           counter_c = counter_c + 1
        END IF
        counter_a = counter_a + 1
     END IF
  END DO sum_loop

  !! Handle case where one was blank
  cleanup_a: DO WHILE (counter_a .LE. SIZE(inner_index_a))
     inner_index_c(counter_c) = inner_index_a(counter_a)
     values_c(counter_c) = values_a(counter_a)*alpha
     counter_a = counter_a + 1
     counter_c = counter_c + 1
  END DO cleanup_a
  cleanup_b: DO WHILE (counter_b .LE. SIZE(inner_index_b))
     inner_index_c(counter_c) = inner_index_b(counter_b)
     values_c(counter_c) = values_b(counter_b)
     counter_b = counter_b + 1
     counter_c = counter_c + 1
  END DO cleanup_b

  total_values_c = counter_c - 1

  END SUBROUTINE AddSparseVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(A,B)
  PURE FUNCTION DotSparseVectors_r(inner_index_a,values_a,inner_index_b, &
       & values_b) RESULT(product)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of values for A.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_a
    !> List of values for B.
    REAL(NTREAL), DIMENSION(:), INTENT(IN) :: values_b
    !> Dot product.
    REAL(NTREAL) :: product
    !! Temporary Variables
    REAL(NTREAL) :: working_value_a, working_value_b

  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b

  counter_a = 1
  counter_b = 1
  product = 0
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        product = product + working_value_a * working_value_b
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_b
        counter_a = counter_a + 1
     END IF
  END DO sum_loop

  END FUNCTION DotSparseVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(A,B)
  PURE FUNCTION DotSparseVectors_c(inner_index_a,values_a,inner_index_b, &
       & values_b) RESULT(product)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of values for A.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN) :: values_a
    !> List of values for B.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN) :: values_b
    !> Dot product.
    COMPLEX(NTCOMPLEX) :: product
    !! Temporary Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b

  counter_a = 1
  counter_b = 1
  product = 0
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        product = product + working_value_a * working_value_b
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_b
        counter_a = counter_a + 1
     END IF
  END DO sum_loop

  END FUNCTION DotSparseVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply Vectors C = A Mult B
  PURE SUBROUTINE PairwiseMultiplyVectors_r(inner_index_a,values_a, &
       & inner_index_b,values_b,inner_index_c,values_c,total_values_c)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices computed for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    REAL(NTREAL), DIMENSION(:), INTENT(IN)  :: values_b
    !> List of values computed for C.
    REAL(NTREAL), DIMENSION(:), INTENT(OUT) :: values_c
    !> This is the total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Temporary Variables
    REAL(NTREAL) :: working_value_a, working_value_b

  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b, counter_c

  counter_a = 1
  counter_b = 1
  counter_c = 1
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        inner_index_c(counter_c) = working_index_a
        values_c(counter_c) = working_value_a * working_value_b
        counter_c = counter_c + 1
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_b
        counter_a = counter_a + 1
     END IF
  END DO sum_loop
  total_values_c = counter_c - 1

  END SUBROUTINE PairwiseMultiplyVectors_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Pairwise Multiply Vectors C = A Mult B
  PURE SUBROUTINE PairwiseMultiplyVectors_c(inner_index_a,values_a, &
       & inner_index_b, values_b,inner_index_c,values_c,total_values_c)
    !> List of indices for A.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_a
    !> List of indices for B.
    INTEGER, DIMENSION(:), INTENT(IN)  :: inner_index_b
    !> List of indices computed for C.
    INTEGER, DIMENSION(:), INTENT(OUT) :: inner_index_c
    !> List of values for A.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_a
    !> List of values for B.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(IN)  :: values_b
    !> This is the total number of values in C.
    COMPLEX(NTCOMPLEX), DIMENSION(:), INTENT(OUT) :: values_c
    !> This is the total number of values in C.
    INTEGER, INTENT(OUT) :: total_values_c
    !! Temporary Variables
    COMPLEX(NTCOMPLEX) :: working_value_a, working_value_b

  INTEGER :: working_index_a, working_index_b
  !! Counter Variables
  INTEGER :: counter_a, counter_b, counter_c

  counter_a = 1
  counter_b = 1
  counter_c = 1
  sum_loop: DO WHILE(counter_a .LE. SIZE(inner_index_a) .AND. counter_b .LE. &
       & SIZE(inner_index_b))
     !! Current inner indices and values
     working_index_a = inner_index_a(counter_a)
     working_index_b = inner_index_b(counter_b)
     working_value_a = values_a(counter_a)
     working_value_b = values_b(counter_b)
     !! Figure out from which vector an insertion will be performed
     IF (working_index_a .EQ. working_index_b) THEN
        inner_index_c(counter_c) = working_index_a
        values_c(counter_c) = working_value_a * working_value_b
        counter_c = counter_c + 1
        counter_a = counter_a + 1
        counter_b = counter_b + 1
     ELSE IF (working_index_a .GT. working_index_b) THEN
        counter_b = counter_b + 1
     ELSE !! implies working_index_b > working_index_b
        counter_a = counter_a + 1
     END IF
  END DO sum_loop
  total_values_c = counter_c - 1

  END SUBROUTINE PairwiseMultiplyVectors_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SVectorModule
