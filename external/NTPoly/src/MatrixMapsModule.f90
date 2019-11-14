

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Simplfiying Per Element Operations on Matrices.
MODULE MatrixMapsModule
  USE DataTypesModule, ONLY : NTREAL, NTCOMPLEX
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, &
       & GetMatrixTripletList, FillMatrixFromTripletList
  USE TripletListModule
  USE TripletModule, ONLY : Triplet_r, Triplet_c
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MapMatrix_psr
  PUBLIC :: MapMatrix_psc
  PUBLIC :: MapTripletList
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MapMatrix_psr
     MODULE PROCEDURE MapMatrix_psr
     MODULE PROCEDURE MapMatrixArray_psr
  END INTERFACE
  INTERFACE MapMatrix_psc
     MODULE PROCEDURE MapMatrix_psc
     MODULE PROCEDURE MapMatrixArray_psc
  END INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MapTripletList
     MODULE PROCEDURE MapTripletList_r
     MODULE PROCEDURE MapTripletList_c
     MODULE PROCEDURE MapTripletListArray_r
     MODULE PROCEDURE MapTripletListArray_c
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a distributed matrix, apply this procedure to each element (real).
  SUBROUTINE MapMatrix_psr(inmat, outmat, proc)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val) RESULT(valid)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !! Local Variables
    TYPE(TripletList_r) :: inlist, outlist

  !! Convert to a triplet list, map the triplet list, fill.
  CALL ConstructEmptyMatrix(outmat, inmat)
  CALL GetMatrixTripletList(inmat, inlist)
  CALL MapTripletList(inlist, outlist, proc, &
       & num_slices_in=inmat%process_grid%num_process_slices, &
       & my_slice_in=inmat%process_grid%my_slice)
  CALL FillMatrixFromTripletList(outmat, outlist, preduplicated_in=.FALSE.)

  !! Cleanup
  CALL DestructTripletList(inlist)
  CALL DestructTripletList(outlist)
  END SUBROUTINE MapMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a distributed matrix, apply this procedure to each element (complex).
  SUBROUTINE MapMatrix_psc(inmat, outmat, proc)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val) RESULT(valid)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !! Local Variables
    TYPE(TripletList_c) :: inlist, outlist

  !! Convert to a triplet list, map the triplet list, fill.
  CALL ConstructEmptyMatrix(outmat, inmat)
  CALL GetMatrixTripletList(inmat, inlist)
  CALL MapTripletList(inlist, outlist, proc, &
       & num_slices_in=inmat%process_grid%num_process_slices, &
       & my_slice_in=inmat%process_grid%my_slice)
  CALL FillMatrixFromTripletList(outmat, outlist, preduplicated_in=.FALSE.)

  !! Cleanup
  CALL DestructTripletList(inlist)
  CALL DestructTripletList(outlist)
  END SUBROUTINE MapMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletList_r(inlist, outlist, proc, num_slices_in, &
       & my_slice_in)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_r), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_r), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val) RESULT(valid)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> How many process slices to do this mapping on (default is 1)
    INTEGER, INTENT(IN), OPTIONAL :: num_slices_in
    !> What process slice this process should compute (default is 0).
    INTEGER, INTENT(IN), OPTIONAL :: my_slice_in
    !! Local Variables
    TYPE(Triplet_r) :: temp

  INTEGER :: II
  LOGICAL :: valid
  INTEGER :: num_slices
  INTEGER :: my_slice

  IF (PRESENT(num_slices_in)) THEN
     num_slices = num_slices_in
     IF (PRESENT(my_slice_in)) THEN
        my_slice = my_slice_in
     ELSE
        my_slice = 0
     END IF
  ELSE
     num_slices = 1
     my_slice = 1
  END IF

  CALL ConstructTripletList(outlist)
  DO II = my_slice+1, inlist%CurrentSize, num_slices
     CALL GetTripletAt(inlist, II, temp)
     valid = proc(temp%index_row, temp%index_column, temp%point_value)
     IF (valid) THEN
        CALL AppendToTripletList(outlist, temp)
     END IF
  END DO
  END SUBROUTINE MapTripletList_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletList_c(inlist, outlist, proc, num_slices_in, &
       & my_slice_in)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_c), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_c), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val) RESULT(valid)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> How many process slices to do this mapping on (default is 1)
    INTEGER, INTENT(IN), OPTIONAL :: num_slices_in
    !> What process slice this process should compute (default is 0).
    INTEGER, INTENT(IN), OPTIONAL :: my_slice_in
    !! Local Variables
    TYPE(Triplet_c) :: temp

  INTEGER :: II
  LOGICAL :: valid
  INTEGER :: num_slices
  INTEGER :: my_slice

  IF (PRESENT(num_slices_in)) THEN
     num_slices = num_slices_in
     IF (PRESENT(my_slice_in)) THEN
        my_slice = my_slice_in
     ELSE
        my_slice = 0
     END IF
  ELSE
     num_slices = 1
     my_slice = 1
  END IF

  CALL ConstructTripletList(outlist)
  DO II = my_slice+1, inlist%CurrentSize, num_slices
     CALL GetTripletAt(inlist, II, temp)
     valid = proc(temp%index_row, temp%index_column, temp%point_value)
     IF (valid) THEN
        CALL AppendToTripletList(outlist, temp)
     END IF
  END DO
  END SUBROUTINE MapTripletList_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a distributed matrix, apply this procedure to each element (real).
  SUBROUTINE MapMatrixArray_psr(inmat, outmat, proc, supp_in)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN) :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN) :: supp_in
    !! Local Variables
    TYPE(TripletList_r) :: inlist, outlist

  !! Convert to a triplet list, map the triplet list, fill.
  CALL ConstructEmptyMatrix(outmat, inmat)
  CALL GetMatrixTripletList(inmat, inlist)
  CALL MapTripletList(inlist, outlist, proc, supp_in=supp_in, &
       & num_slices_in=inmat%process_grid%num_process_slices, &
       & my_slice_in=inmat%process_grid%my_slice)
  CALL FillMatrixFromTripletList(outmat, outlist, preduplicated_in=.FALSE.)

  !! Cleanup
  CALL DestructTripletList(inlist)
  CALL DestructTripletList(outlist)
  END SUBROUTINE MapMatrixArray_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a distributed matrix, apply this procedure to each element (complex).
  SUBROUTINE MapMatrixArray_psc(inmat, outmat, proc, supp_in)
    !> The matrix to apply the procedure to.
    TYPE(Matrix_ps), INTENT(IN) :: inmat
    !> The matrix where each element has had proc called on it.
    TYPE(Matrix_ps), INTENT(INOUT) :: outmat
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN) :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN) :: supp_in
    !! Local Variables
    TYPE(TripletList_c) :: inlist, outlist

  !! Convert to a triplet list, map the triplet list, fill.
  CALL ConstructEmptyMatrix(outmat, inmat)
  CALL GetMatrixTripletList(inmat, inlist)
  CALL MapTripletList(inlist, outlist, proc, supp_in=supp_in, &
       & num_slices_in=inmat%process_grid%num_process_slices, &
       & my_slice_in=inmat%process_grid%my_slice)
  CALL FillMatrixFromTripletList(outmat, outlist, preduplicated_in=.FALSE.)

  !! Cleanup
  CALL DestructTripletList(inlist)
  CALL DestructTripletList(outlist)
  END SUBROUTINE MapMatrixArray_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletListArray_r(inlist, outlist, proc, supp_in, &
       & num_slices_in, my_slice_in)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_r), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_r), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTREAL
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         REAL(KIND=NTREAL), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN) :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    REAL(KIND=NTREAL), DIMENSION(:), INTENT(IN) :: supp_in
    !> How many process slices to do this mapping on (default is 1)
    INTEGER, INTENT(IN), OPTIONAL :: num_slices_in
    !> What process slice this process should compute (default is 0).
    INTEGER, INTENT(IN), OPTIONAL :: my_slice_in
    !! Local Variables
    TYPE(Triplet_r) :: temp

  INTEGER :: II
  LOGICAL :: valid
  INTEGER :: num_slices
  INTEGER :: my_slice

  IF (PRESENT(num_slices_in)) THEN
     num_slices = num_slices_in
     IF (PRESENT(my_slice_in)) THEN
        my_slice = my_slice_in
     ELSE
        my_slice = 0
     END IF
  ELSE
     num_slices = 1
     my_slice = 1
  END IF

  CALL ConstructTripletList(outlist)
  DO II = my_slice+1, inlist%CurrentSize, num_slices
     CALL GetTripletAt(inlist, II, temp)
     valid = proc(temp%index_row, temp%index_column, temp%point_value, supp_in)
     IF (valid) THEN
        CALL AppendToTripletList(outlist, temp)
     END IF
  END DO
  END SUBROUTINE MapTripletListArray_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Given a triplet list, apply this procedure to each element.
  SUBROUTINE MapTripletListArray_c(inlist, outlist, proc, supp_in, &
       & num_slices_in, my_slice_in)
    !> The matrix to apply the procedure to.
    TYPE(TripletList_c), INTENT(IN) :: inlist
    !> The matrix where each element has had proc called on it.
    TYPE(TripletList_c), INTENT(INOUT) :: outlist
    INTERFACE
       !> The procedure to apply to each element.
       FUNCTION proc(row, column, val, supp_in) RESULT(valid)
         USE DataTypesModule, ONLY : NTCOMPLEX
         !> The row value of an element.
         INTEGER, INTENT(INOUT) :: row
         !> The column value of an element.
         INTEGER, INTENT(INOUT) :: column
         !> The actual value of an element.
         COMPLEX(KIND=NTCOMPLEX), INTENT(INOUT) :: val
         !> Any supplementary data you need to pass the map can packed here.
         COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN) :: supp_in
         !> Set this to false to filter an element.
         LOGICAL :: valid
       END FUNCTION proc
    END INTERFACE
    !> Any supplementary data you need to pass the map can packed here.
    COMPLEX(KIND=NTCOMPLEX), DIMENSION(:), INTENT(IN) :: supp_in
    !> How many process slices to do this mapping on (default is 1)
    INTEGER, INTENT(IN), OPTIONAL :: num_slices_in
    !> What process slice this process should compute (default is 0).
    INTEGER, INTENT(IN), OPTIONAL :: my_slice_in
    !! Local Variables
    TYPE(Triplet_c) :: temp

  INTEGER :: II
  LOGICAL :: valid
  INTEGER :: num_slices
  INTEGER :: my_slice

  IF (PRESENT(num_slices_in)) THEN
     num_slices = num_slices_in
     IF (PRESENT(my_slice_in)) THEN
        my_slice = my_slice_in
     ELSE
        my_slice = 0
     END IF
  ELSE
     num_slices = 1
     my_slice = 1
  END IF

  CALL ConstructTripletList(outlist)
  DO II = my_slice+1, inlist%CurrentSize, num_slices
     CALL GetTripletAt(inlist, II, temp)
     valid = proc(temp%index_row, temp%index_column, temp%point_value, supp_in)
     IF (valid) THEN
        CALL AppendToTripletList(outlist, temp)
     END IF
  END DO
  END SUBROUTINE MapTripletListArray_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixMapsModule
