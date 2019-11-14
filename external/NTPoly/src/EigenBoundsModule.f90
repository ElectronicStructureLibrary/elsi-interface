

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module for computing estimates of the bounds of the spectrum of a matrix.
MODULE EigenBoundsModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteElement, &
       & WriteListElement, WriteHeader
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, DotMatrix, &
       & IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, GetMatrixTripletList, FillMatrixFromTripletList
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters, &
       & SolverParameters_init
  USE TripletListModule
  USE TripletModule, ONLY : Triplet_r
  USE NTMPIModule
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GershgorinBounds
  PUBLIC :: PowerBounds
CONTAINS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the minimum and maximum eigenvalue of a matrix.
  !> Uses the Gershgorin theorem.
  SUBROUTINE GershgorinBounds(this,min_value,max_value)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> A lower bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: min_value
    !> An uppder bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !! Local Data
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_min
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: per_column_max
    !! Counters/Temporary
    INTEGER :: counter
    INTEGER :: local_column
    INTEGER :: ierr

    IF (this%is_complex) THEN
  !! Allocate Space For Result
  ALLOCATE(per_column_min(this%local_columns))
  ALLOCATE(per_column_max(this%local_columns))

  !! Compute The Local Contribution
  per_column_min = 0
  per_column_max = 0
  CALL GetMatrixTripletList(this, triplet_list_c)
  DO counter = 1, triplet_list_c%CurrentSize
     local_column = triplet_list_c%data(counter)%index_column - &
          & this%start_column + 1
     IF (triplet_list_c%data(counter)%index_row .EQ. &
          & triplet_list_c%data(counter)%index_column) THEN
        per_column_min(local_column) = per_column_min(local_column) + &
             & REAL(triplet_list_c%data(counter)%point_value,KIND=NTREAL)
        per_column_max(local_column) = per_column_max(local_column) + &
             & REAL(triplet_list_c%data(counter)%point_value,KIND=NTREAL)
     ELSE
        per_column_min(local_column) = per_column_min(local_column) - &
             & ABS(triplet_list_c%data(counter)%point_value)
        per_column_max(local_column) = per_column_max(local_column) + &
             & ABS(triplet_list_c%data(counter)%point_value)
     END IF
  END DO

  !! Sum Along Columns
  CALL MPI_Allreduce(MPI_IN_PLACE,per_column_min,SIZE(per_column_min), &
       & MPINTREAL,MPI_SUM,this%process_grid%column_comm,ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE,per_column_max,SIZE(per_column_max), &
       & MPINTREAL,MPI_SUM,this%process_grid%column_comm,ierr)

  min_value = MINVAL(per_column_min)
  max_value = MAXVAL(per_column_max)

  CALL MPI_Allreduce(MPI_IN_PLACE,min_value,1,MPINTREAL,MPI_MIN, &
       & this%process_grid%row_comm, ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE,max_value,1,MPINTREAL,MPI_MAX, &
       & this%process_grid%row_comm, ierr)

  CALL DestructTripletList(triplet_list_c)
  DEALLOCATE(per_column_min)
  DEALLOCATE(per_column_max)
    ELSE
  !! Allocate Space For Result
  ALLOCATE(per_column_min(this%local_columns))
  ALLOCATE(per_column_max(this%local_columns))

  !! Compute The Local Contribution
  per_column_min = 0
  per_column_max = 0
  CALL GetMatrixTripletList(this, triplet_list_r)
  DO counter = 1, triplet_list_r%CurrentSize
     local_column = triplet_list_r%data(counter)%index_column - &
          & this%start_column + 1
     IF (triplet_list_r%data(counter)%index_row .EQ. &
          & triplet_list_r%data(counter)%index_column) THEN
        per_column_min(local_column) = per_column_min(local_column) + &
             & REAL(triplet_list_r%data(counter)%point_value,KIND=NTREAL)
        per_column_max(local_column) = per_column_max(local_column) + &
             & REAL(triplet_list_r%data(counter)%point_value,KIND=NTREAL)
     ELSE
        per_column_min(local_column) = per_column_min(local_column) - &
             & ABS(triplet_list_r%data(counter)%point_value)
        per_column_max(local_column) = per_column_max(local_column) + &
             & ABS(triplet_list_r%data(counter)%point_value)
     END IF
  END DO

  !! Sum Along Columns
  CALL MPI_Allreduce(MPI_IN_PLACE,per_column_min,SIZE(per_column_min), &
       & MPINTREAL,MPI_SUM,this%process_grid%column_comm,ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE,per_column_max,SIZE(per_column_max), &
       & MPINTREAL,MPI_SUM,this%process_grid%column_comm,ierr)

  min_value = MINVAL(per_column_min)
  max_value = MAXVAL(per_column_max)

  CALL MPI_Allreduce(MPI_IN_PLACE,min_value,1,MPINTREAL,MPI_MIN, &
       & this%process_grid%row_comm, ierr)
  CALL MPI_Allreduce(MPI_IN_PLACE,max_value,1,MPINTREAL,MPI_MAX, &
       & this%process_grid%row_comm, ierr)

  CALL DestructTripletList(triplet_list_r)
  DEALLOCATE(per_column_min)
  DEALLOCATE(per_column_max)
    END IF
  END SUBROUTINE GershgorinBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute a bounds on the maximum eigenvalue of a matrix.
  !> Uses The Power Method.
  SUBROUTINE PowerBounds(this,max_value,solver_parameters_in)
    !> The matrix to compute the min/max of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> An upper bound on the eigenspectrum.
    REAL(NTREAL), INTENT(OUT) :: max_value
    !> The parameters for this calculation.
    TYPE(SolverParameters_t), INTENT(IN), OPTIONAL :: solver_parameters_in
    !! Handling Optional Parameters
    TYPE(SolverParameters_t) :: solver_parameters
    !! Local Data
    TYPE(Matrix_ps) :: vector, vector2, TempMat
    REAL(NTREAL) :: scale_value
    REAL(NTREAL) :: norm_value
    TYPE(TripletList_r) :: temp_list
    TYPE(Triplet_r) :: temp_triplet
    INTEGER :: outer_counter
    TYPE(MatrixMemoryPool_p) :: pool

    !! Optional Parameters
    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_init()
       solver_parameters%max_iterations = 10
    END IF

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Power Bounds Solver")
       CALL EnterSubLog
       CALL PrintParameters(solver_parameters)
    END IF

    !! Diagonal matrices serve as vectors.
    CALL ConstructEmptyMatrix(vector, this)
    CALL ConstructEmptyMatrix(vector2, this)

    !! Guess Vector
    temp_list = ConstructTripletList_r()
    IF (this%process_grid%global_rank .EQ. 0) THEN
       temp_triplet%index_row = 1
       temp_triplet%index_column = 1
       temp_triplet%point_value = 1.0_NTREAL
       CALL AppendToTripletList(temp_list,temp_triplet)
    END IF
    CALL FillMatrixFromTripletList(vector,temp_list)

    !! Iterate
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       IF (solver_parameters%be_verbose .AND. outer_counter .GT. 1) THEN
          CALL WriteListElement(key="Round", value=outer_counter-1)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", value=norm_value)
          CALL ExitSubLog
       END IF

       !! x = Ax
       CALL MatrixMultiply(this,vector,vector2, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
       !! x = x/||x||
       scale_value = 1.0/MatrixNorm(vector2)
       CALL ScaleMatrix(vector2,scale_value)

       !! Check if Converged
       CALL IncrementMatrix(vector2,vector,-1.0_NTREAL)
       norm_value = MatrixNorm(vector)

       CALL CopyMatrix(vector2,vector)

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",value=outer_counter-1)
    END IF

    !! Compute The Largest Eigenvalue
    CALL DotMatrix(vector, vector, scale_value)
    CALL MatrixMultiply(this,vector,vector2, &
         & threshold_in=solver_parameters%threshold, memory_pool_in=pool)
    CALL DotMatrix(vector, vector2, max_value)
    max_value = max_value / scale_value

    IF (solver_parameters%be_verbose) THEN
       CALL WriteElement(key="Max_Eigen_Value",value=max_value)
       CALL ExitSubLog
    END IF

    !! Cleanup
    CALL DestructMatrix(vector)
    CALL DestructMatrix(vector2)
    CALL DestructMatrix(TempMat)
    CALL DestructMatrixMemoryPool(pool)
  END SUBROUTINE PowerBounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE EigenBoundsModule
