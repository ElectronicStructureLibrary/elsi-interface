!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Computing The Square Root of a Matrix.
MODULE SquareRootSolversModule
  USE DataTypesModule, ONLY : NTREAL
  USE EigenBoundsModule, ONLY : GershgorinBounds
  USE LoadBalancerModule, ONLY : PermuteMatrix, UndoPermuteMatrix
  USE LoggingModule, ONLY : EnterSubLog, ExitSubLog, WriteListElement, &
       & WriteHeader, WriteElement, WriteCitation
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & DestructMatrixMemoryPool
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixNorm, &
       & IncrementMatrix, ScaleMatrix
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, FillMatrixIdentity, PrintMatrixInformation
  USE SolverParametersModule, ONLY : SolverParameters_t, PrintParameters
  USE MPI
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Solvers
  PUBLIC :: SquareRoot
  PUBLIC :: InverseSquareRoot
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root of a matrix.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] solver_parameters_in parameters for the solver, optional.
  !! @param[in] order_in polynomial order for calculation (default 5).
  SUBROUTINE SquareRoot(InputMat, OutputMat, solver_parameters_in, order_in)
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(SolverParameters_t),INTENT(in),OPTIONAL :: solver_parameters_in
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    TYPE(SolverParameters_t) :: solver_parameters
    INTEGER :: order

    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .FALSE.,&
            & order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .FALSE.)
    END IF

  END SUBROUTINE SquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the inverse square root of a matrix.
  !! @param[in] InputMat the matrix to compute.
  !! @param[out] OutputMat the resulting matrix.
  !! @param[in] solver_parameters_in parameters for the solver, optional.
  !! @param[in] order_in polynomial order for calculation (default 5).
  SUBROUTINE InverseSquareRoot(InputMat, OutputMat, solver_parameters_in, &
       & order_in)
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(SolverParameters_t),INTENT(in),OPTIONAL :: solver_parameters_in
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    TYPE(SolverParameters_t) :: solver_parameters
    INTEGER :: order

    IF (PRESENT(solver_parameters_in)) THEN
       solver_parameters = solver_parameters_in
    ELSE
       solver_parameters = SolverParameters_t()
    END IF

    IF (PRESENT(order_in)) THEN
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .TRUE.,&
            & order_in)
    ELSE
       CALL SquareRootSelector(InputMat, OutputMat, solver_parameters, .TRUE.)
    END IF

  END SUBROUTINE InverseSquareRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> This routine picks the appropriate solver method
  !! @param[in] InputMat the matrix to compute.
  !! @param[inout] OutputMat the matrix computed.
  !! @param[in] solver_parameters parameters about how to solve.
  !! @param[in] compute_inverse true if we are computing the inverse square root
  !! @param[in] order_in the polynomial degree to use (optional, default=5)
  SUBROUTINE SquareRootSelector(InputMat, OutputMat, solver_parameters, &
       & compute_inverse, order_in)
    TYPE(Matrix_ps), INTENT(in)  :: InputMat
    TYPE(Matrix_ps), INTENT(inout) :: OutputMat
    TYPE(SolverParameters_t),INTENT(in) :: solver_parameters
    LOGICAL, INTENT(IN) :: compute_inverse
    INTEGER, INTENT(IN), OPTIONAL :: order_in
    !! Local Variables
    INTEGER :: order

    IF (PRESENT(order_in)) THEN
       order = order_in
    ELSE
       order = 5
    END IF

    SELECT CASE(order)
    CASE(2)
       CALL NewtonSchultzISROrder2(InputMat, OutputMat, solver_parameters, &
            & compute_inverse)
    CASE DEFAULT
       CALL NewtonSchultzISRTaylor(InputMat, OutputMat, solver_parameters, &
            & order, compute_inverse)
    END SELECT

  END SUBROUTINE SquareRootSelector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root or inverse square root of a matrix.
  !! Based on the Newton-Schultz algorithm presented in: \cite jansik2007linear
  !! @param[in] Mat Matrix.
  !! @param[out] OutMat = Mat^-1/2 or Mat^1/2.
  !! @param[in] solver_parameters parameters for the solver
  !! @param[in] compute_inverse whether to compute the inverse square root.
  SUBROUTINE NewtonSchultzISROrder2(Mat, OutMat, solver_parameters, &
       & compute_inverse)
    !! Parameters
    TYPE(Matrix_ps), INTENT(in)  :: Mat
    TYPE(Matrix_ps), INTENT(inout) :: OutMat
    TYPE(SolverParameters_t), INTENT(in) :: solver_parameters
    LOGICAL, INTENT(in) :: compute_inverse
    !! Local Variables
    REAL(NTREAL) :: lambda
    TYPE(Matrix_ps) :: X_k,T_k,Temp,Identity
    TYPE(Matrix_ps) :: SquareRootMat
    TYPE(Matrix_ps) :: InverseSquareRootMat
    !! Temporary Variables
    REAL(NTREAL) :: e_min, e_max
    REAL(NTREAL) :: max_between
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool1

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteCitation("jansik2007linear")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(X_k, Mat)
    CALL ConstructEmptyMatrix(SquareRootMat, Mat)
    CALL ConstructEmptyMatrix(InverseSquareRootMat, Mat)
    CALL ConstructEmptyMatrix(T_k, Mat)
    CALL ConstructEmptyMatrix(Temp, Mat)
    CALL ConstructEmptyMatrix(Identity, Mat)
    CALL FillMatrixIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(Mat,e_min,e_max)
    max_between = MAX(ABS(e_min),ABS(e_max))
    lambda = 1.0_NTREAL/max_between

    !! Initialize
    CALL FillMatrixIdentity(InverseSquareRootMat)
    CALL CopyMatrix(Mat,SquareRootMat)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat, SquareRootMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(Identity, Identity, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
       CALL PermuteMatrix(InverseSquareRootMat, InverseSquareRootMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k
       CALL MatrixMultiply(SquareRootMat,InverseSquareRootMat,X_k, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL GershgorinBounds(X_k,e_min,e_max)
       max_between = MAX(ABS(e_min),ABS(e_max))
       lambda = 1.0_NTREAL/max_between

       CALL ScaleMatrix(X_k,lambda)

       !! Check if Converged
       CALL CopyMatrix(Identity,Temp)
       CALL IncrementMatrix(X_k,Temp,-1.0_NTREAL)
       norm_value = MatrixNorm(Temp)

       !! Compute T_k
       CALL CopyMatrix(Identity,T_k)
       CALL ScaleMatrix(T_k,3.0_NTREAL)
       CALL IncrementMatrix(X_k,T_k,-1.0_NTREAL)
       CALL ScaleMatrix(T_k,0.5_NTREAL)

       !! Compute Z_k+1
       CALL CopyMatrix(InverseSquareRootMat,Temp)
       CALL MatrixMultiply(Temp,T_k,InverseSquareRootMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL ScaleMatrix(InverseSquareRootMat,SQRT(lambda))

       !! Compute Y_k+1
       CALL CopyMatrix(SquareRootMat, Temp)
       CALL MatrixMultiply(T_k,Temp,SquareRootMat, &
            & threshold_in=solver_parameters%threshold, memory_pool_in=pool1)
       CALL ScaleMatrix(SquareRootMat,SQRT(lambda))

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round", int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence", float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL CopyMatrix(InverseSquareRootMat, OutMat)
    ELSE
       CALL CopyMatrix(SquareRootMat, OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat, OutMat, &
            & solver_parameters%BalancePermutation, memorypool_in=pool1)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructMatrix(Temp)
    CALL DestructMatrix(X_k)
    CALL DestructMatrix(SquareRootMat)
    CALL DestructMatrix(InverseSquareRootMat)
    CALL DestructMatrix(T_k)
    CALL DestructMatrixMemoryPool(pool1)
  END SUBROUTINE NewtonSchultzISROrder2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the square root or inverse square root of a matrix.
  !! Based on the Newton-Schultz algorithm with higher order polynomials.
  !! @param[in] Mat Matrix.
  !! @param[out] OutMat = Mat^-1/2 or Mat^1/2.
  !! @param[in] solver_parameters parameters for the solver
  !! @param[in] compute_inverse whether to compute the inverse square root.
  SUBROUTINE NewtonSchultzISRTaylor(Mat, OutMat, solver_parameters, &
       & taylor_order, compute_inverse)
    !! Parameters
    TYPE(Matrix_ps), INTENT(in)  :: Mat
    TYPE(Matrix_ps), INTENT(inout) :: OutMat
    TYPE(SolverParameters_t), INTENT(in) :: solver_parameters
    INTEGER, INTENT(in) :: taylor_order
    LOGICAL, INTENT(in) :: compute_inverse
    !! Local Variables
    REAL(NTREAL) :: lambda
    REAL(NTREAL) :: aa,bb,cc,dd
    REAL(NTREAL) :: a,b,c,d
    TYPE(Matrix_ps) :: X_k,Temp,Temp2,Identity
    TYPE(Matrix_ps) :: SquareRootMat
    TYPE(Matrix_ps) :: InverseSquareRootMat
    !! Temporary Variables
    REAL(NTREAL) :: e_min,e_max
    REAL(NTREAL) :: max_between
    INTEGER :: outer_counter
    REAL(NTREAL) :: norm_value
    TYPE(MatrixMemoryPool_p) :: pool1

    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Newton Schultz Inverse Square Root")
       CALL EnterSubLog
       CALL WriteCitation("jansik2007linear")
       CALL PrintParameters(solver_parameters)
    END IF

    !! Construct All The Necessary Matrices
    CALL ConstructEmptyMatrix(X_k, Mat)
    CALL ConstructEmptyMatrix(SquareRootMat, Mat)
    CALL ConstructEmptyMatrix(InverseSquareRootMat, Mat)
    CALL ConstructEmptyMatrix(Temp, Mat)
    IF (taylor_order == 5) THEN
       CALL ConstructEmptyMatrix(Temp2, Mat)
    END IF
    CALL ConstructEmptyMatrix(Identity, Mat)
    CALL FillMatrixIdentity(Identity)

    !! Compute the lambda scaling value.
    CALL GershgorinBounds(Mat,e_min,e_max)
    max_between = MAX(ABS(e_min),ABS(e_max))
    lambda = 1.0_NTREAL/max_between

    !! Initialize
    CALL FillMatrixIdentity(InverseSquareRootMat)
    CALL CopyMatrix(Mat,SquareRootMat)
    CALL ScaleMatrix(SquareRootMat,lambda)

    !! Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL PermuteMatrix(SquareRootMat,SquareRootMat, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
       CALL PermuteMatrix(Identity,Identity, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
       CALL PermuteMatrix(InverseSquareRootMat,InverseSquareRootMat, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
    END IF

    !! Iterate.
    IF (solver_parameters%be_verbose) THEN
       CALL WriteHeader("Iterations")
       CALL EnterSubLog
    END IF
    outer_counter = 1
    norm_value = solver_parameters%converge_diff + 1.0_NTREAL
    DO outer_counter = 1,solver_parameters%max_iterations
       !! Compute X_k = Z_k * Y_k - I
       CALL MatrixMultiply(InverseSquareRootMat,SquareRootMat,X_k, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)
       CALL IncrementMatrix(Identity,X_k,-1.0_NTREAL)
       norm_value = MatrixNorm(X_k)

       SELECT CASE(taylor_order)
       CASE(3)
          !! Compute X_k^2
          CALL MatrixMultiply(X_k,X_k,Temp, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

          !! X_k = I - 1/2 X_k + 3/8 X_k^2 + ...
          CALL ScaleMatrix(X_k,-0.5_NTREAL)
          CALL IncrementMatrix(Identity,X_k)
          CALL IncrementMatrix(Temp,X_k,0.375_NTREAL)
       CASE(5)
          !! Compute p(x) = x^4 + A*x^3 + B*x^2 + C*x + D
          !! Scale to make coefficient of x^4 equal to 1
          aa = -40.0_NTREAL/35.0_NTREAL
          bb = 48.0_NTREAL/35.0_NTREAL
          cc = -64.0_NTREAL/35.0_NTREAL
          dd = 128.0_NTREAL/35.0_NTREAL

          !! Knuth's method
          !! p = (z+x+b) * (z+c) + d
          !! z = x * (x+a)
          !! a = (A-1)/2
          !! b = B*(a+1) - C - a*(a+1)*(a+1)
          !! c = B - b - a*(a+1)
          !! d = D - b*c
          a = (aa-1.0_NTREAL)/2.0_NTREAL
          b = bb*(a+1.0_NTREAL)-cc-a*(a+1.0_NTREAL)**2
          c = bb-b-a*(a+1.0_NTREAL)
          d = dd-b*c

          !! Compute Temp = z = x * (x+a)
          CALL MatrixMultiply(X_k,X_k,Temp, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)
          CALL IncrementMatrix(X_k,Temp,a)

          !! Compute Temp2 = z + x + b
          CALL CopyMatrix(Identity,Temp2)
          CALL ScaleMatrix(Temp2,b)
          CALL IncrementMatrix(X_k,Temp2)
          CALL IncrementMatrix(Temp,Temp2)

          !! Compute Temp = z + c
          CALL IncrementMatrix(Identity,Temp,c)

          !! Compute X_k = (z+x+b) * (z+c) + d = Temp2 * Temp + d
          CALL MatrixMultiply(Temp2,Temp,X_k, &
               & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)
          CALL IncrementMatrix(Identity,X_k,d)

          !! Scale back to the target coefficients
          CALL ScaleMatrix(X_k,35.0_NTREAL/128.0_NTREAL)
       END SELECT

       !! Compute Z_k+1 = Z_k * X_k
       CALL CopyMatrix(InverseSquareRootMat,Temp)
       CALL MatrixMultiply(X_k,Temp,InverseSquareRootMat, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

       !! Compute Y_k+1 = X_k * Y_k
       CALL CopyMatrix(SquareRootMat,Temp)
       CALL MatrixMultiply(Temp,X_k,SquareRootMat, &
            & threshold_in=solver_parameters%threshold,memory_pool_in=pool1)

       IF (solver_parameters%be_verbose) THEN
          CALL WriteListElement(key="Round",int_value_in=outer_counter)
          CALL EnterSubLog
          CALL WriteElement(key="Convergence",float_value_in=norm_value)
          CALL ExitSubLog
       END IF

       IF (norm_value .LE. solver_parameters%converge_diff) THEN
          EXIT
       END IF
    END DO
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
       CALL WriteElement(key="Total_Iterations",int_value_in=outer_counter)
       CALL PrintMatrixInformation(InverseSquareRootMat)
    END IF

    IF (compute_inverse) THEN
       CALL ScaleMatrix(InverseSquareRootMat,SQRT(lambda))
       CALL CopyMatrix(InverseSquareRootMat,OutMat)
    ELSE
       CALL ScaleMatrix(SquareRootMat,1.0_NTREAL/SQRT(lambda))
       CALL CopyMatrix(SquareRootMat,OutMat)
    END IF

    !! Undo Load Balancing Step
    IF (solver_parameters%do_load_balancing) THEN
       CALL UndoPermuteMatrix(OutMat,OutMat, &
            & solver_parameters%BalancePermutation,memorypool_in=pool1)
    END IF

    !! Cleanup
    IF (solver_parameters%be_verbose) THEN
       CALL ExitSubLog
    END IF

    CALL DestructMatrix(X_k)
    CALL DestructMatrix(SquareRootMat)
    CALL DestructMatrix(InverseSquareRootMat)
    CALL DestructMatrix(Temp)
    IF (taylor_order == 5) THEN
       CALL DestructMatrix(Temp2)
    END IF
    CALL DestructMatrix(Identity)
    CALL DestructMatrixMemoryPool(pool1)
  END SUBROUTINE NewtonSchultzISRTaylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE SquareRootSolversModule
