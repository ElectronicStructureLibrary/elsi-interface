

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A Module For Performing Distributed Sparse Matrix Algebra Operations.
MODULE PSMatrixAlgebraModule
  USE DataTypesModule, ONLY : NTREAL, MPINTREAL, NTCOMPLEX, MPINTCOMPLEX
  USE GemmTasksModule
  USE PMatrixMemoryPoolModule, ONLY : MatrixMemoryPool_p, &
       & CheckMemoryPoolValidity, DestructMatrixMemoryPool
  USE PSMatrixModule, ONLY : Matrix_ps, ConstructEmptyMatrix, CopyMatrix, &
       & DestructMatrix, ConvertMatrixToComplex, ConjugateMatrix, &
       & MergeMatrixLocalBlocks
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceAndComposeMatrixSizes, &
       & ReduceAndComposeMatrixData, ReduceAndComposeMatrixCleanup, &
       & ReduceANdSumMatrixSizes, ReduceAndSumMatrixData, &
       & ReduceAndSumMatrixCleanup, TestReduceSizeRequest, &
       & TestReduceInnerRequest, TestReduceDataRequest
  USE SMatrixAlgebraModule, ONLY : MatrixMultiply, MatrixGrandSum, &
       & PairwiseMultiplyMatrix, IncrementMatrix, ScaleMatrix, &
       & MatrixColumnNorm
  USE SMatrixModule, ONLY : Matrix_lsr, Matrix_lsc, DestructMatrix, CopyMatrix,&
       & TransposeMatrix, ComposeMatrixColumns, MatrixToTripletList
  USE TimerModule, ONLY : StartTimer, StopTimer
  USE TripletListModule
  USE NTMPIModule
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MatrixSigma
  PUBLIC :: MatrixMultiply
  PUBLIC :: MatrixGrandSum
  PUBLIC :: PairwiseMultiplyMatrix
  PUBLIC :: MatrixNorm
  PUBLIC :: DotMatrix
  PUBLIC :: IncrementMatrix
  PUBLIC :: ScaleMatrix
  PUBLIC :: MatrixTrace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE MatrixSigma
     MODULE PROCEDURE MatrixSigma_ps
  END INTERFACE
  INTERFACE MatrixMultiply
     MODULE PROCEDURE MatrixMultiply_ps
  END INTERFACE
  INTERFACE MatrixGrandSum
     MODULE PROCEDURE MatrixGrandSum_psr
     MODULE PROCEDURE MatrixGrandSum_psc
  END INTERFACE
  INTERFACE PairwiseMultiplyMatrix
     MODULE PROCEDURE PairwiseMultiplyMatrix_ps
  END INTERFACE
  INTERFACE MatrixNorm
     MODULE PROCEDURE MatrixNorm_ps
  END INTERFACE
  INTERFACE DotMatrix
     MODULE PROCEDURE DotMatrix_psr
     MODULE PROCEDURE DotMatrix_psc
  END INTERFACE
  INTERFACE IncrementMatrix
     MODULE PROCEDURE IncrementMatrix_ps
  END INTERFACE
  INTERFACE ScaleMatrix
     MODULE PROCEDURE ScaleMatrix_psr
     MODULE PROCEDURE ScaleMatrix_psc
  END INTERFACE
  INTERFACE MatrixTrace
     MODULE PROCEDURE MatrixTrace_psr
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute sigma for the inversion method.
  !> See \cite ozaki2001efficient for details.
  SUBROUTINE MatrixSigma_ps(this, sigma_value)
    !> The matrix to compute the sigma value of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> Sigma
    REAL(NTREAL), INTENT(OUT) :: sigma_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: column_sigma_contribution
    !! Counters/Temporary
    INTEGER :: inner_counter, outer_counter
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data_c)

  ALLOCATE(column_sigma_contribution(merged_local_data_c%columns))
  column_sigma_contribution = 0
  DO outer_counter = 1, merged_local_data_c%columns
     DO inner_counter = merged_local_data_c%outer_index(outer_counter), &
          & merged_local_data_c%outer_index(outer_counter+1)-1
        column_sigma_contribution(outer_counter) = &
             & column_sigma_contribution(outer_counter) + &
             & ABS(merged_local_data_c%values(inner_counter+1))
     END DO
  END DO
  CALL MPI_Allreduce(MPI_IN_PLACE,column_sigma_contribution,&
       & merged_local_data_c%columns,MPINTREAL,MPI_SUM, &
       & this%process_grid%column_comm, ierr)
  CALL MPI_Allreduce(MAXVAL(column_sigma_contribution),sigma_value,1, &
       & MPINTREAL,MPI_MAX, this%process_grid%row_comm, ierr)
  sigma_value = 1.0d+0/(sigma_value**2)

  DEALLOCATE(column_sigma_contribution)
  CALL DestructMatrix(merged_local_data_c)
    ELSE
  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data_r)

  ALLOCATE(column_sigma_contribution(merged_local_data_r%columns))
  column_sigma_contribution = 0
  DO outer_counter = 1, merged_local_data_r%columns
     DO inner_counter = merged_local_data_r%outer_index(outer_counter), &
          & merged_local_data_r%outer_index(outer_counter+1)-1
        column_sigma_contribution(outer_counter) = &
             & column_sigma_contribution(outer_counter) + &
             & ABS(merged_local_data_r%values(inner_counter+1))
     END DO
  END DO
  CALL MPI_Allreduce(MPI_IN_PLACE,column_sigma_contribution,&
       & merged_local_data_r%columns,MPINTREAL,MPI_SUM, &
       & this%process_grid%column_comm, ierr)
  CALL MPI_Allreduce(MAXVAL(column_sigma_contribution),sigma_value,1, &
       & MPINTREAL,MPI_MAX, this%process_grid%row_comm, ierr)
  sigma_value = 1.0d+0/(sigma_value**2)

  DEALLOCATE(column_sigma_contribution)
  CALL DestructMatrix(merged_local_data_r)
    ENDIF
  END SUBROUTINE MatrixSigma_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Multiply two matrices together, and add to the third.
  !> C := alpha*matA*matB+ beta*matC
  SUBROUTINE MatrixMultiply_ps(matA, matB ,matC, alpha_in, beta_in, &
       & threshold_in, memory_pool_in)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)        :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)        :: matB
    !> matC = alpha*matA*matB + beta*matC
    TYPE(Matrix_ps), INTENT(INOUT)     :: matC
    !> Scales the multiplication
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> Scales matrix we sum on to.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: beta_in
    !> For flushing values to zero.
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !> A memory pool for the calculation.
    TYPE(MatrixMemoryPool_p), OPTIONAL, INTENT(INOUT) :: memory_pool_in
    !! Local Versions of Optional Parameter
    TYPE(Matrix_ps) :: matAConverted
    TYPE(Matrix_ps) :: matBConverted
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: beta
    REAL(NTREAL) :: threshold
    TYPE(MatrixMemoryPool_p) :: memory_pool

    !! Handle the optional parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(beta_in)) THEN
       beta = 0.0
    ELSE
       beta = beta_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0.0
    ELSE
       threshold = threshold_in
    END IF

    !! Setup Memory Pool
    IF (PRESENT(memory_pool_in)) THEN
       IF (matA%is_complex) THEN
          IF (.NOT. CheckMemoryPoolValidity(memory_pool_in, matA)) THEN
             CALL DestructMatrixMemoryPool(memory_pool_in)
             memory_pool_in = MatrixMemoryPool_p(matA)
          END IF
       ELSE
          IF (.NOT. CheckMemoryPoolValidity(memory_pool_in, matB)) THEN
             CALL DestructMatrixMemoryPool(memory_pool_in)
             memory_pool_in = MatrixMemoryPool_p(matB)
          END IF
       END IF
    ELSE
       IF (matA%is_complex) THEN
          memory_pool = MatrixMemoryPool_p(matA)
       ELSE
          memory_pool = MatrixMemoryPool_p(matB)
       END IF
    END IF


    !! Perform Upcasting
    IF (matB%is_complex .AND. .NOT. matA%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, matAConverted)
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_ps_imp(matAConverted, matB, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_ps_imp(matAConverted, matB, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    ELSE IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, matBConverted)
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_ps_imp(matA, matBConverted, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_ps_imp(matA, matBConverted, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    ELSE
       IF (PRESENT(memory_pool_in)) THEN
          CALL MatrixMultiply_ps_imp(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool_in)
       ELSE
          CALL MatrixMultiply_ps_imp(matA, matB, matC, alpha, beta, &
               & threshold, memory_pool)
       END IF
    END IF

    CALL DestructMatrixMemoryPool(memory_pool)
    CALL DestructMatrix(matAConverted)
    CALL DestructMatrix(matBConverted)

  END SUBROUTINE MatrixMultiply_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> The actual implementation of matrix multiply is here. Takes the
  !> same parameters as the standard multiply, but nothing is optional.
  SUBROUTINE MatrixMultiply_ps_imp(matA, matB ,matC, alpha, beta, &
       & threshold, memory_pool)
    !! Parameters
    TYPE(Matrix_ps), INTENT(IN)    :: matA
    TYPE(Matrix_ps), INTENT(IN)    :: matB
    TYPE(Matrix_ps), INTENT(INOUT) :: matC
    REAL(NTREAL), INTENT(IN) :: alpha
    REAL(NTREAL), INTENT(IN) :: beta
    REAL(NTREAL), INTENT(IN) :: threshold
    TYPE(MatrixMemoryPool_p), INTENT(INOUT) :: memory_pool
    TYPE(Matrix_ps) :: matAB
    !! Temporary Matrices
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalRowContribution_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT_r
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution_r
    TYPE(Matrix_lsr), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution_r
    TYPE(Matrix_lsr), DIMENSION(:,:), ALLOCATABLE :: SliceContribution_r
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: AdjacentABlocks_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: LocalRowContribution_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredRowContribution_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredRowContributionT_c
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: TransposedBBlocks_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: LocalColumnContribution_c
    TYPE(Matrix_lsc), DIMENSION(:), ALLOCATABLE :: GatheredColumnContribution_c
    TYPE(Matrix_lsc), DIMENSION(:,:), ALLOCATABLE :: SliceContribution_c
    !! Communication Helpers
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: row_helper
    TYPE(ReduceHelper_t), DIMENSION(:), ALLOCATABLE :: column_helper
    TYPE(ReduceHelper_t), DIMENSION(:,:), ALLOCATABLE :: slice_helper
    !! For Iterating Over Local Blocks
    INTEGER :: II, II2
    INTEGER :: JJ, JJ2
    INTEGER :: duplicate_start_column, duplicate_offset_column
    INTEGER :: duplicate_start_row, duplicate_offset_row
    REAL(NTREAL) :: working_threshold
    !! Scheduling the A work
    INTEGER, DIMENSION(:), ALLOCATABLE :: ATasks
    INTEGER :: ATasks_completed
    !! Scheduling the B work
    INTEGER, DIMENSION(:), ALLOCATABLE :: BTasks
    INTEGER :: BTasks_completed
    !! Scheduling the AB work
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ABTasks
    INTEGER :: ABTasks_completed

    IF (matA%is_complex) THEN
  CALL StartTimer("GEMM")

  !! The threshold needs to be smaller if we are doing a sliced version
  !! because you might flush a value that would be kept in the summed version.
  IF (matA%process_grid%num_process_slices .GT. 1) THEN
     working_threshold = threshold/(matA%process_grid%num_process_slices*1000)
  ELSE
     working_threshold = threshold
  END IF

  !! Construct The Temporary Matrices
  CALL ConstructEmptyMatrix(matAB, matA)

  ALLOCATE(AdjacentABlocks_c(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns/&
       & matAB%process_grid%num_process_slices))
  ALLOCATE(LocalRowContribution_c(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(GatheredRowContribution_c(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(GatheredRowContributionT_c(matAB%process_grid%number_of_blocks_rows))

  ALLOCATE(TransposedBBlocks_c(matAB%process_grid%number_of_blocks_rows/&
       & matAB%process_grid%num_process_slices, &
       & matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(LocalColumnContribution_c(&
       & matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(GatheredColumnContribution_c(&
       & matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(SliceContribution_c(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns))

  !! Helpers
  ALLOCATE(row_helper(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(column_helper(matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(slice_helper(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns))

  !! Construct the task queues
  ALLOCATE(ATasks(matAB%process_grid%number_of_blocks_rows))
  DO II=1,matAB%process_grid%number_of_blocks_rows
     ATasks(II) = LocalGatherA
  END DO
  ALLOCATE(BTasks(matAB%process_grid%number_of_blocks_columns))
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     BTasks(JJ) = LocalGatherB
  END DO
  ALLOCATE(ABTasks(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns))
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     DO II=1,matAB%process_grid%number_of_blocks_rows
        ABTasks(II,JJ) = AwaitingAB
     END DO
  END DO

  !! Setup A Tasks
  duplicate_start_column = matAB%process_grid%my_slice+1
  duplicate_offset_column = matAB%process_grid%num_process_slices

  !! Setup B Tasks
  duplicate_start_row = matAB%process_grid%my_slice+1
  duplicate_offset_row = matAB%process_grid%num_process_slices

  !! Run A Tasks
  ATasks_completed = 0
  BTasks_completed = 0
  ABTasks_completed = 0
  !$OMP PARALLEL
  !$OMP MASTER
  DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
       & BTasks_completed .LT. SIZE(BTasks) .OR. &
       & ABTasks_completed .LT. SIZE(ABTasks))
     DO II=1, matAB%process_grid%number_of_blocks_rows
        SELECT CASE (ATasks(II))
        CASE(LocalGatherA)
           ATasks(II) = TaskRunningA
           !$OMP TASK DEFAULT(SHARED), PRIVATE(JJ2), FIRSTPRIVATE(II)
           !! First Align The Data We Are Working With
           DO JJ2=1, &
                & matAB%process_grid%number_of_blocks_columns/ &
                & matAB%process_grid%num_process_slices
              CALL CopyMatrix(matA%local_data_c(II, &
                   & duplicate_start_column+duplicate_offset_column*(JJ2-1)),&
                   & AdjacentABlocks_c(II,JJ2))
           END DO
           !! Then Do A Local Gather
           CALL ComposeMatrixColumns(AdjacentABlocks_c(II,:), &
                & LocalRowContribution_c(II))
           ATasks(II) = SendSizeA
           !$OMP END TASK
        CASE(SendSizeA)
           !! Then Start A Global Gather
           CALL ReduceAndComposeMatrixSizes(LocalRowContribution_c(II), &
                & matAB%process_grid%blocked_row_comm(II), &
                & GatheredRowContribution_c(II), row_helper(II))
           ATasks(II) = ComposeA
        CASE(ComposeA)
           IF (TestReduceSizeRequest(row_helper(II))) THEN
              CALL ReduceAndComposeMatrixData(LocalRowContribution_c(II), &
                   & matAB%process_grid%blocked_row_comm(II), &
                   & GatheredRowContribution_c(II), row_helper(II))
              ATasks(II) = WaitInnerA
           END IF
        CASE(WaitInnerA)
           IF (TestReduceInnerRequest(row_helper(II))) THEN
              ATasks(II) = WaitDataA
           END IF
        CASE(WaitDataA)
           IF (TestReduceDataRequest(row_helper(II))) THEN
              ATasks(II) = AdjustIndicesA
           END IF
        CASE(AdjustIndicesA)
           ATasks(II) = TaskRunningA
           !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II)
           CALL ReduceAndComposeMatrixCleanup(LocalRowContribution_c(II), &
                & GatheredRowContribution_c(II), row_helper(II))
           CALL TransposeMatrix(GatheredRowContribution_c(II), &
                & GatheredRowContributionT_c(II))
           ATasks(II) = CleanupA
           !$OMP END TASK
        CASE(CleanupA)
           ATasks(II) = FinishedA
           ATasks_completed = ATasks_completed + 1
        END SELECT
     END DO
     !! B Tasks
     DO JJ=1,matAB%process_grid%number_of_blocks_columns
        SELECT CASE (BTasks(JJ))
        CASE(LocalGatherB)
           BTasks(JJ) = TaskRunningB
           !$OMP TASK DEFAULT(SHARED), PRIVATE(II2), FIRSTPRIVATE(JJ)
           !! First Transpose The Data We Are Working With
           DO II2=1, matAB%process_grid%number_of_blocks_rows/&
                & matAB%process_grid%num_process_slices
              CALL TransposeMatrix(matB%local_data_c(duplicate_start_row+&
                   & duplicate_offset_row*(II2-1),JJ), &
                   & TransposedBBlocks_c(II2,JJ))
           END DO
           !! Then Do A Local Gather
           CALL ComposeMatrixColumns(TransposedBBlocks_c(:,JJ), &
                & LocalColumnContribution_c(JJ))
           BTasks(JJ) = SendSizeB
           !$OMP END TASK
        CASE(SendSizeB)
           !! Then A Global Gather
           CALL ReduceAndComposeMatrixSizes(LocalColumnContribution_c(JJ), &
                & matAB%process_grid%blocked_column_comm(JJ), &
                & GatheredColumnContribution_c(JJ), column_helper(JJ))
           BTasks(JJ) = LocalComposeB
        CASE(LocalComposeB)
           IF (TestReduceSizeRequest(column_helper(JJ))) THEN
              CALL ReduceAndComposeMatrixData(LocalColumnContribution_c(JJ),&
                   & matAB%process_grid%blocked_column_comm(JJ), &
                   & GatheredColumnContribution_c(JJ), column_helper(JJ))
              BTasks(JJ) = WaitInnerB
           END IF
        CASE(WaitInnerB)
           IF (TestReduceInnerRequest(column_helper(JJ))) THEN
              BTasks(JJ) = WaitDataB
           END IF
        CASE(WaitDataB)
           IF (TestReduceDataRequest(column_helper(JJ))) THEN
              BTasks(JJ) = AdjustIndicesB
           END IF
        CASE(AdjustIndicesB)
           BTasks(JJ) = TaskRunningB
           !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(JJ)
           CALL ReduceAndComposeMatrixCleanup(LocalColumnContribution_c(JJ), &
                & GatheredColumnContribution_c(JJ), column_helper(JJ))
           BTasks(JJ) = CleanupB
           !$OMP END TASK
        CASE(CleanupB)
           BTasks(JJ) = FinishedB
           BTasks_completed = BTasks_completed + 1
        END SELECT
     END DO
     !! AB Tasks
     DO II=1,matAB%process_grid%number_of_blocks_rows
        DO JJ=1,matAB%process_grid%number_of_blocks_columns
           SELECT CASE(ABTasks(II,JJ))
           CASE (AwaitingAB)
              IF (ATasks(II) .EQ. FinishedA .AND. &
                   & BTasks(JJ) .EQ. FinishedB) THEN
                 ABTasks(II,JJ) = GemmAB
              END IF
           CASE (GemmAB)
              ABTasks(II,JJ) = TaskRunningAB
              !$OMP TASK DEFAULT(shared), FIRSTPRIVATE(II,JJ)
              CALL MatrixMultiply(GatheredRowContributionT_c(II), &
                   & GatheredColumnContribution_c(JJ), &
                   & SliceContribution_c(II,JJ), &
                   & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                   & alpha_in=alpha, threshold_in=working_threshold, &
                   & blocked_memory_pool_in=memory_pool%grid_c(II,JJ))
              !! We can exit early if there is only one process slice
              IF (matAB%process_grid%num_process_slices .EQ. 1) THEN
                 ABTasks(II,JJ) = CleanupAB
                 CALL CopyMatrix(SliceContribution_c(II,JJ), matAB%local_data_c(II,JJ))
              ELSE
                 ABTasks(II,JJ) = SendSizeAB
              END IF
              !$OMP END TASK
           CASE(SendSizeAB)
              CALL ReduceAndSumMatrixSizes(SliceContribution_c(II,JJ),&
                   & matAB%process_grid%blocked_between_slice_comm(II,JJ), &
                   & matAB%local_data_c(II,JJ), slice_helper(II,JJ))
              ABTasks(II,JJ) = GatherAndSumAB
           CASE (GatherAndSumAB)
              IF (TestReduceSizeRequest(slice_helper(II,JJ))) THEN
                 CALL ReduceAndSumMatrixData(SliceContribution_c(II,JJ), &
                      & matAB%local_data_c(II,JJ), &
                      & matAB%process_grid%blocked_between_slice_comm(II,JJ),&
                      & slice_helper(II,JJ))
                 ABTasks(II,JJ) = WaitInnerAB
              END IF
           CASE (WaitInnerAB)
              IF (TestReduceInnerRequest(slice_helper(II,JJ))) THEN
                 ABTasks(II,JJ) = WaitDataAB
              END IF
           CASE (WaitDataAB)
              IF (TestReduceDataRequest(slice_helper(II,JJ))) THEN
                 ABTasks(II,JJ) = LocalSumAB
              END IF
           CASE(LocalSumAB)
              ABTasks(II,JJ) = TaskRunningAB
              !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II,JJ)
              CALL ReduceAndSumMatrixCleanup(SliceContribution_c(II,JJ), &
                   & matAB%local_data_c(II,JJ), threshold, slice_helper(II,JJ))
              ABTasks(II,JJ) = CleanupAB
              !$OMP END TASK
           CASE(CleanupAB)
              ABTasks(II,JJ) = FinishedAB
              ABTasks_completed = ABTasks_completed + 1
           END SELECT
        END DO
     END DO
  END DO
  !$OMP END MASTER
  !$OMP END PARALLEL

  !! Copy to output matrix.
  IF (beta .EQ. 0.0) THEN
     CALL CopyMatrix(matAB,matC)
  ELSE
     CALL ScaleMatrix(MatC,beta)
     CALL IncrementMatrix(MatAB,MatC)
  END IF

  !! Cleanup
  CALL DestructMatrix(matAB)
  DEALLOCATE(row_helper)
  DEALLOCATE(column_helper)
  DEALLOCATE(slice_helper)
  DEALLOCATE(ATasks)
  DEALLOCATE(BTasks)
  DEALLOCATE(ABTasks)

  !! Deallocate Buffers From A
  DO II=1,matAB%process_grid%number_of_blocks_rows
     DO JJ2=1,matAB%process_grid%number_of_blocks_columns/&
          & matAB%process_grid%num_process_slices
        CALL DestructMatrix(AdjacentABlocks_c(II,JJ2))
     END DO
     CALL DestructMatrix(LocalRowContribution_c(II))
     CALL DestructMatrix(GatheredRowContribution_c(II))
  END DO
  DEALLOCATE(AdjacentABlocks_c)
  DEALLOCATE(LocalRowContribution_c)
  DEALLOCATE(GatheredRowContribution_c)
  !! Deallocate Buffers From B
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     DO II2=1,matAB%process_grid%number_of_blocks_rows/&
          & matAB%process_grid%num_process_slices
        CALL DestructMatrix(TransposedBBlocks_c(II2,JJ))
     END DO
     CALL DestructMatrix(LocalColumnContribution_c(JJ))
  END DO
  DEALLOCATE(TransposedBBlocks_c)
  DEALLOCATE(LocalColumnContribution_c)
  !! Deallocate Buffers From Multiplying The Block
  DO II=1,matAB%process_grid%number_of_blocks_rows
     CALL DestructMatrix(GatheredRowContributionT_c(II))
  END DO
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     CALL DestructMatrix(GatheredColumnContribution_c(JJ))
  END DO
  DEALLOCATE(GatheredRowContributionT_c)
  DEALLOCATE(GatheredColumnContribution_c)
  !! Deallocate Buffers From Sum
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     DO II=1,matAB%process_grid%number_of_blocks_rows
        CALL DestructMatrix(SliceContribution_c(II,JJ))
     END DO
  END DO
  DEALLOCATE(SliceContribution_c)

  CALL StopTimer("GEMM")
    ELSE
  CALL StartTimer("GEMM")

  !! The threshold needs to be smaller if we are doing a sliced version
  !! because you might flush a value that would be kept in the summed version.
  IF (matA%process_grid%num_process_slices .GT. 1) THEN
     working_threshold = threshold/(matA%process_grid%num_process_slices*1000)
  ELSE
     working_threshold = threshold
  END IF

  !! Construct The Temporary Matrices
  CALL ConstructEmptyMatrix(matAB, matA)

  ALLOCATE(AdjacentABlocks_r(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns/&
       & matAB%process_grid%num_process_slices))
  ALLOCATE(LocalRowContribution_r(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(GatheredRowContribution_r(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(GatheredRowContributionT_r(matAB%process_grid%number_of_blocks_rows))

  ALLOCATE(TransposedBBlocks_r(matAB%process_grid%number_of_blocks_rows/&
       & matAB%process_grid%num_process_slices, &
       & matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(LocalColumnContribution_r(&
       & matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(GatheredColumnContribution_r(&
       & matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(SliceContribution_r(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns))

  !! Helpers
  ALLOCATE(row_helper(matAB%process_grid%number_of_blocks_rows))
  ALLOCATE(column_helper(matAB%process_grid%number_of_blocks_columns))
  ALLOCATE(slice_helper(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns))

  !! Construct the task queues
  ALLOCATE(ATasks(matAB%process_grid%number_of_blocks_rows))
  DO II=1,matAB%process_grid%number_of_blocks_rows
     ATasks(II) = LocalGatherA
  END DO
  ALLOCATE(BTasks(matAB%process_grid%number_of_blocks_columns))
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     BTasks(JJ) = LocalGatherB
  END DO
  ALLOCATE(ABTasks(matAB%process_grid%number_of_blocks_rows, &
       & matAB%process_grid%number_of_blocks_columns))
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     DO II=1,matAB%process_grid%number_of_blocks_rows
        ABTasks(II,JJ) = AwaitingAB
     END DO
  END DO

  !! Setup A Tasks
  duplicate_start_column = matAB%process_grid%my_slice+1
  duplicate_offset_column = matAB%process_grid%num_process_slices

  !! Setup B Tasks
  duplicate_start_row = matAB%process_grid%my_slice+1
  duplicate_offset_row = matAB%process_grid%num_process_slices

  !! Run A Tasks
  ATasks_completed = 0
  BTasks_completed = 0
  ABTasks_completed = 0
  !$OMP PARALLEL
  !$OMP MASTER
  DO WHILE (ATasks_completed .LT. SIZE(ATasks) .OR. &
       & BTasks_completed .LT. SIZE(BTasks) .OR. &
       & ABTasks_completed .LT. SIZE(ABTasks))
     DO II=1, matAB%process_grid%number_of_blocks_rows
        SELECT CASE (ATasks(II))
        CASE(LocalGatherA)
           ATasks(II) = TaskRunningA
           !$OMP TASK DEFAULT(SHARED), PRIVATE(JJ2), FIRSTPRIVATE(II)
           !! First Align The Data We Are Working With
           DO JJ2=1, &
                & matAB%process_grid%number_of_blocks_columns/ &
                & matAB%process_grid%num_process_slices
              CALL CopyMatrix(matA%local_data_r(II, &
                   & duplicate_start_column+duplicate_offset_column*(JJ2-1)),&
                   & AdjacentABlocks_r(II,JJ2))
           END DO
           !! Then Do A Local Gather
           CALL ComposeMatrixColumns(AdjacentABlocks_r(II,:), &
                & LocalRowContribution_r(II))
           ATasks(II) = SendSizeA
           !$OMP END TASK
        CASE(SendSizeA)
           !! Then Start A Global Gather
           CALL ReduceAndComposeMatrixSizes(LocalRowContribution_r(II), &
                & matAB%process_grid%blocked_row_comm(II), &
                & GatheredRowContribution_r(II), row_helper(II))
           ATasks(II) = ComposeA
        CASE(ComposeA)
           IF (TestReduceSizeRequest(row_helper(II))) THEN
              CALL ReduceAndComposeMatrixData(LocalRowContribution_r(II), &
                   & matAB%process_grid%blocked_row_comm(II), &
                   & GatheredRowContribution_r(II), row_helper(II))
              ATasks(II) = WaitInnerA
           END IF
        CASE(WaitInnerA)
           IF (TestReduceInnerRequest(row_helper(II))) THEN
              ATasks(II) = WaitDataA
           END IF
        CASE(WaitDataA)
           IF (TestReduceDataRequest(row_helper(II))) THEN
              ATasks(II) = AdjustIndicesA
           END IF
        CASE(AdjustIndicesA)
           ATasks(II) = TaskRunningA
           !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II)
           CALL ReduceAndComposeMatrixCleanup(LocalRowContribution_r(II), &
                & GatheredRowContribution_r(II), row_helper(II))
           CALL TransposeMatrix(GatheredRowContribution_r(II), &
                & GatheredRowContributionT_r(II))
           ATasks(II) = CleanupA
           !$OMP END TASK
        CASE(CleanupA)
           ATasks(II) = FinishedA
           ATasks_completed = ATasks_completed + 1
        END SELECT
     END DO
     !! B Tasks
     DO JJ=1,matAB%process_grid%number_of_blocks_columns
        SELECT CASE (BTasks(JJ))
        CASE(LocalGatherB)
           BTasks(JJ) = TaskRunningB
           !$OMP TASK DEFAULT(SHARED), PRIVATE(II2), FIRSTPRIVATE(JJ)
           !! First Transpose The Data We Are Working With
           DO II2=1, matAB%process_grid%number_of_blocks_rows/&
                & matAB%process_grid%num_process_slices
              CALL TransposeMatrix(matB%local_data_r(duplicate_start_row+&
                   & duplicate_offset_row*(II2-1),JJ), &
                   & TransposedBBlocks_r(II2,JJ))
           END DO
           !! Then Do A Local Gather
           CALL ComposeMatrixColumns(TransposedBBlocks_r(:,JJ), &
                & LocalColumnContribution_r(JJ))
           BTasks(JJ) = SendSizeB
           !$OMP END TASK
        CASE(SendSizeB)
           !! Then A Global Gather
           CALL ReduceAndComposeMatrixSizes(LocalColumnContribution_r(JJ), &
                & matAB%process_grid%blocked_column_comm(JJ), &
                & GatheredColumnContribution_r(JJ), column_helper(JJ))
           BTasks(JJ) = LocalComposeB
        CASE(LocalComposeB)
           IF (TestReduceSizeRequest(column_helper(JJ))) THEN
              CALL ReduceAndComposeMatrixData(LocalColumnContribution_r(JJ),&
                   & matAB%process_grid%blocked_column_comm(JJ), &
                   & GatheredColumnContribution_r(JJ), column_helper(JJ))
              BTasks(JJ) = WaitInnerB
           END IF
        CASE(WaitInnerB)
           IF (TestReduceInnerRequest(column_helper(JJ))) THEN
              BTasks(JJ) = WaitDataB
           END IF
        CASE(WaitDataB)
           IF (TestReduceDataRequest(column_helper(JJ))) THEN
              BTasks(JJ) = AdjustIndicesB
           END IF
        CASE(AdjustIndicesB)
           BTasks(JJ) = TaskRunningB
           !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(JJ)
           CALL ReduceAndComposeMatrixCleanup(LocalColumnContribution_r(JJ), &
                & GatheredColumnContribution_r(JJ), column_helper(JJ))
           BTasks(JJ) = CleanupB
           !$OMP END TASK
        CASE(CleanupB)
           BTasks(JJ) = FinishedB
           BTasks_completed = BTasks_completed + 1
        END SELECT
     END DO
     !! AB Tasks
     DO II=1,matAB%process_grid%number_of_blocks_rows
        DO JJ=1,matAB%process_grid%number_of_blocks_columns
           SELECT CASE(ABTasks(II,JJ))
           CASE (AwaitingAB)
              IF (ATasks(II) .EQ. FinishedA .AND. &
                   & BTasks(JJ) .EQ. FinishedB) THEN
                 ABTasks(II,JJ) = GemmAB
              END IF
           CASE (GemmAB)
              ABTasks(II,JJ) = TaskRunningAB
              !$OMP TASK DEFAULT(shared), FIRSTPRIVATE(II,JJ)
              CALL MatrixMultiply(GatheredRowContributionT_r(II), &
                   & GatheredColumnContribution_r(JJ), &
                   & SliceContribution_r(II,JJ), &
                   & IsATransposed_in=.TRUE., IsBTransposed_in=.TRUE., &
                   & alpha_in=alpha, threshold_in=working_threshold, &
                   & blocked_memory_pool_in=memory_pool%grid_r(II,JJ))
              !! We can exit early if there is only one process slice
              IF (matAB%process_grid%num_process_slices .EQ. 1) THEN
                 ABTasks(II,JJ) = CleanupAB
                 CALL CopyMatrix(SliceContribution_r(II,JJ), matAB%local_data_r(II,JJ))
              ELSE
                 ABTasks(II,JJ) = SendSizeAB
              END IF
              !$OMP END TASK
           CASE(SendSizeAB)
              CALL ReduceAndSumMatrixSizes(SliceContribution_r(II,JJ),&
                   & matAB%process_grid%blocked_between_slice_comm(II,JJ), &
                   & matAB%local_data_r(II,JJ), slice_helper(II,JJ))
              ABTasks(II,JJ) = GatherAndSumAB
           CASE (GatherAndSumAB)
              IF (TestReduceSizeRequest(slice_helper(II,JJ))) THEN
                 CALL ReduceAndSumMatrixData(SliceContribution_r(II,JJ), &
                      & matAB%local_data_r(II,JJ), &
                      & matAB%process_grid%blocked_between_slice_comm(II,JJ),&
                      & slice_helper(II,JJ))
                 ABTasks(II,JJ) = WaitInnerAB
              END IF
           CASE (WaitInnerAB)
              IF (TestReduceInnerRequest(slice_helper(II,JJ))) THEN
                 ABTasks(II,JJ) = WaitDataAB
              END IF
           CASE (WaitDataAB)
              IF (TestReduceDataRequest(slice_helper(II,JJ))) THEN
                 ABTasks(II,JJ) = LocalSumAB
              END IF
           CASE(LocalSumAB)
              ABTasks(II,JJ) = TaskRunningAB
              !$OMP TASK DEFAULT(SHARED), FIRSTPRIVATE(II,JJ)
              CALL ReduceAndSumMatrixCleanup(SliceContribution_r(II,JJ), &
                   & matAB%local_data_r(II,JJ), threshold, slice_helper(II,JJ))
              ABTasks(II,JJ) = CleanupAB
              !$OMP END TASK
           CASE(CleanupAB)
              ABTasks(II,JJ) = FinishedAB
              ABTasks_completed = ABTasks_completed + 1
           END SELECT
        END DO
     END DO
  END DO
  !$OMP END MASTER
  !$OMP END PARALLEL

  !! Copy to output matrix.
  IF (beta .EQ. 0.0) THEN
     CALL CopyMatrix(matAB,matC)
  ELSE
     CALL ScaleMatrix(MatC,beta)
     CALL IncrementMatrix(MatAB,MatC)
  END IF

  !! Cleanup
  CALL DestructMatrix(matAB)
  DEALLOCATE(row_helper)
  DEALLOCATE(column_helper)
  DEALLOCATE(slice_helper)
  DEALLOCATE(ATasks)
  DEALLOCATE(BTasks)
  DEALLOCATE(ABTasks)

  !! Deallocate Buffers From A
  DO II=1,matAB%process_grid%number_of_blocks_rows
     DO JJ2=1,matAB%process_grid%number_of_blocks_columns/&
          & matAB%process_grid%num_process_slices
        CALL DestructMatrix(AdjacentABlocks_r(II,JJ2))
     END DO
     CALL DestructMatrix(LocalRowContribution_r(II))
     CALL DestructMatrix(GatheredRowContribution_r(II))
  END DO
  DEALLOCATE(AdjacentABlocks_r)
  DEALLOCATE(LocalRowContribution_r)
  DEALLOCATE(GatheredRowContribution_r)
  !! Deallocate Buffers From B
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     DO II2=1,matAB%process_grid%number_of_blocks_rows/&
          & matAB%process_grid%num_process_slices
        CALL DestructMatrix(TransposedBBlocks_r(II2,JJ))
     END DO
     CALL DestructMatrix(LocalColumnContribution_r(JJ))
  END DO
  DEALLOCATE(TransposedBBlocks_r)
  DEALLOCATE(LocalColumnContribution_r)
  !! Deallocate Buffers From Multiplying The Block
  DO II=1,matAB%process_grid%number_of_blocks_rows
     CALL DestructMatrix(GatheredRowContributionT_r(II))
  END DO
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     CALL DestructMatrix(GatheredColumnContribution_r(JJ))
  END DO
  DEALLOCATE(GatheredRowContributionT_r)
  DEALLOCATE(GatheredColumnContribution_r)
  !! Deallocate Buffers From Sum
  DO JJ=1,matAB%process_grid%number_of_blocks_columns
     DO II=1,matAB%process_grid%number_of_blocks_rows
        CALL DestructMatrix(SliceContribution_r(II,JJ))
     END DO
  END DO
  DEALLOCATE(SliceContribution_r)

  CALL StopTimer("GEMM")
    END IF
  END SUBROUTINE MatrixMultiply_ps_imp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix into a single value.
  SUBROUTINE MatrixGrandSum_psr(this, sum)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: this
    !> The sum of all elements.
    REAL(NTREAL), INTENT(OUT) :: sum
    !! Local Data
    INTEGER :: II, JJ
    REAL(NTREAL) :: temp_r
    COMPLEX(NTCOMPLEX) :: temp_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
  sum = 0
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL MatrixGrandSum(this%local_data_c(II,JJ), temp_c)
        sum = sum + REAL(temp_c, KIND=NTREAL)
     END DO
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTREAL, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)
    ELSE
  sum = 0
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL MatrixGrandSum(this%local_data_r(II,JJ), temp_r)
        sum = sum + REAL(temp_r, KIND=NTREAL)
     END DO
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTREAL, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)
    END IF
  END SUBROUTINE MatrixGrandSum_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sum up the elements in a matrix into a single value.
  SUBROUTINE MatrixGrandSum_psc(this, sum)
    !> The matrix to compute.
    TYPE(Matrix_ps), INTENT(IN)  :: this
    !> The sum of all elements.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: sum
    !! Local Data
    INTEGER :: II, JJ
    REAL(NTREAL) :: temp_r
    COMPLEX(NTCOMPLEX) :: temp_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
  sum = 0
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL MatrixGrandSum(this%local_data_c(II,JJ), temp_c)
        sum = sum + REAL(temp_c, KIND=NTREAL)
     END DO
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTCOMPLEX, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)
    ELSE
  sum = 0
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL MatrixGrandSum(this%local_data_r(II,JJ), temp_r)
        sum = sum + REAL(temp_r, KIND=NTREAL)
     END DO
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, sum, 1, MPINTCOMPLEX, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)
    END IF
  END SUBROUTINE MatrixGrandSum_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Elementwise multiplication. C_ij = A_ij * B_ij.
  !> Also known as a Hadamard product.
  RECURSIVE SUBROUTINE PairwiseMultiplyMatrix_ps(matA, matB, matC)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> matC = MatA mult MatB.
    TYPE(Matrix_ps), INTENT(INOUT)  :: matC
    !! Local Data
    TYPE(Matrix_ps) :: converted_matrix
    INTEGER :: II, JJ

    IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, converted_matrix)
       CALL PairwiseMultiplyMatrix(matA, converted_matrix, matC)
       CALL DestructMatrix(converted_matrix)
    ELSE IF (.NOT. matA%is_complex .AND. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, converted_matrix)
       CALL PairwiseMultiplyMatrix(converted_matrix, matB, matC)
       CALL DestructMatrix(converted_matrix)
    ELSE IF (matA%is_complex .AND. matB%is_complex) THEN
  CALL ConstructEmptyMatrix(matC, matA%actual_matrix_dimension, &
       & matA%process_grid, matA%is_complex)

  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, matA%process_grid%number_of_blocks_columns
     DO II = 1, matA%process_grid%number_of_blocks_rows
        CALL PairwiseMultiplyMatrix(matA%local_data_c(II,JJ), matB%local_data_c(II,JJ), &
             & matC%local_data_c(II,JJ))
     END DO
  END DO
  !$omp end do
  !$omp end parallel
    ELSE
  CALL ConstructEmptyMatrix(matC, matA%actual_matrix_dimension, &
       & matA%process_grid, matA%is_complex)

  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, matA%process_grid%number_of_blocks_columns
     DO II = 1, matA%process_grid%number_of_blocks_rows
        CALL PairwiseMultiplyMatrix(matA%local_data_r(II,JJ), matB%local_data_r(II,JJ), &
             & matC%local_data_r(II,JJ))
     END DO
  END DO
  !$omp end do
  !$omp end parallel
    END IF
  END SUBROUTINE PairwiseMultiplyMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the norm of a distributed sparse matrix along the rows.
  FUNCTION MatrixNorm_ps(this) RESULT(norm_value)
    !> The matrix to compute the norm of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The norm value of the full distributed sparse matrix.
    REAL(NTREAL) :: norm_value
    !! Local Data
    REAL(NTREAL), DIMENSION(:), ALLOCATABLE :: local_norm
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data_c)
  ALLOCATE(local_norm(merged_local_data_c%columns))

  !! Sum Along Columns
  CALL MatrixColumnNorm(merged_local_data_c,local_norm)
  CALL MPI_Allreduce(MPI_IN_PLACE,local_norm,SIZE(local_norm), &
       & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

  !! Find Max Value Amonst Columns
  norm_value = MAXVAL(local_norm)
  CALL MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPINTREAL,MPI_MAX, &
       & this%process_grid%row_comm, ierr)

  CALL DestructMatrix(merged_local_data_c)
  DEALLOCATE(local_norm)
    ELSE
  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data_r)
  ALLOCATE(local_norm(merged_local_data_r%columns))

  !! Sum Along Columns
  CALL MatrixColumnNorm(merged_local_data_r,local_norm)
  CALL MPI_Allreduce(MPI_IN_PLACE,local_norm,SIZE(local_norm), &
       & MPINTREAL, MPI_SUM, this%process_grid%column_comm, ierr)

  !! Find Max Value Amonst Columns
  norm_value = MAXVAL(local_norm)
  CALL MPI_Allreduce(MPI_IN_PLACE,norm_value,1,MPINTREAL,MPI_MAX, &
       & this%process_grid%row_comm, ierr)

  CALL DestructMatrix(merged_local_data_r)
  DEALLOCATE(local_norm)
    END IF
  END FUNCTION MatrixNorm_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !> Note that a dot product is the sum of elementwise multiplication, not
  !> traditional matrix multiplication.
  SUBROUTINE DotMatrix_psr(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> The dot product.
    REAL(NTREAL), INTENT(OUT) :: product

  !! Local Data
  TYPE(Matrix_ps) :: matAH
  TYPE(Matrix_ps) :: matC

  IF (matA%is_complex) THEN
     CALL CopyMatrix(matA, matAH)
     CALL ConjugateMatrix(matAH)
     CALL PairwiseMultiplyMatrix(matAH, matB, matC)
     CALL DestructMatrix(matAH)
  ELSE
     CALL PairwiseMultiplyMatrix(matA,matB,matC)
  END IF

  CALL MatrixGrandSum(matC, product)
  CALL DestructMatrix(matC)
  END SUBROUTINE DotMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> product = dot(Matrix A,Matrix B)
  !> Note that a dot product is the sum of elementwise multiplication, not
  !> traditional matrix multiplication.
  SUBROUTINE DotMatrix_psc(matA, matB, product)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(IN)  :: matB
    !> The dot product.
    COMPLEX(NTCOMPLEX), INTENT(OUT) :: product

  !! Local Data
  TYPE(Matrix_ps) :: matAH
  TYPE(Matrix_ps) :: matC

  IF (matA%is_complex) THEN
     CALL CopyMatrix(matA, matAH)
     CALL ConjugateMatrix(matAH)
     CALL PairwiseMultiplyMatrix(matAH, matB, matC)
     CALL DestructMatrix(matAH)
  ELSE
     CALL PairwiseMultiplyMatrix(matA,matB,matC)
  END IF

  CALL MatrixGrandSum(matC, product)
  CALL DestructMatrix(matC)
  END SUBROUTINE DotMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Matrix B = alpha*Matrix A + Matrix B (AXPY)
  !> This will utilize the sparse vector increment routine.
  RECURSIVE SUBROUTINE IncrementMatrix_ps(matA, matB, alpha_in, threshold_in)
    !> Matrix A.
    TYPE(Matrix_ps), INTENT(IN)  :: matA
    !> Matrix B.
    TYPE(Matrix_ps), INTENT(INOUT)  :: matB
    !> Multiplier (default= 1.0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: alpha_in
    !> For flushing values to zero (default=0).
    REAL(NTREAL), OPTIONAL, INTENT(IN) :: threshold_in
    !! Local Data
    TYPE(Matrix_ps) :: converted_matrix
    REAL(NTREAL) :: alpha
    REAL(NTREAL) :: threshold
    INTEGER :: II, JJ

    !! Optional Parameters
    IF (.NOT. PRESENT(alpha_in)) THEN
       alpha = 1.0d+0
    ELSE
       alpha = alpha_in
    END IF
    IF (.NOT. PRESENT(threshold_in)) THEN
       threshold = 0
    ELSE
       threshold = threshold_in
    END IF

    IF (matA%is_complex .AND. .NOT. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matB, converted_matrix)
       CALL IncrementMatrix(matA, converted_matrix, alpha, threshold)
       CALL CopyMatrix(converted_matrix, matB)
    ELSE IF (.NOT. matA%is_complex .AND. matB%is_complex) THEN
       CALL ConvertMatrixToComplex(matA, converted_matrix)
       CALL IncrementMatrix(converted_matrix, matB, alpha, threshold)
    ELSE IF (matA%is_complex .AND. matB%is_complex) THEN
  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, matA%process_grid%number_of_blocks_columns
     DO II = 1, matA%process_grid%number_of_blocks_rows
        CALL IncrementMatrix(matA%local_data_c(II,JJ), matB%local_data_c(II,JJ), alpha, &
             & threshold)
     END DO
  END DO
  !$omp end do
  !$omp end parallel
    ELSE
  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, matA%process_grid%number_of_blocks_columns
     DO II = 1, matA%process_grid%number_of_blocks_rows
        CALL IncrementMatrix(matA%local_data_r(II,JJ), matB%local_data_r(II,JJ), alpha, &
             & threshold)
     END DO
  END DO
  !$omp end do
  !$omp end parallel
    END IF

    CALL DestructMatrix(converted_matrix)

  END SUBROUTINE IncrementMatrix_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  SUBROUTINE ScaleMatrix_psr(this, constant)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    REAL(NTREAL), INTENT(IN) :: constant
    !! Local Data
    INTEGER :: II, JJ

    IF (this%is_complex) THEN
  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL ScaleMatrix(this%local_data_c(II,JJ),constant)
     END DO
  END DO
  !$omp end do
  !$omp end parallel
    ELSE
  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL ScaleMatrix(this%local_data_r(II,JJ),constant)
     END DO
  END DO
  !$omp end do
  !$omp end parallel
    END IF

  END SUBROUTINE ScaleMatrix_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Will scale a distributed sparse matrix by a constant.
  RECURSIVE SUBROUTINE ScaleMatrix_psc(this, constant)
    !> Matrix to scale.
    TYPE(Matrix_ps), INTENT(INOUT) :: this
    !> A constant scale factor.
    COMPLEX(NTCOMPLEX), INTENT(IN) :: constant
    !! Local Data
    TYPE(Matrix_ps) :: this_c
    INTEGER :: II, JJ

    IF (this%is_complex) THEN
  !$omp parallel
  !$omp do collapse(2)
  DO JJ = 1, this%process_grid%number_of_blocks_columns
     DO II = 1, this%process_grid%number_of_blocks_rows
        CALL ScaleMatrix(this%local_data_c(II,JJ),constant)
     END DO
  END DO
  !$omp end do
  !$omp end parallel
    ELSE
       CALL ConvertMatrixToComplex(this, this_c)
       CALL ScaleMatrix_psc(this_c, constant)
       CALL CopyMatrix(this_c, this)
       CALL DestructMatrix(this_c)
    END IF

  END SUBROUTINE ScaleMatrix_psc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Compute the trace of the matrix.
  SUBROUTINE MatrixTrace_psr(this, trace_value)
    !> The matrix to compute the trace of.
    TYPE(Matrix_ps), INTENT(IN) :: this
    !> The trace value of the full distributed sparse matrix.
    REAL(NTREAL), INTENT(OUT) :: trace_value
    !! Local data
    TYPE(TripletList_r) :: triplet_list_r
    TYPE(TripletList_c) :: triplet_list_c
    !! Counters/Temporary
    INTEGER :: counter
    TYPE(Matrix_lsr) :: merged_local_data_r
    TYPE(Matrix_lsc) :: merged_local_data_c
    INTEGER :: ierr

    IF (this%is_complex) THEN
  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data_c)

  !! Compute The Local Contribution
  trace_value = 0
  CALL MatrixToTripletList(merged_local_data_c, triplet_list_c)
  DO counter = 1, triplet_list_c%CurrentSize
     IF (this%start_row + triplet_list_c%data(counter)%index_row .EQ. &
          & this%start_column + triplet_list_c%data(counter)%index_column) THEN
        trace_value = trace_value + &
             & REAL(triplet_list_c%data(counter)%point_value, NTREAL)
     END IF
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)

  CALL DestructMatrix(merged_local_data_c)
    ELSE
  !! Merge all the local data
  CALL MergeMatrixLocalBlocks(this, merged_local_data_r)

  !! Compute The Local Contribution
  trace_value = 0
  CALL MatrixToTripletList(merged_local_data_r, triplet_list_r)
  DO counter = 1, triplet_list_r%CurrentSize
     IF (this%start_row + triplet_list_r%data(counter)%index_row .EQ. &
          & this%start_column + triplet_list_r%data(counter)%index_column) THEN
        trace_value = trace_value + &
             & REAL(triplet_list_r%data(counter)%point_value, NTREAL)
     END IF
  END DO

  !! Sum Among Process Slice
  CALL MPI_Allreduce(MPI_IN_PLACE, trace_value, 1, MPINTREAL, &
       & MPI_SUM, this%process_grid%within_slice_comm, ierr)

  CALL DestructMatrix(merged_local_data_r)
    END IF
  END SUBROUTINE MatrixTrace_psr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE PSMatrixAlgebraModule
