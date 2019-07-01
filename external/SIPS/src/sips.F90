!  This is a module for simplifying and enhancing the usage of spectrum slicing
!  eigensolver and other eigensolvers provided in SLEPc.
!
!  Copyright (c) 2016-2017, Murat Keceli
!  All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!  POSSIBILITY OF SUCH DAMAGE.
!
MODULE M_SIPS

#include <slepc/finclude/slepc.h>

USE slepceps

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(0.0d0)
    PetscInt           :: nrank
    PetscInt           :: rank
    PetscInt           :: istart
    PetscInt           :: iend
    EPS                :: eps
    Mat                :: math
    Mat                :: mats
    Mat                :: newmath
    Mat                :: submatA
    MPI_Comm           :: matcomm
    PetscErrorCode     :: ierr
    PetscBool          :: flag_update_eps
    PetscBool          :: flag_update_ham
    PetscBool          :: flag_load_ham_ovlp

CONTAINS

    SUBROUTINE sips_initialize()

        IMPLICIT NONE

        CALL SlepcInitialize(PETSC_NULL_CHARACTER,ierr)
        CHKERRQ(ierr)

        CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
        CALL MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)

        flag_update_eps    = .true.
        flag_update_ham    = .true.
        flag_load_ham_ovlp = .true.

    END SUBROUTINE

    SUBROUTINE sips_finalize()

        IMPLICIT NONE

        CALL EPSDestroy(eps,ierr)
        CHKERRQ(ierr)

        CALL MatDestroy(math,ierr)
        CHKERRQ(ierr)

        CALL MatDestroy(mats,ierr)
        CHKERRQ(ierr)

        CALL MatDestroy(newmath,ierr)
        CHKERRQ(ierr)

        CALL MatDestroy(submatA,ierr)
        CHKERRQ(ierr)

        flag_update_eps    = .true.
        flag_update_ham    = .true.
        flag_load_ham_ovlp = .true.

    END SUBROUTINE

    SUBROUTINE sips_set_eps(stdevp)

        IMPLICIT NONE

        PetscInt, INTENT(IN) :: stdevp ! 0: Generalized problem

        ST  :: st
        KSP :: ksp
        PC  :: pc

        CALL EPSCreate(PETSC_COMM_WORLD,eps,ierr)
        CHKERRQ(ierr)

        IF (stdevp == 0) THEN
            CALL EPSSetOperators(eps,math,mats,ierr)
            CHKERRQ(ierr)

            CALL EPSSetProblemType(eps,EPS_GHEP,ierr)
            CHKERRQ(ierr)
        ELSE
            CALL EPSSetOperators(eps,math,PETSC_NULL_MAT,ierr)
            CHKERRQ(ierr)

            CALL EPSSetProblemType(eps,EPS_HEP,ierr)
            CHKERRQ(ierr)
        END IF

        CALL EPSSetWhichEigenpairs(eps,EPS_ALL,ierr)
        CHKERRQ(ierr)

        CALL EPSSetType(eps,EPSKRYLOVSCHUR,ierr)
        CHKERRQ(ierr)

        CALL EPSGetST(eps,st,ierr)
        CHKERRQ(ierr)

        CALL STSetType(st,STSINVERT,ierr)
        CHKERRQ(ierr)

        CALL STSetMatStructure(st,SAME_NONZERO_PATTERN,ierr)
        CHKERRQ(ierr)

        CALL STGetKSP(st,ksp,ierr)
        CHKERRQ(ierr)

        CALL KSPSetType(ksp,KSPPREONLY,ierr)
        CHKERRQ(ierr)

        CALL KSPGetPC(ksp,pc,ierr)
        CHKERRQ(ierr)

        CALL PCSetType(pc,PCCHOLESKY,ierr)
        CHKERRQ(ierr)

        CALL EPSKrylovSchurSetDetectZeros(eps,PETSC_TRUE,ierr)
        CHKERRQ(ierr)

        CALL EPSSetPurify(eps,PETSC_FALSE,ierr)
        CHKERRQ(ierr)

        CALL PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
        CHKERRQ(ierr)

        ! Add MUMPS options (currently no better way of setting this):
        ! Turn off ScaLAPACK in inertia computation
        CALL PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                 "-mat_mumps_icntl_13 1",ierr)
        CHKERRQ(ierr)

        ! Increase workspace with a percentage (50, 100 or more)
        CALL PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                 "-mat_mumps_icntl_14 100",ierr)
        CHKERRQ(ierr)

        ! Detect null pivots (a shift is equal to an eigenvalue)
        CALL PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                 "-mat_mumps_icntl_24 1",ierr)
        CHKERRQ(ierr)

        ! Choose PT-SCOTCH as reordering method
        CALL PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                 "-mat_mumps_icntl_29 1",ierr)
        CHKERRQ(ierr)

        ! Tolerance for null pivot detection
        CALL PetscOptionsInsertString(PETSC_NULL_OPTIONS,&
                 "-mat_mumps_cntl_3 1e-10",ierr)
        CHKERRQ(ierr)

    END SUBROUTINE

    ! Updates SLEPc eps object with new subcomm matrix A.
    SUBROUTINE sips_update_eps(nsub)

        IMPLICIT NONE

        PetscInt, INTENT(IN) :: nsub

        IF (flag_update_eps) THEN
            CALL EPSKrylovSchurGetSubcommMats(eps,submatA,PETSC_NULL_MAT,ierr)
            CHKERRQ(ierr)

            CALL PetscObjectGetComm(submatA,matcomm,ierr)
            CHKERRQ(ierr)

            CALL MatCreateRedundantMatrix(newmath,nsub,matcomm,&
                     MAT_INITIAL_MATRIX,submatA,ierr)
            CHKERRQ(ierr)
        ELSE
            CALL MatCreateRedundantMatrix(newmath,nsub,matcomm,&
                     MAT_REUSE_MATRIX,submatA,ierr)
            CHKERRQ(ierr)
        END IF

        CALL EPSKrylovSchurUpdateSubcommMats(eps,0.0_dp,1.0_dp,submatA,1.0_dp,&
                 0.0_dp,PETSC_NULL_MAT,SAME_NONZERO_PATTERN,PETSC_FALSE,ierr)
        CHKERRQ(ierr)

        flag_update_eps = .false.

    END SUBROUTINE

    SUBROUTINE sips_solve_eps(nconv)

        IMPLICIT NONE

        PetscInt, INTENT(OUT) :: nconv

        CALL EPSSetUp(eps,ierr)
        CHKERRQ(ierr)

        CALL EPSSolve(eps,ierr)
        CHKERRQ(ierr)

        CALL EPSGetConverged(eps,nconv,ierr)
        CHKERRQ(ierr)

    END SUBROUTINE

    SUBROUTINE sips_get_eigenvalues(nev,evals)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)  :: nev
        PetscReal, INTENT(OUT) :: evals(nev)

        PetscReal :: eig
        PetscInt  :: i

        DO i = 1,nev
            CALL EPSGetEigenvalue(eps,i-1,eig,PETSC_NULL_SCALAR,ierr)
            CHKERRQ(ierr)

            evals(i) = eig
        END DO

    END SUBROUTINE

    SUBROUTINE sips_load_ham(ncol_g,ncol_l,nnz_l,col_idx,row_ptr,ham_val)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: ncol_g            ! Global size of H
        PetscInt,  INTENT(IN)    :: ncol_l            ! Local size
        PetscInt,  INTENT(IN)    :: nnz_l             ! Local non-zeros
        PetscInt,  INTENT(INOUT) :: col_idx(nnz_l)    ! Column index
        PetscInt,  INTENT(INOUT) :: row_ptr(ncol_l+1) ! Row pointer
        PetscReal, INTENT(IN)    :: ham_val(nnz_l)    ! Non-zero values

        PetscInt :: i
        PetscInt :: k
        PetscInt :: nnz_this_row
        PetscInt :: which_row(1)

        ! Index conversion
        col_idx = col_idx-1
        row_ptr = row_ptr-1

        ! Create PETSC matrix
        CALL MatCreate(PETSC_COMM_WORLD,math,ierr)
        CHKERRQ(ierr)

        CALL MatSetSizes(math,ncol_l,ncol_l,ncol_g,ncol_g,ierr)
        CHKERRQ(ierr)

        CALL MatSetFromOptions(math,ierr)
        CHKERRQ(ierr)

        CALL MatSetUp(math,ierr)
        CHKERRQ(ierr)

        CALL MatGetOwnershipRange(math,istart,iend,ierr)
        CHKERRQ(ierr)

        ! Set values
        k = 1
        DO i = istart,iend-1
            nnz_this_row = row_ptr(k+1)-row_ptr(k)
            which_row(1) = i

            CALL MatSetValues(math,1,which_row,nnz_this_row,&
                     col_idx(row_ptr(k)+1:row_ptr(k+1)),&
                     ham_val(row_ptr(k)+1:row_ptr(k+1)),INSERT_VALUES,ierr)
            CHKERRQ(ierr)

            k = k+1
        END DO

        ! Assemble matrix
        CALL MatAssemblyBegin(math,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

        ! Index conversion
        col_idx = col_idx+1
        row_ptr = row_ptr+1

        ! Assemble matrix
        CALL MatAssemblyEnd(math,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

    END SUBROUTINE

    SUBROUTINE sips_update_ham(ncol_g,ncol_l,nnz_l,col_idx,row_ptr,ham_val)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: ncol_g            ! Global size of H
        PetscInt,  INTENT(IN)    :: ncol_l            ! Local size
        PetscInt,  INTENT(IN)    :: nnz_l             ! Local non-zeros
        PetscInt,  INTENT(INOUT) :: col_idx(nnz_l)    ! Column index
        PetscInt,  INTENT(INOUT) :: row_ptr(ncol_l+1) ! Row pointer
        PetscReal, INTENT(IN)    :: ham_val(nnz_l)    ! Non-zero values

        PetscInt :: i
        PetscInt :: k
        PetscInt :: nnz_this_row
        PetscInt :: which_row(1)

        ! Index conversion
        col_idx = col_idx-1
        row_ptr = row_ptr-1

        ! Create a new ham matrix from the old one
        IF (flag_update_ham) THEN
            CALL MatDuplicate(math,MAT_COPY_VALUES,newmath,ierr)
            CHKERRQ(ierr)

            CALL MatSetUp(newmath,ierr)
            CHKERRQ(ierr)
        END IF

        ! Set values
        k = 1
        DO i = istart,iend-1
            nnz_this_row = row_ptr(k+1)-row_ptr(k)
            which_row(1) = i

            CALL MatSetValues(newmath,1,which_row,nnz_this_row,&
                     col_idx(row_ptr(k)+1:row_ptr(k+1)),&
                     ham_val(row_ptr(k)+1:row_ptr(k+1)),INSERT_VALUES,ierr)
            CHKERRQ(ierr)

            k = k+1
        END DO

        ! Assemble matrix
        CALL MatAssemblyBegin(newmath,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

        ! Index conversion
        col_idx = col_idx+1
        row_ptr = row_ptr+1

        ! Assemble matrix
        CALL MatAssemblyEnd(newmath,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

        flag_update_ham = .false.

    END SUBROUTINE

    SUBROUTINE sips_load_ham_ovlp(ncol_g,ncol_l,nnz_l,col_idx,row_ptr,ham_val,&
                   ovlp_val)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: ncol_g            ! Global size of H
        PetscInt,  INTENT(IN)    :: ncol_l            ! Local size
        PetscInt,  INTENT(IN)    :: nnz_l             ! Local non-zeros
        PetscInt,  INTENT(INOUT) :: col_idx(nnz_l)    ! Column index
        PetscInt,  INTENT(INOUT) :: row_ptr(ncol_l+1) ! Row pointer
        PetscReal, INTENT(IN)    :: ham_val(nnz_l)    ! Non-zero values
        PetscReal, INTENT(IN)    :: ovlp_val(nnz_l)   ! Non-zero values

        PetscInt :: i
        PetscInt :: k
        PetscInt :: nnz_this_row
        PetscInt :: which_row(1)

        ! Index conversion
        col_idx = col_idx-1
        row_ptr = row_ptr-1

        ! Create PETSc matrix
        IF (flag_load_ham_ovlp) THEN
            CALL MatCreate(PETSC_COMM_WORLD,mats,ierr)
            CHKERRQ(ierr)

            CALL MatSetSizes(mats,ncol_l,ncol_l,ncol_g,ncol_g,ierr)
            CHKERRQ(ierr)

            CALL MatSetFromOptions(mats,ierr)
            CHKERRQ(ierr)

            CALL MatSetUp(mats,ierr)
            CHKERRQ(ierr)
        END If

        CALL MatGetOwnershipRange(mats,istart,iend,ierr)
        CHKERRQ(ierr)

        ! Set values
        k = 1
        DO i = istart,iend-1
            nnz_this_row = row_ptr(k+1)-row_ptr(k)
            which_row(1) = i

            CALL MatSetValues(mats,1,which_row,nnz_this_row,&
                     col_idx(row_ptr(k)+1:row_ptr(k+1)),&
                     ovlp_val(row_ptr(k)+1:row_ptr(k+1)),INSERT_VALUES,ierr)
            CHKERRQ(ierr)

            k = k+1
        END DO

        ! Assemble matrix
        CALL MatAssemblyBegin(mats,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

        ! Assemble matrix
        CALL MatAssemblyEnd(mats,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

        ! Create ham matrix from ovlp
        IF (flag_load_ham_ovlp) THEN
            CALL MatDuplicate(mats,MAT_COPY_VALUES,math,ierr)
            CHKERRQ(ierr)

            CALL MatSetUp(math,ierr)
            CHKERRQ(ierr)
        END IF

        ! Set values
        k = 1
        DO i = istart,iend-1
            nnz_this_row = row_ptr(k+1)-row_ptr(k)
            which_row(1) = i

            CALL MatSetValues(math,1,which_row,nnz_this_row,&
                     col_idx(row_ptr(k)+1:row_ptr(k+1)),&
                     ham_val(row_ptr(k)+1:row_ptr(k+1)),INSERT_VALUES,ierr)
            CHKERRQ(ierr)

            k = k+1
        END DO

        ! Assemble matrix
        CALL MatAssemblyBegin(math,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

        ! Index conversion
        col_idx = col_idx+1
        row_ptr = row_ptr+1

        CALL MatAssemblyEnd(math,MAT_FINAL_ASSEMBLY,ierr)
        CHKERRQ(ierr)

        flag_load_ham_ovlp = .false.

    END SUBROUTINE

    SUBROUTINE sips_set_slices(nsub,subs)

        IMPLICIT NONE

        PetscInt,  INTENT(IN) :: nsub
        PetscReal, INTENT(IN) :: subs(nsub+1)

        CALL EPSSetInterval(eps,subs(1),subs(nsub+1),ierr)
        CHKERRQ(ierr)

        CALL EPSKrylovSchurSetPartitions(eps,nsub,ierr)
        CHKERRQ(ierr)

        CALL EPSKrylovSchurSetSubintervals(eps,subs,ierr)
        CHKERRQ(ierr)

    END SUBROUTINE

    SUBROUTINE sips_get_inertias(nsub,subs,inertias)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: nsub
        PetscReal, INTENT(INOUT) :: subs(nsub+1)
        PetscInt,  INTENT(OUT)   :: inertias(nsub+1)

        PetscInt :: nshift

        CALL sips_set_slices(nsub,subs)

        CALL EPSSetUp(eps,ierr)
        CHKERRQ(ierr)

        CALL EPSKrylovSchurGetInertias(eps,nshift,subs,inertias,ierr)
        CHKERRQ(ierr)

    END SUBROUTINE

    SUBROUTINE sips_get_eigenvectors(nev,lrow,evecs)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: nev
        PetscInt,  INTENT(IN)    :: lrow
        PetscReal, INTENT(INOUT) :: evecs(lrow,nev)

        Vec                :: xr
        PetscInt           :: i
        PetscReal, POINTER :: vec_tmp(:)

        CALL MatCreateVecs(math,xr,PETSC_NULL_VEC,ierr)
        CHKERRQ(ierr)

        DO i = 1,nev
            CALL EPSGetEigenvector(eps,i-1,xr,PETSC_NULL_VEC,ierr)
            CHKERRQ(ierr)

            CALL VecGetArrayReadF90(xr,vec_tmp,ierr)
            CHKERRQ(ierr)

            evecs(1:lrow,i) = vec_tmp(1:lrow)

            CALL VecRestoreArrayReadF90(xr,vec_tmp,ierr)
            CHKERRQ(ierr)
        END DO

        CALL VecDestroy(xr,ierr)
        CHKERRQ(ierr)

    END SUBROUTINE

    SUBROUTINE sips_get_slices(algr,nev,nsub,buf,subbuf,evals,subs)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)  :: algr         ! Slicing algorithm
        PetscInt,  INTENT(IN)  :: nev          ! Number of eigenvalues
        PetscInt,  INTENT(IN)  :: nsub         ! Number of slices
        PetscReal, INTENT(IN)  :: buf          ! Buffer for global interval
        PetscReal, INTENT(IN)  :: subbuf       ! Buffer for subintervals
        PetscReal, INTENT(IN)  :: evals(nev)   ! Eigenvalues
        PetscReal, INTENT(OUT) :: subs(nsub+1) ! Subintervals

        PetscInt :: i
        PetscInt :: ngap
        PetscInt :: newid
        PetscInt :: ids(nsub)

        PetscReal, PARAMETER :: dtol = 0.00001d0 ! Evals within this value
                                                 ! always in one cluster
        PetscReal, PARAMETER :: gap = 0.005d0    ! Evals with a gap > this value
                                                 ! separated into two clusters

        IF (nsub == 1) THEN
            subs(1) = evals(1)-buf
            subs(2) = evals(nev)+buf

            RETURN
        END IF

        SELECT CASE (algr)
        CASE (0) ! Equally spaced
            subs(1:nsub+1) = get_linspace(evals(1)-buf,evals(nev)+buf,nsub+1)
        CASE (2) ! Equally populated
            ids            = 0
            ids(1:nsub)    = get_cluster_ids_uniform_population(nev,evals,nsub)
            subs(1:nsub+1) = get_subs_from_cluster_ids(buf,subbuf,nsub,&
                                 ids(1:nsub),nev,evals)
        CASE (4) ! Gap sensitive equally populated
            ids  = 0
            ngap = get_number_of_gaps(nev,evals,gap)

            IF (ngap > INT(nsub/2)) THEN
                ngap = INT(nsub/2)
            END IF

            ids(1:nsub-ngap) = get_cluster_ids_uniform_population(nev,evals,&
                                   nsub-ngap)

            DO i = 1,ngap
                newid            = get_gap_id(nev,nsub-ngap+i-1,evals,&
                                       ids(1:nsub-ngap+i-1))
                ids(nsub-ngap+i) = newid
            END DO

            CALL sort_array(nsub,ids)
            CALL fix_cluster_ids(nev,nsub,evals,ids(1:nsub),dtol)

            subs(1:nsub+1) = get_subs_from_cluster_ids(buf,subbuf,nsub,&
                                 ids(1:nsub),nev,evals)
        CASE DEFAULT
            STOP "SLEPc-SIPs: Unknown slicing method."
        END SELECT

    END SUBROUTINE

    FUNCTION get_linspace(first,last,n) RESULT(x)

        IMPLICIT NONE

        PetscReal, INTENT(IN) :: first
        PetscReal, INTENT(IN) :: last
        PetscInt,  INTENT(IN) :: n
        PetscReal             :: x(n)

        PetscReal :: step
        PetscInt  :: i

        step = (last-first)/(n-1)

        DO i = 1,n
            x(i) = first+(i-1)*step
        END DO

    END FUNCTION

    ! Given a 1D array of eigenvalues, return the indices of k clusters, with
    ! uniform population distribution.
    ! Indices correspond to lowest number in the cluster.
    FUNCTION get_cluster_ids_uniform_population(nev,evals,k) RESULT(ids)

        IMPLICIT NONE

        PetscInt,  INTENT(IN) :: nev
        PetscReal, INTENT(IN) :: evals(nev)
        PetscInt,  INTENT(IN) :: k

        PetscInt :: ids(k)
        PetscInt :: d
        PetscInt :: i
        PetscInt :: j
        PetscInt :: r

        d      = nev/k
        r      = MOD(nev,k)
        ids(1) = 1

        j = 1
        DO i = 1,k-1
            IF (j < r) THEN
                ids(i+1) = i*d+j+1
                j        = j+1
            ELSE
                ids(i+1) = i*d+r+1
            END IF
        END DO

    END FUNCTION

    ! Returns subintervals for given cluster ids and evals.
    FUNCTION get_subs_from_cluster_ids(buffer,subbuffer,nsub,ids,nev,evals) &
                 result(subs)

        IMPLICIT NONE

        PetscReal, INTENT(IN) :: buffer
        PetscReal, INTENT(IN) :: subbuffer
        PetscInt,  INTENT(IN) :: nsub
        PetscInt,  INTENT(IN) :: ids(nsub)
        PetscInt,  INTENT(IN) :: nev
        PetscReal, INTENT(IN) :: evals(nev)
        PetscReal             :: subs(nsub+1)

        PetscInt :: i

        subs(1)      = evals(1)-buffer
        subs(nsub+1) = evals(nev)+buffer

        DO i = 2,nsub
            subs(i) = evals(ids(i))-subbuffer
        END DO

    END FUNCTION

    ! Returns the number of gaps for given evals and a gap size.
    FUNCTION get_number_of_gaps(nev,evals,gap) RESULT(ngap)

        IMPLICIT NONE

        PetscInt,  INTENT(IN) :: nev
        PetscReal, INTENT(IN) :: evals(nev)
        PetscReal, INTENT(IN) :: gap
        PetscInt              :: ngap

        PetscInt  :: i
        PetscReal :: dist

        ngap = 0

        DO i = 1,nev-1
            dist = evals(i+1)-evals(i)

            IF (dist > gap) THEN
                ngap = ngap+1
            END IF
        END DO

    END FUNCTION

    FUNCTION get_gap_id(nev,nid,evals,ids) RESULT(id)

        IMPLICIT NONE

        PetscInt,  INTENT(IN) :: nev
        PetscInt,  INTENT(IN) :: nid
        PetscReal, INTENT(IN) :: evals(nev)
        PetscInt,  INTENT(IN) :: ids(nid)
        PetscInt              :: id

        PetscInt :: i
        PetscInt :: idid

        PetscReal, ALLOCATABLE :: dists(:)

        ALLOCATE(dists(nev-1))

        DO i = 1,nev-1
            dists(i) = evals(i+1)-evals(i)
        END DO

        id = MAXLOC(dists,DIM=1)+1

        DO WHILE (.TRUE.)
            idid = get_index(nid,ids,id)

            IF (idid > 0) THEN
                id = MAXLOC(dists,DIM=1,MASK=(dists<dists(id-1)))+1
            ELSE
                EXIT
            END IF
        END DO

    END FUNCTION

    ! Retuns the index of a target value in a given array.
    ! Returns the index of the first match.
    ! Returns 0 if not found.
    FUNCTION get_index(n,x,val) RESULT(id)

        IMPLICIT NONE

        PetscInt, INTENT(IN) :: n
        PetscInt, INTENT(IN) :: x(n)
        PetscInt, INTENT(IN) :: val
        PetscInt             :: id

        PetscInt :: i

        id = 0

        DO i = 1,n
            IF (x(i) == val) THEN
                id = i
                EXIT
            END IF
        END DO

    END FUNCTION

    ! Sorts the given array in ascending order.
    SUBROUTINE sort_array(n,x)

        IMPLICIT NONE

        PetscInt, INTENT(IN)    :: n
        PetscInt, INTENT(INOUT) :: x(n)

        PetscInt :: i
        PetscInt :: tmp
        PetscInt :: minid
        PetscInt :: curmin

        curmin = MINVAL(x,DIM=1)

        DO i = 1,n
            minid    = MINLOC(x(i:n),DIM=1,MASK=(x(i:n)>=curmin))+i-1
            curmin   = x(minid)
            tmp      = x(i)
            x(i)     = x(minid)
            x(minid) = tmp
        END DO

    END SUBROUTINE

    ! If any almost degenerate (within dtol) evals are in two clusters, smaller
    ! evals are merged into the cluster of the larger evals.
    SUBROUTINE fix_cluster_ids(nev,nid,evals,ids,dtol)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: nev
        PetscInt,  INTENT(IN)    :: nid
        PetscReal, INTENT(IN)    :: evals(nev)
        PetscInt,  INTENT(INOUT) :: ids(nid)
        PetscReal, INTENT(IN)    :: dtol

        PetscInt :: i
        PetscInt :: oldids(nid)

        oldids = 0

        DO WHILE (ANY(ids-oldids /= 0))
            oldids = ids

            DO i = 1,nid-1
                IF (evals(ids(i+1))-evals(ids(i+1)-1) < dtol) THEN
                    ids(i+1) = ids(i+1)-1
                END IF
            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE inertias_to_eigenvalues(nsub,nev,subs,inertias,evals)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)  :: nsub
        PetscInt,  INTENT(IN)  :: nev
        PetscReal, INTENT(IN)  :: subs(nsub+1)
        PetscInt,  INTENT(IN)  :: inertias(nsub+1)
        PetscReal, INTENT(OUT) :: evals(nev)

        PetscInt  :: i
        PetscInt  :: j
        PetscInt  :: k
        PetscInt  :: ne
        PetscInt  :: last
        PetscReal :: sw
        PetscReal :: step

        k    = 0
        last = nsub

        DO i = 1,nsub
            ne = inertias(i+1)-inertias(i)
            sw = subs(i+1)-subs(i)

            IF (ne == 0) THEN
                CYCLE
            ELSE IF (ne == 1) THEN
                k        = k+1
                evals(k) = subs(i)
            ELSE
                step = sw/REAL(ne)

                DO j = 1,ne
                    k = k+1

                    IF (k >= nev) THEN
                        last = i
                        EXIT
                    END IF

                    evals(k) = subs(i)+(j-1)*step
                END DO
            END IF
            IF (k >= nev) THEN
                last = i
                EXIT
            END IF
        END DO

        evals(nev) = subs(last+1)

    END SUBROUTINE

    SUBROUTINE sips_get_slices_from_inertias(nev,nsub,inertias,subs)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: nev              ! Number of eigenvalues
        PetscInt,  INTENT(IN)    :: nsub             ! Number of slices
        PetscInt,  INTENT(IN)    :: inertias(nsub+1) ! Inertias
        PetscReal, INTENT(INOUT) :: subs(nsub+1)     ! Subintervals

        PetscReal :: evals(nev)

        CALL inertias_to_eigenvalues(nsub,nev,subs,inertias,evals)

        CALL sips_get_slices(2,nev,nsub,0.d0,0.d0,evals,subs)

    END SUBROUTINE

    SUBROUTINE sips_get_dm(ncol_l,nnz_l,col_idx,row_ptr,nev,occ,dm_val)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: ncol_l            ! Local size
        PetscInt,  INTENT(IN)    :: nnz_l             ! Local non-zeros
        PetscInt,  INTENT(INOUT) :: col_idx(nnz_l)    ! Column index
        PetscInt,  INTENT(INOUT) :: row_ptr(ncol_l+1) ! Row pointer
        PetscInt,  INTENT(IN)    :: nev               ! Number of eigenvalues
        PetscReal, INTENT(IN)    :: occ(nev)          ! Occupation numbers
        PetscReal, INTENT(OUT)   :: dm_val(nnz_l)     ! Non-zero values

        Vec                :: xr
        Vec                :: xrseq
        VecScatter         :: ctx
        PetscInt           :: i
        PetscInt           :: j
        PetscInt           :: k
        PetscInt           :: i_val
        PetscInt           :: nnz_this_row
        PetscReal          :: tmp
        PetscReal, POINTER :: vec_tmp(:)

        ! Index conversion
        col_idx = col_idx-1
        row_ptr = row_ptr-1

        CALL MatCreateVecs(math,xr,PETSC_NULL_VEC,ierr)
        CHKERRQ(ierr)

        CALL VecScatterCreateToAll(xr,ctx,xrseq,ierr)
        CHKERRQ(ierr)

        dm_val = 0.0_dp

        DO i = 1,nev
            CALL EPSGetEigenvector(eps,i-1,xr,PETSC_NULL_VEC,ierr)
            CHKERRQ(ierr)

            CALL VecScatterBegin(ctx,xr,xrseq,INSERT_VALUES,SCATTER_FORWARD,&
                     ierr)
            CHKERRQ(ierr)

            CALL VecScatterEnd(ctx,xr,xrseq,INSERT_VALUES,SCATTER_FORWARD,ierr)
            CHKERRQ(ierr)

            CALL VecGetArrayReadF90(xrseq,vec_tmp,ierr)
            CHKERRQ(ierr)

            i_val = 0

            DO j = istart+1,iend
                nnz_this_row = row_ptr(j-istart+1)-row_ptr(j-istart)

                DO k = 1,nnz_this_row
                    i_val         = i_val+1
                    tmp           = occ(i)*vec_tmp(j)*vec_tmp(col_idx(i_val)+1)
                    dm_val(i_val) = dm_val(i_val)+tmp
                END DO
            END DO

            CALL VecRestoreArrayReadF90(xrseq,vec_tmp,ierr)
            CHKERRQ(ierr)
        END DO

        CALL VecDestroy(xr,ierr)
        CHKERRQ(ierr)

        CALL VecDestroy(xrseq,ierr)
        CHKERRQ(ierr)

        CALL VecScatterDestroy(ctx,ierr)
        CHKERRQ(ierr)

        ! Index conversion
        col_idx = col_idx+1
        row_ptr = row_ptr+1

    END SUBROUTINE

    SUBROUTINE sips_get_edm(ncol_l,nnz_l,col_idx,row_ptr,nev,occ,edm_val)

        IMPLICIT NONE

        PetscInt,  INTENT(IN)    :: ncol_l            ! Local size
        PetscInt,  INTENT(IN)    :: nnz_l             ! Local non-zeros
        PetscInt,  INTENT(INOUT) :: col_idx(nnz_l)    ! Column index
        PetscInt,  INTENT(INOUT) :: row_ptr(ncol_l+1) ! Row pointer
        PetscInt,  INTENT(IN)    :: nev               ! Number of eigenvalues
        PetscReal, INTENT(IN)    :: occ(nev)          ! Occupation numbers
        PetscReal, INTENT(OUT)   :: edm_val(nnz_l)    ! Non-zero values

        Vec                :: xr
        Vec                :: xrseq
        VecScatter         :: ctx
        PetscInt           :: i
        PetscInt           :: j
        PetscInt           :: k
        PetscInt           :: i_val
        PetscInt           :: nnz_this_row
        PetscReal          :: tmp
        PetscReal          :: eval
        PetscReal, POINTER :: vec_tmp(:)

        ! Index conversion
        col_idx = col_idx-1
        row_ptr = row_ptr-1

        CALL MatCreateVecs(math,xr,PETSC_NULL_VEC,ierr)
        CHKERRQ(ierr)

        CALL VecScatterCreateToAll(xr,ctx,xrseq,ierr)
        CHKERRQ(ierr)

        edm_val = 0.0_dp

        DO i = 1,nev
            CALL EPSGetEigenpair(eps,i-1,eval,PETSC_NULL_SCALAR,xr,&
                     PETSC_NULL_VEC,ierr)
            CHKERRQ(ierr)

            CALL VecScatterBegin(ctx,xr,xrseq,INSERT_VALUES,SCATTER_FORWARD,&
                     ierr)
            CHKERRQ(ierr)

            CALL VecScatterEnd(ctx,xr,xrseq,INSERT_VALUES,SCATTER_FORWARD,ierr)
            CHKERRQ(ierr)

            CALL VecGetArrayReadF90(xrseq,vec_tmp,ierr)
            CHKERRQ(ierr)

            i_val = 0

            DO j = istart+1,iend
                nnz_this_row = row_ptr(j-istart+1)-row_ptr(j-istart)

                DO k = 1,nnz_this_row
                    i_val          = i_val+1
                    tmp            = eval*occ(i)*vec_tmp(j)*&
                                         vec_tmp(col_idx(i_val)+1)
                    edm_val(i_val) = edm_val(i_val)+tmp
                END DO
            END DO

            CALL VecRestoreArrayReadF90(xrseq,vec_tmp,ierr)
            CHKERRQ(ierr)
        END DO

        ! Index conversion
        col_idx = col_idx+1
        row_ptr = row_ptr+1

        CALL VecDestroy(xr,ierr)
        CHKERRQ(ierr)

        CALL VecDestroy(xrseq,ierr)
        CHKERRQ(ierr)

        CALL VecScatterDestroy(ctx,ierr)
        CHKERRQ(ierr)

    END SUBROUTINE

END MODULE M_SIPS
