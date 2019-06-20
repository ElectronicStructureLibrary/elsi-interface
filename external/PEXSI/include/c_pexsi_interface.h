/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.

   Author: Lin Lin

   This file is part of PEXSI. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   (1) Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
   (2) Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
   (3) Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   You are under no obligation whatsoever to provide any bug fixes, patches, or
   upgrades to the features, functionality or performance of the source code
   ("Enhancements") to anyone; however, if you choose to make your Enhancements
   available either publicly, or directly to Lawrence Berkeley National
   Laboratory, without imposing a separate written license agreement for such
   Enhancements, then you hereby grant the following license: a non-exclusive,
   royalty-free perpetual license to install, use, modify, prepare derivative
   works, incorporate into other computer software, distribute, and sublicense
   such enhancements or derivative works thereof, in binary and source code form.
*/
/**
 * @file c_pexsi_interface.h
 * @brief Interface subroutines of %PEXSI that can be called by C.
 *
 * @date Original:      2013-01-31
 * @date Revision:      2014-03-07  Second generation interface.
 */
#ifndef _PEXSI_C_PEXSI_INTERFACE_H_
#define _PEXSI_C_PEXSI_INTERFACE_H_
#include <mpi.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"{
#endif

/**
 * @brief Read the sizes of a DistSparseMatrix in formatted form (txt)
 * for allocating memory in C.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[out] size (global) Number of rows and columns of the matrix.
 * @param[out] nnz (global) Total number of nonzeros.
 * @param[out] nnzLocal (local) Number of local nonzeros.
 * @param[out] numColLocal (local) Number of local columns.
 * @param[in]  comm (global) MPI communicator.
 */
void ReadDistSparseMatrixFormattedHeadInterface (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Comm comm );

/**
 * @brief Reading the data of a formatted DistSparseMatrix.
 *
 * This routine assumes that the arrays have been allocated outside this
 * subroutine.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[in] size (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros.
 * @param[in] nnzLocal (local) Number of local nonzeros.
 * @param[in] numColLocal (local) Number of local columns.
 * @param[out] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[out] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[out] nzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values in CSC format.
 * @param[in]  comm (global) MPI communicator.
 */
void ReadDistSparseMatrixFormattedInterface(
		char*     filename,
		int       size,
		int       nnz,
		int       nnzLocal,
		int       numColLocal,
		int*      colptrLocal,
		int*      rowindLocal,
		double*   nzvalLocal,
		MPI_Comm  comm );


/**
 * @brief Read the sizes of a DistSparseMatrix in unformatted form
 * (csc) for allocating memory in C.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[out] size (global) Number of rows and columns of the matrix.
 * @param[out] nnz (global) Total number of nonzeros.
 * @param[out] nnzLocal (local) Number of local nonzeros.
 * @param[out] numColLocal (local) Number of local columns.
 * @param[in]  comm (global) MPI communicator.
 */
void ReadDistSparseMatrixHeadInterface (
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Comm comm );

/**
 * @brief Actual reading the data of a DistSparseMatrix using MPI-IO,
 * assuming that the arrays have been allocated outside this
 * subroutine.
 *
 * This routine can be much faster than reading a DistSparseMatrix
 * sequentially, especially compared to the version using formatted
 * input @ref ReadDistSparseMatrixFormattedInterface.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[in] size (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros.
 * @param[in] nnzLocal (local) Number of local nonzeros.
 * @param[in] numColLocal (local) Number of local columns.
 * @param[out] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[out] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[out] nzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values in CSC format.
 * @param[in]  comm (global) MPI communicator.
 */
void ParaReadDistSparseMatrixInterface (
    char*     filename,
    int       size,
    int       nnz,
    int       nnzLocal,
    int       numColLocal,
    int*      colptrLocal,
    int*      rowindLocal,
    double*   nzvalLocal,
    MPI_Comm  comm );



// *********************************************************************
// The following routines belong to the second version of the interface
// *********************************************************************

/**
 * @brief A handle for holding the internal %PEXSI data structure.
 *
 * @note This handle can be also used with FORTRAN, with the `INTEGER*8`
 * data structure, or the `INTEGER(C_INTPTR_T)` structure if
 * ISO_C_BINDING is used.
 */
typedef intptr_t  PPEXSIPlan;

/**
 * @struct PPEXSIOptions
 * @brief Structure for the input parameters in DFT calculations.
 */
typedef struct {
    /**
     * @brief  NumSpin, default is 2.0
     */
    double        spin;
    /**
     * @brief  Temperature, in the same unit as H
     */
    double        temperature;
    /**
     * @brief  Spectral gap. **Note** This can be set to be 0 in most cases.
     */
    double        gap;
    /**
     * @brief  An upper bound for the spectral radius of \f$S^{-1} H\f$.
     */
    double        deltaE;
    /**
     * @brief  Number of terms in the pole expansion.
     */
    int           numPole;
    /**
     * @brief  Whether inertia counting is used at the very beginning.
     */
    int           isInertiaCount;
    /**
     * @brief  Maximum number of %PEXSI iterations after each inertia
     * counting procedure.
     */
    int           maxPEXSIIter;
    /**
     * @brief  Initial guess of lower bound for mu.
     */
    double        muMin0;
    /**
     * @brief  Initial guess of upper bound for mu.
     */
    double        muMax0;
     /**
     * @brief  Initial guess for mu (for the solver) (AG)
     */
    double        mu0;
    /**
     * @brief  Stopping criterion in terms of the chemical potential
     * for the inertia counting procedure.
     */
    double        muInertiaTolerance;
    /**
     * @brief  If the chemical potential is not in the initial interval,
     * the interval is expanded by muInertiaExpansion.
     */
    double        muInertiaExpansion;
    /**
     * @brief  Safe guard criterion in terms of the chemical potential
     * to reinvoke the inertia counting procedure.
     */
    double        muPEXSISafeGuard;
    /**
     * @brief  Stopping criterion of the %PEXSI iteration in terms of the
     * number of electrons compared to numElectronExact.
     */
    double        numElectronPEXSITolerance;
    /**
     * @brief  matrixType (global) Type of input H and S matrices.
     * - = 0   : Real symmetric (default)
     * - = 1   : General complex matrices (not implemented yet)
     */
    int           matrixType;
    /**
     * @brief  Whether to perform symbolic factorization.
     */
    int           isSymbolicFactorize;
    /**
     * @brief  Whether to construct PSelInv communication pattern.
     */
    int           isConstructCommPattern;

    /**
     * @brief  Solver used to do the factorization prior to the selected
     * inversion.
     * - = 0   : SuperLU_DIST.
     * - = 1   : symPACK (For symmetric matrices only).
     */
    int           solver;

    /**
     * @brief  Storage space used by the Selected Inversion algorithm for symmetric matrices.
     * - = 0   : Non symmetric storage.
     * - = 1   : Symmetric storage (lower memory usage).
     */
    int           symmetricStorage;


    /**
     * @brief  Ordering strategy for factorization and selected
     * inversion. When SuperLU is used:
     * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
     *   option in SuperLU_DIST).
     * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
     *   option in SuperLU_DIST).
     * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
     *   option in SuperLU_DIST).
     * When symPACK is used:
     * - = 0   : Parallel ordering using PT-SCOTCH.
     * - = 1   : Sequential ordering using SCOTCH.
     * - = 2   : Multiple minimum degree ordering.
     * - = 3   : Approximate minimum degree ordering.
     * - = 4   : Parallel ordering using PARMETIS.
     * - = 5   : Sequential ordering using METIS.
     */
    int           ordering;
    /**
     * @brief  row permutation strategy for factorization and selected
     * inversion.
     * - = 0   : No row permutation (NOROWPERM
     *   option in SuperLU_DIST).
     * - = 1   : Make diagonal entry larger than off diagonal ( LargeDiag
     *   option in SuperLU_DIST).
     */
    int           rowOrdering;
    /**
     * @brief  Number of processors for PARMETIS/PT-SCOTCH.  Only used
     * if the ordering == 0.
     */
    int           npSymbFact;
    /**
     * @brief  Matrix structure.
     * - = 0   : Unsymmetric matrix
     * - = 1   : Symmetric matrix (default).
     */
    int           symmetric;
    /**
     * @brief  Transpose.
     * - = 0   : Factor non transposed matrix (default).
     * - = 1   : Factor transposed matrix.
     */
    int           transpose;
    /**
     * @brief  The pole expansion method to be used.
     * - = 1   : Cauchy Contour Integral method used.
     * - = 2   : Moussa optimized method.
     */
    int           method;
    /**
     * @brief  The point parallelizaion of PEXSI.
     * - = 2  : Recommend two points parallelization
     */
    int           nPoints;
    /**
     * @brief  The driver version of the PEXSI.
     * - = 1.0.0 : Latest version is 1.0.0
     */
    //char         driverVersion[10] ;//= "1.0.0";
    /**
     * @brief  The level of output information.
     * - = 0   : No output.
     * - = 1   : Basic output (default)
     * - = 2   : Detailed output.
     */
    int           verbosity;
} PPEXSIOptions;


/**
 * @brief Set the default options for DFT driver.
 *
 * All default values assume the input unit (for H) is Rydberg.
 *
 * @param[in] options (global) Pointer to the options containing input
 * parameters for the driver.
 */
void PPEXSISetDefaultOptions(
    PPEXSIOptions*   options );


/**
 * @brief Initialize the %PEXSI plan.
 *
 * In %PEXSI, a matrix is generally referred to as a "pole". The
 * factorization and selected inversion procedure for a pole is computed
 * in parallel using `numProcRow * numProcCol` processors.
 *
 * When only selected inversion (PSelInv) is used, it is recommended to
 * set the mpisize of the communicator `comm` to be
 * `numProcRow * numProcCol`.
 *
 * When %PEXSI is used to evaluate a large number of inverse matrices
 * such as in the electronic structure calculation, mpisize should be
 * `numPole*numProcRow * numProcCol`, where `numPole` inverse matrices
 * can be processed in parallel.
 *
 *
 * @param[in] comm  (global) Communicator used for the entire %PEXSI procedure.  The
 * size of this communicator should be a multiple of npPerPole =
 * numProcRow * numProcCol, which is the total number of processors used
 * by each individual pole.
 * @param[in] numProcRow (global) Number of processors in the row communication group
 * for each pole.
 * @param[in] numProcCol (global) Number of processors in the column communication group
 * for each pole.
 * @param[in] outputFileIndex (local) The index for the %PEXSI output file.  For
 * instance, if this index is 1, then the corresponding processor will
 * output to the file `logPEXSI1`.  If the index is negative, then no output file is given.
 * **Note** Each processor must output to a **different** file.  By
 * default, outputFileIndex can be set as mpirank.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 *
 * @return (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 */
PPEXSIPlan PPEXSIPlanInitialize(
    MPI_Comm      comm,
    int           numProcRow,
    int           numProcCol,
    int           outputFileIndex,
    int*          info );


/**
 * @brief Load the real H and S matrices into the %PEXSI
 * internal data structure. H and S can be either real symmetric or real
 * unsymmetric.
 *
 * @note Only input from the processors associated with the first pole
 * is required. The information will be broadcast to the other
 * processors in the communicator.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] nrows (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros of H.
 * @param[in] nnzLocal (local) Number of local nonzeros of H.
 * @param[in] numColLocal (local) Number of local columns for H.
 * @param[in] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[in] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[in] HnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of H in CSC format.
 * @param[in] isSIdentity (global) Whether S is an identity matrix.
 * If so, the variable SnzvalLocal is omitted.
 * @param[in] SnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * value of S in CSC format.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSILoadRealHSMatrix(
    PPEXSIPlan    plan,
    PPEXSIOptions options,
    int           nrows,
    int           nnz,
    int           nnzLocal,
    int           numColLocal,
    int*          colptrLocal,
    int*          rowindLocal,
    double*       HnzvalLocal,
    int           isSIdentity,
    double*       SnzvalLocal,
    int*          info );

/**
 * @brief Load the complex H and S matrices into the %PEXSI internal
 * data structure. H and S can be either complex symmetric or complex unsymmetric.
 *
 * @note Only input from the processors associated with the first pole
 * is required. The information will be broadcast to the other
 * processors in the communicator.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] nrows (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros of H.
 * @param[in] nnzLocal (local) Number of local nonzeros of H.
 * @param[in] numColLocal (local) Number of local columns for H.
 * @param[in] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[in] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[in] HnzvalLocal (local) Dimension: 2*nnzLocal. Local nonzero
 * values of H in CSC format.
 * @param[in] isSIdentity (global) Whether S is an identity matrix.
 * If so, the variable SnzvalLocal is omitted.
 * @param[in] SnzvalLocal (local) Dimension: 2*nnzLocal. Local nonzero
 * value of S in CSC format.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSILoadComplexHSMatrix(
    PPEXSIPlan    plan,
    PPEXSIOptions options,
    int           nrows,
    int           nnz,
    int           nnzLocal,
    int           numColLocal,
    int*          colptrLocal,
    int*          rowindLocal,
    double*       HnzvalLocal,
    int           isSIdentity,
    double*       SnzvalLocal,
    int*          info );



/**
 * @brief Separately perform symbolic factorization to prepare
 * factorization and selected inversion for real arithmetic matrices.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISymbolicFactorizeRealSymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info );

/**
 * @brief Separately perform symbolic factorization to prepare
 * factorization and selected inversion for complex arithmetic matrices.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISymbolicFactorizeComplexSymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info );


/**
 * @brief Directly compute the negative inertia at a set of shifts for real matrices.
 *
 * This can be used as an "expert" interface for solving KSDFT with
 * user-implemented heuristics strategies.
 *
 * @note Only input from the processors associated with the first pole
 * is required. The information will be broadcast to the other
 * processors in the communicator.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] numShift (global) Number of shifts.
 * @param[in] shiftList (global) The list of shifts. Size: numShift
 * @param[out] inertiaList (global) The list of inertia counts (in
 * double precision but are of integer values). Size: numShift
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIInertiaCountRealMatrix(
    /* Input parameters */
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int               numShift,
    double*           shiftList,
    /* Output parameters */
    double*           inertiaList,
    int*              info );

/**
 * @brief Directly compute the negative inertia at a set of shifts for complex matrices.
 *
 * This can be used as an "expert" interface for solving KSDFT with
 * user-implemented heuristics strategies.
 *
 * @note Only input from the processors associated with the first pole
 * is required. The information will be broadcast to the other
 * processors in the communicator.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] numShift (global) Number of shifts.
 * @param[in] shiftList (global) The list of shifts. Size: numShift
 * @param[out] inertiaList (global) The list of inertia counts (in
 * double precision but are of integer values). Size: numShift
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIInertiaCountComplexMatrix(
    /* Input parameters */
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int               numShift,
    double*           shiftList,
    /* Output parameters */
    double*           inertiaList,
    int*              info );



/**
 * @brief Separately perform symbolic factorization to prepare
 * factorization and selected inversion for real arithmetic matrices.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] AnzvalLocal non zero values required to compute row permutation
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISymbolicFactorizeRealUnsymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    int*              info );

/**
 * @brief Separately perform symbolic factorization to prepare
 * factorization and selected inversion for complex arithmetic matrices.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] AnzvalLocal non zero values required to compute row permutation
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISymbolicFactorizeComplexUnsymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    int*              info );


/**
 * @brief Directly compute the negative inertia at a set of shifts.
 *
 * This can be used as an "expert" interface for solving KSDFT with
 * user-implemented heuristics strategies.
 *
 * @note Only input from the processors associated with the first pole
 * is required. The information will be broadcast to the other
 * processors in the communicator.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] numShift (global) Number of shifts.
 * @param[in] shiftList (global) The list of shifts. Size: numShift
 * @param[out] inertiaList (global) The list of inertia counts (in
 * double precision but are of integer values). Size: numShift
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
/*
void PPEXSIInertiaCountRealUnsymmetricMatrix(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int               numShift,
    double*           shiftList,
    double*           inertiaList,
    int*              info );
*/












/**
 * @brief Compute the density matrices and number of electrons for a
 * given chemical potential.
 *
 * This can be used as an "expert" interface for solving KSDFT with
 * user-implemented heuristics strategies.
 *
 * **Input/Output**
 *
 * The input parameter options are controlled through the structure
 * PPEXSIOptions.  The default value can be obtained through
 * PPEXSISetDefaultOptions.

 * The input H and S matrices should be given by loading functions
 * (currently it is PPEXSILoadRealSymmetricHSMatrix).  The output
 * matrices should be obtained from retrieving functions (currently it
 * is PPEXSIRetrieveRealSymmetricDFTMatrix).
 *
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] numElectronExact (global) Exact number of electrons, i.e.
 * \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * - In the case that convergence is reached within maxPEXSIIter steps, the
 * value of muPEXSI is the last mu used to achieve accuracy within
 * numElectronPEXSITolerance.
 * - In the case that convergence is not reached within maxPEXSIIter steps,
 * and the update from Newton's iteration does not exceed
 * muPEXSISafeGuard, the value of muPEXSI is the last mu plus the update
 * from Newton's iteration.
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at the last step.
 * **Note** In the case that convergence is not reached within maxPEXSIIter steps,
 * and numElectron does not correspond to the number of electrons
 * evaluated at muPEXSI.
 * @param[out] muMinInertia (global) Lower bound for mu after the last
 * inertia count procedure.
 * @param[out] muMaxInertia (global) Upper bound for mu after the last
 * inertia count procedure.
 * @param[out] numTotalInertiaIter (global) Number of total inertia
 * counting procedure.
 * @param[out] numTotalPEXSIIter (global) Number of total %PEXSI
 * evaluation procedure.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSICalculateFermiOperatorReal(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            mu,
    double            numElectronExact,
    double*           numElectronPEXSI,
    double*           numElectronDrvMuPEXSI,
    int*              info );

/**
 * @brief Compute the density matrices and number of electrons for a
 * given chemical potential.
 *
 * This can be used as an "expert" interface for solving KSDFT with
 * user-implemented heuristics strategies.
 *
 * **Input/Output**
 *
 * The input parameter options are controlled through the structure
 * PPEXSIOptions.  The default value can be obtained through
 * PPEXSISetDefaultOptions.

 * The input H and S matrices should be given by loading functions
 * (currently it is PPEXSILoadRealSymmetricHSMatrix).  The output
 * matrices should be obtained from retrieving functions (currently it
 * is PPEXSIRetrieveRealSymmetricDFTMatrix).
 *
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] numElectronExact (global) Exact number of electrons, i.e.
 * \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * - In the case that convergence is reached within maxPEXSIIter steps, the
 * value of muPEXSI is the last mu used to achieve accuracy within
 * numElectronPEXSITolerance.
 * - In the case that convergence is not reached within maxPEXSIIter steps,
 * and the update from Newton's iteration does not exceed
 * muPEXSISafeGuard, the value of muPEXSI is the last mu plus the update
 * from Newton's iteration.
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at the last step.
 * **Note** In the case that convergence is not reached within maxPEXSIIter steps,
 * and numElectron does not correspond to the number of electrons
 * evaluated at muPEXSI.
 * @param[out] muMinInertia (global) Lower bound for mu after the last
 * inertia count procedure.
 * @param[out] muMaxInertia (global) Upper bound for mu after the last
 * inertia count procedure.
 * @param[out] numTotalInertiaIter (global) Number of total inertia
 * counting procedure.
 * @param[out] numTotalPEXSIIter (global) Number of total %PEXSI
 * evaluation procedure.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSICalculateFermiOperatorReal3(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            numElectronExact,
    double*           mu,
    double*           numElectronPEXSI,
    int*              info );





/**
 * @brief Compute the density matrices and number of electrons for a
 * given chemical potential. Here H and S are complex Hermitian matrices
 * (even though S is real in electronic structure calculation)
 *
 * This can be used as an "expert" interface for solving KSDFT with
 * user-implemented heuristics strategies.
 *
 * **Input/Output**
 *
 * The input parameter options are controlled through the structure
 * PPEXSIOptions.  The default value can be obtained through
 * PPEXSISetDefaultOptions.

 * The input H and S matrices should be given by loading functions
 * (currently it is PPEXSILoadRealSymmetricHSMatrix).  The output
 * matrices should be obtained from retrieving functions (currently it
 * is PPEXSIRetrieveRealSymmetricDFTMatrix).
 *
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] numElectronExact (global) Exact number of electrons, i.e.
 * \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * - In the case that convergence is reached within maxPEXSIIter steps, the
 * value of muPEXSI is the last mu used to achieve accuracy within
 * numElectronPEXSITolerance.
 * - In the case that convergence is not reached within maxPEXSIIter steps,
 * and the update from Newton's iteration does not exceed
 * muPEXSISafeGuard, the value of muPEXSI is the last mu plus the update
 * from Newton's iteration.
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at the last step.
 * **Note** In the case that convergence is not reached within maxPEXSIIter steps,
 * and numElectron does not correspond to the number of electrons
 * evaluated at muPEXSI.
 * @param[out] muMinInertia (global) Lower bound for mu after the last
 * inertia count procedure.
 * @param[out] muMaxInertia (global) Upper bound for mu after the last
 * inertia count procedure.
 * @param[out] numTotalInertiaIter (global) Number of total inertia
 * counting procedure.
 * @param[out] numTotalPEXSIIter (global) Number of total %PEXSI
 * evaluation procedure.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSICalculateFermiOperatorComplex(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            mu,
    double            numElectronExact,
    double*           numElectronPEXSI,
    double*           numElectronDrvMuPEXSI,
    int*              info );


/**
 * @brief Simplified driver interface for computing the selected
 * elements of a real symmetric matrix.
 *
 * @note The computation is only performed using the group of processors
 * corresponding to the first pole.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] AnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of A in CSC format.
 * @param[out] AinvnzvalLocal (local) Dimension: nnzLocal.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISelInvRealSymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info );


/**
 * @brief Simplified driver interface for computing the selected
 * elements of a complex symmetric matrix.
 *
 * @note The computation is only performed using the group of processors
 * corresponding to the first pole.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] AnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of A in CSC format.
 * - Use 2 double for one complex number. This
 * ensures the compatibility with FORTRAN.
 * - Real part: AnzvalLocal[2*k]. Imag part: AnzvalLocal[2*k+1].
 * @param[out] AinvnzvalLocal (local) Dimension: 2*nnzLocal.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISelInvComplexSymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info );


/**
 * @brief Simplified driver interface for computing the selected
 * elements of a real unsymmetric matrix.
 *
 * @note The computation is only performed using the group of processors
 * corresponding to the first pole.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] AnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of A in CSC format.
 * @param[out] AinvnzvalLocal (local) Dimension: nnzLocal.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISelInvRealUnsymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info );


/**
 * @brief Simplified driver interface for computing the selected
 * elements of a complex unsymmetric matrix.
 *
 * @note The computation is only performed using the group of processors
 * corresponding to the first pole.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] AnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of A in CSC format.
 * - Use 2 double for one complex number. This
 * ensures the compatibility with FORTRAN.
 * - Real part: AnzvalLocal[2*k]. Imag part: AnzvalLocal[2*k+1].
 * @param[out] AinvnzvalLocal (local) Dimension: 2*nnzLocal.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSISelInvComplexUnsymmetricMatrix (
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*           AnzvalLocal,
    double*           AinvnzvalLocal,
    int*              info );



/**
 * @brief Simplified driver for solving Kohn-Sham DFT.
 *
 * This function contains both the inertia counting step for estimating
 * the chemical potential, and the Newton's iteration for updating the
 * chemical potential.  Heuristics are built into this routine.  Expert
 * users and developers can modify this routine to obtain better
 * heuristics.  The implementation of this function contains **all** the
 * heuristics for a DFT solver.
 *
 * The input parameter options are controlled through the structure
 * PPEXSIOptions.  The default value can be obtained through
 * PPEXSISetDefaultOptions.
 *
 *
 * **Basic strategy of the heuristics**
 *
 * - If isInertiaCount == 1, then the inertia counting procedure is
 * invoked until the chemical potential interval is refined from size
 * (muMax0-muMin0) to muInertiaTolerance, or the estimated band gap is
 * larger than muInertiaTolerance.
 * - If the change of step in Newton's iteration is larger than
 * muPEXSISafeGuard, the the inertia counting procedure is invoked
 * again, starting from (muMin0, muMax0).  If Newton's iteration fails
 * again, the subroutine returns error message with info = 1.
 * - The number of shifts in the inertia count is always automatically
 * chosen to be (mpisize / npPerPole). This minimizes the wall clock
 * time of the inertia counting procedure.
 *
 * **Complex Hermitian case**
 *
 * This file should work for both real symmetric and complex Hermitian H
 * and S matrices.  However, the complex Hermitian case require
 * asymmetric PSelInv which will be in the future work.
 *
 * **Input/Output**
 *
 * The input H and S matrices should be given by loading functions
 * (currently it is PPEXSILoadRealSymmetricHSMatrix).  The output
 * matrices should be obtained from retrieving functions (currently it
 * is PPEXSIRetrieveRealSymmetricDFTMatrix).
 *
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] numElectronExact (global) Exact number of electrons, i.e.
 * \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * - In the case that convergence is reached within maxPEXSIIter steps, the
 * value of muPEXSI is the last mu used to achieve accuracy within
 * numElectronPEXSITolerance.
 * - In the case that convergence is not reached within maxPEXSIIter steps,
 * and the update from Newton's iteration does not exceed
 * muPEXSISafeGuard, the value of muPEXSI is the last mu plus the update
 * from Newton's iteration.
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at the last step.
 * **Note** In the case that convergence is not reached within maxPEXSIIter steps,
 * and numElectron does not correspond to the number of electrons
 * evaluated at muPEXSI.
 * @param[out] muMinInertia (global) Lower bound for mu after the last
 * inertia count procedure.
 * @param[out] muMaxInertia (global) Upper bound for mu after the last
 * inertia count procedure.
 * @param[out] numTotalInertiaIter (global) Number of total inertia
 * counting procedure.
 * @param[out] numTotalPEXSIIter (global) Number of total %PEXSI
 * evaluation procedure.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIDFTDriver(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            numElectronExact,
    double*           muPEXSI,
    double*           numElectronPEXSI,
    double*           muMinInertia,
    double*           muMaxInertia,
    int*              numTotalInertiaIter,
    int*              numTotalPEXSIIter,
    int*              info );

/**
 * @brief Simplified driver for solving Kohn-Sham DFT. Version 2. This
 * is deprecated.
 *
 * This function contains both the inertia counting step for estimating
 * the chemical potential, and strategy for updating the chemical
 * potential WITHOUT Newton's iteration.  Heuristics are built into this
 * routine.  Expert users and developers can modify this routine to
 * obtain better heuristics.  The implementation of this function
 * contains **all** the heuristics for a DFT solver.
 *
 * The input parameter options are controlled through the structure
 * PPEXSIOptions.  The default value can be obtained through
 * PPEXSISetDefaultOptions.
 *
 * **Basic strategy of the heuristics**
 *
 * - If isInertiaCount == 1, then the inertia counting procedure is
 * invoked until the chemical potential interval is refined from size
 * (muMax0-muMin0) to muInertiaTolerance, or the estimated band gap is
 * larger than muInertiaTolerance.
 * - The number of shifts in the inertia count is always automatically
 * chosen to be (mpisize / npPerPole). This minimizes the wall clock
 * time of the inertia counting procedure.
 *
 * **Input/Output**
 *
 * The input H and S matrices should be given by loading functions
 * (currently it is PPEXSILoadRealSymmetricHSMatrix).  The output
 * matrices should be obtained from retrieving functions (currently it
 * is PPEXSIRetrieveRealSymmetricDFTMatrix).
 *
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] numElectronExact (global) Exact number of electrons, i.e.
 * \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * - In the case that convergence is reached within maxPEXSIIter steps, the
 * value of muPEXSI is the last mu used to achieve accuracy within
 * numElectronPEXSITolerance.
 * - In the case that convergence is not reached within maxPEXSIIter steps,
 * and the update from Newton's iteration does not exceed
 * muPEXSISafeGuard, the value of muPEXSI is the last mu plus the update
 * from Newton's iteration.
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at the last step.
 * **Note** In the case that convergence is not reached within maxPEXSIIter steps,
 * and numElectron does not correspond to the number of electrons
 * evaluated at muPEXSI.
 * @param[out] muMinInertia (global) Lower bound for mu after the last
 * inertia count procedure.
 * @param[out] muMaxInertia (global) Upper bound for mu after the last
 * inertia count procedure.
 * @param[out] numTotalInertiaIter (global) Number of total inertia
 * counting procedure.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIDFTDriver2_Deprecate(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double            numElectronExact,
    double*           muPEXSI,
    double*           numElectronPEXSI,
    double*           muMinInertia,
    double*           muMaxInertia,
    int*              numTotalInertiaIter,
    int*              info );

/**
 * @brief Simplified driver for solving Kohn-Sham DFT. Version 2.
 *
 * This function contains both the inertia counting step for estimating
 * the chemical potential, and strategy for updating the chemical
 * potential WITHOUT Newton's iteration.  Heuristics are built into this
 * routine.  Expert users and developers can modify this routine to
 * obtain better heuristics.  The implementation of this function
 * contains **all** the heuristics for a DFT solver.
 *
 * The input parameter options are controlled through the structure
 * PPEXSIOptions.  The default value can be obtained through
 * PPEXSISetDefaultOptions.
 *
 * **Basic strategy of the heuristics**
 *
 * - If isInertiaCount == 1, then the inertia counting procedure is
 * invoked until the chemical potential interval is refined from size
 * (muMax0-muMin0) to muInertiaTolerance, or the estimated band gap is
 * larger than muInertiaTolerance.
 * - The number of shifts in the inertia count is always automatically
 * chosen to be (mpisize / npPerPole). This minimizes the wall clock
 * time of the inertia counting procedure.
 *
 * **Input/Output**
 *
 * The input H and S matrices should be given by loading functions
 * (currently it is PPEXSILoadRealSymmetricHSMatrix).  The output
 * matrices should be obtained from retrieving functions (currently it
 * is PPEXSIRetrieveRealSymmetricDFTMatrix).
 *
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[in] numElectronExact (global) Exact number of electrons, i.e.
 * \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * - In the case that convergence is reached within maxPEXSIIter steps, the
 * value of muPEXSI is the last mu used to achieve accuracy within
 * numElectronPEXSITolerance.
 * - In the case that convergence is not reached within maxPEXSIIter steps,
 * and the update from Newton's iteration does not exceed
 * muPEXSISafeGuard, the value of muPEXSI is the last mu plus the update
 * from Newton's iteration.
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at the last step.
 * **Note** In the case that convergence is not reached within maxPEXSIIter steps,
 * and numElectron does not correspond to the number of electrons
 * evaluated at muPEXSI.
 * @param[out] muMinInertia (global) Lower bound for mu after the last
 * inertia count procedure.
 * @param[out] muMaxInertia (global) Upper bound for mu after the last
 * inertia count procedure.
 * @param[out] numTotalInertiaIter (global) Number of total inertia
 * counting procedure.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIDFTDriver2(
    PPEXSIPlan        plan,
    PPEXSIOptions*    options,
    double            numElectronExact,
    //int               method,
    //int               nPoints,
    double*           muPEXSI,
    double*           numElectronPEXSI,
    int*              numTotalInertiaIter,
    int*              info );


/**
 * @brief Retrieve the output matrices after running PPEXSIDFTDriver for real input matrices.
 * this is only used for the PEXSI method = 1
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIRetrieveRealDFTMatrix(
    PPEXSIPlan        plan,
    double*      DMnzvalLocal,
    double*     EDMnzvalLocal,
    double*     FDMnzvalLocal,
    double*     totalEnergyH,
    double*     totalEnergyS,
    double*     totalFreeEnergy,
    int*              info );


/**
 * @brief Retrieve the output matrices after running PPEXSIDFTDriver for complex input matrices.
 * this is only used for the PEXSI method = 1
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out]  DMnzvalLocal (local)  Dimension: 2*nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: 2*nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: 2*nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIRetrieveComplexDFTMatrix(
    PPEXSIPlan        plan,
    double*      DMnzvalLocal,
    double*     EDMnzvalLocal,
    double*     FDMnzvalLocal,
    double*     totalEnergyH,
    double*     totalEnergyS,
    double*     totalFreeEnergy,
    int*              info );


/**
 * @brief Release the memory used by %PEXSI.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIPlanFinalize(
    PPEXSIPlan plan,
    int*       info );

/**
 * @brief Pole expansion for density matrix (DM).
 *
 * @param[out] zshift  (global) Quadrature point (shift).
 * @param[out] zweight (global) Quadrature weight.
 * @param[in]  Npole   (global) Number of poles.
 * @param[in]  temp    (global) Temperature (unit: au).
 * @param[in]  gap     (global) Energy gap (unit: au).
 * @param[in]  deltaE  (global) Spectral radius (unit: au).
 * @param[in]  mu      (global) Chemical potential (unit: au).
 */
void PPEXSIGetPoleDM(
    double*       zshift,
    double*       zweight,
    int           Npole,
    double        temp,
    double        gap,
    double        deltaE,
    double        mu );

/**
 * @brief Pole expansion for energy density matrix (EDM).
 *
 * @param[out] zshift  (global) Quadrature point (shift).
 * @param[out] zweight (global) Quadrature weight.
 * @param[in]  Npole   (global) Number of poles.
 * @param[in]  temp    (global) Temperature (unit: au).
 * @param[in]  gap     (global) Energy gap (unit: au).
 * @param[in]  deltaE  (global) Spectral radius (unit: au).
 * @param[in]  mu      (global) Chemical potential (unit: au).
 */
void PPEXSIGetPoleEDM(
    double*       zshift,
    double*       zweight,
    int           Npole,
    double        temp,
    double        gap,
    double        deltaE,
    double        mu );

/**
 * @brief Pole expansion for free energy density matrix (FDM).
 *
 * @param[out] zshift  (global) Quadrature point (shift).
 * @param[out] zweight (global) Quadrature weight.
 * @param[in]  Npole   (global) Number of poles.
 * @param[in]  temp    (global) Temperature (unit: au).
 * @param[in]  gap     (global) Energy gap (unit: au).
 * @param[in]  deltaE  (global) Spectral radius (unit: au).
 * @param[in]  mu      (global) Chemical potential (unit: au).
 */
void PPEXSIGetPoleFDM(
    double*       zshift,
    double*       zweight,
    int           Npole,
    double        temp,
    double        gap,
    double        deltaE,
    double        mu );

/**
 * @brief Pole expansion for correct the complex EDM matrix
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 *
 */
void PPEXSICalculateEDMCorrectionComplex(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info );

/**
 * @brief Pole expansion for correct the complex EDM matrix
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 *
 */
void PPEXSICalculateEDMCorrectionReal(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    int*              info );

/**
 * @brief interpolate the Density matrix(DM) and Energy Density Matrix EDM
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 *
 */
void PPEXSIInterpolateDMReal(
    PPEXSIPlan        plan,
    PPEXSIOptions*    options,
    double            numElectronExact,
    double            numElectronPEXSI,
    double *          NeVec,
    double *          muPEXSI,
    int*              info );

/**
 * @brief interpolate the complex Density matrix(DM) and Energy Density Matrix EDM
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] options (global) Other input parameters for the DFT driver.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 *
 */
void PPEXSIInterpolateDMComplex(
    PPEXSIPlan        plan,
    PPEXSIOptions*    options,
    double            numElectronExact,
    double            numElectronPEXSI,
    double *          NeVec,
    double *          muPEXSI,
    int*              info );

/**
 * @brief Retrieve the output DM matrices after running PPEXSIDFTDriver for real input matrices.
 * this is only used for the PEXSI method = 2
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out]  totalEnergyH(local)  H*DM energy
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIRetrieveRealDM(
    PPEXSIPlan  plan,
    double*     DMnzvalLocal,
    double*     totalEnergyH,
    int*        info );

/**
 * @brief Retrieve the output DM matrices after running PPEXSIDFTDriver for real input matrices.
 * this is only used for the PEXSI method = 2, S inverse will be computed in this procedure when
 * it is not identity.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out]  totalEnergyS(local)  Tr[S*EDM] energy
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIRetrieveRealEDM(
    PPEXSIPlan  plan,
    PPEXSIOptions     options,
    double*     EDMnzvalLocal,
    double*     totalEnergyS,
    int*        info );

/**
 * @brief Retrieve the output FDM matrices after running PPEXSIDFTDriver for real input matrices.
 * this is only used for the PEXSI method = 1 and 3
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out]  totalEnergyF(local)  Free energy
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
void PPEXSIRetrieveRealFDM(
    PPEXSIPlan  plan,
    PPEXSIOptions     options,
    double*     FDMnzvalLocal,
    double*     totalEnergyF,
    int*        info );
 */


/**
 * @brief Retrieve the output matrices after running PPEXSIDFTDriver for complex input matrices.
 * this is only used for the PEXSI method = 2
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out]  DMnzvalLocal (local)  Dimension: 2*nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIRetrieveComplexDM(
    PPEXSIPlan        plan,
    double*      DMnzvalLocal,
    double*     totalEnergyH,
    int*              info );

/**
 * @brief Retrieve the output matrices after running PPEXSIDFTDriver for complex input matrices.
 * this is only used for the PEXSI method = 2, S inverse will be computed in this procedure when
 * it is not identity.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[out] EDMnzvalLocal (local)  Dimension: 2*nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] info (local) whether the current processor returns the correct information.
 * - = 0: successful exit.
 * - > 0: unsuccessful.
 */
void PPEXSIRetrieveComplexEDM(
    PPEXSIPlan        plan,
    PPEXSIOptions     options,
    double*     EDMnzvalLocal,
    double*     totalEnergyS,
    int*              info );



#ifdef __cplusplus
}// extern "C"
#endif

#endif // _PEXSI_C_PEXSI_INTERFACE_H_
