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
 * @file c_pexsi_interface_old.h
 * @brief Interface subroutines of %PEXSI that can be called by C.
 *
 * The FORTRAN interface are wrappers of the C interface.  The header
 * file is not needed here, and the interface routines are given
 * directly in interface.cpp.
 * 
 * @date 2013-01-31
 */
#ifndef _PEXSI_C_PEXSI_INTERFACE_H_ 
#define _PEXSI_C_PEXSI_INTERFACE_H_
#include <mpi.h>

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



/**
 * @brief Driver interface for computing the selected  elements of a
 * matrix (in complex arithmietic).  
 *
 * Evaluate the selected elements of 
 * \f$(H - z S)^{-1}\f$ for a shift \f$z\in\mathbb{C}\f$.
 * 
 * **Note** 
 * - H and S are always real symmetric matrices, are saved in the
 *   compressed sparse column (CSC) format, and have the same sparsity
 *   pattern.
 * - Complex arithmetic operation is used for the factorization
 *   and selected inversion even if z is real.
 * - The input of S is optional, and the shift z can be 0.
 *
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
 * @param[in] zShift  (global) Dimension: 2. Use 2 real numbers to
 * denote a complex shift. 
 * - Real part: zShift[0]
 * - Imag part: zShift[1]    
 * @param[in] ordering (global) Ordering strategy for factorization and selected
 * inversion.  
 * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
 *   option in SuperLU_DIST).
 * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
 *   option in SuperLU_DIST).
 * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
 *   option in SuperLU_DIST).
 * @param[in] npSymbFact (global) Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
 * @param[in] comm (global) MPI communicator.
 * @param[out] AinvnzvalLocal (local) Dimension: 2*nnzLocal. Local nonzero
 * values of the selected elements of \f$(H - z S)^{-1}\f$.
 * - Use 2 double for one complex number. This ensures the compatibility with FORTRAN.
 * - Real part: AinvnzvalLocal[2*k]. Imag part: AinvnzvalLocal[2*k+1].
 * @param[out] info 
 * - = 0: successful exit.  
 * - > 0: unsuccessful.
 */
void PPEXSISelInvInterface ( 
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       HnzvalLocal,                  
    int           isSIdentity,                  
    double*       SnzvalLocal,                  
    double*       zShift,                       
    int           ordering,                
    int           npSymbFact,                   
    MPI_Comm	    comm,                         
    double*       AinvnzvalLocal,
    int*          info
    );

/**
 * @brief Driver interface for estimating the chemical potential using
 * inertia counting procedure.
 *
 * The main purpose of this routine is to obtain two set of intervals.
 *
 * The first interval is [muMinInertia, muMaxInertia]. This should be a
 * tight interval so that 
 *   \f[
 *     N_e(\mathrm{muMinInertia}) \le \mathrm{numElectronExact} 
 *     \le N_e(\mathrm{muMaxInertia}).
 *   \f]
 * This interval can be used as a bound interval for @ref PPEXSISolveInterface.
 *
 * The second interval is [muLowerEdge, muUpperEdge]. This interval is
 * defined by
 *   \f[
 *     N_e(\mathrm{muLowerEdge}) = \mathrm{numElectronExact} - eps,
 *   \f]
 *   \f[
 *     N_e(\mathrm{muUpppeEdge}) = \mathrm{numElectronExact} + eps.
 *   \f]
 *   eps is a small number currently fixed to be 0.1.
 *   For metals (\f$N'_e(\mu_{\mathrm{exact}}) > 0\f$), it should be
 *   expected that muLowerEdge is very close to muUpperEdge. For
 *   insulators, (\f$N'_e(\mu_{\mathrm{exact}}) \approx 0\f$),
 *   we have muUpperEdge > muLowerEdge, with the difference
 *   characterizing the band gap.
 *
 *   For either metal or insulator, 
 *   \f[
 *       \frac{1}{2} ( \mathrm{muLowerEdge} + \mathrm{muUpperEdge} ) 
 *   \f]
 *   is a good guess for the chemical potential to be provided to 
 *   @ref PPEXSISolveInterface.
 *
 * These intervals are refined for in an iterative procedure, until
 * \f[
 *   N_e(\mathrm{muMaxInertia}) - N_e(\mathrm{muMinInertia}) <
 *   numElectronTolerance.
 * \f]
 *
 * The counting of the number of electrons (eigenvalues) is obtained by
 * counting the number of negative eigenvalues for the matrix \f$H-\mu_l
 * S\f$ for a series of \f$\{\mu_l\}\f$. The algorithm makes use of
 * Sylvester's law of inertia.  
 *
 * Finite temperature effect is also taken into account when computing
 * the number of electrons.
 * 
 * **Note** 
 * - H and S are always real symmetric matrices, are saved in the
 *   compressed sparse column (CSC) format, and have the same sparsity
 *   pattern.
 * - This algorithm uses real arithmetic LU decomposition.
 * - The input of S is optional, and the shift z can be 0.
 *
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
 * @param[in] temperature (global) Temperature, in the same unit as H
 * @param[in] numElectronExact (global) Exact number of electrons, i.e. \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[in] muMin0   (global) Initial guess of lower bound for mu
 * @param[in] muMax0   (global) Initial guess of upper bound for mu
 * @param[in] numPole  (global) Number of shifts (chemical potentials) evaluated at each iteration.
 * @param[in] maxIter  (global) Maximum number of iterations.
 * @param[in] numElectronTolerance (global) Stopping criterion of the
 * iterations.
 * @param[in] ordering (global) Ordering strategy for factorization and selected
 * inversion.  
 * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
 *   option in SuperLU_DIST).
 * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
 *   option in SuperLU_DIST).
 * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
 *   option in SuperLU_DIST).
 * @param[in] npSymbFact (global) Number of processors for
 * PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
 * @param[in] comm (global) MPI communicator.
 * @param[out] muMinInertia (global) Lower bound for mu after the inertia count.
 * @param[out] muMaxInertia (global) Upper bound for mu after the inertia count.
 * @param[out] muLowerEdge  (global) muLowerEdge satisfies
 * N_e(muLowerEdge) = numElectronExact - eps.  
 * @param[out] muUpperEdge  (global) muUpperEdge satisfies
 * N_e(muUpperEdge) = numElectronExact + eps.  
 * @param[out] numIter      (global) Number of iterations used before
 * reaching the stopping criterion.
 * @param[out] muList       (global) Dimension: numIter. The list of shifts after the final
 * iteration.
 * @param[out] numElectronList (global) Dimension: numIter. The number of electrons at the
 * each shift. Finite temperature effect is taken into account. For
 * obtaining the number of electrons without the finite temperature
 * effect, see @ref PPEXSIRawInertiaCountInterface.
 * @param[out] info 
 * - = 0: successful exit.  
 * - > 0: unsuccessful.
 */
void PPEXSIInertiaCountInterface(
    int           nrows,  
    int           nnz,    
    int           nnzLocal, 
    int           numColLocal,
    int*          colptrLocal,
    int*          rowindLocal,
    double*       HnzvalLocal,
    int           isSIdentity,
    double*       SnzvalLocal,
    double        temperature,
    double        numElectronExact,
    double        muMin0,          
    double        muMax0,         
    int           numPole,       
    int           maxIter,      
    double        numElectronTolerance,
    int           ordering,           
    int           npPerPole,         
    int           npSymbFact,       
    MPI_Comm	    comm,            
    double*       muMinInertia,                
    double*       muMaxInertia,               
    double*       muLowerEdge,               
    double*       muUpperEdge,              
    int*          numIter,                 
    double*       muList,                 
    double*       numElectronList,       
    int*          info                  
    );


/** 
 * @brief Directly compute the negative inertia at a set of shifts.
 * 
 *  This is a simplified version of @ref PPEXSIInertiaCountInterface,
 *  without the including the finite temperature effect, or
 *  iterative refinement of the chemical potential.  This routine can be
 *  used for evaluating density of states after @ref PPEXSISolveInterface 
 *  has reached convergence.
 * 
 *  See @ref PPEXSIInertiaCountInterface for explanation of
 *  parameters.
 */
void PPEXSIRawInertiaCountInterface(
		// Input parameters
		int           nrows,  
		int           nnz,    
		int           nnzLocal, 
		int           numColLocal, 
		int*          colptrLocal, 
		int*          rowindLocal, 
		double*       HnzvalLocal, 
		int           isSIdentity, 
		double*       SnzvalLocal,
		double        muMin0,    
		double        muMax0,   
		int           numPole, 
		int           ordering,
		int           npPerPole, 
		int           npSymbFact,
		MPI_Comm	    comm,     
		double*       muList, 
		int   *       numElectronList,   
		int*          info              
		);


/** 
 * @brief Main driver routine for solving Kohn-Sham density functional
 * theory using %PEXSI.
 * 
 * Given a good initial guess of the chemical potential and the bound
 * for the chemical potential, this routine
 * uses the Newton method to find the converged chemical potential for
 * a given number of electrons at certain temperature.
 * 
 * After convergence, this routine also returns the density matrix for
 * evaluating the electron density.  It also returns the free energy
 * density matrix for evaluating the (Helmholtz or Mermin) free energy,
 * and the energy density matrix for evaluating the atomic force.
 *
 * **Note**
 * - H and S are always real symmetric matrices, are saved in the
 *   compressed sparse column (CSC) format, and have the same sparsity
 *   pattern.
 * - The input of S is optional, and the shift z can be 0.
 * - The initial guess for the chemical potential mu0, and the lower / 
 *   upper bounds muMin0 / muMax0 should be obtained from 
 *   @ref PPEXSIInertiaCountInterface, or from heuristics using
 *   information from previous steps.
 * - In the current version of the interface routines, the number of 
 * processors used for each selected inversion process must be a square 
 * number.
 *
 * @todo The estimation of deltaE should be estimated in a more
 * automatic way later, using a few steps of Lanczos.
 *
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
 * @param[in] temperature (global) Temperature, in the same unit as H
 * @param[in] numElectronExact (global) Exact number of electrons, i.e. \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[in] mu0      (global) Initial guess for mu.
 * @param[in] muMin0   (global) Initial guess of lower bound for mu.
 * @param[in] muMax0   (global) Initial guess of upper bound for mu.
 * @param[in] gap      (global) A lower bound for the band gap.
 * Currently it is recommended to set this to 0, which works for all
 * systems.
 * @param[in] deltaE   (global) An upper bound for the spectral radius
 * of \f$S^{-1} H\f$.  
 * @param[in] numPole  (global) Number of poles.
 * @param[in] maxIter  (global) Maximum number of iterations.
 * @param[in] numElectronTolerance (global) Stopping criterion of the
 * iterations.
 * @param[in] ordering (global) Ordering strategy for factorization and selected
 * inversion.  
 * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
 *   option in SuperLU_DIST).
 * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
 *   option in SuperLU_DIST).
 * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
 *   option in SuperLU_DIST).
 * @param[in] npPerPole  (global) Number of processors per pole
 * (ppp).  
 * @param[in] npSymbFact (global) Number of processors for
 * PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
 * @param[in] comm (global) MPI communicator.
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at muPEXSI after the last iteration.
 * @param[out] muMinPEXSI   (global) Lower bound for mu after the last
 * iteration.
 * @param[out] muMaxPEXSI   (global) Upper bound for mu after the last
 * iteration.
 * @param[out] numIter      (global) Number of iterations used before
 * reaching the stopping criterion.
 * @param[out] muList       (global) Dimension: numIter. The history of chemical potential.
 * @param[out] numElectronList (global) Dimension: numIter. The history of the number of
 * electrons evaluated at each mu in muList.
 * @param[out] numElectronDrvList (global) Dimension: numIter. The history of the derivative
 * of the number of electrons with respect to mu, evaluated at each mu
 * in muList. This is used in the Newton iteration.
 * @param[out] info 
 * - = 0: successful exit.  
 * - > 0: unsuccessful.
 */ 
void PPEXSISolveInterface (
		int           nrows,                        
	  int           nnz,                          
		int           nnzLocal,                     
		int           numColLocal,                  
		int*          colptrLocal,                  
		int*          rowindLocal,                  
		double*       HnzvalLocal,                  
		int           isSIdentity,                  
		double*       SnzvalLocal,                 
		double        temperature,                 
		double        numElectronExact,            
		double        mu0,                         
		double        muMin0,                      
		double        muMax0,                      
		double        gap,                        
		double        deltaE,                    
		int           numPole,                  
		int           maxIter,                 
		double        numElectronTolerance,   
		int           ordering,               
	  int           npPerPole,             
		int           npSymbFact,           
	  MPI_Comm	    comm,                 
		double*      DMnzvalLocal,                  
		double*     EDMnzvalLocal,                 
		double*     FDMnzvalLocal,                
		double*       muPEXSI,                   
		double*       numElectronPEXSI,         
		double*       muMinPEXSI,              
		double*       muMaxPEXSI,             
		int*          numIter,               
		double*       muList,               
		double*       numElectronList,     
		double*       numElectronDrvList, 
		int*          info                
		);


/**
 *  @brief Driver interface for computing the local density of states.
 *
 *  The local density of states (LDOS) can be evaluated as
 *  \f[
 *     n(r;E) = \lim_{\eta->0+} \mathrm{Im} \frac{1}{\pi} 
 *       \langle r\vert (H-(E+i \eta)I)^{-1} \vert r \rangle.
 *  \f]
 *  For nonorthogonal basis functions, this routine returns
 *  of the matrix
 *  \f[
 *     \frac{1}{\pi} \mathrm{Im} ( H - (E+i \eta) S )^{-1},
 *  \f]
 *  and the LDOS in the real space can be constructed in the same way
 *  that the density is constructed.
 *
 * **Note** 
 * - H and S are always real symmetric matrices, are saved in the
 *   compressed sparse column (CSC) format, and have the same sparsity
 *   pattern.
 * - The input of S is optional, and the shift z can be 0.
 *
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
 * @param[in] Energy      (global) Real part of the shift.
 * @param[in] eta         (global) Imaginary part of the shift, or the
 * "broadening" parameter \f$\eta\f$.
 * @param[in] ordering (global) Ordering strategy for factorization and selected
 * inversion.  
 * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
 *   option in SuperLU_DIST).
 * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
 *   option in SuperLU_DIST).
 * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
 *   option in SuperLU_DIST).
 * @param[in] npSymbFact (global) Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
 * @param[in] comm (global) MPI communicator.
 * @param[out] localDOSnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of the selected elements of \f$\frac{1}{\pi} \mathrm{Im} ( H - (E+i \eta) S )^{-1}\f$.
 * @param[out] info 
 * - = 0: successful exit.  
 * - > 0: unsuccessful.
 */    
void PPEXSILocalDOSInterface (
		int           nrows, 
	  int           nnz,   
		int           nnzLocal, 
		int           numColLocal, 
		int*          colptrLocal, 
		int*          rowindLocal,
		double*       HnzvalLocal,
		int           isSIdentity, 
		double*       SnzvalLocal, 
		double        Energy,      
		double        eta,         
		int           ordering,   
		int           npSymbFact,
	  MPI_Comm	    comm,        
		double*       localDOSnzvalLocal, 
		int*          info               
		);


/**
 * @brief Simplified driver interface for computing the selected
 * elements of a real symmetric symmetric.
 *
 * @param[in] nrows (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros of H.
 * @param[in] nnzLocal (local) Number of local nonzeros of H.
 * @param[in] numColLocal (local) Number of local columns for H.
 * @param[in] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[in] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[in] AnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of A in CSC format.  
 * @param[in] ordering (global) Ordering strategy for factorization and selected
 * inversion.  
 * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
 *   option in SuperLU_DIST).
 * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
 *   option in SuperLU_DIST).
 * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
 *   option in SuperLU_DIST).
 * @param[in] npSymbFact (global) Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
 * @param[in] comm (global) MPI communicator. NOTE: mpisize should be
 * equal to nprow * npcol
 * @param[in] nprow (global) Number of row processors.
 * @param[in] npcol (global) Number of column processors.  
 * @param[out] AinvnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of the selected elements of \f$A^{-1}\f$.
 * @param[out] info 
 * - = 0: successful exit.  
 * - > 0: unsuccessful.
 */
void PSelInvRealSymmetricInterface ( 
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       AnzvalLocal,                  
    int           ordering,                
    int           npSymbFact,                   
    MPI_Comm	    comm,
    int           nprow,
    int           npcol,
    double*       AinvnzvalLocal,
    int*          info
    );


/**
 * @brief Simplified driver interface for computing the selected
 * elements of a complex symmetric symmetric.
 *
 * @param[in] nrows (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros of H.
 * @param[in] nnzLocal (local) Number of local nonzeros of H.
 * @param[in] numColLocal (local) Number of local columns for H.
 * @param[in] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[in] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[in] AnzvalLocal (local) Dimension: 2*nnzLocal. Local nonzero
 * values of A in CSC format.  
 * - Use 2 double for one complex number. This
 * ensures the compatibility with FORTRAN.  
 * - Real part: AnzvalLocal[2*k]. Imag part: AnzvalLocal[2*k+1].
 * @param[in] ordering (global) Ordering strategy for factorization and selected
 * inversion.  
 * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
 *   option in SuperLU_DIST).
 * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
 *   option in SuperLU_DIST).
 * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
 *   option in SuperLU_DIST).
 * @param[in] npSymbFact (global) Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
 * @param[in] comm (global) MPI communicator. NOTE: mpisize should be
 * equal to nprow * npcol
 * @param[in] nprow (global) Number of row processors.
 * @param[in] npcol (global) Number of column processors.  
 * @param[out] AinvnzvalLocal (local) Dimension: 2*nnzLocal. Local nonzero
 * values of the selected elements of \f$A^{-1}\f$.
 * - Use 2 double for one complex number. This ensures the compatibility with FORTRAN.
 * - Real part: AinvnzvalLocal[2*k]. Imag part: AinvnzvalLocal[2*k+1].
 * @param[out] info 
 * - = 0: successful exit.  
 * - > 0: unsuccessful.
 */
void PSelInvComplexSymmetricInterface ( 
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       AnzvalLocal,                  
    int           ordering,                
    int           npSymbFact,                   
    MPI_Comm	    comm,
    int           nprow,
    int           npcol,
    double*       AinvnzvalLocal,
    int*          info
    );


#ifdef __cplusplus
}// extern "C"
#endif

#endif // _PEXSI_C_PEXSI_INTERFACE_H_

