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
/// @file ppexsi.hpp
/// @brief Main class for parallel %PEXSI.
/// @date Original:      2012-11-20  Initially started.
/// @date Revision:      2014-03-09  Second generation interface.
#ifndef _PPEXSI_HPP_
#define _PPEXSI_HPP_
#include "pexsi/environment.hpp"
#include "pexsi/sparse_matrix.hpp"
#include "pexsi/NumVec.hpp"
#include "pexsi/utility.hpp"
#include "pexsi/pole.hpp"
#include "pexsi/mpi_interf.hpp"
#include "pexsi/SuperLUGrid.hpp"
#include "pexsi/superlu_dist_interf.hpp"
#include "pexsi/pselinv.hpp"
//#include "pexsi/ngchol_interf.hpp"
//#include "pexsi/c_pexsi_interface.h"

namespace PEXSI{

	/// @class PPEXSIData
	///
	/// @brief Main class for parallel %PEXSI.
	///
	class PPEXSIData{
	private:
		// *********************************************************************
		// Computational variables
		// *********************************************************************

		std::vector<Complex>  zshift_;      // Complex shift for the pole expansion
		std::vector<Complex>  zweightRho_;  // Complex weight for the pole expansion for density
		std::vector<Complex>  zweightRhoDrvMu_;  // Complex weight for the pole expansion for derivative of the Fermi-Dirac with respect to the chemical potential
		std::vector<Complex>  zweightRhoDrvT_;   // Complex weight for the pole expansion for derivative of the Fermi-Dirac with respect to the temperature T (1/beta, in au)
		std::vector<Complex>  zweightHelmholtz_;  // Complex shift for the pole expansion for Helmholtz free energy
		std::vector<Complex>  zweightForce_;  // Complex weight for the pole expansion for force

    // Outer layer communicator. Also used for distributing the
    // DistSparseMatrix.  Each DistSparseMatrix is replicated in the row
    // (numPoleGroup) direction of gridPole.
		const GridType*           gridPole_;          
		const GridType*           gridSelInv_;        // Inner layer communicator for SelInv
	
    // Inner layer communicator for SuperLU factorization
   	const SuperLUGrid<Real>*       gridSuperLUReal_;           
   	const SuperLUGrid<Complex>*    gridSuperLUComplex_;           


    DistSparseMatrix<Real>     HRealMat_;
    DistSparseMatrix<Real>     SRealMat_;

    DistSparseMatrix<Real>     shiftRealMat_;
    DistSparseMatrix<Complex>  shiftComplexMat_;

    DistSparseMatrix<Real>     shiftInvRealMat_;
    DistSparseMatrix<Complex>  shiftInvComplexMat_;

		DistSparseMatrix<Real>     rhoRealMat_;                   // Density matrix 
	  DistSparseMatrix<Real>     rhoDrvMuRealMat_;              // Derivative of the Fermi-Dirac with respect to mu
		DistSparseMatrix<Real>     rhoDrvTRealMat_;               // Derivative of the Fermi-Dirac with respect to T
		DistSparseMatrix<Real>     freeEnergyDensityRealMat_;     // Helmholtz free energy density matrix
		DistSparseMatrix<Real>     energyDensityRealMat_;         // Energy density matrix for computing the Pulay force

    // SuperLUMatrix and PMatrix structures These structures are saved
    // to avoid repetitive symbolic factorization process, and saved in
    // pointer form because of the constructors.
    SuperLUMatrix<Real>*       luRealMat_;
    SuperLUMatrix<Complex>*    luComplexMat_;

    SuperLUOptions             luOpt_;

    PMatrix<Real>*             PMRealMat_;
    PMatrix<Complex>*          PMComplexMat_;

    // Whether the matrices have been loaded into HRealMat_ and
    // SRealMat_
    bool                       isMatrixLoaded_;
    // Whether the matrices (luMat and PMat) have obtained symbolic
    // information
    bool                       isRealSymmetricSymbolicFactorized_;
    bool                       isComplexSymmetricSymbolicFactorized_;
    // Supernode partition for the real matrix
    SuperNodeType              superReal_;             
    // Supernode partition for the complex matrix
    SuperNodeType              superComplex_;             

		// Saves all the indices of diagonal elements in H, so that
		// H.nzvalLocal(diagIdxLocal_[j]) are diagonal elements for all j.
		// This is manly used when S is implicitly given as an identity matrix.
		std::vector<Int>           diagIdxLocal_;    

    // Energy computed from Tr[H*DM]
    Real                       totalEnergyH_;
    // Energy computed from Tr[S*EDM]
    Real                       totalEnergyS_;
    // Free energy 
    Real                       totalFreeEnergy_;


		// *********************************************************************
		// Saved variables for nonlinear iterations
		// *********************************************************************

	public:
    PPEXSIData(
        MPI_Comm   comm,
        Int        numProcRow, 
        Int        numProcCol, 
        Int        outputFileIndex );

		~PPEXSIData();

    void LoadRealSymmetricMatrix(
        Int           nrows,                        
        Int           nnz,                          
        Int           nnzLocal,                     
        Int           numColLocal,                  
        Int*          colptrLocal,                  
        Int*          rowindLocal,                  
        Real*         HnzvalLocal,                  
        Int           isSIdentity,                  
        Real*         SnzvalLocal,
        Int           verbosity );


    
    /// @brief Symbolically factorize the loaded matrices for real
    /// arithmetic factorization and selected inversion.
    ///
    /// The symbolic information is saved internally at luRealMat_ and
    /// PMRealMat_.
    ///
		/// @param[in] ColPerm   Permutation method used for SuperLU_DIST
		///
		/// @param[in] numProcSymbFact Number of processors used for parallel
		/// symbolic factorization and PARMETIS/PT-SCOTCH.
    /// @param[in] verbosity The level of output information.
    /// - = 0   : No output.
    /// - = 1   : Basic output (default)
    /// - = 2   : Detailed output.
    void SymbolicFactorizeRealSymmetricMatrix(
				std::string                    ColPerm,
				Int                            numProcSymbFact,
        Int                            verbosity );
 
   
    /// @brief Symbolically factorize the loaded matrices for complex
    /// arithmetic factorization and selected inversion.
    ///
    /// The symbolic information is saved internally at luComplexMat_ and
    /// PMComplexMat_.
    ///
		/// @param[in] ColPerm   Permutation method used for SuperLU_DIST
		///
		/// @param[in] numProcSymbFact Number of processors used for parallel
		/// symbolic factorization and PARMETIS/PT-SCOTCH.
    /// @param[in] verbosity The level of output information.
    /// - = 0   : No output.
    /// - = 1   : Basic output (default)
    /// - = 2   : Detailed output.
    void SymbolicFactorizeComplexSymmetricMatrix(
				std::string                    ColPerm,
				Int                            numProcSymbFact,
        Int                            verbosity );

    void SelInvRealSymmetricMatrix(
          double*           AnzvalLocal,                  
          Int               verbosity,
          double*           AinvnzvalLocal );

    void SelInvComplexSymmetricMatrix(
          double*           AnzvalLocal,                  
          Int               verbosity,
          double*           AinvnzvalLocal );

		/// @brief Compute the negative inertia (the number of eigenvalues
		/// below a shift) for real symmetric matrices.  The factorization
    /// uses real arithemetic factorization routine.
		///
		/// This subroutine computes the negative inertia of the matrix
		///
		/// I = H - shift * S
		///
		/// where I is the same as the number of eigenvalues lambda for
		///
		/// H x = lambda S x
		///
		/// with lambda < shift according to the Sylvester's law of inertia.
		///
		/// @param[in]  shiftVec Shift vectors.
		/// @param[out] inertiaVec Negative inertia count, the same size as
		/// shiftVec.
		/// @param[in] HMat Hamiltonian matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix.
		/// @param[in] SMat Overlap matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix.
		///
		/// **Note**: If SMat.size == 0, SMat is treated as an identity matrix.
		/// 
    /// @param[in] verbosity The level of output information.
    /// - = 0   : No output.
    /// - = 1   : Basic output (default)
    /// - = 2   : Detailed output.
		void CalculateNegativeInertiaReal(
				const std::vector<Real>&       shiftVec, 
				std::vector<Real>&             inertiaVec,
        Int                            verbosity );


    /// @brief Compute the Fermi operator for a given chemical
    /// potential for real symmetric matrices.
		///
    /// This routine also computes the single particle density matrix,
    /// the Helmholtz free energy density matrix, and the energy density
    /// matrix (for computing the Pulay force) simultaneously.   These
    /// matrices can be called later via member functions DensityMatrix,
    /// FreeEnergyDensityMatrix, EnergyDensityMatrix.
		///
		/// @param[in] numPole Number of poles for the pole expansion
		///	@param[in] temperature  Temperature
		/// @param[in] gap Band gap
		/// @param[in] deltaE Upperbound of the spectrum width
		/// @param[in] mu Initial guess of chemical potential.
		/// @param[in] numElectronExact  Exact number of electrons.
    /// @param[in] numElectronTolerance  Tolerance for the number of
    /// electrons. This is just used to discard some poles in the pole
    /// expansion.
    /// @param[in] verbosity The level of output information.
    /// - = 0   : No output.
    /// - = 1   : Basic output (default)
    /// - = 2   : Detailed output.
    /// @param[out] numElectron The number of electron calculated at mu.
    /// @param[out] numElectronDrvMu The derivative of the number of
    /// electron calculated with respect to the chemical potential at mu.
		void CalculateFermiOperatorReal(
				Int   numPole, 
				Real  temperature,
				Real  gap,
				Real  deltaE,
				Real  mu,
        Real  numElectronExact, 
        Real  numElectronTolerance,
        Int   verbosity,
        Real& numElectron,
        Real& numElectronDrvMu );


    /// @brief Main driver for solving KSDFT.
    void DFTDriver(
        Real       numElectronExact,
        Real       temperature,
        Real       gap,
        Real       deltaE,
        Int        numPole, 
        Int        isInertiaCount,
        Int        maxPEXSIIter,
        Real       muMin0,
        Real       muMax0,
        Real       mu0,
        Real       muInertiaTolerance,
        Real       muInertiaExpansion,
        Real       muPEXSISafeGuard,
        Real       numElectronPEXSITolerance,
        Int        matrixType,
        Int        isSymbolicFactorize,
        Int        ordering,
        Int        numProcSymbFact,
        Int        verbosity,
        Real&      muPEXSI,                   
        Real&      numElectronPEXSI,         
        Real&      muMinInertia,              
        Real&      muMaxInertia,             
        Int&       numTotalInertiaIter,   
        Int&       numTotalPEXSIIter );


    // *********************************************************************
    // Access data
    // *********************************************************************

    const GridType*  GridPole() const {return gridPole_;}

    /// @brief Density matrix.
    ///
    /// Can be used to estimate the number of electrons by Tr[DM*S]
    /// or the band energy via Tr[DM*H]
    const DistSparseMatrix<Real>&   RhoRealMat() const {return rhoRealMat_;}

    /// @brief Energy density matrix.  
    ///
    /// Can be used to estimate the total band energy via Tr[EDM*S] or
    /// the force, including the Hellman-Feynman force and the Pulay
    /// force. 
    const DistSparseMatrix<Real>&   EnergyDensityRealMat() const {return energyDensityRealMat_;}

    /// @brief Total Helmholtz free energy matrix (band energy part only).  
    ///
    /// The Helmholtz free energy is computed by Tr[rho_f*H].
		///
		/// For more information see 
		/// Alavi, A., Kohanoff, J., Parrinello, M., & Frenkel, D. (1994). Ab
		/// initio molecular dynamics with excited electrons. Physical review
		/// letters, 73(19), 2599–2602. 
    const DistSparseMatrix<Real>&   FreeEnergyDensityRealMat() const {return freeEnergyDensityRealMat_;}

    Real   TotalEnergyH() const {return totalEnergyH_;}
    
    Real   TotalEnergyS() const {return totalEnergyS_;}

    Real   TotalFreeEnergy() const {return totalFreeEnergy_;}


	}; // PPEXSIData


} // namespace PEXSI
#endif // _PPEXSI_HPP_
