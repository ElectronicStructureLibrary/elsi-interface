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
/// @file ppexsi_old.cpp
/// @brief Implementation of the parallel %PEXSI.
/// @date 2012-11-20  Initially started.
/// @date 2014-03-09  Revise for the new interface.
///
///
#include "ppexsi.hpp"

namespace PEXSI{

PPEXSIData::PPEXSIData	( const GridType* g, Int nprow, Int npcol ): gridPole_(g)
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::PPEXSIData");
#endif
	if( nprow != npcol ){
		#ifdef USE_ABORT
abort();
#endif
throw std::runtime_error( "PSelInv only allows to use square grid." );
	}
	if( gridPole_->numProcCol != nprow * npcol ){
		#ifdef USE_ABORT
abort();
#endif
throw std::runtime_error( "The number of processors numProcCol do not match nprow * npcol." );
	}
	zGridSuperLU_  = new SuperLUGrid<Complex>( gridPole_->rowComm, nprow, npcol );
	dGridSuperLU_  = new SuperLUGrid<Real>( gridPole_->rowComm, nprow, npcol );
	gridSelInv_   = new GridType( gridPole_->rowComm, nprow, npcol );

  
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::PPEXSIData  ----- 


PPEXSIData::~PPEXSIData	(  )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::~PPEXSIData");
#endif
	if( zGridSuperLU_ != NULL ){
		delete zGridSuperLU_;
	}
	
	if( dGridSuperLU_ != NULL ){
		delete dGridSuperLU_;
	}
	if( gridSelInv_ != NULL ){
		delete gridSelInv_;
	}
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::~PPEXSIData  ----- 


Real PPEXSIData::CalculateChemicalPotentialNewtonFD ( 
			const Int iter, 
			const Real numElectronExact, 
			const Real numElectronTolerance, 
			const std::vector<Real>& muList,
			const std::vector<Real>& numElectronList )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateChemicalPotentialNewtonFD");
#endif
  // FIXME Magic number here
	Real  muMin = -5.0, muMax = 5.0, muMinStep = 0.01;;
	Real  muNew;

	if( iter == 0 ){
		if( numElectronExact > numElectronList[iter] ){
			muNew = muList[iter] + muMinStep;
		}
		else{
			muNew = muList[iter] - muMinStep;
		}
	}
	else{
		if( std::abs(numElectronList[iter] -  numElectronList[iter-1])
		    < numElectronTolerance ){
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "The number of electrons did not change." << std::endl;
#endif
			if( numElectronExact > numElectronList[iter] ){
				muNew = muList[iter] + muMinStep;
			}
			else{
				muNew = muList[iter] - muMinStep;
			}
		}
		else {
			muNew = muList[iter] + (muList[iter] - muList[iter-1]) / 
				( numElectronList[iter] - numElectronList[iter-1] ) *
				( numElectronExact - numElectronList[iter] );
			if( muNew < muMin || muNew > muMax ){
#if ( _DEBUGlevel_ >= 0 )
				statusOFS << "muNew = " << muNew << " is out of bound ["
					<< muMin << ", " << muMax << "]" << std::endl;
#endif
				if( numElectronExact > numElectronList[iter] ){
					muNew = muList[iter] + muMinStep;
				}
				else{
					muNew = muList[iter] - muMinStep;
				}
			}
		} // if ( numElectron changed )
	} // if (iter == 0)

#ifndef _RELEASE_
	PopCallStack();
#endif

	return muNew;
} 		// -----  end of method PPEXSIData::CalculateChemicalPotentialNewtonFD  ----- 

Real PPEXSIData::CalculateChemicalPotentialNewtonBisection ( 
			const Real numElectronExact, 
			const Real numElectron,
			const Real numElectronDrvMu,
			const Real mu,
			const Real muMin,
			const Real muMax )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateChemicalPotentialNewtonBisection");
#endif
	Real  muNew, muNewton;

	muNewton = mu - ( numElectron - numElectronExact ) / numElectronDrvMu;

	if( muNewton < muMin || muNewton > muMax ){
		// Bisection method
		if( numElectron < numElectronExact )
			muNew = 0.5 * ( mu + muMax );
		else
			muNew = 0.5 * ( mu + muMin );
	}
	else
		muNew = muNewton;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return muNew;
} 		// -----  end of method PPEXSIData::CalculateChemicalPotentialNewtonBisection  ----- 

// Main subroutine for the electronic structure calculation
void PPEXSIData::Solve( 
		Int  numPole, 
		Real temperature,
		Real numElectronExact,
		Real gap,
		Real deltaE,
		Real& mu,
		Real& muMin,
		Real& muMax,
		const DistSparseMatrix<Real>&  HMat,
		const DistSparseMatrix<Real>&  SMat,
		Int  muMaxIter,
		Real numElectronTolerance,
		std::string         ColPerm,
		Int                 numProcSymbFact,
		bool isFreeEnergyDensityMatrix, 
		bool isEnergyDensityMatrix,
		bool isDerivativeTMatrix,
		std::vector<Real>&	muList,
		std::vector<Real>&  numElectronList,
		std::vector<Real>&  numElectronDrvMuList,
		bool&               isConverged
		){
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::Solve");
#endif


	// *********************************************************************
	// Check the input parameters
	// *********************************************************************
	if( numPole % 2 != 0 ){
		#ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "Must be even number of poles!" );
	}

	// TODO Check H and S have the same pattern

	// TODO Check H and S agree with gridPole_


	// *********************************************************************
	// Initialize
	// *********************************************************************
	muList.clear();
	numElectronList.clear();
	numElectronDrvMuList.clear();
	isConverged = false;

	DistSparseMatrix<Complex>  AMat;              // A = H - z * S
	// rename for convenience
	DistSparseMatrix<Real>& rhoMat       = rhoMat_;     
	DistSparseMatrix<Real>& rhoDrvMuMat  = rhoDrvMuMat_;
	DistSparseMatrix<Real>& rhoDrvTMat   = rhoDrvTMat_;
	DistSparseMatrix<Real>& hmzMat       = freeEnergyDensityMat_;
	DistSparseMatrix<Real>& frcMat       = energyDensityMat_;

	// Get the diagonal indices for H and save it n diagIdxLocal_
	{
		Int numColLocal      = HMat.colptrLocal.m() - 1;
		Int numColLocalFirst = HMat.size / gridSelInv_->mpisize;
		Int firstCol         = gridSelInv_->mpirank * numColLocalFirst;
		
		diagIdxLocal_.clear();

		for( Int j = 0; j < numColLocal; j++ ){
			Int jcol = firstCol + j + 1;
			for( Int i = HMat.colptrLocal(j)-1; 
				 	 i < HMat.colptrLocal(j+1)-1; i++ ){
				Int irow = HMat.rowindLocal(i);
				if( irow == jcol ){
					diagIdxLocal_.push_back( i );
				}
			}
		} // for (j)
	}


	// Copy the pattern
	CopyPattern( HMat, AMat );
	CopyPattern( HMat, rhoMat );
	CopyPattern( HMat, rhoDrvMuMat );
	if( isFreeEnergyDensityMatrix )
		CopyPattern( HMat, hmzMat );
	if( isEnergyDensityMatrix )
		CopyPattern( HMat, frcMat );
	if( isDerivativeTMatrix )
		CopyPattern( HMat, rhoDrvTMat );

	SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value

#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "AMat.nnzLocal = " << AMat.nnzLocal << std::endl;
	statusOFS << "AMat.nnz      = " << AMat.Nnz()    << std::endl;
#endif

	SuperLUOptions   luOpt;

	luOpt.ColPerm = ColPerm;
	luOpt.numProcSymbFact = numProcSymbFact;
	// TODO Introduce maxPipelineDepth as an adjustable parameter when needed.
	luOpt.maxPipelineDepth = -1;

	SuperLUMatrix<Complex>    luMat( *zGridSuperLU_, luOpt );  // SuperLU matrix.

	// *********************************************************************
	// Symbolic factorization.  
	// Each numPoleGroup perform independently
	// *********************************************************************
	Real timeSta, timeEnd;
	GetTime( timeSta );
	luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
	GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "AMat is converted to SuperMatrix." << std::endl;
	statusOFS << "Time for SuperMatrix conversion is " <<
		timeEnd - timeSta << " [s]" << std::endl << std::endl;
#endif
	GetTime( timeSta );
	luMat.SymbolicFactorize();
	GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "Symbolic factorization is finished." << std::endl;
	statusOFS << "Time for symbolic factorization is " <<
		timeEnd - timeSta << " [s]" << std::endl << std::endl;
#endif
	luMat.SymbolicToSuperNode( super_ );
	luMat.DestroyAOnly();


	// Compute the number of nonzeros from PMatrix
	{
		PMatrix<Scalar> PMloc( gridSelInv_, &super_ , &luOpt); // A^{-1} in PMatrix format
		luMat.LUstructToPMatrix( PMloc );
#if ( _DEBUGlevel_ >= 0 )
		Int nnzLocal = PMloc.NnzLocal();
		statusOFS << "Number of local nonzeros (L+U) = " << nnzLocal << std::endl;
		LongInt nnz  = PMloc.Nnz();
		statusOFS << "Number of nonzeros (L+U)       = " << nnz << std::endl;
#endif
	}

#if ( _DEBUGlevel_ >= 1 )
	statusOFS << "perm: "    << std::endl << super_.perm     << std::endl;
	statusOFS << "permInv: " << std::endl << super_.permInv  << std::endl;
	statusOFS << "superIdx:" << std::endl << super_.superIdx << std::endl;
	statusOFS << "superPtr:" << std::endl << super_.superPtr << std::endl; 
#endif

	Real muNow = mu;
	Real numElectronNow;
	Real numElectronDrvMuNow;

	Real timeMuSta, timeMuEnd;

	// Iteration with chemical potentials
	for(Int iter = 0; iter < muMaxIter; iter++){
		GetTime( timeMuSta );

#if ( _DEBUGlevel_ >= 0 )
		{
			std::ostringstream msg;
			msg << "Iteration " << iter << ", mu = " << muNow;
			PrintBlock( statusOFS, msg.str() );
		}
#endif

		// Reinitialize the variables
		SetValue( rhoMat.nzvalLocal, 0.0 );
		SetValue( rhoDrvMuMat.nzvalLocal, 0.0 );
		if( isFreeEnergyDensityMatrix )
			SetValue( hmzMat.nzvalLocal, 0.0 );
		if( isEnergyDensityMatrix )
			SetValue( frcMat.nzvalLocal, 0.0 );
		if( isDerivativeTMatrix )
			SetValue( rhoDrvTMat.nzvalLocal, 0.0 );


		// Refine the pole expansion until 
		// 1) The number of poles used is numPole
		// 2) The worst case bound for the relative error of the pole
		// expansion using numPole poles is less than
		//
		//   numElectronTolerance / numElectronExact  
		
		// numPoleInput is the number of poles to be given to other parts of
		// the pole expansion, which is larger than or equal to numPole.
		Int numPoleInput;
		// poleIdx is of size numPole.  Only poles with index in poleIdx is
		// used for computation. The rest of the poles are discarded
		// according to tolerance criterion
		//
		//   numElectronTolerance / numElectronExact / numPole
		std::vector<Int>  poleIdx(numPole);
		{
			// Setup a grid from (muNow - deltaE, muNow + deltaE), and measure
			// the error bound in the L^infty sense on this grid.
			//
			// fdGrid:      Exact Fermi-Dirac function evaluated on xGrid
			// fdPoleGrid:  Fermi-Dirac function using pole expansion
			// evaluated on the grid.
			Int numX = 10000;
			std::vector<Real>    xGrid( numX );
			std::vector<Real>    fdGrid( numX );

			Real x0 = muNow - deltaE;
			Real x1 = muNow + deltaE;
			Real h  = (x1 - x0) / (numX - 1);
			Real ez;
			for( Int i = 0; i < numX; i++ ){
				xGrid[i]  = x0 + i * h;
				if( xGrid[i] - muNow >= 0 ){
					ez = std::exp(- (xGrid[i] - muNow) / temperature );
					fdGrid[i] = 2.0 * ez / (1.0 + ez);
				}
				else{
					ez = std::exp((xGrid[i] - muNow) / temperature );
					fdGrid[i] = 2.0 / (1.0 + ez);
				}
			}


			numPoleInput = numPole;
			Real tol;
			Int  numPoleSignificant;
				
			Int poleIter = 0;
			do{
				// If the number of significant poles is less than numPole,
				// increase numPoleInput by 2 at a time and redo the
				// computation.
				if( poleIter > 0 )
					numPoleInput += 2;

				zshift_.resize( numPoleInput );
				zweightRho_.resize( numPoleInput );
				GetPoleDensity( &zshift_[0], &zweightRho_[0],
						numPoleInput, temperature, gap, deltaE, muNow ); 

				std::vector<Real>  maxMagPole(numPoleInput);
				for( Int l = 0; l < numPoleInput; l++ ){
					maxMagPole[l] = 0.0;
				}
					
				// Compute the approximation due to pole expansion, as well as
				// the maximum magnitude of each pole
				Complex cpxmag;
				Real    mag;
				numPoleSignificant = 0;
				tol = numElectronTolerance / numElectronExact / numPoleInput;
				for( Int l = 0; l < numPoleInput; l++ ){
					for( Int i = 0; i < numX; i++ ){
						cpxmag = zweightRho_[l] / ( xGrid[i] - zshift_[l] );
						mag    = cpxmag.imag();
						maxMagPole[l] = ( maxMagPole[l] >= mag ) ?  maxMagPole[l] : mag;
					}
					if( maxMagPole[l] > tol ){
						numPoleSignificant++;
					}	
				} // for (l)

				// Pick the most significant numPole poles and update poleIdx
				// Sort in DESCENDING order
				std::vector<Int>  sortIdx( numPoleInput );
				for( Int i = 0; i < sortIdx.size(); i++ ){
					sortIdx[i]      = i;
				}
				std::sort( sortIdx.begin(), sortIdx.end(), 
						IndexComp<std::vector<Real>& >( maxMagPole ) ) ;
				std::reverse( sortIdx.begin(), sortIdx.end() );

				for( Int l = 0; l < numPole; l++ ){
					poleIdx[l]      = sortIdx[l];
				}
				

				// Update poleIter
				poleIter++; 
			} while( numPoleSignificant < numPole );


			// Compute the relative error 
			std::vector<Real>  fdPoleGrid( numX );
			Real errRel, errorRelMax;
			Real errAbs, errorAbsMax;
			Real errorTotal;
			Real errEPS = 1e-1;

			errorRelMax = 0.0;
			errorAbsMax = 0.0; 
			for( Int i = 0; i < numX; i++ ){
				fdPoleGrid[i] = 0.0;
				for( Int lidx = 0; lidx < numPoleInput; lidx++ ){
					Int l = lidx;
					Complex cpxmag = zweightRho_[l] / ( xGrid[i] - zshift_[l] );
					fdPoleGrid[i] += cpxmag.imag();
				}
				errAbs = std::abs( fdPoleGrid[i] - fdGrid[i] );
				errorAbsMax = ( errorAbsMax >= errAbs ) ? errorAbsMax : errAbs;

				if( std::abs(fdGrid[i]) > errEPS ){
					errRel = std::abs( fdPoleGrid[i] - fdGrid[i] ) / ( std::abs( fdGrid[i] ) );
					errorRelMax = ( errorRelMax >= errRel ) ? errorRelMax : errRel;
				}
			}

			errorTotal = errorRelMax * numElectronExact + errorAbsMax * HMat.size;

#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "Pole expansion indicates that the error "
				<< "of numElectron is bounded by "
				<< errorTotal << std::endl;
			statusOFS << "Required accuracy: numElectronTolerance is "
				<< numElectronTolerance << std::endl << std::endl;

			if( errorTotal > numElectronTolerance ){
				statusOFS << "WARNING!!! " 
					<< "Pole expansion may not be accurate enough to reach numElectronTolerance. " << std::endl
					<< "Try to increase numPole or increase numElectronTolerance." << std::endl << std::endl;

			}
			statusOFS << "numPoleInput =" << numPoleInput << std::endl;
			statusOFS << "numPoleSignificant = " << numPoleSignificant << std::endl;
#endif
		}


		// Initialize the number of electrons
		numElectronNow  = 0.0;

		//Initialize the pole expansion
		zweightRhoDrvMu_.resize( numPoleInput );

		GetPoleDensityDrvMu( &zshift_[0], &zweightRhoDrvMu_[0],
				numPoleInput, temperature, gap, deltaE, muNow ); 

		if( isFreeEnergyDensityMatrix ){
			std::vector<Complex>  zshiftTmp( numPoleInput );
			zweightHelmholtz_.resize( numPoleInput );
			GetPoleHelmholtz( &zshiftTmp[0], &zweightHelmholtz_[0], 
				numPoleInput, temperature, gap, deltaE, muNow ); 
		}

		if( isEnergyDensityMatrix ){
			std::vector<Complex>  zshiftTmp( numPoleInput );
			zweightForce_.resize( numPoleInput );
			GetPoleForce( &zshiftTmp[0], &zweightForce_[0],
				numPoleInput, temperature, gap, deltaE, muNow ); 
		}

		if( isDerivativeTMatrix ){
			std::vector<Complex>  zshiftTmp( numPoleInput );
			zweightRhoDrvT_.resize( numPoleInput );
			GetPoleDensityDrvT( &zshiftTmp[0], &zweightRhoDrvT_[0],
				numPoleInput, temperature, gap, deltaE, muNow ); 
		}

#if ( _DEBUGlevel_ >= 1 )
		statusOFS << "zshift" << std::endl << zshift_ << std::endl;
		statusOFS << "zweightRho" << std::endl << zweightRho_ << std::endl;
#endif

		// for each pole, perform LDLT factoriation and selected inversion
		Real timePoleSta, timePoleEnd;

		Int numPoleComputed = 0;
		for(Int lidx = 0; lidx < numPole; lidx++){
			if( MYROW( gridPole_ ) == PROW( lidx, gridPole_ ) ){

				Int l = poleIdx[lidx];

				GetTime( timePoleSta );
#if ( _DEBUGlevel_ >= 0 )
				statusOFS << "Pole " << lidx << " processing..." << std::endl;
				statusOFS << "zshift           = " << zshift_[l] << std::endl;
				statusOFS	<< "zweightRho       = " << zweightRho_[l] << std::endl;
				statusOFS	<< "zweightRhoDrvMu  = " << zweightRhoDrvMu_[l] << std::endl;
				if( isFreeEnergyDensityMatrix )
          statusOFS << "zweightHelmholtz = " << zweightHelmholtz_[l] << std::endl;
				if( isEnergyDensityMatrix )
          statusOFS << "zweightForce     = " << zweightForce_[l] << std::endl;
				if( isDerivativeTMatrix )
					statusOFS << "zweightRhoDrvT   = " << zweightRhoDrvT_[l] << std::endl;
#endif
				{
					numPoleComputed++;

					if( SMat.size != 0 ){
						// S is not an identity matrix
						for( Int i = 0; i < HMat.nnzLocal; i++ ){
							AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift_[l] * SMat.nzvalLocal(i);
						}
					}
					else{
						// S is an identity matrix
						for( Int i = 0; i < HMat.nnzLocal; i++ ){
							AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
						}

						for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
							AMat.nzvalLocal( diagIdxLocal_[i] ) -= zshift_[l];
						}
					} // if (SMat.size != 0 )


					// *********************************************************************
					// Factorization
					// *********************************************************************
					// Important: the distribution in pzsymbfact is going to mess up the
					// A matrix.  Recompute the matrix A here.
#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
#endif
					luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
#endif

					Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

					GetTime( timeTotalFactorizationSta );

					// Data redistribution
#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "Before Distribute." << std::endl;
#endif
					luMat.Distribute();
#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "After Distribute." << std::endl;
#endif

					// Numerical factorization
#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "Before NumericalFactorize." << std::endl;
#endif
					luMat.NumericalFactorize();
#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "After NumericalFactorize." << std::endl;
#endif
					luMat.DestroyAOnly();

					GetTime( timeTotalFactorizationEnd );

#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
#endif

					// *********************************************************************
					// Selected inversion
					// *********************************************************************
					Real timeTotalSelInvSta, timeTotalSelInvEnd;
					GetTime( timeTotalSelInvSta );

					PMatrix<Scalar> PMloc( gridSelInv_, &super_, &luOpt ); // A^{-1} in PMatrix format

					luMat.LUstructToPMatrix( PMloc );
					
          // P2p communication version
          PMloc.ConstructCommunicationPattern();

          // Collective communication version
//          PMloc.ConstructCommunicationPattern_Collectives();

					PMloc.PreSelInv();

					// Main subroutine for selected inversion
          //
          // P2p communication version
          PMloc.SelInv();

          // Collective communication version
//          PMloc.SelInv_Collectives();

					GetTime( timeTotalSelInvEnd );

#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "Time for total selected inversion is " <<
						timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;
#endif

					// *********************************************************************
					// Postprocessing
					// *********************************************************************

					DistSparseMatrix<Complex>  AinvMat;       // A^{-1} in DistSparseMatrix format

					Real timePostProcessingSta, timePostProcessingEnd;

					GetTime( timePostProcessingSta );

					PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );

#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "rhoMat.nnzLocal = " << rhoMat.nnzLocal << std::endl;
					statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
#endif


					// Update the density matrix. The following lines are equivalent to
					//
					//				for( Int i = 0; i < rhoMat.nnzLocal; i++ ){
					//					rhoMat.nzvalLocal(i) += 
					//						zweightRho_[l].real() * AinvMat.nzvalLocal(i).imag() + 
					//						zweightRho_[l].imag() * AinvMat.nzvalLocal(i).real();
					//				}
					// 
					// But done more cache-efficiently with blas.
					Real* AinvMatRealPtr = (Real*)AinvMat.nzvalLocal.Data();
					Real* AinvMatImagPtr = AinvMatRealPtr + 1;
					blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].real(), AinvMatImagPtr, 2, 
							rhoMat.nzvalLocal.Data(), 1 );
					blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].imag(), AinvMatRealPtr, 2,
							rhoMat.nzvalLocal.Data(), 1 );

					// Derivative of the Fermi-Dirac with respect to mu
					blas::Axpy( rhoDrvMuMat.nnzLocal, zweightRhoDrvMu_[l].real(), AinvMatImagPtr, 2, 
							rhoDrvMuMat.nzvalLocal.Data(), 1 );
					blas::Axpy( rhoDrvMuMat.nnzLocal, zweightRhoDrvMu_[l].imag(), AinvMatRealPtr, 2,
							rhoDrvMuMat.nzvalLocal.Data(), 1 );

					// Free energy density matrix
					if( isFreeEnergyDensityMatrix ){
						blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].real(), AinvMatImagPtr, 2,
								hmzMat.nzvalLocal.Data(), 1 );
						blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].imag(), AinvMatRealPtr, 2,
								hmzMat.nzvalLocal.Data(), 1 );
					}

					// Energy density matrix
					if( isEnergyDensityMatrix ){
						blas::Axpy( frcMat.nnzLocal, zweightForce_[l].real(), AinvMatImagPtr, 2,
								frcMat.nzvalLocal.Data(), 1 );
						blas::Axpy( frcMat.nnzLocal, zweightForce_[l].imag(), AinvMatRealPtr, 2, 
								frcMat.nzvalLocal.Data(), 1 );
					}

					// Derivative of the Fermi-Dirac with respect to T
					if( isDerivativeTMatrix ){
						blas::Axpy( rhoDrvTMat.nnzLocal, zweightRhoDrvT_[l].real(), AinvMatImagPtr, 2, 
								rhoDrvTMat.nzvalLocal.Data(), 1 );
						blas::Axpy( rhoDrvTMat.nnzLocal, zweightRhoDrvT_[l].imag(), AinvMatRealPtr, 2,
								rhoDrvTMat.nzvalLocal.Data(), 1 );
					}

					// Update the free energy density matrix and energy density matrix similarly
					GetTime( timePostProcessingEnd );

#if ( _DEBUGlevel_ >= 0 )
					statusOFS << "Time for postprocessing is " <<
						timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
#endif

				}
				GetTime( timePoleEnd );

#if ( _DEBUGlevel_ >= 0 )
				statusOFS << "Time for pole " << lidx << " is " <<
					timePoleEnd - timePoleSta << " [s]" << std::endl << std::endl;
#endif

			} // if I am in charge of this pole
#if ( _DEBUGlevel_ >= 0 )
			// Output the number of electrons at each step for debugging,
			// if there is no parallelization among poles.
			// This debug mode is currently only available if SMat is not
			// implicitly given by  an identity matrix 
			if( gridPole_->numProcRow == 1 && SMat.size != 0 ){
				Real numElecLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
							1, rhoMat_.nzvalLocal.Data(), 1 );

				Real numElec;
				mpi::Allreduce( &numElecLocal, &numElec, 1, MPI_SUM, gridPole_->comm ); 

#if ( _DEBUGlevel_ >= 0 )
				statusOFS << std::endl << "No parallelization of poles, output the number of electrons up to this pole." << std::endl;
				statusOFS << "numElecLocal = " << numElecLocal << std::endl;
				statusOFS << "numElecTotal = " << numElec << std::endl << std::endl;
#endif
				
			}
#endif
		} // for(lidx)

		// Reduce the density matrix across the processor rows in gridPole_
		{
			DblNumVec nzvalRhoMatLocal = rhoMat.nzvalLocal;
			SetValue( rhoMat.nzvalLocal, 0.0 );

			mpi::Allreduce( nzvalRhoMatLocal.Data(), rhoMat.nzvalLocal.Data(),
					rhoMat.nnzLocal, MPI_SUM, gridPole_->colComm );
		}

		// Reduce the derivative of density matrix with respect to mu across
		// the processor rows in gridPole_ 
		{
			DblNumVec nzvalRhoDrvMuMatLocal = rhoDrvMuMat.nzvalLocal;
			SetValue( rhoDrvMuMat.nzvalLocal, 0.0 );

			mpi::Allreduce( nzvalRhoDrvMuMatLocal.Data(), rhoDrvMuMat.nzvalLocal.Data(),
					rhoDrvMuMat.nnzLocal, MPI_SUM, gridPole_->colComm );
		}

		// Reduce the free energy density matrix across the processor rows in gridPole_ 
		if( isFreeEnergyDensityMatrix ){
			DblNumVec nzvalHmzMatLocal = hmzMat.nzvalLocal;
			SetValue( hmzMat.nzvalLocal, 0.0 );

			mpi::Allreduce( nzvalHmzMatLocal.Data(), hmzMat.nzvalLocal.Data(),
					hmzMat.nnzLocal, MPI_SUM, gridPole_->colComm );
		}

		// Reduce the energy density matrix across the processor rows in gridPole_ 
		if( isEnergyDensityMatrix ){
			DblNumVec nzvalFrcMatLocal = frcMat.nzvalLocal;
			SetValue( frcMat.nzvalLocal, 0.0 );

			mpi::Allreduce( nzvalFrcMatLocal.Data(), frcMat.nzvalLocal.Data(),
					frcMat.nnzLocal, MPI_SUM, gridPole_->colComm );
		}

		// Reduce the derivative of density matrix with respect to T across
		// the processor rows in gridPole_ 
		if( isDerivativeTMatrix ){
			DblNumVec nzvalRhoDrvTMatLocal = rhoDrvTMat.nzvalLocal;
			SetValue( rhoDrvTMat.nzvalLocal, 0.0 );

			mpi::Allreduce( nzvalRhoDrvTMatLocal.Data(), rhoDrvTMat.nzvalLocal.Data(),
					rhoDrvTMat.nnzLocal, MPI_SUM, gridPole_->colComm );
		}

		// All processors groups compute the number of electrons, and total
		// energy, and optimally the helmholtz free energy

		numElectronNow = CalculateNumElectron( SMat );
		
		numElectronDrvMuNow = CalculateNumElectronDrvMu( SMat );

    Real totalEnergy = CalculateTotalEnergy( HMat );

#if ( _DEBUGlevel_ >= 0 )
		statusOFS<< "muNow = " << muNow << std::endl; 
		statusOFS<< "numElectronNow = " << numElectronNow << std::endl;
#endif

		Real totalFreeEnergy;
		if( isFreeEnergyDensityMatrix )
			totalFreeEnergy = CalculateFreeEnergy( SMat ) + muNow * numElectronNow;
		Real numElectronDrvT;
		if( isDerivativeTMatrix )
			numElectronDrvT =  CalculateNumElectronDrvT( SMat );

		muList.push_back(muNow);
		numElectronList.push_back( numElectronNow );
		numElectronDrvMuList.push_back( numElectronDrvMuNow );


		GetTime( timeMuEnd );

#if ( _DEBUGlevel_ >= 0 )
		statusOFS << std::endl << "Time for mu iteration " << iter << " is " <<
			timeMuEnd - timeMuSta << " [s]" << std::endl << std::endl;
#endif


#if ( _DEBUGlevel_ >= 0 )
		// Output status 
		statusOFS << std::endl;
		Print( statusOFS, "mu                          = ", muList[iter] );
		Print( statusOFS, "muMin                       = ", muMin );
		Print( statusOFS, "muMax                       = ", muMax ); 
		Print( statusOFS, "Number of poles computed    = ", numPoleComputed );
		Print( statusOFS, "Computed number of electron = ", numElectronList[iter] );
		Print( statusOFS, "d Ne / d mu                 = ", numElectronDrvMuList[iter] );
		Print( statusOFS, "Exact number of electron    = ", numElectronExact );
		Print( statusOFS, "Total energy                = ", totalEnergy );
		if( isFreeEnergyDensityMatrix )
			Print( statusOFS, "Total free energy           = ", totalFreeEnergy );
		if( isDerivativeTMatrix )
			Print( statusOFS, "d Ne / d T                  = ", numElectronDrvT );
#endif

		// Check convergence
		if( std::abs( numElectronExact - numElectronList[iter] ) <
				numElectronTolerance ){
			isConverged = true;
			break;
		}

		// Update the chemical potential

		muNow = CalculateChemicalPotentialNewtonBisection( 
				numElectronExact, numElectronList[iter], 
				numElectronDrvMuList[iter], muList[iter], muMin, muMax );
	
		// Update the bisection interval
		if( numElectronList[iter] < numElectronExact )
			muMin = muList[iter];
		else
			muMax = muList[iter];

	} // for ( iteration of the chemical potential )

	// The final chemical potential for which the number of electrons is
	// NOT computed yet.
	mu = muNow;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::Solve----- 



Real
PPEXSIData::CalculateNumElectron	( const DistSparseMatrix<Real>& SMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateNumElectron");
#endif
	Real numElecLocal = 0.0, numElec = 0.0;
	
	// TODO Check SMat and rhoMat has the same sparsity if SMat is not
	// implicitly given by an identity matrix.

	if( SMat.size != 0 ){
		// S is not an identity matrix
		numElecLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
				1, rhoMat_.nzvalLocal.Data(), 1 );
	}
	else{
		// S is an identity matrix
		DblNumVec& nzval = rhoMat_.nzvalLocal;
		for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
			numElecLocal += nzval(diagIdxLocal_[i]);
		}

	} // if ( SMat.size != 0 )
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "numElecLocal = " << numElecLocal << std::endl;
#endif

	mpi::Allreduce( &numElecLocal, &numElec, 1, MPI_SUM, rhoMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return numElec;
} 		// -----  end of method PPEXSIData::CalculateNumElectron  ----- 

Real
PPEXSIData::CalculateNumElectronDrvMu	( const DistSparseMatrix<Real>& SMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateNumElectronDrvMu");
#endif
	Real numElecDrvLocal = 0.0, numElecDrv = 0.0;
	
	// TODO Check SMat and rhoDrvMuMat has the same sparsity


	if( SMat.size != 0 ){
		// S is not an identity matrix
		numElecDrvLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
				1, rhoDrvMuMat_.nzvalLocal.Data(), 1 );
	}
	else{
		// S is an identity matrix
		DblNumVec& nzval = rhoDrvMuMat_.nzvalLocal;
		for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
			numElecDrvLocal += nzval(diagIdxLocal_[i]);
		}
	}

#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "numElecDrvLocal = " << numElecDrvLocal << std::endl;
#endif

	mpi::Allreduce( &numElecDrvLocal, &numElecDrv, 1, MPI_SUM, rhoDrvMuMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return numElecDrv;
} 		// -----  end of method PPEXSIData::CalculateNumElectronDrvMu  ----- 

Real
PPEXSIData::CalculateNumElectronDrvT	( const DistSparseMatrix<Real>& SMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateNumElectronDrvT");
#endif
	Real numElecDrvLocal = 0.0, numElecDrv = 0.0;
	
	// TODO Check SMat and rhoDrvTMat has the same sparsity

	if( SMat.size != 0 ){
		// S is not an identity matrix
		numElecDrvLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
				1, rhoDrvTMat_.nzvalLocal.Data(), 1 );
	}
	else{
		// S is an identity matrix
		DblNumVec& nzval = rhoDrvTMat_.nzvalLocal;
		for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
			numElecDrvLocal += nzval(diagIdxLocal_[i]);
		}
	}
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "numElecDrvLocal = " << numElecDrvLocal << std::endl;
#endif

	mpi::Allreduce( &numElecDrvLocal, &numElecDrv, 1, MPI_SUM, rhoDrvTMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return numElecDrv;
} 		// -----  end of method PPEXSIData::CalculateNumElectronDrvT  ----- 


Real
PPEXSIData::CalculateTotalEnergy	( const DistSparseMatrix<Real>& HMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateTotalEnergy");
#endif
	
	Real totalEnergyLocal = 0.0, totalEnergy = 0.0;
	
	// TODO Check HMat and rhoMat has the same sparsity

	totalEnergyLocal = blas::Dot( HMat.nnzLocal, HMat.nzvalLocal.Data(),
			1, rhoMat_.nzvalLocal.Data(), 1 );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "TotalEnergyLocal = " << totalEnergyLocal << std::endl;
#endif

	mpi::Allreduce( &totalEnergyLocal, &totalEnergy, 1, MPI_SUM, rhoMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return totalEnergy;
} 		// -----  end of method PPEXSIData::CalculateTotalEnergy  ----- 

Real 
PPEXSIData::CalculateFreeEnergy	( const DistSparseMatrix<Real>& SMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateFreeEnergy");
#endif
	
	Real totalFreeEnergyLocal = 0.0, totalFreeEnergy = 0.0;

	// TODO Check SMat and freeEnergyDensityMat_ has the same sparsity

	if( SMat.size != 0 ){
		// S is not an identity matrix
		totalFreeEnergyLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
				1, freeEnergyDensityMat_.nzvalLocal.Data(), 1 );
	}
	else{
		// S is an identity matrix
		DblNumVec& nzval = freeEnergyDensityMat_.nzvalLocal;
		for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
		  totalFreeEnergyLocal += nzval(diagIdxLocal_[i]);
		}
	} // if ( SMat.size != 0 )

#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "TotalFreeEnergyLocal = " << totalFreeEnergyLocal << std::endl;
#endif

	mpi::Allreduce( &totalFreeEnergyLocal, &totalFreeEnergy, 1, MPI_SUM, 
			freeEnergyDensityMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return totalFreeEnergy;
} 		// -----  end of method PPEXSIData::CalculateFreeEnergy  ----- 


Real
PPEXSIData::CalculateForce	( 
		const DistSparseMatrix<Real>& HDerivativeMat,  
		const DistSparseMatrix<Real>& SDerivativeMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateForce");
#endif

	Real totalForceLocal = 0.0, totalForce = 0.0;

	// TODO Check HDerivativeMat, SDerivativeMat, rhoMat and
	// energyDensityMat_ has the same sparsity pattern

	totalForceLocal = - blas::Dot( HDerivativeMat.nnzLocal,
			HDerivativeMat.nzvalLocal.Data(), 1, rhoMat_.nzvalLocal.Data(), 1
			);

	if( SDerivativeMat.size != 0 ){
		// If S is not an identity matrix, compute the Pulay force
		totalForceLocal += blas::Dot( SDerivativeMat.nnzLocal,
				SDerivativeMat.nzvalLocal.Data(), 1,
				energyDensityMat_.nzvalLocal.Data(), 1 );
	}

#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "TotalForceLocal = " << totalForceLocal << std::endl;
#endif

	mpi::Allreduce( &totalForceLocal, &totalForce, 1, MPI_SUM, rhoMat_.comm );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return totalForce;
} 		// -----  end of method PPEXSIData::CalculateForce  ----- 



Real
PPEXSIData::EstimateZeroTemperatureChemicalPotential	( 
		Real temperature,
		Real mu,
	 	const DistSparseMatrix<Real>& SMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::EstimateZeroTemperatureChemicalPotential");
#endif
 
  Real numElecDrvMu, numElecDrvT, muDrvT;
  Real beta;
	Real mu0;

	beta = 1.0 / (temperature);

  numElecDrvMu = CalculateNumElectronDrvMu( SMat );
  
  numElecDrvT  = CalculateNumElectronDrvT ( SMat );

	muDrvT       = -numElecDrvT / numElecDrvMu;

	mu0 = mu - 0.5 * muDrvT / beta;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return mu0;
} 		// -----  end of method PPEXSIData::EstimateZeroTemperatureChemicalPotential  ----- 

void PPEXSIData::CalculateNegativeInertia( 
		const std::vector<Real>&       shiftVec, 
		std::vector<Real>&             inertiaVec,
		const DistSparseMatrix<Real>&  HMat,
		const DistSparseMatrix<Real>&  SMat,
		std::string                    ColPerm,
		Int                            numProcSymbFact
		){
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateNegativeInertia");
#endif



	// *********************************************************************
	// Initialize
	// *********************************************************************

	DistSparseMatrix<Complex>  AMat;              // A = H - z * S

	// Get the diagonal indices for H and save it n diagIdxLocal_
	{
		Int numColLocal      = HMat.colptrLocal.m() - 1;
		Int numColLocalFirst = HMat.size / gridSelInv_->mpisize;
		Int firstCol         = gridSelInv_->mpirank * numColLocalFirst;
		
		diagIdxLocal_.clear();

		for( Int j = 0; j < numColLocal; j++ ){
			Int jcol = firstCol + j + 1;
			for( Int i = HMat.colptrLocal(j)-1; 
				 	 i < HMat.colptrLocal(j+1)-1; i++ ){
				Int irow = HMat.rowindLocal(i);
				if( irow == jcol ){
					diagIdxLocal_.push_back( i );
				}
			}
		} // for (j)
	}


	// Copy the pattern
	CopyPattern( HMat, AMat );

	SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value


	SuperLUOptions   luOpt;

	luOpt.ColPerm = ColPerm;
	luOpt.numProcSymbFact = numProcSymbFact;

	SuperLUMatrix<Complex>    luMat( *zGridSuperLU_, luOpt );  // SuperLU matrix.

	// *********************************************************************
	// Symbolic factorization.  
	// Each numPoleGroup perform independently
	// *********************************************************************
	Real timeSta, timeEnd;
	GetTime( timeSta );
	luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
	GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "AMat is converted to SuperMatrix." << std::endl;
	statusOFS << "Time for SuperMatrix conversion is " <<
		timeEnd - timeSta << " [s]" << std::endl << std::endl;
#endif
	GetTime( timeSta );
	luMat.SymbolicFactorize();
	GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "Symbolic factorization is finished." << std::endl;
	statusOFS << "Time for symbolic factorization is " <<
		timeEnd - timeSta << " [s]" << std::endl << std::endl;
#endif
	luMat.SymbolicToSuperNode( super_ );
	luMat.DestroyAOnly();


	// Compute the number of nonzeros from PMatrix
	{
		PMatrix<Scalar> PMloc( gridSelInv_, &super_ , &luOpt); // A^{-1} in PMatrix format
		luMat.LUstructToPMatrix( PMloc );
#if ( _DEBUGlevel_ >= 0 )
		Int nnzLocal = PMloc.NnzLocal();
		statusOFS << "Number of local nonzeros (L+U) = " << nnzLocal << std::endl;
		LongInt nnz  = PMloc.Nnz();
		statusOFS << "Number of nonzeros (L+U)       = " << nnz << std::endl;
#endif
	}

#if ( _DEBUGlevel_ >= 1 )
	statusOFS << "perm: "    << std::endl << super_.perm     << std::endl;
	statusOFS << "permInv: " << std::endl << super_.permInv  << std::endl;
	statusOFS << "superIdx:" << std::endl << super_.superIdx << std::endl;
	statusOFS << "superPtr:" << std::endl << super_.superPtr << std::endl; 
#endif


	Real timeShiftSta, timeShiftEnd;
	
	Int numShift = shiftVec.size();
	std::vector<Real>  inertiaVecLocal(numShift);
	inertiaVec.resize(numShift);
	for(Int l = 0; l < numShift; l++){
		inertiaVecLocal[l] = 0.0;
		inertiaVec[l]      = 0.0;
	}

	for(Int l = 0; l < numShift; l++){
		if( MYROW( gridPole_ ) == PROW( l, gridPole_ ) ){

			GetTime( timeShiftSta );

#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "Shift " << l << " = " << shiftVec[l] 
				<< " processing..." << std::endl;
#endif

			if( SMat.size != 0 ){
				// S is not an identity matrix
				for( Int i = 0; i < HMat.nnzLocal; i++ ){
					AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - shiftVec[l] * SMat.nzvalLocal(i);
				}
			}
			else{
				// S is an identity matrix
				for( Int i = 0; i < HMat.nnzLocal; i++ ){
					AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
				}

				for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
					AMat.nzvalLocal( diagIdxLocal_[i] ) -= shiftVec[l];
				}
			} // if (SMat.size != 0 )


			// *********************************************************************
			// Factorization
			// *********************************************************************
			// Important: the distribution in pzsymbfact is going to mess up the
			// A matrix.  Recompute the matrix A here.
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "Before DistSparseMatrixToSuperMatrixNRloc." << std::endl;
#endif
			luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "After DistSparseMatrixToSuperMatrixNRloc." << std::endl;
#endif

			Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

			GetTime( timeTotalFactorizationSta );

			// Data redistribution
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "Before Distribute." << std::endl;
#endif
			luMat.Distribute();
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "After Distribute." << std::endl;
#endif

			// Numerical factorization
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "Before NumericalFactorize." << std::endl;
#endif
			luMat.NumericalFactorize();
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "After NumericalFactorize." << std::endl;
#endif
			luMat.DestroyAOnly();

			GetTime( timeTotalFactorizationEnd );

#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 
#endif

			// *********************************************************************
			// Compute inertia
			// *********************************************************************
			Real timeInertiaSta, timeInertiaEnd;
			GetTime( timeInertiaSta );

			PMatrix<Scalar> PMloc( gridSelInv_, &super_, &luOpt ); // A^{-1} in PMatrix format

			luMat.LUstructToPMatrix( PMloc );

			// Compute the negative inertia of the matrix.
			PMloc.GetNegativeInertia( inertiaVecLocal[l] );

			GetTime( timeInertiaEnd );

#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "Time for computing the inertia is " <<
				timeInertiaEnd  - timeInertiaSta << " [s]" << std::endl;
#endif


		} // if I am in charge of this shift
	} // for(l)

	// Collect all the negative inertia together
	mpi::Allreduce( &inertiaVecLocal[0], &inertiaVec[0], numShift, 
			MPI_SUM, gridPole_->colComm );


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::CalculateNegativeInertia ----- 


void
PPEXSIData::EstimateSpectralRadius(
    Int                            method,
    const DistSparseMatrix<Real>&  HMat,
    const DistSparseMatrix<Real>&  SMat,
    const NumVec<Real>&            v0,
    Real                           tol,
    Int                            maxNumIter,
    Int&                           numIter,
    Real&                          sigma ){
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::EstimateSpectralRadius");
#endif
  Int n = HMat.size;

  bool isSIdentity = ( SMat.size == 0 );

  if( HMat.size != SMat.size && (!isSIdentity) ){
		std::ostringstream msg;
		msg << std::endl
			<< "The sizes of H and S do not match." << std::endl
      << "H.size = " << HMat.size << std::endl
      << "S.size = " << SMat.size << std::endl;
		#ifdef USE_ABORT
abort();
#endif
throw std::runtime_error( msg.str().c_str() );
  }

  // Locally Optimal Conjugate Gradient method
  if( method == 0 ){
    // Q = [X W P]
    // X: The current vector
    // W: (Preconditioned) residual
    // P: Conjugate direction
    // A: The same as H
    // B: The same as S
    DblNumMat  Q, AQ, BQ; 
    Q.Resize( n, 3 );
    AQ.Resize( n, 3 );
    BQ.Resize( n, 3 );
    
    DblNumVec X( n, false, Q.VecData(0) );
    DblNumVec W( n, false, Q.VecData(1) );
    DblNumVec P( n, false, Q.VecData(2) );
    DblNumVec AX( n, false, AQ.VecData(0) );
    DblNumVec AW( n, false, AQ.VecData(1) );
    DblNumVec AP( n, false, AQ.VecData(2) );
    DblNumVec BX( n, false, BQ.VecData(0) );
    DblNumVec BW( n, false, BQ.VecData(1) );
    DblNumVec BP( n, false, BQ.VecData(2) );

    // Residual
    DblNumVec R( n );

    if( v0.m() != 0 ){
      // Use v0
      if( v0.m() != HMat.size ){
        std::ostringstream msg;
        msg << std::endl
          << "The sizes of H and v0 do not match." << std::endl
          << "H.size = " << HMat.size << std::endl
          << "v0.size = " << v0.m() << std::endl;
        #ifdef USE_ABORT
abort();
#endif
throw std::runtime_error( msg.str().c_str() );
      }
      blas::Copy( n, v0.Data(), 1, X.Data(), 1 );
    }
    else{
      // Start from random initial vector
      UniformRandom( X );
    }

    Real XNorm = std::sqrt( Energy(X) );
    blas::Scal( n, 1.0 / XNorm, X.Data(), 1 );

    Int numConv = 0;
    Int sizeQ;

    Real lambda = 0.0, lambdaP = 0.0;
    Real RNorm, PNorm;
    DblNumMat QAQ(3, 3), QBQ(3, 3);
    DblNumVec ev(3);
    DblNumVec ef(3);
    SetValue( QAQ, 0.0 );
    SetValue( QBQ, 0.0 );
    SetValue( ev, 0.0 );
    SetValue( ef, 0.0 );

    Int iter;
    for( iter = 1; iter <= maxNumIter; iter++ ){
      DistSparseMatMultGlobalVec( 
          1.0, HMat, X, 0.0, AX );
      if( !isSIdentity ){
        DistSparseMatMultGlobalVec( 
            1.0, SMat, X, 0.0, BX );
      }
      else{
        blas::Copy( n, X.Data(), 1, BX.Data(), 1 );
      }
      
      QAQ(0,0) = blas::Dot( n, X.Data(), 1, AX.Data(), 1 );
      QBQ(0,0) = blas::Dot( n, X.Data(), 1, BX.Data(), 1 );

      lambdaP = lambda;
      lambda = QAQ(0,0) / QBQ(0,0);
      // Residual
      for( Int i = 0; i < n; i++ )
        R(i) = AX(i) - BX(i) *lambda;
      
      RNorm = std::sqrt( Energy(R) );

#if ( _DEBUGlevel_ >= 2 )
      statusOFS << "R = " << R << std::endl;
      statusOFS << "X = " << X << std::endl;
      statusOFS << "Energy(R) = " << Energy(R) << std::endl;
      statusOFS << "RNorm = " << RNorm << std::endl;
#endif
      
      if( RNorm < tol ){
        numConv = 1;
      }
      if( iter > 1 ){
        if( std::abs( lambdaP - lambda ) < tol * std::abs( lambda ) ){
          numConv = 1;
        }
      }

      if( numConv >= 1 ){
        break;
      }


      // No preconditioner applied to the residual
      for( Int i = 0; i < n; i++ )
        W(i) = R(i) / RNorm;

      
      DistSparseMatMultGlobalVec( 
          1.0, HMat, W, 0.0, AW );
      if( !isSIdentity ){
        DistSparseMatMultGlobalVec( 
            1.0, SMat, W, 0.0, BW );
      }
      else{
        blas::Copy( n, W.Data(), 1, BW.Data(), 1 );
      }

      if( iter == 1 )
        sizeQ = 2;
      else
        sizeQ = 3;

      blas::Gemm( 'T', 'N', sizeQ, sizeQ, n, 
          1.0, Q.Data(), n, AQ.Data(), n, 0.0, QAQ.Data(), 3 );

      blas::Gemm( 'T', 'N', sizeQ, sizeQ, n, 
          1.0, Q.Data(), n, BQ.Data(), n, 0.0, QBQ.Data(), 3 );

#if ( _DEBUGlevel_ >= 2 )
      statusOFS << "QAQ = " << QAQ << std::endl;
      statusOFS << "QBQ = " << QBQ << std::endl;
#endif

      // ( Q^T A Q ) V = ( Q^T B Q ) V D
      lapack::Sygvd( 1, 'V', 'L', sizeQ, QAQ.Data(), 3, 
          QBQ.Data(), 3, ev.Data() );

      // ef is the eigenvector corresponding to the largest eigenvalue
      blas::Copy( sizeQ, QAQ.VecData(sizeQ-1), 1, ef.Data(), 1 );
      
#if ( _DEBUGlevel_ >= 2 )
      statusOFS << "ev = " << ev << std::endl;
      statusOFS << "ef = " << ef << std::endl;
#endif


      if( iter > 1 ){
        //  X =  Q * ef
        // AX = AQ * ef
        // BX = BQ * ef
        //  P =  W * ef[1] +  P * ef[2]
        // AP = AW * ef[1] + AP * ef[2]
        // BP = BW * ef[1] + BP * ef[2]

        for( Int i = 0; i < n; i++ ){
          X[i] = ef[0] * X[i] + ef[1] * W[i] + ef[2] * P[i];
          AX[i] = ef[0] * AX[i] + ef[1] * AW[i] + ef[2] * AP[i]; 
          BX[i] = ef[0] * BX[i] + ef[1] * BW[i] + ef[2] * BP[i];

          P[i] = ef[1] * W[i] + ef[2] * P[i];
          AP[i] = ef[1] * AW[i] + ef[2] * AP[i]; 
          BP[i] = ef[1] * BW[i] + ef[2] * BP[i];
        }
        PNorm = std::sqrt( Energy( P ) );
        blas::Scal( n, 1.0 / PNorm, P.Data(), 1 );
        blas::Scal( n, 1.0 / PNorm, AP.Data(), 1 );
        blas::Scal( n, 1.0 / PNorm, BP.Data(), 1 );
      }
      else{
        //  X =  Q * ef
        // AX = AQ * ef
        // BX = BQ * ef
        // P  = W;
        // AP = AW;
        // BP = BW;
        for( Int i = 0; i < n; i++ ){
          X[i] = ef[0] * X[i] + ef[1] * W[i];
          AX[i] = ef[0] * AX[i] + ef[1] * AW[i]; 
          BX[i] = ef[0] * BX[i] + ef[1] * BW[i];

          P[i] = W[i];
          AP[i] = AW[i]; 
          BP[i] = BW[i];
        }
      }
    } // for (iter)

    // Final computation of the eigenvalue
    DistSparseMatMultGlobalVec( 
        1.0, HMat, X, 0.0, AX );
    if( !isSIdentity ){
      DistSparseMatMultGlobalVec( 
          1.0, SMat, X, 0.0, BX );
    }
    else{
      blas::Copy( n, X.Data(), 1, BX.Data(), 1 );
    }

    QAQ(0,0) = blas::Dot( n, X.Data(), 1, AX.Data(), 1 );
    QBQ(0,0) = blas::Dot( n, X.Data(), 1, BX.Data(), 1 );

    // Output variable
    numIter = (iter <= maxNumIter) ? iter : maxNumIter;
    sigma = QAQ(0,0) / QBQ(0,0);
  }
  else{
    #ifdef USE_ABORT
abort();
#endif
throw std::runtime_error( "Method != 0 is not yet supported." );
  }


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function PPEXSIData::EstimateSpectralRadius  ----- 
} //  namespace PEXSI
