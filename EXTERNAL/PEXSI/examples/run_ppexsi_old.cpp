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
/// @file run_ppexsi.cpp
/// @brief Test for the PPEXSI module using SuperLU and PSelInv.
/// @date 2012-12-24
#include "ppexsi.hpp"

using namespace PEXSI;
using namespace std;

// FIXME
extern "C" { 
	double seekeig_(int *, int *, double *, double *, double *);
}


void Usage(){
  std::cout 
    << "WARNING: This file is out of date.  Do not use it for any serious purpose!" << std::endl
		<< "run_ppexsi -temp [temp] -mu0 [mu0] -muMin [muMin] -muMax [muMax] -numel [numel] -numPole [numPole] -deltaE [deltaE] -gap [gap] -H [Hfile] -S [Sfile] -npPerPole [npPole] -npSymbFact [npsymbfact] -colperm [colperm] -muiter [muiter]" << std::endl 
		<< "temp:    Temperature (unit: K)" << std::endl
		<< "mu0:     Initial guess for chemical potential" << std::endl
		<< "muMin:   Lower bound for chemical potential" << std::endl
		<< "muMax:   Upper bound for chemical potential" << std::endl
		<< "numel:   Exact number of electrons (spin-restricted)" << std::endl
		<< "numPole: Number of poles." << std::endl
		<< "deltaE:  guess for the width of the spectrum of H-mu S" << std::endl
		<< "gap:     guess for the distance betweeen the spectrum and mu" << std::endl
		<< "H: Hamiltonian matrix " << std::endl
		<< "S: Overlap     matrix. if omitted, the overlap matrix is treated as an identity matrix implicitly." << std::endl
		<< "npPerPole: number of processors used for each pole" << std::endl
		<< "npSymbFact: number of processors used for parallel symbolic factorization" << std::endl
		<< "colperm: permutation method (for SuperLU_DIST)" << std::endl
	  << "muiter:  number of iterations for the chemical potential" << std::endl
		<< "formatted: whether the input of H/S matrices are formatted (1) or unformatted (csc format, 0)" << std::endl;
}

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


	if( argc < 25 || argc%2 == 0 ) {
		if( mpirank == 0 ) Usage();
		MPI_Finalize();
		return 0;
	}
			

	
	try{
		stringstream  ss;
		ss << "logPEXSI" << mpirank;
		statusOFS.open( ss.str().c_str() );


		// *********************************************************************
		// Input parameter
		// *********************************************************************
		std::map<std::string,std::string> options;
		OptionsCreate(argc, argv, options);
		
		Real numElectronTolerance = 1e-2;

		Real temperature;
		if( options.find("-temp") != options.end() ){
			temperature = std::atof(options["-temp"].c_str());
		}
		else{
      throw std::logic_error("temp must be provided.");
		}

		Real mu0;
		if( options.find("-mu0") != options.end() ){
			mu0 = std::atof(options["-mu0"].c_str());
		}
		else{
      throw std::logic_error("mu0 must be provided.");
		}

		Real muMin;
		if( options.find("-muMin") != options.end() ){
			muMin = std::atof(options["-muMin"].c_str());
		}
		else{
      throw std::logic_error("muMin must be provided.");
		}

		Real muMax;
		if( options.find("-muMax") != options.end() ){
			muMax = std::atof(options["-muMax"].c_str());
		}
		else{
      throw std::logic_error("muMax must be provided.");
		}

		Real numElectronExact;
    if( options.find("-numel") != options.end() ){
			numElectronExact = std::atof(options["-numel"].c_str());
		}
		else{
      throw std::logic_error("numel must be provided.");
		}
		
		Int numPole;
    if( options.find("-numPole") != options.end() ){
			numPole = std::atof(options["-numPole"].c_str());
		}
		else{
      throw std::logic_error("numPole must be provided.");
		}


		Real deltaE;
    if( options.find("-deltaE") != options.end() ){
			deltaE = std::atof(options["-deltaE"].c_str());
		}
		else{
      throw std::logic_error("deltaE must be provided.");
		}

		Real gap;
    if( options.find("-gap") != options.end() ){
			gap = std::atof(options["-gap"].c_str());
		}
		else{
      throw std::logic_error("gap must be provided.");
		}


		Int npPerPole;
    if( options.find("-npPerPole") != options.end() ){
			npPerPole = std::atoi(options["-npPerPole"].c_str());
		}
		else{
      throw std::logic_error("npPerPole must be provided.");
		}

		Int npSymbFact;
    if( options.find("-npSymbFact") != options.end() ){
			npSymbFact = std::atoi(options["-npSymbFact"].c_str());
		}
		else{
      throw std::logic_error("npSymbFact must be provided.");
		}

		Int isFormatted;
		if( options.find("-formatted") != options.end() ){
			isFormatted = std::atoi(options["-formatted"].c_str());
		}
		else{
			isFormatted = 0; // Binary file
		}
   
		std::string Hfile, Sfile;                   
		if( options.find("-H") != options.end() ){ 
			Hfile = options["-H"];
		}
		else{
      throw std::logic_error("Hfile must be provided.");
		}

		if( options.find("-S") != options.end() ){ 
			Sfile = options["-S"];
		}
		else{
			statusOFS << "-S option is not given. " 
				<< "Treat the overlap matrix as an identity matrix." 
				<< std::endl << std::endl;
		}

		std::string ColPerm;
		if( options.find("-colperm") != options.end() ){ 
			ColPerm = options["-colperm"];
		}
		else{
      throw std::logic_error("colperm must be provided.");
		}

		bool isFreeEnergyDensityMatrix = true;

		bool isEnergyDensityMatrix     = true;

		bool isDerivativeTMatrix       = true;

		bool isInertiaCount            = true;

		// *********************************************************************
		// Check the input parameters
		// *********************************************************************
		if( mpisize % npPerPole != 0 ){
			throw std::logic_error( "mpisize cannot be divided evenly by npPerPole." );
		}

	  Int nprow = iround( std::sqrt( (double)npPerPole) );
		Int npcol = npPerPole / nprow;
		if( npPerPole != nprow * npcol || nprow != npcol ){
			throw std::runtime_error( "npPerPole must be a perfect square due to the current implementation of PSelInv." );
		}

		Int muMaxIter;
    if( options.find("-muiter") != options.end() ){
			muMaxIter = std::atoi(options["-muiter"].c_str());
		}
		else{
      throw std::logic_error("muiter must be provided.");
		}


		// Initialize
		GridType gridPole( MPI_COMM_WORLD, mpisize / npPerPole, npPerPole );
		PPEXSIData pexsi( &gridPole, nprow, npcol );

		// *********************************************************************
		// Read input matrix
		// *********************************************************************


		DistSparseMatrix<Real> HMat;

		Real timeSta, timeEnd;

		// The first processor row gridPole reads first, and then
		// broadcast the information to other row processors in gridPole.
		// This can be improved later by MPI_Bcast directly.

		// HMat
		std::vector<char> sstr;
		Int sizeStm;
		if( MYROW( &gridPole ) == 0 ){
			std::stringstream sstm;


			if( isFormatted )
				ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, gridPole.rowComm ); 
			else
				ReadDistSparseMatrix( Hfile.c_str(), HMat, gridPole.rowComm ); 

			serialize( HMat.size, sstm, NO_MASK );
			serialize( HMat.nnz,  sstm, NO_MASK );
			serialize( HMat.nnzLocal, sstm, NO_MASK );
			serialize( HMat.colptrLocal, sstm, NO_MASK );
			serialize( HMat.rowindLocal, sstm, NO_MASK );
			serialize( HMat.nzvalLocal,  sstm, NO_MASK );
			sstr.resize( Size( sstm ) );
		  sstm.read( &sstr[0], sstr.size() ); 	
			sizeStm = sstr.size();
		}


		MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );
		statusOFS << "sizeStm = " << sizeStm << std::endl;

		if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

		MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

		if( MYROW( &gridPole ) != 0 ){
			std::stringstream sstm;
			sstm.write( &sstr[0], sizeStm );
			deserialize( HMat.size, sstm, NO_MASK );
			deserialize( HMat.nnz,  sstm, NO_MASK );
			deserialize( HMat.nnzLocal, sstm, NO_MASK );
			deserialize( HMat.colptrLocal, sstm, NO_MASK );
			deserialize( HMat.rowindLocal, sstm, NO_MASK );
			deserialize( HMat.nzvalLocal,  sstm, NO_MASK );
		}
		statusOFS << "size     = " << HMat.size     << std::endl;
		statusOFS << "nnzLocal = " << HMat.nnzLocal << std::endl;
		statusOFS << "nnz      = " << HMat.Nnz()    << std::endl;
		// Communicator
		HMat.comm = gridPole.rowComm;

		sstr.clear();

		// SMat
		DistSparseMatrix<Real> SMat;
		if( Sfile.empty() ){
			// Set the size to be zero.  This will tell PPEXSI.Solve to treat
			// the overlap matrix as an identity matrix implicitly.
			SMat.size = 0;  
		}
		else{
			// SMat is given directly
			if( MYROW( &gridPole ) == 0 ){
				std::stringstream sstm;
				if( isFormatted )
					ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, gridPole.rowComm ); 
				else
					ReadDistSparseMatrix( Sfile.c_str(), SMat, gridPole.rowComm ); 

				serialize( SMat.size, sstm, NO_MASK );
				serialize( SMat.nnz,  sstm, NO_MASK );
				serialize( SMat.nnzLocal, sstm, NO_MASK );
				serialize( SMat.colptrLocal, sstm, NO_MASK );
				serialize( SMat.rowindLocal, sstm, NO_MASK );
				serialize( SMat.nzvalLocal,  sstm, NO_MASK );
				sstr.resize( Size( sstm ) );
				sstm.read( &sstr[0], sstr.size() ); 	
				sizeStm = sstr.size();
			}

			MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );

			if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

			MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

			if( MYROW( &gridPole ) != 0 ){
				std::stringstream sstm;
				sstm.write( &sstr[0], sizeStm );
				deserialize( SMat.size, sstm, NO_MASK );
				deserialize( SMat.nnz,  sstm, NO_MASK );
				deserialize( SMat.nnzLocal, sstm, NO_MASK );
				deserialize( SMat.colptrLocal, sstm, NO_MASK );
				deserialize( SMat.rowindLocal, sstm, NO_MASK );
				deserialize( SMat.nzvalLocal,  sstm, NO_MASK );
			}
			// Communicator
			SMat.comm = gridPole.rowComm;

			sstr.clear();
		} // if (Sfile.empty())


		Print(statusOFS, "mu0                    = ", mu0);
		Print(statusOFS, "muMin                  = ", muMin);
		Print(statusOFS, "muMax                  = ", muMax); 
		Print(statusOFS, "numElectronExact       = ", numElectronExact);
		Print(statusOFS, "deltaE                 = ", deltaE);
		Print(statusOFS, "gap                    = ", gap);
		Print(statusOFS, "temperature            = ", temperature);
		Print(statusOFS, "numPole                = ", numPole);
		Print(statusOFS, "numElectronTolerance   = ", numElectronTolerance);
		Print(statusOFS, "ColPerm                = ", ColPerm );
		Print(statusOFS, "muMaxIter              = ", muMaxIter);
		Print(statusOFS, "mpisize                = ", mpisize );
		Print(statusOFS, "npPerPole              = ", npPerPole );
		Print(statusOFS, "npSymbFact              = ", npSymbFact );
		Print(statusOFS, "isFreeEnergyMatrix     = ", isFreeEnergyDensityMatrix );
		Print(statusOFS, "isEnergyMatrix         = ", isEnergyDensityMatrix ); 
		Print(statusOFS, "isInertiaCount         = ", isInertiaCount ); 


		// *********************************************************************
		// Solve
		// *********************************************************************
		std::vector<Real>  muList;
		std::vector<Real>  numElectronList;
		std::vector<Real>  numElectronDrvList;
	  bool isConverged;	
		
		Real timeSolveSta, timeSolveEnd;

		Real mu;
		if( isInertiaCount ){
			// Inertia count overwrites the initial input mu, as well as the
			// bound muMin and muMax.

			// It keeps refining till the accuracy of the number of electrons
			// becomes less than 1.

			// Number of inertia counts is the same as the number of poles
			Int  numShift = numPole;
			std::vector<Real>  shiftVec( numShift );
			std::vector<Real>  inertiaVec( numShift );

			for( Int l = 0; l < numShift; l++ ){
				shiftVec[l] = muMin + l * (muMax - muMin) / (numShift-1);
			}

			GetTime( timeSta );
			pexsi.CalculateNegativeInertiaReal( 
					shiftVec,
					inertiaVec,
					HMat,
					SMat,
					ColPerm,
				  npSymbFact );

			GetTime( timeEnd );

			for( Int l = 0; l < numShift; l++ ){
				// Inertia is multiplied by 2.0 to reflect the doubly occupied
				// orbitals.
				inertiaVec[l] *= 2.0;
				statusOFS << std::setiosflags(std::ios::left) 
					<< std::setw(LENGTH_VAR_NAME) << "Shift = "
					<< std::setw(LENGTH_VAR_DATA) << shiftVec[l]
					<< std::setw(LENGTH_VAR_NAME) << "Inertia = "
					<< std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
					<< std::endl << std::endl;
			}

			Print( statusOFS, "Time for inertia count = ", 
					timeEnd - timeSta );
			{
				double muStart = (muMin  + muMax)/2.0;
				int nelec = numElectronExact;
				mu = seekeig_(&nelec, &numShift, &shiftVec[0],&inertiaVec[0], &muStart); 
			}

		// Just for numerical stability
		const Real EPS = 1e-6;
		std::vector<Real>::iterator vi0, vi1;
		vi0 = std::lower_bound( inertiaVec.begin(), inertiaVec.end(), 
				numElectronExact-EPS );
		vi1 = std::upper_bound( inertiaVec.begin(), inertiaVec.end(), 
				numElectronExact+EPS );
		if( vi0 == inertiaVec.begin() || vi0 == inertiaVec.end() ){
			throw std::runtime_error("Increase the range of [muMin, muMax].");
		}
		if( vi1 == inertiaVec.begin() || vi1 == inertiaVec.end() ){
			throw std::runtime_error("Increase the range of [muMin, muMax].");
		}
		Int idx0 = vi0 - inertiaVec.begin() - 1;
		Int idx1 = vi1 - inertiaVec.begin();
		// Adjust muMin, muMax by a safer bound taking into account the
		// temperature effect.  
		muMin = shiftVec[idx0] - temperature / au2K;
		muMax = shiftVec[idx1] + temperature / au2K;

			statusOFS << std::endl << "After the inertia count," << std::endl;
			Print( statusOFS, "mu      = ", mu );
			Print( statusOFS, "muMin   = ", muMin );
			Print( statusOFS, "muMax   = ", muMax );
			statusOFS << std::endl;

		}
		else{
			// Otherwise, use the initial input mu
			mu = mu0;
		}

		GetTime( timeSolveSta );
		pexsi.Solve( 
				numPole,
				temperature,
				numElectronExact,
				gap,
				deltaE,
				mu,
				muMin,
				muMax,
				HMat,
				SMat,
				muMaxIter,
				numElectronTolerance,
				ColPerm,
				npSymbFact,
				isFreeEnergyDensityMatrix,
				isEnergyDensityMatrix,
				isDerivativeTMatrix,
				muList,
				numElectronList,
				numElectronDrvList,
			  isConverged	);

		GetTime( timeSolveEnd );

		PrintBlock( statusOFS, "Solve finished." );
		if( isConverged ){
			statusOFS << "PEXSI has converged with " << muList.size() << 
				" iterations" << std::endl;
		}
		else {
			statusOFS << "PEXSI did not converge with " << muList.size() << 
				" iterations" << std::endl;
		}
		Print( statusOFS, "mu                   = ", mu );
		Print( statusOFS, "muMin                = ", muMin ); 
		Print( statusOFS, "muMax                = ", muMax ); 
		Print( statusOFS, "Total time for PEXSI = ", 
				timeSolveEnd - timeSolveSta );


		
		// *********************************************************************
		// Solve for the second
		// Using temperature expansion for the chemical potential
		// *********************************************************************

		if(0){
			Real muZeroT = pexsi.EstimateZeroTemperatureChemicalPotential(
					temperature,
					*muList.rbegin(),
					SMat );

			PrintBlock( statusOFS, "Second calculation at low temperature." );

			Print( statusOFS, "mu (T=0 estimate)    = ", 
					muZeroT );

			mu0   = muZeroT;
			// Safe choice of bound
			muMin = mu0 - temperature / au2K;
			muMax = mu0 + temperature / au2K;

			// New temperature
			temperature = 300.0;


			Print(statusOFS, "mu0                    = ", mu0);
			Print(statusOFS, "muMin                  = ", muMin);
			Print(statusOFS, "muMax                  = ", muMax); 
			Print(statusOFS, "numElectronExact       = ", numElectronExact);
			Print(statusOFS, "deltaE                 = ", deltaE);
			Print(statusOFS, "gap                    = ", gap);
			Print(statusOFS, "temperature            = ", temperature);
			Print(statusOFS, "numPole                = ", numPole);
			Print(statusOFS, "numElectronTolerance   = ", numElectronTolerance);
			Print(statusOFS, "ColPerm                = ", ColPerm );
			Print(statusOFS, "muMaxIter              = ", muMaxIter);
			Print(statusOFS, "mpisize                = ", mpisize );
			Print(statusOFS, "npPerPole              = ", npPerPole );
			Print(statusOFS, "isFreeEnergyMatrix     = ", isFreeEnergyDensityMatrix );
			Print(statusOFS, "isEnergyMatrix         = ", isEnergyDensityMatrix ); 


			GetTime( timeSolveSta );

			pexsi.Solve( 
					numPole,
					temperature,
					numElectronExact,
					gap,
					deltaE,
					mu0,
					muMin,
					muMax,
					HMat,
					SMat,
					muMaxIter,
					numElectronTolerance,
					ColPerm,
					npSymbFact,
					isFreeEnergyDensityMatrix,
					isEnergyDensityMatrix,
					isDerivativeTMatrix,
					muList,
					numElectronList,
					numElectronDrvList,
					isConverged	);

			GetTime( timeSolveEnd );

			PrintBlock( statusOFS, "Solve finished." );
			if( isConverged ){
				statusOFS << "PEXSI has converged with " << muList.size() << 
					" iterations" << std::endl;
			}
			else {
				statusOFS << "PEXSI did not converge with " << muList.size() << 
					" iterations" << std::endl;
			}
			Print( statusOFS, "mu                   = ", 
					*muList.rbegin() );
			Print( statusOFS, "Total time for PEXSI = ", 
					timeSolveEnd - timeSolveSta );

		}

		statusOFS.close();
	}
	catch( std::exception& e )
	{
		statusOFS << std::endl << " ERROR!!! Proc " << mpirank << " caught exception with message: "
			<< e.what() << std::endl;
		statusOFS.close();
		statusOFS << std::endl << " ERROR!!! Proc " << mpirank << " caught exception with message: "
			<< e.what() << std::endl;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}
	
	MPI_Finalize();

	return 0;
}
