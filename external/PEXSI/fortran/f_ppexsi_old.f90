!	 Copyright (c) 2012 The Regents of the University of California,
!	 through Lawrence Berkeley National Laboratory.  
!
!  Author: Lin Lin
!	 
!  This file is part of PEXSI. All rights reserved.
!
!	 Redistribution and use in source and binary forms, with or without
!	 modification, are permitted provided that the following conditions are met:
!
!	 (1) Redistributions of source code must retain the above copyright notice, this
!	 list of conditions and the following disclaimer.
!	 (2) Redistributions in binary form must reproduce the above copyright notice,
!	 this list of conditions and the following disclaimer in the documentation
!	 and/or other materials provided with the distribution.
!	 (3) Neither the name of the University of California, Lawrence Berkeley
!	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
!	 be used to endorse or promote products derived from this software without
!	 specific prior written permission.
!
!	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
!	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!	 You are under no obligation whatsoever to provide any bug fixes, patches, or
!	 upgrades to the features, functionality or performance of the source code
!	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
!	 available either publicly, or directly to Lawrence Berkeley National
!	 Laboratory, without imposing a separate written license agreement for such
!	 Enhancements, then you hereby grant the following license: a non-exclusive,
!	 royalty-free perpetual license to install, use, modify, prepare derivative
!	 works, incorporate into other computer software, distribute, and sublicense
!	 such enhancements or derivative works thereof, in binary and source code form.
!
!> @file f_ppexsi.f90
!> @brief Test the new FORTRAN interface for PPEXSI.
!> @date 2013-04-10
program f_ppexsi
implicit none
include 'mpif.h'

integer :: nrows, nnz, nnzLocal, numColLocal
integer, allocatable, dimension(:) ::  colptrLocal, rowindLocal
integer, allocatable, dimension(:) ::  inertiaListInt
double precision, allocatable, dimension(:) :: &
	HnzvalLocal, SnzvalLocal, DMnzvalLocal, EDMnzvalLocal, &
	FDMnzvalLocal, muList, numElectronList, numElectronDrvList,&
	shiftList, inertiaList, localDOSnzvalLocal
double precision:: Energy, eta
integer :: numPole
double precision :: temperature, numElectronExact, numElectron,&
	gap, deltaE
double precision ::   muMin0, muMax0, muInertia, muMinInertia, muMaxInertia,&
	muLowerEdge, muUpperEdge, muPEXSI, muMinPEXSI, muMaxPEXSI

integer:: inertiaMaxIter, inertiaIter, muMaxIter, muIter
integer:: ordering
integer:: isInertiaCount
double precision :: PEXSINumElectronTolerance, &
	inertiaNumElectronTolerance
integer:: npPerPole, nprow, npcol, npSymbFact
integer :: mpirank, mpisize, ierr
double precision:: timeSta, timeEnd
character*32 :: Hfile, Sfile
integer:: isFormatted
integer:: isSIdentity
! Communicator for reading the matrix, with size npPerPole
integer:: readComm
integer:: isProcRead
integer:: info
integer:: i


namelist/InputVars/ &
	temperature         ,&
	numElectronExact    ,&
	numPole             ,&
	gap                 ,&
	deltaE              ,&
	muMin0              ,&
	muMax0              ,&
	inertiaMaxIter      ,&
	muMaxIter           ,&
	inertiaNumElectronTolerance ,&
	PEXSINumElectronTolerance   ,&
	npPerPole           ,&
	npSymbFact          ,&
	Hfile               ,&
	Sfile               ,&
  isFormatted         ,&
	isSIdentity         ,&
	ordering            


call mpi_init( ierr )
call mpi_comm_rank( MPI_COMM_WORLD, mpirank, ierr )
call mpi_comm_size( MPI_COMM_WORLD, mpisize, ierr )


! Below is the data for a DNA matrix
! Temperature should be in the same unit as the H matrix. Here it is Rydberg.
temperature      = 0.0019d0

numElectronExact = 2442.0d0
numPole          = 40
gap              = 0.0d0
! deltaE is in theory the spectrum width, but in practice can be much smaller
! than | E_max - mu |.  It is found that deltaE that is slightly bigger
! than  | E_min - mu | is usually good enough.
deltaE           = 30.0d0
! Initial searching interval for the chemical potential for inertia count.
muMin0           = -1.0d0
muMax0           =  0.0d0
! Maximum number of iterations for computing the inertia
inertiaMaxIter   = 1
! Maximum number of iterations for PEXSI iteration
muMaxIter        = 3
! Stop inertia count if Ne(muMax) - Ne(muMin) < inertiaNumElectronTolerance
inertiaNumElectronTolerance = 100
! Stop mu-iteration if numElectronTolerance is < numElectronTolerance.
PEXSINumElectronTolerance = 1d-2
! Number of processors used for each pole. At the moment use mpisize.
! Later can be changed to 
npPerPole        = 4
! Number of processors used for paralle symbolic factorization. This is only
! relevant if PARMETIS/PT-SCOTCH is used.
npSymbFact       = 1

Hfile            = "H.matrix"
! Empty Sfile means the overlap matrix is identity
Sfile            = "S.matrix"
! Whether S is an identity matrix
isSIdentity      = 1
! Whether the H and S files are formatted
isFormatted      = 1

! Ordering 
!   0   : PARMETIS
!   1   : METIS_AT_PLUS_A
!   2   : MMD_AT_PLUS_A
ordering         = 2

! Read the actual input parameters from standard input
if( mpirank  == 0 ) then
	read(unit=*, nml=InputVars)
endif

! Broadcast the input parameters

call MPI_BCAST( temperature,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( numElectronExact, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( numPole, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( gap, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( deltaE, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( muMin0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( muMax0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( inertiaMaxIter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( muMaxIter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 
call MPI_BCAST( inertiaNumElectronTolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( PEXSINumElectronTolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( npPerPole, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 
call MPI_BCAST( npSymbFact, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 
call MPI_BCAST( Hfile, len(Hfile), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( Sfile, len(Sfile), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
call MPI_BCAST( isFormatted, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 
call MPI_BCAST( isSIdentity, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 
call MPI_BCAST( ordering, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr ) 


! Read and compute the size/local size of the arrays 
! The conversion of the string to fit the C format is important.

! Split the processors to read matrix
if( mpirank < npPerPole ) then
	isProcRead = 1
else
	isProcRead = 0
endif

call mpi_comm_split( MPI_COMM_WORLD, isProcRead, mpirank, readComm, ierr )


allocate( muList( muMaxIter ) )
allocate( numElectronList( muMaxIter ) )
allocate( numElectronDrvList( muMaxIter ) )

allocate( shiftList( numPole ) )
allocate( inertiaList( numPole ) )

if( isProcRead == 1 ) then
	write(*,*) "Proc ", mpirank, " is reading file..."
  if( isFormatted == 1 ) then
		call f_read_distsparsematrix_formatted_head( &
			trim(Hfile)//char(0),&
			nrows,&
			nnz,&
			nnzLocal,&
			numColLocal,&
			readComm )
	else
		call f_read_distsparsematrix_head( &
			trim(Hfile)//char(0),&
			nrows,&
			nnz,&
			nnzLocal,&
			numColLocal,&
			readComm )
	endif

	if( mpirank .eq. 0 ) then
		write(*,*) "Matrix size (local data on proc 0):" 
		write(*,*) "size = ", nrows
		write(*,*) "nnz  = ", nnz
		write(*,*) "nnzLocal = ", nnzLocal
		write(*,*) "numColLocal = ", numColLocal
	endif

	! Allocate memory
	allocate( colptrLocal( numColLocal + 1 ) )
	allocate( rowindLocal( nnzLocal ) )
	allocate( HnzvalLocal( nnzLocal ) )
	allocate( SnzvalLocal( nnzLocal ) ) 
	allocate( DMnzvalLocal( nnzLocal ) ) 
	allocate( EDMnzvalLocal( nnzLocal ) ) 
	allocate( FDMnzvalLocal( nnzLocal ) ) 
	allocate( localDOSnzvalLocal( nnzLocal ) )

	timeSta = mpi_wtime()

	if( isFormatted == 1 ) then
		call f_read_distsparsematrix_formatted (&
			trim(Hfile)//char(0),&
			nrows,&
			nnz,&
			nnzLocal,&
			numColLocal,&
			colptrLocal,&
			rowindLocal,&
			HnzvalLocal,&
			readComm )
	else
		call f_para_read_distsparsematrix (&
			trim(Hfile)//char(0),&
			nrows,&
			nnz,&
			nnzLocal,&
			numColLocal,&
			colptrLocal,&
			rowindLocal,&
			HnzvalLocal,&
			readComm )
	endif

	if( isSIdentity == 0 ) then
		if( isFormatted == 1 ) then
			call f_read_distsparsematrix_formatted (&
				trim(Sfile)//char(0),&
				nrows,&
				nnz,&
				nnzLocal,&
				numColLocal,&
				colptrLocal,&
				rowindLocal,&
				SnzvalLocal,&
				readComm )	
		else
			call f_para_read_distsparsematrix (&
				trim(Sfile)//char(0),&
				nrows,&
				nnz,&
				nnzLocal,&
				numColLocal,&
				colptrLocal,&
				rowindLocal,&
				SnzvalLocal,&
				readComm )	
		endif
	endif

	timeEnd = mpi_wtime()

endif

call mpi_barrier( MPI_COMM_WORLD, ierr )

if( mpirank == 0 ) then
	write(*,*) "Time for reading H/S matrices is ", &
		timeEnd - timeSta, " [s]"
endif

! Step 1. Estimate the range of chemical potential

call f_ppexsi_inertiacount_interface(&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
	rowindLocal,&
	HnzvalLocal,&
	isSIdentity,&
	SnzvalLocal,&
	temperature,&
	numElectronExact,&
	muMin0,&
	muMax0,&
	numPole,&
	inertiaMaxIter,&
	inertiaNumElectronTolerance,&
	ordering,&
	npPerPole,&
	npSymbFact,&
	MPI_COMM_WORLD,&
	muMinInertia,&
	muMaxInertia,&
	muLowerEdge,&
	muUpperEdge,&
	inertiaIter,&
	shiftList,&
	inertiaList,&
	info)


if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif

muInertia = (muLowerEdge + muUpperEdge)/2.0;

if( mpirank == 0 ) then
	write(*,*) "The computed finite temperature inertia = "
	do i = 1, numPole
		write(*,*) "Shift = ", shiftList(i), "inertia = ", inertiaList(i)
	enddo
endif


! Step 2. Solve KSDFT using PEXSISolve

call f_ppexsi_solve_interface(&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
  rowindLocal,&
	HnzvalLocal,&
	isSIdentity,&
	SnzvalLocal,&
	temperature,&
	numElectronExact,&
	muInertia,&
	muMinInertia,&
	muMaxInertia,&
	gap,&
	deltaE,&
	numPole,&
	muMaxIter,&
	PEXSINumElectronTolerance,&
	ordering,&
	npPerPole,&
	npSymbFact,&
	MPI_COMM_WORLD,&
	DMnzvalLocal,&
	EDMnzvalLocal,&
	FDMnzvalLocal,&
	muPEXSI,&
	numElectron,&
	muMinPEXSI,&
	muMaxPEXSI,&
	muIter,&
	muList,&
	numElectronList,&
	numElectronDrvList,&
	info)

if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif

if( mpirank == 0 ) then
	write(*,*) "PEXSI Solve finished."
endif

! Step 3. Post processing

! Compute the "raw inertia" after the calculation

if( mpirank == 0 ) then
	write(*,*) 
	write(*,*)  "Compute the raw inertia after the calculation."
endif


allocate( inertiaListInt( numPole ) ) 

call f_ppexsi_raw_inertiacount_interface(&
	nrows,&
	nnz,&
	nnzLocal,&
	numColLocal,&
	colptrLocal,&
	rowindLocal,&
	HnzvalLocal,&
	isSIdentity,&
	SnzvalLocal,&
	muMinInertia,&
	muMaxInertia,&
	numPole,&
	ordering,&
	npPerPole,&
	npSymbFact,&
	MPI_COMM_WORLD,&
	shiftList,&
	inertiaListInt,&
	info)

if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif


if( mpirank == 0 ) then
	write(*,*) "The computed zero temperature RAW inertia = "
	do i = 1, numPole
		write(*,*) "Shift = ", shiftList(i), "inertia = ", inertiaListInt(i)
	enddo
endif


! Compute the local density of states
! Only the first pole group participates in the computation of the selected
! inversion for a single shift.
if( isProcRead == 1 ) then
	Energy = 1.0;
	eta    = 0.001;
	call f_ppexsi_localdos_interface(&
		nrows,&
		nnz,&
		nnzLocal,&
		numColLocal,&
		colptrLocal,&
		rowindLocal,&
		HnzvalLocal,&
		isSIdentity,&
		SnzvalLocal,&
		Energy,&
		eta,&
		ordering,&
		npSymbFact,&
		readComm,&
		localDOSnzvalLocal,&
		info)

	if( info .ne. 0 ) then
		call mpi_finalize( ierr )
		call exit(info)
	endif

end if

call mpi_barrier( MPI_COMM_WORLD, ierr )


! Final output.

if( mpirank == 0 ) then
  write(*, *)
	write(*, *) "After inertia count,"
	write(*, *) "muMinInertia  = ", muMinInertia
	write(*, *) "muMaxInertia  = ", muMaxInertia
	write(*, *) "muLowerEdge   = ", muLowerEdge
	write(*, *) "muUpperEdge   = ", muUpperEdge
	write(*, *) "inertiaIter   = ", inertiaIter
	! write(*, *) "shift ,           inertia count"
	! do i = 1, numPole 
		! write(*,*) shiftList(i), inertiaList(i)
	! enddo
	write(*, *) 
	write(*, *) "After PEXSI,"
	write(*, *) "muPEXSI       = ", muPEXSI
	write(*, *) "numElectron   = ", numElectron
	write(*, *) "muMinPEXSI    = ", muMinPEXSI
	write(*, *) "muMaxPEXSI    = ", muMaxPEXSI
	write(*, *) "muIter        = ", muIter
endif

deallocate( muList )
deallocate( numElectronList )
deallocate( numElectronDrvList )

deallocate( shiftList )
deallocate( inertiaList )
deallocate( inertiaListInt )

if( isProcRead == 1 ) then
	deallocate( colptrLocal )
	deallocate( rowindLocal )
	deallocate( HnzvalLocal )
	deallocate( SnzvalLocal )
	deallocate( DMnzvalLocal )
	deallocate( EDMnzvalLocal )
	deallocate( FDMnzvalLocal )
	deallocate( localDOSnzvalLocal );
endif


call mpi_comm_free( readComm, ierr )
call mpi_finalize( ierr )

end program f_ppexsi

