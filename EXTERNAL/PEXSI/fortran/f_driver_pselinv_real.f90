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
!> @file f_driver_pselinv_real.f90
!> @brief FORTRAN version of the driver for using PSelInv for real symmetric
!> matrices.
!> @date 2014-04-01
program f_driver_pselinv_real
use f_ppexsi_interface
use iso_c_binding
implicit none
include 'mpif.h'

integer(c_int) :: nrows, nnz, nnzLocal, numColLocal
integer(c_int), allocatable, dimension(:) ::  colptrLocal, rowindLocal
real(c_double), allocatable, dimension(:) ::  &
  HnzvalLocal, SnzvalLocal, AnzvalLocal, AinvnzvalLocal
integer(c_int):: nprow, npcol, npSymbFact, outputFileIndex
integer :: mpirank, mpisize, ierr
double precision:: timeSta, timeEnd
character*32 :: Hfile
integer(c_int):: info
integer(c_intptr_t) :: plan
type(f_ppexsi_options) :: options

integer:: i, j 
integer:: numColLocalFirst, firstCol
integer:: irow, jcol


call mpi_init( ierr )
call mpi_comm_rank( MPI_COMM_WORLD, mpirank, ierr )
call mpi_comm_size( MPI_COMM_WORLD, mpisize, ierr )

Hfile            = "lap2dr.matrix"

nprow = 1
npcol = 1

call f_read_distsparsematrix_formatted_head( &
  trim(Hfile)//char(0),&
  nrows,&
  nnz,&
  nnzLocal,&
  numColLocal,&
  MPI_COMM_WORLD )

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
allocate( AnzvalLocal( nnzLocal ) )
allocate( AinvnzvalLocal( nnzLocal ) )

call f_read_distsparsematrix_formatted (&
  trim(Hfile)//char(0),&
  nrows,&
  nnz,&
  nnzLocal,&
  numColLocal,&
  colptrLocal,&
  rowindLocal,&
  HnzvalLocal,&
  MPI_COMM_WORLD )

do i = 1, nnzLocal
  AnzvalLocal(i) = HnzvalLocal(i)
enddo

if( mpirank .eq. 0 ) then
  write(*,*) "Finish reading the matrix" 
endif

! Each processor outputs information
outputFileIndex = mpirank

plan = f_ppexsi_plan_initialize(&
  MPI_COMM_WORLD,&
  nprow,&
  npcol,&
  outputFileIndex,&
  info )

if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif

call f_ppexsi_set_default_options(&
  options )

if( mpirank .eq. 0 ) then
  write(*,*) "Finish initializing the PPEXSI plan" 
endif

call f_ppexsi_load_real_symmetric_hs_matrix(&
      plan,&       
      options,&
      nrows,&
      nnz,&
      nnzLocal,&
      numColLocal,&
      colptrLocal,& 
      rowindLocal,&
      HnzvalLocal,&
      1,&
      SnzvalLocal,&
      info ) 

if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif


if( mpirank == 0 ) then
  write(*,*)  "Finish setting up the matrix."
endif

call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&
  plan,&
  options,&
  info)

if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif


if( mpirank == 0 ) then
  write(*,*)  "Finish symbolic factorization."
endif

call f_ppexsi_selinv_real_symmetric_matrix(&
  plan,&
  options,&
  AnzvalLocal,&
  AinvnzvalLocal,&
  info)

if( info .ne. 0 ) then
	call mpi_finalize( ierr )
	call exit(info)
endif


if( mpirank == 0 ) then
  write(*,*)  "Finish selected inversion."
endif


! The first processor output the diagonal elements in natural order
if( mpirank == 0 ) then
  numColLocalFirst = nrows / mpisize;
  firstCol         = mpirank * numColLocalFirst;
  do j = 1, numColLocal
    jcol = firstCol + j
    do i = colptrLocal(j), colptrLocal(j+1)
      irow = rowindLocal(i)
      if( irow == jcol ) then
        write(*,"(A,2I5,A,F15.10)") &
          "Ainv[", irow, irow, "]= ", AinvnzvalLocal(i)
      endif
    enddo
  enddo
endif

call f_ppexsi_plan_finalize( plan, info )

call mpi_finalize( ierr )

deallocate( colptrLocal )
deallocate( rowindLocal )
deallocate( HnzvalLocal )
deallocate( AnzvalLocal )
deallocate( AinvnzvalLocal )


end program f_driver_pselinv_real

