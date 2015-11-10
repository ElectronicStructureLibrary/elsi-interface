!Copyright (c) 2015, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!> 
!! This program illustratrates how to read and solve an eigenvalue problem 
!! created by the ELSI interface
!!

program read_and_solve

  use iso_c_binding
  use ELSI

  implicit none
  include 'mpif.h'

  ! This is the ELSI test suite
  ! First we will test the writing and reading of a matrix to a file

  integer :: blocksize = 16
  real*8  :: n_electrons

  integer :: myid, i_eigenvector
  real*8  :: e_tot


  !  Pharse command line argumnents, if given
   INTEGER*4 :: iargc
   character*16 arg1

   if (iargc() == 1) then
      call getarg(1, arg1)
      read(arg1, *) n_electrons
   endif


  ! First: the parallel treatment
  call elsi_initialize_mpi()
  
  ! Second set some ELSI specifications
  call elsi_set_method(PEXSI)
  call elsi_set_mode(REAL_VALUES)
  
  ! Third define the problem
  call elsi_initialize_problem_from_file("elsi_eigenvalue_problem.hdf5",&
        blocksize, blocksize)
  
  ! Initialize problem distribution
  call elsi_initialize_blacs()

  ! Read eigenvalue problem
  call elsi_allocate_matrices()
  call elsi_read_ev_problem("elsi_eigenvalue_problem.hdf5")

  ! Solve the eigenvalue problem
  call elsi_solve_ev_problem(n_electrons)
  
  call elsi_get_total_energy(e_tot)

  call elsi_get_myid(myid)
  if (myid == 0) then 
     write (*,'(A,E19.12)') "total energy : ", e_tot
  end if

  ! elsi shutdown
  call elsi_finalize()

end program
