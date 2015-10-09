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
!! This program illustrates how to read and write an eigenvalue problem 
!! created by the ELSI Interface
!! 

program read_and_write

  use iso_c_binding
  use ELSI

  implicit none

  ! First the parallel treatment
  call elsi_initialize_mpi()
  
  ! Second set some ELSI specifications
  call elsi_set_method(ELPA)
  call elsi_set_mode(REAL_VALUES)
  
  ! Third define the problem
  call elsi_initialize_problem_from_file("elsi_eigenvalue_problem.hdf5")
  
  ! Forth distribute the problem
  call elsi_initialize_blacs()

  ! Read eigenvalue problem
  call elsi_allocate_matrices()
  call elsi_read_ev_problem("elsi_eigenvalue_problem.hdf5")

  ! Write eigenvalue problem to another file
  call elsi_write_ev_problem("elsi_eigenvalue_problem_out.hdf5")

  ! elsi shutdown
  call elsi_finalize()

end program
