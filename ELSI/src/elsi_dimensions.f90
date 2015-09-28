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
!! This module contains variables accessible in ELSI and related modules
!!

module ELSI_DIMENSIONS

  use iso_c_binding

  implicit none

  !> Global Matrix Dimensions
  integer :: n_g_rank       !< Rank of eigenvalue problem 
  
  !> Global Matrix Blocks
  integer :: n_b_rows       !< Number of rows 
  integer :: n_b_cols       !< Number of columns 

  !> Processor grid
  integer :: n_p_rows       !< Number of processor rows 
  integer :: n_p_cols       !< Number of processor columns
  
  !> Local Matrix Dimensions
  integer :: n_l_rows        !< Number of rows for this process
  integer :: n_l_cols        !< Number of columns for this process

  !> MPI variables
  integer :: myid            !< local process id
  integer :: n_procs         !< number of mpi processes
  integer :: mpi_comm_global !< global mpi communicator
  integer :: mpierr          !< mpi error handler
  logical :: external_mpi    !< Has somebody else initialzed mpi

  !> BLACS variables
  integer :: blacs_ctxt     !< local blacs context
  integer :: sc_desc(9)     !< local blacs context
  integer :: mpi_comm_row   !< row    communicatior
  integer :: mpi_comm_col   !< column communicatior
  integer :: my_p_row       !< process row    position
  integer :: my_p_col       !< process column position
  integer :: blacs_info     !< Info code for blacs related calls 
  logical :: external_blacs !< Has somebody else initialized blacs? 

  !> HDF5 variables
  integer :: h5err          !< HDF5 error code

  !> Overlap
  logical :: overlap_is_unity = .True. !< Is the overlap unity
  integer :: n_eigenvectors            !< Number of eigenvectors to be calculated

  !> ELPA variables
  real*8, parameter :: elpa_step_switch = 1.d0 !< This parameter sets the threshold when to switch from ELPA2 to ELPA1 
  
  contains

  !EMPTY

end module ELSI_DIMENSIONS
