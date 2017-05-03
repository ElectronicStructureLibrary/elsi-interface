! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This program tests ELSI density matrix solver.
!!
program test_dm_real_sparse

   use ELSI_PRECISION, only : dp
   use ELSI

   implicit none

   include "mpif.h"

   character(128) :: arg1

   integer :: n_proc,myid
   integer :: mpi_comm_global,mpierr
   integer :: matrix_size
   integer :: n_states
   integer :: solver
   integer :: global_nnz,local_nnz,local_col
   integer :: n_pexsi_poles
   integer, allocatable :: row_idx(:),col_ptr(:)

   real(kind=dp) :: n_electrons
   real(kind=dp) :: e_test
   real(kind=dp) :: t1,t2
   real(kind=dp), allocatable :: H(:),S(:),D(:)

   ! NOTE: The energy computed by PEXSI, e_test, is only accurate if using
   !       enough poles. In this test program, the number of poles is set to
   !       half of the number of MPI tasks. Therefore, converged result is
   !       expected with at least 80 MPI tasks, i.e. 40 PEXSI poles.

   ! Initialize MPI
   call MPI_Init(mpierr)
   mpi_comm_global = MPI_COMM_WORLD
   call MPI_Comm_size(mpi_comm_global,n_proc,mpierr)
   call MPI_Comm_rank(mpi_comm_global,myid,mpierr)

   if(mod(n_proc,4) /= 0) then
      if(myid == 0) then
         write(*,'("  ################################################")')
         write(*,'("  ##  This test program must be launched with   ##")')
         write(*,'("  ##  the number of MPI tasks being a multiple  ##")')
         write(*,'("  ##  of 4, i.e. 4, 8, 12, ...                  ##")')
         write(*,'("  ################################################")')
         call MPI_Abort(mpi_comm_global,0,mpierr)
         stop
      endif
   endif

   ! Read command line arguments
   if(COMMAND_ARGUMENT_COUNT() == 1) then
      call GET_COMMAND_ARGUMENT(1,arg1)
      read(arg1,*) solver
      if(solver /= 3) then
         if(myid == 0) then
            write(*,'("  ################################################")')
            write(*,'("  ##  Wrong choice of solver!!                  ##")')
            write(*,'("  ##  Please choose:                            ##")')
            write(*,'("  ##  PEXSI = 3                                 ##")')
            write(*,'("  ################################################")')
            call MPI_Abort(mpi_comm_global,0,mpierr)
            stop
         endif
      endif
   else
      solver = 3
   endif

   if(myid == 0) then
      write(*,'("  ################################")')
      write(*,'("  ##     ELSI TEST PROGRAMS     ##")')
      write(*,'("  ################################")')
      write(*,*)
      if(solver == 3) then
         write(*,'("  This test program performs the following computational steps:")')
         write(*,*)
         write(*,'("  1) Generates Hamiltonian and overlap matrices;")')
         write(*,'("  2) Computes the density matrix with pole expansion and selected")')
         write(*,'("     inversion method.")')
         write(*,*)
         write(*,'("  Now start testing  elsi_dm_real_sparse + PEXSI")')
      endif
      write(*,*)
   endif

   ! Generate test matrices by hand
   global_nnz = 44
   local_nnz = 22
   local_col = 5

   allocate(H(local_nnz))
   allocate(S(local_nnz))
   allocate(D(local_nnz))
   allocate(row_idx(local_nnz))
   allocate(col_ptr(local_col+1))

   ! Only first pole needs the data
   if(myid == 0) then
      row_idx = (/1,2,4,6,7,9,1,2,4,6,7,9,3,8,1,2,4,6,7,9,5,10/)
      col_ptr = (/1,7,13,15,21,23/)
      H = (/  -0.388670241778891_dp,-4.311480451396187e-2_dp,   -0.172349864644966_dp,&
              -0.412610292189610_dp,   -0.143358918306738_dp,    0.170711531858969_dp,&
           -4.311480451396187e-2_dp,    0.409475020776789_dp,-3.277717121147698e-2_dp,&
              -0.143358918306738_dp,-6.487842673409908e-3_dp, 4.623415594693785e-2_dp,&
               0.819652099421313_dp,    0.140028009400052_dp,   -0.172349864644966_dp,&
           -3.277717121147698e-2_dp,    0.622390391755544_dp,   -0.170711531858969_dp,&
           -4.623415594693791e-2_dp,   -0.257530556898739_dp,    0.819652099421314_dp,&
               0.140028009400052_dp/)
      S = (/   0.999983934950546_dp, 4.016953945433239e-5_dp,-4.480947479392433e-6_dp,&
               0.773518727907150_dp,   -0.228105494172661_dp,   -0.444961515210702_dp,&
            4.016953945433239e-5_dp,    0.999915572079161_dp, 1.859665487610805e-5_dp,&
              -0.228105494172661_dp,    0.477261561338787_dp,   -0.413062386431313_dp,&
               0.999999144962076_dp,    0.582454761626444_dp,-4.480947479392433e-6_dp,&
            1.859665487610805e-5_dp,     1.00000367037239_dp,    0.444961515210702_dp,&
               0.413062386431313_dp, 1.260959980250462e-2_dp,    0.999999144962076_dp,&
               0.582454761626444_dp/)
   elseif(myid == 1) then
      row_idx = (/1,2,4,6,7,9,1,2,4,6,7,9,3,8,1,2,4,6,7,9,5,10/)
      col_ptr = (/1,7,13,15,21,23/)
      H = (/  -0.412610292189610_dp,   -0.143358918306738_dp,   -0.170711531858969_dp,&
              -0.388670241778891_dp,-4.311480451396186e-2_dp,    0.172349864644966_dp,&
              -0.143358918306738_dp,-6.487842673409908e-3_dp,-4.623415594693791e-2_dp,&
           -4.311480451396186e-2_dp,    0.409475020776790_dp, 3.277717121147693e-2_dp,&
               0.140028009400052_dp,    0.819652099421313_dp,    0.170711531858969_dp,&
            4.623415594693785e-2_dp,   -0.257530556898739_dp,    0.172349864644966_dp,&
            3.277717121147693e-2_dp,    0.622390391755544_dp,    0.140028009400052_dp,&
               0.819652099421314_dp/)
      S = (/   0.773518727907150_dp,   -0.228105494172661_dp,    0.444961515210702_dp,&
               0.999983934950546_dp, 4.016953945417670e-5_dp, 4.480947479398607e-6_dp,&
              -0.228105494172661_dp,    0.477261561338787_dp,    0.413062386431313_dp,&
            4.016953945417670e-5_dp,    0.999915572079161_dp,-1.859665487617809e-5_dp,&
               0.582454761626444_dp,    0.999999144962076_dp,   -0.444961515210702_dp,&
              -0.413062386431313_dp, 1.260959980250462e-2_dp, 4.480947479398607e-6_dp,&
           -1.859665487617809e-5_dp,     1.00000367037239_dp,    0.582454761626444_dp,&
               0.999999144962076_dp/)
   else
      row_idx = 0
      col_ptr = 0
      H = 0.0_dp
      S = 0.0_dp
   endif

   if(myid == 0) then
      write(*,'("  Finished test matrices generation")')
   endif

   ! Initialize ELSI
   matrix_size = 10
   n_electrons = 2.0_dp
   n_states = 1
   n_pexsi_poles = n_proc/2

   call elsi_init(solver,1,1,matrix_size,n_electrons,n_states)
   call elsi_set_mpi(mpi_comm_global)
   call elsi_set_csc(global_nnz,local_nnz,local_col,row_idx,col_ptr)

   ! Only 2 PEXSI poles for quick test
   call elsi_customize_pexsi(temperature=4e-4_dp,&
                             n_poles=n_pexsi_poles,&
                             max_iteration=1,&
                             mu_min=-2.0_dp,&
                             mu_max=2.0_dp,&
                             mu_inertia_expansion=1.0_dp,&
                             n_electron_accuracy=1e-5_dp)

   ! Customize ELSI
   call elsi_customize(print_detail=.true.)

   t1 = MPI_Wtime()

   ! Solve problem
   call elsi_dm_real_sparse(H,S,D,e_test)

   t2 = MPI_Wtime()

   if(myid == 0) then
      write(*,'("  Finished test program")')
      write(*,'("  | Total computation time :",F10.3,"s")') t2-t1
      write(*,*)
   endif

   ! Finalize ELSI
   call elsi_finalize()

   deallocate(H)
   deallocate(S)
   deallocate(D)
   deallocate(row_idx)
   deallocate(col_ptr)

   call MPI_Finalize(mpierr)

end program
