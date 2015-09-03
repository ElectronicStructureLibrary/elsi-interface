! ELPA2 -- 2-stage solver for ELPA
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

module ELPA2

! Version 1.1.2, 2011-02-21

  USE ELPA1
  use elpa_pdgeqrf
#ifdef GPU_VERSION
  use cuda_routines
  use cuda_c_kernel
  use iso_c_binding
#endif

  implicit none

  PRIVATE ! By default, all routines contained are private

  ! The following routines are public:

  public :: solve_evp_real_2stage
  public :: solve_evp_complex_2stage

  public :: bandred_real
  public :: tridiag_band_real
  public :: trans_ev_tridi_to_band_real
  public :: trans_ev_band_to_full_real

  public :: bandred_complex
  public :: tridiag_band_complex
  public :: trans_ev_tridi_to_band_complex
  public :: trans_ev_band_to_full_complex

  private :: validate_input

  integer, public :: which_qr_decomposition = 1     ! defines, which QR-decomposition algorithm will be used
                                                    ! 0 for unblocked
                                                    ! 1 for blocked (maxrank: nblk)
!-------------------------------------------------------------------------------

  ! The following array contains the Householder vectors of the
  ! transformation band -> tridiagonal.
  ! It is allocated and set in tridiag_band_real and used in
  ! trans_ev_tridi_to_band_real.
  ! It must be deallocated by the user after trans_ev_tridi_to_band_real!

  real*8, allocatable :: hh_trans_real(:,:)
  complex*16, allocatable :: hh_trans_complex(:,:)

!-------------------------------------------------------------------------------

include 'mpif.h'

contains

subroutine solve_evp_real_2stage(na, nev, a, lda, ev, q, ldq, sdq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_real_2stage: Solves the real eigenvalue problem with a 2 stage approach
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, sdq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
   real*8, intent(inout) :: a(lda,*), ev(na), q(ldq,sdq)

   integer my_pe, n_pes, my_prow, my_pcol, np_rows, np_cols, mpierr
   integer nbw, num_blocks
   real*8, allocatable :: tmat(:,:,:), e(:)
   real*8 ttt0, ttt1, ttts

   ! First, check if the input is valid
   if (.not. validate_input(na, nev, nblk)) return

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32

!   nbw = (31/nblk+1)*nblk
   nbw = nblk

   num_blocks = (na-1)/nbw + 1

   allocate(tmat(nbw,nbw,num_blocks))

   ! Reduction full -> band

   ttt0 = MPI_Wtime()
   ttts = ttt0
   call bandred_real(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time bandred_real               :',ttt1-ttt0

   ! Reduction band -> tridiagonal

   allocate(e(na))

   ttt0 = MPI_Wtime()
   call tridiag_band_real(na, nbw, nblk, a, lda, ev, e, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time tridiag_band_real          :',ttt1-ttt0

   call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)
 
   ttt1 = MPI_Wtime()
   time_evp_fwd = ttt1-ttts

   ! Solve tridiagonal system

   ttt0 = MPI_Wtime()
   call solve_tridi(na, nev, ev, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time solve_tridi                :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0
   ttts = ttt1

   deallocate(e)

   ! Backtransform stage 1

   ttt0 = MPI_Wtime()
   call trans_ev_tridi_to_band_real(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_tridi_to_band_real:',ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   deallocate(hh_trans_real)

   ! Backtransform stage 2

   ttt0 = MPI_Wtime()
   call trans_ev_band_to_full_real(na, nev, nblk, nbw, a, lda, tmat, q, ldq, sdq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_band_to_full_real :',ttt1-ttt0
   time_evp_back = ttt1-ttts

   deallocate(tmat)
!print*,'EV1:',ev(1),  
1 format(a,f10.3) 
end subroutine solve_evp_real_2stage

#ifdef GPU_VERSION
#warning "Compiling ELPA GPU Version."
!! Use this one for the GPU
!!real datatype
#include "bandred_real_gpu_v2.f90" 
#include "tridiag_band_real_gpu_v1.f90"
#include "ev_tridi_band_gpu_v3.f90"
#include "ev_band_full_gpu_v1.f90"

!complex datatype
#include "bandred_complex_gpu_v2.f90"
#include "tridiag_band_complex.f90"
#include "ev_tridi_band_complex_gpu_v3.f90"
#include "ev_band_full_complex_gpu_v1.f90"

#else 
#warning "Compiling ELPA CPU Version"
#warning "Set the -DGPU_VERSION flag for the CPU verison"
!real datatype
#include "bandred_real.f90"
#include "tridiag_band_real.f90"
#include "ev_tridi_band.f90"
#include "ev_band_full.f90"

!complex datatype
#include "bandred_complex.f90"
#include "tridiag_band_complex.f90"
#include "ev_tridi_complex.f90"
#include "ev_band_complex.f90"
#endif

!-------------------------------------------------------------------------------

!---BELOW ARE THE ROUTINES FOR THE COMPLEX CASE

!-------------------------------------------------------------------------------




subroutine solve_evp_complex_2stage(na, nev, a, lda, ev, q, ldq, sdq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)

!-------------------------------------------------------------------------------
!  solve_evp_complex_2stage: Solves the complex eigenvalue problem with a 2 stage approach
!
!  Parameters
!
!  na          Order of matrix a
!
!  nev         Number of eigenvalues needed
!
!  a(lda,*)    Distributed matrix for which eigenvalues are to be computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in scalapack).
!              Destroyed on exit (upper and lower half).
!
!  lda         Leading dimension of a
!
!  ev(na)      On output: eigenvalues of a, every processor gets the complete set
!
!  q(ldq,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  ldq         Leading dimension of q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm_all
!              MPI-Communicator for the total processor set
!
!-------------------------------------------------------------------------------

   implicit none

   integer, intent(in) :: na, nev, lda, ldq, sdq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
   complex*16, intent(inout) :: a(lda,*), q(ldq,sdq)
   real*8, intent(inout) :: ev(na)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows, l_cols_nev, nbw, num_blocks
   complex*16, allocatable :: tmat(:,:,:)
   real*8, allocatable :: q_real(:,:), e(:)
   real*8 ttt0, ttt1, ttts

   ! First, check if the input is valid

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32

   nbw = (31/nblk+1)*nblk

   num_blocks = (na-1)/nbw + 1

   allocate(tmat(nbw,nbw,num_blocks))

   ! Reduction full -> band

   ttt0 = MPI_Wtime()
   ttts = ttt0
   call bandred_complex(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time bandred_complex               :',ttt1-ttt0

   ! Reduction band -> tridiagonal

   allocate(e(na))

   ttt0 = MPI_Wtime()
   call tridiag_band_complex(na, nbw, nblk, a, lda, ev, e, mpi_comm_rows, mpi_comm_cols, mpi_comm_all)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time tridiag_band_complex          :',ttt1-ttt0

   call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ttt1 = MPI_Wtime()
   time_evp_fwd = ttt1-ttts

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q
   l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

   allocate(q_real(l_rows,l_cols))

   ! Solve tridiagonal system

   ttt0 = MPI_Wtime()
   call solve_tridi(na, nev, ev, e, q_real, ubound(q_real,1), nblk, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times)  &
      print 1,'Time solve_tridi                   :',ttt1-ttt0
   time_evp_solve = ttt1-ttt0
   ttts = ttt1

   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   deallocate(e, q_real)

   ! Backtransform stage 1

   ttt0 = MPI_Wtime()
   call trans_ev_tridi_to_band_complex(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_tridi_to_band_complex:',ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   deallocate(hh_trans_complex)

   ! Backtransform stage 2

   ttt0 = MPI_Wtime()
   call trans_ev_band_to_full_complex(na, nev, nblk, nbw, a, lda, tmat, q, ldq, sdq, mpi_comm_rows, mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print 1,'Time trans_ev_band_to_full_complex :',ttt1-ttt0
   time_evp_back = ttt1-ttts

   deallocate(tmat)

1  format(a,f10.3)

end subroutine solve_evp_complex_2stage


!-------------------------------------------------------------------------------

! Check if the given input is valid
logical function validate_input(na, nev, nblk)

    integer, intent(in) :: na, nev, nblk

    if (na .le. 0) then
        print *, 'Fatal error: The matrix dimension must be at least 1'
        validate_input = .false.
        return
    endif

    if ((nev .le. 0) .or. (nev .gt. na)) then
        print *, 'Fatal error: The number of EVs must be between 1 and the matrix dimension'
        validate_input = .false.
        return
    endif

    if ((nblk .lt. 1) .or. (nblk .gt. 128)) then
        print *, 'Fatal error: The SCALAPACK block size must be between 1 and 128'
        validate_input = .false.
        return
    endif

    validate_input = .true.
    return
end function

end module ELPA2
