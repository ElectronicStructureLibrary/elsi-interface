!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!      Informatik,
!    - Technische Universität München, Lehrstuhl für Informatik mit
!      Schwerpunkt Wissenschaftliches Rechnen ,
!    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!      and
!    - IBM Deutschland GmbH
!
!    This particular source code file contains additions, changes and
!    enhancements authored by Intel Corporation which is not part of
!    the ELPA consortium.
!
!    More information can be found here:
!    http://elpa.mpcdf.mpg.de/
!
!    ELPA is free software: you can redistribute it and/or modify
!    it under the terms of the version 3 of the license of the
!    GNU Lesser General Public License as published by the Free
!    Software Foundation.
!
!    ELPA is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
!
! ELPA1 -- Faster replacements for ScaLAPACK symmetric eigenvalue routines
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

!> \brief Fortran module which provides the routines to use the one-stage ELPA solver
module ELPA1

  use, intrinsic :: iso_c_binding
  use elpa_utilities
  use elpa1_auxiliary

  implicit none

  private

  public :: get_elpa_row_col_comms
  public :: get_elpa_communicators
  public :: elpa_get_communicators

  public :: solve_evp_real
  public :: elpa_solve_evp_real_1stage_double

  public :: solve_evp_real_1stage
  public :: solve_evp_real_1stage_double

  public :: solve_evp_complex
  public :: elpa_solve_evp_complex_1stage_double
  public :: solve_evp_complex_1stage
  public :: solve_evp_complex_1stage_double

! imported from elpa1_auxilliary

  public :: elpa_mult_at_b_real_double
  public :: mult_at_b_real

  public :: elpa_mult_ah_b_complex_double
  public :: mult_ah_b_complex

  public :: elpa_invert_trm_real_double
  public :: invert_trm_real

  public :: elpa_invert_trm_complex_double
  public :: invert_trm_complex

  public :: elpa_cholesky_real_double
  public :: cholesky_real

  public :: elpa_cholesky_complex_double
  public :: cholesky_complex

  public :: elpa_solve_tridi_double

! Timing results, set by every call to solve_evp_xxx

  real(kind=c_double), public :: time_evp_fwd
  real(kind=c_double), public :: time_evp_solve
  real(kind=c_double), public :: time_evp_back

  logical, public :: elpa_print_times = .false.

  interface get_elpa_row_col_comms
    module procedure elpa_get_communicators
  end interface

  interface get_elpa_communicators
    module procedure elpa_get_communicators
  end interface

  interface solve_evp_real
    module procedure elpa_solve_evp_real_1stage_double
  end interface

  interface solve_evp_real_1stage
    module procedure elpa_solve_evp_real_1stage_double
  end interface

  interface solve_evp_real_1stage_double
    module procedure elpa_solve_evp_real_1stage_double
  end interface

  interface solve_evp_complex
    module procedure elpa_solve_evp_complex_1stage_double
  end interface

  interface solve_evp_complex_1stage
    module procedure elpa_solve_evp_complex_1stage_double
  end interface

  interface solve_evp_complex_1stage_double
    module procedure elpa_solve_evp_complex_1stage_double
  end interface

contains

!-------------------------------------------------------------------------------

! All ELPA routines need MPI communicators for communicating within
! rows or columns of processes, these are set here.
! mpi_comm_rows/mpi_comm_cols can be free'd with MPI_Comm_free if not used any more.
!
!  Parameters
!
!> \param  mpi_comm_global   Global communicator for the calculations (in)
!>
!> \param  my_prow           Row coordinate of the calling process in the process grid (in)
!>
!> \param  my_pcol           Column coordinate of the calling process in the process grid (in)
!>
!> \param  mpi_comm_rows     Communicator for communicating within rows of processes (out)
!>
!> \param  mpi_comm_cols     Communicator for communicating within columns of processes (out)
!> \result mpierr            integer error value of mpi_comm_split function

function elpa_get_communicators(mpi_comm_global, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols) result(mpierr)

  use elpa_mpi
  use, intrinsic :: iso_c_binding

  implicit none

  integer(kind=c_int), intent(in) :: mpi_comm_global, my_prow, my_pcol
  integer(kind=c_int), intent(out) :: mpi_comm_rows, mpi_comm_cols

  integer(kind=c_int) :: mpierr

! mpi_comm_rows is used for communicating WITHIN rows, i.e. all processes
! having the same column coordinate share one mpi_comm_rows.
! So the "color" for splitting is my_pcol and the "key" is my row coordinate.
! Analogous for mpi_comm_cols

  call MPI_Comm_split(mpi_comm_global,my_pcol,my_prow,mpi_comm_rows,mpierr)
  call MPI_Comm_split(mpi_comm_global,my_prow,my_pcol,mpi_comm_cols,mpierr)

end function elpa_get_communicators

!> \brief elpa_solve_evp_real_1stage_double: Fortran function to solve the real double-precision eigenvalue problem with 1-stage solver
!>
!  Parameters
!
!> \param  na                   Order of matrix a
!>
!> \param  nev                  Number of eigenvalues needed.
!>                              The smallest nev eigenvalues/eigenvectors are calculated.
!>
!> \param  a(lda,matrixCols)    Distributed matrix for which eigenvalues are to be computed.
!>                              Distribution is like in Scalapack.
!>                              The full matrix must be set (not only one half like in scalapack).
!>                              Destroyed on exit (upper and lower half).
!>
!>  \param lda                  Leading dimension of a
!>
!>  \param ev(na)               On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)    On output: Eigenvectors of a
!>                              Distribution is like in Scalapack.
!>                              Must be always dimensioned to the full size (corresponding to (na,na))
!>                              even if only a part of the eigenvalues is needed.
!>
!>  \param ldq                  Leading dimension of q
!>
!>  \param nblk                 blocksize of cyclic distribution, must be the same in both directions!
!>
!>  \param mpi_comm_rows        MPI-Communicator for rows
!>  \param mpi_comm_cols        MPI-Communicator for columns
!>  \param mpi_comm_all         global MPI communicator
!>
!>  \result                     success

function elpa_solve_evp_real_1stage_double(na, nev, a, lda, ev, q, ldq, nblk, &
  mpi_comm_rows, mpi_comm_cols, mpi_comm_all) result(success)

  use precision
  use, intrinsic :: iso_c_binding
  use elpa_mpi
  use elpa1_compute

  implicit none

  integer(kind=c_int), intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, &
                                     mpi_comm_cols, mpi_comm_all
  real(kind=c_double), intent(out) :: ev(na)

  real(kind=c_double), intent(inout) :: a(lda,*)
  real(kind=c_double), intent(out) :: q(ldq,*)

  logical :: success

  integer(kind=c_int) :: my_pe, n_pes, my_prow, my_pcol, mpierr
  real(kind=c_double), allocatable :: e(:), tau(:)
  real(kind=c_double) :: ttt0, ttt1 ! MPI_WTIME always needs double

  call MPI_Comm_rank(mpi_comm_all,my_pe,mpierr)
  call MPI_Comm_size(mpi_comm_all,n_pes,mpierr)

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)

  success = .true.

  allocate(e(na), tau(na))

  ttt0 = MPI_Wtime()

  call tridiag_real_double(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, ev, e, tau)

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time full => tridiagonal    :",ttt1-ttt0," s"
  endif

  ttt0 = MPI_Wtime()

  call solve_tridi_double(na, nev, ev, e, q, ldq, nblk, mpi_comm_rows, &
       mpi_comm_cols, success)

  if (.not.(success)) return

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"
  endif

  ttt0 = MPI_Wtime()

  call trans_ev_real_double(na, nev, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

  ttt1 = MPI_Wtime()

  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time ev tridiagonal => full :",ttt1-ttt0," s"
  endif

  deallocate(e, tau)

end function elpa_solve_evp_real_1stage_double

!> \brief elpa_solve_evp_complex_1stage_double: Fortran function to solve the complex double-precision eigenvalue problem with 1-stage solver
!>
!  Parameters
!
!> \param  na                   Order of matrix a
!>
!> \param  nev                  Number of eigenvalues needed.
!>                              The smallest nev eigenvalues/eigenvectors are calculated.
!>
!> \param  a(lda,matrixCols)    Distributed matrix for which eigenvalues are to be computed.
!>                              Distribution is like in Scalapack.
!>                              The full matrix must be set (not only one half like in scalapack).
!>                              Destroyed on exit (upper and lower half).
!>
!>  \param lda                  Leading dimension of a
!>
!>  \param ev(na)               On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)    On output: Eigenvectors of a
!>                              Distribution is like in Scalapack.
!>                              Must be always dimensioned to the full size (corresponding to (na,na))
!>                              even if only a part of the eigenvalues is needed.
!>
!>  \param ldq                  Leading dimension of q
!>
!>  \param nblk                 blocksize of cyclic distribution, must be the same in both directions!
!>
!>  \param mpi_comm_rows        MPI-Communicator for rows
!>  \param mpi_comm_cols        MPI-Communicator for columns
!>  \param mpi_comm_all         global MPI Communicator
!>
!>  \result                     success

function elpa_solve_evp_complex_1stage_double(na, nev, a, lda, ev, q, ldq, nblk, &
  mpi_comm_rows, mpi_comm_cols, mpi_comm_all) result(success)

  use precision
  use, intrinsic :: iso_c_binding
  use elpa_mpi
  use elpa1_compute

  implicit none

  integer(kind=c_int), intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all

  complex(kind=c_double), intent(inout) :: a(lda,*)
  complex(kind=c_double), intent(out) :: q(ldq,*)

  real(kind=c_double), intent(out) :: ev(na)

  logical :: success

  integer(kind=c_int) :: my_pe, n_pes, my_prow, my_pcol, np_rows, np_cols, mpierr
  integer(kind=c_int) :: l_rows, l_cols, l_cols_nev
  real(kind=c_double), allocatable :: q_real(:,:), e(:)
  complex(kind=c_double), allocatable :: tau(:)
  real(kind=c_double) :: ttt0, ttt1 ! MPI_WTIME always needs double

  call MPI_Comm_rank(mpi_comm_all,my_pe,mpierr)
  call MPI_Comm_size(mpi_comm_all,n_pes,mpierr)

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

  success = .true.

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

  l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

  allocate(e(na), tau(na))
  allocate(q_real(l_rows,l_cols))

  ttt0 = MPI_Wtime()

  call tridiag_complex_double(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols, ev, e, tau)

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time full => tridiagonal    :",ttt1-ttt0," s"
  endif

  ttt0 = MPI_Wtime()

  call solve_tridi_double(na, nev, ev, e, q_real, l_rows, nblk, mpi_comm_rows, &
       mpi_comm_cols, success)

  if (.not.(success)) return

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"
  endif

  ttt0 = MPI_Wtime()
  q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

  call trans_ev_complex_double(na, nev, a, lda, tau, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time ev tridiagonal => full :",ttt1-ttt0," s"
  endif

  deallocate(q_real)
  deallocate(e, tau)

end function elpa_solve_evp_complex_1stage_double

end module ELPA1
