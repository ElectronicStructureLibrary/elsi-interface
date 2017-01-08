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

module ELPA1
  use, intrinsic :: iso_c_binding, only : c_double
  use elpa_utilities
  use elpa1_auxiliary

  implicit none

  private

  public :: get_elpa_row_col_comms           !< old, deprecated interface: Sets MPI row/col communicators DO NOT USE
  public :: get_elpa_communicators           !< Sets MPI row/col communicators

  public :: solve_evp_real                   !< old, deprecated interface: Driver routine for real double-precision eigenvalue problem DO NOT USE
  public :: solve_evp_real_1stage            !< Driver routine for real double-precision eigenvalue problem
  public :: solve_evp_real_1stage_double     !< Driver routine for real double-precision eigenvalue problem
  public :: solve_evp_complex                !< old, deprecated interface:  Driver routine for complex double-precision eigenvalue problem DO NOT USE
  public :: solve_evp_complex_1stage         !< Driver routine for complex double-precision eigenvalue problem
  public :: solve_evp_complex_1stage_double  !< Driver routine for complex double-precision eigenvalue problem

  ! imported from elpa1_auxilliary
  public :: elpa_mult_at_b_real_double       !< Multiply double-precision real matrices A**T * B
  public :: mult_at_b_real                   !< old, deprecated interface to multiply double-precision real matrices A**T * B  DO NOT USE

  public :: elpa_mult_ah_b_complex_double    !< Multiply double-precision complex matrices A**H * B
  public :: mult_ah_b_complex                !< old, deprecated interface to multiply double-preicion complex matrices A**H * B  DO NOT USE

  public :: elpa_invert_trm_real_double      !< Invert double-precision real triangular matrix
  public :: invert_trm_real                  !< old, deprecated interface to invert double-precision real triangular matrix  DO NOT USE

  public :: elpa_invert_trm_complex_double   !< Invert double-precision complex triangular matrix
  public :: invert_trm_complex               !< old, deprecated interface to invert double-precision complex triangular matrix  DO NOT USE

  public :: elpa_cholesky_real_double        !< Cholesky factorization of a double-precision real matrix
  public :: cholesky_real                    !< old, deprecated interface to do Cholesky factorization of a double-precision real matrix  DO NOT USE

  public :: elpa_cholesky_complex_double     !< Cholesky factorization of a double-precision complex matrix
  public :: cholesky_complex                 !< old, deprecated interface to do Cholesky factorization of a double-precision complex matrix  DO NOT USE

  public :: elpa_solve_tridi_double          !< Solve a double-precision tridiagonal eigensystem with divide and conquer method

  logical, public :: elpa_print_times = .false. !< Set elpa_print_times to .true. for explicit timing outputs

  interface get_elpa_row_col_comms
    module procedure get_elpa_communicators
  end interface

  interface solve_evp_real
    module procedure solve_evp_real_1stage_double
  end interface

  interface solve_evp_real_1stage
    module procedure solve_evp_real_1stage_double
  end interface

  interface solve_evp_complex
    module procedure solve_evp_complex_1stage_double
  end interface

  interface solve_evp_complex_1stage
    module procedure solve_evp_complex_1stage_double
  end interface

contains

!-------------------------------------------------------------------------------
!> \param  mpi_comm_global   Global communicator for the calculations (in)
!>
!> \param  my_prow           Row coordinate of the calling process in the process grid (in)
!>
!> \param  my_pcol           Column coordinate of the calling process in the process grid (in)
!>
!> \param  mpi_comm_rows     Communicator for communicating within rows of processes (out)
!>
!> \param  mpi_comm_cols     Communicator for communicating within columns of processes (out)
!>
!> \result mpierr            integer error value of mpi_comm_split function
!-------------------------------------------------------------------------------

function get_elpa_communicators(mpi_comm_global,my_prow,my_pcol,mpi_comm_rows,&
                                mpi_comm_cols) result(mpierr)

   use precision
   use elpa_mpi

   implicit none

   integer(kind=ik), intent(in)  :: mpi_comm_global,my_prow,my_pcol
   integer(kind=ik), intent(out) :: mpi_comm_rows,mpi_comm_cols

   integer(kind=ik)              :: mpierr

   ! mpi_comm_rows is used for communicating WITHIN rows, i.e. all processes
   ! having the same column coordinate share one mpi_comm_rows
   ! So the "color" for splitting is my_pcol and the "key" is my row coordinate
   ! Analogous for mpi_comm_cols
   call mpi_comm_split(mpi_comm_global,my_pcol,my_prow,mpi_comm_rows,mpierr)
   call mpi_comm_split(mpi_comm_global,my_prow,my_pcol,mpi_comm_cols,mpierr)

end function get_elpa_communicators

!-------------------------------------------------------------------------------
!> \param  na               Order of matrix a
!>
!> \param  nev              Number of eigenvalues needed
!>                          The smallest nev eigenvalues/eigenvectors are calculated
!>
!> \param  a(lda,*)         Distributed matrix for which eigenvalues are to be computed
!>                          Distribution is like in ScaLAPACK
!>                          The full matrix must be set
!>                          Destroyed on exit (upper and lower half)
!>
!>  \param lda              Leading dimension of a
!>
!>  \param ev(na)           Output: eigenvalues of a
!>                          All processors get the complete set
!>
!>  \param q(ldq,*)         Output: Eigenvectors of a
!>                          Distribution is like in ScaLAPACK
!>                          Must be always dimensioned to the full size (na,na)
!>                          even if only a part of the eigenvalues is needed
!>
!>  \param ldq              Leading dimension of q
!>
!>  \param nblk             Blocksize of cyclic distribution
!>                          Must be the same in both directions
!>
!>  \param matrixCols       Distributed number of matrix columns
!>
!>  \param mpi_comm_rows    MPI-Communicator for rows
!>  \param mpi_comm_cols    MPI-Communicator for columns
!>
!>  \result                 success
!-------------------------------------------------------------------------------

function solve_evp_real_1stage_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
                                      mpi_comm_rows,mpi_comm_cols) result(success)
   use precision
   use iso_c_binding
   use elpa_mpi
   use elpa1_compute

   implicit none

   integer(kind=ik), intent(in) :: na,nev,lda,ldq,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols
!   real(kind=REAL_DATATYPE)     :: a(lda,matrixCols),ev(na),q(ldq,matrixCols)
   real(kind=rk8)               :: a(lda,*),ev(na),q(ldq,*)

   integer(kind=ik)             :: my_prow,my_pcol,mpierr
   real(kind=rk8), allocatable  :: e(:),tau(:)
   real(kind=c_double)          :: ttt0,ttt1
   logical                      :: success
   logical, save                :: firstCall = .true.
   logical                      :: wantDebug

   success   = .true.
   firstCall = .false.
   wantDebug = .false.

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)

   allocate(e(na))
   allocate(tau(na))

   ttt0 = MPI_Wtime()
   call tridiag_real_double(na,a,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols,ev,e,tau)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print *,"  Time full ==> tridiagonal    :",ttt1-ttt0

   ttt0 = MPI_Wtime()
   call solve_tridi_double(na,nev,ev,e,q,ldq,nblk,matrixCols,mpi_comm_rows,&
                           mpi_comm_cols,wantDebug,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print *,"  Time solve tridiagonal       :",ttt1-ttt0

   ttt0 = MPI_Wtime()
   call trans_ev_real_double(na,nev,a,lda,tau,q,ldq,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print *,"  Time ev tridiagonal ==> full :",ttt1-ttt0

   deallocate(e)
   deallocate(tau)

end function solve_evp_real_1stage_double

!-------------------------------------------------------------------------------
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
!>  \param matrixCols           distributed number of matrix columns
!>
!>  \param mpi_comm_rows        MPI-Communicator for rows
!>  \param mpi_comm_cols        MPI-Communicator for columns
!>
!>  \result                     success
!-------------------------------------------------------------------------------

function solve_evp_complex_1stage_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
                                         mpi_comm_rows,mpi_comm_cols) result(success)

   use precision
   use iso_c_binding
   use elpa_mpi
   use elpa1_compute

   implicit none

   integer(kind=ik), intent(in) :: na,nev,lda,ldq,nblk,matrixCols
   integer(kind=ik), intent(in) :: mpi_comm_rows,mpi_comm_cols
   complex(kind=ck8)            :: a(lda,*),q(ldq,*)
   real(kind=rk8)               :: ev(na)

   integer(kind=ik)               :: my_prow,my_pcol,np_rows,np_cols,mpierr
   integer(kind=ik)               :: l_rows,l_cols,l_cols_nev
   real(kind=rk8), allocatable    :: q_real(:,:),e(:)
   complex(kind=ck8), allocatable :: tau(:)
   real(kind=c_double)            :: ttt0,ttt1
   logical                        :: success
   logical, save                  :: firstCall = .true.
   logical                        :: wantDebug

   success   = .true.
   wantDebug = .false.
   firstCall = .false.

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   l_rows = local_index(na,my_prow,np_rows,nblk,-1) ! Local rows of a and q
   l_cols = local_index(na,my_pcol,np_cols,nblk,-1) ! Local columns of q
   l_cols_nev = local_index(nev,my_pcol,np_cols,nblk,-1) ! Local columns corresponding to nev

   allocate(e(na))
   allocate(tau(na))
   allocate(q_real(l_rows,l_cols))

   ttt0 = MPI_Wtime()
   call tridiag_complex_double(na,a,lda,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols,ev,e,tau)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print *,"  Time full ==> tridiagonal    :",ttt1-ttt0

   ttt0 = MPI_Wtime()
   call solve_tridi_double(na,nev,ev,e,q_real,l_rows,nblk,matrixCols,mpi_comm_rows,&
                           mpi_comm_cols,wantDebug,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print *,"  Time solve tridiagonal       :",ttt1-ttt0

   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   ttt0 = MPI_Wtime()
   call trans_ev_complex_double(na,nev,a,lda,tau,q,ldq,nblk,matrixCols,mpi_comm_rows,mpi_comm_cols)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      print *,"  Time ev tridiagonal ==> full :",ttt1-ttt0

   deallocate(q_real)
   deallocate(e)
   deallocate(tau)

end function solve_evp_complex_1stage_double

end module ELPA1
