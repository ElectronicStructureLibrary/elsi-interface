!   This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), fomerly known as
!      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
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

! ELPA2 -- 2-stage solver for ELPA
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

!> \brief Fortran module which provides the routines to use the 2-stage ELPA solver
module ELPA2

  use elpa_utilities
  use elpa1, only : elpa_print_times, time_evp_back, time_evp_fwd, time_evp_solve
  use elpa2_utilities

  implicit none

  private

  public :: solve_evp_real_2stage_double
  public :: solve_evp_complex_2stage_double
  public :: elpa_solve_evp_real_2stage_double
  public :: elpa_solve_evp_complex_2stage_double

  interface solve_evp_real_2stage
    module procedure solve_evp_real_2stage_double
  end interface

  interface elpa_solve_evp_real_2stage_double
    module procedure solve_evp_real_2stage_double
  end interface

  interface solve_evp_complex_2stage
    module procedure solve_evp_complex_2stage_double
  end interface

  interface elpa_solve_evp_complex_2stage_double
    module procedure solve_evp_complex_2stage_double
  end interface

  contains
!-------------------------------------------------------------------------------
!>  \brief solve_evp_real_2stage_double: Fortran function to solve the double-precision real eigenvalue problem with a 2 stage approach
!>
!>  Parameters
!>
!>  \param na                                   Order of matrix a
!>
!>  \param nev                                  Number of eigenvalues needed
!>
!>  \param a(lda,matrixCols)                    Distributed matrix for which eigenvalues are to be computed.
!>                                              Distribution is like in Scalapack.
!>                                              The full matrix must be set (not only one half like in scalapack).
!>                                              Destroyed on exit (upper and lower half).
!>
!>  \param lda                                  Leading dimension of a
!>
!>  \param ev(na)                               On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)                    On output: Eigenvectors of a
!>                                              Distribution is like in Scalapack.
!>                                              Must be always dimensioned to the full size (corresponding to (na,na))
!>                                              even if only a part of the eigenvalues is needed.
!>
!>  \param ldq                                  Leading dimension of q
!>
!>  \param nblk                                 blocksize of cyclic distribution, must be the same in both directions!
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------

function solve_evp_real_2stage_double(na,nev,a,lda,ev,q,ldq,nblk,&
  mpi_comm_rows,mpi_comm_cols,mpi_comm_all) result(success)

  use elpa1_compute
  use elpa2_compute
  use elpa_mpi
  use, intrinsic :: iso_c_binding

  implicit none

  integer(kind=c_int) :: THIS_REAL_ELPA_KERNEL

  integer(kind=c_int), intent(in) :: na, nev, lda, ldq, mpi_comm_rows, &
                                     mpi_comm_cols, mpi_comm_all
  integer(kind=c_int), intent(in) :: nblk
  real(kind=c_double), intent(inout) :: ev(na)

  real(kind=c_double), intent(inout) :: a(lda,*), q(ldq,*)

  real(kind=c_double), allocatable :: hh_trans_real(:,:)

  integer(kind=c_int) :: my_pe, n_pes, my_prow, my_pcol, np_rows, np_cols, mpierr
  integer(kind=c_int) :: nbw, num_blocks
  real(kind=c_double), allocatable :: tmat(:,:,:), e(:)
  real(kind=c_double) :: ttt0, ttt1
  logical :: success
  integer(kind=c_int) :: istat
  character(200) :: errorMessage

  call MPI_Comm_rank(mpi_comm_all,my_pe,mpierr)
  call MPI_Comm_size(mpi_comm_all,n_pes,mpierr)
  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

  success = .true.

  THIS_REAL_ELPA_KERNEL = DEFAULT_REAL_ELPA_KERNEL

! Bandwidth must be a multiple of nblk
! Set to a value >= 32
  nbw = (63/nblk+1)*nblk
  num_blocks = (na-1)/nbw+1

  allocate(tmat(nbw,nbw,num_blocks),stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_real_2stage: error when allocating tmat "//errorMessage
    stop
  endif

! Reduction full -> band
  ttt0 = MPI_Wtime()

  call bandred_real_double(na,a,lda,nblk,nbw,mpi_comm_rows,mpi_comm_cols,tmat,success)

  if(.not.(success)) return
  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time full => band           :",ttt1-ttt0," s"
  endif

! Reduction band -> tridiagonal
  allocate(e(na),stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_real_2stage: error when allocating e "//errorMessage
    stop
  endif

  ttt0 = MPI_Wtime()

  call tridiag_band_real_double(na,nbw,nblk,a,lda,ev,e,hh_trans_real,&
       mpi_comm_rows,mpi_comm_cols,mpi_comm_all)

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time band => tridiagonal    :",ttt1-ttt0," s"
  endif

  call MPI_Bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
  call MPI_Bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

! Solve tridiagonal system
  ttt0 = MPI_Wtime()

  call solve_tridi_double(na,nev,ev,e,q,ldq,nblk,mpi_comm_rows,mpi_comm_cols,success)

  if(.not.(success)) return
  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"
  endif

  deallocate(e,stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_real_2stage: error when deallocating e "//errorMessage
    stop
  endif

! Backtransform stage 1
  ttt0 = MPI_Wtime()

  call trans_ev_tridi_to_band_real_double(na,nev,nblk,nbw,q,ldq,hh_trans_real,&
       mpi_comm_rows,mpi_comm_cols,success,THIS_REAL_ELPA_KERNEL)

  if(.not.(success)) return
  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time ev tridiagonal => band :",ttt1-ttt0," s"
  endif

! We can now deallocate the stored householder vectors
  deallocate(hh_trans_real,stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_real_2stage: error when deallocating hh_trans_real "//errorMessage
    stop
  endif

! Backtransform stage 2
  ttt0 = MPI_Wtime()

  call trans_ev_band_to_full_real_double(na,nev,nblk,nbw,a,lda,tmat,q,ldq,&
       mpi_comm_rows,mpi_comm_cols)

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time ev band => full        :",ttt1-ttt0," s"
  endif

  deallocate(tmat,stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_real_2stage: error when deallocating tmat"//errorMessage
    stop
  endif

end function solve_evp_real_2stage_double

!>  \brief solve_evp_complex_2stage_double: Fortran function to solve the double-precision complex eigenvalue problem with a 2 stage approach
!>
!>  Parameters
!>
!>  \param na                                   Order of matrix a
!>
!>  \param nev                                  Number of eigenvalues needed
!>
!>  \param a(lda,matrixCols)                    Distributed matrix for which eigenvalues are to be computed.
!>                                              Distribution is like in Scalapack.
!>                                              The full matrix must be set (not only one half like in scalapack).
!>                                              Destroyed on exit (upper and lower half).
!>
!>  \param lda                                  Leading dimension of a
!>
!>  \param ev(na)                               On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)                    On output: Eigenvectors of a
!>                                              Distribution is like in Scalapack.
!>                                              Must be always dimensioned to the full size (corresponding to (na,na))
!>                                              even if only a part of the eigenvalues is needed.
!>
!>  \param ldq                                  Leading dimension of q
!>
!>  \param nblk                                 blocksize of cyclic distribution, must be the same in both directions!
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------

function solve_evp_complex_2stage_double(na,nev,a,lda,ev,q,ldq,nblk,&
  mpi_comm_rows,mpi_comm_cols,mpi_comm_all) result(success)

  use elpa1_compute
  use elpa2_compute
  use elpa_mpi
  use, intrinsic :: iso_c_binding

  implicit none

  integer(kind=c_int) :: THIS_COMPLEX_ELPA_KERNEL
  integer(kind=c_int), intent(in) :: na, nev, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  real(kind=c_double), intent(inout) :: ev(na)

  complex(kind=c_double), intent(inout) :: a(lda,*), q(ldq,*)

  complex(kind=c_double), allocatable :: hh_trans_complex(:,:)

  integer(kind=c_int) :: my_prow, my_pcol, np_rows, np_cols, mpierr, my_pe, n_pes
  integer(kind=c_int) :: l_cols, l_rows, l_cols_nev, nbw, num_blocks
  complex(kind=c_double), allocatable :: tmat(:,:,:)
  real(kind=c_double), allocatable :: q_real(:,:), e(:)
  real(kind=c_double) :: ttt0, ttt1

  logical :: success
  integer(kind=c_int) :: istat
  character(200) :: errorMessage

  call MPI_Comm_rank(mpi_comm_all,my_pe,mpierr)
  call MPI_Comm_size(mpi_comm_all,n_pes,mpierr)
  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

  success = .true.

  THIS_COMPLEX_ELPA_KERNEL = DEFAULT_COMPLEX_ELPA_KERNEL

! Choose bandwidth, must be a multiple of nblk, set to a value >= 32
  nbw = (31/nblk+1)*nblk
  num_blocks = (na-1)/nbw + 1

  allocate(tmat(nbw,nbw,num_blocks),stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_complex_2stage: error when allocating tmat"//errorMessage
    stop
  endif

! Reduction full -> band
  ttt0 = MPI_Wtime()

  call bandred_complex_double(na,a,lda,nblk,nbw,mpi_comm_rows,mpi_comm_cols,tmat,success)

  if(.not.(success)) then
    return
  endif

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time full => band           :",ttt1-ttt0," s"
  endif

! Reduction band -> tridiagonal
  allocate(e(na),stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_complex_2stage: error when allocating e"//errorMessage
    stop
  endif

  ttt0 = MPI_Wtime()

  call tridiag_band_complex_double(na,nbw,nblk,a,lda,ev,e,hh_trans_complex,&
       mpi_comm_rows,mpi_comm_cols,mpi_comm_all)

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time band => tridiagonal    :",ttt1-ttt0," s"
  endif

  call MPI_Bcast(ev,na,mpi_real8,0,mpi_comm_all,mpierr)
  call MPI_Bcast(e,na,mpi_real8,0,mpi_comm_all,mpierr)

  l_rows = local_index(na,my_prow,np_rows,nblk,-1) ! Local rows of a and q
  l_cols = local_index(na,my_pcol,np_cols,nblk,-1) ! Local columns of q
  l_cols_nev = local_index(nev,my_pcol,np_cols,nblk,-1) ! Local columns corresponding to nev

  allocate(q_real(l_rows,l_cols),stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_complex_2stage: error when allocating q_real"//errorMessage
    stop
  endif

! Solve tridiagonal system
  ttt0 = MPI_Wtime()

  call solve_tridi_double(na,nev,ev,e,q_real,ubound(q_real,dim=1),nblk,&
       mpi_comm_rows,mpi_comm_cols,success)

  if(.not.(success)) return

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"
  endif

  q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

  deallocate(e,q_real,stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_complex_2stage: error when deallocating e, q_real"//errorMessage
    stop
  endif

! Backtransform stage 1
  ttt0 = MPI_Wtime()

  call trans_ev_tridi_to_band_complex_double(na,nev,nblk,nbw,q,ldq,hh_trans_complex,&
       mpi_comm_rows,mpi_comm_cols,success,THIS_COMPLEX_ELPA_KERNEL)

  if(.not.(success)) return
  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time ev tridiagonal => band :",ttt1-ttt0," s"
  endif

! We can now deallocate the stored householder vectors
  deallocate(hh_trans_complex,stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_complex_2stage: error when deallocating hh_trans_complex"//errorMessage
    stop
  endif

! Backtransform stage 2
  ttt0 = MPI_Wtime()

  call trans_ev_band_to_full_complex_double(na,nev,nblk,nbw,a,lda,tmat,q,ldq,&
       mpi_comm_rows,mpi_comm_cols)

  ttt1 = MPI_Wtime()
  if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) then
    write(use_unit,"(A,F10.3,A)") "  | Time ev band => full        :",ttt1-ttt0," s"
  endif

  deallocate(tmat,stat=istat,errmsg=errorMessage)
  if(istat .ne. 0) then
    write(error_unit,*) "solve_evp_complex_2stage: error when deallocating tmat "//errorMessage
    stop
  endif

end function solve_evp_complex_2stage_double

end module ELPA2
