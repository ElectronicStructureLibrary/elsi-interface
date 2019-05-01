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

! Version 1.1.2, 2011-02-21

  use elpa_utilities
  use elpa1, only : elpa_print_times, time_evp_back, time_evp_fwd, time_evp_solve
  use elpa2_utilities

  implicit none

  PRIVATE ! By default, all routines contained are private

! The following routines are public:

  public :: solve_evp_real_2stage_double               !< old, deprecated interface: Driver routine for real double-precision eigenvalue problem. will be deleted at some point
  public :: solve_evp_complex_2stage_double            !< old, deprecated interface: Driver routine for complex double-precision eigenvalue problem. will be deleted at some point
  public :: elpa_solve_evp_real_2stage_double          !< Driver routine for real double-precision 2-stage eigenvalue problem
  public :: elpa_solve_evp_complex_2stage_double       !< Driver routine for complex double-precision 2-stage eigenvalue problem

!-------------------------------------------------------------------------------
!>  \brief solve_evp_real_2stage: Old, deprecated interface for elpa_solve_evp_real_2stage_double
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
!>  \param matrixCols                           local columns of matrix a and q
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \param THIS_REAL_ELPA_KERNEL_API (optional) specify used ELPA2 kernel via API
!>
!>  \param useQR (optional)                     use QR decomposition
!>  \param useGPU (optional)                    decide whether to use GPUs or not

!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------
  interface solve_evp_real_2stage
    module procedure solve_evp_real_2stage_double
  end interface

!-------------------------------------------------------------------------------
!>  \brief elpa_solve_evp_real_2stage_double: Fortran function to solve the real double-precision eigenvalue problem with a 2 stage approach. This is called by "elpa_solve_evp_real_double"
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
!>  \param matrixCols                           local columns of matrix a and q
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \param THIS_REAL_ELPA_KERNEL_API (optional) specify used ELPA2 kernel via API
!>
!>  \param useQR (optional)                     use QR decomposition
!>  \param useGPU (optional)                    decide whether to use GPUs or not
!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------
  interface elpa_solve_evp_real_2stage_double
    module procedure solve_evp_real_2stage_double
  end interface

!-------------------------------------------------------------------------------
!>  \brief solve_evp_complex_2stage: Old, deprecated interface for elpa_solve_evp_complex_2stage_double
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
!>  \param matrixCols                           local columns of matrix a and q
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \param THIS_REAL_ELPA_KERNEL_API (optional) specify used ELPA2 kernel via API
!>
!>  \param useGPU (optional)                    decide whether to use GPUs or not
!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------
  interface solve_evp_complex_2stage
    module procedure solve_evp_complex_2stage_double
  end interface

!-------------------------------------------------------------------------------
!>  \brief elpa_solve_evp_complex_2stage_double: Fortran function to solve the complex double-precision eigenvalue problem with a 2 stage approach. This is called by "elpa_solve_evp_complex_double"
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
!>  \param matrixCols                           local columns of matrix a and q
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \param THIS_REAL_ELPA_KERNEL_API (optional) specify used ELPA2 kernel via API
!>
!>  \param useGPU (optional)                    decide whether to use GPUs or not
!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------
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
!>  \param matrixCols                           local columns of matrix a and q
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \param THIS_REAL_ELPA_KERNEL_API (optional) specify used ELPA2 kernel via API
!>
!>  \param useQR (optional)                     use QR decomposition
!>  \param useGPU (optional)                    decide whether to use GPUs or not
!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------

function solve_evp_real_2stage_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
            mpi_comm_rows,mpi_comm_cols,mpi_comm_all,THIS_REAL_ELPA_KERNEL_API,&
            useQR,useGPU) result(success)

   use elpa1_compute
   use elpa2_compute
   use elpa_mpi
   use, intrinsic :: iso_c_binding

   implicit none

   logical, intent(in), optional             :: useQR, useGPU
   logical                                   :: useQRActual
   integer(kind=c_int), intent(in), optional :: THIS_REAL_ELPA_KERNEL_API
   integer(kind=c_int)                       :: THIS_REAL_ELPA_KERNEL

   integer(kind=c_int), intent(in)           :: na, nev, lda, ldq, matrixCols, mpi_comm_rows, &
                                                mpi_comm_cols, mpi_comm_all
   integer(kind=c_int), intent(in)           :: nblk
   real(kind=c_double), intent(inout)        :: ev(na)

   real(kind=c_double), intent(inout)        :: a(lda,*), q(ldq,*)

   real(kind=c_double), allocatable          :: hh_trans_real(:,:)

   integer(kind=c_int)                       :: my_pe, n_pes, my_prow, my_pcol, np_rows, np_cols, mpierr
   integer(kind=c_int)                       :: nbw, num_blocks
   real(kind=c_double), allocatable          :: tmat(:,:,:), e(:)
   integer(kind=c_intptr_t)                  :: tmat_dev, q_dev, a_dev
   real(kind=c_double)                       :: ttt0, ttt1
   logical                                   :: success
   logical                                   :: wantDebug
   integer(kind=c_int)                       :: istat
   character(200)                            :: errorMessage
   logical                                   :: do_useGPU, do_useGPU_4

    call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
    call mpi_comm_size(mpi_comm_all,n_pes,mpierr)
    call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
    call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
    call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

    wantDebug   = .false.
    success     = .true.
    useQRActual = .false.
    do_useGPU   = .false.
    do_useGPU_4 = .false.

    THIS_REAL_ELPA_KERNEL = get_actual_real_kernel()

    if(THIS_REAL_ELPA_KERNEL .eq. REAL_ELPA_KERNEL_GPU) then
       if(nblk .ne. 128) then
          print *, " ERROR: ELPA2 GPU kernel requires blocksize = 128. Exiting.."
          stop
       endif

! 4th step
       do_useGPU_4 = .true.
    endif

! Bandwidth must be a multiple of nblk
! Set to a value >= 32
    if(do_useGPU) then
       nbw = nblk
    else
       nbw = (63/nblk+1)*nblk
    endif
    num_blocks = (na-1)/nbw+1

    allocate(tmat(nbw,nbw,num_blocks),stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_real_2stage: error when allocating tmat "//errorMessage
       stop
    endif

! Reduction full -> band
    ttt0 = MPI_Wtime()

    call bandred_real_double(na,a,a_dev,lda,nblk,nbw,matrixCols,num_blocks,&
            mpi_comm_rows,mpi_comm_cols,tmat,tmat_dev,wantDebug,do_useGPU,success,useQRActual)

    if(.not.(success)) return
    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time full => band           :",ttt1-ttt0," s"

! Reduction band -> tridiagonal
    allocate(e(na),stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_real_2stage: error when allocating e "//errorMessage
       stop
    endif

    ttt0 = MPI_Wtime()

    call tridiag_band_real_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_real,&
            mpi_comm_rows,mpi_comm_cols,mpi_comm_all)

    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time band => tridiagonal    :",ttt1-ttt0," s"

    call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
    call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

! Solve tridiagonal system
    ttt0 = MPI_Wtime()

    call solve_tridi_double(na,nev,ev,e,q,ldq,nblk,matrixCols,mpi_comm_rows,&
            mpi_comm_cols,wantDebug,success)

    if(.not.(success)) return
    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"

    deallocate(e,stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_real_2stage: error when deallocating e "//errorMessage
       stop
    endif

! Backtransform stage 1
    ttt0 = MPI_Wtime()

    call trans_ev_tridi_to_band_real_double(na,nev,nblk,nbw,q,q_dev,ldq,matrixCols,hh_trans_real,&
            mpi_comm_rows,mpi_comm_cols,wantDebug,do_useGPU_4,success,THIS_REAL_ELPA_KERNEL)

    if(.not.(success)) return
    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time ev tridiagonal => band :",ttt1-ttt0," s"

! We can now deallocate the stored householder vectors
    deallocate(hh_trans_real,stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_real_2stage: error when deallocating hh_trans_real "//errorMessage
       stop
    endif

! Backtransform stage 2
    ttt0 = MPI_Wtime()

    call trans_ev_band_to_full_real_double(na,nev,nblk,nbw,a,a_dev,lda,tmat,tmat_dev,q,q_dev,ldq,&
            matrixCols,num_blocks,mpi_comm_rows,mpi_comm_cols,do_useGPU,useQRActual)

    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time ev band => full        :",ttt1-ttt0," s"

    deallocate(tmat,stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_real_2stage: error when deallocating tmat"//errorMessage
       stop
    endif

    if(my_prow==0 .and. my_pcol==0) then
       if(do_useGPU) then
          if(do_useGPU_4) then
             print *,"  GPU has been used for the 1st, 4th, 5th steps of ELPA2"
          else
             print *,"  GPU has been used for the 1st, 5th steps of ELPA2"
          endif
       endif
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
!>  \param matrixCols                           local columns of matrix a and q
!>
!>  \param mpi_comm_rows                        MPI communicator for rows
!>  \param mpi_comm_cols                        MPI communicator for columns
!>  \param mpi_comm_all                         MPI communicator for the total processor set
!>
!>  \param THIS_REAL_ELPA_KERNEL_API (optional) specify used ELPA2 kernel via API
!>  \param useGPU (optional)                    decide whether GPUs should be used or not
!>
!>
!>  \result success                             logical, false if error occured
!-------------------------------------------------------------------------------

function solve_evp_complex_2stage_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
            mpi_comm_rows,mpi_comm_cols,mpi_comm_all,THIS_COMPLEX_ELPA_KERNEL_API,&
            useGPU) result(success)

    use elpa1_compute
    use elpa2_compute
    use elpa_mpi
    use, intrinsic :: iso_c_binding

    implicit none

    logical, intent(in), optional             :: useGPU
    integer(kind=c_int), intent(in), optional :: THIS_COMPLEX_ELPA_KERNEL_API
    integer(kind=c_int)                       :: THIS_COMPLEX_ELPA_KERNEL
    integer(kind=c_int), intent(in)           :: na, nev, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
    real(kind=c_double), intent(inout)        :: ev(na)

    complex(kind=c_double), intent(inout)     :: a(lda,*), q(ldq,*)

    complex(kind=c_double), allocatable       :: hh_trans_complex(:,:)

    integer(kind=c_int)                       :: my_prow, my_pcol, np_rows, np_cols, mpierr, my_pe, n_pes
    integer(kind=c_int)                       :: l_cols, l_rows, l_cols_nev, nbw, num_blocks
    complex(kind=c_double), allocatable       :: tmat(:,:,:)
    real(kind=c_double), allocatable          :: q_real(:,:), e(:)
    real(kind=c_double)                       :: ttt0, ttt1

    logical                                   :: success, wantDebug
    integer(kind=c_int)                       :: istat
    character(200)                            :: errorMessage
    logical                                   :: do_useGPU, do_useGPU_4

    call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
    call mpi_comm_size(mpi_comm_all,n_pes,mpierr)
    call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
    call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
    call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

    do_useGPU   = .false.
    do_useGPU_4 = .false.
    wantDebug   = .false.
    success     = .true.

    THIS_COMPLEX_ELPA_KERNEL = get_actual_complex_kernel()

    if(THIS_COMPLEX_ELPA_KERNEL .eq. COMPLEX_ELPA_KERNEL_GPU) then
       if(nblk .ne. 128) then
          print *, " ERROR: ELPA2 GPU kernel requires blocksize = 128. Exiting.."
          stop
       endif

! 4th step
       do_useGPU_4 = .true.
    endif

! Choose bandwidth, must be a multiple of nblk, set to a value >= 32
    nbw = (31/nblk+1)*nblk
    num_blocks = (na-1)/nbw + 1

    allocate(tmat(nbw,nbw,num_blocks),stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_complex_2stage: error when allocating tmat"//errorMessage
       stop
    endif

! Reduction full -> band
    ttt0 = MPI_Wtime()

    call bandred_complex_double(na,a,lda,nblk,nbw,matrixCols,num_blocks,&
            mpi_comm_rows,mpi_comm_cols,tmat,wantDebug,do_useGPU,success)

    if(.not.(success)) then

       return
    endif
    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time full => band           :",ttt1-ttt0," s"

! Reduction band -> tridiagonal
    allocate(e(na),stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_complex_2stage: error when allocating e"//errorMessage
       stop
    endif

    ttt0 = MPI_Wtime()

    call tridiag_band_complex_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_complex,&
            mpi_comm_rows,mpi_comm_cols,mpi_comm_all)

    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time band => tridiagonal    :",ttt1-ttt0," s"

    call mpi_bcast(ev,na,mpi_real8,0,mpi_comm_all,mpierr)
    call mpi_bcast(e,na,mpi_real8,0,mpi_comm_all,mpierr)

    l_rows = local_index(na,my_prow,np_rows,nblk,-1) ! Local rows of a and q
    l_cols = local_index(na,my_pcol,np_cols,nblk,-1) ! Local columns of q
    l_cols_nev = local_index(nev,my_pcol,np_cols,nblk,-1) ! Local columns corresponding to nev

    allocate(q_real(l_rows,l_cols),stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_complex_2stage: error when allocating q_real"//errorMessage
       stop
    endif

! Solve tridiagonal system
    ttt0 = MPI_Wtime()

    call solve_tridi_double(na,nev,ev,e,q_real,ubound(q_real,dim=1),nblk,matrixCols,&
            mpi_comm_rows,mpi_comm_cols,wantDebug,success)

    if(.not.(success)) return

    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times)  &
       write(use_unit,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"

    q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

    deallocate(e,q_real,stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_complex_2stage: error when deallocating e, q_real"//errorMessage
       stop
    endif

! Backtransform stage 1
    ttt0 = MPI_Wtime()

    call trans_ev_tridi_to_band_complex_double(na,nev,nblk,nbw,q,ldq,matrixCols,hh_trans_complex,&
            mpi_comm_rows,mpi_comm_cols,wantDebug,do_useGPU_4,success,THIS_COMPLEX_ELPA_KERNEL)

    if(.not.(success)) return
    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time ev tridiagonal => band :",ttt1-ttt0," s"

! We can now deallocate the stored householder vectors
    deallocate(hh_trans_complex,stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_complex_2stage: error when deallocating hh_trans_complex"//errorMessage
       stop
    endif

! Backtransform stage 2
    ttt0 = MPI_Wtime()

    call trans_ev_band_to_full_complex_double(na,nev,nblk,nbw,a,lda,tmat,q,ldq,matrixCols,num_blocks,&
            mpi_comm_rows,mpi_comm_cols,do_useGPU)

    ttt1 = MPI_Wtime()
    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
       write(use_unit,"(A,F10.3,A)") "  | Time ev band => full        :",ttt1-ttt0," s"

    deallocate(tmat,stat=istat,errmsg=errorMessage)
    if(istat .ne. 0) then
       print *,"solve_evp_complex_2stage: error when deallocating tmat "//errorMessage
       stop
    endif

    if(my_prow==0 .and. my_pcol==0) then
       if(do_useGPU) then
          if(do_useGPU_4) then
             print *,"  GPU has been used for the 1st, 4th, 5th steps of ELPA2"
          else
             print *,"  GPU has been used for the 1st, 5th steps of ELPA2"
          endif
       endif
    endif

end function solve_evp_complex_2stage_double

end module ELPA2
