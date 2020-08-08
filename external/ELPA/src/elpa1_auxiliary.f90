!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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
!
! This file has been rewritten by A. Marek, MPCDF

!> \brief Fortran module which provides helper routines for matrix calculations
module elpa1_auxiliary_impl
   use elpa_utilities
   use elpa_abstract_impl
   use elpa1_compute, solve_tridi_&
   &double&
   &_private_impl => solve_tridi_&
   &double&
   &_impl
   use elpa1_compute, solve_tridi_&
   &single&
   &_private_impl => solve_tridi_&
   &single&
   &_impl
   use elpa_mpi
   use precision
   use elpa_omp
   use, intrinsic :: iso_c_binding
   use cuda_functions
   use mod_check_for_gpu

   implicit none

   public :: elpa_mult_at_b_real_double_impl      !< Multiply double-precision real matrices A**T * B

   public :: elpa_mult_ah_b_complex_double_impl   !< Multiply double-precision complex matrices A**H * B

   public :: elpa_invert_trm_real_double_impl    !< Invert double-precision real triangular matrix

   public :: elpa_invert_trm_complex_double_impl  !< Invert double-precision complex triangular matrix

   public :: elpa_cholesky_real_double_impl       !< Cholesky factorization of a double-precision real matrix

   public :: elpa_cholesky_complex_double_impl    !< Cholesky factorization of a double-precision complex matrix

   public :: elpa_solve_tridi_double_impl         !< Solve tridiagonal eigensystem for a double-precision matrix with divide and conquer method

   public :: elpa_cholesky_real_single_impl       !< Cholesky factorization of a single-precision real matrix
   public :: elpa_invert_trm_real_single_impl     !< Invert single-precision real triangular matrix
   public :: elpa_mult_at_b_real_single_impl      !< Multiply single-precision real matrices A**T * B
   public :: elpa_solve_tridi_single_impl         !< Solve tridiagonal eigensystem for a single-precision matrix with divide and conquer method

   public :: elpa_cholesky_complex_single_impl    !< Cholesky factorization of a single-precision complex matrix
   public :: elpa_invert_trm_complex_single_impl  !< Invert single-precision complex triangular matrix
   public :: elpa_mult_ah_b_complex_single_impl   !< Multiply single-precision complex matrices A**H * B

contains

   function elpa_cholesky_real_double_impl (obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none
      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: rck = C_DOUBLE
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      real(kind=rck)      :: a(obj%local_nrows,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)              :: n, nc, i, info
      integer(kind=BLAS_KIND)       :: infoBLAS
      integer(kind=ik)              :: lcs, lce, lrs, lre
      integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

      real(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
      logical                       :: wantDebug
      logical                       :: success
      integer(kind=ik)              :: istat, debug, error
      character(200)                :: errorMessage
      integer(kind=ik)              :: nrThreads

      call obj%timer%start("elpa_cholesky_&
      &real&
      &_&
      &double&
      &")

      nrThreads=1

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error )
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug",debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for debug settings. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI, kind=c_int)
      np_rows = int(np_rowsMPI, kind=c_int)
      my_pcol = int(my_pcolMPI, kind=c_int)
      np_cols = int(np_colsMPI, kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

      tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
      tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

      l_rows_tile = tile_size/np_rows ! local rows of a tile
      l_cols_tile = tile_size/np_cols ! local cols of a tile

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp1", 149,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp2", 152,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatr", 158,  istat,  errorMessage)

      allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatc", 161,  istat,  errorMessage)

      tmatr = 0
      tmatc = 0

      do n = 1, na, nblk
         ! Calculate first local row and column of the still remaining matrix
         ! on the local processor

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

         if (n+nblk > na) then

            ! This is the last step, just do a Cholesky-Factorization
            ! of the remaining block

            if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call DPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND), infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &real&

                  &: Error in dpotrf: ",info
                  success = .false.
                  return
               endif

            endif

            exit ! Loop

         endif

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then

               ! The process owning the upper left remaining block does the
               ! Cholesky-Factorization of this block
               call obj%timer%start("blas")

               call DPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND) , infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &real&

                  &: Error in dpotrf 2: ",info
                  success = .false.
                  return
               endif

               nc = 0
               do i=1,nblk
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")

            call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
               MPI_REAL8,         &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")

            nc = 0
            do i=1,nblk
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call DTRSM('L', 'U', 'T', 'N', int(nblk,kind=BLAS_KIND),  &
               int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
               int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
            call obj%timer%stop("blas")
         endif

         do i=1,nblk

            if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)

            call obj%timer%start("mpi_communication")
            if (l_cols-l_colx+1>0) &
               call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_REAL8, &
               int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")
         enddo
         ! this has to be checked since it was changed substantially when doing type safe
         call elpa_transpose_vectors_&
         &real&
         &_&
         &double &
            (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
            tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
            n, na, nblk, nblk, nrThreads)

         do i=0,(na-1)/tile_size
            lcs = max(l_colx,i*l_cols_tile+1)
            lce = min(l_cols,(i+1)*l_cols_tile)
            lrs = l_rowx
            lre = min(l_rows,(i+1)*l_rows_tile)
            if (lce<lcs .or. lre<lrs) cycle
            call obj%timer%start("blas")
            call DGEMM('N', 'T', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
               int(nblk,kind=BLAS_KIND), -ONE,  &
               tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
               int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
               ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")

         enddo

      enddo

      deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 319,  istat,  errorMessage)

      ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

      do i=1,na
         if (my_pcol==pcol(i, nblk, np_cols)) then
            ! column i is on local processor
            l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
            l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
            a(l_row1:l_rows,l_col1) = 0
         endif
      enddo

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_cholesky_&
      &real&
      &_&
      &double&
      &")

   end function elpa_cholesky_real_double_impl

   function elpa_cholesky_real_single_impl(obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none
      integer, parameter :: rk = C_FLOAT
      integer, parameter :: rck = C_FLOAT
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      real(kind=rck)      :: a(obj%local_nrows,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)              :: n, nc, i, info
      integer(kind=BLAS_KIND)       :: infoBLAS
      integer(kind=ik)              :: lcs, lce, lrs, lre
      integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

      real(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
      logical                       :: wantDebug
      logical                       :: success
      integer(kind=ik)              :: istat, debug, error
      character(200)                :: errorMessage
      integer(kind=ik)              :: nrThreads

      call obj%timer%start("elpa_cholesky_&
      &real&
      &_&
      &single&
      &")

      nrThreads=1

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error )
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug",debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for debug settings. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI, kind=c_int)
      np_rows = int(np_rowsMPI, kind=c_int)
      my_pcol = int(my_pcolMPI, kind=c_int)
      np_cols = int(np_colsMPI, kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

      tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
      tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

      l_rows_tile = tile_size/np_rows ! local rows of a tile
      l_cols_tile = tile_size/np_cols ! local cols of a tile

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp1", 149,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp2", 152,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatr", 158,  istat,  errorMessage)

      allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatc", 161,  istat,  errorMessage)

      tmatr = 0
      tmatc = 0

      do n = 1, na, nblk
         ! Calculate first local row and column of the still remaining matrix
         ! on the local processor

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

         if (n+nblk > na) then

            ! This is the last step, just do a Cholesky-Factorization
            ! of the remaining block

            if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call SPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND), infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &real&

                  &: Error in dpotrf: ",info
                  success = .false.
                  return
               endif

            endif

            exit ! Loop

         endif

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then

               ! The process owning the upper left remaining block does the
               ! Cholesky-Factorization of this block
               call obj%timer%start("blas")

               call SPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND) , infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &real&

                  &: Error in dpotrf 2: ",info
                  success = .false.
                  return
               endif

               nc = 0
               do i=1,nblk
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")

            call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
               MPI_REAL4,         &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")

            nc = 0
            do i=1,nblk
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call STRSM('L', 'U', 'T', 'N', int(nblk,kind=BLAS_KIND),  &
               int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
               int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
            call obj%timer%stop("blas")
         endif

         do i=1,nblk

            if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)

            call obj%timer%start("mpi_communication")
            if (l_cols-l_colx+1>0) &
               call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_REAL4, &
               int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")
         enddo
         ! this has to be checked since it was changed substantially when doing type safe
         call elpa_transpose_vectors_&
         &real&
         &_&
         &single &
            (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
            tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
            n, na, nblk, nblk, nrThreads)

         do i=0,(na-1)/tile_size
            lcs = max(l_colx,i*l_cols_tile+1)
            lce = min(l_cols,(i+1)*l_cols_tile)
            lrs = l_rowx
            lre = min(l_rows,(i+1)*l_rows_tile)
            if (lce<lcs .or. lre<lrs) cycle
            call obj%timer%start("blas")
            call SGEMM('N', 'T', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
               int(nblk,kind=BLAS_KIND), -ONE,  &
               tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
               int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
               ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")

         enddo

      enddo

      deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 319,  istat,  errorMessage)

      ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

      do i=1,na
         if (my_pcol==pcol(i, nblk, np_cols)) then
            ! column i is on local processor
            l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
            l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
            a(l_row1:l_rows,l_col1) = 0
         endif
      enddo

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_cholesky_&
      &real&
      &_&
      &single&
      &")

   end function elpa_cholesky_real_single_impl

!> \brief  elpa_invert_trm_real_double: Inverts a double-precision real upper triangular matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_invert_trm_real_double_impl(obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none
      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: rck = C_DOUBLE
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      real(kind=rck)     :: a(obj%local_nrows,*)

      integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)             :: n, nc, i, info, ns, nb
      integer(kind=BLAS_KIND)      :: infoBLAS
      real(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
      logical                      :: wantDebug
      logical                      :: success
      integer(kind=ik)             :: istat, debug, error
      character(200)               :: errorMessage

      call obj%timer%start("elpa_invert_trm_&
      &real&
      &_&
      &double&
      &")

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for debug. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .true.
      endif
      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp1", 133,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp2", 136,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat1", 142,  istat,  errorMessage)

      allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat2", 145,  istat,  errorMessage)

      tmat1 = 0
      tmat2 = 0

      ns = ((na-1)/nblk)*nblk + 1

      do n = ns,1,-nblk

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         nb = nblk
         if (na-n+1 < nblk) nb = na-n+1

         l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call DTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                  infoBLAS)
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
                  &real&

                  &: Error in DTRTRI"

                  success = .false.
                  call obj%timer%stop("elpa_invert_trm_&
                  &real&
                  &_&
                  &double&
                  &")
                  return
               endif

               nc = 0
               do i=1,nb
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_REAL8,       &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            nc = 0
            do i=1,nb
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call DTRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
               tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
            if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

         endif

         if (l_row1>1) then
            if (my_pcol==pcol(n, nblk, np_cols)) then
               tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
               a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
            endif

            do i=1,nb
               call obj%timer%start("mpi_communication")
               call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_REAL8, &
                  int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

               call obj%timer%stop("mpi_communication")
            enddo
         endif
         call obj%timer%start("mpi_communication")
         if (l_cols-l_col1+1>0) &
            call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_REAL8, &
            int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

         call obj%timer%stop("mpi_communication")

         call obj%timer%start("blas")
         if (l_row1>1 .and. l_cols-l_col1+1>0) &
            call DGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
            int(nb,kind=BLAS_KIND), -ONE, &
            tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
            int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
            a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

         call obj%timer%stop("blas")

      enddo

      deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 260,  istat,  errorMessage)

      call obj%timer%stop("elpa_invert_trm_&
      &real&
      &_&
      &double&
      &")
   end function elpa_invert_trm_real_double_impl

!> \brief  elpa_invert_trm_real_single_impl: Inverts a single-precision real upper triangular matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure

   function elpa_invert_trm_real_single_impl(obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none
      integer, parameter :: rk = C_FLOAT
      integer, parameter :: rck = C_FLOAT
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      real(kind=rck)     :: a(obj%local_nrows,*)

      integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)             :: n, nc, i, info, ns, nb
      integer(kind=BLAS_KIND)      :: infoBLAS
      real(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
      logical                      :: wantDebug
      logical                      :: success
      integer(kind=ik)             :: istat, debug, error
      character(200)               :: errorMessage

      call obj%timer%start("elpa_invert_trm_&
      &real&
      &_&
      &single&
      &")

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for debug. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .true.
      endif
      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp1", 133,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp2", 136,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat1", 142,  istat,  errorMessage)

      allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat2", 145,  istat,  errorMessage)

      tmat1 = 0
      tmat2 = 0

      ns = ((na-1)/nblk)*nblk + 1

      do n = ns,1,-nblk

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         nb = nblk
         if (na-n+1 < nblk) nb = na-n+1

         l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call STRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                  infoBLAS)
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
                  &real&

                  &: Error in DTRTRI"

                  success = .false.
                  call obj%timer%stop("elpa_invert_trm_&
                  &real&
                  &_&
                  &single&
                  &")
                  return
               endif

               nc = 0
               do i=1,nb
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_REAL4,       &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            nc = 0
            do i=1,nb
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call STRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
               tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
            if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

         endif

         if (l_row1>1) then
            if (my_pcol==pcol(n, nblk, np_cols)) then
               tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
               a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
            endif

            do i=1,nb
               call obj%timer%start("mpi_communication")
               call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_REAL4, &
                  int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

               call obj%timer%stop("mpi_communication")
            enddo
         endif
         call obj%timer%start("mpi_communication")
         if (l_cols-l_col1+1>0) &
            call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_REAL4, &
            int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

         call obj%timer%stop("mpi_communication")

         call obj%timer%start("blas")
         if (l_row1>1 .and. l_cols-l_col1+1>0) &
            call SGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
            int(nb,kind=BLAS_KIND), -ONE, &
            tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
            int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
            a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

         call obj%timer%stop("blas")

      enddo

      deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 260,  istat,  errorMessage)

      call obj%timer%stop("elpa_invert_trm_&
      &real&
      &_&
      &single&
      &")
   end function elpa_invert_trm_real_single_impl

!> \brief  elpa_cholesky_complex_double_impl: Cholesky factorization of a double-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_complex_double_impl(obj, a) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none

      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: ck = C_DOUBLE_COMPLEX
      integer, parameter :: rck = C_DOUBLE_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      complex(kind=rck)      :: a(obj%local_nrows,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)              :: n, nc, i, info
      integer(kind=BLAS_KIND)       :: infoBLAS
      integer(kind=ik)              :: lcs, lce, lrs, lre
      integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

      complex(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
      logical                       :: wantDebug
      logical                       :: success
      integer(kind=ik)              :: istat, debug, error
      character(200)                :: errorMessage
      integer(kind=ik)              :: nrThreads

      call obj%timer%start("elpa_cholesky_&
      &complex&
      &_&
      &double&
      &")

      nrThreads=1

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error )
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug",debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for debug settings. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI, kind=c_int)
      np_rows = int(np_rowsMPI, kind=c_int)
      my_pcol = int(my_pcolMPI, kind=c_int)
      np_cols = int(np_colsMPI, kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

      tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
      tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

      l_rows_tile = tile_size/np_rows ! local rows of a tile
      l_cols_tile = tile_size/np_cols ! local cols of a tile

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp1", 149,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp2", 152,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatr", 158,  istat,  errorMessage)

      allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatc", 161,  istat,  errorMessage)

      tmatr = 0
      tmatc = 0

      do n = 1, na, nblk
         ! Calculate first local row and column of the still remaining matrix
         ! on the local processor

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

         if (n+nblk > na) then

            ! This is the last step, just do a Cholesky-Factorization
            ! of the remaining block

            if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call ZPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND), infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &complex&

                  &: Error in zpotrf: ",info
                  success = .false.
                  return
               endif

            endif

            exit ! Loop

         endif

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then

               ! The process owning the upper left remaining block does the
               ! Cholesky-Factorization of this block
               call obj%timer%start("blas")

               call ZPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND) , infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &complex&

                  &: Error in zpotrf 2: ",info

                  success = .false.
                  return
               endif

               nc = 0
               do i=1,nblk
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")

            call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
               MPI_DOUBLE_COMPLEX,      &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")

            nc = 0
            do i=1,nblk
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call ZTRSM('L', 'U', 'C', 'N', int(nblk,kind=BLAS_KIND),  &
               int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
               int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
            call obj%timer%stop("blas")
         endif

         do i=1,nblk

            if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))

            call obj%timer%start("mpi_communication")
            if (l_cols-l_colx+1>0) &
               call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
               int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")
         enddo
         ! this has to be checked since it was changed substantially when doing type safe
         call elpa_transpose_vectors_&
         &complex&
         &_&
         &double &
            (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
            tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
            n, na, nblk, nblk, nrThreads)

         do i=0,(na-1)/tile_size
            lcs = max(l_colx,i*l_cols_tile+1)
            lce = min(l_cols,(i+1)*l_cols_tile)
            lrs = l_rowx
            lre = min(l_rows,(i+1)*l_rows_tile)
            if (lce<lcs .or. lre<lrs) cycle
            call obj%timer%start("blas")
            call ZGEMM('N', 'C', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
               int(nblk,kind=BLAS_KIND), -ONE,  &
               tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
               int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
               ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")

         enddo

      enddo

      deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 319,  istat,  errorMessage)

      ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

      do i=1,na
         if (my_pcol==pcol(i, nblk, np_cols)) then
            ! column i is on local processor
            l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
            l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
            a(l_row1:l_rows,l_col1) = 0
         endif
      enddo

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_cholesky_&
      &complex&
      &_&
      &double&
      &")

   end function elpa_cholesky_complex_double_impl

!> \brief  elpa_cholesky_complex_single_impl: Cholesky factorization of a single-precision complex hermitian matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_cholesky_complex_single_impl(obj, a) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none

      integer, parameter :: rk = C_FLOAT
      integer, parameter :: ck = C_FLOAT_COMPLEX
      integer, parameter :: rck = C_FLOAT_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)              :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      complex(kind=rck)      :: a(obj%local_nrows,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)              :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)              :: n, nc, i, info
      integer(kind=BLAS_KIND)       :: infoBLAS
      integer(kind=ik)              :: lcs, lce, lrs, lre
      integer(kind=ik)              :: tile_size, l_rows_tile, l_cols_tile

      complex(kind=rck), allocatable    :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)
      logical                       :: wantDebug
      logical                       :: success
      integer(kind=ik)              :: istat, debug, error
      character(200)                :: errorMessage
      integer(kind=ik)              :: nrThreads

      call obj%timer%start("elpa_cholesky_&
      &complex&
      &_&
      &single&
      &")

      nrThreads=1

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error )
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug",debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for debug settings. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI, kind=c_int)
      np_rows = int(np_rowsMPI, kind=c_int)
      my_pcol = int(my_pcolMPI, kind=c_int)
      np_cols = int(np_colsMPI, kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

      tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
      tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

      l_rows_tile = tile_size/np_rows ! local rows of a tile
      l_cols_tile = tile_size/np_cols ! local cols of a tile

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp1", 149,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmp2", 152,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatr", 158,  istat,  errorMessage)

      allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_cholesky: tmatc", 161,  istat,  errorMessage)

      tmatr = 0
      tmatc = 0

      do n = 1, na, nblk
         ! Calculate first local row and column of the still remaining matrix
         ! on the local processor

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         l_rowx = local_index(n+nblk, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nblk, my_pcol, np_cols, nblk, +1)

         if (n+nblk > na) then

            ! This is the last step, just do a Cholesky-Factorization
            ! of the remaining block

            if (my_prow==prow(n, nblk, np_rows) .and. my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call CPOTRF('U', int(na-n+1,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND), infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &complex&

                  &: Error in zpotrf: ",info
                  success = .false.
                  return
               endif

            endif

            exit ! Loop

         endif

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then

               ! The process owning the upper left remaining block does the
               ! Cholesky-Factorization of this block
               call obj%timer%start("blas")

               call CPOTRF('U', int(nblk,kind=BLAS_KIND), a(l_row1,l_col1), &
                  int(matrixRows,kind=BLAS_KIND) , infoBLAS )
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_cholesky_&
                  &complex&

                  &: Error in zpotrf 2: ",info

                  success = .false.
                  return
               endif

               nc = 0
               do i=1,nblk
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")

            call MPI_Bcast(tmp1, int(nblk*(nblk+1)/2,kind=MPI_KIND),      &
               MPI_COMPLEX,      &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")

            nc = 0
            do i=1,nblk
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call CTRSM('L', 'U', 'C', 'N', int(nblk,kind=BLAS_KIND),  &
               int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, tmp2, &
               int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND) )
            call obj%timer%stop("blas")
         endif

         do i=1,nblk

            if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))

            call obj%timer%start("mpi_communication")
            if (l_cols-l_colx+1>0) &
               call MPI_Bcast(tmatc(l_colx,i), int(l_cols-l_colx+1,kind=MPI_KIND), MPI_COMPLEX, &
               int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

            call obj%timer%stop("mpi_communication")
         enddo
         ! this has to be checked since it was changed substantially when doing type safe
         call elpa_transpose_vectors_&
         &complex&
         &_&
         &single &
            (obj, tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
            tmatr, ubound(tmatr,dim=1), mpi_comm_rows, &
            n, na, nblk, nblk, nrThreads)

         do i=0,(na-1)/tile_size
            lcs = max(l_colx,i*l_cols_tile+1)
            lce = min(l_cols,(i+1)*l_cols_tile)
            lrs = l_rowx
            lre = min(l_rows,(i+1)*l_rows_tile)
            if (lce<lcs .or. lre<lrs) cycle
            call obj%timer%start("blas")
            call CGEMM('N', 'C', int(lre-lrs+1,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
               int(nblk,kind=BLAS_KIND), -ONE,  &
               tmatr(lrs,1), int(ubound(tmatr,dim=1),kind=BLAS_KIND), tmatc(lcs,1), &
               int(ubound(tmatc,dim=1),kind=BLAS_KIND), &
               ONE, a(lrs,lcs), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")

         enddo

      enddo

      deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_cholesky: tmp1, tmp2, tmatr, tmatc", 319,  istat,  errorMessage)

      ! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

      do i=1,na
         if (my_pcol==pcol(i, nblk, np_cols)) then
            ! column i is on local processor
            l_col1 = local_index(i  , my_pcol, np_cols, nblk, +1) ! local column number
            l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
            a(l_row1:l_rows,l_col1) = 0
         endif
      enddo

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_cholesky_&
      &complex&
      &_&
      &single&
      &")

   end function elpa_cholesky_complex_single_impl

!> \brief  elpa_invert_trm_complex_double_impl: Inverts a double-precision complex upper triangular matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_invert_trm_complex_double_impl(obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none

      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: ck = C_DOUBLE_COMPLEX
      integer, parameter :: rck = C_DOUBLE_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      complex(kind=rck)     :: a(obj%local_nrows,*)

      integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)             :: n, nc, i, info, ns, nb
      integer(kind=BLAS_KIND)      :: infoBLAS
      complex(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
      logical                      :: wantDebug
      logical                      :: success
      integer(kind=ik)             :: istat, debug, error
      character(200)               :: errorMessage

      call obj%timer%start("elpa_invert_trm_&
      &complex&
      &_&
      &double&
      &")

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for debug. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .true.
      endif
      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp1", 133,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp2", 136,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat1", 142,  istat,  errorMessage)

      allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat2", 145,  istat,  errorMessage)

      tmat1 = 0
      tmat2 = 0

      ns = ((na-1)/nblk)*nblk + 1

      do n = ns,1,-nblk

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         nb = nblk
         if (na-n+1 < nblk) nb = na-n+1

         l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call ZTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                  infoBLAS)
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
                  &complex&

                  &: Error in ZTRTRI"

                  success = .false.
                  call obj%timer%stop("elpa_invert_trm_&
                  &complex&
                  &_&
                  &double&
                  &")
                  return
               endif

               nc = 0
               do i=1,nb
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,       &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            nc = 0
            do i=1,nb
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call ZTRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
               tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
            if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

         endif

         if (l_row1>1) then
            if (my_pcol==pcol(n, nblk, np_cols)) then
               tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
               a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
            endif

            do i=1,nb
               call obj%timer%start("mpi_communication")
               call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                  int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

               call obj%timer%stop("mpi_communication")
            enddo
         endif
         call obj%timer%start("mpi_communication")
         if (l_cols-l_col1+1>0) &
            call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
            int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

         call obj%timer%stop("mpi_communication")

         call obj%timer%start("blas")
         if (l_row1>1 .and. l_cols-l_col1+1>0) &
            call ZGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
            int(nb,kind=BLAS_KIND), -ONE, &
            tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
            int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
            a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

         call obj%timer%stop("blas")

      enddo

      deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 260,  istat,  errorMessage)

      call obj%timer%stop("elpa_invert_trm_&
      &complex&
      &_&
      &double&
      &")
   end function elpa_invert_trm_complex_double_impl

!> \brief  elpa_invert_trm_complex_single_impl: Inverts a single-precision complex upper triangular matrix
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%local_nrows   Leading dimension of a
!> \param     - obj%local_ncols   local columns of matrix a
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param  a(lda,matrixCols)      Distributed matrix which should be inverted
!>                                Distribution is like in Scalapack.
!>                                Only upper triangle needs to be set.
!>                                The lower triangle is not referenced.
!> \result succes                 logical, reports success or failure
   function elpa_invert_trm_complex_single_impl(obj, a) result(success)
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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

      implicit none

      integer, parameter :: rk = C_FLOAT
      integer, parameter :: ck = C_FLOAT_COMPLEX
      integer, parameter :: rck = C_FLOAT_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)             :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      complex(kind=rck)     :: a(obj%local_nrows,*)

      integer(kind=ik)             :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)       :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik)             :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
      integer(kind=ik)             :: n, nc, i, info, ns, nb
      integer(kind=BLAS_KIND)      :: infoBLAS
      complex(kind=rck), allocatable   :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)
      logical                      :: wantDebug
      logical                      :: success
      integer(kind=ik)             :: istat, debug, error
      character(200)               :: errorMessage

      call obj%timer%start("elpa_invert_trm_&
      &complex&
      &_&
      &single&
      &")

      na         = obj%na
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Error getting option for debug. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .true.
      endif
      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

      allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp1", 133,  istat,  errorMessage)

      allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmp2", 136,  istat,  errorMessage)

      tmp1 = 0
      tmp2 = 0

      allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat1", 142,  istat,  errorMessage)

      allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_invert_trm: tmat2", 145,  istat,  errorMessage)

      tmat1 = 0
      tmat2 = 0

      ns = ((na-1)/nblk)*nblk + 1

      do n = ns,1,-nblk

         l_row1 = local_index(n, my_prow, np_rows, nblk, +1)
         l_col1 = local_index(n, my_pcol, np_cols, nblk, +1)

         nb = nblk
         if (na-n+1 < nblk) nb = na-n+1

         l_rowx = local_index(n+nb, my_prow, np_rows, nblk, +1)
         l_colx = local_index(n+nb, my_pcol, np_cols, nblk, +1)

         if (my_prow==prow(n, nblk, np_rows)) then

            if (my_pcol==pcol(n, nblk, np_cols)) then
               call obj%timer%start("blas")

               call CTRTRI('U', 'N', int(nb,kind=BLAS_KIND), a(l_row1,l_col1), int(matrixRows,kind=BLAS_KIND), &
                  infoBLAS)
               info = int(infoBLAS,kind=ik)
               call obj%timer%stop("blas")

               if (info/=0) then
                  if (wantDebug) write(error_unit,*) "elpa_invert_trm_&
                  &complex&

                  &: Error in ZTRTRI"

                  success = .false.
                  call obj%timer%stop("elpa_invert_trm_&
                  &complex&
                  &_&
                  &single&
                  &")
                  return
               endif

               nc = 0
               do i=1,nb
                  tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
                  nc = nc+i
               enddo
            endif
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(tmp1, int(nb*(nb+1)/2,kind=MPI_KIND), MPI_COMPLEX,       &
               int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            nc = 0
            do i=1,nb
               tmp2(1:i,i) = tmp1(nc+1:nc+i)
               nc = nc+i
            enddo

            call obj%timer%start("blas")
            if (l_cols-l_colx+1>0) &
               call CTRMM('L', 'U', 'N', 'N', int(nb,kind=BLAS_KIND), int(l_cols-l_colx+1,kind=BLAS_KIND), ONE, &
               tmp2, int(ubound(tmp2,dim=1),kind=BLAS_KIND), a(l_row1,l_colx), int(matrixRows,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            if (l_colx<=l_cols)   tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
            if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

         endif

         if (l_row1>1) then
            if (my_pcol==pcol(n, nblk, np_cols)) then
               tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
               a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
            endif

            do i=1,nb
               call obj%timer%start("mpi_communication")
               call MPI_Bcast(tmat1(1,i), int(l_row1-1,kind=MPI_KIND), MPI_COMPLEX, &
                  int(pcol(n, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)

               call obj%timer%stop("mpi_communication")
            enddo
         endif
         call obj%timer%start("mpi_communication")
         if (l_cols-l_col1+1>0) &
            call MPI_Bcast(tmat2(1,l_col1), int((l_cols-l_col1+1)*nblk,kind=MPI_KIND), MPI_COMPLEX, &
            int(prow(n, nblk, np_rows),kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)

         call obj%timer%stop("mpi_communication")

         call obj%timer%start("blas")
         if (l_row1>1 .and. l_cols-l_col1+1>0) &
            call CGEMM('N', 'N', int(l_row1-1,kind=BLAS_KIND), int(l_cols-l_col1+1,kind=BLAS_KIND), &
            int(nb,kind=BLAS_KIND), -ONE, &
            tmat1, int(ubound(tmat1,dim=1),kind=BLAS_KIND), tmat2(1,l_col1), &
            int(ubound(tmat2,dim=1),kind=BLAS_KIND), ONE, &
            a(1,l_col1), int(matrixRows,kind=BLAS_KIND) )

         call obj%timer%stop("blas")

      enddo

      deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_invert_trm: tmp1, tmp2, tmat1, tmat2", 260,  istat,  errorMessage)

      call obj%timer%stop("elpa_invert_trm_&
      &complex&
      &_&
      &single&
      &")
   end function elpa_invert_trm_complex_single_impl

   function elpa_mult_at_b_real_double_impl(obj, uplo_a, uplo_c, ncb, a, b, ldb, ldbCols, &
      c, ldc, ldcCols) result(success)
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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
!
! Author: A. Marek, MPCDF

      implicit none

      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: rck = C_DOUBLE
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj

      character*1                   :: uplo_a, uplo_c

      integer(kind=ik), intent(in)  :: ldb, ldbCols, ldc, ldcCols
      integer(kind=ik)              :: na, ncb
      real(kind=rck)       :: a(obj%local_nrows,*), b(ldb,*), c(ldc,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
      integer(kind=MPI_KIND)        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=MPI_KIND)        :: mpierr, myidMPI
      integer(kind=ik)              :: l_cols, l_rows, l_rows_np
      integer(kind=ik)              :: np, n, nb, nblk_mult, lrs, lre, lcs, lce
      integer(kind=ik)              :: gcol_min, gcol, goff
      integer(kind=ik)              :: nstor, nr_done, noff, np_bc, n_aux_bc, nvals
      integer(kind=ik), allocatable :: lrs_save(:), lre_save(:)

      logical                       :: a_lower, a_upper, c_lower, c_upper
      real(kind=rck), pointer :: aux_mat(:,:), tmp1(:,:)
      real(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage
      character(20)                 :: gpuString
      logical                       :: success
      logical                       :: successCUDA
      logical                       :: useGPU
      integer(kind=c_int)           :: gpu, numGPU
      integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
      integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
      integer(kind=c_intptr_t)      :: aux_dev, a_dev, b_dev, tmp1_dev
      type(c_ptr)                   :: aux_host, tmp1_host
      integer(kind=c_intptr_t)      :: num
      integer(kind=c_intptr_t)      :: aux_off, a_off, b_off
      integer(kind=ik)              :: multiply_at_a
      integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
      &double&
      &_&
      &real

      success = .true.

      ! GPU settings
      call obj%get("gpu", gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif

      useGPU = (gpu == 1)

      if(useGPU) then
         gpuString = "_gpu"
      else
         gpuString = ""
      endif

      call obj%timer%start("elpa_mult_at_b_&
      &real&
      &_&
      &double&
      &"//gpuString)

      na          = obj%na
      nblk        = obj%nblk
      matrixRows  = obj%local_nrows
      matrixCols  = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_parent",mpi_comm_all,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_parent. Aborting..."
         stop
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_all, kind=MPI_KIND) ,myidMPI ,mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      myid    = int(myidMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")

      call obj%get("multiply_at_a",multiply_at_a,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for multiply_at_a. Aborting..."
         stop
      endif

      if (multiply_at_a == 0) then ! C = A^T * B
         l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
         l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b
      else ! C = A^T * A
         l_rows = local_index(ncb, my_prow, np_rows, nblk, -1) ! Local rows of a
         l_cols = local_index(na,  my_pcol, np_cols, nblk, -1) ! Local cols of a
      endif

      ! Block factor for matrix multiplications, must be a multiple of nblk

      if (na/np_rows <= 256) then
         nblk_mult = (31/nblk+1)*nblk
      else
         nblk_mult = (63/nblk+1)*nblk
      endif

      if (useGPU) then
         call obj%timer%start("check_for_gpu")
         if (check_for_gpu(myid,numGPU)) then
            ! set the neccessary parameters
            cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
            cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
            cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
            cudaHostRegisterPortable = cuda_hostRegisterPortable()
            cudaHostRegisterMapped   = cuda_hostRegisterMapped()
         else
            print *,"GPUs are requested but not detected! Aborting..."
            success = .false.
            return
         endif
         call obj%timer%stop("check_for_gpu")

         if (multiply_at_a == 0) then ! C = A^T * B
            ! copy b to b_dev
            num = ldb*ldbCols*size_of_datatype
            successCUDA = cuda_malloc(b_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_b: b_dev", 212,  successCUDA)

            successCUDA = cuda_host_register(int(loc(b),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_b: b", 217,  successCUDA)

            successCUDA = cuda_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_b: b to b_dev", 221,  successCUDA)
         else ! C = A^T * A
            ! copy a to a_dev
            num = matrixRows*matrixCols*size_of_datatype
            successCUDA = cuda_malloc(a_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_a: a_dev", 226,  successCUDA)

            successCUDA = cuda_host_register(int(loc(a),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_a: a", 231,  successCUDA)

            successCUDA = cuda_memcpy(a_dev,int(loc(a),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_a: a to a_dev", 235,  successCUDA)
         endif

         num = l_rows*nblk_mult*size_of_datatype
         successCUDA = cuda_malloc_host(aux_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: aux_host", 240,  successCUDA)

         call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

         successCUDA = cuda_malloc(aux_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: aux_dev", 245,  successCUDA)

         num = nblk_mult*l_cols*size_of_datatype
         successCUDA = cuda_malloc_host(tmp1_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: tmp1_host", 249,  successCUDA)

         call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

         successCUDA = cuda_malloc(tmp1_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: tmp1_dev", 254,  successCUDA)
      else ! useGPU
         allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa_mult_at_b: aux_mat", 257,  istat,  errorMessage)
      endif ! useGPU

      allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: aux_bc", 261,  istat,  errorMessage)

      allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lrs_save", 264,  istat,  errorMessage)

      allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lre_save", 267,  istat,  errorMessage)

      a_lower = .false.
      a_upper = .false.
      c_lower = .false.
      c_upper = .false.

      if (multiply_at_a == 0) then ! C = A^T * B
         if (uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
         if (uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
      endif
      if (uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
      if (uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

      ! Build up the result matrix by processor rows

      do np = 0, np_rows-1

         ! In this turn, procs of row np assemble the result

         l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

         nr_done = 0 ! Number of rows done
         aux_mat = 0
         nstor = 0   ! Number of columns stored in aux_mat

         ! Loop over the blocks on row np

         do nb=0,(l_rows_np-1)/nblk

            goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

            ! Get the processor column which owns this block (A is transposed, so we need the column)
            ! and the offset in blocks within this column.
            ! The corresponding block column in A is then broadcast to all for multiplication with B

            np_bc = MOD(goff,np_cols)
            noff = goff/np_cols
            n_aux_bc = 0

            ! Gather up the complete block column of A on the owner

            do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

               gcol = goff*nblk + n ! global column corresponding to n
               if (nstor==0 .and. n==1) gcol_min = gcol

               lrs = 1       ! 1st local row number for broadcast
               lre = l_rows  ! last local row number for broadcast
               if (a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  if (my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
                  n_aux_bc = n_aux_bc + nvals
               endif

               lrs_save(n) = lrs
               lre_save(n) = lre

            enddo

            ! Broadcast block column
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(aux_bc, int(n_aux_bc,kind=MPI_KIND),    &
               MPI_REAL8,  &
               int(np_bc,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            ! Insert what we got in aux_mat

            n_aux_bc = 0
            do n = 1, min(l_rows_np-nb*nblk,nblk)
               nstor = nstor+1
               lrs = lrs_save(n)
               lre = lre_save(n)
               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
                  n_aux_bc = n_aux_bc + nvals
               endif
            enddo

            ! If we got nblk_mult columns in aux_mat or this is the last block
            ! do the matrix multiplication

            if (nstor==nblk_mult .or. nb*nblk+nblk >= l_rows_np) then

               lrs = 1       ! 1st local row number for multiply
               lre = l_rows  ! last local row number for multiply
               if (a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               lcs = 1       ! 1st local col number for multiply
               lce = l_cols  ! last local col number for multiply
               if (c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
               if (c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

               if (lcs<=lce) then
                  allocate(tmp1(nstor,lcs:lce), tmp2(nstor,lcs:lce), stat=istat, errmsg=errorMessage)
                  call check_alloc("elpa_mult_at_b_&
                  &real ", "tmp1", istat, errorMessage)

                  if (lrs<=lre) then
                     if (useGPU) then
                        num = l_rows*nblk_mult*size_of_datatype
                        successCUDA = cuda_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                           num, cudaMemcpyHostToDevice)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: aux_mat to aux_dev", 377,  successCUDA)

                        call obj%timer%start("cublas")

                        aux_off = (lrs-1)*size_of_datatype
                        if (multiply_at_a == 0) then ! C = A^T * B
                           b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

                           call cublas_DGEMM('T', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                              tmp1_dev, nstor)
                        else ! C = A^T * A
                           a_off = ((lcs-1)*matrixRows+lrs-1)*size_of_datatype

                           call cublas_DGEMM('T', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, a_dev+a_off, matrixRows, ZERO, &
                              tmp1_dev, nstor)
                        endif

                        call obj%timer%stop("cublas")

                        num = nstor*(lce-lcs+1)*size_of_datatype
                        successCUDA = cuda_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                           tmp1_dev, num, cudaMemcpyDeviceToHost)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: tmp1_dev to tmp1", 401,  successCUDA)
                     else ! useGPU
                        call obj%timer%start("blas")

                        if (multiply_at_a == 0) then ! C = A^T * B
                           call DGEMM('T', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        else
                           call DGEMM('T', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              a(lrs,lcs), int(matrixRows,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        endif

                        call obj%timer%stop("blas")
                     endif ! useGPU
                  else
                     tmp1 = 0
                  endif

                  ! Sum up the results and send to processor row np
                  call obj%timer%start("mpi_communication")
                  call mpi_reduce(tmp1, tmp2, int(nstor*(lce-lcs+1),kind=MPI_KIND),  MPI_REAL8, &
                     MPI_SUM, int(np,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
                  call obj%timer%stop("mpi_communication")
                  ! Put the result into C
                  if (my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

                  deallocate(tmp1,tmp2, stat=istat, errmsg=errorMessage)
                  call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 442,  istat,  errorMessage)
               endif

               nr_done = nr_done+nstor
               nstor=0
               aux_mat(:,:)=0
            endif
         enddo
      enddo

      if (useGPU) then
         if (multiply_at_a == 0) then ! C = A^T * B
            successCUDA = cuda_free(b_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: b_dev", 455,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(b),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: b", 458,  successCUDA)
         else ! C = A^T * A
            successCUDA = cuda_free(a_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: a_dev", 461,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(a),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: a", 464,  successCUDA)
         endif

         nullify(aux_mat)
         nullify(tmp1)

         successCUDA = cuda_free_host(aux_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: aux_host", 471,  successCUDA)

         successCUDA = cuda_free(aux_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: aux_dev", 474,  successCUDA)

         successCUDA = cuda_free_host(tmp1_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: tmp1_host", 477,  successCUDA)

         successCUDA = cuda_free(tmp1_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: tmp1_dev", 480,  successCUDA)
      else ! useGPU
         deallocate(aux_mat, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa_mult_at_b: aux_mat", 483,  istat,  errorMessage)
      endif ! useGPU

      deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 487,  istat,  errorMessage)

      call obj%timer%stop("elpa_mult_at_b_&
      &real&
      &_&
      &double&
      &"//gpuString)
   end function elpa_mult_at_b_real_double_impl

!> \brief  elpa_mult_at_b_real_single_impl: Performs C : = A**T * B
!>         where   A is a square matrix (obj%na,obj%na) which is optionally upper or lower triangular
!>                 B is a (obj%na,ncb) matrix
!>                 C is a (obj%na,ncb) matrix where optionally only the upper or lower
!>                   triangle may be computed
!> \details

!> \param  uplo_a               'U' if A is upper triangular
!>                              'L' if A is lower triangular
!>                              anything else if A is a full matrix
!>                              Please note: This pertains to the original A (as set in the calling program)
!>                                           whereas the transpose of A is used for calculations
!>                              If uplo_a is 'U' or 'L', the other triangle is not used at all,
!>                              i.e. it may contain arbitrary numbers
!> \param uplo_c                'U' if only the upper diagonal part of C is needed
!>                              'L' if only the upper diagonal part of C is needed
!>                              anything else if the full matrix C is needed
!>                              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!>                                            written to a certain extent, i.e. one shouldn't rely on the content there!
!> \param na                    Number of rows/columns of A, number of rows of B and C
!> \param ncb                   Number of columns  of B and C
!> \param a                     matrix a
!> \param obj%local_nrows       leading dimension of matrix a, set with class method obj%set("local_nrows",value)
!> \param b                     matrix b
!> \param ldb                   leading dimension of matrix b
!> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \param c                     matrix c
!> \param ldc                   leading dimension of matrix c
!> \result success

   function elpa_mult_at_b_real_single_impl(obj, uplo_a, uplo_c, ncb, a, b, ldb, ldbCols, &
      c, ldc, ldcCols) result(success)

!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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
!
! Author: A. Marek, MPCDF

      implicit none

      integer, parameter :: rk = C_FLOAT
      integer, parameter :: rck = C_FLOAT
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj

      character*1                   :: uplo_a, uplo_c

      integer(kind=ik), intent(in)  :: ldb, ldbCols, ldc, ldcCols
      integer(kind=ik)              :: na, ncb
      real(kind=rck)       :: a(obj%local_nrows,*), b(ldb,*), c(ldc,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
      integer(kind=MPI_KIND)        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=MPI_KIND)        :: mpierr, myidMPI
      integer(kind=ik)              :: l_cols, l_rows, l_rows_np
      integer(kind=ik)              :: np, n, nb, nblk_mult, lrs, lre, lcs, lce
      integer(kind=ik)              :: gcol_min, gcol, goff
      integer(kind=ik)              :: nstor, nr_done, noff, np_bc, n_aux_bc, nvals
      integer(kind=ik), allocatable :: lrs_save(:), lre_save(:)

      logical                       :: a_lower, a_upper, c_lower, c_upper
      real(kind=rck), pointer :: aux_mat(:,:), tmp1(:,:)
      real(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage
      character(20)                 :: gpuString
      logical                       :: success
      logical                       :: successCUDA
      logical                       :: useGPU
      integer(kind=c_int)           :: gpu, numGPU
      integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
      integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
      integer(kind=c_intptr_t)      :: aux_dev, a_dev, b_dev, tmp1_dev
      type(c_ptr)                   :: aux_host, tmp1_host
      integer(kind=c_intptr_t)      :: num
      integer(kind=c_intptr_t)      :: aux_off, a_off, b_off
      integer(kind=ik)              :: multiply_at_a
      integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
      &single&
      &_&
      &real

      success = .true.

      ! GPU settings
      call obj%get("gpu", gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif

      useGPU = (gpu == 1)

      if(useGPU) then
         gpuString = "_gpu"
      else
         gpuString = ""
      endif

      call obj%timer%start("elpa_mult_at_b_&
      &real&
      &_&
      &single&
      &"//gpuString)

      na          = obj%na
      nblk        = obj%nblk
      matrixRows  = obj%local_nrows
      matrixCols  = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_parent",mpi_comm_all,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_parent. Aborting..."
         stop
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_all, kind=MPI_KIND) ,myidMPI ,mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      myid    = int(myidMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")

      call obj%get("multiply_at_a",multiply_at_a,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for multiply_at_a. Aborting..."
         stop
      endif

      if (multiply_at_a == 0) then ! C = A^T * B
         l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
         l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b
      else ! C = A^T * A
         l_rows = local_index(ncb, my_prow, np_rows, nblk, -1) ! Local rows of a
         l_cols = local_index(na,  my_pcol, np_cols, nblk, -1) ! Local cols of a
      endif

      ! Block factor for matrix multiplications, must be a multiple of nblk

      if (na/np_rows <= 256) then
         nblk_mult = (31/nblk+1)*nblk
      else
         nblk_mult = (63/nblk+1)*nblk
      endif

      if (useGPU) then
         call obj%timer%start("check_for_gpu")
         if (check_for_gpu(myid,numGPU)) then
            ! set the neccessary parameters
            cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
            cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
            cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
            cudaHostRegisterPortable = cuda_hostRegisterPortable()
            cudaHostRegisterMapped   = cuda_hostRegisterMapped()
         else
            print *,"GPUs are requested but not detected! Aborting..."
            success = .false.
            return
         endif
         call obj%timer%stop("check_for_gpu")

         if (multiply_at_a == 0) then ! C = A^T * B
            ! copy b to b_dev
            num = ldb*ldbCols*size_of_datatype
            successCUDA = cuda_malloc(b_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_b: b_dev", 212,  successCUDA)

            successCUDA = cuda_host_register(int(loc(b),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_b: b", 217,  successCUDA)

            successCUDA = cuda_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_b: b to b_dev", 221,  successCUDA)
         else ! C = A^T * A
            ! copy a to a_dev
            num = matrixRows*matrixCols*size_of_datatype
            successCUDA = cuda_malloc(a_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_a: a_dev", 226,  successCUDA)

            successCUDA = cuda_host_register(int(loc(a),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_a: a", 231,  successCUDA)

            successCUDA = cuda_memcpy(a_dev,int(loc(a),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_a: a to a_dev", 235,  successCUDA)
         endif

         num = l_rows*nblk_mult*size_of_datatype
         successCUDA = cuda_malloc_host(aux_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: aux_host", 240,  successCUDA)

         call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

         successCUDA = cuda_malloc(aux_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: aux_dev", 245,  successCUDA)

         num = nblk_mult*l_cols*size_of_datatype
         successCUDA = cuda_malloc_host(tmp1_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: tmp1_host", 249,  successCUDA)

         call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

         successCUDA = cuda_malloc(tmp1_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: tmp1_dev", 254,  successCUDA)
      else ! useGPU
         allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa_mult_at_b: aux_mat", 257,  istat,  errorMessage)
      endif ! useGPU

      allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: aux_bc", 261,  istat,  errorMessage)

      allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lrs_save", 264,  istat,  errorMessage)

      allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lre_save", 267,  istat,  errorMessage)

      a_lower = .false.
      a_upper = .false.
      c_lower = .false.
      c_upper = .false.

      if (multiply_at_a == 0) then ! C = A^T * B
         if (uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
         if (uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
      endif
      if (uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
      if (uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

      ! Build up the result matrix by processor rows

      do np = 0, np_rows-1

         ! In this turn, procs of row np assemble the result

         l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

         nr_done = 0 ! Number of rows done
         aux_mat = 0
         nstor = 0   ! Number of columns stored in aux_mat

         ! Loop over the blocks on row np

         do nb=0,(l_rows_np-1)/nblk

            goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

            ! Get the processor column which owns this block (A is transposed, so we need the column)
            ! and the offset in blocks within this column.
            ! The corresponding block column in A is then broadcast to all for multiplication with B

            np_bc = MOD(goff,np_cols)
            noff = goff/np_cols
            n_aux_bc = 0

            ! Gather up the complete block column of A on the owner

            do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

               gcol = goff*nblk + n ! global column corresponding to n
               if (nstor==0 .and. n==1) gcol_min = gcol

               lrs = 1       ! 1st local row number for broadcast
               lre = l_rows  ! last local row number for broadcast
               if (a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  if (my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
                  n_aux_bc = n_aux_bc + nvals
               endif

               lrs_save(n) = lrs
               lre_save(n) = lre

            enddo

            ! Broadcast block column
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(aux_bc, int(n_aux_bc,kind=MPI_KIND),    &
               MPI_REAL4,  &
               int(np_bc,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            ! Insert what we got in aux_mat

            n_aux_bc = 0
            do n = 1, min(l_rows_np-nb*nblk,nblk)
               nstor = nstor+1
               lrs = lrs_save(n)
               lre = lre_save(n)
               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
                  n_aux_bc = n_aux_bc + nvals
               endif
            enddo

            ! If we got nblk_mult columns in aux_mat or this is the last block
            ! do the matrix multiplication

            if (nstor==nblk_mult .or. nb*nblk+nblk >= l_rows_np) then

               lrs = 1       ! 1st local row number for multiply
               lre = l_rows  ! last local row number for multiply
               if (a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               lcs = 1       ! 1st local col number for multiply
               lce = l_cols  ! last local col number for multiply
               if (c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
               if (c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

               if (lcs<=lce) then
                  allocate(tmp1(nstor,lcs:lce), tmp2(nstor,lcs:lce), stat=istat, errmsg=errorMessage)
                  call check_alloc("elpa_mult_at_b_&
                  &real ", "tmp1", istat, errorMessage)

                  if (lrs<=lre) then
                     if (useGPU) then
                        num = l_rows*nblk_mult*size_of_datatype
                        successCUDA = cuda_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                           num, cudaMemcpyHostToDevice)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: aux_mat to aux_dev", 377,  successCUDA)

                        call obj%timer%start("cublas")

                        aux_off = (lrs-1)*size_of_datatype
                        if (multiply_at_a == 0) then ! C = A^T * B
                           b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

                           call cublas_SGEMM('T', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                              tmp1_dev, nstor)
                        else ! C = A^T * A
                           a_off = ((lcs-1)*matrixRows+lrs-1)*size_of_datatype

                           call cublas_SGEMM('T', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, a_dev+a_off, matrixRows, ZERO, &
                              tmp1_dev, nstor)
                        endif

                        call obj%timer%stop("cublas")

                        num = nstor*(lce-lcs+1)*size_of_datatype
                        successCUDA = cuda_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                           tmp1_dev, num, cudaMemcpyDeviceToHost)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: tmp1_dev to tmp1", 401,  successCUDA)
                     else ! useGPU
                        call obj%timer%start("blas")

                        if (multiply_at_a == 0) then ! C = A^T * B
                           call SGEMM('T', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        else
                           call SGEMM('T', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              a(lrs,lcs), int(matrixRows,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        endif

                        call obj%timer%stop("blas")
                     endif ! useGPU
                  else
                     tmp1 = 0
                  endif

                  ! Sum up the results and send to processor row np
                  call obj%timer%start("mpi_communication")
                  call mpi_reduce(tmp1, tmp2, int(nstor*(lce-lcs+1),kind=MPI_KIND),  MPI_REAL4, &
                     MPI_SUM, int(np,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
                  call obj%timer%stop("mpi_communication")
                  ! Put the result into C
                  if (my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

                  deallocate(tmp1,tmp2, stat=istat, errmsg=errorMessage)
                  call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 442,  istat,  errorMessage)
               endif

               nr_done = nr_done+nstor
               nstor=0
               aux_mat(:,:)=0
            endif
         enddo
      enddo

      if (useGPU) then
         if (multiply_at_a == 0) then ! C = A^T * B
            successCUDA = cuda_free(b_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: b_dev", 455,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(b),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: b", 458,  successCUDA)
         else ! C = A^T * A
            successCUDA = cuda_free(a_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: a_dev", 461,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(a),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: a", 464,  successCUDA)
         endif

         nullify(aux_mat)
         nullify(tmp1)

         successCUDA = cuda_free_host(aux_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: aux_host", 471,  successCUDA)

         successCUDA = cuda_free(aux_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: aux_dev", 474,  successCUDA)

         successCUDA = cuda_free_host(tmp1_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: tmp1_host", 477,  successCUDA)

         successCUDA = cuda_free(tmp1_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: tmp1_dev", 480,  successCUDA)
      else ! useGPU
         deallocate(aux_mat, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa_mult_at_b: aux_mat", 483,  istat,  errorMessage)
      endif ! useGPU

      deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 487,  istat,  errorMessage)

      call obj%timer%stop("elpa_mult_at_b_&
      &real&
      &_&
      &single&
      &"//gpuString)

   end function elpa_mult_at_b_real_single_impl

!> \brief  elpa_mult_ah_b_complex_double_impl: Performs C : = A**H * B
!>         where   A is a square matrix (obj%na,obj%na) which is optionally upper or lower triangular
!>                 B is a (obj%na,ncb) matrix
!>                 C is a (obj%na,ncb) matrix where optionally only the upper or lower
!>                   triangle may be computed
!> \details
!>
!> \param  uplo_a               'U' if A is upper triangular
!>                              'L' if A is lower triangular
!>                              anything else if A is a full matrix
!>                              Please note: This pertains to the original A (as set in the calling program)
!>                                           whereas the transpose of A is used for calculations
!>                              If uplo_a is 'U' or 'L', the other triangle is not used at all,
!>                              i.e. it may contain arbitrary numbers
!> \param uplo_c                'U' if only the upper diagonal part of C is needed
!>                              'L' if only the upper diagonal part of C is needed
!>                              anything else if the full matrix C is needed
!>                              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!>                                            written to a certain extent, i.e. one shouldn't rely on the content there!
!> \param na                    Number of rows/columns of A, number of rows of B and C
!> \param ncb                   Number of columns  of B and C
!> \param a                     matrix a
!> \param obj%local_ncols       leading dimension of matrix a, set with class method obj%set("local_nrows",value)
!> \param ldaCols               columns of matrix a
!> \param b                     matrix b
!> \param ldb                   leading dimension of matrix b
!> \param ldbCols               columns of matrix b
!> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \param c                     matrix c
!> \param ldc                   leading dimension of matrix c
!> \result success

   function elpa_mult_ah_b_complex_double_impl(obj, uplo_a, uplo_c, ncb, a, b, ldb, ldbCols, &
      c, ldc, ldcCols) result(success)
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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
!
! Author: A. Marek, MPCDF

      implicit none

      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: ck = C_DOUBLE_COMPLEX
      integer, parameter :: rck = C_DOUBLE_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout) :: obj

      character*1                   :: uplo_a, uplo_c

      integer(kind=ik), intent(in)  :: ldb, ldbCols, ldc, ldcCols
      integer(kind=ik)              :: na, ncb
      complex(kind=rck)       :: a(obj%local_nrows,*), b(ldb,*), c(ldc,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
      integer(kind=MPI_KIND)        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=MPI_KIND)        :: mpierr, myidMPI
      integer(kind=ik)              :: l_cols, l_rows, l_rows_np
      integer(kind=ik)              :: np, n, nb, nblk_mult, lrs, lre, lcs, lce
      integer(kind=ik)              :: gcol_min, gcol, goff
      integer(kind=ik)              :: nstor, nr_done, noff, np_bc, n_aux_bc, nvals
      integer(kind=ik), allocatable :: lrs_save(:), lre_save(:)

      logical                       :: a_lower, a_upper, c_lower, c_upper
      complex(kind=rck), pointer :: aux_mat(:,:), tmp1(:,:)
      complex(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage
      character(20)                 :: gpuString
      logical                       :: success
      logical                       :: successCUDA
      logical                       :: useGPU
      integer(kind=c_int)           :: gpu, numGPU
      integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
      integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
      integer(kind=c_intptr_t)      :: aux_dev, a_dev, b_dev, tmp1_dev
      type(c_ptr)                   :: aux_host, tmp1_host
      integer(kind=c_intptr_t)      :: num
      integer(kind=c_intptr_t)      :: aux_off, a_off, b_off
      integer(kind=ik)              :: multiply_at_a
      integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
      &double&
      &_&
      &complex

      success = .true.

      ! GPU settings
      call obj%get("gpu", gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif

      useGPU = (gpu == 1)

      if(useGPU) then
         gpuString = "_gpu"
      else
         gpuString = ""
      endif

      call obj%timer%start("elpa_mult_at_b_&
      &complex&
      &_&
      &double&
      &"//gpuString)

      na          = obj%na
      nblk        = obj%nblk
      matrixRows  = obj%local_nrows
      matrixCols  = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_parent",mpi_comm_all,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_parent. Aborting..."
         stop
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_all, kind=MPI_KIND) ,myidMPI ,mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      myid    = int(myidMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")

      call obj%get("multiply_at_a",multiply_at_a,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for multiply_at_a. Aborting..."
         stop
      endif

      if (multiply_at_a == 0) then ! C = A^T * B
         l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
         l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b
      else ! C = A^T * A
         l_rows = local_index(ncb, my_prow, np_rows, nblk, -1) ! Local rows of a
         l_cols = local_index(na,  my_pcol, np_cols, nblk, -1) ! Local cols of a
      endif

      ! Block factor for matrix multiplications, must be a multiple of nblk

      if (na/np_rows <= 256) then
         nblk_mult = (31/nblk+1)*nblk
      else
         nblk_mult = (63/nblk+1)*nblk
      endif

      if (useGPU) then
         call obj%timer%start("check_for_gpu")
         if (check_for_gpu(myid,numGPU)) then
            ! set the neccessary parameters
            cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
            cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
            cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
            cudaHostRegisterPortable = cuda_hostRegisterPortable()
            cudaHostRegisterMapped   = cuda_hostRegisterMapped()
         else
            print *,"GPUs are requested but not detected! Aborting..."
            success = .false.
            return
         endif
         call obj%timer%stop("check_for_gpu")

         if (multiply_at_a == 0) then ! C = A^T * B
            ! copy b to b_dev
            num = ldb*ldbCols*size_of_datatype
            successCUDA = cuda_malloc(b_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_b: b_dev", 212,  successCUDA)

            successCUDA = cuda_host_register(int(loc(b),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_b: b", 217,  successCUDA)

            successCUDA = cuda_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_b: b to b_dev", 221,  successCUDA)
         else ! C = A^T * A
            ! copy a to a_dev
            num = matrixRows*matrixCols*size_of_datatype
            successCUDA = cuda_malloc(a_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_a: a_dev", 226,  successCUDA)

            successCUDA = cuda_host_register(int(loc(a),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_a: a", 231,  successCUDA)

            successCUDA = cuda_memcpy(a_dev,int(loc(a),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_a: a to a_dev", 235,  successCUDA)
         endif

         num = l_rows*nblk_mult*size_of_datatype
         successCUDA = cuda_malloc_host(aux_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: aux_host", 240,  successCUDA)

         call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

         successCUDA = cuda_malloc(aux_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: aux_dev", 245,  successCUDA)

         num = nblk_mult*l_cols*size_of_datatype
         successCUDA = cuda_malloc_host(tmp1_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: tmp1_host", 249,  successCUDA)

         call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

         successCUDA = cuda_malloc(tmp1_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: tmp1_dev", 254,  successCUDA)
      else ! useGPU
         allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa_mult_at_b: aux_mat", 257,  istat,  errorMessage)
      endif ! useGPU

      allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: aux_bc", 261,  istat,  errorMessage)

      allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lrs_save", 264,  istat,  errorMessage)

      allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lre_save", 267,  istat,  errorMessage)

      a_lower = .false.
      a_upper = .false.
      c_lower = .false.
      c_upper = .false.

      if (multiply_at_a == 0) then ! C = A^T * B
         if (uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
         if (uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
      endif
      if (uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
      if (uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

      ! Build up the result matrix by processor rows

      do np = 0, np_rows-1

         ! In this turn, procs of row np assemble the result

         l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

         nr_done = 0 ! Number of rows done
         aux_mat = 0
         nstor = 0   ! Number of columns stored in aux_mat

         ! Loop over the blocks on row np

         do nb=0,(l_rows_np-1)/nblk

            goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

            ! Get the processor column which owns this block (A is transposed, so we need the column)
            ! and the offset in blocks within this column.
            ! The corresponding block column in A is then broadcast to all for multiplication with B

            np_bc = MOD(goff,np_cols)
            noff = goff/np_cols
            n_aux_bc = 0

            ! Gather up the complete block column of A on the owner

            do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

               gcol = goff*nblk + n ! global column corresponding to n
               if (nstor==0 .and. n==1) gcol_min = gcol

               lrs = 1       ! 1st local row number for broadcast
               lre = l_rows  ! last local row number for broadcast
               if (a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  if (my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
                  n_aux_bc = n_aux_bc + nvals
               endif

               lrs_save(n) = lrs
               lre_save(n) = lre

            enddo

            ! Broadcast block column
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(aux_bc, int(n_aux_bc,kind=MPI_KIND),    &
               MPI_DOUBLE_COMPLEX,  &
               int(np_bc,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            ! Insert what we got in aux_mat

            n_aux_bc = 0
            do n = 1, min(l_rows_np-nb*nblk,nblk)
               nstor = nstor+1
               lrs = lrs_save(n)
               lre = lre_save(n)
               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
                  n_aux_bc = n_aux_bc + nvals
               endif
            enddo

            ! If we got nblk_mult columns in aux_mat or this is the last block
            ! do the matrix multiplication

            if (nstor==nblk_mult .or. nb*nblk+nblk >= l_rows_np) then

               lrs = 1       ! 1st local row number for multiply
               lre = l_rows  ! last local row number for multiply
               if (a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               lcs = 1       ! 1st local col number for multiply
               lce = l_cols  ! last local col number for multiply
               if (c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
               if (c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

               if (lcs<=lce) then
                  allocate(tmp1(nstor,lcs:lce), tmp2(nstor,lcs:lce), stat=istat, errmsg=errorMessage)
                  call check_alloc("elpa_mult_at_b_&
                  &complex ", "tmp1", istat, errorMessage)

                  if (lrs<=lre) then
                     if (useGPU) then
                        num = l_rows*nblk_mult*size_of_datatype
                        successCUDA = cuda_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                           num, cudaMemcpyHostToDevice)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: aux_mat to aux_dev", 377,  successCUDA)

                        call obj%timer%start("cublas")

                        aux_off = (lrs-1)*size_of_datatype
                        if (multiply_at_a == 0) then ! C = A^T * B
                           b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

                           call cublas_ZGEMM('C', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                              tmp1_dev, nstor)
                        else ! C = A^T * A
                           a_off = ((lcs-1)*matrixRows+lrs-1)*size_of_datatype

                           call cublas_ZGEMM('C', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, a_dev+a_off, matrixRows, ZERO, &
                              tmp1_dev, nstor)
                        endif

                        call obj%timer%stop("cublas")

                        num = nstor*(lce-lcs+1)*size_of_datatype
                        successCUDA = cuda_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                           tmp1_dev, num, cudaMemcpyDeviceToHost)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: tmp1_dev to tmp1", 401,  successCUDA)
                     else ! useGPU
                        call obj%timer%start("blas")

                        if (multiply_at_a == 0) then ! C = A^T * B
                           call ZGEMM('C', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        else
                           call ZGEMM('C', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              a(lrs,lcs), int(matrixRows,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        endif

                        call obj%timer%stop("blas")
                     endif ! useGPU
                  else
                     tmp1 = 0
                  endif

                  ! Sum up the results and send to processor row np
                  call obj%timer%start("mpi_communication")
                  call mpi_reduce(tmp1, tmp2, int(nstor*(lce-lcs+1),kind=MPI_KIND),  MPI_DOUBLE_COMPLEX, &
                     MPI_SUM, int(np,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
                  call obj%timer%stop("mpi_communication")
                  ! Put the result into C
                  if (my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

                  deallocate(tmp1,tmp2, stat=istat, errmsg=errorMessage)
                  call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 442,  istat,  errorMessage)
               endif

               nr_done = nr_done+nstor
               nstor=0
               aux_mat(:,:)=0
            endif
         enddo
      enddo

      if (useGPU) then
         if (multiply_at_a == 0) then ! C = A^T * B
            successCUDA = cuda_free(b_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: b_dev", 455,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(b),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: b", 458,  successCUDA)
         else ! C = A^T * A
            successCUDA = cuda_free(a_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: a_dev", 461,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(a),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: a", 464,  successCUDA)
         endif

         nullify(aux_mat)
         nullify(tmp1)

         successCUDA = cuda_free_host(aux_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: aux_host", 471,  successCUDA)

         successCUDA = cuda_free(aux_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: aux_dev", 474,  successCUDA)

         successCUDA = cuda_free_host(tmp1_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: tmp1_host", 477,  successCUDA)

         successCUDA = cuda_free(tmp1_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: tmp1_dev", 480,  successCUDA)
      else ! useGPU
         deallocate(aux_mat, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa_mult_at_b: aux_mat", 483,  istat,  errorMessage)
      endif ! useGPU

      deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 487,  istat,  errorMessage)

      call obj%timer%stop("elpa_mult_at_b_&
      &complex&
      &_&
      &double&
      &"//gpuString)

   end function elpa_mult_ah_b_complex_double_impl

!> \brief  elpa_mult_ah_b_complex_single_impl: Performs C : = A**H * B
!>         where   A is a square matrix (obj%na,obj%na) which is optionally upper or lower triangular
!>                 B is a (obj%na,ncb) matrix
!>                 C is a (obj%na,ncb) matrix where optionally only the upper or lower
!>                   triangle may be computed
!> \details
!>
!> \param  uplo_a               'U' if A is upper triangular
!>                              'L' if A is lower triangular
!>                              anything else if A is a full matrix
!>                              Please note: This pertains to the original A (as set in the calling program)
!>                                           whereas the transpose of A is used for calculations
!>                              If uplo_a is 'U' or 'L', the other triangle is not used at all,
!>                              i.e. it may contain arbitrary numbers
!> \param uplo_c                'U' if only the upper diagonal part of C is needed
!>                              'L' if only the upper diagonal part of C is needed
!>                              anything else if the full matrix C is needed
!>                              Please note: Even when uplo_c is 'U' or 'L', the other triangle may be
!>                                            written to a certain extent, i.e. one shouldn't rely on the content there!
!> \param na                    Number of rows/columns of A, number of rows of B and C
!> \param ncb                   Number of columns  of B and C
!> \param a                     matrix a
!> \param lda                   leading dimension of matrix a
!> \param ldaCols               columns of matrix a
!> \param b                     matrix b
!> \param ldb                   leading dimension of matrix b
!> \param ldbCols               columns of matrix b
!> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \param c                     matrix c
!> \param ldc                   leading dimension of matrix c
!> \result success

   function elpa_mult_ah_b_complex_single_impl(obj, uplo_a, uplo_c, ncb, a, b, ldb, ldbCols, &
      c, ldc, ldcCols) result(success)

!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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
!
! Author: A. Marek, MPCDF

      implicit none

      integer, parameter :: rk = C_FLOAT
      integer, parameter :: ck = C_FLOAT_COMPLEX
      integer, parameter :: rck = C_FLOAT_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout) :: obj

      character*1                   :: uplo_a, uplo_c

      integer(kind=ik), intent(in)  :: ldb, ldbCols, ldc, ldcCols
      integer(kind=ik)              :: na, ncb
      complex(kind=rck)       :: a(obj%local_nrows,*), b(ldb,*), c(ldc,*)
      integer(kind=ik)              :: my_prow, my_pcol, np_rows, np_cols, myid
      integer(kind=MPI_KIND)        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=MPI_KIND)        :: mpierr, myidMPI
      integer(kind=ik)              :: l_cols, l_rows, l_rows_np
      integer(kind=ik)              :: np, n, nb, nblk_mult, lrs, lre, lcs, lce
      integer(kind=ik)              :: gcol_min, gcol, goff
      integer(kind=ik)              :: nstor, nr_done, noff, np_bc, n_aux_bc, nvals
      integer(kind=ik), allocatable :: lrs_save(:), lre_save(:)

      logical                       :: a_lower, a_upper, c_lower, c_upper
      complex(kind=rck), pointer :: aux_mat(:,:), tmp1(:,:)
      complex(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage
      character(20)                 :: gpuString
      logical                       :: success
      logical                       :: successCUDA
      logical                       :: useGPU
      integer(kind=c_int)           :: gpu, numGPU
      integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
      integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
      integer(kind=c_intptr_t)      :: aux_dev, a_dev, b_dev, tmp1_dev
      type(c_ptr)                   :: aux_host, tmp1_host
      integer(kind=c_intptr_t)      :: num
      integer(kind=c_intptr_t)      :: aux_off, a_off, b_off
      integer(kind=ik)              :: multiply_at_a
      integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
      &single&
      &_&
      &complex

      success = .true.

      ! GPU settings
      call obj%get("gpu", gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif

      useGPU = (gpu == 1)

      if(useGPU) then
         gpuString = "_gpu"
      else
         gpuString = ""
      endif

      call obj%timer%start("elpa_mult_at_b_&
      &complex&
      &_&
      &single&
      &"//gpuString)

      na          = obj%na
      nblk        = obj%nblk
      matrixRows  = obj%local_nrows
      matrixCols  = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_parent",mpi_comm_all,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_parent. Aborting..."
         stop
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)
      call mpi_comm_rank(int(mpi_comm_all, kind=MPI_KIND) ,myidMPI ,mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)
      myid    = int(myidMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")

      call obj%get("multiply_at_a",multiply_at_a,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for multiply_at_a. Aborting..."
         stop
      endif

      if (multiply_at_a == 0) then ! C = A^T * B
         l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
         l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b
      else ! C = A^T * A
         l_rows = local_index(ncb, my_prow, np_rows, nblk, -1) ! Local rows of a
         l_cols = local_index(na,  my_pcol, np_cols, nblk, -1) ! Local cols of a
      endif

      ! Block factor for matrix multiplications, must be a multiple of nblk

      if (na/np_rows <= 256) then
         nblk_mult = (31/nblk+1)*nblk
      else
         nblk_mult = (63/nblk+1)*nblk
      endif

      if (useGPU) then
         call obj%timer%start("check_for_gpu")
         if (check_for_gpu(myid,numGPU)) then
            ! set the neccessary parameters
            cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
            cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
            cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
            cudaHostRegisterPortable = cuda_hostRegisterPortable()
            cudaHostRegisterMapped   = cuda_hostRegisterMapped()
         else
            print *,"GPUs are requested but not detected! Aborting..."
            success = .false.
            return
         endif
         call obj%timer%stop("check_for_gpu")

         if (multiply_at_a == 0) then ! C = A^T * B
            ! copy b to b_dev
            num = ldb*ldbCols*size_of_datatype
            successCUDA = cuda_malloc(b_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_b: b_dev", 212,  successCUDA)

            successCUDA = cuda_host_register(int(loc(b),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_b: b", 217,  successCUDA)

            successCUDA = cuda_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_b: b to b_dev", 221,  successCUDA)
         else ! C = A^T * A
            ! copy a to a_dev
            num = matrixRows*matrixCols*size_of_datatype
            successCUDA = cuda_malloc(a_dev,num)
            call check_alloc_CUDA_f("elpa_mult_at_a: a_dev", 226,  successCUDA)

            successCUDA = cuda_host_register(int(loc(a),kind=c_intptr_t),num,&
               cudaHostRegisterDefault)

            call check_host_register_CUDA_f("elpa_mult_at_a: a", 231,  successCUDA)

            successCUDA = cuda_memcpy(a_dev,int(loc(a),kind=c_intptr_t),num,&
               cudaMemcpyHostToDevice)
            call check_memcpy_CUDA_f("elpa_mult_at_a: a to a_dev", 235,  successCUDA)
         endif

         num = l_rows*nblk_mult*size_of_datatype
         successCUDA = cuda_malloc_host(aux_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: aux_host", 240,  successCUDA)

         call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

         successCUDA = cuda_malloc(aux_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: aux_dev", 245,  successCUDA)

         num = nblk_mult*l_cols*size_of_datatype
         successCUDA = cuda_malloc_host(tmp1_host,num)
         call check_host_alloc_CUDA_f("elpa_mult_at_b: tmp1_host", 249,  successCUDA)

         call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

         successCUDA = cuda_malloc(tmp1_dev,num)
         call check_alloc_CUDA_f("elpa_mult_at_b: tmp1_dev", 254,  successCUDA)
      else ! useGPU
         allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa_mult_at_b: aux_mat", 257,  istat,  errorMessage)
      endif ! useGPU

      allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: aux_bc", 261,  istat,  errorMessage)

      allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lrs_save", 264,  istat,  errorMessage)

      allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa_mult_at_b: lre_save", 267,  istat,  errorMessage)

      a_lower = .false.
      a_upper = .false.
      c_lower = .false.
      c_upper = .false.

      if (multiply_at_a == 0) then ! C = A^T * B
         if (uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
         if (uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
      endif
      if (uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
      if (uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

      ! Build up the result matrix by processor rows

      do np = 0, np_rows-1

         ! In this turn, procs of row np assemble the result

         l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

         nr_done = 0 ! Number of rows done
         aux_mat = 0
         nstor = 0   ! Number of columns stored in aux_mat

         ! Loop over the blocks on row np

         do nb=0,(l_rows_np-1)/nblk

            goff  = nb*np_rows + np ! Global offset in blocks corresponding to nb

            ! Get the processor column which owns this block (A is transposed, so we need the column)
            ! and the offset in blocks within this column.
            ! The corresponding block column in A is then broadcast to all for multiplication with B

            np_bc = MOD(goff,np_cols)
            noff = goff/np_cols
            n_aux_bc = 0

            ! Gather up the complete block column of A on the owner

            do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

               gcol = goff*nblk + n ! global column corresponding to n
               if (nstor==0 .and. n==1) gcol_min = gcol

               lrs = 1       ! 1st local row number for broadcast
               lre = l_rows  ! last local row number for broadcast
               if (a_lower) lrs = local_index(gcol, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  if (my_pcol == np_bc) aux_bc(n_aux_bc+1:n_aux_bc+nvals) = a(lrs:lre,noff*nblk+n)
                  n_aux_bc = n_aux_bc + nvals
               endif

               lrs_save(n) = lrs
               lre_save(n) = lre

            enddo

            ! Broadcast block column
            call obj%timer%start("mpi_communication")
            call MPI_Bcast(aux_bc, int(n_aux_bc,kind=MPI_KIND),    &
               MPI_COMPLEX,  &
               int(np_bc,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
            ! Insert what we got in aux_mat

            n_aux_bc = 0
            do n = 1, min(l_rows_np-nb*nblk,nblk)
               nstor = nstor+1
               lrs = lrs_save(n)
               lre = lre_save(n)
               if (lrs<=lre) then
                  nvals = lre-lrs+1
                  aux_mat(lrs:lre,nstor) = aux_bc(n_aux_bc+1:n_aux_bc+nvals)
                  n_aux_bc = n_aux_bc + nvals
               endif
            enddo

            ! If we got nblk_mult columns in aux_mat or this is the last block
            ! do the matrix multiplication

            if (nstor==nblk_mult .or. nb*nblk+nblk >= l_rows_np) then

               lrs = 1       ! 1st local row number for multiply
               lre = l_rows  ! last local row number for multiply
               if (a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
               if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

               lcs = 1       ! 1st local col number for multiply
               lce = l_cols  ! last local col number for multiply
               if (c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
               if (c_lower) lce = MIN(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

               if (lcs<=lce) then
                  allocate(tmp1(nstor,lcs:lce), tmp2(nstor,lcs:lce), stat=istat, errmsg=errorMessage)
                  call check_alloc("elpa_mult_at_b_&
                  &complex ", "tmp1", istat, errorMessage)

                  if (lrs<=lre) then
                     if (useGPU) then
                        num = l_rows*nblk_mult*size_of_datatype
                        successCUDA = cuda_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                           num, cudaMemcpyHostToDevice)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: aux_mat to aux_dev", 377,  successCUDA)

                        call obj%timer%start("cublas")

                        aux_off = (lrs-1)*size_of_datatype
                        if (multiply_at_a == 0) then ! C = A^T * B
                           b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

                           call cublas_CGEMM('C', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                              tmp1_dev, nstor)
                        else ! C = A^T * A
                           a_off = ((lcs-1)*matrixRows+lrs-1)*size_of_datatype

                           call cublas_CGEMM('C', 'N', nstor, lce-lcs+1, &
                              lre-lrs+1, ONE, aux_dev+aux_off, l_rows, a_dev+a_off, matrixRows, ZERO, &
                              tmp1_dev, nstor)
                        endif

                        call obj%timer%stop("cublas")

                        num = nstor*(lce-lcs+1)*size_of_datatype
                        successCUDA = cuda_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                           tmp1_dev, num, cudaMemcpyDeviceToHost)
                        call check_memcpy_CUDA_f("elpa_mult_at_b: tmp1_dev to tmp1", 401,  successCUDA)
                     else ! useGPU
                        call obj%timer%start("blas")

                        if (multiply_at_a == 0) then ! C = A^T * B
                           call CGEMM('C', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        else
                           call CGEMM('C', 'N', int(nstor,kind=BLAS_KIND), &
                              int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                              ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                              a(lrs,lcs), int(matrixRows,kind=BLAS_KIND), ZERO, tmp1, &
                              int(nstor,kind=BLAS_KIND))
                        endif

                        call obj%timer%stop("blas")
                     endif ! useGPU
                  else
                     tmp1 = 0
                  endif

                  ! Sum up the results and send to processor row np
                  call obj%timer%start("mpi_communication")
                  call mpi_reduce(tmp1, tmp2, int(nstor*(lce-lcs+1),kind=MPI_KIND),  MPI_COMPLEX, &
                     MPI_SUM, int(np,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), mpierr)
                  call obj%timer%stop("mpi_communication")
                  ! Put the result into C
                  if (my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

                  deallocate(tmp1,tmp2, stat=istat, errmsg=errorMessage)
                  call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 442,  istat,  errorMessage)
               endif

               nr_done = nr_done+nstor
               nstor=0
               aux_mat(:,:)=0
            endif
         enddo
      enddo

      if (useGPU) then
         if (multiply_at_a == 0) then ! C = A^T * B
            successCUDA = cuda_free(b_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: b_dev", 455,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(b),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: b", 458,  successCUDA)
         else ! C = A^T * A
            successCUDA = cuda_free(a_dev)
            call check_dealloc_CUDA_f("elpa_mult_a_b: a_dev", 461,  successCUDA)

            successCUDA = cuda_host_unregister(int(loc(a),kind=c_intptr_t))
            call check_host_unregister_CUDA_f("elpa_mult_a_b: a", 464,  successCUDA)
         endif

         nullify(aux_mat)
         nullify(tmp1)

         successCUDA = cuda_free_host(aux_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: aux_host", 471,  successCUDA)

         successCUDA = cuda_free(aux_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: aux_dev", 474,  successCUDA)

         successCUDA = cuda_free_host(tmp1_host)
         call check_host_dealloc_CUDA_f("elpa_mult_a_b: tmp1_host", 477,  successCUDA)

         successCUDA = cuda_free(tmp1_dev)
         call check_dealloc_CUDA_f("elpa_mult_a_b: tmp1_dev", 480,  successCUDA)
      else ! useGPU
         deallocate(aux_mat, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa_mult_at_b: aux_mat", 483,  istat,  errorMessage)
      endif ! useGPU

      deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 487,  istat,  errorMessage)

      call obj%timer%stop("elpa_mult_at_b_&
      &complex&
      &_&
      &single&
      &"//gpuString)

   end function elpa_mult_ah_b_complex_single_impl

!> \brief  elpa_solve_tridi_double_impl: Solve tridiagonal eigensystem for a double-precision matrix with divide and conquer method
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%nev           number of eigenvalues/vectors to be computed
!> \param     - obj%local_nrows   Leading dimension of q
!> \param     - obj%local_ncols   local columns of matrix q
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param d                       array d(na) on input diagonal elements of tridiagonal matrix, on
!>                                output the eigenvalues in ascending order
!> \param e                       array e(na) on input subdiagonal elements of matrix, on exit destroyed
!> \param q                       on exit : matrix q(ldq,matrixCols) contains the eigenvectors
!> \result succes                 logical, reports success or failure
   function elpa_solve_tridi_double_impl(obj, d, e, q) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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
!
! Author: A. Marek, MPCDF

      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)         :: na, nev, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      real(kind=rk8) :: d(obj%na), e(obj%na)
      real(kind=rk8) :: q(obj%local_nrows,*)

      logical                  :: wantDebug
      logical                  :: success

      integer                  :: debug, error
      integer                  :: nrThreads

      call obj%timer%start("elpa_solve_tridi_public_&
      &real&
      &_&
      &double&
      &")
      na         = obj%na
      nev        = obj%nev
      nblk       = obj%nblk
      matrixRows = obj%local_nrows
      matrixCols = obj%local_ncols

      nrThreads=1

      call obj%get("mpi_comm_rows", mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols", mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug",debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for debug. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .false.
      endif
      success = .false.

      call solve_tridi_&
      &double&
      &_private_impl(obj, na, nev, d, e, q, matrixRows, nblk, matrixCols, &
         mpi_comm_rows, mpi_comm_cols,.false., wantDebug, success, &
         nrThreads)

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_solve_tridi_public_&
      &real&
      &_&
      &double&
      &")

   end function

!> \brief  elpa_solve_tridi_single_impl: Solve tridiagonal eigensystem for a single-precision matrix with divide and conquer method
!> \details
!> \param  obj                    elpa_t object contains:
!> \param     - obj%na            Order of matrix
!> \param     - obj%nev           number of eigenvalues/vectors to be computed
!> \param     - obj%local_nrows   Leading dimension of q
!> \param     - obj%local_ncols   local columns of matrix q
!> \param     - obj%nblk          blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows MPI communicator for rows
!> \param     - obj%mpi_comm_cols MPI communicator for columns
!> \param     - obj%wantDebug     logical, more debug information on failure
!> \param d                       array d(na) on input diagonal elements of tridiagonal matrix, on
!>                                output the eigenvalues in ascending order
!> \param e                       array e(na) on input subdiagonal elements of matrix, on exit destroyed
!> \param q                       on exit : matrix q(ldq,matrixCols) contains the eigenvectors
!> \result succes                 logical, reports success or failure
   function elpa_solve_tridi_single_impl(obj, d, e, q) result(success)

!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), formerly known as
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
!
! Author: A. Marek, MPCDF

      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)         :: na, nev, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      real(kind=rk4) :: d(obj%na), e(obj%na)
      real(kind=rk4) :: q(obj%local_nrows,*)

      logical                  :: wantDebug
      logical                  :: success

      integer                  :: debug, error
      integer                  :: nrThreads

      call obj%timer%start("elpa_solve_tridi_public_&
      &real&
      &_&
      &single&
      &")
      na         = obj%na
      nev        = obj%nev
      nblk       = obj%nblk
      matrixRows = obj%local_nrows
      matrixCols = obj%local_ncols

      nrThreads=1

      call obj%get("mpi_comm_rows", mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_rows. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols", mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for mpi_comm_cols. Aborting..."
         stop
      endif

      call obj%get("debug",debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for debug. Aborting..."
         stop
      endif
      if (debug == 1) then
         wantDebug = .true.
      else
         wantDebug = .false.
      endif
      success = .false.

      call solve_tridi_&
      &single&
      &_private_impl(obj, na, nev, d, e, q, matrixRows, nblk, matrixCols, &
         mpi_comm_rows, mpi_comm_cols,.false., wantDebug, success, &
         nrThreads)

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_solve_tridi_public_&
      &real&
      &_&
      &single&
      &")

   end function

end module elpa1_auxiliary_impl

