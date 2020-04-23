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

!> \mainpage
!> Eigenvalue SoLvers for Petaflop-Applications (ELPA)
!> \par
!> http://elpa.mpcdf.mpg.de
!>
!> \par
!>    The ELPA library was originally created by the ELPA consortium,
!>    consisting of the following organizations:
!>
!>    - Max Planck Computing and Data Facility (MPCDF) formerly known as
!>      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
!>    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
!>      Informatik,
!>    - Technische Universität München, Lehrstuhl für Informatik mit
!>      Schwerpunkt Wissenschaftliches Rechnen ,
!>    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
!>    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
!>      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
!>      and
!>    - IBM Deutschland GmbH
!>
!>   Some parts and enhancements of ELPA have been contributed and authored
!>   by the Intel Corporation which is not part of the ELPA consortium.
!>
!>   Contributions to the ELPA source have been authored by (in alphabetical order):
!>
!> \author T. Auckenthaler, Volker Blum, A. Heinecke, L. Huedepohl, R. Johanni, Werner Jürgens, and A. Marek

!> \brief Fortran module which provides the routines to use the one-stage ELPA solver
module elpa1_impl
   use, intrinsic :: iso_c_binding
   use elpa_utilities
   use elpa1_auxiliary_impl

   implicit none

   ! The following routines are public:
   private

   public :: elpa_solve_evp_real_1stage_double_impl    !< Driver routine for real double-precision 1-stage eigenvalue problem

   public :: elpa_solve_evp_real_1stage_single_impl    !< Driver routine for real single-precision 1-stage eigenvalue problem

   public :: elpa_solve_evp_complex_1stage_double_impl !< Driver routine for complex 1-stage eigenvalue problem
   public :: elpa_solve_evp_complex_1stage_single_impl !< Driver routine for complex 1-stage eigenvalue problem

   ! imported from elpa1_auxilliary

   public :: elpa_mult_at_b_real_double_impl       !< Multiply double-precision real matrices A**T * B

   public :: elpa_mult_ah_b_complex_double_impl    !< Multiply double-precision complex matrices A**H * B

   public :: elpa_invert_trm_real_double_impl      !< Invert double-precision real triangular matrix

   public :: elpa_invert_trm_complex_double_impl   !< Invert double-precision complex triangular matrix

   public :: elpa_cholesky_real_double_impl        !< Cholesky factorization of a double-precision real matrix

   public :: elpa_cholesky_complex_double_impl     !< Cholesky factorization of a double-precision complex matrix

   public :: elpa_solve_tridi_double_impl          !< Solve a double-precision tridiagonal eigensystem with divide and conquer method

   public :: elpa_mult_at_b_real_single_impl       !< Multiply single-precision real matrices A**T * B
   public :: elpa_invert_trm_real_single_impl      !< Invert single-precision real triangular matrix
   public :: elpa_cholesky_real_single_impl        !< Cholesky factorization of a single-precision real matrix
   public :: elpa_solve_tridi_single_impl          !< Solve a single-precision tridiagonal eigensystem with divide and conquer method

   public :: elpa_mult_ah_b_complex_single_impl    !< Multiply single-precision complex matrices A**H * B
   public :: elpa_invert_trm_complex_single_impl   !< Invert single-precision complex triangular matrix
   public :: elpa_cholesky_complex_single_impl     !< Cholesky factorization of a single-precision complex matrix

contains

!> \brief elpa_solve_evp_real_1stage_double_impl: Fortran function to solve the real double-precision eigenvalue problem with 1-stage solver
!>
!> \details
!> \param  obj                      elpa_t object contains:
!> \param     - obj%na              Order of matrix
!> \param     - obj%nev             number of eigenvalues/vectors to be computed
!>                                  The smallest nev eigenvalues/eigenvectors are calculated.
!> \param     - obj%local_nrows     Leading dimension of a
!> \param     - obj%local_ncols     local columns of matrix q
!> \param     - obj%nblk            blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows   MPI communicator for rows
!> \param     - obj%mpi_comm_cols   MPI communicator for columns
!> \param     - obj%mpi_comm_parent MPI communicator for columns
!> \param     - obj%gpu             use GPU version (1 or 0)
!>
!> \param  a(lda,matrixCols)        Distributed matrix for which eigenvalues are to be computed.
!>                                  Distribution is like in Scalapack.
!>                                  The full matrix must be set (not only one half like in scalapack).
!>                                  Destroyed on exit (upper and lower half).
!>
!>  \param ev(na)                   On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)        On output: Eigenvectors of a
!>                                  Distribution is like in Scalapack.
!>                                  Must be always dimensioned to the full size (corresponding to (na,na))
!>                                  even if only a part of the eigenvalues is needed.
!>
!>
!>  \result                       success

   function elpa_solve_evp_&
   &real&
   &_1stage_&
   &double&
   &_impl (obj, &
      a, &
      ev, &
      q) result(success)
      use precision
      use cuda_functions
      use mod_check_for_gpu
      use, intrinsic :: iso_c_binding
      use elpa_abstract_impl
      use elpa_mpi
      use elpa1_compute
      use elpa_omp

      implicit none
      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: rck = C_DOUBLE
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout)                         :: obj
      real(kind=rk8), intent(out)                              :: ev(obj%na)

      real(kind=rck), intent(inout), target                     :: a(obj%local_nrows,*)
      real(kind=rck), optional,target,intent(out)               :: q(obj%local_nrows,*)

      real(kind=c_double), allocatable           :: tau(:)
      real(kind=c_double), allocatable, target   :: q_dummy(:,:)
      real(kind=c_double), pointer               :: q_actual(:,:)

      integer(kind=c_int)                             :: l_cols, l_rows, l_cols_nev, np_rows, np_cols
      integer(kind=MPI_KIND)                          :: np_rowsMPI, np_colsMPI

      logical                                         :: useGPU
      integer(kind=c_int)                             :: skewsymmetric
      logical                                         :: isSkewsymmetric
      logical                                         :: success

      logical                                         :: do_useGPU, do_useGPU_tridiag, &
         do_useGPU_solve_tridi, do_useGPU_trans_ev
      integer(kind=ik)                                :: numberOfGPUDevices

      integer(kind=c_int)                             :: my_pe, n_pes, my_prow, my_pcol
      integer(kind=MPI_KIND)                          :: mpierr, my_peMPI, n_pesMPI, my_prowMPI, my_pcolMPI
      real(kind=c_double), allocatable         :: e(:)
      logical                                         :: wantDebug
      integer(kind=c_int)                             :: istat, debug, gpu
      character(200)                                  :: errorMessage
      integer(kind=ik)                                :: na, nev, nblk, matrixCols, &
         mpi_comm_rows, mpi_comm_cols,        &
         mpi_comm_all, check_pd, i, error, matrixRows

      logical                                         :: do_tridiag, do_solve, do_trans_ev
      integer(kind=ik)                                :: nrThreads
      integer(kind=ik)                                :: global_index

      logical                                         :: reDistributeMatrix, doRedistributeMatrix

      call obj%timer%start("elpa_solve_evp_&
      &real&
      &_1stage_&
      &double&
      &")

      reDistributeMatrix = .false.

      matrixRows = obj%local_nrows
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_parent", mpi_comm_all, error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), my_peMPI, mpierr)
      my_pe = int(my_peMPI,kind=c_int)

      nrThreads = 1

      success = .true.

      if (present(q)) then
         obj%eigenvalues_only = .false.
      else
         obj%eigenvalues_only = .true.
      endif

      na         = obj%na
      nev        = obj%nev
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      ! special case na = 1
      if (na .eq. 1) then
         ev(1) = a(1,1)
         if (.not.(obj%eigenvalues_only)) then
            q(1,1) = ONE
         endif

         ! restore original OpenMP settings
         call obj%timer%stop("elpa_solve_evp_&
         &real&
         &_1stage_&
         &double&
         &")
         return
      endif

      if (nev == 0) then
         nev = 1
         obj%eigenvalues_only = .true.
      endif

      call obj%get("gpu",gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif
      if (gpu .eq. 1) then
         useGPU =.true.
      else
         useGPU = .false.
      endif

      call obj%get("is_skewsymmetric",skewsymmetric,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for skewsymmetric. Aborting..."
         stop
      endif

      isSkewsymmetric = (skewsymmetric == 1)

      call obj%timer%start("mpi_communication")

      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)

      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      np_rows = int(np_rowsMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem setting option for debug. Aborting..."
         stop
      endif
      wantDebug = debug == 1
      do_useGPU = .false.


      if (useGPU) then
         call obj%timer%start("check_for_gpu")

         if (check_for_gpu(my_pe,numberOfGPUDevices, wantDebug=wantDebug)) then
            do_useGPU = .true.
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
      endif

      do_useGPU_tridiag = do_useGPU
      do_useGPU_solve_tridi = do_useGPU
      do_useGPU_trans_ev = do_useGPU
      ! only if we want (and can) use GPU in general, look what are the
      ! requirements for individual routines. Implicitly they are all set to 1, so
      ! unles specified otherwise by the user, GPU versions of all individual
      ! routines should be used
      if(do_useGPU) then
         call obj%get("gpu_tridiag", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_tridiag. Aborting..."
            stop
         endif
         do_useGPU_tridiag = (gpu == 1)

         call obj%get("gpu_solve_tridi", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_solve_tridi. Aborting..."
            stop
         endif
         do_useGPU_solve_tridi = (gpu == 1)

         call obj%get("gpu_trans_ev", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_trans_ev. Aborting..."
            stop
         endif
         do_useGPU_trans_ev = (gpu == 1)
      endif
      ! for elpa1 the easy thing is, that the individual phases of the algorithm
      ! do not share any data on the GPU.

      ! allocate a dummy q_intern, if eigenvectors should not be commputed and thus q is NOT present
      if (.not.(obj%eigenvalues_only)) then
         q_actual => q(1:matrixRows,1:matrixCols)
      else
         allocate(q_dummy(1:matrixRows,1:matrixCols), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa1_template: q_dummy", 386,  istat,  errorMessage)
         q_actual => q_dummy
      endif

      allocate(e(na), tau(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa1_template: e, tau", 404,  istat,  errorMessage)

      ! start the computations
      ! as default do all three steps (this might change at some point)
      do_tridiag  = .true.
      do_solve    = .true.
      do_trans_ev = .true.

      if (do_tridiag) then
         call obj%timer%start("forward")

         call tridiag_&
         &real&
         &_&
         &double&
         & (obj, na, a, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, ev, e, tau, do_useGPU_tridiag, wantDebug, &
            nrThreads)

         call obj%timer%stop("forward")
      endif  !do_tridiag

      if (do_solve) then
         call obj%timer%start("solve")
         call solve_tridi_&
         &double&
         & (obj, na, nev, ev, e,  &
            q_actual, matrixRows,          &
            nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_solve_tridi, wantDebug, success, nrThreads)

         call obj%timer%stop("solve")
         if (.not.(success)) return
      endif !do_solve

      if (obj%eigenvalues_only) then
         do_trans_ev = .false.
      else
         call obj%get("check_pd",check_pd,error)
         if (error .ne. ELPA_OK) then
            print *,"Problem setting option for check_pd. Aborting..."
            stop
         endif
         if (check_pd .eq. 1) then
            check_pd = 0
            do i = 1, na
               if (ev(i) .gt. 1e-11_rk8) then
                  check_pd = check_pd + 1
               endif
            enddo
            if (check_pd .lt. na) then
               ! not positiv definite => eigenvectors needed
               do_trans_ev = .true.
            else
               do_trans_ev = .false.
            endif
         endif ! check_pd
      endif ! eigenvalues_only

      if (do_trans_ev) then
         ! q must be given thats why from here on we can use q and not q_actual
         if (isSkewsymmetric) then
            ! Extra transformation step for skew-symmetric matrix. Multiplication with diagonal complex matrix D.
            ! This makes the eigenvectors complex.
            ! For now real part of eigenvectors is generated in first half of q, imaginary part in second part.
            q(1:matrixRows, matrixCols+1:2*matrixCols) = 0.0
            do i = 1, matrixRows
!        global_index = indxl2g(i, nblk, my_prow, 0, np_rows)
               global_index = np_rows*nblk*((i-1)/nblk) + MOD(i-1,nblk) + MOD(np_rows+my_prow-0, np_rows)*nblk + 1
               if (mod(global_index-1,4) .eq. 0) then
                  ! do nothing
               end if
               if (mod(global_index-1,4) .eq. 1) then
                  q(i,matrixCols+1:2*matrixCols) = q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
               if (mod(global_index-1,4) .eq. 2) then
                  q(i,1:matrixCols) = -q(i,1:matrixCols)
               end if
               if (mod(global_index-1,4) .eq. 3) then
                  q(i,matrixCols+1:2*matrixCols) = -q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
            end do
         endif

         call obj%timer%start("back")
         ! In the skew-symmetric case this transforms the real part
         call trans_ev_&
         &real&
         &_&
         &double&
         & (obj, na, nev, a, matrixRows, tau, q, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         if (isSkewsymmetric) then
            ! Transform imaginary part
            ! Transformation of real and imaginary part could also be one call of trans_ev_tridi acting on the n x 2n matrix.
            call trans_ev_&
            &real&
            &_&
            &double&
            & (obj, na, nev, a, matrixRows, tau, q(1:matrixRows, matrixCols+1:2*matrixCols), matrixRows, nblk, matrixCols, &
               mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         endif

         call obj%timer%stop("back")
      endif ! do_trans_ev

      deallocate(e, tau, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa1_template: e, tau", 558,  istat,  errorMessage)

      if (obj%eigenvalues_only) then
         deallocate(q_dummy, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa1_template: q_dummy", 562,  istat,  errorMessage)
      endif

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_solve_evp_&
      &real&
      &_1stage_&
      &double&
      &")
   end function

!> \brief elpa_solve_evp_real_1stage_single_impl: Fortran function to solve the real single-precision eigenvalue problem with 1-stage solver
!> \details
!> \param  obj                      elpa_t object contains:
!> \param     - obj%na              Order of matrix
!> \param     - obj%nev             number of eigenvalues/vectors to be computed
!>                                  The smallest nev eigenvalues/eigenvectors are calculated.
!> \param     - obj%local_nrows     Leading dimension of a
!> \param     - obj%local_ncols     local columns of matrix q
!> \param     - obj%nblk            blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows   MPI communicator for rows
!> \param     - obj%mpi_comm_cols   MPI communicator for columns
!> \param     - obj%mpi_comm_parent MPI communicator for columns
!> \param     - obj%gpu             use GPU version (1 or 0)
!>
!> \param  a(lda,matrixCols)        Distributed matrix for which eigenvalues are to be computed.
!>                                  Distribution is like in Scalapack.
!>                                  The full matrix must be set (not only one half like in scalapack).
!>                                  Destroyed on exit (upper and lower half).
!>
!>  \param ev(na)                   On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)        On output: Eigenvectors of a
!>                                  Distribution is like in Scalapack.
!>                                  Must be always dimensioned to the full size (corresponding to (na,na))
!>                                  even if only a part of the eigenvalues is needed.
!>
!>
!>  \result                       success

   function elpa_solve_evp_&
   &real&
   &_1stage_&
   &single&
   &_impl (obj, &
      a, &
      ev, &
      q) result(success)
      use precision
      use cuda_functions
      use mod_check_for_gpu
      use, intrinsic :: iso_c_binding
      use elpa_abstract_impl
      use elpa_mpi
      use elpa1_compute
      use elpa_omp

      implicit none
      integer, parameter :: rk = C_FLOAT
      integer, parameter :: rck = C_FLOAT
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout)                         :: obj
      real(kind=rk4), intent(out)                              :: ev(obj%na)

      real(kind=rck), intent(inout), target                     :: a(obj%local_nrows,*)
      real(kind=rck), optional,target,intent(out)               :: q(obj%local_nrows,*)

      real(kind=c_float), allocatable           :: tau(:)
      real(kind=c_float), allocatable, target   :: q_dummy(:,:)
      real(kind=c_float), pointer               :: q_actual(:,:)

      integer(kind=c_int)                             :: l_cols, l_rows, l_cols_nev, np_rows, np_cols
      integer(kind=MPI_KIND)                          :: np_rowsMPI, np_colsMPI

      logical                                         :: useGPU
      integer(kind=c_int)                             :: skewsymmetric
      logical                                         :: isSkewsymmetric
      logical                                         :: success

      logical                                         :: do_useGPU, do_useGPU_tridiag, &
         do_useGPU_solve_tridi, do_useGPU_trans_ev
      integer(kind=ik)                                :: numberOfGPUDevices

      integer(kind=c_int)                             :: my_pe, n_pes, my_prow, my_pcol
      integer(kind=MPI_KIND)                          :: mpierr, my_peMPI, n_pesMPI, my_prowMPI, my_pcolMPI
      real(kind=c_float), allocatable         :: e(:)
      logical                                         :: wantDebug
      integer(kind=c_int)                             :: istat, debug, gpu
      character(200)                                  :: errorMessage
      integer(kind=ik)                                :: na, nev, nblk, matrixCols, &
         mpi_comm_rows, mpi_comm_cols,        &
         mpi_comm_all, check_pd, i, error, matrixRows

      logical                                         :: do_tridiag, do_solve, do_trans_ev
      integer(kind=ik)                                :: nrThreads
      integer(kind=ik)                                :: global_index

      logical                                         :: reDistributeMatrix, doRedistributeMatrix

      call obj%timer%start("elpa_solve_evp_&
      &real&
      &_1stage_&
      &single&
      &")

      reDistributeMatrix = .false.

      matrixRows = obj%local_nrows
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_parent", mpi_comm_all, error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), my_peMPI, mpierr)
      my_pe = int(my_peMPI,kind=c_int)

      nrThreads = 1

      success = .true.

      if (present(q)) then
         obj%eigenvalues_only = .false.
      else
         obj%eigenvalues_only = .true.
      endif

      na         = obj%na
      nev        = obj%nev
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      ! special case na = 1
      if (na .eq. 1) then
         ev(1) = a(1,1)
         if (.not.(obj%eigenvalues_only)) then
            q(1,1) = ONE
         endif

         ! restore original OpenMP settings
         call obj%timer%stop("elpa_solve_evp_&
         &real&
         &_1stage_&
         &single&
         &")
         return
      endif

      if (nev == 0) then
         nev = 1
         obj%eigenvalues_only = .true.
      endif

      call obj%get("gpu",gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif
      if (gpu .eq. 1) then
         useGPU =.true.
      else
         useGPU = .false.
      endif

      call obj%get("is_skewsymmetric",skewsymmetric,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for skewsymmetric. Aborting..."
         stop
      endif

      isSkewsymmetric = (skewsymmetric == 1)

      call obj%timer%start("mpi_communication")

      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)

      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      np_rows = int(np_rowsMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem setting option for debug. Aborting..."
         stop
      endif
      wantDebug = debug == 1
      do_useGPU = .false.


      if (useGPU) then
         call obj%timer%start("check_for_gpu")

         if (check_for_gpu(my_pe,numberOfGPUDevices, wantDebug=wantDebug)) then
            do_useGPU = .true.
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
      endif

      do_useGPU_tridiag = do_useGPU
      do_useGPU_solve_tridi = do_useGPU
      do_useGPU_trans_ev = do_useGPU
      ! only if we want (and can) use GPU in general, look what are the
      ! requirements for individual routines. Implicitly they are all set to 1, so
      ! unles specified otherwise by the user, GPU versions of all individual
      ! routines should be used
      if(do_useGPU) then
         call obj%get("gpu_tridiag", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_tridiag. Aborting..."
            stop
         endif
         do_useGPU_tridiag = (gpu == 1)

         call obj%get("gpu_solve_tridi", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_solve_tridi. Aborting..."
            stop
         endif
         do_useGPU_solve_tridi = (gpu == 1)

         call obj%get("gpu_trans_ev", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_trans_ev. Aborting..."
            stop
         endif
         do_useGPU_trans_ev = (gpu == 1)
      endif
      ! for elpa1 the easy thing is, that the individual phases of the algorithm
      ! do not share any data on the GPU.

      ! allocate a dummy q_intern, if eigenvectors should not be commputed and thus q is NOT present
      if (.not.(obj%eigenvalues_only)) then
         q_actual => q(1:matrixRows,1:matrixCols)
      else
         allocate(q_dummy(1:matrixRows,1:matrixCols), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa1_template: q_dummy", 386,  istat,  errorMessage)
         q_actual => q_dummy
      endif

      allocate(e(na), tau(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa1_template: e, tau", 404,  istat,  errorMessage)

      ! start the computations
      ! as default do all three steps (this might change at some point)
      do_tridiag  = .true.
      do_solve    = .true.
      do_trans_ev = .true.

      if (do_tridiag) then
         call obj%timer%start("forward")

         call tridiag_&
         &real&
         &_&
         &single&
         & (obj, na, a, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, ev, e, tau, do_useGPU_tridiag, wantDebug, &
            nrThreads)

         call obj%timer%stop("forward")
      endif  !do_tridiag

      if (do_solve) then
         call obj%timer%start("solve")
         call solve_tridi_&
         &single&
         & (obj, na, nev, ev, e,  &
            q_actual, matrixRows,          &
            nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_solve_tridi, wantDebug, success, nrThreads)

         call obj%timer%stop("solve")
         if (.not.(success)) return
      endif !do_solve

      if (obj%eigenvalues_only) then
         do_trans_ev = .false.
      else
         call obj%get("check_pd",check_pd,error)
         if (error .ne. ELPA_OK) then
            print *,"Problem setting option for check_pd. Aborting..."
            stop
         endif
         if (check_pd .eq. 1) then
            check_pd = 0
            do i = 1, na
               if (ev(i) .gt. 1e-4_rk4) then
                  check_pd = check_pd + 1
               endif
            enddo
            if (check_pd .lt. na) then
               ! not positiv definite => eigenvectors needed
               do_trans_ev = .true.
            else
               do_trans_ev = .false.
            endif
         endif ! check_pd
      endif ! eigenvalues_only

      if (do_trans_ev) then
         ! q must be given thats why from here on we can use q and not q_actual
         if (isSkewsymmetric) then
            ! Extra transformation step for skew-symmetric matrix. Multiplication with diagonal complex matrix D.
            ! This makes the eigenvectors complex.
            ! For now real part of eigenvectors is generated in first half of q, imaginary part in second part.
            q(1:matrixRows, matrixCols+1:2*matrixCols) = 0.0
            do i = 1, matrixRows
!        global_index = indxl2g(i, nblk, my_prow, 0, np_rows)
               global_index = np_rows*nblk*((i-1)/nblk) + MOD(i-1,nblk) + MOD(np_rows+my_prow-0, np_rows)*nblk + 1
               if (mod(global_index-1,4) .eq. 0) then
                  ! do nothing
               end if
               if (mod(global_index-1,4) .eq. 1) then
                  q(i,matrixCols+1:2*matrixCols) = q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
               if (mod(global_index-1,4) .eq. 2) then
                  q(i,1:matrixCols) = -q(i,1:matrixCols)
               end if
               if (mod(global_index-1,4) .eq. 3) then
                  q(i,matrixCols+1:2*matrixCols) = -q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
            end do
         endif

         call obj%timer%start("back")
         ! In the skew-symmetric case this transforms the real part
         call trans_ev_&
         &real&
         &_&
         &single&
         & (obj, na, nev, a, matrixRows, tau, q, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         if (isSkewsymmetric) then
            ! Transform imaginary part
            ! Transformation of real and imaginary part could also be one call of trans_ev_tridi acting on the n x 2n matrix.
            call trans_ev_&
            &real&
            &_&
            &single&
            & (obj, na, nev, a, matrixRows, tau, q(1:matrixRows, matrixCols+1:2*matrixCols), matrixRows, nblk, matrixCols, &
               mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         endif

         call obj%timer%stop("back")
      endif ! do_trans_ev

      deallocate(e, tau, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa1_template: e, tau", 558,  istat,  errorMessage)

      if (obj%eigenvalues_only) then
         deallocate(q_dummy, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa1_template: q_dummy", 562,  istat,  errorMessage)
      endif

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_solve_evp_&
      &real&
      &_1stage_&
      &single&
      &")
   end function

!> \brief elpa_solve_evp_complex_1stage_double_impl: Fortran function to solve the complex double-precision eigenvalue problem with 1-stage solver
!> \details
!> \param  obj                      elpa_t object contains:
!> \param     - obj%na              Order of matrix
!> \param     - obj%nev             number of eigenvalues/vectors to be computed
!>                                  The smallest nev eigenvalues/eigenvectors are calculated.
!> \param     - obj%local_nrows     Leading dimension of a
!> \param     - obj%local_ncols     local columns of matrix q
!> \param     - obj%nblk            blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows   MPI communicator for rows
!> \param     - obj%mpi_comm_cols   MPI communicator for columns
!> \param     - obj%mpi_comm_parent MPI communicator for columns
!> \param     - obj%gpu             use GPU version (1 or 0)
!>
!> \param  a(lda,matrixCols)        Distributed matrix for which eigenvalues are to be computed.
!>                                  Distribution is like in Scalapack.
!>                                  The full matrix must be set (not only one half like in scalapack).
!>                                  Destroyed on exit (upper and lower half).
!>
!>  \param ev(na)                   On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)        On output: Eigenvectors of a
!>                                  Distribution is like in Scalapack.
!>                                  Must be always dimensioned to the full size (corresponding to (na,na))
!>                                  even if only a part of the eigenvalues is needed.
!>
!>
!>  \result                       success

   function elpa_solve_evp_&
   &complex&
   &_1stage_&
   &double&
   &_impl (obj, &
      a, &
      ev, &
      q) result(success)
      use precision
      use cuda_functions
      use mod_check_for_gpu
      use, intrinsic :: iso_c_binding
      use elpa_abstract_impl
      use elpa_mpi
      use elpa1_compute
      use elpa_omp

      implicit none

      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: ck = C_DOUBLE_COMPLEX
      integer, parameter :: rck = C_DOUBLE_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout)                         :: obj
      real(kind=RK8), intent(out)                              :: ev(obj%na)

      complex(kind=rck), intent(inout), target                     :: a(obj%local_nrows,*)
      complex(kind=rck), optional,target,intent(out)               :: q(obj%local_nrows,*)

      real(kind=RK8), allocatable             :: q_real(:,:)
      complex(kind=c_double), allocatable        :: tau(:)
      complex(kind=c_double), allocatable,target :: q_dummy(:,:)
      complex(kind=c_double), pointer            :: q_actual(:,:)

      integer(kind=c_int)                             :: l_cols, l_rows, l_cols_nev, np_rows, np_cols
      integer(kind=MPI_KIND)                          :: np_rowsMPI, np_colsMPI

      logical                                         :: useGPU
      integer(kind=c_int)                             :: skewsymmetric
      logical                                         :: isSkewsymmetric
      logical                                         :: success

      logical                                         :: do_useGPU, do_useGPU_tridiag, &
         do_useGPU_solve_tridi, do_useGPU_trans_ev
      integer(kind=ik)                                :: numberOfGPUDevices

      integer(kind=c_int)                             :: my_pe, n_pes, my_prow, my_pcol
      integer(kind=MPI_KIND)                          :: mpierr, my_peMPI, n_pesMPI, my_prowMPI, my_pcolMPI
      real(kind=c_double), allocatable         :: e(:)
      logical                                         :: wantDebug
      integer(kind=c_int)                             :: istat, debug, gpu
      character(200)                                  :: errorMessage
      integer(kind=ik)                                :: na, nev, nblk, matrixCols, &
         mpi_comm_rows, mpi_comm_cols,        &
         mpi_comm_all, check_pd, i, error, matrixRows

      logical                                         :: do_tridiag, do_solve, do_trans_ev
      integer(kind=ik)                                :: nrThreads
      integer(kind=ik)                                :: global_index

      logical                                         :: reDistributeMatrix, doRedistributeMatrix

      call obj%timer%start("elpa_solve_evp_&
      &complex&
      &_1stage_&
      &double&
      &")

      reDistributeMatrix = .false.

      matrixRows = obj%local_nrows
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_parent", mpi_comm_all, error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), my_peMPI, mpierr)
      my_pe = int(my_peMPI,kind=c_int)

      nrThreads = 1

      success = .true.

      if (present(q)) then
         obj%eigenvalues_only = .false.
      else
         obj%eigenvalues_only = .true.
      endif

      na         = obj%na
      nev        = obj%nev
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      ! special case na = 1
      if (na .eq. 1) then
         ev(1) = real(a(1,1))
         if (.not.(obj%eigenvalues_only)) then
            q(1,1) = ONE
         endif

         ! restore original OpenMP settings
         call obj%timer%stop("elpa_solve_evp_&
         &complex&
         &_1stage_&
         &double&
         &")
         return
      endif

      if (nev == 0) then
         nev = 1
         obj%eigenvalues_only = .true.
      endif

      call obj%get("gpu",gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif
      if (gpu .eq. 1) then
         useGPU =.true.
      else
         useGPU = .false.
      endif

      call obj%get("is_skewsymmetric",skewsymmetric,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for skewsymmetric. Aborting..."
         stop
      endif

      isSkewsymmetric = (skewsymmetric == 1)

      call obj%timer%start("mpi_communication")

      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)

      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      np_rows = int(np_rowsMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem setting option for debug. Aborting..."
         stop
      endif
      wantDebug = debug == 1
      do_useGPU = .false.


      if (useGPU) then
         call obj%timer%start("check_for_gpu")

         if (check_for_gpu(my_pe,numberOfGPUDevices, wantDebug=wantDebug)) then
            do_useGPU = .true.
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
      endif

      do_useGPU_tridiag = do_useGPU
      do_useGPU_solve_tridi = do_useGPU
      do_useGPU_trans_ev = do_useGPU
      ! only if we want (and can) use GPU in general, look what are the
      ! requirements for individual routines. Implicitly they are all set to 1, so
      ! unles specified otherwise by the user, GPU versions of all individual
      ! routines should be used
      if(do_useGPU) then
         call obj%get("gpu_tridiag", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_tridiag. Aborting..."
            stop
         endif
         do_useGPU_tridiag = (gpu == 1)

         call obj%get("gpu_solve_tridi", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_solve_tridi. Aborting..."
            stop
         endif
         do_useGPU_solve_tridi = (gpu == 1)

         call obj%get("gpu_trans_ev", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_trans_ev. Aborting..."
            stop
         endif
         do_useGPU_trans_ev = (gpu == 1)
      endif
      ! for elpa1 the easy thing is, that the individual phases of the algorithm
      ! do not share any data on the GPU.

      ! allocate a dummy q_intern, if eigenvectors should not be commputed and thus q is NOT present
      if (.not.(obj%eigenvalues_only)) then
         q_actual => q(1:matrixRows,1:matrixCols)
      else
         allocate(q_dummy(1:matrixRows,1:matrixCols), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa1_template: q_dummy", 386,  istat,  errorMessage)
         q_actual => q_dummy
      endif

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

      l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

      allocate(q_real(l_rows,l_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa1_template: q_real", 401,  istat,  errorMessage)
      allocate(e(na), tau(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa1_template: e, tau", 404,  istat,  errorMessage)

      ! start the computations
      ! as default do all three steps (this might change at some point)
      do_tridiag  = .true.
      do_solve    = .true.
      do_trans_ev = .true.

      if (do_tridiag) then
         call obj%timer%start("forward")

         call tridiag_&
         &complex&
         &_&
         &double&
         & (obj, na, a, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, ev, e, tau, do_useGPU_tridiag, wantDebug, &
            nrThreads)

         call obj%timer%stop("forward")
      endif  !do_tridiag

      if (do_solve) then
         call obj%timer%start("solve")
         call solve_tridi_&
         &double&
         & (obj, na, nev, ev, e,  &
            q_real, l_rows,  &
            nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_solve_tridi, wantDebug, success, nrThreads)

         call obj%timer%stop("solve")
         if (.not.(success)) return
      endif !do_solve

      if (obj%eigenvalues_only) then
         do_trans_ev = .false.
      else
         call obj%get("check_pd",check_pd,error)
         if (error .ne. ELPA_OK) then
            print *,"Problem setting option for check_pd. Aborting..."
            stop
         endif
         if (check_pd .eq. 1) then
            check_pd = 0
            do i = 1, na
               if (ev(i) .gt. 1e-11_rk8) then
                  check_pd = check_pd + 1
               endif
            enddo
            if (check_pd .lt. na) then
               ! not positiv definite => eigenvectors needed
               do_trans_ev = .true.
            else
               do_trans_ev = .false.
            endif
         endif ! check_pd
      endif ! eigenvalues_only

      if (do_trans_ev) then
         ! q must be given thats why from here on we can use q and not q_actual
         q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)
         if (isSkewsymmetric) then
            ! Extra transformation step for skew-symmetric matrix. Multiplication with diagonal complex matrix D.
            ! This makes the eigenvectors complex.
            ! For now real part of eigenvectors is generated in first half of q, imaginary part in second part.
            q(1:matrixRows, matrixCols+1:2*matrixCols) = 0.0
            do i = 1, matrixRows
!        global_index = indxl2g(i, nblk, my_prow, 0, np_rows)
               global_index = np_rows*nblk*((i-1)/nblk) + MOD(i-1,nblk) + MOD(np_rows+my_prow-0, np_rows)*nblk + 1
               if (mod(global_index-1,4) .eq. 0) then
                  ! do nothing
               end if
               if (mod(global_index-1,4) .eq. 1) then
                  q(i,matrixCols+1:2*matrixCols) = q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
               if (mod(global_index-1,4) .eq. 2) then
                  q(i,1:matrixCols) = -q(i,1:matrixCols)
               end if
               if (mod(global_index-1,4) .eq. 3) then
                  q(i,matrixCols+1:2*matrixCols) = -q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
            end do
         endif

         call obj%timer%start("back")
         ! In the skew-symmetric case this transforms the real part
         call trans_ev_&
         &complex&
         &_&
         &double&
         & (obj, na, nev, a, matrixRows, tau, q, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         if (isSkewsymmetric) then
            ! Transform imaginary part
            ! Transformation of real and imaginary part could also be one call of trans_ev_tridi acting on the n x 2n matrix.
            call trans_ev_&
            &complex&
            &_&
            &double&
            & (obj, na, nev, a, matrixRows, tau, q(1:matrixRows, matrixCols+1:2*matrixCols), matrixRows, nblk, matrixCols, &
               mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         endif

         call obj%timer%stop("back")
      endif ! do_trans_ev

      deallocate(q_real, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa1_template: q_real", 554,  istat,  errorMessage)

      deallocate(e, tau, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa1_template: e, tau", 558,  istat,  errorMessage)

      if (obj%eigenvalues_only) then
         deallocate(q_dummy, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa1_template: q_dummy", 562,  istat,  errorMessage)
      endif

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_solve_evp_&
      &complex&
      &_1stage_&
      &double&
      &")
   end function

!> \brief elpa_solve_evp_complex_1stage_single_impl: Fortran function to solve the complex single-precision eigenvalue problem with 1-stage solver
!> \details
!> \param  obj                      elpa_t object contains:
!> \param     - obj%na              Order of matrix
!> \param     - obj%nev             number of eigenvalues/vectors to be computed
!>                                  The smallest nev eigenvalues/eigenvectors are calculated.
!> \param     - obj%local_nrows     Leading dimension of a
!> \param     - obj%local_ncols     local columns of matrix q
!> \param     - obj%nblk            blocksize of cyclic distribution, must be the same in both directions!
!> \param     - obj%mpi_comm_rows   MPI communicator for rows
!> \param     - obj%mpi_comm_cols   MPI communicator for columns
!> \param     - obj%mpi_comm_parent MPI communicator for columns
!> \param     - obj%gpu             use GPU version (1 or 0)
!>
!> \param  a(lda,matrixCols)        Distributed matrix for which eigenvalues are to be computed.
!>                                  Distribution is like in Scalapack.
!>                                  The full matrix must be set (not only one half like in scalapack).
!>                                  Destroyed on exit (upper and lower half).
!>
!>  \param ev(na)                   On output: eigenvalues of a, every processor gets the complete set
!>
!>  \param q(ldq,matrixCols)        On output: Eigenvectors of a
!>                                  Distribution is like in Scalapack.
!>                                  Must be always dimensioned to the full size (corresponding to (na,na))
!>                                  even if only a part of the eigenvalues is needed.
!>
!>
!>  \result                       success

   function elpa_solve_evp_&
   &complex&
   &_1stage_&
   &single&
   &_impl (obj, &
      a, &
      ev, &
      q) result(success)
      use precision
      use cuda_functions
      use mod_check_for_gpu
      use, intrinsic :: iso_c_binding
      use elpa_abstract_impl
      use elpa_mpi
      use elpa1_compute
      use elpa_omp

      implicit none

      integer, parameter :: rk = C_FLOAT
      integer, parameter :: ck = C_FLOAT_COMPLEX
      integer, parameter :: rck = C_FLOAT_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
      class(elpa_abstract_impl_t), intent(inout)                         :: obj
      real(kind=RK4), intent(out)                              :: ev(obj%na)

      complex(kind=rck), intent(inout), target                     :: a(obj%local_nrows,*)
      complex(kind=rck), optional,target,intent(out)               :: q(obj%local_nrows,*)

      real(kind=RK4), allocatable             :: q_real(:,:)
      complex(kind=c_float), allocatable        :: tau(:)
      complex(kind=c_float), allocatable,target :: q_dummy(:,:)
      complex(kind=c_float), pointer            :: q_actual(:,:)

      integer(kind=c_int)                             :: l_cols, l_rows, l_cols_nev, np_rows, np_cols
      integer(kind=MPI_KIND)                          :: np_rowsMPI, np_colsMPI

      logical                                         :: useGPU
      integer(kind=c_int)                             :: skewsymmetric
      logical                                         :: isSkewsymmetric
      logical                                         :: success

      logical                                         :: do_useGPU, do_useGPU_tridiag, &
         do_useGPU_solve_tridi, do_useGPU_trans_ev
      integer(kind=ik)                                :: numberOfGPUDevices

      integer(kind=c_int)                             :: my_pe, n_pes, my_prow, my_pcol
      integer(kind=MPI_KIND)                          :: mpierr, my_peMPI, n_pesMPI, my_prowMPI, my_pcolMPI
      real(kind=c_float), allocatable         :: e(:)
      logical                                         :: wantDebug
      integer(kind=c_int)                             :: istat, debug, gpu
      character(200)                                  :: errorMessage
      integer(kind=ik)                                :: na, nev, nblk, matrixCols, &
         mpi_comm_rows, mpi_comm_cols,        &
         mpi_comm_all, check_pd, i, error, matrixRows

      logical                                         :: do_tridiag, do_solve, do_trans_ev
      integer(kind=ik)                                :: nrThreads
      integer(kind=ik)                                :: global_index

      logical                                         :: reDistributeMatrix, doRedistributeMatrix

      call obj%timer%start("elpa_solve_evp_&
      &complex&
      &_1stage_&
      &single&
      &")

      reDistributeMatrix = .false.

      matrixRows = obj%local_nrows
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_parent", mpi_comm_all, error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND), my_peMPI, mpierr)
      my_pe = int(my_peMPI,kind=c_int)

      nrThreads = 1

      success = .true.

      if (present(q)) then
         obj%eigenvalues_only = .false.
      else
         obj%eigenvalues_only = .true.
      endif

      na         = obj%na
      nev        = obj%nev
      matrixRows = obj%local_nrows
      nblk       = obj%nblk
      matrixCols = obj%local_ncols

      call obj%get("mpi_comm_rows",mpi_comm_rows,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif
      call obj%get("mpi_comm_cols",mpi_comm_cols,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option. Aborting..."
         stop
      endif

      ! special case na = 1
      if (na .eq. 1) then
         ev(1) = real(a(1,1))
         if (.not.(obj%eigenvalues_only)) then
            q(1,1) = ONE
         endif

         ! restore original OpenMP settings
         call obj%timer%stop("elpa_solve_evp_&
         &complex&
         &_1stage_&
         &single&
         &")
         return
      endif

      if (nev == 0) then
         nev = 1
         obj%eigenvalues_only = .true.
      endif

      call obj%get("gpu",gpu,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for gpu. Aborting..."
         stop
      endif
      if (gpu .eq. 1) then
         useGPU =.true.
      else
         useGPU = .false.
      endif

      call obj%get("is_skewsymmetric",skewsymmetric,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem getting option for skewsymmetric. Aborting..."
         stop
      endif

      isSkewsymmetric = (skewsymmetric == 1)

      call obj%timer%start("mpi_communication")

      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)

      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

      np_rows = int(np_rowsMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      call obj%get("debug", debug,error)
      if (error .ne. ELPA_OK) then
         print *,"Problem setting option for debug. Aborting..."
         stop
      endif
      wantDebug = debug == 1
      do_useGPU = .false.


      if (useGPU) then
         call obj%timer%start("check_for_gpu")

         if (check_for_gpu(my_pe,numberOfGPUDevices, wantDebug=wantDebug)) then
            do_useGPU = .true.
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
      endif

      do_useGPU_tridiag = do_useGPU
      do_useGPU_solve_tridi = do_useGPU
      do_useGPU_trans_ev = do_useGPU
      ! only if we want (and can) use GPU in general, look what are the
      ! requirements for individual routines. Implicitly they are all set to 1, so
      ! unles specified otherwise by the user, GPU versions of all individual
      ! routines should be used
      if(do_useGPU) then
         call obj%get("gpu_tridiag", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_tridiag. Aborting..."
            stop
         endif
         do_useGPU_tridiag = (gpu == 1)

         call obj%get("gpu_solve_tridi", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_solve_tridi. Aborting..."
            stop
         endif
         do_useGPU_solve_tridi = (gpu == 1)

         call obj%get("gpu_trans_ev", gpu, error)
         if (error .ne. ELPA_OK) then
            print *,"Problem getting option for gpu_trans_ev. Aborting..."
            stop
         endif
         do_useGPU_trans_ev = (gpu == 1)
      endif
      ! for elpa1 the easy thing is, that the individual phases of the algorithm
      ! do not share any data on the GPU.

      ! allocate a dummy q_intern, if eigenvectors should not be commputed and thus q is NOT present
      if (.not.(obj%eigenvalues_only)) then
         q_actual => q(1:matrixRows,1:matrixCols)
      else
         allocate(q_dummy(1:matrixRows,1:matrixCols), stat=istat, errmsg=errorMessage)
         call check_allocate_f("elpa1_template: q_dummy", 386,  istat,  errorMessage)
         q_actual => q_dummy
      endif

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

      l_cols_nev = local_index(nev, my_pcol, np_cols, nblk, -1) ! Local columns corresponding to nev

      allocate(q_real(l_rows,l_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa1_template: q_real", 401,  istat,  errorMessage)
      allocate(e(na), tau(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("elpa1_template: e, tau", 404,  istat,  errorMessage)

      ! start the computations
      ! as default do all three steps (this might change at some point)
      do_tridiag  = .true.
      do_solve    = .true.
      do_trans_ev = .true.

      if (do_tridiag) then
         call obj%timer%start("forward")

         call tridiag_&
         &complex&
         &_&
         &single&
         & (obj, na, a, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, ev, e, tau, do_useGPU_tridiag, wantDebug, &
            nrThreads)

         call obj%timer%stop("forward")
      endif  !do_tridiag

      if (do_solve) then
         call obj%timer%start("solve")
         call solve_tridi_&
         &single&
         & (obj, na, nev, ev, e,  &
            q_real, l_rows,  &
            nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_solve_tridi, wantDebug, success, nrThreads)

         call obj%timer%stop("solve")
         if (.not.(success)) return
      endif !do_solve

      if (obj%eigenvalues_only) then
         do_trans_ev = .false.
      else
         call obj%get("check_pd",check_pd,error)
         if (error .ne. ELPA_OK) then
            print *,"Problem setting option for check_pd. Aborting..."
            stop
         endif
         if (check_pd .eq. 1) then
            check_pd = 0
            do i = 1, na
               if (ev(i) .gt. 1e-4_rk4) then
                  check_pd = check_pd + 1
               endif
            enddo
            if (check_pd .lt. na) then
               ! not positiv definite => eigenvectors needed
               do_trans_ev = .true.
            else
               do_trans_ev = .false.
            endif
         endif ! check_pd
      endif ! eigenvalues_only

      if (do_trans_ev) then
         ! q must be given thats why from here on we can use q and not q_actual
         q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)
         if (isSkewsymmetric) then
            ! Extra transformation step for skew-symmetric matrix. Multiplication with diagonal complex matrix D.
            ! This makes the eigenvectors complex.
            ! For now real part of eigenvectors is generated in first half of q, imaginary part in second part.
            q(1:matrixRows, matrixCols+1:2*matrixCols) = 0.0
            do i = 1, matrixRows
!        global_index = indxl2g(i, nblk, my_prow, 0, np_rows)
               global_index = np_rows*nblk*((i-1)/nblk) + MOD(i-1,nblk) + MOD(np_rows+my_prow-0, np_rows)*nblk + 1
               if (mod(global_index-1,4) .eq. 0) then
                  ! do nothing
               end if
               if (mod(global_index-1,4) .eq. 1) then
                  q(i,matrixCols+1:2*matrixCols) = q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
               if (mod(global_index-1,4) .eq. 2) then
                  q(i,1:matrixCols) = -q(i,1:matrixCols)
               end if
               if (mod(global_index-1,4) .eq. 3) then
                  q(i,matrixCols+1:2*matrixCols) = -q(i,1:matrixCols)
                  q(i,1:matrixCols) = 0
               end if
            end do
         endif

         call obj%timer%start("back")
         ! In the skew-symmetric case this transforms the real part
         call trans_ev_&
         &complex&
         &_&
         &single&
         & (obj, na, nev, a, matrixRows, tau, q, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         if (isSkewsymmetric) then
            ! Transform imaginary part
            ! Transformation of real and imaginary part could also be one call of trans_ev_tridi acting on the n x 2n matrix.
            call trans_ev_&
            &complex&
            &_&
            &single&
            & (obj, na, nev, a, matrixRows, tau, q(1:matrixRows, matrixCols+1:2*matrixCols), matrixRows, nblk, matrixCols, &
               mpi_comm_rows, mpi_comm_cols, do_useGPU_trans_ev)
         endif

         call obj%timer%stop("back")
      endif ! do_trans_ev

      deallocate(q_real, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa1_template: q_real", 554,  istat,  errorMessage)

      deallocate(e, tau, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("elpa1_template: e, tau", 558,  istat,  errorMessage)

      if (obj%eigenvalues_only) then
         deallocate(q_dummy, stat=istat, errmsg=errorMessage)
         call check_deallocate_f("elpa1_template: q_dummy", 562,  istat,  errorMessage)
      endif

      ! restore original OpenMP settings

      call obj%timer%stop("elpa_solve_evp_&
      &complex&
      &_1stage_&
      &single&
      &")
   end function

end module ELPA1_impl
