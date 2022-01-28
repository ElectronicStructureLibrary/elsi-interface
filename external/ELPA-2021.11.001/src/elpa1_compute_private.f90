










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

!> \brief Fortran module which contains the source of ELPA 1stage
module elpa1_compute
  use elpa_utilities
  use elpa_mpi
  implicit none

  PRIVATE ! set default to private

  public :: tridiag_real_double               ! Transform real symmetric matrix to tridiagonal form
  public :: tridiag_real
  public :: trans_ev_real_double              ! Transform real eigenvectors of a tridiagonal matrix back
  public :: trans_ev_real

  !public :: solve_tridi_double
  public :: solve_tridi_double_impl

  interface tridiag_real
    module procedure tridiag_real_double
  end interface

  interface trans_ev_real
    module procedure trans_ev_real_double
  end interface

  public :: tridiag_real_single        ! Transform real single-precision symmetric matrix to tridiagonal form
  public :: trans_ev_real_single       ! Transform real  single-precision eigenvectors of a tridiagonal matrix back
  !public :: solve_tridi_single
  public :: solve_tridi_single_impl

  public :: tridiag_complex_double            ! Transform complex hermitian matrix to tridiagonal form
  public :: tridiag_complex
  public :: trans_ev_complex_double           ! Transform eigenvectors of a tridiagonal matrix back
  public :: trans_ev_complex

  interface tridiag_complex
    module procedure tridiag_complex_double
  end interface

  interface trans_ev_complex
    module procedure trans_ev_complex_double
  end interface

  public :: tridiag_complex_single     ! Transform complex single-precision hermitian matrix to tridiagonal form
  public :: trans_ev_complex_single    ! Transform complex single-precision eigenvectors of a tridiagonal matrix back

  public :: hh_transform_real_double
  public :: hh_transform_real
  public :: elpa_reduce_add_vectors_real_double
  public :: elpa_reduce_add_vectors_real
  public :: elpa_transpose_vectors_real_double
  public :: elpa_transpose_vectors_ss_real_double
  public :: elpa_transpose_vectors_real
  public :: elpa_transpose_vectors_ss_real

  interface hh_transform_real
    module procedure hh_transform_real_double
  end interface

  interface elpa_reduce_add_vectors_real
    module procedure elpa_reduce_add_vectors_real_double
  end interface

  interface elpa_transpose_vectors_real
    module procedure elpa_transpose_vectors_real_double
  end interface
  
  interface elpa_transpose_vectors_ss_real
    module procedure elpa_transpose_vectors_ss_real_double
  end interface

  public :: hh_transform_real_single
  public :: elpa_reduce_add_vectors_real_single
  public :: elpa_transpose_vectors_real_single
  public :: elpa_transpose_vectors_ss_real_single

  public :: hh_transform_complex_double
  public :: hh_transform_complex
  public :: elpa_reduce_add_vectors_complex_double
  public :: elpa_reduce_add_vectors_complex
  public :: elpa_transpose_vectors_complex_double
  public :: elpa_transpose_vectors_ss_complex_double
  public :: elpa_transpose_vectors_complex
  public :: elpa_transpose_vectors_ss_complex

  interface hh_transform_complex
    module procedure hh_transform_complex_double
  end interface

  interface elpa_reduce_add_vectors_complex
    module procedure elpa_reduce_add_vectors_complex_double
  end interface

  interface elpa_transpose_vectors_complex
    module procedure elpa_transpose_vectors_complex_double
  end interface
  
  interface elpa_transpose_vectors_ss_complex
    module procedure elpa_transpose_vectors_ss_complex_double
  end interface

  public :: hh_transform_complex_single
  public :: elpa_reduce_add_vectors_complex_single
  public :: elpa_transpose_vectors_complex_single
  public :: elpa_transpose_vectors_ss_complex_single

  contains

! real double precision first






















!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&real&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL8,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL8,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_double" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&real&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL8,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL8,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_double" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&real&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  real(kind=c_double), intent(in)    :: vmat_s(ld_s,nvc)
  real(kind=c_double), intent(inout) :: vmat_t(ld_t,nvc)

  real(kind=c_double), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: auxstride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_double" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    auxstride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * auxstride

      call obj%timer%start("mpi_communication")

      if (k>0) call mpi_reduce(aux1, aux2, k, &
                    MPI_REAL8,  &
                    MPI_SUM, int(ipt,kind=MPI_KIND), int(comm_t,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")


      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif

    endif

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_double" &
  )
end subroutine


! single precision
























!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&real&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL4,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL4,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &real&
  &" // &
  &"_single" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&real&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  real(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  real(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL4,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_REAL4,    &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &real&
  &" // &
  &"_single" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&real&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  real(kind=c_float), intent(in)    :: vmat_s(ld_s,nvc)
  real(kind=c_float), intent(inout) :: vmat_t(ld_t,nvc)

  real(kind=c_float), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: auxstride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_single" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    auxstride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * auxstride

      call obj%timer%start("mpi_communication")

      if (k>0) call mpi_reduce(aux1, aux2, k, &
                    MPI_REAL4,  &
                    MPI_SUM, int(ipt,kind=MPI_KIND), int(comm_t,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")


      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif

    endif

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &real&
  &" // &
  &"_single" &
  )
end subroutine



! double precision























!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&complex&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_DOUBLE_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_DOUBLE_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_double" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&complex&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_double), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_double), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_double), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_double" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_DOUBLE_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_DOUBLE_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_double" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&complex&
&_&
&double &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  complex(kind=c_double), intent(in)    :: vmat_s(ld_s,nvc)
  complex(kind=c_double), intent(inout) :: vmat_t(ld_t,nvc)

  complex(kind=c_double), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: auxstride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_double" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    auxstride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * auxstride

      call obj%timer%start("mpi_communication")

      if (k>0) call mpi_reduce(aux1, aux2, k, &
                    MPI_DOUBLE_COMPLEX, &
                    MPI_SUM, int(ipt,kind=MPI_KIND), int(comm_t,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")


      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif

    endif

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_double" &
  )
end subroutine

























!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_&
&complex&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_&
  &complex&
  &" // &
  &"_single" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)




subroutine elpa_transpose_vectors_ss_&
&complex&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvs, nvr, nvc, nblk, nrThreads, comm_s_isRows)

!-------------------------------------------------------------------------------
! This routine transposes an array of vectors which are distributed in
! communicator comm_s into its transposed form distributed in communicator comm_t.
! There must be an identical copy of vmat_s in every communicator comm_s.
! After this routine, there is an identical copy of vmat_t in every communicator comm_t.
!
! vmat_s    original array of vectors
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors in transposed form
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvs       global index where to start in vmat_s/vmat_t
!           Please note: this is kind of a hint, some values before nvs will be
!           accessed in vmat_s/put into vmat_t
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_abstract_impl
  use elpa_mpi

  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik), intent(in)                      :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=c_float), intent(in)   :: vmat_s(ld_s,nvc)
  complex(kind=c_float), intent(inout):: vmat_t(ld_t,nvc)

  complex(kind=c_float), allocatable  :: aux(:)
  integer(kind=ik)                                  :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                            :: mypsMPI, myptMPI, npsMPI, nptMPI
  integer(kind=ik)                                  :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                            :: mpierr
  integer(kind=ik)                                  :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik)                                  :: auxstride
  integer(kind=ik), intent(in)                      :: nrThreads
  integer(kind=ik)                                  :: istat
  character(200)                                    :: errorMessage

  integer(kind=MPI_KIND)                            :: bcast_request1
  logical                                           :: useNonBlockingCollectives
  logical                                           :: useNonBlockingCollectivesRows
  logical                                           :: useNonBlockingCollectivesCols
  logical, intent(in)                               :: comm_s_isRows
  integer(kind=c_int)                               :: non_blocking_collectives_rows, error, &
                                                       non_blocking_collectives_cols

  call obj%timer%start("&
          &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_single" &
  )

  call obj%get("nbc_row_transpose_vectors", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_rows in transpose_vectors. Aborting..."
    stop
  endif

  call obj%get("nbc_col_transpose_vectors", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for_cols in transpose_vectors. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  if (comm_s_isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND),mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND),npsMPI ,mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND),myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND),nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)


  call obj%timer%stop("mpi_communication")
  ! The basic idea of this routine is that for every block (in the block cyclic
  ! distribution), the processor within comm_t which owns the diagonal
  ! broadcasts its values of vmat_s to all processors within comm_t.
  ! Of course this has not to be done for every block separately, since
  ! the communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  ! Get the number of blocks to be skipped at the begin.
  ! This must be a multiple of lcm_s_t (else it is getting complicated),
  ! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_transpose_vectors: aux", 187,  istat,  errorMessage)
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if (mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
!      if(nblks_comm==0) cycle
      if (nblks_comm .ne. 0) then
        if (myps == ips) then
!          k = 0
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!              k = k+nblk
            enddo
          enddo
        endif


        if (useNonBlockingCollectives) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND), bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(aux, int(nblks_comm*nblk*nvc,kind=MPI_KIND),    &
                      MPI_COMPLEX, &
                      int(ips,kind=MPI_KIND), int(comm_s,kind=MPI_KIND),  mpierr)
          call obj%timer%stop("mpi_communication")
        endif

!        k = 0
        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = - aux(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif
    endif

  enddo
  deallocate(aux, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_transpose_vectors: aux", 286,  istat,  errorMessage)

  call obj%timer%stop("&
  &elpa_transpose_vectors_ss_&
  &complex&
  &" // &
  &"_single" &
  )

end subroutine





!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine elpa_reduce_add_vectors_&
&complex&
&_&
&single &
(obj, vmat_s, ld_s, comm_s, vmat_t, ld_t, comm_t, nvr, nvc, nblk, nrThreads)

!-------------------------------------------------------------------------------
! This routine does a reduce of all vectors in vmat_s over the communicator comm_t.
! The result of the reduce is gathered on the processors owning the diagonal
! and added to the array of vectors vmat_t (which is distributed over comm_t).
!
! Opposed to elpa_transpose_vectors, there is NO identical copy of vmat_s
! in the different members within vmat_t (else a reduce wouldn't be necessary).
! After this routine, an allreduce of vmat_t has to be done.
!
! vmat_s    array of vectors to be reduced and added
! ld_s      leading dimension of vmat_s
! comm_s    communicator over which vmat_s is distributed
! vmat_t    array of vectors to which vmat_s is added
! ld_t      leading dimension of vmat_t
! comm_t    communicator over which vmat_t is distributed
! nvr       global length of vmat_s/vmat_t
! nvc       number of columns in vmat_s/vmat_t
! nblk      block size of block cyclic distribution
!
!-------------------------------------------------------------------------------

  use precision
  use elpa_mpi
  use elpa_abstract_impl
  implicit none

  class(elpa_abstract_impl_t), intent(inout)         :: obj
  integer(kind=ik), intent(in)                       :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  complex(kind=c_float), intent(in)    :: vmat_s(ld_s,nvc)
  complex(kind=c_float), intent(inout) :: vmat_t(ld_t,nvc)

  complex(kind=c_float), allocatable   :: aux1(:), aux2(:)
  integer(kind=ik)                                   :: myps, mypt, nps, npt
  integer(kind=MPI_KIND)                             :: mypsMPI, npsMPI, myptMPI, nptMPI
  integer(kind=ik)                                   :: n, lc, k, i, ips, ipt, ns, nl
  integer(kind=MPI_KIND)                             :: mpierr
  integer(kind=ik)                                   :: lcm_s_t, nblks_tot
  integer(kind=ik)                                   :: auxstride
  integer(kind=ik), intent(in)                       :: nrThreads
  integer(kind=ik)                                   :: istat
  character(200)                                     :: errorMessage

  call obj%timer%start("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_single" &
  )

  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(comm_s,kind=MPI_KIND), mypsMPI, mpierr)
  call mpi_comm_size(int(comm_s,kind=MPI_KIND), npsMPI,  mpierr)
  call mpi_comm_rank(int(comm_t,kind=MPI_KIND), myptMPI, mpierr)
  call mpi_comm_size(int(comm_t,kind=MPI_KIND), nptMPI ,mpierr)
  myps = int(mypsMPI,kind=c_int)
  nps = int(npsMPI,kind=c_int)
  mypt = int(myptMPI,kind=c_int)
  npt = int(nptMPI,kind=c_int)

  call obj%timer%stop("mpi_communication")

  ! Look to elpa_transpose_vectors for the basic idea!

  ! The communictation pattern repeats in the global matrix after
  ! the least common multiple of (nps,npt) blocks

  lcm_s_t   = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux1", 129,  istat,  errorMessage)

  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_reduce_add: aux2", 132,  istat,  errorMessage)
  aux1(:) = 0
  aux2(:) = 0
  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    auxstride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if (myps == ips) then

!      k = 0
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
!          k = k+nblk
        enddo
      enddo

      k = nvc * auxstride

      call obj%timer%start("mpi_communication")

      if (k>0) call mpi_reduce(aux1, aux2, k, &
                    MPI_COMPLEX, &
                    MPI_SUM, int(ipt,kind=MPI_KIND), int(comm_t,kind=MPI_KIND), mpierr)

      call obj%timer%stop("mpi_communication")


      if (mypt == ipt) then
!        k = 0
        do lc=1,nvc
          do  i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
!            k = k+nblk
          enddo
        enddo
      endif

    endif

  enddo

  deallocate(aux1, aux2, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_reduce_add: aux1, aux2", 218,  istat,  errorMessage)

  call obj%timer%stop("elpa_reduce_add_vectors_&
  &complex&
  &" // &
  &"_single" &
  )
end subroutine




! real double precision





















!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)









!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &real&
  &_&
  &double &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads, isSkewsymmetric)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: rck = C_DOUBLE
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  real(kind=rck), intent(out)          :: tau(na)
  real(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik), parameter                   :: max_stored_uv = 32
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  !integer(kind=ik)                             :: l_cols, l_rows
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs

  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: a_offset

  integer(kind=ik), intent(in)                  :: max_threads

  real(kind=rk)                                 :: vnorm2
  real(kind=rck)                       :: vav, x, aux(2*max_stored_uv), aux1(2), aux2(2), vrl, xf

  integer(kind=c_intptr_t)                      :: num
  real(kind=rck), allocatable          :: tmp(:)
  real(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  real(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  real(kind=rck), pointer             :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  real(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &real
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_&
  &real&
  &" // &
  "_double" // &
  gpuString )

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for min_tile_size. Aborting..."
    stop
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows
  !
  allocate(tmp(MAX(max_local_rows,max_local_cols)), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
  &real ", "tmp", istat, errorMessage)

  ! allocate v_row 1 element longer to allow store and broadcast tau together with it
  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &real ", "v_row", istat, errorMessage)
  
      allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
       &real ", "v_col", istat, errorMessage)
  
      allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &real ", "u_col", istat, errorMessage)
  
      allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &real ", "u_row", istat, errorMessage)
    else

      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 346,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 351,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 356,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))

      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 361,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))

      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_roes", 367,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 372,  successGPU)

      num = na * 8
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 381,  successGPU)

      num = na * 8
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 390,  successGPU)
    endif
  else
    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &real ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_row", istat, errorMessage)
      
  endif


  tmp = 0
  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU .and. .not.(useIntelGPU) ) then
     successGPU = gpu_malloc(v_row_dev, max_local_rows * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_row_dev", 429,  successGPU)

     successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

     call check_alloc_GPU_f("tridiag: u_row_dev", 433,  successGPU)

     successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_col_dev", 436,  successGPU)

     successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: u_col_dev", 439,  successGPU)

     successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 442,  successGPU)

     successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 445,  successGPU)
  endif !useGPU

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) &
  d_vec(na) = a_mat(l_rows,l_cols)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! allocate memmory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 476,  successGPU)

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("tridiag: a_mat", 480,  successGPU)

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 484,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1

    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

      ! copy l_cols + 1 column of A to v_row
      if (useGPU) then
        if (useIntelGPU) then
          v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
        else
          a_offset = l_cols * matrixRows * size_of_datatype
          ! we use v_row on the host at the moment! successGPU = gpu_memcpy(v_row_dev, a_dev + a_offset, 
          ! (l_rows)*size_of_PRECISION_real, gpuMemcpyDeviceToDevice)

          successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), &
                                    a_dev + a_offset, (l_rows)* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag a_dev 1", 519,  successGPU)
        endif
      else
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        if (useIntelGPU) then
                !print *,"intel phase aaaaaaaaaaaaaaaaaaaaaaaaaa"
          if (wantDebug) call obj%timer%start("mkl_offload")
          aux(1:2*n_stored_vecs) = uv_stored_cols(l_cols+1,1:2*n_stored_vecs)



          if (wantDebug) call obj%timer%stop("mkl_offload")
        else
          if (wantDebug) call obj%timer%start("blas")
          call DGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            uv_stored_cols(l_cols+1,1), &
                            int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND), &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif
      endif

      if (my_prow == prow(istep-1, nblk, np_rows)) then
        aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1))
        aux1(2) = v_row(l_rows)
      else
        aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
        aux1(2) = 0.
      endif

      if (useNonBlockingCollectivesRows) then
        if (wantDebug) call obj%timer%start("mpi_communication_non_blocking")
        call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL8, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
        call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication_non_blocking")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL8, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif

      vnorm2 = aux2(1)
      vrl    = aux2(2)

      ! Householder transformation
      call hh_transform_real_&
               &double &
               (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      ! Scale v_row and store Householder Vector for back transformation

      v_row(1:l_rows) = v_row(1:l_rows) * xf
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        v_row(l_rows) = 1.

        ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
        e_vec(istep-1) = vrl
      endif

      ! store Householder Vector for back transformation
      a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)

      ! add tau after the end of actuall v_row, to be broadcasted with it
      v_row(l_rows+1) = tau(istep)
    endif !(my_pcol == pcol(istep, nblk, np_cols))

!          call prmat(na, useGpu, a_mat, a_dev, matrixRows, matrixCols, nblk, my_prow, my_pcol, np_rows, np_cols, "HH vec stored",  na - istep + 1)

    if (useNonBlockingCollectivesCols) then
      if (wantDebug) call obj%timer%start("mpi_nbc_communication")
      ! Broadcast the Householder Vector (and tau) along columns
      call mpi_ibcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_REAL8,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   bcast_request1, mpierr)
      call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
      if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
    else
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_bcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_REAL8,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")
    endif

    !recover tau, which has been broadcasted together with v_row
    tau(istep) =  v_row(l_rows+1)

    ! Transpose Householder Vector v_row -> v_col
    call elpa_transpose_vectors_&
        &real&
        &_&
        &double &
              (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
               1, istep-1, 1, nblk, max_threads, .true.)

    ! Calculate u = (A + VU**T + UV**T)*v

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    u_col(1:l_cols) = 0
    u_row(1:l_rows) = 0
    if (l_rows > 0 .and. l_cols> 0 ) then
     if (useGPU .and. .not.(useIntelGPU)) then
       successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_col_dev", 693,  successGPU)

       successGPU = gpu_memset(u_row_dev, 0, l_rows * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_row_dev", 696,  successGPU)

       successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                     l_cols * size_of_datatype, gpuMemcpyHostToDevice)

       call check_memcpy_GPU_f("tridiag: v_col_dev", 701,  successGPU)

       successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                 l_rows * size_of_datatype, gpuMemcpyHostToDevice)
       call check_memcpy_GPU_f("tridiag: v_row_dev", 705,  successGPU)
     endif ! useGPU

     !if (useIntelGPU) then
     !  ! needed later when we can do explicit memcopy
     !endif

     do i= 0, (istep-2)/tile_size
       l_col_beg = i*l_cols_per_tile+1
       l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
       if (l_col_end < l_col_beg) cycle
       do j = 0, i
         l_row_beg = j*l_rows_per_tile+1
         l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
         if (l_row_end < l_row_beg) cycle

         ! multiplication by blocks is efficient only for CPU
         ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
         ! CPU implementation) or by one large matrix Vector multiply
         if (.not. useGPU) then
           if (wantDebug) call obj%timer%start("blas")
           call DGEMV('T',  &
                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                       ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                       v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                       ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

           if (i/=j) then
             if (isSkewsymmetric) then
               call DGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                   -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)

             else
               call DGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                   ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)
             endif
           endif
           if (wantDebug) call obj%timer%stop("blas")
         endif ! not useGPU

            enddo  ! j=0,i
          enddo  ! i=0,(istep-2)/tile_size

          if (useGPU) then
            if (mat_vec_as_one_block) then
              if (useIntelGPU) then
                 if (wantDebug) call obj%timer%start("mkl_offload")

                if (wantDebug) call obj%timer%stop("mkl_offload")

              else
                ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
                ! this requires altering of the algorithm when later explicitly updating the matrix
                ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
                if (wantDebug) call obj%timer%start("gpublas")
                call gpublas_DGEMV('T', l_rows,l_cols,  &
                                          ONE, a_dev, matrixRows,                   &
                                          v_row_dev , 1,                          &
                                          ONE, u_col_dev, 1)

       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_DGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + a_offset, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
                if (wantDebug) call obj%timer%stop("gpublas")
              endif
            else  ! mat_vec_as_one_block
              !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
              !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")

                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype

                    call gpublas_DGEMV('T', &
                                l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                ONE, a_dev + a_offset, matrixRows,  &
                                v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                                ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1)
                endif
              enddo

              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,i*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")


                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype
                    if (isSkewsymmetric) then
                       call gpublas_DGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   -ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                    else
                       call gpublas_DGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                   endif
                endif
              enddo
            end if !multiplication as one block / per stripes

            if (.not.(useIntelGPU)) then
              successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), &
                          u_col_dev, l_cols * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_col_dev 1", 954,  successGPU)

              successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), &
                          u_row_dev, l_rows * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_row_dev 1", 958,  successGPU)
            endif
            !if (useIntelGPU) then
            !  
            !endif

          endif ! useGPU


         ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
         if (n_stored_vecs > 0) then
           if (wantDebug) call obj%timer%start("blas")
           call DGEMV('T',     &
                               int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                               v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)

           call DGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                               aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
           if (wantDebug) call obj%timer%stop("blas")
         endif

       endif  ! (l_rows>0 .and. l_cols>0)

       ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
       ! on the processors containing the diagonal
       ! This is only necessary if u_row has been calculated, i.e. if the
       ! global tile size is smaller than the global remaining matrix

       if (tile_size < istep-1) then

         call elpa_reduce_add_vectors_&
         &real&
         &_&
         &double &
         (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
         mpi_comm_cols, istep-1, 1, nblk, max_threads)

       endif

       ! Sum up all the u_col(:) parts, transpose u_col -> u_row

       if (l_cols>0) then
         tmp(1:l_cols) = u_col(1:l_cols)
         if (useNonBlockingCollectivesRows) then
           if (wantDebug) call obj%timer%start("mpi_nbc_communication")
           call mpi_iallreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_REAL8,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
           call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
         else
           if (wantDebug) call obj%timer%start("mpi_communication")
           call mpi_allreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_REAL8,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
           if (wantDebug) call obj%timer%stop("mpi_communication")
         endif
       endif
       if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
          &real&
          &_&
          &double &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       else
          call elpa_transpose_vectors_&
          &real&
          &_&
          &double &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       endif

       ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
       x = 0
       if (l_cols>0)  &
       x = dot_product(v_col(1:l_cols),u_col(1:l_cols))

       if (useNonBlockingCollectivesCols) then
         if (wantDebug) call obj%timer%start("mpi_nbc_communication")
         call mpi_iallreduce(x, vav, 1_MPI_KIND, MPI_REAL8, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
               allreduce_request3, mpierr)
         call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
       else
         if (wantDebug) call obj%timer%start("mpi_communication")
         call mpi_allreduce(x, vav, 1_MPI_KIND, MPI_REAL8, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
                mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")
       endif

       ! store u and v in the matrices U and V
       ! these matrices are stored combined in one here

       do j=1,l_rows
         vu_stored_rows(j,2*n_stored_vecs+1) = tau(istep)*v_row(j)
         vu_stored_rows(j,2*n_stored_vecs+2) = 0.5*tau(istep)*vav*v_row(j) - u_row(j)
       enddo
       do j=1,l_cols
         uv_stored_cols(j,2*n_stored_vecs+1) = 0.5*tau(istep)*vav*v_col(j) - u_col(j)
         uv_stored_cols(j,2*n_stored_vecs+2) = tau(istep)*v_col(j)
       enddo

       ! We have calculated another Hauseholder Vector, number of implicitly stored increased
       n_stored_vecs = n_stored_vecs+1

       ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
       if (n_stored_vecs == max_stored_uv .or. istep == 3) then

         if (useGPU .and. .not.(useIntelGPU)) then
           successGPU = gpu_memcpy(vu_stored_rows_dev, int(loc(vu_stored_rows(1,1)),kind=c_intptr_t), &
                                     max_local_rows * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_rows_dev", 1109,  successGPU)

           successGPU = gpu_memcpy(uv_stored_cols_dev, int(loc(uv_stored_cols(1,1)),kind=c_intptr_t), &
                                     max_local_cols * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_cols_dev", 1114,  successGPU)
         endif
          !if (useIntelGPU) then
          !  ! needed later when we can do explicit offloads
          !endif

         do i = 0, (istep-2)/tile_size
           ! go over tiles above (or on) the diagonal
           l_col_beg = i*l_cols_per_tile+1
           l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
           l_row_beg = 1
           l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
           if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) &
           cycle


           if (useGPU) then
             if (.not. mat_vec_as_one_block) then
               if (useIntelGPU) then
                  if (wantDebug) call obj%timer%start("mkl_offload")

                  call DGEMM('N', 'T',     &
                                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), &
                                       int(l_col_end-l_col_beg+1,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),     &
                                       ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv),                    &
                                       int(max_local_rows,kind=BLAS_KIND), uv_stored_cols(l_col_beg,1),                    &
                                       int(max_local_cols,kind=BLAS_KIND), ONE, a_mat(l_row_beg,l_col_beg),                &
                                       int(matrixRows,kind=BLAS_KIND))
                  if (wantDebug) call obj%timer%stop("mkl_offload")

               else
                 ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
                 ! we than use the same calls as for CPU version
                 if (wantDebug) call obj%timer%start("gpublas")
                 call gpublas_DGEMM('N', 'T',     &
                                         l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                         ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                         size_of_datatype,  &
                                         max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                         size_of_datatype,  &
                                         max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                         size_of_datatype , matrixRows)
                 if (wantDebug) call obj%timer%stop("gpublas")
               endif
             endif
           else !useGPU
             if (wantDebug) call obj%timer%start("blas")
             call DGEMM('N', 'T',                &
                                  int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                  int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                  ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                  int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                  uv_stored_cols(l_col_beg,1), &
                                  int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                  ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
             if (wantDebug) call obj%timer%stop("blas")
           endif !useGPU
         enddo

         if (useGPU) then
           if (mat_vec_as_one_block) then
             if (useIntelGPU) then
                if (wantDebug) call obj%timer%start("mkl_offload")
                call DGEMM('N', 'T', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                                    int(2*n_stored_vecs, kind=BLAS_KIND), ONE,  &
                                    vu_stored_rows, int(max_local_rows,kind=BLAS_KIND), &
                                    uv_stored_cols, int(max_local_cols,kind=BLAS_KIND),  &
                                    ONE, a_mat, int(matrixRows,kind=BLAS_KIND))
                if (wantDebug) call obj%timer%stop("mkl_offload")
             else
               !update whole (remaining) part of matrix, including tiles below diagonal
               !we can do that in one large cublas call
               if (wantDebug) call obj%timer%start("gpublas")
               call gpublas_DGEMM('N', 'T', l_rows, l_cols, 2*n_stored_vecs,   &
                                         ONE, vu_stored_rows_dev, max_local_rows, &
                                         uv_stored_cols_dev, max_local_cols,  &
                                         ONE, a_dev, matrixRows)
               if (wantDebug) call obj%timer%stop("gpublas")
             endif
           endif
         endif

         n_stored_vecs = 0
       endif

       if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then
         if (useGPU) then
           if (useIntelGPU) then
                       ! if (useIntelGPU) then
          ! needed at a later time when we can do explcit mem copys
          ! endif

           else
             !a_mat(l_rows,l_cols) = a_dev(l_rows,l_cols)
              a_offset = ((l_rows - 1) + matrixRows * (l_cols - 1)) * size_of_datatype

              successGPU = gpu_memcpy(int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), a_dev + a_offset, &
                                      1 *  size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: a_dev 3", 1239,  successGPU)
           endif
         endif
         if (n_stored_vecs > 0) then
           a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) &
                       + dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs),uv_stored_cols(l_cols,1:2*n_stored_vecs))
         end if
         if (isSkewsymmetric) then
           d_vec(istep-1) = 0.0_rk
         else
           d_vec(istep-1) = a_mat(l_rows,l_cols)
         endif

         if (useGPU) then
           if (useIntelGPU) then
          ! if (useIntelGPU) then
          ! needed at a later time when we can expicit mem copy
          ! endif
           else
             !a_dev(l_rows,l_cols) = a_mat(l_rows,l_cols)
             !successGPU = gpu_threadsynchronize()
             !call check_memcpy_GPU_f("tridiag: a_dev 4a5a", 1265,  successGPU)

             successGPU = gpu_memcpy(a_dev + a_offset, int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), &
                                     int(1 * size_of_datatype, kind=c_intptr_t), gpuMemcpyHostToDevice)
             call check_memcpy_GPU_f("tridiag: a_dev 4", 1269,  successGPU)
           endif
         endif
       endif

     enddo ! main cycle over istep=na,3,-1


  ! Store e_vec(1)

  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(2, nblk, np_cols)) then
    if (useGPU) then
      if (useIntelGPU) then
        e_vec(1) = a_mat(1,l_cols) ! use last l_cols value of loop above
      else
        successGPU = gpu_memcpy(int(loc(e_vec(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                                1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 7", 1364,  successGPU)
      endif
    else !useGPU
      e_vec(1) = a_mat(1,l_cols) ! use last l_cols value of loop above
    endif !useGPU
  endif

  ! Store d_vec(1)
  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols)) then
    if(useGPU) then
      if (useIntelGPU) then
        if (isSkewsymmetric) then
          d_vec(1) = 0.0_rk
        else
          d_vec(1) = a_mat(1,1)
        endif
      else
        successGPU = gpu_memcpy(int(loc(d_vec(1)),kind=c_intptr_t), a_dev, 1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 8", 1382,  successGPU)
      endif
    else !useGPU
      if (isSkewsymmetric) then
        d_vec(1) = 0.0_rk
      else
        d_vec(1) = a_mat(1,1)
      endif
    endif !useGPU
  endif

  deallocate(tmp, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp", 1395,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 1400,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 1403,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 1406,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 1409,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 1412,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 1415,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 1418,  successGPU)
  endif
  ! if (useIntelGPU) then
  ! needed at a later time when we can do explicit frees
  ! endif


  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 1428,  istat,  errorMessage)

  if (useNonBlockingCollectivesRows) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request4, mpierr)
    call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6,mpierr)
    call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif
  if (useNonBlockingCollectivesCols) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request5, mpierr)
    call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request7, mpierr)
    call mpi_wait(allreduce_request7, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif

  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 1477,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
           deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
     call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1482,  istat,  errorMessage)
    else
      successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: a_mat", 1485,  successGPU)

      successGPU = gpu_free_host(v_row_host)
      call check_host_dealloc_GPU_f("tridiag: v_row_host", 1488,  successGPU)
      nullify(v_row)

      successGPU = gpu_free_host(v_col_host)
      call check_host_dealloc_GPU_f("tridiag: v_col_host", 1492,  successGPU)
      nullify(v_col)

      successGPU = gpu_free_host(u_col_host)
      call check_host_dealloc_GPU_f("tridiag: u_col_host", 1496,  successGPU)
      nullify(u_col)

      successGPU = gpu_free_host(u_row_host)
      call check_host_dealloc_GPU_f("tridiag: u_row_host", 1500,  successGPU)
      nullify(u_row)

      successGPU = gpu_host_unregister(int(loc(uv_stored_cols),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: uv_stored_cols", 1504,  successGPU)

      successGPU = gpu_host_unregister(int(loc(vu_stored_rows),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: vu_stored_rows", 1507,  successGPU)

      successGPU = gpu_host_unregister(int(loc(e_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: e_vec", 1510,  successGPU)

      successGPU = gpu_host_unregister(int(loc(d_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: d_vec", 1513,  successGPU)
    endif
  else
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1517,  istat,  errorMessage)
  endif

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 1521,  istat,  errorMessage)

  call obj%timer%stop("tridiag_&
  &real&
  &" // &
  "_double" // &
  gpuString )

end subroutine tridiag_&
&real&
&_&
&double





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>

subroutine trans_ev_&
&real&
&_&
&double &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: rck = C_DOUBLE
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  real(kind=rck), intent(in)           :: tau(na)

  real(kind=rck), intent(inout)        :: a_mat(lda,*)
  real(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  real(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  real(kind=rck), pointer              :: tmp1(:), tmp2(:)
  real(kind=rck), allocatable          :: h1(:), h2(:)
  real(kind=rck), pointer              :: tmat(:,:)
  real(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &real
  integer(kind=ik)                              :: error
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk

  if (useIntelGPU) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp2", istat, errorMessage)
  endif

  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif


  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&real&
    !&", "hvm1", istat, errorMessage)
    num = (max_local_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(hvm1_host,num)
    call check_alloc_GPU_f("trans_ev: hvm1_host", 304,  successGPU)
    call c_f_pointer(hvm1_host,hvm1,(/(max_local_rows*max_stored_rows)/))

    num = (max_stored_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmat_host,num)
    call check_alloc_GPU_f("trans_ev: tmat_host", 309,  successGPU)
    call c_f_pointer(tmat_host,tmat,(/max_stored_rows,max_stored_rows/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_alloc_GPU_f("trans_ev: tmp1_host", 314,  successGPU)
    call c_f_pointer(tmp1_host,tmp1,(/(max_local_cols*max_stored_rows)/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp2_host,num)
    call check_alloc_GPU_f("trans_ev: tmp2_host", 319,  successGPU)
    call c_f_pointer(tmp2_host,tmp2,(/(max_local_cols*max_stored_rows)/))

    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 323,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 326,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 329,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 333,  successGPU)

    successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev: q_mat", 337,  successGPU)

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 341,  successGPU)
  endif  ! useGPU

  !if (useIntelGPU) then
  !  ! needed at a later time when we can do explicit copys
  !endif


  do istep = 1, na, blockStep
    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo

    if (nb>0) then
      if (useNonBlockingCollectivesCols) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(hvb, int(nb,kind=MPI_KIND), MPI_REAL8 , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
       else
        call obj%timer%start("mpi_communication")
        call mpi_bcast(hvb, int(nb,kind=MPI_KIND), MPI_REAL8 , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif
    endif

    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) &
      call DSYRK('U', 'T',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo
      if (nc>0) then
        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif
      endif
      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call DTRMV('L', 'T' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -h2(nc+1:nc+n)  &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU .and. .not.(useIntelGPU)) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 475,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 480,  successGPU)
      endif

      !if (useIntelGPU) then
      !  ! needed later when we can do explicit copys
      !endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            call obj%timer%start("mkl_offload")
            call obj%timer%stop("mkl_offload")

          else
            call obj%timer%start("gpublas")
            call gpublas_DGEMM('T', 'N',   &
                                     nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                     q_dev, ldq, ZERO, tmp_dev, nstor)
            call obj%timer%stop("gpublas")
          endif
        else ! useGPU

          call obj%timer%start("blas")
          call DGEMM('T', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (useIntelGPU) then
            tmp1(1:l_cols*nstor) = 0
          else
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 532,  successGPU)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0


      if (useGPU .and. .not.(useIntelGPU)) then
        ! In the legacy GPU version, this allreduce was ommited. But probably it has to be done for GPU + MPI
        ! todo: does it need to be copied whole? Wouldn't be a part sufficient?
        successGPU = gpu_memcpy(int(loc(tmp1(1)),kind=c_intptr_t), tmp_dev,  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("trans_ev", 547,  successGPU)
      endif

      !if (useIntelGPU) then
      !   ! needed later
      !endif

      if (useNonBlockingCollectivesRows) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_iallreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
        call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
      else
        call obj%timer%start("mpi_communication")
        call mpi_allreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif

      if (useGPU .and. .not.(useIntelGPU)) then
        ! copy back tmp2 - after reduction...
        successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2(1)),kind=c_intptr_t),  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
      endif ! useGPU


      !if (useIntelGPU) then
      !   ! needed later
      !endif


      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            ! tmp2 = tmat * tmp2
            call obj%timer%start("mkl_offload")

            call obj%timer%stop("mkl_offload")
          else
            call obj%timer%start("gpublas")
            call gpublas_DTRMM('L', 'L', 'N', 'N',     &
                                     nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                     tmp_dev, nstor)

            call gpublas_DGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                     -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                     ONE, q_dev, ldq)
            call obj%timer%stop("gpublas")
          endif
        else !useGPU
          ! tmp2 = tmat * tmp2
          call obj%timer%start("blas")
          call DTRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                             ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND))
          !q_mat = q_mat - hvm*tmp2
          call DGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(nstor,kind=BLAS_KIND),   &
                              -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND), &
                              ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

  enddo ! istep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 693,  istat,  errorMessage)

  if (useGPU) then

    if (useIntelGPU) then
      deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_&       &MATH_DATATYPE&       &: tmat, tmp1, tmp2", 701,  istat,  errorMessage)
    else

      !q_mat = q_dev
      successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                    q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev", 707,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev: q_mat", 710,  successGPU)

      successGPU = gpu_free_host(hvm1_host)
      call check_host_dealloc_GPU_f("trans_ev: hvm1_host", 713,  successGPU)
      nullify(hvm1)

      successGPU = gpu_free_host(tmat_host)
      call check_host_dealloc_GPU_f("trans_ev: tmat_host", 717,  successGPU)
      nullify(tmat)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp1_host", 721,  successGPU)
      nullify(tmp1)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp2_host", 725,  successGPU)
      nullify(tmp2)

      !deallocate(hvm1, stat=istat, errmsg=errorMessage)
      !if (istat .ne. 0) then
      !  print *,"trans_ev_&
      !  &real&
      !  &: error when deallocating hvm1 "//errorMessage
      !  stop 1
      !endif

      !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev", 738,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev", 741,  successGPU)

      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev", 744,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev", 747,  successGPU)
    endif
  else
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 753,  istat,  errorMessage)
  endif


  call obj%timer%stop("trans_ev_&
  &real&
  &" // &
  &"_double" // &
  gpuString )

end subroutine trans_ev_&
&real&
&_&
&double

! now comes a dirty hack:
! the file elpa1_solve_tridi_real_template.F90 must be included twice
! for the legacy and for the new API. In the new API, however, some routines
! must be named "..._impl"

!#include "elpa1_solve_tridi_real_template.F90"




!cannot use "../src/elpa1/../solve_tridi/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine solve_tridi_&
&double_impl &
    ( obj, na, nev, d, e, q, ldq, nblk, matrixCols, mpi_comm_all, mpi_comm_rows, &
                                           mpi_comm_cols, useGPU, wantDebug, success, max_threads )

      use precision
      use elpa_abstract_impl
      use merge_recursive
      use merge_systems
      use elpa_mpi
      use ELPA_utilities
      use distribute_global_column
      use elpa_mpi
      implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: rck = C_DOUBLE
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik), intent(in)               :: na, nev, ldq, nblk, matrixCols, &
                                                    mpi_comm_all, mpi_comm_rows, mpi_comm_cols
      real(kind=rk8), intent(inout)    :: d(na), e(na)
      real(kind=rk8), intent(inout)    :: q(ldq,*)
      logical, intent(in)                        :: useGPU, wantDebug
      logical, intent(out)                       :: success

      integer(kind=ik)                           :: i, j, n, np, nc, nev1, l_cols, l_rows
      integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik), allocatable              :: limits(:), l_col(:), p_col(:), l_col_bc(:), p_col_bc(:)

      integer(kind=ik)                           :: istat
      character(200)                             :: errorMessage
      character(20)                              :: gpuString
      integer(kind=ik), intent(in)               :: max_threads

      if(useGPU) then
        gpuString = "_gpu"
      else
        gpuString = ""
      endif

      call obj%timer%start("solve_tridi" // "_double" // gpuString)

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

      ! Set Q to 0
      q(1:l_rows, 1:l_cols) = 0.0_rk

      ! Get the limits of the subdivisons, each subdivison has as many cols
      ! as fit on the respective processor column

      allocate(limits(0:np_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: limits", 128,  istat,  errorMessage)

      limits(0) = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1) ! number of columns on proc column np

        ! Check for the case that a column has have zero width.
        ! This is not supported!
        ! Scalapack supports it but delivers no results for these columns,
        ! which is rather annoying
        if (nc==0) then
          call obj%timer%stop("solve_tridi" // "_double")
          if (wantDebug) write(error_unit,*) 'ELPA1_solve_tridi: ERROR: Problem contains processor column with zero width'
          success = .false.
          return
        endif
        limits(np+1) = limits(np) + nc
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,np_cols-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      ! Solve sub problems on processsor columns

      nc = limits(my_pcol) ! column after which my problem starts

      if (np_cols>1) then
        nev1 = l_cols ! all eigenvectors are needed
      else
        nev1 = MIN(nev,l_cols)
      endif
      call solve_tridi_col_&
           &double_impl &
             (obj, l_cols, nev1, nc, d(nc+1), e(nc+1), q, ldq, nblk,  &
                        matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads)
      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_double" // gpuString)
        return
      endif
      ! If there is only 1 processor column, we are done

      if (np_cols==1) then
        deallocate(limits, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi: limits", 176,  istat,  errorMessage)

        call obj%timer%stop("solve_tridi" // "_double" // gpuString)
        return
      endif

      ! Set index arrays for Q columns

      ! Dense distribution scheme:

      allocate(l_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col", 187,  istat,  errorMessage)

      allocate(p_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col", 190,  istat,  errorMessage)

      n = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1)
        do i=1,nc
          n = n+1
          l_col(n) = i
          p_col(n) = np
        enddo
      enddo

      ! Block cyclic distribution scheme, only nev columns are set:

      allocate(l_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col_bc", 205,  istat,  errorMessage)

      allocate(p_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col_bc", 208,  istat,  errorMessage)

      p_col_bc(:) = -1
      l_col_bc(:) = -1

      do i = 0, na-1, nblk*np_cols
        do j = 0, np_cols-1
          do n = 1, nblk
            if (i+j*nblk+n <= MIN(nev,na)) then
              p_col_bc(i+j*nblk+n) = j
              l_col_bc(i+j*nblk+n) = i/np_cols + n
             endif
           enddo
         enddo
      enddo

      ! Recursively merge sub problems
      call merge_recursive_&
           &double &
           (obj, 0, np_cols, ldq, matrixCols, nblk, &
           l_col, p_col, l_col_bc, p_col_bc, limits, &
           np_cols, na, q, d, e, &
           mpi_comm_all, mpi_comm_rows, mpi_comm_cols,&
           useGPU, wantDebug, success, max_threads)

      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_double" // gpuString)
        return
      endif

      deallocate(limits,l_col,p_col,l_col_bc,p_col_bc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi: limits, l_col, p_col, l_col_bc, p_col_bc", 239,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi" // "_double" // gpuString)
      return

    end subroutine solve_tridi_&
        &double_impl

    subroutine solve_tridi_col_&
    &double_impl &
      ( obj, na, nev, nqoff, d, e, q, ldq, nblk, matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads )

   ! Solves the symmetric, tridiagonal eigenvalue problem on one processor column
   ! with the divide and conquer method.
   ! Works best if the number of processor rows is a power of 2!
      use precision
      use elpa_abstract_impl
      use elpa_mpi
      use merge_systems
      use ELPA_utilities
      use distribute_global_column
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik)              :: na, nev, nqoff, ldq, nblk, matrixCols, mpi_comm_rows
      real(kind=rk8)      :: d(na), e(na)
      real(kind=rk8)      :: q(ldq,*)

      integer(kind=ik), parameter   :: min_submatrix_size = 16 ! Minimum size of the submatrices to be used

      real(kind=rk8), allocatable    :: qmat1(:,:), qmat2(:,:)
      integer(kind=ik)              :: i, n, np
      integer(kind=ik)              :: ndiv, noff, nmid, nlen, max_size
      integer(kind=ik)              :: my_prow, np_rows
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, np_rowsMPI

      integer(kind=ik), allocatable :: limits(:), l_col(:), p_col_i(:), p_col_o(:)
      logical, intent(in)           :: useGPU, wantDebug
      logical, intent(out)          :: success
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage

      integer(kind=ik), intent(in)  :: max_threads

      integer(kind=MPI_KIND)        :: bcast_request1, bcast_request2
      logical                       :: useNonBlockingCollectivesRows
      integer(kind=c_int)           :: non_blocking_collectives, error

      call obj%timer%start("solve_tridi_col" // "_double")

      call obj%get("nbc_row_solve_tridi", non_blocking_collectives, error)
      if (error .ne. ELPA_OK) then
        print *,"Problem setting option for non blocking collectives for rows in solve_tridi. Aborting..."
        stop
      endif

      if (non_blocking_collectives .eq. 1) then
        useNonBlockingCollectivesRows = .true.
      else
        useNonBlockingCollectivesRows = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.
      ! Calculate the number of subdivisions needed.

      n = na
      ndiv = 1
      do while(2*ndiv<=np_rows .and. n>2*min_submatrix_size)
        n = ((n+3)/4)*2 ! the bigger one of the two halves, we want EVEN boundaries
        ndiv = ndiv*2
      enddo

      ! If there is only 1 processor row and not all eigenvectors are needed
      ! and the matrix size is big enough, then use 2 subdivisions
      ! so that merge_systems is called once and only the needed
      ! eigenvectors are calculated for the final problem.

      if (np_rows==1 .and. nev<na .and. na>2*min_submatrix_size) ndiv = 2

      allocate(limits(0:ndiv), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits", 442,  istat,  errorMessage)

      limits(0) = 0
      limits(ndiv) = na

      n = ndiv
      do while(n>1)
        n = n/2 ! n is always a power of 2
        do i=0,ndiv-1,2*n
          ! We want to have even boundaries (for cache line alignments)
          limits(i+n) = limits(i) + ((limits(i+2*n)-limits(i)+3)/4)*2
        enddo
      enddo

      ! Calculate the maximum size of a subproblem

      max_size = 0
      do i=1,ndiv
        max_size = MAX(max_size,limits(i)-limits(i-1))
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,ndiv-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      if (np_rows==1)    then

        ! For 1 processor row there may be 1 or 2 subdivisions
        do n=0,ndiv-1
          noff = limits(n)        ! Start of subproblem
          nlen = limits(n+1)-noff ! Size of subproblem

          call solve_tridi_single_problem_&
          &double_impl &
                                  (obj, nlen,d(noff+1),e(noff+1), &
                                    q(nqoff+noff+1,noff+1),ubound(q,dim=1), wantDebug, success)

          if (.not.(success)) return
        enddo

      else

        ! Solve sub problems in parallel with solve_tridi_single
        ! There is at maximum 1 subproblem per processor

        allocate(qmat1(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1", 492,  istat,  errorMessage)

        allocate(qmat2(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat2", 495,  istat,  errorMessage)

        qmat1 = 0 ! Make sure that all elements are defined

        if (my_prow < ndiv) then

          noff = limits(my_prow)        ! Start of subproblem
          nlen = limits(my_prow+1)-noff ! Size of subproblem
          call solve_tridi_single_problem_&
          &double_impl &
                                    (obj, nlen,d(noff+1),e(noff+1),qmat1, &
                                    ubound(qmat1,dim=1), wantDebug, success)

          if (.not.(success)) return
        endif

        ! Fill eigenvectors in qmat1 into global matrix q

        do np = 0, ndiv-1

          noff = limits(np)
          nlen = limits(np+1)-noff
          if (useNonBlockingCollectivesRows) then
            call obj%timer%start("mpi_nbc_communication")
            call mpi_ibcast(d(noff+1), int(nlen,kind=MPI_KIND), MPI_REAL8, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), bcast_request1, mpierr)
            call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)

            qmat2 = qmat1
            call mpi_ibcast(qmat2, int(max_size*max_size,kind=MPI_KIND), MPI_REAL8, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), bcast_request2, mpierr)
            call mpi_wait(bcast_request2, MPI_STATUS_IGNORE, mpierr)
            call obj%timer%stop("mpi_nbc_communication")
          else
            call obj%timer%start("mpi_communication")
            call mpi_bcast(d(noff+1), int(nlen,kind=MPI_KIND), MPI_REAL8, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)

            qmat2 = qmat1
            call mpi_bcast(qmat2, int(max_size*max_size,kind=MPI_KIND), MPI_REAL8, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
          endif
          do i=1,nlen

            call distribute_global_column_&
            &double &
                     (obj, qmat2(1,i), q(1,noff+i), nqoff+noff, nlen, my_prow, np_rows, nblk)
          enddo

        enddo

        deallocate(qmat1, qmat2, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1, qmat2", 558,  istat,  errorMessage)

      endif

      ! Allocate and set index arrays l_col and p_col

      allocate(l_col(na), p_col_i(na),  p_col_o(na), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: l_col, p_col_i, p_col_o", 565,  istat,  errorMessage)

      do i=1,na
        l_col(i) = i
        p_col_i(i) = 0
        p_col_o(i) = 0
      enddo

      ! Merge subproblems

      n = 1
      do while(n<ndiv) ! if ndiv==1, the problem was solved by single call to solve_tridi_single

        do i=0,ndiv-1,2*n

          noff = limits(i)
          nmid = limits(i+n) - noff
          nlen = limits(i+2*n) - noff

          if (nlen == na) then
            ! Last merge, set p_col_o=-1 for unneeded (output) eigenvectors
            p_col_o(nev+1:na) = -1
          endif
          call merge_systems_&
          &double &
                              (obj, nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, nqoff+noff, nblk, &
                               matrixCols, int(mpi_comm_rows,kind=ik), int(mpi_comm_self,kind=ik), &
                               l_col(noff+1), p_col_i(noff+1), &
                               l_col(noff+1), p_col_o(noff+1), 0, 1, useGPU, wantDebug, success, max_threads)
          if (.not.(success)) return

        enddo

        n = 2*n

      enddo

      deallocate(limits, l_col, p_col_i, p_col_o, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits, l_col, p_col_i, p_col_o", 603,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi_col" // "_double")

    end subroutine solve_tridi_col_&
    &double_impl

    subroutine solve_tridi_single_problem_&
    &double_impl &
    (obj, nlen, d, e, q, ldq, wantDebug, success)

   ! Solves the symmetric, tridiagonal eigenvalue problem on a single processor.
   ! Takes precautions if DSTEDC fails or if the eigenvalues are not ordered correctly.
     use precision
     use elpa_abstract_impl
     use elpa_blas_interfaces
     use ELPA_utilities
     implicit none
     class(elpa_abstract_impl_t), intent(inout) :: obj
     integer(kind=ik)                         :: nlen, ldq
     real(kind=rk8)                 :: d(nlen), e(nlen), q(ldq,nlen)

     real(kind=rk8), allocatable    :: work(:), qtmp(:), ds(:), es(:)
     real(kind=rk8)                 :: dtmp

     integer(kind=ik)              :: i, j, lwork, liwork, info
     integer(kind=BLAS_KIND)       :: infoBLAS
     integer(kind=ik), allocatable :: iwork(:)

     logical, intent(in)           :: wantDebug
     logical, intent(out)          :: success
      integer(kind=ik)             :: istat
      character(200)               :: errorMessage

     call obj%timer%start("solve_tridi_single" // "_double")

     success = .true.
     allocate(ds(nlen), es(nlen), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: ds, es", 641,  istat,  errorMessage)

     ! Save d and e for the case that dstedc fails

     ds(:) = d(:)
     es(:) = e(:)

     ! First try dstedc, this is normally faster but it may fail sometimes (why???)

     lwork = 1 + 4*nlen + nlen**2
     liwork =  3 + 5*nlen
     allocate(work(lwork), iwork(liwork), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: work, iwork", 653,  istat,  errorMessage)
     call obj%timer%start("blas")
     call DSTEDC('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND),    &
                          work, int(lwork,kind=BLAS_KIND), int(iwork,kind=BLAS_KIND), int(liwork,kind=BLAS_KIND), &
                          infoBLAS)
     info = int(infoBLAS,kind=ik)
     call obj%timer%stop("blas")

     if (info /= 0) then

       ! DSTEDC failed, try DSTEQR. The workspace is enough for DSTEQR.

       write(error_unit,'(a,i8,a)') 'Warning: Lapack routine DSTEDC failed, info= ',info,', Trying DSTEQR!'

       d(:) = ds(:)
       e(:) = es(:)
       call obj%timer%start("blas")
       call DSTEQR('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND), work, infoBLAS )
       info = int(infoBLAS,kind=ik)
       call obj%timer%stop("blas")

       ! If DSTEQR fails also, we don't know what to do further ...

       if (info /= 0) then
         if (wantDebug) &
           write(error_unit,'(a,i8,a)') 'ELPA1_solve_tridi_single: ERROR: Lapack routine DSTEQR failed, info= ',info,', Aborting!'
           success = .false.
           return
         endif
       end if

       deallocate(work,iwork,ds,es, stat=istat, errmsg=errorMessage)
       call check_deallocate_f("solve_tridi_single: work, iwork, ds, es", 685,  istat,  errorMessage)

      ! Check if eigenvalues are monotonically increasing
      ! This seems to be not always the case  (in the IBM implementation of dstedc ???)

      do i=1,nlen-1
        if (d(i+1)<d(i)) then
          if (abs(d(i+1) - d(i)) / abs(d(i+1) + d(i)) > 1e-14_rk8) then
            write(error_unit,'(a,i8,2g25.16)') '***WARNING: Monotony error dste**:',i+1,d(i),d(i+1)
          else
            write(error_unit,'(a,i8,2g25.16)') 'Info: Monotony error dste{dc,qr}:',i+1,d(i),d(i+1)
            write(error_unit,'(a)') 'The eigenvalues from a lapack call are not sorted to machine precision.'
            write(error_unit,'(a)') 'In this extent, this is completely harmless.'
            write(error_unit,'(a)') 'Still, we keep this info message just in case.'
          end if
          allocate(qtmp(nlen), stat=istat, errmsg=errorMessage)
          call check_allocate_f("solve_tridi_single: qtmp", 705,  istat,  errorMessage)

          dtmp = d(i+1)
          qtmp(1:nlen) = q(1:nlen,i+1)
          do j=i,1,-1
            if (dtmp<d(j)) then
              d(j+1)        = d(j)
              q(1:nlen,j+1) = q(1:nlen,j)
            else
              exit ! Loop
            endif
          enddo
          d(j+1)        = dtmp
          q(1:nlen,j+1) = qtmp(1:nlen)
          deallocate(qtmp, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("solve_tridi_single: qtmp", 720,  istat,  errorMessage)

       endif
     enddo
     call obj%timer%stop("solve_tridi_single" // "_double")

    end subroutine solve_tridi_single_problem_&
    &double_impl

!#include "elpa1_merge_systems_real_template.F90"










subroutine hh_transform_real_&
&double &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine DLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
  implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: rck = C_DOUBLE
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  real(kind=rk), intent(inout)       :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  real(kind=rk), intent(out)         :: xf, tau

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &real&
     	      &" // &
                   &"_double" )


  if ( XNORM_SQ==0.0_rk ) then

    if ( ALPHA>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHA**2 + XNORM_SQ ), ALPHA )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHA = XNORM_SQ / ALPHA

      TAU = ALPHA / BETA
      ALPHA = -ALPHA
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &real&
  &" // &
  &"_double" )

end subroutine hh_transform_real_&
    &double




! real single precision

























!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)









!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &real&
  &_&
  &single &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads, isSkewsymmetric)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_FLOAT
  integer, parameter :: rck = C_FLOAT
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  real(kind=rck), intent(out)          :: tau(na)
  real(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik), parameter                   :: max_stored_uv = 32
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  !integer(kind=ik)                             :: l_cols, l_rows
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs

  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: a_offset

  integer(kind=ik), intent(in)                  :: max_threads

  real(kind=rk)                                 :: vnorm2
  real(kind=rck)                       :: vav, x, aux(2*max_stored_uv), aux1(2), aux2(2), vrl, xf

  integer(kind=c_intptr_t)                      :: num
  real(kind=rck), allocatable          :: tmp(:)
  real(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  real(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  real(kind=rck), pointer             :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  real(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &real
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_&
  &real&
  &" // &
  "_single" // &
  gpuString )

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for min_tile_size. Aborting..."
    stop
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows
  !
  allocate(tmp(MAX(max_local_rows,max_local_cols)), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
  &real ", "tmp", istat, errorMessage)

  ! allocate v_row 1 element longer to allow store and broadcast tau together with it
  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &real ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &real ", "v_row", istat, errorMessage)
  
      allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
       &real ", "v_col", istat, errorMessage)
  
      allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &real ", "u_col", istat, errorMessage)
  
      allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &real ", "u_row", istat, errorMessage)
    else

      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 346,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 351,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 356,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))

      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 361,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))

      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_roes", 367,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 372,  successGPU)

      num = na * 4
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 381,  successGPU)

      num = na * 4
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 390,  successGPU)
    endif
  else
    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &real ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &real ", "u_row", istat, errorMessage)
      
  endif


  tmp = 0
  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU .and. .not.(useIntelGPU) ) then
     successGPU = gpu_malloc(v_row_dev, max_local_rows * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_row_dev", 429,  successGPU)

     successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

     call check_alloc_GPU_f("tridiag: u_row_dev", 433,  successGPU)

     successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_col_dev", 436,  successGPU)

     successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: u_col_dev", 439,  successGPU)

     successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 442,  successGPU)

     successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 445,  successGPU)
  endif !useGPU

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) &
  d_vec(na) = a_mat(l_rows,l_cols)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! allocate memmory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 476,  successGPU)

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("tridiag: a_mat", 480,  successGPU)

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 484,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1

    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

      ! copy l_cols + 1 column of A to v_row
      if (useGPU) then
        if (useIntelGPU) then
          v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
        else
          a_offset = l_cols * matrixRows * size_of_datatype
          ! we use v_row on the host at the moment! successGPU = gpu_memcpy(v_row_dev, a_dev + a_offset, 
          ! (l_rows)*size_of_PRECISION_real, gpuMemcpyDeviceToDevice)

          successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), &
                                    a_dev + a_offset, (l_rows)* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag a_dev 1", 519,  successGPU)
        endif
      else
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        if (useIntelGPU) then
                !print *,"intel phase aaaaaaaaaaaaaaaaaaaaaaaaaa"
          if (wantDebug) call obj%timer%start("mkl_offload")
          aux(1:2*n_stored_vecs) = uv_stored_cols(l_cols+1,1:2*n_stored_vecs)



          if (wantDebug) call obj%timer%stop("mkl_offload")
        else
          if (wantDebug) call obj%timer%start("blas")
          call SGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            uv_stored_cols(l_cols+1,1), &
                            int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND), &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif
      endif

      if (my_prow == prow(istep-1, nblk, np_rows)) then
        aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1))
        aux1(2) = v_row(l_rows)
      else
        aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
        aux1(2) = 0.
      endif

      if (useNonBlockingCollectivesRows) then
        if (wantDebug) call obj%timer%start("mpi_communication_non_blocking")
        call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL4, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
        call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication_non_blocking")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL4, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif

      vnorm2 = aux2(1)
      vrl    = aux2(2)

      ! Householder transformation
      call hh_transform_real_&
               &single &
               (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      ! Scale v_row and store Householder Vector for back transformation

      v_row(1:l_rows) = v_row(1:l_rows) * xf
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        v_row(l_rows) = 1.

        ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
        e_vec(istep-1) = vrl
      endif

      ! store Householder Vector for back transformation
      a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)

      ! add tau after the end of actuall v_row, to be broadcasted with it
      v_row(l_rows+1) = tau(istep)
    endif !(my_pcol == pcol(istep, nblk, np_cols))

!          

    if (useNonBlockingCollectivesCols) then
      if (wantDebug) call obj%timer%start("mpi_nbc_communication")
      ! Broadcast the Householder Vector (and tau) along columns
      call mpi_ibcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_REAL4,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   bcast_request1, mpierr)
      call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
      if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
    else
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_bcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_REAL4,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")
    endif

    !recover tau, which has been broadcasted together with v_row
    tau(istep) =  v_row(l_rows+1)

    ! Transpose Householder Vector v_row -> v_col
    call elpa_transpose_vectors_&
        &real&
        &_&
        &single &
              (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
               1, istep-1, 1, nblk, max_threads, .true.)

    ! Calculate u = (A + VU**T + UV**T)*v

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    u_col(1:l_cols) = 0
    u_row(1:l_rows) = 0
    if (l_rows > 0 .and. l_cols> 0 ) then
     if (useGPU .and. .not.(useIntelGPU)) then
       successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_col_dev", 693,  successGPU)

       successGPU = gpu_memset(u_row_dev, 0, l_rows * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_row_dev", 696,  successGPU)

       successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                     l_cols * size_of_datatype, gpuMemcpyHostToDevice)

       call check_memcpy_GPU_f("tridiag: v_col_dev", 701,  successGPU)

       successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                 l_rows * size_of_datatype, gpuMemcpyHostToDevice)
       call check_memcpy_GPU_f("tridiag: v_row_dev", 705,  successGPU)
     endif ! useGPU

     !if (useIntelGPU) then
     !  ! needed later when we can do explicit memcopy
     !endif

     do i= 0, (istep-2)/tile_size
       l_col_beg = i*l_cols_per_tile+1
       l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
       if (l_col_end < l_col_beg) cycle
       do j = 0, i
         l_row_beg = j*l_rows_per_tile+1
         l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
         if (l_row_end < l_row_beg) cycle

         ! multiplication by blocks is efficient only for CPU
         ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
         ! CPU implementation) or by one large matrix Vector multiply
         if (.not. useGPU) then
           if (wantDebug) call obj%timer%start("blas")
           call SGEMV('T',  &
                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                       ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                       v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                       ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

           if (i/=j) then
             if (isSkewsymmetric) then
               call SGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                   -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)

             else
               call SGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                   ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)
             endif
           endif
           if (wantDebug) call obj%timer%stop("blas")
         endif ! not useGPU

            enddo  ! j=0,i
          enddo  ! i=0,(istep-2)/tile_size

          if (useGPU) then
            if (mat_vec_as_one_block) then
              if (useIntelGPU) then
                 if (wantDebug) call obj%timer%start("mkl_offload")

                if (wantDebug) call obj%timer%stop("mkl_offload")

              else
                ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
                ! this requires altering of the algorithm when later explicitly updating the matrix
                ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
                if (wantDebug) call obj%timer%start("gpublas")
                call gpublas_SGEMV('T', l_rows,l_cols,  &
                                          ONE, a_dev, matrixRows,                   &
                                          v_row_dev , 1,                          &
                                          ONE, u_col_dev, 1)

       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_SGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + a_offset, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
                if (wantDebug) call obj%timer%stop("gpublas")
              endif
            else  ! mat_vec_as_one_block
              !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
              !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")

                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype

                    call gpublas_SGEMV('T', &
                                l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                ONE, a_dev + a_offset, matrixRows,  &
                                v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                                ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1)
                endif
              enddo

              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,i*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")


                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype
                    if (isSkewsymmetric) then
                       call gpublas_SGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   -ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                    else
                       call gpublas_SGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                   endif
                endif
              enddo
            end if !multiplication as one block / per stripes

            if (.not.(useIntelGPU)) then
              successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), &
                          u_col_dev, l_cols * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_col_dev 1", 954,  successGPU)

              successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), &
                          u_row_dev, l_rows * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_row_dev 1", 958,  successGPU)
            endif
            !if (useIntelGPU) then
            !  
            !endif

          endif ! useGPU


         ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
         if (n_stored_vecs > 0) then
           if (wantDebug) call obj%timer%start("blas")
           call SGEMV('T',     &
                               int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                               v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)

           call SGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                               aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
           if (wantDebug) call obj%timer%stop("blas")
         endif

       endif  ! (l_rows>0 .and. l_cols>0)

       ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
       ! on the processors containing the diagonal
       ! This is only necessary if u_row has been calculated, i.e. if the
       ! global tile size is smaller than the global remaining matrix

       if (tile_size < istep-1) then

         call elpa_reduce_add_vectors_&
         &real&
         &_&
         &single &
         (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
         mpi_comm_cols, istep-1, 1, nblk, max_threads)

       endif

       ! Sum up all the u_col(:) parts, transpose u_col -> u_row

       if (l_cols>0) then
         tmp(1:l_cols) = u_col(1:l_cols)
         if (useNonBlockingCollectivesRows) then
           if (wantDebug) call obj%timer%start("mpi_nbc_communication")
           call mpi_iallreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_REAL4,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
           call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
         else
           if (wantDebug) call obj%timer%start("mpi_communication")
           call mpi_allreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_REAL4,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
           if (wantDebug) call obj%timer%stop("mpi_communication")
         endif
       endif
       if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
          &real&
          &_&
          &single &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       else
          call elpa_transpose_vectors_&
          &real&
          &_&
          &single &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       endif

       ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
       x = 0
       if (l_cols>0)  &
       x = dot_product(v_col(1:l_cols),u_col(1:l_cols))

       if (useNonBlockingCollectivesCols) then
         if (wantDebug) call obj%timer%start("mpi_nbc_communication")
         call mpi_iallreduce(x, vav, 1_MPI_KIND, MPI_REAL4, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
               allreduce_request3, mpierr)
         call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
       else
         if (wantDebug) call obj%timer%start("mpi_communication")
         call mpi_allreduce(x, vav, 1_MPI_KIND, MPI_REAL4, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
                mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")
       endif

       ! store u and v in the matrices U and V
       ! these matrices are stored combined in one here

       do j=1,l_rows
         vu_stored_rows(j,2*n_stored_vecs+1) = tau(istep)*v_row(j)
         vu_stored_rows(j,2*n_stored_vecs+2) = 0.5*tau(istep)*vav*v_row(j) - u_row(j)
       enddo
       do j=1,l_cols
         uv_stored_cols(j,2*n_stored_vecs+1) = 0.5*tau(istep)*vav*v_col(j) - u_col(j)
         uv_stored_cols(j,2*n_stored_vecs+2) = tau(istep)*v_col(j)
       enddo

       ! We have calculated another Hauseholder Vector, number of implicitly stored increased
       n_stored_vecs = n_stored_vecs+1

       ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
       if (n_stored_vecs == max_stored_uv .or. istep == 3) then

         if (useGPU .and. .not.(useIntelGPU)) then
           successGPU = gpu_memcpy(vu_stored_rows_dev, int(loc(vu_stored_rows(1,1)),kind=c_intptr_t), &
                                     max_local_rows * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_rows_dev", 1109,  successGPU)

           successGPU = gpu_memcpy(uv_stored_cols_dev, int(loc(uv_stored_cols(1,1)),kind=c_intptr_t), &
                                     max_local_cols * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_cols_dev", 1114,  successGPU)
         endif
          !if (useIntelGPU) then
          !  ! needed later when we can do explicit offloads
          !endif

         do i = 0, (istep-2)/tile_size
           ! go over tiles above (or on) the diagonal
           l_col_beg = i*l_cols_per_tile+1
           l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
           l_row_beg = 1
           l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
           if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) &
           cycle


           if (useGPU) then
             if (.not. mat_vec_as_one_block) then
               if (useIntelGPU) then
                  if (wantDebug) call obj%timer%start("mkl_offload")

                  call SGEMM('N', 'T',     &
                                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), &
                                       int(l_col_end-l_col_beg+1,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),     &
                                       ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv),                    &
                                       int(max_local_rows,kind=BLAS_KIND), uv_stored_cols(l_col_beg,1),                    &
                                       int(max_local_cols,kind=BLAS_KIND), ONE, a_mat(l_row_beg,l_col_beg),                &
                                       int(matrixRows,kind=BLAS_KIND))
                  if (wantDebug) call obj%timer%stop("mkl_offload")

               else
                 ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
                 ! we than use the same calls as for CPU version
                 if (wantDebug) call obj%timer%start("gpublas")
                 call gpublas_SGEMM('N', 'T',     &
                                         l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                         ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                         size_of_datatype,  &
                                         max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                         size_of_datatype,  &
                                         max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                         size_of_datatype , matrixRows)
                 if (wantDebug) call obj%timer%stop("gpublas")
               endif
             endif
           else !useGPU
             if (wantDebug) call obj%timer%start("blas")
             call SGEMM('N', 'T',                &
                                  int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                  int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                  ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                  int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                  uv_stored_cols(l_col_beg,1), &
                                  int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                  ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
             if (wantDebug) call obj%timer%stop("blas")
           endif !useGPU
         enddo

         if (useGPU) then
           if (mat_vec_as_one_block) then
             if (useIntelGPU) then
                if (wantDebug) call obj%timer%start("mkl_offload")
                call SGEMM('N', 'T', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                                    int(2*n_stored_vecs, kind=BLAS_KIND), ONE,  &
                                    vu_stored_rows, int(max_local_rows,kind=BLAS_KIND), &
                                    uv_stored_cols, int(max_local_cols,kind=BLAS_KIND),  &
                                    ONE, a_mat, int(matrixRows,kind=BLAS_KIND))
                if (wantDebug) call obj%timer%stop("mkl_offload")
             else
               !update whole (remaining) part of matrix, including tiles below diagonal
               !we can do that in one large cublas call
               if (wantDebug) call obj%timer%start("gpublas")
               call gpublas_SGEMM('N', 'T', l_rows, l_cols, 2*n_stored_vecs,   &
                                         ONE, vu_stored_rows_dev, max_local_rows, &
                                         uv_stored_cols_dev, max_local_cols,  &
                                         ONE, a_dev, matrixRows)
               if (wantDebug) call obj%timer%stop("gpublas")
             endif
           endif
         endif

         n_stored_vecs = 0
       endif

       if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then
         if (useGPU) then
           if (useIntelGPU) then
                       ! if (useIntelGPU) then
          ! needed at a later time when we can do explcit mem copys
          ! endif

           else
             !a_mat(l_rows,l_cols) = a_dev(l_rows,l_cols)
              a_offset = ((l_rows - 1) + matrixRows * (l_cols - 1)) * size_of_datatype

              successGPU = gpu_memcpy(int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), a_dev + a_offset, &
                                      1 *  size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: a_dev 3", 1239,  successGPU)
           endif
         endif
         if (n_stored_vecs > 0) then
           a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) &
                       + dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs),uv_stored_cols(l_cols,1:2*n_stored_vecs))
         end if
         if (isSkewsymmetric) then
           d_vec(istep-1) = 0.0_rk
         else
           d_vec(istep-1) = a_mat(l_rows,l_cols)
         endif

         if (useGPU) then
           if (useIntelGPU) then
          ! if (useIntelGPU) then
          ! needed at a later time when we can expicit mem copy
          ! endif
           else
             !a_dev(l_rows,l_cols) = a_mat(l_rows,l_cols)
             !successGPU = gpu_threadsynchronize()
             !call check_memcpy_GPU_f("tridiag: a_dev 4a5a", 1265,  successGPU)

             successGPU = gpu_memcpy(a_dev + a_offset, int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), &
                                     int(1 * size_of_datatype, kind=c_intptr_t), gpuMemcpyHostToDevice)
             call check_memcpy_GPU_f("tridiag: a_dev 4", 1269,  successGPU)
           endif
         endif
       endif

     enddo ! main cycle over istep=na,3,-1


  ! Store e_vec(1)

  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(2, nblk, np_cols)) then
    if (useGPU) then
      if (useIntelGPU) then
        e_vec(1) = a_mat(1,l_cols) ! use last l_cols value of loop above
      else
        successGPU = gpu_memcpy(int(loc(e_vec(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                                1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 7", 1364,  successGPU)
      endif
    else !useGPU
      e_vec(1) = a_mat(1,l_cols) ! use last l_cols value of loop above
    endif !useGPU
  endif

  ! Store d_vec(1)
  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols)) then
    if(useGPU) then
      if (useIntelGPU) then
        if (isSkewsymmetric) then
          d_vec(1) = 0.0_rk
        else
          d_vec(1) = a_mat(1,1)
        endif
      else
        successGPU = gpu_memcpy(int(loc(d_vec(1)),kind=c_intptr_t), a_dev, 1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 8", 1382,  successGPU)
      endif
    else !useGPU
      if (isSkewsymmetric) then
        d_vec(1) = 0.0_rk
      else
        d_vec(1) = a_mat(1,1)
      endif
    endif !useGPU
  endif

  deallocate(tmp, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp", 1395,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 1400,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 1403,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 1406,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 1409,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 1412,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 1415,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 1418,  successGPU)
  endif
  ! if (useIntelGPU) then
  ! needed at a later time when we can do explicit frees
  ! endif


  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 1428,  istat,  errorMessage)

  if (useNonBlockingCollectivesRows) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request4, mpierr)
    call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6,mpierr)
    call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif
  if (useNonBlockingCollectivesCols) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request5, mpierr)
    call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request7, mpierr)
    call mpi_wait(allreduce_request7, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif

  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 1477,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
           deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
     call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1482,  istat,  errorMessage)
    else
      successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: a_mat", 1485,  successGPU)

      successGPU = gpu_free_host(v_row_host)
      call check_host_dealloc_GPU_f("tridiag: v_row_host", 1488,  successGPU)
      nullify(v_row)

      successGPU = gpu_free_host(v_col_host)
      call check_host_dealloc_GPU_f("tridiag: v_col_host", 1492,  successGPU)
      nullify(v_col)

      successGPU = gpu_free_host(u_col_host)
      call check_host_dealloc_GPU_f("tridiag: u_col_host", 1496,  successGPU)
      nullify(u_col)

      successGPU = gpu_free_host(u_row_host)
      call check_host_dealloc_GPU_f("tridiag: u_row_host", 1500,  successGPU)
      nullify(u_row)

      successGPU = gpu_host_unregister(int(loc(uv_stored_cols),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: uv_stored_cols", 1504,  successGPU)

      successGPU = gpu_host_unregister(int(loc(vu_stored_rows),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: vu_stored_rows", 1507,  successGPU)

      successGPU = gpu_host_unregister(int(loc(e_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: e_vec", 1510,  successGPU)

      successGPU = gpu_host_unregister(int(loc(d_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: d_vec", 1513,  successGPU)
    endif
  else
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1517,  istat,  errorMessage)
  endif

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 1521,  istat,  errorMessage)

  call obj%timer%stop("tridiag_&
  &real&
  &" // &
  "_single" // &
  gpuString )

end subroutine tridiag_&
&real&
&_&
&single





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>

subroutine trans_ev_&
&real&
&_&
&single &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_FLOAT
  integer, parameter :: rck = C_FLOAT
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  real(kind=rck), intent(in)           :: tau(na)

  real(kind=rck), intent(inout)        :: a_mat(lda,*)
  real(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  real(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  real(kind=rck), pointer              :: tmp1(:), tmp2(:)
  real(kind=rck), allocatable          :: h1(:), h2(:)
  real(kind=rck), pointer              :: tmat(:,:)
  real(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &real
  integer(kind=ik)                              :: error
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk

  if (useIntelGPU) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp2", istat, errorMessage)
  endif

  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &real&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &real&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif


  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&real&
    !&", "hvm1", istat, errorMessage)
    num = (max_local_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(hvm1_host,num)
    call check_alloc_GPU_f("trans_ev: hvm1_host", 304,  successGPU)
    call c_f_pointer(hvm1_host,hvm1,(/(max_local_rows*max_stored_rows)/))

    num = (max_stored_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmat_host,num)
    call check_alloc_GPU_f("trans_ev: tmat_host", 309,  successGPU)
    call c_f_pointer(tmat_host,tmat,(/max_stored_rows,max_stored_rows/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_alloc_GPU_f("trans_ev: tmp1_host", 314,  successGPU)
    call c_f_pointer(tmp1_host,tmp1,(/(max_local_cols*max_stored_rows)/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp2_host,num)
    call check_alloc_GPU_f("trans_ev: tmp2_host", 319,  successGPU)
    call c_f_pointer(tmp2_host,tmp2,(/(max_local_cols*max_stored_rows)/))

    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 323,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 326,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 329,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 333,  successGPU)

    successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev: q_mat", 337,  successGPU)

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 341,  successGPU)
  endif  ! useGPU

  !if (useIntelGPU) then
  !  ! needed at a later time when we can do explicit copys
  !endif


  do istep = 1, na, blockStep
    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo

    if (nb>0) then
      if (useNonBlockingCollectivesCols) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(hvb, int(nb,kind=MPI_KIND), MPI_REAL4 , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
       else
        call obj%timer%start("mpi_communication")
        call mpi_bcast(hvb, int(nb,kind=MPI_KIND), MPI_REAL4 , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif
    endif

    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) &
      call SSYRK('U', 'T',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo
      if (nc>0) then
        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif
      endif
      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call STRMV('L', 'T' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -h2(nc+1:nc+n)  &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU .and. .not.(useIntelGPU)) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 475,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 480,  successGPU)
      endif

      !if (useIntelGPU) then
      !  ! needed later when we can do explicit copys
      !endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            call obj%timer%start("mkl_offload")
            call obj%timer%stop("mkl_offload")

          else
            call obj%timer%start("gpublas")
            call gpublas_SGEMM('T', 'N',   &
                                     nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                     q_dev, ldq, ZERO, tmp_dev, nstor)
            call obj%timer%stop("gpublas")
          endif
        else ! useGPU

          call obj%timer%start("blas")
          call SGEMM('T', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (useIntelGPU) then
            tmp1(1:l_cols*nstor) = 0
          else
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 532,  successGPU)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0


      if (useGPU .and. .not.(useIntelGPU)) then
        ! In the legacy GPU version, this allreduce was ommited. But probably it has to be done for GPU + MPI
        ! todo: does it need to be copied whole? Wouldn't be a part sufficient?
        successGPU = gpu_memcpy(int(loc(tmp1(1)),kind=c_intptr_t), tmp_dev,  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("trans_ev", 547,  successGPU)
      endif

      !if (useIntelGPU) then
      !   ! needed later
      !endif

      if (useNonBlockingCollectivesRows) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_iallreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
        call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
      else
        call obj%timer%start("mpi_communication")
        call mpi_allreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif

      if (useGPU .and. .not.(useIntelGPU)) then
        ! copy back tmp2 - after reduction...
        successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2(1)),kind=c_intptr_t),  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
      endif ! useGPU


      !if (useIntelGPU) then
      !   ! needed later
      !endif


      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            ! tmp2 = tmat * tmp2
            call obj%timer%start("mkl_offload")

            call obj%timer%stop("mkl_offload")
          else
            call obj%timer%start("gpublas")
            call gpublas_STRMM('L', 'L', 'N', 'N',     &
                                     nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                     tmp_dev, nstor)

            call gpublas_SGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                     -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                     ONE, q_dev, ldq)
            call obj%timer%stop("gpublas")
          endif
        else !useGPU
          ! tmp2 = tmat * tmp2
          call obj%timer%start("blas")
          call STRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                             ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND))
          !q_mat = q_mat - hvm*tmp2
          call SGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(nstor,kind=BLAS_KIND),   &
                              -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND), &
                              ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

  enddo ! istep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 693,  istat,  errorMessage)

  if (useGPU) then

    if (useIntelGPU) then
      deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_&       &MATH_DATATYPE&       &: tmat, tmp1, tmp2", 701,  istat,  errorMessage)
    else

      !q_mat = q_dev
      successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                    q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev", 707,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev: q_mat", 710,  successGPU)

      successGPU = gpu_free_host(hvm1_host)
      call check_host_dealloc_GPU_f("trans_ev: hvm1_host", 713,  successGPU)
      nullify(hvm1)

      successGPU = gpu_free_host(tmat_host)
      call check_host_dealloc_GPU_f("trans_ev: tmat_host", 717,  successGPU)
      nullify(tmat)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp1_host", 721,  successGPU)
      nullify(tmp1)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp2_host", 725,  successGPU)
      nullify(tmp2)

      !deallocate(hvm1, stat=istat, errmsg=errorMessage)
      !if (istat .ne. 0) then
      !  print *,"trans_ev_&
      !  &real&
      !  &: error when deallocating hvm1 "//errorMessage
      !  stop 1
      !endif

      !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev", 738,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev", 741,  successGPU)

      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev", 744,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev", 747,  successGPU)
    endif
  else
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 753,  istat,  errorMessage)
  endif


  call obj%timer%stop("trans_ev_&
  &real&
  &" // &
  &"_single" // &
  gpuString )

end subroutine trans_ev_&
&real&
&_&
&single

! now comes a dirty hack:
! the file elpa1_solve_tridi_real_template.F90 must be included twice
! for the legacy and for the new API. In the new API, however, some routines
! must be named "..._impl"

!#include "elpa1_solve_tridi_real_template.F90"




!cannot use "../src/elpa1/../solve_tridi/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


subroutine solve_tridi_&
&single_impl &
    ( obj, na, nev, d, e, q, ldq, nblk, matrixCols, mpi_comm_all, mpi_comm_rows, &
                                           mpi_comm_cols, useGPU, wantDebug, success, max_threads )

      use precision
      use elpa_abstract_impl
      use merge_recursive
      use merge_systems
      use elpa_mpi
      use ELPA_utilities
      use distribute_global_column
      use elpa_mpi
      implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_FLOAT
  integer, parameter :: rck = C_FLOAT
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik), intent(in)               :: na, nev, ldq, nblk, matrixCols, &
                                                    mpi_comm_all, mpi_comm_rows, mpi_comm_cols
      real(kind=rk4), intent(inout)    :: d(na), e(na)
      real(kind=rk4), intent(inout)    :: q(ldq,*)
      logical, intent(in)                        :: useGPU, wantDebug
      logical, intent(out)                       :: success

      integer(kind=ik)                           :: i, j, n, np, nc, nev1, l_cols, l_rows
      integer(kind=ik)                           :: my_prow, my_pcol, np_rows, np_cols
      integer(kind=MPI_KIND)                     :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
      integer(kind=ik), allocatable              :: limits(:), l_col(:), p_col(:), l_col_bc(:), p_col_bc(:)

      integer(kind=ik)                           :: istat
      character(200)                             :: errorMessage
      character(20)                              :: gpuString
      integer(kind=ik), intent(in)               :: max_threads

      if(useGPU) then
        gpuString = "_gpu"
      else
        gpuString = ""
      endif

      call obj%timer%start("solve_tridi" // "_single" // gpuString)

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      success = .true.

      l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
      l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

      ! Set Q to 0
      q(1:l_rows, 1:l_cols) = 0.0_rk

      ! Get the limits of the subdivisons, each subdivison has as many cols
      ! as fit on the respective processor column

      allocate(limits(0:np_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: limits", 128,  istat,  errorMessage)

      limits(0) = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1) ! number of columns on proc column np

        ! Check for the case that a column has have zero width.
        ! This is not supported!
        ! Scalapack supports it but delivers no results for these columns,
        ! which is rather annoying
        if (nc==0) then
          call obj%timer%stop("solve_tridi" // "_single")
          if (wantDebug) write(error_unit,*) 'ELPA1_solve_tridi: ERROR: Problem contains processor column with zero width'
          success = .false.
          return
        endif
        limits(np+1) = limits(np) + nc
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,np_cols-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      ! Solve sub problems on processsor columns

      nc = limits(my_pcol) ! column after which my problem starts

      if (np_cols>1) then
        nev1 = l_cols ! all eigenvectors are needed
      else
        nev1 = MIN(nev,l_cols)
      endif
      call solve_tridi_col_&
           &single_impl &
             (obj, l_cols, nev1, nc, d(nc+1), e(nc+1), q, ldq, nblk,  &
                        matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads)
      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_single" // gpuString)
        return
      endif
      ! If there is only 1 processor column, we are done

      if (np_cols==1) then
        deallocate(limits, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi: limits", 176,  istat,  errorMessage)

        call obj%timer%stop("solve_tridi" // "_single" // gpuString)
        return
      endif

      ! Set index arrays for Q columns

      ! Dense distribution scheme:

      allocate(l_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col", 187,  istat,  errorMessage)

      allocate(p_col(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col", 190,  istat,  errorMessage)

      n = 0
      do np=0,np_cols-1
        nc = local_index(na, np, np_cols, nblk, -1)
        do i=1,nc
          n = n+1
          l_col(n) = i
          p_col(n) = np
        enddo
      enddo

      ! Block cyclic distribution scheme, only nev columns are set:

      allocate(l_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: l_col_bc", 205,  istat,  errorMessage)

      allocate(p_col_bc(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("solve_tridi: p_col_bc", 208,  istat,  errorMessage)

      p_col_bc(:) = -1
      l_col_bc(:) = -1

      do i = 0, na-1, nblk*np_cols
        do j = 0, np_cols-1
          do n = 1, nblk
            if (i+j*nblk+n <= MIN(nev,na)) then
              p_col_bc(i+j*nblk+n) = j
              l_col_bc(i+j*nblk+n) = i/np_cols + n
             endif
           enddo
         enddo
      enddo

      ! Recursively merge sub problems
      call merge_recursive_&
           &single &
           (obj, 0, np_cols, ldq, matrixCols, nblk, &
           l_col, p_col, l_col_bc, p_col_bc, limits, &
           np_cols, na, q, d, e, &
           mpi_comm_all, mpi_comm_rows, mpi_comm_cols,&
           useGPU, wantDebug, success, max_threads)

      if (.not.(success)) then
        call obj%timer%stop("solve_tridi" // "_single" // gpuString)
        return
      endif

      deallocate(limits,l_col,p_col,l_col_bc,p_col_bc, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi: limits, l_col, p_col, l_col_bc, p_col_bc", 239,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi" // "_single" // gpuString)
      return

    end subroutine solve_tridi_&
        &single_impl

    subroutine solve_tridi_col_&
    &single_impl &
      ( obj, na, nev, nqoff, d, e, q, ldq, nblk, matrixCols, mpi_comm_rows, useGPU, wantDebug, success, max_threads )

   ! Solves the symmetric, tridiagonal eigenvalue problem on one processor column
   ! with the divide and conquer method.
   ! Works best if the number of processor rows is a power of 2!
      use precision
      use elpa_abstract_impl
      use elpa_mpi
      use merge_systems
      use ELPA_utilities
      use distribute_global_column
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik)              :: na, nev, nqoff, ldq, nblk, matrixCols, mpi_comm_rows
      real(kind=rk4)      :: d(na), e(na)
      real(kind=rk4)      :: q(ldq,*)

      integer(kind=ik), parameter   :: min_submatrix_size = 16 ! Minimum size of the submatrices to be used

      real(kind=rk4), allocatable    :: qmat1(:,:), qmat2(:,:)
      integer(kind=ik)              :: i, n, np
      integer(kind=ik)              :: ndiv, noff, nmid, nlen, max_size
      integer(kind=ik)              :: my_prow, np_rows
      integer(kind=MPI_KIND)        :: mpierr, my_prowMPI, np_rowsMPI

      integer(kind=ik), allocatable :: limits(:), l_col(:), p_col_i(:), p_col_o(:)
      logical, intent(in)           :: useGPU, wantDebug
      logical, intent(out)          :: success
      integer(kind=ik)              :: istat
      character(200)                :: errorMessage

      integer(kind=ik), intent(in)  :: max_threads

      integer(kind=MPI_KIND)        :: bcast_request1, bcast_request2
      logical                       :: useNonBlockingCollectivesRows
      integer(kind=c_int)           :: non_blocking_collectives, error

      call obj%timer%start("solve_tridi_col" // "_single")

      call obj%get("nbc_row_solve_tridi", non_blocking_collectives, error)
      if (error .ne. ELPA_OK) then
        print *,"Problem setting option for non blocking collectives for rows in solve_tridi. Aborting..."
        stop
      endif

      if (non_blocking_collectives .eq. 1) then
        useNonBlockingCollectivesRows = .true.
      else
        useNonBlockingCollectivesRows = .false.
      endif

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")
      success = .true.
      ! Calculate the number of subdivisions needed.

      n = na
      ndiv = 1
      do while(2*ndiv<=np_rows .and. n>2*min_submatrix_size)
        n = ((n+3)/4)*2 ! the bigger one of the two halves, we want EVEN boundaries
        ndiv = ndiv*2
      enddo

      ! If there is only 1 processor row and not all eigenvectors are needed
      ! and the matrix size is big enough, then use 2 subdivisions
      ! so that merge_systems is called once and only the needed
      ! eigenvectors are calculated for the final problem.

      if (np_rows==1 .and. nev<na .and. na>2*min_submatrix_size) ndiv = 2

      allocate(limits(0:ndiv), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits", 442,  istat,  errorMessage)

      limits(0) = 0
      limits(ndiv) = na

      n = ndiv
      do while(n>1)
        n = n/2 ! n is always a power of 2
        do i=0,ndiv-1,2*n
          ! We want to have even boundaries (for cache line alignments)
          limits(i+n) = limits(i) + ((limits(i+2*n)-limits(i)+3)/4)*2
        enddo
      enddo

      ! Calculate the maximum size of a subproblem

      max_size = 0
      do i=1,ndiv
        max_size = MAX(max_size,limits(i)-limits(i-1))
      enddo

      ! Subdivide matrix by subtracting rank 1 modifications

      do i=1,ndiv-1
        n = limits(i)
        d(n) = d(n)-abs(e(n))
        d(n+1) = d(n+1)-abs(e(n))
      enddo

      if (np_rows==1)    then

        ! For 1 processor row there may be 1 or 2 subdivisions
        do n=0,ndiv-1
          noff = limits(n)        ! Start of subproblem
          nlen = limits(n+1)-noff ! Size of subproblem

          call solve_tridi_single_problem_&
          &single_impl &
                                  (obj, nlen,d(noff+1),e(noff+1), &
                                    q(nqoff+noff+1,noff+1),ubound(q,dim=1), wantDebug, success)

          if (.not.(success)) return
        enddo

      else

        ! Solve sub problems in parallel with solve_tridi_single
        ! There is at maximum 1 subproblem per processor

        allocate(qmat1(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1", 492,  istat,  errorMessage)

        allocate(qmat2(max_size,max_size), stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat2", 495,  istat,  errorMessage)

        qmat1 = 0 ! Make sure that all elements are defined

        if (my_prow < ndiv) then

          noff = limits(my_prow)        ! Start of subproblem
          nlen = limits(my_prow+1)-noff ! Size of subproblem
          call solve_tridi_single_problem_&
          &single_impl &
                                    (obj, nlen,d(noff+1),e(noff+1),qmat1, &
                                    ubound(qmat1,dim=1), wantDebug, success)

          if (.not.(success)) return
        endif

        ! Fill eigenvectors in qmat1 into global matrix q

        do np = 0, ndiv-1

          noff = limits(np)
          nlen = limits(np+1)-noff
          if (useNonBlockingCollectivesRows) then
            call obj%timer%start("mpi_nbc_communication")
            call mpi_ibcast(d(noff+1), int(nlen,kind=MPI_KIND), MPI_REAL4, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), bcast_request1, mpierr)
            call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)

            qmat2 = qmat1
            call mpi_ibcast(qmat2, int(max_size*max_size,kind=MPI_KIND), MPI_REAL4, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), bcast_request2, mpierr)
            call mpi_wait(bcast_request2, MPI_STATUS_IGNORE, mpierr)
            call obj%timer%stop("mpi_nbc_communication")
          else
            call obj%timer%start("mpi_communication")
            call mpi_bcast(d(noff+1), int(nlen,kind=MPI_KIND), MPI_REAL4, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)

            qmat2 = qmat1
            call mpi_bcast(qmat2, int(max_size*max_size,kind=MPI_KIND), MPI_REAL4, int(np,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")
          endif
          do i=1,nlen

            call distribute_global_column_&
            &single &
                     (obj, qmat2(1,i), q(1,noff+i), nqoff+noff, nlen, my_prow, np_rows, nblk)
          enddo

        enddo

        deallocate(qmat1, qmat2, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("solve_tridi_col: qmat1, qmat2", 558,  istat,  errorMessage)

      endif

      ! Allocate and set index arrays l_col and p_col

      allocate(l_col(na), p_col_i(na),  p_col_o(na), stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: l_col, p_col_i, p_col_o", 565,  istat,  errorMessage)

      do i=1,na
        l_col(i) = i
        p_col_i(i) = 0
        p_col_o(i) = 0
      enddo

      ! Merge subproblems

      n = 1
      do while(n<ndiv) ! if ndiv==1, the problem was solved by single call to solve_tridi_single

        do i=0,ndiv-1,2*n

          noff = limits(i)
          nmid = limits(i+n) - noff
          nlen = limits(i+2*n) - noff

          if (nlen == na) then
            ! Last merge, set p_col_o=-1 for unneeded (output) eigenvectors
            p_col_o(nev+1:na) = -1
          endif
          call merge_systems_&
          &single &
                              (obj, nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, nqoff+noff, nblk, &
                               matrixCols, int(mpi_comm_rows,kind=ik), int(mpi_comm_self,kind=ik), &
                               l_col(noff+1), p_col_i(noff+1), &
                               l_col(noff+1), p_col_o(noff+1), 0, 1, useGPU, wantDebug, success, max_threads)
          if (.not.(success)) return

        enddo

        n = 2*n

      enddo

      deallocate(limits, l_col, p_col_i, p_col_o, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("solve_tridi_col: limits, l_col, p_col_i, p_col_o", 603,  istat,  errorMessage)

      call obj%timer%stop("solve_tridi_col" // "_single")

    end subroutine solve_tridi_col_&
    &single_impl

    subroutine solve_tridi_single_problem_&
    &single_impl &
    (obj, nlen, d, e, q, ldq, wantDebug, success)

   ! Solves the symmetric, tridiagonal eigenvalue problem on a single processor.
   ! Takes precautions if DSTEDC fails or if the eigenvalues are not ordered correctly.
     use precision
     use elpa_abstract_impl
     use elpa_blas_interfaces
     use ELPA_utilities
     implicit none
     class(elpa_abstract_impl_t), intent(inout) :: obj
     integer(kind=ik)                         :: nlen, ldq
     real(kind=rk4)                 :: d(nlen), e(nlen), q(ldq,nlen)

     real(kind=rk4), allocatable    :: work(:), qtmp(:), ds(:), es(:)
     real(kind=rk4)                 :: dtmp

     integer(kind=ik)              :: i, j, lwork, liwork, info
     integer(kind=BLAS_KIND)       :: infoBLAS
     integer(kind=ik), allocatable :: iwork(:)

     logical, intent(in)           :: wantDebug
     logical, intent(out)          :: success
      integer(kind=ik)             :: istat
      character(200)               :: errorMessage

     call obj%timer%start("solve_tridi_single" // "_single")

     success = .true.
     allocate(ds(nlen), es(nlen), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: ds, es", 641,  istat,  errorMessage)

     ! Save d and e for the case that dstedc fails

     ds(:) = d(:)
     es(:) = e(:)

     ! First try dstedc, this is normally faster but it may fail sometimes (why???)

     lwork = 1 + 4*nlen + nlen**2
     liwork =  3 + 5*nlen
     allocate(work(lwork), iwork(liwork), stat=istat, errmsg=errorMessage)
     call check_allocate_f("solve_tridi_single: work, iwork", 653,  istat,  errorMessage)
     call obj%timer%start("blas")
     call SSTEDC('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND),    &
                          work, int(lwork,kind=BLAS_KIND), int(iwork,kind=BLAS_KIND), int(liwork,kind=BLAS_KIND), &
                          infoBLAS)
     info = int(infoBLAS,kind=ik)
     call obj%timer%stop("blas")

     if (info /= 0) then

       ! DSTEDC failed, try DSTEQR. The workspace is enough for DSTEQR.

       write(error_unit,'(a,i8,a)') 'Warning: Lapack routine DSTEDC failed, info= ',info,', Trying DSTEQR!'

       d(:) = ds(:)
       e(:) = es(:)
       call obj%timer%start("blas")
       call SSTEQR('I', int(nlen,kind=BLAS_KIND), d, e, q, int(ldq,kind=BLAS_KIND), work, infoBLAS )
       info = int(infoBLAS,kind=ik)
       call obj%timer%stop("blas")

       ! If DSTEQR fails also, we don't know what to do further ...

       if (info /= 0) then
         if (wantDebug) &
           write(error_unit,'(a,i8,a)') 'ELPA1_solve_tridi_single: ERROR: Lapack routine DSTEQR failed, info= ',info,', Aborting!'
           success = .false.
           return
         endif
       end if

       deallocate(work,iwork,ds,es, stat=istat, errmsg=errorMessage)
       call check_deallocate_f("solve_tridi_single: work, iwork, ds, es", 685,  istat,  errorMessage)

      ! Check if eigenvalues are monotonically increasing
      ! This seems to be not always the case  (in the IBM implementation of dstedc ???)

      do i=1,nlen-1
        if (d(i+1)<d(i)) then
          if (abs(d(i+1) - d(i)) / abs(d(i+1) + d(i)) > 1e-14_rk4) then
            write(error_unit,'(a,i8,2g25.16)') '***WARNING: Monotony error dste**:',i+1,d(i),d(i+1)
          else
            write(error_unit,'(a,i8,2g25.16)') 'Info: Monotony error dste{dc,qr}:',i+1,d(i),d(i+1)
            write(error_unit,'(a)') 'The eigenvalues from a lapack call are not sorted to machine precision.'
            write(error_unit,'(a)') 'In this extent, this is completely harmless.'
            write(error_unit,'(a)') 'Still, we keep this info message just in case.'
          end if
          allocate(qtmp(nlen), stat=istat, errmsg=errorMessage)
          call check_allocate_f("solve_tridi_single: qtmp", 705,  istat,  errorMessage)

          dtmp = d(i+1)
          qtmp(1:nlen) = q(1:nlen,i+1)
          do j=i,1,-1
            if (dtmp<d(j)) then
              d(j+1)        = d(j)
              q(1:nlen,j+1) = q(1:nlen,j)
            else
              exit ! Loop
            endif
          enddo
          d(j+1)        = dtmp
          q(1:nlen,j+1) = qtmp(1:nlen)
          deallocate(qtmp, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("solve_tridi_single: qtmp", 720,  istat,  errorMessage)

       endif
     enddo
     call obj%timer%stop("solve_tridi_single" // "_single")

    end subroutine solve_tridi_single_problem_&
    &single_impl

!#include "elpa1_merge_systems_real_template.F90"










subroutine hh_transform_real_&
&single &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine DLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
  implicit none
!    Copyright 2011, A. Marek
!
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
  integer, parameter :: rk = C_FLOAT
  integer, parameter :: rck = C_FLOAT
  real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  real(kind=rk), intent(inout)       :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  real(kind=rk), intent(out)         :: xf, tau

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &real&
     	      &" // &
                   &"_single" )


  if ( XNORM_SQ==0.0_rk ) then

    if ( ALPHA>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHA**2 + XNORM_SQ ), ALPHA )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHA = XNORM_SQ / ALPHA

      TAU = ALPHA / BETA
      ALPHA = -ALPHA
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &real&
  &" // &
  &"_single" )

end subroutine hh_transform_real_&
    &single




! complex double precision

































!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &complex&
  &_&
  &double &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads, isSkewsymmetric)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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

  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: ck = C_DOUBLE_COMPLEX
  integer, parameter :: rck = C_DOUBLE_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  complex(kind=rck), intent(out)          :: tau(na)
  complex(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik), parameter                   :: max_stored_uv = 32
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  !integer(kind=ik)                             :: l_cols, l_rows
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs

  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: a_offset

  integer(kind=ik), intent(in)                  :: max_threads

  real(kind=rk)                                 :: vnorm2
  complex(kind=rck)                       :: vav, x, aux(2*max_stored_uv), aux1(2), aux2(2), vrl, xf
  complex(kind=rck)                             :: aux3(1)

  integer(kind=c_intptr_t)                      :: num
  complex(kind=rck), allocatable          :: tmp(:)
  complex(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  complex(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  complex(kind=rck), pointer             :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  complex(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &complex
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_&
  &complex&
  &" // &
  "_double" // &
  gpuString )

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for min_tile_size. Aborting..."
    stop
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows
  !
  allocate(tmp(MAX(max_local_rows,max_local_cols)), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
  &complex ", "tmp", istat, errorMessage)

  ! allocate v_row 1 element longer to allow store and broadcast tau together with it
  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &complex ", "v_row", istat, errorMessage)
  
      allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
       &complex ", "v_col", istat, errorMessage)
  
      allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &complex ", "u_col", istat, errorMessage)
  
      allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &complex ", "u_row", istat, errorMessage)
    else

      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 346,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 351,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 356,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))

      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 361,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))

      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_roes", 367,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 372,  successGPU)

      num = na * 8
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 381,  successGPU)

      num = na * 8
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 390,  successGPU)
    endif
  else
    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &complex ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_row", istat, errorMessage)
      
  endif


  tmp = 0
  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU .and. .not.(useIntelGPU) ) then
     successGPU = gpu_malloc(v_row_dev, max_local_rows * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_row_dev", 429,  successGPU)

     successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

     call check_alloc_GPU_f("tridiag: u_row_dev", 433,  successGPU)

     successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_col_dev", 436,  successGPU)

     successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: u_col_dev", 439,  successGPU)

     successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 442,  successGPU)

     successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 445,  successGPU)
  endif !useGPU

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) &
  d_vec(na) = real(a_mat(l_rows,l_cols), kind=rk)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! allocate memmory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 476,  successGPU)

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("tridiag: a_mat", 480,  successGPU)

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 484,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1

    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

      ! copy l_cols + 1 column of A to v_row
      if (useGPU) then
        if (useIntelGPU) then
          v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
        else
          a_offset = l_cols * matrixRows * size_of_datatype
          ! we use v_row on the host at the moment! successGPU = gpu_memcpy(v_row_dev, a_dev + a_offset, 
          ! (l_rows)*size_of_PRECISION_real, gpuMemcpyDeviceToDevice)

          successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), &
                                    a_dev + a_offset, (l_rows)* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag a_dev 1", 519,  successGPU)
        endif
      else
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        aux(1:2*n_stored_vecs) = conjg(uv_stored_cols(l_cols+1,1:2*n_stored_vecs))
        if (useIntelGPU) then
                !print *,"intel phase aaaaaaaaaaaaaaaaaaaaaaaaaa"
          if (wantDebug) call obj%timer%start("mkl_offload")



          if (wantDebug) call obj%timer%stop("mkl_offload")
        else
          if (wantDebug) call obj%timer%start("blas")
          call ZGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            aux, 1_BLAS_KIND,  &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif
      endif

      if (my_prow == prow(istep-1, nblk, np_rows)) then
        aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1))
        aux1(2) = v_row(l_rows)
      else
        aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
        aux1(2) = 0.
      endif

      if (useNonBlockingCollectivesRows) then
        if (wantDebug) call obj%timer%start("mpi_communication_non_blocking")
        call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_DOUBLE_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
        call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication_non_blocking")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_DOUBLE_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif

      vnorm2 = real(aux2(1),kind=rk)
      vrl    = aux2(2)

      ! Householder transformation
      call hh_transform_complex_&
               &double &
               (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      ! Scale v_row and store Householder Vector for back transformation

      v_row(1:l_rows) = v_row(1:l_rows) * xf
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        v_row(l_rows) = 1.

        ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
        e_vec(istep-1) = real(vrl,kind=rk)
      endif

      ! store Householder Vector for back transformation
      a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)

      ! add tau after the end of actuall v_row, to be broadcasted with it
      v_row(l_rows+1) = tau(istep)
    endif !(my_pcol == pcol(istep, nblk, np_cols))

!          

    if (useNonBlockingCollectivesCols) then
      if (wantDebug) call obj%timer%start("mpi_nbc_communication")
      ! Broadcast the Householder Vector (and tau) along columns
      call mpi_ibcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   bcast_request1, mpierr)
      call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
      if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
    else
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_bcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")
    endif

    !recover tau, which has been broadcasted together with v_row
    tau(istep) =  v_row(l_rows+1)

    ! Transpose Householder Vector v_row -> v_col
    call elpa_transpose_vectors_&
        &complex&
        &_&
        &double &
              (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
               1, istep-1, 1, nblk, max_threads, .true.)

    ! Calculate u = (A + VU**T + UV**T)*v

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    u_col(1:l_cols) = 0
    u_row(1:l_rows) = 0
    if (l_rows > 0 .and. l_cols> 0 ) then
     if (useGPU .and. .not.(useIntelGPU)) then
       successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_col_dev", 693,  successGPU)

       successGPU = gpu_memset(u_row_dev, 0, l_rows * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_row_dev", 696,  successGPU)

       successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                     l_cols * size_of_datatype, gpuMemcpyHostToDevice)

       call check_memcpy_GPU_f("tridiag: v_col_dev", 701,  successGPU)

       successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                 l_rows * size_of_datatype, gpuMemcpyHostToDevice)
       call check_memcpy_GPU_f("tridiag: v_row_dev", 705,  successGPU)
     endif ! useGPU

     !if (useIntelGPU) then
     !  ! needed later when we can do explicit memcopy
     !endif

     do i= 0, (istep-2)/tile_size
       l_col_beg = i*l_cols_per_tile+1
       l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
       if (l_col_end < l_col_beg) cycle
       do j = 0, i
         l_row_beg = j*l_rows_per_tile+1
         l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
         if (l_row_end < l_row_beg) cycle

         ! multiplication by blocks is efficient only for CPU
         ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
         ! CPU implementation) or by one large matrix Vector multiply
         if (.not. useGPU) then
           if (wantDebug) call obj%timer%start("blas")
           call ZGEMV('C',  &
                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                       ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                       v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                       ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

           if (i/=j) then
             if (isSkewsymmetric) then
               call ZGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                   -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)

             else
               call ZGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                   ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)
             endif
           endif
           if (wantDebug) call obj%timer%stop("blas")
         endif ! not useGPU

            enddo  ! j=0,i
          enddo  ! i=0,(istep-2)/tile_size

          if (useGPU) then
            if (mat_vec_as_one_block) then
              if (useIntelGPU) then
                 if (wantDebug) call obj%timer%start("mkl_offload")

                if (wantDebug) call obj%timer%stop("mkl_offload")

              else
                ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
                ! this requires altering of the algorithm when later explicitly updating the matrix
                ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
                if (wantDebug) call obj%timer%start("gpublas")
                call gpublas_ZGEMV('C', l_rows,l_cols,  &
                                          ONE, a_dev, matrixRows,                   &
                                          v_row_dev , 1,                          &
                                          ONE, u_col_dev, 1)

       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_ZGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + a_offset, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
                if (wantDebug) call obj%timer%stop("gpublas")
              endif
            else  ! mat_vec_as_one_block
              !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
              !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")

                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype

                    call gpublas_ZGEMV('C', &
                                l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                ONE, a_dev + a_offset, matrixRows,  &
                                v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                                ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1)
                endif
              enddo

              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,i*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")


                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype
                    if (isSkewsymmetric) then
                       call gpublas_ZGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   -ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                    else
                       call gpublas_ZGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                   endif
                endif
              enddo
            end if !multiplication as one block / per stripes

            if (.not.(useIntelGPU)) then
              successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), &
                          u_col_dev, l_cols * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_col_dev 1", 954,  successGPU)

              successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), &
                          u_row_dev, l_rows * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_row_dev 1", 958,  successGPU)
            endif
            !if (useIntelGPU) then
            !  
            !endif

          endif ! useGPU


         ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
         if (n_stored_vecs > 0) then
           if (wantDebug) call obj%timer%start("blas")
           call ZGEMV('C',     &
                               int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                               v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)

           call ZGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                               aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
           if (wantDebug) call obj%timer%stop("blas")
         endif

       endif  ! (l_rows>0 .and. l_cols>0)

       ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
       ! on the processors containing the diagonal
       ! This is only necessary if u_row has been calculated, i.e. if the
       ! global tile size is smaller than the global remaining matrix

       if (tile_size < istep-1) then

         call elpa_reduce_add_vectors_&
         &complex&
         &_&
         &double &
         (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
         mpi_comm_cols, istep-1, 1, nblk, max_threads)

       endif

       ! Sum up all the u_col(:) parts, transpose u_col -> u_row

       if (l_cols>0) then
         tmp(1:l_cols) = u_col(1:l_cols)
         if (useNonBlockingCollectivesRows) then
           if (wantDebug) call obj%timer%start("mpi_nbc_communication")
           call mpi_iallreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
           call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
         else
           if (wantDebug) call obj%timer%start("mpi_communication")
           call mpi_allreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
           if (wantDebug) call obj%timer%stop("mpi_communication")
         endif
       endif
       if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
          &complex&
          &_&
          &double &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       else
          call elpa_transpose_vectors_&
          &complex&
          &_&
          &double &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       endif

       ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
       x = 0
       if (l_cols>0)  &
       x = dot_product(v_col(1:l_cols),u_col(1:l_cols))

       if (useNonBlockingCollectivesCols) then
         if (wantDebug) call obj%timer%start("mpi_nbc_communication")
         call mpi_iallreduce(x, vav, 1_MPI_KIND, MPI_DOUBLE_COMPLEX, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
               allreduce_request3, mpierr)
         call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
       else
         if (wantDebug) call obj%timer%start("mpi_communication")
         call mpi_allreduce(x, vav, 1_MPI_KIND, MPI_DOUBLE_COMPLEX, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
                mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")
       endif

       ! store u and v in the matrices U and V
       ! these matrices are stored combined in one here

       do j=1,l_rows
         vu_stored_rows(j,2*n_stored_vecs+1) = conjg(tau(istep))*v_row(j)
         vu_stored_rows(j,2*n_stored_vecs+2) = 0.5*conjg(tau(istep))*vav*v_row(j) - u_row(j)
       enddo
       do j=1,l_cols
         uv_stored_cols(j,2*n_stored_vecs+1) = 0.5*conjg(tau(istep))*vav*v_col(j) - u_col(j)
         uv_stored_cols(j,2*n_stored_vecs+2) = conjg(tau(istep))*v_col(j)
       enddo

       ! We have calculated another Hauseholder Vector, number of implicitly stored increased
       n_stored_vecs = n_stored_vecs+1

       ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
       if (n_stored_vecs == max_stored_uv .or. istep == 3) then

         if (useGPU .and. .not.(useIntelGPU)) then
           successGPU = gpu_memcpy(vu_stored_rows_dev, int(loc(vu_stored_rows(1,1)),kind=c_intptr_t), &
                                     max_local_rows * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_rows_dev", 1109,  successGPU)

           successGPU = gpu_memcpy(uv_stored_cols_dev, int(loc(uv_stored_cols(1,1)),kind=c_intptr_t), &
                                     max_local_cols * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_cols_dev", 1114,  successGPU)
         endif
          !if (useIntelGPU) then
          !  ! needed later when we can do explicit offloads
          !endif

         do i = 0, (istep-2)/tile_size
           ! go over tiles above (or on) the diagonal
           l_col_beg = i*l_cols_per_tile+1
           l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
           l_row_beg = 1
           l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
           if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) &
           cycle


           if (useGPU) then
             if (.not. mat_vec_as_one_block) then
               if (useIntelGPU) then
                  if (wantDebug) call obj%timer%start("mkl_offload")

                  call ZGEMM('N', 'C',     &
                                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), &
                                       int(l_col_end-l_col_beg+1,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),     &
                                       ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv),                    &
                                       int(max_local_rows,kind=BLAS_KIND), uv_stored_cols(l_col_beg,1),                    &
                                       int(max_local_cols,kind=BLAS_KIND), ONE, a_mat(l_row_beg,l_col_beg),                &
                                       int(matrixRows,kind=BLAS_KIND))
                  if (wantDebug) call obj%timer%stop("mkl_offload")

               else
                 ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
                 ! we than use the same calls as for CPU version
                 if (wantDebug) call obj%timer%start("gpublas")
                 call gpublas_ZGEMM('N', 'C',     &
                                         l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                         ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                         size_of_datatype,  &
                                         max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                         size_of_datatype,  &
                                         max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                         size_of_datatype , matrixRows)
                 if (wantDebug) call obj%timer%stop("gpublas")
               endif
             endif
           else !useGPU
             if (wantDebug) call obj%timer%start("blas")
             call ZGEMM('N', 'C',                &
                                  int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                  int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                  ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                  int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                  uv_stored_cols(l_col_beg,1), &
                                  int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                  ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
             if (wantDebug) call obj%timer%stop("blas")
           endif !useGPU
         enddo

         if (useGPU) then
           if (mat_vec_as_one_block) then
             if (useIntelGPU) then
                if (wantDebug) call obj%timer%start("mkl_offload")
                call ZGEMM('N', 'C', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                                    int(2*n_stored_vecs, kind=BLAS_KIND), ONE,  &
                                    vu_stored_rows, int(max_local_rows,kind=BLAS_KIND), &
                                    uv_stored_cols, int(max_local_cols,kind=BLAS_KIND),  &
                                    ONE, a_mat, int(matrixRows,kind=BLAS_KIND))
                if (wantDebug) call obj%timer%stop("mkl_offload")
             else
               !update whole (remaining) part of matrix, including tiles below diagonal
               !we can do that in one large cublas call
               if (wantDebug) call obj%timer%start("gpublas")
               call gpublas_ZGEMM('N', 'C', l_rows, l_cols, 2*n_stored_vecs,   &
                                         ONE, vu_stored_rows_dev, max_local_rows, &
                                         uv_stored_cols_dev, max_local_cols,  &
                                         ONE, a_dev, matrixRows)
               if (wantDebug) call obj%timer%stop("gpublas")
             endif
           endif
         endif

         n_stored_vecs = 0
       endif

       if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then
         if (useGPU) then
           if (useIntelGPU) then
                       ! if (useIntelGPU) then
          ! needed at a later time when we can do explcit mem copys
          ! endif

           else
             !a_mat(l_rows,l_cols) = a_dev(l_rows,l_cols)
              a_offset = ((l_rows - 1) + matrixRows * (l_cols - 1)) * size_of_datatype

              successGPU = gpu_memcpy(int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), a_dev + a_offset, &
                                      1 *  size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: a_dev 3", 1239,  successGPU)
           endif
         endif
         if (n_stored_vecs > 0) then
           a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) &
                       + dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs),uv_stored_cols(l_cols,1:2*n_stored_vecs))
         end if
         d_vec(istep-1) = real(a_mat(l_rows,l_cols),kind=rk)

         if (useGPU) then
           if (useIntelGPU) then
          ! if (useIntelGPU) then
          ! needed at a later time when we can expicit mem copy
          ! endif
           else
             !a_dev(l_rows,l_cols) = a_mat(l_rows,l_cols)
             !successGPU = gpu_threadsynchronize()
             !call check_memcpy_GPU_f("tridiag: a_dev 4a5a", 1265,  successGPU)

             successGPU = gpu_memcpy(a_dev + a_offset, int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), &
                                     int(1 * size_of_datatype, kind=c_intptr_t), gpuMemcpyHostToDevice)
             call check_memcpy_GPU_f("tridiag: a_dev 4", 1269,  successGPU)
           endif
         endif
       endif

     enddo ! main cycle over istep=na,3,-1

     ! Store e_vec(1) and d_vec(1)

     if (my_pcol==pcol(2, nblk, np_cols)) then
      if (my_prow==prow(1, nblk, np_rows)) then
       ! We use last l_cols value of loop above
       if (useGPU) then
         if (useIntelGPU) then
            vrl = a_mat(1,l_cols)
         else
           successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                                   1 * size_of_datatype, gpuMemcpyDeviceToHost)
           call check_memcpy_GPU_f("tridiag: a_dev 5", 1288,  successGPU)
           vrl = aux3(1)
         endif
       else !useGPU
         vrl = a_mat(1,l_cols)
       endif !useGPU
       call hh_transform_complex_&
       &double &
       (obj, vrl, 0.0_rk, xf, tau(2), wantDebug)
       e_vec(1) = real(vrl,kind=rk)
       a_mat(1,l_cols) = 1. ! for consistency only
     endif
     if (useNonBlockingCollectivesRows) then
       if (wantDebug) call obj%timer%start("mpi_nbc_communication")
       call mpi_ibcast(tau(2), 1_MPI_KIND, MPI_DOUBLE_COMPLEX, int(prow(1, nblk, np_rows),kind=MPI_KIND), &
                   int(mpi_comm_rows,kind=MPI_KIND), bcast_request2, mpierr)
       call mpi_wait(bcast_request2, MPI_STATUS_IGNORE, mpierr)
       if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
     else
       if (wantDebug) call obj%timer%start("mpi_communication")
       call mpi_bcast(tau(2), 1_MPI_KIND, MPI_DOUBLE_COMPLEX, int(prow(1, nblk, np_rows),kind=MPI_KIND), &
                   int(mpi_comm_rows,kind=MPI_KIND),  mpierr)
       if (wantDebug) call obj%timer%stop("mpi_communication")
     endif

   endif

   if (useNonBlockingCollectivesCols) then
     if (wantDebug) call obj%timer%start("mpi_nbc_communication")
     call mpi_ibcast(tau(2), 1_MPI_KIND, MPI_DOUBLE_COMPLEX, int(pcol(2, nblk, np_cols),kind=MPI_KIND), &
                  int(mpi_comm_cols,kind=MPI_KIND), bcast_request3, mpierr)
     call mpi_wait(bcast_request3, MPI_STATUS_IGNORE, mpierr)
     if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
   else
     if (wantDebug) call obj%timer%start("mpi_communication")
     call mpi_bcast(tau(2), 1_MPI_KIND, MPI_DOUBLE_COMPLEX, int(pcol(2, nblk, np_cols),kind=MPI_KIND), &
                  int(mpi_comm_cols,kind=MPI_KIND), mpierr)
     if (wantDebug) call obj%timer%stop("mpi_communication")
   endif

  if (my_prow == prow(1, nblk, np_rows) .and. my_pcol == pcol(1, nblk, np_cols))  then
    if (useGPU) then
      if (useIntelGPU) then
        d_vec(1) = DREAL(a_mat(1,1))
      else
        successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev, &
                               1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 6", 1344,  successGPU)
        d_vec(1) = DREAL(aux3(1))
      endif
    else !useGPU
      d_vec(1) = DREAL(a_mat(1,1))
    endif !useGPU
  endif



  deallocate(tmp, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp", 1395,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 1400,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 1403,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 1406,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 1409,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 1412,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 1415,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 1418,  successGPU)
  endif
  ! if (useIntelGPU) then
  ! needed at a later time when we can do explicit frees
  ! endif


  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 1428,  istat,  errorMessage)

  if (useNonBlockingCollectivesRows) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request4, mpierr)
    call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6,mpierr)
    call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif
  if (useNonBlockingCollectivesCols) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request5, mpierr)
    call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request7, mpierr)
    call mpi_wait(allreduce_request7, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif

  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 1477,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
           deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
     call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1482,  istat,  errorMessage)
    else
      successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: a_mat", 1485,  successGPU)

      successGPU = gpu_free_host(v_row_host)
      call check_host_dealloc_GPU_f("tridiag: v_row_host", 1488,  successGPU)
      nullify(v_row)

      successGPU = gpu_free_host(v_col_host)
      call check_host_dealloc_GPU_f("tridiag: v_col_host", 1492,  successGPU)
      nullify(v_col)

      successGPU = gpu_free_host(u_col_host)
      call check_host_dealloc_GPU_f("tridiag: u_col_host", 1496,  successGPU)
      nullify(u_col)

      successGPU = gpu_free_host(u_row_host)
      call check_host_dealloc_GPU_f("tridiag: u_row_host", 1500,  successGPU)
      nullify(u_row)

      successGPU = gpu_host_unregister(int(loc(uv_stored_cols),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: uv_stored_cols", 1504,  successGPU)

      successGPU = gpu_host_unregister(int(loc(vu_stored_rows),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: vu_stored_rows", 1507,  successGPU)

      successGPU = gpu_host_unregister(int(loc(e_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: e_vec", 1510,  successGPU)

      successGPU = gpu_host_unregister(int(loc(d_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: d_vec", 1513,  successGPU)
    endif
  else
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1517,  istat,  errorMessage)
  endif

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 1521,  istat,  errorMessage)

  call obj%timer%stop("tridiag_&
  &complex&
  &" // &
  "_double" // &
  gpuString )

end subroutine tridiag_&
&complex&
&_&
&double





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>

subroutine trans_ev_&
&complex&
&_&
&double &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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

  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: ck = C_DOUBLE_COMPLEX
  integer, parameter :: rck = C_DOUBLE_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  complex(kind=rck), intent(in)           :: tau(na)

  complex(kind=rck), intent(inout)        :: a_mat(lda,*)
  complex(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  complex(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  complex(kind=rck), pointer              :: tmp1(:), tmp2(:)
  complex(kind=rck), allocatable          :: h1(:), h2(:)
  complex(kind=rck), pointer              :: tmat(:,:)
  complex(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_&
                                                                      &complex
  integer(kind=ik)                              :: error
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk

  if (useIntelGPU) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp2", istat, errorMessage)
  endif

  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif

  ! In the complex case tau(2) /= 0
  if (my_prow == prow(1, nblk, np_rows)) then
    q_mat(1,1:l_cols) = q_mat(1,1:l_cols)*(ONE-tau(2))
  endif

  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&complex&
    !&", "hvm1", istat, errorMessage)
    num = (max_local_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(hvm1_host,num)
    call check_alloc_GPU_f("trans_ev: hvm1_host", 304,  successGPU)
    call c_f_pointer(hvm1_host,hvm1,(/(max_local_rows*max_stored_rows)/))

    num = (max_stored_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmat_host,num)
    call check_alloc_GPU_f("trans_ev: tmat_host", 309,  successGPU)
    call c_f_pointer(tmat_host,tmat,(/max_stored_rows,max_stored_rows/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_alloc_GPU_f("trans_ev: tmp1_host", 314,  successGPU)
    call c_f_pointer(tmp1_host,tmp1,(/(max_local_cols*max_stored_rows)/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp2_host,num)
    call check_alloc_GPU_f("trans_ev: tmp2_host", 319,  successGPU)
    call c_f_pointer(tmp2_host,tmp2,(/(max_local_cols*max_stored_rows)/))

    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 323,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 326,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 329,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 333,  successGPU)

    successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev: q_mat", 337,  successGPU)

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 341,  successGPU)
  endif  ! useGPU

  !if (useIntelGPU) then
  !  ! needed at a later time when we can do explicit copys
  !endif


  do istep = 1, na, blockStep
    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo

    if (nb>0) then
      if (useNonBlockingCollectivesCols) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(hvb, int(nb,kind=MPI_KIND), MPI_DOUBLE_COMPLEX , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
       else
        call obj%timer%start("mpi_communication")
        call mpi_bcast(hvb, int(nb,kind=MPI_KIND), MPI_DOUBLE_COMPLEX , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif
    endif

    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) &
      call ZHERK('U', 'C',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo
      if (nc>0) then
        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif
      endif
      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call ZTRMV('L', 'C' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -conjg(h2(nc+1:nc+n)) &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU .and. .not.(useIntelGPU)) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 475,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 480,  successGPU)
      endif

      !if (useIntelGPU) then
      !  ! needed later when we can do explicit copys
      !endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            call obj%timer%start("mkl_offload")
            call obj%timer%stop("mkl_offload")

          else
            call obj%timer%start("gpublas")
            call gpublas_ZGEMM('C', 'N',   &
                                     nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                     q_dev, ldq, ZERO, tmp_dev, nstor)
            call obj%timer%stop("gpublas")
          endif
        else ! useGPU

          call obj%timer%start("blas")
          call ZGEMM('C', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (useIntelGPU) then
            tmp1(1:l_cols*nstor) = 0
          else
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 532,  successGPU)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0


      if (useGPU .and. .not.(useIntelGPU)) then
        ! In the legacy GPU version, this allreduce was ommited. But probably it has to be done for GPU + MPI
        ! todo: does it need to be copied whole? Wouldn't be a part sufficient?
        successGPU = gpu_memcpy(int(loc(tmp1(1)),kind=c_intptr_t), tmp_dev,  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("trans_ev", 547,  successGPU)
      endif

      !if (useIntelGPU) then
      !   ! needed later
      !endif

      if (useNonBlockingCollectivesRows) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_iallreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
        call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
      else
        call obj%timer%start("mpi_communication")
        call mpi_allreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif

      if (useGPU .and. .not.(useIntelGPU)) then
        ! copy back tmp2 - after reduction...
        successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2(1)),kind=c_intptr_t),  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
      endif ! useGPU


      !if (useIntelGPU) then
      !   ! needed later
      !endif


      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            ! tmp2 = tmat * tmp2
            call obj%timer%start("mkl_offload")

            call obj%timer%stop("mkl_offload")
          else
            call obj%timer%start("gpublas")
            call gpublas_ZTRMM('L', 'L', 'N', 'N',     &
                                     nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                     tmp_dev, nstor)

            call gpublas_ZGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                     -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                     ONE, q_dev, ldq)
            call obj%timer%stop("gpublas")
          endif
        else !useGPU
          ! tmp2 = tmat * tmp2
          call obj%timer%start("blas")
          call ZTRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                             ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND))
          !q_mat = q_mat - hvm*tmp2
          call ZGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(nstor,kind=BLAS_KIND),   &
                              -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND), &
                              ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

  enddo ! istep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 693,  istat,  errorMessage)

  if (useGPU) then

    if (useIntelGPU) then
      deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_&       &MATH_DATATYPE&       &: tmat, tmp1, tmp2", 701,  istat,  errorMessage)
    else

      !q_mat = q_dev
      successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                    q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev", 707,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev: q_mat", 710,  successGPU)

      successGPU = gpu_free_host(hvm1_host)
      call check_host_dealloc_GPU_f("trans_ev: hvm1_host", 713,  successGPU)
      nullify(hvm1)

      successGPU = gpu_free_host(tmat_host)
      call check_host_dealloc_GPU_f("trans_ev: tmat_host", 717,  successGPU)
      nullify(tmat)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp1_host", 721,  successGPU)
      nullify(tmp1)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp2_host", 725,  successGPU)
      nullify(tmp2)

      !deallocate(hvm1, stat=istat, errmsg=errorMessage)
      !if (istat .ne. 0) then
      !  print *,"trans_ev_&
      !  &complex&
      !  &: error when deallocating hvm1 "//errorMessage
      !  stop 1
      !endif

      !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev", 738,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev", 741,  successGPU)

      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev", 744,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev", 747,  successGPU)
    endif
  else
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 753,  istat,  errorMessage)
  endif


  call obj%timer%stop("trans_ev_&
  &complex&
  &" // &
  &"_double" // &
  gpuString )

end subroutine trans_ev_&
&complex&
&_&
&double







subroutine hh_transform_complex_&
&double &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine ZLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
  implicit none
!    Copyright 2011, A. Marek
!
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

  integer, parameter :: rk = C_DOUBLE
  integer, parameter :: ck = C_DOUBLE_COMPLEX
  integer, parameter :: rck = C_DOUBLE_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  complex(kind=ck), intent(inout) :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  complex(kind=ck), intent(out)   :: xf, tau
  real(kind=rk)                      :: ALPHR, ALPHI

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &complex&
     	      &" // &
                   &"_double" )

  ALPHR = real( ALPHA, kind=rk )
  ALPHI = DIMAG( ALPHA )

  if ( XNORM_SQ==0.0_rk .AND. ALPHI==0.0_rk ) then

    if ( ALPHR>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHR**2 + ALPHI**2 + XNORM_SQ ), ALPHR )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHR = ALPHI * (ALPHI/real( ALPHA , kind=rk))
      ALPHR = ALPHR + XNORM_SQ/real( ALPHA, kind=rk )

      TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA )
      ALPHA = DCMPLX( -ALPHR, ALPHI )
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &complex&
  &" // &
  &"_double" )

    end subroutine hh_transform_complex_&
    &double




! complex single precision
































!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param obj	      object of elpa_type
!> \param na          Order of matrix
!>
!> \param a_mat(matrixRows,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param matrixRows         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param d_vec(na)       Diagonal elements (returned), identical on all processors
!>
!> \param e_vec(na)       Off-Diagonal elements (returned), identical on all processors
!>
!> \param tau(na)     Factors for the Householder vectors (returned), needed for back transformation
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!> \param wantDebug   if true more debug information
!>
subroutine tridiag_&
  &complex&
  &_&
  &single &
  (obj, na, a_mat, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau, useGPU, wantDebug, &
   max_threads, isSkewsymmetric)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use matrix_plot
  use elpa_omp
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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

  integer, parameter :: rk = C_FLOAT
  integer, parameter :: ck = C_FLOAT_COMPLEX
  integer, parameter :: rck = C_FLOAT_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  logical, intent(in)                           :: useGPU, wantDebug
  logical, intent(in)                           :: isSkewsymmetric

  complex(kind=rck), intent(out)          :: tau(na)
  complex(kind=rck), intent(inout)        :: a_mat(matrixRows,*)
  real(kind=rk), intent(out)                    :: d_vec(na)
  real(kind=rk), intent(out)                    :: e_vec(na)
  integer(kind=ik), parameter                   :: max_stored_uv = 32
  logical,          parameter                   :: mat_vec_as_one_block = .true.

  ! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=MPI_KIND)                        :: mpierr
  integer(kind=ik)                              :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, &
                                                   max_local_cols
  ! updated after each istep (in the main cycle) to contain number of
  ! local columns and rows of the remaining part of the matrix
  !integer(kind=ik)                             :: l_cols, l_rows
  integer(kind=ik)                              :: l_cols, l_rows
  integer(kind=ik)                              :: n_stored_vecs

  integer(kind=C_intptr_T)                      :: a_dev, v_row_dev, v_col_dev, u_row_dev, u_col_dev, vu_stored_rows_dev, &
                                                   uv_stored_cols_dev
  logical                                       :: successGPU

  integer(kind=ik)                              :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik)                              :: tile_size, l_rows_per_tile, l_cols_per_tile
  integer(kind=c_intptr_t)                      :: a_offset

  integer(kind=ik), intent(in)                  :: max_threads

  real(kind=rk)                                 :: vnorm2
  complex(kind=rck)                       :: vav, x, aux(2*max_stored_uv), aux1(2), aux2(2), vrl, xf
  complex(kind=rck)                             :: aux3(1)

  integer(kind=c_intptr_t)                      :: num
  complex(kind=rck), allocatable          :: tmp(:)
  complex(kind=rck), pointer              :: v_row(:), & ! used to store calculated Householder Vector
                                                   v_col(:)   ! the same Vector, but transposed 
  complex(kind=rck), pointer              :: u_col(:), u_row(:)

  ! the following two matrices store pairs of vectors v and u calculated in each step
  ! at most max_stored_uv Vector pairs are stored, than the matrix A_i is explicitli updated
  ! u and v are stored both in row and Vector forms
  ! pattern: v1,u1,v2,u2,v3,u3,....
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  complex(kind=rck), pointer             :: vu_stored_rows(:,:)
  ! pattern: u1,v1,u2,v2,u3,v3,....
  complex(kind=rck), allocatable         :: uv_stored_cols(:,:)


  type(c_ptr)                                   :: v_row_host, v_col_host
  type(c_ptr)                                   :: u_row_host, u_col_host
  type(c_ptr)                                   :: vu_stored_rows_host, uv_stored_cols_host
  real(kind=rk), allocatable                    :: tmp_real(:)
  integer(kind=ik)                              :: min_tile_size, error
  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString
  integer(kind=ik)                              :: nblockEnd
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &complex
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, bcast_request2, bcast_request3
  integer(kind=MPI_KIND)                        :: allreduce_request1, allreduce_request2, allreduce_request3
  integer(kind=MPI_KIND)                        :: allreduce_request4, allreduce_request5, allreduce_request6, &
                                                   allreduce_request7
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_&
  &complex&
  &" // &
  "_single" // &
  gpuString )

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa1_full_to_tridi", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridiag. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_full_to_tridi", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridiag. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND), my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND), np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND), my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND), np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  ! seems that tile is a square submatrix, consisting by several blocks
  ! it is a smallest possible square submatrix, where blocks being distributed among
  ! processors are "aligned" in both rows and columns
  !  -----------------
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- ...
  ! | 1 4 | 1 4 | 1 4 | ...
  ! | 2 5 | 2 5 | 2 5 | ...
  ! | 3 6 | 3 6 | 3 6 | ...
  !  ----------------- .
  !   : :   : :   : :    .
  !   : :   : :   : :      .
  !
  ! this is a tile, where each number represents block, assigned to a processor with the shown number
  ! size of this small block is nblk
  ! Image is for situation with 6 processors, 3 processor rows and 2 columns
  ! tile_size is thus nblk * 6
  !
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size

  ! make tile_size a smallest possible multiple of previously defined tile size, such that it is
  ! larger or equal to min_tile_size
  ! min_tile_size has been originally hardcoded as 128 * max(np_rows, np_cols), so it is now the implicit value
  ! it can, however, be set by the user
  call obj%get("min_tile_size", min_tile_size ,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for min_tile_size. Aborting..."
    stop
  endif
  if(min_tile_size == 0) then
    ! not set by the user, use the default value
    min_tile_size = 128*max(np_rows, np_cols)
  endif
  tile_size = ((min_tile_size-1)/tile_size+1)*tile_size

  nblockEnd = 3

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  ! localy owned submatrix has size at most max_local_rows x max_local_cols at each processor
  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  ! allocate memmory for vectors
  ! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  ! todo: if something has length max_local_rows, it is actually a column, no?
  ! todo: probably one should read it as v_row = Vector v distributed among rows
  !
  allocate(tmp(MAX(max_local_rows,max_local_cols)), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
  &complex ", "tmp", istat, errorMessage)

  ! allocate v_row 1 element longer to allow store and broadcast tau together with it
  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "uv_stored_cols", istat, errorMessage)

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_&
       &complex ", "vu_stored_rows", istat, errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &complex ", "v_row", istat, errorMessage)
  
      allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
       &complex ", "v_col", istat, errorMessage)
  
      allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &complex ", "u_col", istat, errorMessage)
  
      allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
      call check_alloc("tridiag_&
      &complex ", "u_row", istat, errorMessage)
    else

      num = (max_local_rows+1) * size_of_datatype
      successGPU = gpu_malloc_host(v_row_host, num)
      call check_host_alloc_GPU_f("tridiag: v_row_host", 346,  successGPU)
      call c_f_pointer(v_row_host,v_row,(/(max_local_rows+1)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(v_col_host,num)
      call check_host_alloc_GPU_f("tridiag: v_col_host", 351,  successGPU)
      call c_f_pointer(v_col_host,v_col,(/(max_local_cols)/))

      num = (max_local_cols) * size_of_datatype
      successGPU = gpu_malloc_host(u_col_host,num)
      call check_host_alloc_GPU_f("tridiag: u_col_host", 356,  successGPU)
      call c_f_pointer(u_col_host,u_col,(/(max_local_cols)/))

      num = (max_local_rows) * size_of_datatype
      successGPU = gpu_malloc_host(u_row_host,num)
      call check_host_alloc_GPU_f("tridiag: u_row_host", 361,  successGPU)
      call c_f_pointer(u_row_host,u_row,(/(max_local_rows)/))

      num = (max_local_rows * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(vu_stored_rows),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: vu_stored_roes", 367,  successGPU)

      num = (max_local_cols * 2*max_stored_uv) * size_of_datatype
      successGPU = gpu_host_register(int(loc(uv_stored_cols),kind=c_intptr_t),num,&
                    gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: uv_stored_cols", 372,  successGPU)

      num = na * 4
      successGPU = gpu_host_register(int(loc(e_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: e_vec", 381,  successGPU)

      num = na * 4
      successGPU = gpu_host_register(int(loc(d_vec),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
      call check_host_register_GPU_f("tridiag: d_vec", 390,  successGPU)
    endif
  else
    allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "v_row", istat, errorMessage)

    allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
     &complex ", "v_col", istat, errorMessage)

    allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_col", istat, errorMessage)

    allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("tridiag_&
    &complex ", "u_row", istat, errorMessage)
      
  endif


  tmp = 0
  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  if (useGPU .and. .not.(useIntelGPU) ) then
     successGPU = gpu_malloc(v_row_dev, max_local_rows * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_row_dev", 429,  successGPU)

     successGPU = gpu_malloc(u_row_dev, max_local_rows * size_of_datatype)

     call check_alloc_GPU_f("tridiag: u_row_dev", 433,  successGPU)

     successGPU = gpu_malloc(v_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: v_col_dev", 436,  successGPU)

     successGPU = gpu_malloc(u_col_dev, max_local_cols * size_of_datatype)
     call check_alloc_GPU_f("tridiag: u_col_dev", 439,  successGPU)

     successGPU = gpu_malloc(vu_stored_rows_dev, max_local_rows * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 442,  successGPU)

     successGPU = gpu_malloc(uv_stored_cols_dev, max_local_cols * 2 * max_stored_uv * size_of_datatype)
     call check_alloc_GPU_f("tridiag: vu_stored_rows_dev", 445,  successGPU)
  endif !useGPU

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat

  if (my_prow == prow(na, nblk, np_rows) .and. my_pcol == pcol(na, nblk, np_cols)) &
  d_vec(na) = real(a_mat(l_rows,l_cols), kind=rk)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! allocate memmory for matrix A on the device and than copy the matrix

    num = matrixRows * matrixCols * size_of_datatype

    successGPU = gpu_malloc(a_dev, num)
    call check_alloc_GPU_f("tridiag: a_dev", 476,  successGPU)

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("tridiag: a_mat", 480,  successGPU)

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat(1,1)),kind=c_intptr_t), &
                              num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("tridiag: a_dev", 484,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif

  ! main cycle of tridiagonalization
  ! in each step, 1 Householder Vector is calculated
  do istep = na, nblockEnd ,-1

    ! Calculate number of local rows and columns of the still remaining matrix
    ! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

    ! Calculate Vector for Householder transformation on all procs
    ! owning column istep

    if (my_pcol == pcol(istep, nblk, np_cols)) then

      ! Get Vector to be transformed; distribute last element and norm of
      ! remaining elements to all procs in current column

      ! copy l_cols + 1 column of A to v_row
      if (useGPU) then
        if (useIntelGPU) then
          v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
        else
          a_offset = l_cols * matrixRows * size_of_datatype
          ! we use v_row on the host at the moment! successGPU = gpu_memcpy(v_row_dev, a_dev + a_offset, 
          ! (l_rows)*size_of_PRECISION_real, gpuMemcpyDeviceToDevice)

          successGPU = gpu_memcpy(int(loc(v_row),kind=c_intptr_t), &
                                    a_dev + a_offset, (l_rows)* size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("tridiag a_dev 1", 519,  successGPU)
        endif
      else
        v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)
      endif

      if (n_stored_vecs > 0 .and. l_rows > 0) then
        aux(1:2*n_stored_vecs) = conjg(uv_stored_cols(l_cols+1,1:2*n_stored_vecs))
        if (useIntelGPU) then
                !print *,"intel phase aaaaaaaaaaaaaaaaaaaaaaaaaa"
          if (wantDebug) call obj%timer%start("mkl_offload")



          if (wantDebug) call obj%timer%stop("mkl_offload")
        else
          if (wantDebug) call obj%timer%start("blas")
          call CGEMV('N',   &
                            int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND), &
                            ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND), &
                            aux, 1_BLAS_KIND,  &
                            ONE, v_row, 1_BLAS_KIND)
          if (wantDebug) call obj%timer%stop("blas")
        endif
      endif

      if (my_prow == prow(istep-1, nblk, np_rows)) then
        aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1))
        aux1(2) = v_row(l_rows)
      else
        aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
        aux1(2) = 0.
      endif

      if (useNonBlockingCollectivesRows) then
        if (wantDebug) call obj%timer%start("mpi_communication_non_blocking")
        call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
        call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication_non_blocking")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif

      vnorm2 = real(aux2(1),kind=rk)
      vrl    = aux2(2)

      ! Householder transformation
      call hh_transform_complex_&
               &single &
               (obj, vrl, vnorm2, xf, tau(istep), wantDebug)
      ! Scale v_row and store Householder Vector for back transformation

      v_row(1:l_rows) = v_row(1:l_rows) * xf
      if (my_prow == prow(istep-1, nblk, np_rows)) then
        v_row(l_rows) = 1.

        ! vrl is newly computed off-diagonal element of the final tridiagonal matrix
        e_vec(istep-1) = real(vrl,kind=rk)
      endif

      ! store Householder Vector for back transformation
      a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)

      ! add tau after the end of actuall v_row, to be broadcasted with it
      v_row(l_rows+1) = tau(istep)
    endif !(my_pcol == pcol(istep, nblk, np_cols))

!          

    if (useNonBlockingCollectivesCols) then
      if (wantDebug) call obj%timer%start("mpi_nbc_communication")
      ! Broadcast the Householder Vector (and tau) along columns
      call mpi_ibcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_COMPLEX,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   bcast_request1, mpierr)
      call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
      if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
    else
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_bcast(v_row, int(l_rows+1,kind=MPI_KIND), MPI_COMPLEX,    &
                   int(pcol(istep, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                   mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")
    endif

    !recover tau, which has been broadcasted together with v_row
    tau(istep) =  v_row(l_rows+1)

    ! Transpose Householder Vector v_row -> v_col
    call elpa_transpose_vectors_&
        &complex&
        &_&
        &single &
              (obj, v_row, ubound(v_row,dim=1), mpi_comm_rows, v_col, ubound(v_col,dim=1), mpi_comm_cols, &
               1, istep-1, 1, nblk, max_threads, .true.)

    ! Calculate u = (A + VU**T + UV**T)*v

    ! For cache efficiency, we use only the upper half of the matrix tiles for this,
    ! thus the result is partly in u_col(:) and partly in u_row(:)

    u_col(1:l_cols) = 0
    u_row(1:l_rows) = 0
    if (l_rows > 0 .and. l_cols> 0 ) then
     if (useGPU .and. .not.(useIntelGPU)) then
       successGPU = gpu_memset(u_col_dev, 0, l_cols * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_col_dev", 693,  successGPU)

       successGPU = gpu_memset(u_row_dev, 0, l_rows * size_of_datatype)
       call check_memcpy_GPU_f("tridiag: u_row_dev", 696,  successGPU)

       successGPU = gpu_memcpy(v_col_dev, int(loc(v_col(1)),kind=c_intptr_t), &
                     l_cols * size_of_datatype, gpuMemcpyHostToDevice)

       call check_memcpy_GPU_f("tridiag: v_col_dev", 701,  successGPU)

       successGPU = gpu_memcpy(v_row_dev, int(loc(v_row(1)),kind=c_intptr_t), &
                                 l_rows * size_of_datatype, gpuMemcpyHostToDevice)
       call check_memcpy_GPU_f("tridiag: v_row_dev", 705,  successGPU)
     endif ! useGPU

     !if (useIntelGPU) then
     !  ! needed later when we can do explicit memcopy
     !endif

     do i= 0, (istep-2)/tile_size
       l_col_beg = i*l_cols_per_tile+1
       l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
       if (l_col_end < l_col_beg) cycle
       do j = 0, i
         l_row_beg = j*l_rows_per_tile+1
         l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
         if (l_row_end < l_row_beg) cycle

         ! multiplication by blocks is efficient only for CPU
         ! for GPU we introduced 2 other ways, either by stripes (more simmilar to the original
         ! CPU implementation) or by one large matrix Vector multiply
         if (.not. useGPU) then
           if (wantDebug) call obj%timer%start("blas")
           call CGEMV('C',  &
                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                       ONE, a_mat(l_row_beg, l_col_beg), int(matrixRows,kind=BLAS_KIND),         &
                       v_row(l_row_beg:max_local_rows+1), 1_BLAS_KIND,                           &
                       ONE, u_col(l_col_beg:max_local_cols), 1_BLAS_KIND)

           if (i/=j) then
             if (isSkewsymmetric) then
               call CGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                   -ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)

             else
               call CGEMV('N',int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND),  &
                                   ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND),               &
                                   v_col(l_col_beg:max_local_cols), 1_BLAS_KIND, ONE, u_row(l_row_beg:max_local_rows), &
                                   1_BLAS_KIND)
             endif
           endif
           if (wantDebug) call obj%timer%stop("blas")
         endif ! not useGPU

            enddo  ! j=0,i
          enddo  ! i=0,(istep-2)/tile_size

          if (useGPU) then
            if (mat_vec_as_one_block) then
              if (useIntelGPU) then
                 if (wantDebug) call obj%timer%start("mkl_offload")

                if (wantDebug) call obj%timer%stop("mkl_offload")

              else
                ! Unlike for CPU, we (for each MPI thread) do just one large mat-vec multiplication
                ! this requires altering of the algorithm when later explicitly updating the matrix
                ! after max_stored_uv is reached : we need to update all tiles, not only those above diagonal
                if (wantDebug) call obj%timer%start("gpublas")
                call gpublas_CGEMV('C', l_rows,l_cols,  &
                                          ONE, a_dev, matrixRows,                   &
                                          v_row_dev , 1,                          &
                                          ONE, u_col_dev, 1)

       ! todo: try with non transposed!!!
!                 if(i/=j) then
!                   call gpublas_CGEMV('N', l_row_end-l_row_beg+1,l_col_end-l_col_beg+1,  &
!                                             ONE, a_dev + a_offset, matrixRows,                        &
!                                             v_col_dev + (l_col_beg - 1) *                      &
!                                             size_of_datatype, 1,                          &
!                                             ONE, u_row_dev + (l_row_beg - 1) *                 &
!                                             size_of_datatype, 1)
!                 endif
                if (wantDebug) call obj%timer%stop("gpublas")
              endif
            else  ! mat_vec_as_one_block
              !perform multiplication by stripes - it is faster than by blocks, since we call cublas with
              !larger matrices. In general, however, this algorithm is very simmilar to the one with CPU
              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")

                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype

                    call gpublas_CGEMV('C', &
                                l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                ONE, a_dev + a_offset, matrixRows,  &
                                v_row_dev + (l_row_beg - 1) * size_of_datatype, 1,  &
                                ONE, u_col_dev + (l_col_beg - 1) * size_of_datatype, 1)
                endif
              enddo

              do i=0,(istep-2)/tile_size
                  l_col_beg = i*l_cols_per_tile+1
                  l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
                  if(l_col_end<l_col_beg) cycle

                  l_row_beg = 1
                  l_row_end = min(l_rows,i*l_rows_per_tile)
                  
                  if (useIntelGPU) then
                    if (wantDebug) call obj%timer%start("mkl_offload")
                    if (wantDebug) call obj%timer%stop("mkl_offload")


                  else
                    a_offset = ((l_row_beg-1) + (l_col_beg - 1) * matrixRows) * &
                            size_of_datatype
                    if (isSkewsymmetric) then
                       call gpublas_CGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   -ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                    else
                       call gpublas_CGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                                   ONE, a_dev + a_offset, matrixRows, &
                                   v_col_dev + (l_col_beg - 1) * size_of_datatype,1, &
                                   ONE, u_row_dev + (l_row_beg - 1) * size_of_datatype, 1)
                   endif
                endif
              enddo
            end if !multiplication as one block / per stripes

            if (.not.(useIntelGPU)) then
              successGPU = gpu_memcpy(int(loc(u_col(1)),kind=c_intptr_t), &
                          u_col_dev, l_cols * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_col_dev 1", 954,  successGPU)

              successGPU = gpu_memcpy(int(loc(u_row(1)),kind=c_intptr_t), &
                          u_row_dev, l_rows * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: u_row_dev 1", 958,  successGPU)
            endif
            !if (useIntelGPU) then
            !  
            !endif

          endif ! useGPU


         ! second calculate (VU**T + UV**T)*v part of (A + VU**T + UV**T)*v
         if (n_stored_vecs > 0) then
           if (wantDebug) call obj%timer%start("blas")
           call CGEMV('C',     &
                               int(l_rows,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, vu_stored_rows, int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                               v_row,  1_BLAS_KIND, ZERO, aux, 1_BLAS_KIND)

           call CGEMV('N', int(l_cols,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),   &
                               ONE, uv_stored_cols, int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),   &
                               aux, 1_BLAS_KIND, ONE, u_col,  1_BLAS_KIND)
           if (wantDebug) call obj%timer%stop("blas")
         endif

       endif  ! (l_rows>0 .and. l_cols>0)

       ! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
       ! on the processors containing the diagonal
       ! This is only necessary if u_row has been calculated, i.e. if the
       ! global tile size is smaller than the global remaining matrix

       if (tile_size < istep-1) then

         call elpa_reduce_add_vectors_&
         &complex&
         &_&
         &single &
         (obj, u_row, ubound(u_row,dim=1), mpi_comm_rows, u_col, ubound(u_col,dim=1), &
         mpi_comm_cols, istep-1, 1, nblk, max_threads)

       endif

       ! Sum up all the u_col(:) parts, transpose u_col -> u_row

       if (l_cols>0) then
         tmp(1:l_cols) = u_col(1:l_cols)
         if (useNonBlockingCollectivesRows) then
           if (wantDebug) call obj%timer%start("mpi_nbc_communication")
           call mpi_iallreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_COMPLEX,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
           call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
         else
           if (wantDebug) call obj%timer%start("mpi_communication")
           call mpi_allreduce(tmp, u_col, int(l_cols,kind=MPI_KIND), MPI_COMPLEX,    &
                            MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
           if (wantDebug) call obj%timer%stop("mpi_communication")
         endif
       endif
       if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
          &complex&
          &_&
          &single &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       else
          call elpa_transpose_vectors_&
          &complex&
          &_&
          &single &
          (obj, u_col, ubound(u_col,dim=1), mpi_comm_cols, u_row, ubound(u_row,dim=1), &
           mpi_comm_rows, 1, istep-1, 1, nblk, max_threads, .false.)
       endif

       ! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )
       x = 0
       if (l_cols>0)  &
       x = dot_product(v_col(1:l_cols),u_col(1:l_cols))

       if (useNonBlockingCollectivesCols) then
         if (wantDebug) call obj%timer%start("mpi_nbc_communication")
         call mpi_iallreduce(x, vav, 1_MPI_KIND, MPI_COMPLEX, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
               allreduce_request3, mpierr)
         call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
       else
         if (wantDebug) call obj%timer%start("mpi_communication")
         call mpi_allreduce(x, vav, 1_MPI_KIND, MPI_COMPLEX, MPI_SUM, int(mpi_comm_cols,kind=MPI_KIND), &
                mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")
       endif

       ! store u and v in the matrices U and V
       ! these matrices are stored combined in one here

       do j=1,l_rows
         vu_stored_rows(j,2*n_stored_vecs+1) = conjg(tau(istep))*v_row(j)
         vu_stored_rows(j,2*n_stored_vecs+2) = 0.5*conjg(tau(istep))*vav*v_row(j) - u_row(j)
       enddo
       do j=1,l_cols
         uv_stored_cols(j,2*n_stored_vecs+1) = 0.5*conjg(tau(istep))*vav*v_col(j) - u_col(j)
         uv_stored_cols(j,2*n_stored_vecs+2) = conjg(tau(istep))*v_col(j)
       enddo

       ! We have calculated another Hauseholder Vector, number of implicitly stored increased
       n_stored_vecs = n_stored_vecs+1

       ! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
       if (n_stored_vecs == max_stored_uv .or. istep == 3) then

         if (useGPU .and. .not.(useIntelGPU)) then
           successGPU = gpu_memcpy(vu_stored_rows_dev, int(loc(vu_stored_rows(1,1)),kind=c_intptr_t), &
                                     max_local_rows * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_rows_dev", 1109,  successGPU)

           successGPU = gpu_memcpy(uv_stored_cols_dev, int(loc(uv_stored_cols(1,1)),kind=c_intptr_t), &
                                     max_local_cols * 2 * max_stored_uv *          &
                                     size_of_datatype, gpuMemcpyHostToDevice)
           call check_memcpy_GPU_f("tridiag: uv_stored_cols_dev", 1114,  successGPU)
         endif
          !if (useIntelGPU) then
          !  ! needed later when we can do explicit offloads
          !endif

         do i = 0, (istep-2)/tile_size
           ! go over tiles above (or on) the diagonal
           l_col_beg = i*l_cols_per_tile+1
           l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
           l_row_beg = 1
           l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
           if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) &
           cycle


           if (useGPU) then
             if (.not. mat_vec_as_one_block) then
               if (useIntelGPU) then
                  if (wantDebug) call obj%timer%start("mkl_offload")

                  call CGEMM('N', 'C',     &
                                       int(l_row_end-l_row_beg+1,kind=BLAS_KIND), &
                                       int(l_col_end-l_col_beg+1,kind=BLAS_KIND), int(2*n_stored_vecs,kind=BLAS_KIND),     &
                                       ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv),                    &
                                       int(max_local_rows,kind=BLAS_KIND), uv_stored_cols(l_col_beg,1),                    &
                                       int(max_local_cols,kind=BLAS_KIND), ONE, a_mat(l_row_beg,l_col_beg),                &
                                       int(matrixRows,kind=BLAS_KIND))
                  if (wantDebug) call obj%timer%stop("mkl_offload")

               else
                 ! if using mat-vec multiply by stripes, it is enough to update tiles above (or on) the diagonal only
                 ! we than use the same calls as for CPU version
                 if (wantDebug) call obj%timer%start("gpublas")
                 call gpublas_CGEMM('N', 'C',     &
                                         l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs,                      &
                                         ONE, vu_stored_rows_dev + (l_row_beg - 1) *                                         &
                                         size_of_datatype,  &
                                         max_local_rows, uv_stored_cols_dev + (l_col_beg - 1) *                              &
                                         size_of_datatype,  &
                                         max_local_cols, ONE, a_dev + ((l_row_beg - 1) + (l_col_beg - 1) * matrixRows) *     &
                                         size_of_datatype , matrixRows)
                 if (wantDebug) call obj%timer%stop("gpublas")
               endif
             endif
           else !useGPU
             if (wantDebug) call obj%timer%start("blas")
             call CGEMM('N', 'C',                &
                                  int(l_row_end-l_row_beg+1,kind=BLAS_KIND), int(l_col_end-l_col_beg+1,kind=BLAS_KIND), &
                                  int(2*n_stored_vecs,kind=BLAS_KIND),    &
                                  ONE, vu_stored_rows(l_row_beg:max_local_rows,1:2*max_stored_uv), &
                                  int(ubound(vu_stored_rows,dim=1),kind=BLAS_KIND),   &
                                  uv_stored_cols(l_col_beg,1), &
                                  int(ubound(uv_stored_cols,dim=1),kind=BLAS_KIND),        &
                                  ONE, a_mat(l_row_beg,l_col_beg), int(matrixRows,kind=BLAS_KIND))
             if (wantDebug) call obj%timer%stop("blas")
           endif !useGPU
         enddo

         if (useGPU) then
           if (mat_vec_as_one_block) then
             if (useIntelGPU) then
                if (wantDebug) call obj%timer%start("mkl_offload")
                call CGEMM('N', 'C', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                                    int(2*n_stored_vecs, kind=BLAS_KIND), ONE,  &
                                    vu_stored_rows, int(max_local_rows,kind=BLAS_KIND), &
                                    uv_stored_cols, int(max_local_cols,kind=BLAS_KIND),  &
                                    ONE, a_mat, int(matrixRows,kind=BLAS_KIND))
                if (wantDebug) call obj%timer%stop("mkl_offload")
             else
               !update whole (remaining) part of matrix, including tiles below diagonal
               !we can do that in one large cublas call
               if (wantDebug) call obj%timer%start("gpublas")
               call gpublas_CGEMM('N', 'C', l_rows, l_cols, 2*n_stored_vecs,   &
                                         ONE, vu_stored_rows_dev, max_local_rows, &
                                         uv_stored_cols_dev, max_local_cols,  &
                                         ONE, a_dev, matrixRows)
               if (wantDebug) call obj%timer%stop("gpublas")
             endif
           endif
         endif

         n_stored_vecs = 0
       endif

       if (my_prow == prow(istep-1, nblk, np_rows) .and. my_pcol == pcol(istep-1, nblk, np_cols)) then
         if (useGPU) then
           if (useIntelGPU) then
                       ! if (useIntelGPU) then
          ! needed at a later time when we can do explcit mem copys
          ! endif

           else
             !a_mat(l_rows,l_cols) = a_dev(l_rows,l_cols)
              a_offset = ((l_rows - 1) + matrixRows * (l_cols - 1)) * size_of_datatype

              successGPU = gpu_memcpy(int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), a_dev + a_offset, &
                                      1 *  size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("tridiag: a_dev 3", 1239,  successGPU)
           endif
         endif
         if (n_stored_vecs > 0) then
           a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) &
                       + dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs),uv_stored_cols(l_cols,1:2*n_stored_vecs))
         end if
         d_vec(istep-1) = real(a_mat(l_rows,l_cols),kind=rk)

         if (useGPU) then
           if (useIntelGPU) then
          ! if (useIntelGPU) then
          ! needed at a later time when we can expicit mem copy
          ! endif
           else
             !a_dev(l_rows,l_cols) = a_mat(l_rows,l_cols)
             !successGPU = gpu_threadsynchronize()
             !call check_memcpy_GPU_f("tridiag: a_dev 4a5a", 1265,  successGPU)

             successGPU = gpu_memcpy(a_dev + a_offset, int(loc(a_mat(l_rows, l_cols)),kind=c_intptr_t), &
                                     int(1 * size_of_datatype, kind=c_intptr_t), gpuMemcpyHostToDevice)
             call check_memcpy_GPU_f("tridiag: a_dev 4", 1269,  successGPU)
           endif
         endif
       endif

     enddo ! main cycle over istep=na,3,-1

     ! Store e_vec(1) and d_vec(1)

     if (my_pcol==pcol(2, nblk, np_cols)) then
      if (my_prow==prow(1, nblk, np_rows)) then
       ! We use last l_cols value of loop above
       if (useGPU) then
         if (useIntelGPU) then
            vrl = a_mat(1,l_cols)
         else
           successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev + (matrixRows * (l_cols - 1)) * size_of_datatype, &
                                   1 * size_of_datatype, gpuMemcpyDeviceToHost)
           call check_memcpy_GPU_f("tridiag: a_dev 5", 1288,  successGPU)
           vrl = aux3(1)
         endif
       else !useGPU
         vrl = a_mat(1,l_cols)
       endif !useGPU
       call hh_transform_complex_&
       &single &
       (obj, vrl, 0.0_rk, xf, tau(2), wantDebug)
       e_vec(1) = real(vrl,kind=rk)
       a_mat(1,l_cols) = 1. ! for consistency only
     endif
     if (useNonBlockingCollectivesRows) then
       if (wantDebug) call obj%timer%start("mpi_nbc_communication")
       call mpi_ibcast(tau(2), 1_MPI_KIND, MPI_COMPLEX, int(prow(1, nblk, np_rows),kind=MPI_KIND), &
                   int(mpi_comm_rows,kind=MPI_KIND), bcast_request2, mpierr)
       call mpi_wait(bcast_request2, MPI_STATUS_IGNORE, mpierr)
       if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
     else
       if (wantDebug) call obj%timer%start("mpi_communication")
       call mpi_bcast(tau(2), 1_MPI_KIND, MPI_COMPLEX, int(prow(1, nblk, np_rows),kind=MPI_KIND), &
                   int(mpi_comm_rows,kind=MPI_KIND),  mpierr)
       if (wantDebug) call obj%timer%stop("mpi_communication")
     endif

   endif

   if (useNonBlockingCollectivesCols) then
     if (wantDebug) call obj%timer%start("mpi_nbc_communication")
     call mpi_ibcast(tau(2), 1_MPI_KIND, MPI_COMPLEX, int(pcol(2, nblk, np_cols),kind=MPI_KIND), &
                  int(mpi_comm_cols,kind=MPI_KIND), bcast_request3, mpierr)
     call mpi_wait(bcast_request3, MPI_STATUS_IGNORE, mpierr)
     if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
   else
     if (wantDebug) call obj%timer%start("mpi_communication")
     call mpi_bcast(tau(2), 1_MPI_KIND, MPI_COMPLEX, int(pcol(2, nblk, np_cols),kind=MPI_KIND), &
                  int(mpi_comm_cols,kind=MPI_KIND), mpierr)
     if (wantDebug) call obj%timer%stop("mpi_communication")
   endif

  if (my_prow == prow(1, nblk, np_rows) .and. my_pcol == pcol(1, nblk, np_cols))  then
    if (useGPU) then
      if (useIntelGPU) then
        d_vec(1) = REAL(a_mat(1,1))
      else
        successGPU = gpu_memcpy(int(loc(aux3(1)),kind=c_intptr_t), a_dev, &
                               1 * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("tridiag: a_dev 6", 1344,  successGPU)
        d_vec(1) = REAL(aux3(1))
      endif
    else !useGPU
      d_vec(1) = REAL(a_mat(1,1))
    endif !useGPU
  endif



  deallocate(tmp, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp", 1395,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: should we leave a_mat on the device for further use?
    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("tridiag: a_dev 9", 1400,  successGPU)

    successGPU = gpu_free(v_row_dev)
    call check_dealloc_GPU_f("tridiag: v_row_dev", 1403,  successGPU)

    successGPU = gpu_free(u_row_dev)
    call check_dealloc_GPU_f("tridiag: (u_row_dev", 1406,  successGPU)

    successGPU = gpu_free(v_col_dev)
    call check_dealloc_GPU_f("tridiag: v_col_dev", 1409,  successGPU)

    successGPU = gpu_free(u_col_dev)
    call check_dealloc_GPU_f("tridiag: u_col_dev ", 1412,  successGPU)

    successGPU = gpu_free(vu_stored_rows_dev)
    call check_dealloc_GPU_f("tridiag: vu_stored_rows_dev ", 1415,  successGPU)

    successGPU = gpu_free(uv_stored_cols_dev)
    call check_dealloc_GPU_f("tridiag:uv_stored_cols_dev ", 1418,  successGPU)
  endif
  ! if (useIntelGPU) then
  ! needed at a later time when we can do explicit frees
  ! endif


  ! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag: tmp_real", 1428,  istat,  errorMessage)

  if (useNonBlockingCollectivesRows) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request4, mpierr)
    call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6,mpierr)
    call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif
  if (useNonBlockingCollectivesCols) then
    if (wantDebug) call obj%timer%start("mpi_nbc_communication")
    tmp_real = d_vec
    call mpi_iallreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request5, mpierr)
    call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

    tmp_real = e_vec
    call mpi_iallreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), allreduce_request7, mpierr)
    call mpi_wait(allreduce_request7, MPI_STATUS_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
  else
    if (wantDebug) call obj%timer%start("mpi_communication")
    tmp_real = d_vec
    call mpi_allreduce(tmp_real, d_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    tmp_real = e_vec
    call mpi_allreduce(tmp_real, e_vec, int(na,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_cols,kind=MPI_KIND), mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif

  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: tmp_real", 1477,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
           deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
     call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1482,  istat,  errorMessage)
    else
      successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: a_mat", 1485,  successGPU)

      successGPU = gpu_free_host(v_row_host)
      call check_host_dealloc_GPU_f("tridiag: v_row_host", 1488,  successGPU)
      nullify(v_row)

      successGPU = gpu_free_host(v_col_host)
      call check_host_dealloc_GPU_f("tridiag: v_col_host", 1492,  successGPU)
      nullify(v_col)

      successGPU = gpu_free_host(u_col_host)
      call check_host_dealloc_GPU_f("tridiag: u_col_host", 1496,  successGPU)
      nullify(u_col)

      successGPU = gpu_free_host(u_row_host)
      call check_host_dealloc_GPU_f("tridiag: u_row_host", 1500,  successGPU)
      nullify(u_row)

      successGPU = gpu_host_unregister(int(loc(uv_stored_cols),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: uv_stored_cols", 1504,  successGPU)

      successGPU = gpu_host_unregister(int(loc(vu_stored_rows),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: vu_stored_rows", 1507,  successGPU)

      successGPU = gpu_host_unregister(int(loc(e_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: e_vec", 1510,  successGPU)

      successGPU = gpu_host_unregister(int(loc(d_vec),kind=c_intptr_t))
      call check_host_unregister_GPU_f("tridiag: d_vec", 1513,  successGPU)
    endif
  else
    deallocate(v_row, v_col, u_row, u_col, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("tridiag: v_row, v_col, u_row, u_col", 1517,  istat,  errorMessage)
  endif

  deallocate(vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag: vu_stored_rows, uv_stored_cols", 1521,  istat,  errorMessage)

  call obj%timer%stop("tridiag_&
  &complex&
  &" // &
  "_single" // &
  gpuString )

end subroutine tridiag_&
&complex&
&_&
&single





!> \brief Transforms the eigenvectors of a tridiagonal matrix back
!>                     to the eigenvectors of the original matrix
!>                     (like Scalapack Routine PDORMTR)
!>
!  Parameters
!
!> \param na          Order of matrix a_mat, number of rows of matrix q_mat
!>
!> \param nqc         Number of columns of matrix q_mat
!>
!> \param a_mat(lda,matrixCols)  Matrix containing the Householder vectors (i.e. matrix a after tridiag_real)
!>                           Distribution is like in Scalapack.
!>
!> \param lda         Leading dimension of a_mat
!>
!> \param tau(na)     Factors of the Householder vectors
!>
!> \param q_mat           On input: Eigenvectors of tridiagonal matrix
!>                    On output: Transformed eigenvectors
!>                    Distribution is like in Scalapack.
!>
!> \param ldq         Leading dimension of q_mat
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
!>
!> \param matrixCols  local columns of matrix a_mat and q_mat
!>
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
!> \param useGPU      If true,  GPU version of the subroutine will be used
!>

subroutine trans_ev_&
&complex&
&_&
&single &
(obj, na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, useGPU)
  use, intrinsic :: iso_c_binding
  use precision
  use elpa_abstract_impl
  use elpa_blas_interfaces
  use elpa_gpu

  implicit none
!    Copyright 2011, A. Marek
!
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

  integer, parameter :: rk = C_FLOAT
  integer, parameter :: ck = C_FLOAT_COMPLEX
  integer, parameter :: rck = C_FLOAT_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  integer(kind=ik), intent(in)                  :: na, nqc, lda, ldq, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
  complex(kind=rck), intent(in)           :: tau(na)

  complex(kind=rck), intent(inout)        :: a_mat(lda,*)
  complex(kind=rck), intent(inout)        :: q_mat(ldq,*)
  logical, intent(in)                           :: useGPU
  integer(kind=ik)                              :: max_stored_rows, max_stored_rows_fac

  integer(kind=ik)                              :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                        :: mpierr, my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                              :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik)                              :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik)                              :: istep, n, nc, ic, ics, ice, nb, cur_pcol
  integer(kind=ik)                              :: hvn_ubnd, hvm_ubnd

  complex(kind=rck), allocatable          :: hvb(:), hvm(:,:)
  complex(kind=rck), pointer              :: tmp1(:), tmp2(:)
  complex(kind=rck), allocatable          :: h1(:), h2(:)
  complex(kind=rck), pointer              :: tmat(:,:)
  complex(kind=rck), pointer              :: hvm1(:)
  type(c_ptr)                                   :: tmp1_host, tmp2_host
  type(c_ptr)                                   :: hvm1_host, tmat_host

  integer(kind=ik)                              :: istat
  character(200)                                :: errorMessage
  character(20)                                 :: gpuString

  integer(kind=c_intptr_t)                      :: num
  integer(kind=C_intptr_T)                      :: q_dev, tmp_dev, hvm_dev, tmat_dev

  integer(kind=ik)                              :: blockStep
  logical                                       :: successGPU
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_&
                                                                      &complex
  integer(kind=ik)                              :: error
  logical                                       :: useIntelGPU
  integer(kind=MPI_KIND)                        :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                       :: useNonBlockingCollectivesCols
  logical                                       :: useNonBlockingCollectivesRows
  integer(kind=c_int)                           :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

  call obj%get("nbc_row_elpa1_tridi_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa1_tridi_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa1_tridi_to_full. Aborting..."
    stop
  endif

  if (non_blocking_collectives_rows .eq. 1) then
    useNonBlockingCollectivesRows = .true.
  else
    useNonBlockingCollectivesRows = .false.
  endif

  if (non_blocking_collectives_cols .eq. 1) then
    useNonBlockingCollectivesCols = .true.
  else
    useNonBlockingCollectivesCols = .false.
  endif

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

  my_prow = int(my_prowMPI, kind=c_int)
  np_rows = int(np_rowsMPI, kind=c_int)
  my_pcol = int(my_pcolMPI, kind=c_int)
  np_cols = int(np_colsMPI, kind=c_int)
  call obj%timer%stop("mpi_communication")

  call obj%get("max_stored_rows",max_stored_rows_fac, error)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (max_stored_rows_fac/nblk+1)*nblk

  if (useIntelGPU) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp2", istat, errorMessage)
  endif

  if (.not.(useGPU)) then
    allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmat", istat, errorMessage)

    allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp1", istat, errorMessage)

    allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
    call check_alloc("trans_ev_&
    &complex&
    &", "tmp2", istat, errorMessage)
  endif

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "h2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_&
  &complex&
  &", "hvm", istat, errorMessage)

  hvm = 0   ! Must be set to 0 !!!
  hvb = 0   ! Safety only
  blockStep = nblk

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0
  if (useGPU) then
    hvn_ubnd = 0
  endif

  ! In the complex case tau(2) /= 0
  if (my_prow == prow(1, nblk, np_rows)) then
    q_mat(1,1:l_cols) = q_mat(1,1:l_cols)*(ONE-tau(2))
  endif

  if (useGPU .and. .not.(useIntelGPU)) then
    ! todo: this is used only for copying hmv to device.. it should be possible to go without it
    !allocate(hvm1(max_local_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
    !call check_alloc("trans_ev_&
    !&complex&
    !&", "hvm1", istat, errorMessage)
    num = (max_local_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(hvm1_host,num)
    call check_alloc_GPU_f("trans_ev: hvm1_host", 304,  successGPU)
    call c_f_pointer(hvm1_host,hvm1,(/(max_local_rows*max_stored_rows)/))

    num = (max_stored_rows*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmat_host,num)
    call check_alloc_GPU_f("trans_ev: tmat_host", 309,  successGPU)
    call c_f_pointer(tmat_host,tmat,(/max_stored_rows,max_stored_rows/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_alloc_GPU_f("trans_ev: tmp1_host", 314,  successGPU)
    call c_f_pointer(tmp1_host,tmp1,(/(max_local_cols*max_stored_rows)/))

    num = (max_local_cols*max_stored_rows) * size_of_datatype
    successGPU = gpu_malloc_host(tmp2_host,num)
    call check_alloc_GPU_f("trans_ev: tmp2_host", 319,  successGPU)
    call c_f_pointer(tmp2_host,tmp2,(/(max_local_cols*max_stored_rows)/))

    successGPU = gpu_malloc(tmat_dev, max_stored_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 323,  successGPU)

    successGPU = gpu_malloc(hvm_dev, max_local_rows * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 326,  successGPU)

    successGPU = gpu_malloc(tmp_dev, max_local_cols * max_stored_rows * size_of_datatype)
    call check_alloc_GPU_f("trans_ev", 329,  successGPU)

    num = ldq * matrixCols * size_of_datatype
    successGPU = gpu_malloc(q_dev, num)
    call check_alloc_GPU_f("trans_ev", 333,  successGPU)

    successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev: q_mat", 337,  successGPU)

    successGPU = gpu_memcpy(q_dev, int(loc(q_mat(1,1)),kind=c_intptr_t), &
                  num, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("trans_ev", 341,  successGPU)
  endif  ! useGPU

  !if (useIntelGPU) then
  !  ! needed at a later time when we can do explicit copys
  !endif


  do istep = 1, na, blockStep
    ics = MAX(istep,3)
    ice = MIN(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic = ics, ice

      l_colh = local_index(ic  , my_pcol, np_cols, nblk, -1) ! Column of Householder Vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector


      if (my_pcol == cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow == prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo

    if (nb>0) then
      if (useNonBlockingCollectivesCols) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(hvb, int(nb,kind=MPI_KIND), MPI_COMPLEX , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
       else
        call obj%timer%start("mpi_communication")
        call mpi_bcast(hvb, int(nb,kind=MPI_KIND), MPI_COMPLEX , int(cur_pcol,kind=MPI_KIND), &
                   int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif
    endif

    nb = 0
    do ic = ics, ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder Vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      if (useGPU) then
        hvm_ubnd = l_rows
      endif
      nstor = nstor+1
      nb = nb+l_rows
    enddo

    ! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk > max_stored_rows .or. istep+nblk > na .or. (na/np_rows <= 256 .and. nstor >= 32)) then

      ! Calculate scalar products of stored vectors.
      ! This can be done in different ways, we use dsyrk or zherk

      tmat = 0
      call obj%timer%start("blas")
      if (l_rows>0) &
      call CHERK('U', 'C',   &
                         int(nstor,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                         hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), ZERO, tmat, int(max_stored_rows,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      nc = 0
      do n = 1, nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo
      if (nc>0) then
        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce( h1, h2, int(nc,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                                   int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif
      endif
      ! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n = 1, nstor-1
        call obj%timer%start("blas")
        call CTRMV('L', 'C' , 'N', int(n,kind=BLAS_KIND), tmat, &
                            int(max_stored_rows,kind=BLAS_KIND), h2(nc+1), 1_BLAS_KIND)
        call obj%timer%stop("blas")

        tmat(n+1,1:n) = &
        -conjg(h2(nc+1:nc+n)) &
        *tau(ice-nstor+n+1)

        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

      if (useGPU .and. .not.(useIntelGPU)) then
        ! todo: is this reshape really neccessary?
        hvm1(1:hvm_ubnd*nstor) = reshape(hvm(1:hvm_ubnd,1:nstor), (/ hvm_ubnd*nstor /))

        !hvm_dev(1:hvm_ubnd*nstor) = hvm1(1:hvm_ubnd*nstor)
        successGPU = gpu_memcpy(hvm_dev, int(loc(hvm1(1)),kind=c_intptr_t),   &
                      hvm_ubnd * nstor * size_of_datatype, gpuMemcpyHostToDevice)

        call check_memcpy_GPU_f("trans_ev", 475,  successGPU)

        !tmat_dev = tmat
        successGPU = gpu_memcpy(tmat_dev, int(loc(tmat(1,1)),kind=c_intptr_t),   &
                      max_stored_rows * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 480,  successGPU)
      endif

      !if (useIntelGPU) then
      !  ! needed later when we can do explicit copys
      !endif


      ! Q = Q - V * T * V**T * Q

      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            call obj%timer%start("mkl_offload")
            call obj%timer%stop("mkl_offload")

          else
            call obj%timer%start("gpublas")
            call gpublas_CGEMM('C', 'N',   &
                                     nstor, l_cols, l_rows, ONE, hvm_dev, hvm_ubnd,  &
                                     q_dev, ldq, ZERO, tmp_dev, nstor)
            call obj%timer%stop("gpublas")
          endif
        else ! useGPU

          call obj%timer%start("blas")
          call CGEMM('C', 'N',  &
                              int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(l_rows,kind=BLAS_KIND), ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), &
                              q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, int(nstor,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU

      else !l_rows>0

        if (useGPU) then
          if (useIntelGPU) then
            tmp1(1:l_cols*nstor) = 0
          else
            successGPU = gpu_memset(tmp_dev, 0, l_cols * nstor * size_of_datatype)
            call check_memcpy_GPU_f("trans_ev", 532,  successGPU)
          endif
        else
          tmp1(1:l_cols*nstor) = 0
        endif
      endif  !l_rows>0


      if (useGPU .and. .not.(useIntelGPU)) then
        ! In the legacy GPU version, this allreduce was ommited. But probably it has to be done for GPU + MPI
        ! todo: does it need to be copied whole? Wouldn't be a part sufficient?
        successGPU = gpu_memcpy(int(loc(tmp1(1)),kind=c_intptr_t), tmp_dev,  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("trans_ev", 547,  successGPU)
      endif

      !if (useIntelGPU) then
      !   ! needed later
      !endif

      if (useNonBlockingCollectivesRows) then
        call obj%timer%start("mpi_nbc_communication")
        call mpi_iallreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
        call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
        call obj%timer%stop("mpi_nbc_communication")
      else
        call obj%timer%start("mpi_communication")
        call mpi_allreduce(tmp1, tmp2, int(nstor*l_cols,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                         int(mpi_comm_rows,kind=MPI_KIND), mpierr)
        call obj%timer%stop("mpi_communication")
      endif

      if (useGPU .and. .not.(useIntelGPU)) then
        ! copy back tmp2 - after reduction...
        successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2(1)),kind=c_intptr_t),  &
                      max_local_cols * max_stored_rows * size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev", 590,  successGPU)
      endif ! useGPU


      !if (useIntelGPU) then
      !   ! needed later
      !endif


      if (l_rows>0) then
        if (useGPU) then
          if (useIntelGPU) then
            ! tmp2 = tmat * tmp2
            call obj%timer%start("mkl_offload")

            call obj%timer%stop("mkl_offload")
          else
            call obj%timer%start("gpublas")
            call gpublas_CTRMM('L', 'L', 'N', 'N',     &
                                     nstor, l_cols, ONE, tmat_dev, max_stored_rows,  &
                                     tmp_dev, nstor)

            call gpublas_CGEMM('N', 'N' ,l_rows ,l_cols ,nstor,  &
                                     -ONE, hvm_dev, hvm_ubnd, tmp_dev, nstor,   &
                                     ONE, q_dev, ldq)
            call obj%timer%stop("gpublas")
          endif
        else !useGPU
          ! tmp2 = tmat * tmp2
          call obj%timer%start("blas")
          call CTRMM('L', 'L', 'N', 'N', int(nstor,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND),   &
                             ONE, tmat, int(max_stored_rows,kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND))
          !q_mat = q_mat - hvm*tmp2
          call CGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(nstor,kind=BLAS_KIND),   &
                              -ONE, hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(nstor,kind=BLAS_KIND), &
                              ONE, q_mat, int(ldq,kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU
      endif  ! l_rows>0
      nstor = 0
    endif  ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

  enddo ! istep

  deallocate(h1, h2, hvb, hvm, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: h1, h2, hvb, hvm", 693,  istat,  errorMessage)

  if (useGPU) then

    if (useIntelGPU) then
      deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_&       &MATH_DATATYPE&       &: tmat, tmp1, tmp2", 701,  istat,  errorMessage)
    else

      !q_mat = q_dev
      successGPU = gpu_memcpy(int(loc(q_mat(1,1)),kind=c_intptr_t), &
                    q_dev, ldq * matrixCols * size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev", 707,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev: q_mat", 710,  successGPU)

      successGPU = gpu_free_host(hvm1_host)
      call check_host_dealloc_GPU_f("trans_ev: hvm1_host", 713,  successGPU)
      nullify(hvm1)

      successGPU = gpu_free_host(tmat_host)
      call check_host_dealloc_GPU_f("trans_ev: tmat_host", 717,  successGPU)
      nullify(tmat)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp1_host", 721,  successGPU)
      nullify(tmp1)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev: tmp2_host", 725,  successGPU)
      nullify(tmp2)

      !deallocate(hvm1, stat=istat, errmsg=errorMessage)
      !if (istat .ne. 0) then
      !  print *,"trans_ev_&
      !  &complex&
      !  &: error when deallocating hvm1 "//errorMessage
      !  stop 1
      !endif

      !deallocate(q_dev, tmp_dev, hvm_dev, tmat_dev)
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev", 738,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev", 741,  successGPU)

      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev", 744,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev", 747,  successGPU)
    endif
  else
    deallocate(tmat, tmp1, tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_&     &MATH_DATATYPE&     &: tmat, tmp1, tmp2", 753,  istat,  errorMessage)
  endif


  call obj%timer%stop("trans_ev_&
  &complex&
  &" // &
  &"_single" // &
  gpuString )

end subroutine trans_ev_&
&complex&
&_&
&single







subroutine hh_transform_complex_&
&single &
(obj, alpha, xnorm_sq, xf, tau, wantDebug)
  ! Similar to LAPACK routine ZLARFP, but uses ||x||**2 instead of x(:)
  ! and returns the factor xf by which x has to be scaled.
  ! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
  ! since this would be expensive for the parallel implementation.
  use precision
  use elpa_abstract_impl
  implicit none
!    Copyright 2011, A. Marek
!
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

  integer, parameter :: rk = C_FLOAT
  integer, parameter :: ck = C_FLOAT_COMPLEX
  integer, parameter :: rck = C_FLOAT_COMPLEX
  complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)
  class(elpa_abstract_impl_t), intent(inout)    :: obj
  logical, intent(in)                           :: wantDebug
  complex(kind=ck), intent(inout) :: alpha
  real(kind=rk), intent(in)          :: xnorm_sq
  complex(kind=ck), intent(out)   :: xf, tau
  real(kind=rk)                      :: ALPHR, ALPHI

  real(kind=rk)                      :: BETA

  if (wantDebug) call obj%timer%start("hh_transform_&
                   &complex&
     	      &" // &
                   &"_single" )

  ALPHR = real( ALPHA, kind=rk )
  ALPHI = AIMAG( ALPHA )

  if ( XNORM_SQ==0.0_rk .AND. ALPHI==0.0_rk ) then

    if ( ALPHR>=0.0_rk ) then
      TAU = 0.0_rk
    else
      TAU = 2.0_rk
      ALPHA = -ALPHA
    endif
    XF = 0.0_rk

  else

    BETA = SIGN( SQRT( ALPHR**2 + ALPHI**2 + XNORM_SQ ), ALPHR )
    ALPHA = ALPHA + BETA
    IF ( BETA<0 ) THEN
      BETA = -BETA
      TAU  = -ALPHA / BETA
    ELSE
      ALPHR = ALPHI * (ALPHI/real( ALPHA , kind=rk))
      ALPHR = ALPHR + XNORM_SQ/real( ALPHA, kind=rk )

      TAU = CMPLX( ALPHR/BETA, -ALPHI/BETA )
      ALPHA = CMPLX( -ALPHR, ALPHI )
    END IF
    XF = 1.0_rk/ALPHA
    ALPHA = BETA
  endif

  if (wantDebug) call obj%timer%stop("hh_transform_&
  &complex&
  &" // &
  &"_single" )

    end subroutine hh_transform_complex_&
    &single





end module elpa1_compute
