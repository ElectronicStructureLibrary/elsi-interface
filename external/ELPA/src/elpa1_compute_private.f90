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
module ELPA1_COMPUTE

  use elpa_utilities
  use elpa_mpi

  implicit none

  private

  public :: tridiag_real_double
  public :: tridiag_real
  public :: trans_ev_real_double
  public :: trans_ev_real

  public :: solve_tridi_double

  interface tridiag_real
    module procedure tridiag_real_double
  end interface

  interface trans_ev_real
    module procedure trans_ev_real_double
  end interface

  public :: tridiag_complex_double
  public :: tridiag_complex
  public :: trans_ev_complex_double
  public :: trans_ev_complex

  interface tridiag_complex
    module procedure tridiag_complex_double
  end interface

  interface trans_ev_complex
    module procedure trans_ev_complex_double
  end interface

  public :: hh_transform_real_double
  public :: hh_transform_real
  public :: elpa_reduce_add_vectors_real_double
  public :: elpa_reduce_add_vectors_real
  public :: elpa_transpose_vectors_real_double
  public :: elpa_transpose_vectors_real

  interface hh_transform_real
    module procedure hh_transform_real_double
  end interface

  interface elpa_reduce_add_vectors_real
    module procedure elpa_reduce_add_vectors_real_double
  end interface

  interface elpa_transpose_vectors_real
    module procedure elpa_transpose_vectors_real_double
  end interface

  public :: hh_transform_complex_double
  public :: hh_transform_complex
  public :: elpa_reduce_add_vectors_complex_double
  public :: elpa_reduce_add_vectors_complex
  public :: elpa_transpose_vectors_complex_double
  public :: elpa_transpose_vectors_complex

  interface hh_transform_complex
    module procedure hh_transform_complex_double
  end interface

  interface elpa_reduce_add_vectors_complex
    module procedure elpa_reduce_add_vectors_complex_double
  end interface

  interface elpa_transpose_vectors_complex
    module procedure elpa_transpose_vectors_complex_double
  end interface

  contains

! real double precision first

subroutine elpa_transpose_vectors_real_double(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvs,nvr,nvc,nblk)

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
  use elpa_mpi

  implicit none

  integer(kind=ik), intent(in) :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  real(kind=rk8), intent(in) :: vmat_s(ld_s,nvc)
  real(kind=rk8), intent(inout) :: vmat_t(ld_t,nvc)

  real(kind=rk8), allocatable :: aux(:)
  integer(kind=ik) :: myps, mypt, nps, npt
  integer(kind=ik) :: n, lc, k, i, ips, ipt, ns, nl, mpierr
  integer(kind=ik) :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik) :: auxstride

  call MPI_Comm_rank(comm_s,myps,mpierr)
  call MPI_Comm_size(comm_s,nps ,mpierr)
  call MPI_Comm_rank(comm_t,mypt,mpierr)
  call MPI_Comm_size(comm_t,npt ,mpierr)

! The basic idea of this routine is that for every block (in the block cyclic
! distribution), the processor within comm_t which owns the diagonal
! broadcasts its values of vmat_s to all processors within comm_t.
! Of course this has not to be done for every block separately, since
! the communictation pattern repeats in the global matrix after
! the least common multiple of (nps,npt) blocks

  lcm_s_t = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

! Get the number of blocks to be skipped at the begin.
! This must be a multiple of lcm_s_t (else it is getting complicated),
! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))

  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if(mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
      if (nblks_comm .ne. 0) then
        if(myps == ips) then

          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
               k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
               ns = (i/nps)*nblk ! local start of block i
               nl = min(nvr-i*nblk,nblk) ! length
               aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
            enddo
          enddo
        endif

        call MPI_Bcast(aux, nblks_comm*nblk*nvc, MPI_REAL8, ips, comm_s, mpierr)

        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
             k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
             ns = (i/npt)*nblk ! local start of block i
             nl = min(nvr-i*nblk,nblk) ! length
             vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
          enddo
        enddo
      endif
    endif

  enddo

  deallocate(aux)

end subroutine

subroutine elpa_reduce_add_vectors_real_double(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvr,nvc,nblk)

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

  implicit none

  integer(kind=ik), intent(in) :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  real(kind=rk8), intent(in) :: vmat_s(ld_s,nvc)
  real(kind=rk8), intent(inout) :: vmat_t(ld_t,nvc)

  real(kind=rk8), allocatable :: aux1(:), aux2(:)
  integer(kind=ik) :: myps, mypt, nps, npt
  integer(kind=ik) :: n, lc, k, i, ips, ipt, ns, nl, mpierr
  integer(kind=ik) :: lcm_s_t, nblks_tot
  integer(kind=ik) :: auxstride

  call MPI_Comm_rank(comm_s,myps,mpierr)
  call MPI_Comm_size(comm_s,nps ,mpierr)
  call MPI_Comm_rank(comm_t,mypt,mpierr)
  call MPI_Comm_size(comm_t,npt ,mpierr)

! Look to elpa_transpose_vectors for the basic idea!

! The communictation pattern repeats in the global matrix after
! the least common multiple of (nps,npt) blocks

  lcm_s_t = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
  aux1(:) = 0
  aux2(:) = 0

  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    auxstride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if(myps == ips) then
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
        enddo
      enddo

      k = nvc * auxstride

      if(k>0) call MPI_Reduce(aux1, aux2, k, MPI_REAL8, MPI_SUM, ipt, comm_t, mpierr)

      if (mypt == ipt) then
        do lc=1,nvc
          do i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
          enddo
        enddo
      endif

    endif

  enddo

  deallocate(aux1)
  deallocate(aux2)

end subroutine

! double precision

subroutine elpa_transpose_vectors_complex_double(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvs,nvr,nvc,nblk)

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
  use elpa_mpi

  implicit none

  integer(kind=ik), intent(in) :: ld_s, comm_s, ld_t, comm_t, nvs, nvr, nvc, nblk
  complex(kind=ck8), intent(in) :: vmat_s(ld_s,nvc)
  complex(kind=ck8), intent(inout) :: vmat_t(ld_t,nvc)

  complex(kind=ck8), allocatable :: aux(:)
  integer(kind=ik) :: myps, mypt, nps, npt
  integer(kind=ik) :: n, lc, k, i, ips, ipt, ns, nl, mpierr
  integer(kind=ik) :: lcm_s_t, nblks_tot, nblks_comm, nblks_skip
  integer(kind=ik) :: auxstride

  call MPI_Comm_rank(comm_s,myps,mpierr)
  call MPI_Comm_size(comm_s,nps ,mpierr)
  call MPI_Comm_rank(comm_t,mypt,mpierr)
  call MPI_Comm_size(comm_t,npt ,mpierr)

! The basic idea of this routine is that for every block (in the block cyclic
! distribution), the processor within comm_t which owns the diagonal
! broadcasts its values of vmat_s to all processors within comm_t.
! Of course this has not to be done for every block separately, since
! the communictation pattern repeats in the global matrix after
! the least common multiple of (nps,npt) blocks

  lcm_s_t = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

! Get the number of blocks to be skipped at the begin.
! This must be a multiple of lcm_s_t (else it is getting complicated),
! thus some elements before nvs will be accessed/set.

  nblks_skip = ((nvs-1)/(nblk*lcm_s_t))*lcm_s_t

  allocate(aux( ((nblks_tot-nblks_skip+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))

  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    if(mypt == ipt) then

      nblks_comm = (nblks_tot-nblks_skip-n+lcm_s_t-1)/lcm_s_t
      auxstride = nblk * nblks_comm
      if (nblks_comm .ne. 0) then
        if(myps == ips) then
          do lc=1,nvc
            do i = nblks_skip+n, nblks_tot-1, lcm_s_t
              k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
              ns = (i/nps)*nblk ! local start of block i
              nl = min(nvr-i*nblk,nblk) ! length
              aux(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
            enddo
          enddo
        endif

        call MPI_Bcast(aux, nblks_comm*nblk*nvc, MPI_DOUBLE_COMPLEX, ips, comm_s, mpierr)

        do lc=1,nvc
          do i = nblks_skip+n, nblks_tot-1, lcm_s_t
            k = (i - nblks_skip - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = aux(k+1:k+nl)
          enddo
        enddo
      endif
    endif
  enddo

  deallocate(aux)

end subroutine

subroutine elpa_reduce_add_vectors_complex_double(vmat_s,ld_s,comm_s,vmat_t,ld_t,comm_t,nvr,nvc,nblk)

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

  implicit none

  integer(kind=ik), intent(in) :: ld_s, comm_s, ld_t, comm_t, nvr, nvc, nblk
  complex(kind=ck8), intent(in) :: vmat_s(ld_s,nvc)
  complex(kind=ck8), intent(inout) :: vmat_t(ld_t,nvc)

  complex(kind=ck8), allocatable :: aux1(:), aux2(:)
  integer(kind=ik) :: myps, mypt, nps, npt
  integer(kind=ik) :: n, lc, k, i, ips, ipt, ns, nl, mpierr
  integer(kind=ik) :: lcm_s_t, nblks_tot
  integer(kind=ik) :: auxstride

  call MPI_Comm_rank(comm_s,myps,mpierr)
  call MPI_Comm_size(comm_s,nps ,mpierr)
  call MPI_Comm_rank(comm_t,mypt,mpierr)
  call MPI_Comm_size(comm_t,npt ,mpierr)

! Look to elpa_transpose_vectors for the basic idea!

! The communictation pattern repeats in the global matrix after
! the least common multiple of (nps,npt) blocks

  lcm_s_t = least_common_multiple(nps,npt) ! least common multiple of nps, npt

  nblks_tot = (nvr+nblk-1)/nblk ! number of blocks corresponding to nvr

  allocate(aux1( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
  allocate(aux2( ((nblks_tot+lcm_s_t-1)/lcm_s_t) * nblk * nvc ))
  aux1(:) = 0
  aux2(:) = 0

  do n = 0, lcm_s_t-1

    ips = mod(n,nps)
    ipt = mod(n,npt)

    auxstride = nblk * ((nblks_tot - n + lcm_s_t - 1)/lcm_s_t)

    if(myps == ips) then
      do lc=1,nvc
        do i = n, nblks_tot-1, lcm_s_t
          k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
          ns = (i/nps)*nblk ! local start of block i
          nl = min(nvr-i*nblk,nblk) ! length
          aux1(k+1:k+nl) = vmat_s(ns+1:ns+nl,lc)
        enddo
      enddo

      k = nvc * auxstride

      if(k>0) call MPI_Reduce(aux1, aux2, k, MPI_DOUBLE_COMPLEX, MPI_SUM, ipt, comm_t, mpierr)

      if (mypt == ipt) then
        do lc=1,nvc
          do i = n, nblks_tot-1, lcm_s_t
            k = (i - n)/lcm_s_t * nblk + (lc - 1) * auxstride
            ns = (i/npt)*nblk ! local start of block i
            nl = min(nvr-i*nblk,nblk) ! length
            vmat_t(ns+1:ns+nl,lc) = vmat_t(ns+1:ns+nl,lc) + aux2(k+1:k+nl)
          enddo
        enddo
      endif
    endif
  enddo

  deallocate(aux1)
  deallocate(aux2)

end subroutine

! real double precision

!cannot use "elpa1_compute_real_template.X90" because filename with path can be too long for gfortran (max line length)

!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param na          Order of matrix
!>
!> \param a_mat(lda,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param lda         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
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
subroutine tridiag_real_double(na, a_mat, lda, nblk, mpi_comm_rows, mpi_comm_cols, d_vec, e_vec, tau)

  use, intrinsic :: iso_c_binding
  use precision

  implicit none

  integer(kind=ik), intent(in) :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols

  real(kind=rk8), intent(out) :: d_vec(na), e_vec(na), tau(na)

  real(kind=rk8), intent(inout) :: a_mat(lda,*)

  integer(kind=ik), parameter :: max_stored_uv = 32

! id in processor row and column and total numbers of processor rows and columns
  integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=ik) :: mpierr

  integer(kind=ik) :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, max_local_cols

! updated after each istep (in the main cycle) to contain number of
! local columns and rows of the remaining part of the matrix
!integer(kind=ik) :: l_cols, l_rows
  integer(kind=ik) :: l_cols, l_rows

  integer(kind=ik) :: n_stored_vecs
  integer(kind=ik) :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik) :: tile_size, l_rows_per_tile, l_cols_per_tile

  real(kind=rk8) :: vav, vnorm2, x, aux(2*max_stored_uv), aux1(2), aux2(2), vrl, xf

  real(kind=rk8), allocatable :: tmp(:)
  real(kind=rk8), allocatable :: v_row(:) ! used to store calculated Householder vector
  real(kind=rk8), allocatable :: v_col(:) ! the same vector, but transposed - differently distributed among MPI tasks
  real(kind=rk8), allocatable :: u_row(:)
  real(kind=rk8), allocatable :: u_col(:)

! the following two matrices store pairs of vectors v and u calculated in each step
! at most max_stored_uv vector pairs are stored, than the matrix A_i is explicitli updated
! u and v are stored both in row and vector forms
! pattern: v1,u1,v2,u2,v3,u3,....
! todo: It is little bit confusing, I think, that variables _row actually store columns and vice versa
  real(kind=rk8), allocatable :: vu_stored_rows(:,:)
! pattern: u1,v1,u2,v2,u3,v3,....
  real(kind=rk8), allocatable :: uv_stored_cols(:,:)

  integer(kind=ik) :: istat
  character(200) :: errorMessage

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)
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
  tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

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
! todo: probably one should read it as v_row = vector v distributed among rows
!
  allocate(tmp(max(max_local_rows,max_local_cols)), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_real", "tmp", istat, errorMessage)

! allocate v_row 1 element longer to allow store and broadcast tau together with it
  allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_real", "v_row", istat, errorMessage)

  allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_real", "u_row", istat, errorMessage)

  allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_real", "v_col", istat, errorMessage)

  allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_real", "u_col", istat, errorMessage)

  tmp = 0
  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_real", "vu_stored_rows", istat, errorMessage)

  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_real", "uv_stored_cols", istat, errorMessage)

  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat
  if (my_prow==prow(na, nblk, np_rows) .and. my_pcol==pcol(na, nblk, np_cols)) &
    d_vec(na) = a_mat(l_rows,l_cols)

! main cycle of tridiagonalization
! in each step, 1 Householder vector is calculated
  do istep=na,3,-1

! Calculate number of local rows and columns of the still remaining matrix
! on the local processor
    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

! Calculate vector for Householder transformation on all procs
! owning column istep

    if (my_pcol==pcol(istep, nblk, np_cols)) then

! Get vector to be transformed; distribute last element and norm of
! remaining elements to all procs in current column

! copy l_cols + 1 column of A to v_row
      v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)

      if(n_stored_vecs>0 .and. l_rows>0) then
        call DGEMV('N', l_rows, 2*n_stored_vecs, &
             1.0_rk8, vu_stored_rows, ubound(vu_stored_rows,dim=1), &
             uv_stored_cols(l_cols+1,1), ubound(uv_stored_cols,dim=1), &
             1.0_rk8, v_row, 1)
      endif

      if(my_prow==prow(istep-1, nblk, np_rows)) then
        aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1))
        aux1(2) = v_row(l_rows)
      else
        aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
        aux1(2) = 0.
      endif

      call MPI_Allreduce(aux1, aux2, 2, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)

      vnorm2 = aux2(1)
      vrl = aux2(2)

! Householder transformation
      call hh_transform_real_double(vrl, vnorm2, xf, tau(istep))
! Scale v_row and store Householder vector for back transformation

      v_row(1:l_rows) = v_row(1:l_rows) * xf
      if (my_prow==prow(istep-1, nblk, np_rows)) then
        v_row(l_rows) = 1.

! vrl is newly computed off-diagonal element of the final tridiagonal matrix
        e_vec(istep-1) = vrl
      endif

! store Householder vector for back transformation
      a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)

! add tau after the end of actuall v_row, to be broadcasted with it
      v_row(l_rows+1) = tau(istep)
    endif !(my_pcol==pcol(istep, nblk, np_cols))

! Broadcast the Householder vector (and tau) along columns
    call MPI_Bcast(v_row, l_rows+1, MPI_REAL8, pcol(istep, nblk, np_cols), mpi_comm_cols, mpierr)

!recover tau, which has been broadcasted together with v_row
    tau(istep) = v_row(l_rows+1)

! Transpose Householder vector v_row -> v_col
    call elpa_transpose_vectors_real_double(v_row, ubound(v_row,dim=1), mpi_comm_rows, &
         v_col, ubound(v_col,dim=1), mpi_comm_cols, &
         1, istep-1, 1, nblk)

! Calculate u = (A + VU**T + UV**T)*v

! For cache efficiency, we use only the upper half of the matrix tiles for this,
! thus the result is partly in u_col(:) and partly in u_row(:)

    u_col(1:l_cols) = 0
    u_row(1:l_rows) = 0
    if (l_rows>0 .and. l_cols>0) then

      do i=0,(istep-2)/tile_size
        l_col_beg = i*l_cols_per_tile+1
        l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
        if (l_col_end<l_col_beg) cycle
        do j=0,i
          l_row_beg = j*l_rows_per_tile+1
          l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
          if (l_row_end<l_row_beg) cycle

          call DGEMV('T', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
               1.0_rk8, a_mat(l_row_beg, l_col_beg), lda, &
               v_row(l_row_beg), 1, &
               1.0_rk8, u_col(l_col_beg), 1)

          if (i/=j) then
            call DGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                 1.0_rk8, a_mat(l_row_beg,l_col_beg), lda, &
                 v_col(l_col_beg), 1, &
                 1.0_rk8, u_row(l_row_beg), 1)
          endif

        enddo ! j=0,i
      enddo ! i=0,(istep-2)/tile_size

      if (n_stored_vecs>0) then
        call DGEMV('T', l_rows, 2*n_stored_vecs, &
             1.0_rk8, vu_stored_rows, ubound(vu_stored_rows,dim=1), &
             v_row, 1, 0.0_rk8, aux, 1)
        call DGEMV('N', l_cols, 2*n_stored_vecs, &
             1.0_rk8, uv_stored_cols, ubound(uv_stored_cols,dim=1), &
             aux, 1, 1.0_rk8, u_col, 1)
      endif

    endif ! (l_rows>0 .and. l_cols>0)

! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
! on the processors containing the diagonal
! This is only necessary if u_row has been calculated, i.e. if the
! global tile size is smaller than the global remaining matrix

    if (tile_size < istep-1) then
      call elpa_reduce_add_vectors_real_double (u_row, ubound(u_row,dim=1), mpi_comm_rows, &
           u_col, ubound(u_col,dim=1), mpi_comm_cols, istep-1, 1, nblk)
    endif

! Sum up all the u_col(:) parts, transpose u_col -> u_row

    if (l_cols>0) then
      tmp(1:l_cols) = u_col(1:l_cols)

      call MPI_Allreduce(tmp, u_col, l_cols, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)

    endif

    call elpa_transpose_vectors_real_double (u_col, ubound(u_col,dim=1), mpi_comm_cols, &
         u_row, ubound(u_row,dim=1), mpi_comm_rows, 1, istep-1, 1, nblk)
! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )

    x = 0
    if (l_cols>0) &
      x = dot_product(v_col(1:l_cols),u_col(1:l_cols))

    call MPI_Allreduce(x, vav, 1, MPI_REAL8, MPI_SUM, mpi_comm_cols, mpierr)

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

! We have calculated another Hauseholder vector, number of implicitly stored increased
    n_stored_vecs = n_stored_vecs+1

! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
    if (n_stored_vecs==max_stored_uv .or. istep==3) then

      do i=0,(istep-2)/tile_size
        l_col_beg = i*l_cols_per_tile+1
        l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
        l_row_beg = 1
        l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
        if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) &
          cycle

        call DGEMM('N', 'T', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs, &
             1.0_rk8, vu_stored_rows(l_row_beg,1), ubound(vu_stored_rows,dim=1), &
             uv_stored_cols(l_col_beg,1), ubound(uv_stored_cols,dim=1), &
             1.0_rk8, a_mat(l_row_beg,l_col_beg), lda)
      enddo

      n_stored_vecs = 0

    endif

    if (my_prow==prow(istep-1, nblk, np_rows) .and. my_pcol==pcol(istep-1, nblk, np_cols)) then
      if (n_stored_vecs>0) then
        a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) &
          + dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs),uv_stored_cols(l_cols,1:2*n_stored_vecs))
      end if
      d_vec(istep-1) = a_mat(l_rows,l_cols)

    endif

  enddo ! main cycle over istep=na,3,-1

! Store e_vec(1)
  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(2, nblk, np_cols)) then
    e_vec(1) = a_mat(1,l_cols) ! use last l_cols value of loop above
  endif

! Store d_vec(1)
  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols)) then
    d_vec(1) = a_mat(1,1)
  endif

  deallocate(tmp, v_row, u_row, v_col, u_col, vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "tridiag_real: error when deallocating uv_stored_cols "//errorMessage
    stop
  endif

! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp(na), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "tridiag_real: error when allocating tmp "//errorMessage
    stop
  endif

  tmp = d_vec
  call MPI_Allreduce(tmp, d_vec, na, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)
  tmp = d_vec
  call MPI_Allreduce(tmp, d_vec, na, MPI_REAL8, MPI_SUM, mpi_comm_cols, mpierr)
  tmp = e_vec
  call MPI_Allreduce(tmp, e_vec, na, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)
  tmp = e_vec
  call MPI_Allreduce(tmp, e_vec, na, MPI_REAL8, MPI_SUM, mpi_comm_cols, mpierr)

  deallocate(tmp, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "tridiag_real: error when deallocating tmp "//errorMessage
    stop
  endif

end subroutine tridiag_real_double

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
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
subroutine trans_ev_real_double(na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

  use, intrinsic :: iso_c_binding
  use precision

  implicit none

  integer(kind=ik), intent(in) :: na, nqc, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols
  real(kind=rk8), intent(in) :: tau(na)

  real(kind=rk8), intent(inout) :: a_mat(lda,*), q_mat(ldq,*)

  integer(kind=ik) :: max_stored_rows

  integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
  integer(kind=ik) :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik) :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik) :: istep, n, nc, ic, ics, ice, nb, cur_pcol

  real(kind=rk8), allocatable :: tmp1(:), tmp2(:), hvb(:), hvm(:,:)
  real(kind=rk8), allocatable :: tmat(:,:), h1(:), h2(:)
  integer(kind=ik) :: istat
  character(200) :: errorMessage

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)
  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1 ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (63/nblk+1)*nblk

  allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_real", "tmat", istat, errorMessage)

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_real", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_real", "h2", istat, errorMessage)

  allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_real", "tmp1", istat, errorMessage)

  allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_real", "tmp2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_real", "hvn", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_real", "hvm", istat, errorMessage)

  hvm = 0 ! Must be set to 0 !!!
  hvb = 0 ! Safety only

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0

  do istep=1,na,nblk
    ics = max(istep,3)
    ice = min(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic=ics,ice

      l_colh = local_index(ic , my_pcol, np_cols, nblk, -1) ! Column of Householder vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector

      if (my_pcol==cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow==prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo

    if (nb>0) &
      call MPI_Bcast(hvb, nb, MPI_REAL8, cur_pcol, mpi_comm_cols, mpierr)

    nb = 0
    do ic=ics,ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      nstor = nstor+1
      nb = nb+l_rows
    enddo

! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32)) then

! Calculate scalar products of stored vectors.
! This can be done in different ways, we use dsyrk

      tmat = 0
      if (l_rows>0) &
        call DSYRK('U', 'T', nstor, l_rows, &
             1.0_rk8, hvm, ubound(hvm,dim=1), &
             0.0_rk8, tmat, max_stored_rows)

      nc = 0
      do n=1,nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo

      if (nc>0) call MPI_Allreduce( h1, h2, nc, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)

! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n=1,nstor-1
        call DTRMV('L', 'T', 'N', n, &
             tmat, max_stored_rows, &
             h2(nc+1), 1)
        tmat(n+1,1:n) = -h2(nc+1:nc+n)*tau(ice-nstor+n+1)
        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

! Q = Q - V * T * V**T * Q

      if (l_rows>0) then
        call DGEMM('T', 'N', nstor, l_cols, l_rows, &
             1.0_rk8, hvm, ubound(hvm,dim=1), &
             q_mat, ldq, &
             0.0_rk8, tmp1, nstor)

      else !l_rows>0
        tmp1(1:l_cols*nstor) = 0
      endif !l_rows>0

      call MPI_Allreduce(tmp1, tmp2, nstor*l_cols, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)

! copy back tmp2 - after reduction...

      if (l_rows>0) then
! tmp2 = tmat * tmp2
        call DTRMM('L', 'L', 'N', 'N', nstor, l_cols, &
             1.0_rk8, tmat, max_stored_rows, &
             tmp2, nstor)
!q_mat = q_mat - hvm*tmp2
        call DGEMM('N', 'N', l_rows, l_cols, nstor, &
             -1.0_rk8, hvm, ubound(hvm,dim=1), &
             tmp2, nstor, &
             1.0_rk8, q_mat, ldq)

      endif ! l_rows>0
      nstor = 0
    endif ! (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32))

  enddo

  deallocate(tmat, h1, h2, tmp1, tmp2, hvb, hvm, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "trans_ev_real: error when deallocating hvm "//errorMessage
    stop
  endif

end subroutine trans_ev_real_double

subroutine solve_tridi_double( na, nev, d, e, q, ldq, nblk, mpi_comm_rows, &
  mpi_comm_cols, success )

  use precision

  implicit none

  integer(kind=ik), intent(in) :: na, nev, ldq, nblk, mpi_comm_rows, mpi_comm_cols
  real(kind=rk8), intent(inout) :: d(na), e(na)

  real(kind=rk8), intent(inout) :: q(ldq,*)

  logical, intent(out) :: success

  integer(kind=ik) :: i, j, n, np, nc, nev1, l_cols, l_rows
  integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr

  integer(kind=ik), allocatable :: limits(:), l_col(:), p_col(:), l_col_bc(:), p_col_bc(:)

  integer(kind=ik) :: istat
  character(200) :: errorMessage

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

  success = .true.

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and q
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of q

! Set Q to 0
  q(1:l_rows, 1:l_cols) = 0.0_rk8

! Get the limits of the subdivisons, each subdivison has as many cols
! as fit on the respective processor column

  allocate(limits(0:np_cols), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi: error when allocating limits "//errorMessage
    stop
  endif

  limits(0) = 0
  do np=0,np_cols-1
    nc = local_index(na, np, np_cols, nblk, -1) ! number of columns on proc column np

! Check for the case that a column has have zero width.
! This is not supported!
! Scalapack supports it but delivers no results for these columns,
! which is rather annoying
    if (nc==0) then
      write(error_unit,*) 'ELPA1_solve_tridi: ERROR: Problem contains processor column with zero width'
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
    nev1 = min(nev,l_cols)
  endif
  call solve_tridi_col_double(l_cols, nev1, nc, d(nc+1), e(nc+1), q, ldq, nblk, &
       mpi_comm_rows, success)
  if (.not.(success)) then
    return
  endif
! If there is only 1 processor column, we are done

  if (np_cols==1) then
    deallocate(limits, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "solve_tridi: error when deallocating limits "//errorMessage
      stop
    endif

    return
  endif

! Set index arrays for Q columns

! Dense distribution scheme:

  allocate(l_col(na), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi: error when allocating l_col "//errorMessage
    stop
  endif

  allocate(p_col(na), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi: error when allocating p_col "//errorMessage
    stop
  endif

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
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi: error when allocating l_col_bc "//errorMessage
    stop
  endif

  allocate(p_col_bc(na), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi: error when allocating p_col_bc "//errorMessage
    stop
  endif

  p_col_bc(:) = -1
  l_col_bc(:) = -1

  do i = 0, na-1, nblk*np_cols
    do j = 0, np_cols-1
      do n = 1, nblk
        if (i+j*nblk+n <= min(nev,na)) then
          p_col_bc(i+j*nblk+n) = j
          l_col_bc(i+j*nblk+n) = i/np_cols + n
        endif
      enddo
    enddo
  enddo

! Recursively merge sub problems
  call merge_recursive_double(0, np_cols, success)
  if (.not.(success)) then
    return
  endif

  deallocate(limits,l_col,p_col,l_col_bc,p_col_bc, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi: error when deallocating l_col "//errorMessage
    stop
  endif

  return

  contains
  recursive subroutine merge_recursive_double(np_off, nprocs, success)

    use precision

    implicit none

! noff is always a multiple of nblk_ev
! nlen-noff is always > nblk_ev

    integer(kind=ik) :: np_off, nprocs
    integer(kind=ik) :: np1, np2, noff, nlen, nmid, n

    logical, intent(out) :: success

    success = .true.

    if (nprocs<=1) then
! Safety check only
      write(error_unit,*) "ELPA1_merge_recursive: INTERNAL error merge_recursive: nprocs=",nprocs
      success = .false.
      return
    endif
! Split problem into 2 subproblems of size np1 / np2

    np1 = nprocs/2
    np2 = nprocs-np1

    if (np1 > 1) call merge_recursive_double(np_off, np1, success)
    if (.not.(success)) return
    if (np2 > 1) call merge_recursive_double(np_off+np1, np2, success)
    if (.not.(success)) return

    noff = limits(np_off)
    nmid = limits(np_off+np1) - noff
    nlen = limits(np_off+nprocs) - noff

    if (my_pcol==np_off) then
      do n=np_off+np1,np_off+nprocs-1
        call MPI_Send(d(noff+1), nmid, MPI_REAL8, n, 1, mpi_comm_cols, mpierr)
      enddo
    endif

    if (my_pcol>=np_off+np1 .and. my_pcol<np_off+nprocs) then

      call MPI_Recv(d(noff+1), nmid, MPI_REAL8, np_off, 1, mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

    endif

    if (my_pcol==np_off+np1) then
      do n=np_off,np_off+np1-1

        call MPI_Send(d(noff+nmid+1), nlen-nmid, MPI_REAL8, n, 1, mpi_comm_cols, mpierr)

      enddo
    endif
    if (my_pcol>=np_off .and. my_pcol<np_off+np1) then

      call MPI_Recv(d(noff+nmid+1), nlen-nmid, MPI_REAL8, np_off+np1, 1,mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

    endif
    if (nprocs == np_cols) then

! Last merge, result distribution must be block cyclic, noff==0,
! p_col_bc is set so that only nev eigenvalues are calculated
      call merge_systems_double(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, noff, &
           nblk, mpi_comm_rows, mpi_comm_cols, l_col, p_col, &
           l_col_bc, p_col_bc, np_off, nprocs, success )
      if (.not.(success)) return
    else
! Not last merge, leave dense column distribution
      call merge_systems_double(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, noff, &
           nblk, mpi_comm_rows, mpi_comm_cols, l_col(noff+1), p_col(noff+1), &
           l_col(noff+1), p_col(noff+1), np_off, nprocs, success )
      if (.not.(success)) return
    endif

  end subroutine merge_recursive_double

end subroutine solve_tridi_double

subroutine solve_tridi_col_double( na, nev, nqoff, d, e, q, ldq, nblk, mpi_comm_rows, success )

! Solves the symmetric, tridiagonal eigenvalue problem on one processor column
! with the divide and conquer method.
! Works best if the number of processor rows is a power of 2!

  use precision

  implicit none

  integer(kind=ik) :: na, nev, nqoff, ldq, nblk, mpi_comm_rows
  real(kind=rk8) :: d(na), e(na)

  real(kind=rk8) :: q(ldq,*)

  integer(kind=ik), parameter :: min_submatrix_size = 16 ! Minimum size of the submatrices to be used

  real(kind=rk8), allocatable :: qmat1(:,:), qmat2(:,:)
  integer(kind=ik) :: i, n, np
  integer(kind=ik) :: ndiv, noff, nmid, nlen, max_size
  integer(kind=ik) :: my_prow, np_rows, mpierr

  integer(kind=ik), allocatable :: limits(:), l_col(:), p_col_i(:), p_col_o(:)
  logical, intent(out) :: success
  integer(kind=ik) :: istat
  character(200) :: errorMessage

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
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
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi_col: error when allocating limits "//errorMessage
    stop
  endif

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
    max_size = max(max_size,limits(i)-limits(i-1))
  enddo

! Subdivide matrix by subtracting rank 1 modifications

  do i=1,ndiv-1
    n = limits(i)
    d(n) = d(n)-abs(e(n))
    d(n+1) = d(n+1)-abs(e(n))
  enddo

  if (np_rows==1) then

! For 1 processor row there may be 1 or 2 subdivisions
    do n=0,ndiv-1
      noff = limits(n) ! Start of subproblem
      nlen = limits(n+1)-noff ! Size of subproblem

      call solve_tridi_single_problem_double(nlen,d(noff+1),e(noff+1), &
           q(nqoff+noff+1,noff+1),ubound(q,dim=1), success)

      if (.not.(success)) return
    enddo

  else

! Solve sub problems in parallel with solve_tridi_single
! There is at maximum 1 subproblem per processor

    allocate(qmat1(max_size,max_size), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "solve_tridi_col: error when allocating qmat1 "//errorMessage
      stop
    endif

    allocate(qmat2(max_size,max_size), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "solve_tridi_col: error when allocating qmat2 "//errorMessage
      stop
    endif

    qmat1 = 0 ! Make sure that all elements are defined

    if (my_prow < ndiv) then

      noff = limits(my_prow) ! Start of subproblem
      nlen = limits(my_prow+1)-noff ! Size of subproblem
      call solve_tridi_single_problem_double(nlen,d(noff+1),e(noff+1),qmat1, &
           ubound(qmat1,dim=1), success)

      if (.not.(success)) return
    endif

! Fill eigenvectors in qmat1 into global matrix q

    do np = 0, ndiv-1

      noff = limits(np)
      nlen = limits(np+1)-noff

      call MPI_Bcast(d(noff+1), nlen, MPI_REAL8, np, mpi_comm_rows, mpierr)
      qmat2 = qmat1
      call MPI_Bcast(qmat2, max_size*max_size, MPI_REAL8, np, mpi_comm_rows, mpierr)

      do i=1,nlen

        call distribute_global_column_double(qmat2(1,i), q(1,noff+i), nqoff+noff, nlen, my_prow, np_rows, nblk)

      enddo

    enddo

    deallocate(qmat1, qmat2, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "solve_tridi_col: error when deallocating qmat2 "//errorMessage
      stop
    endif

  endif

! Allocate and set index arrays l_col and p_col

  allocate(l_col(na), p_col_i(na), p_col_o(na), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi_col: error when allocating l_col "//errorMessage
    stop
  endif

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
      call merge_systems_double(nlen, nmid, d(noff+1), e(noff+nmid), q, ldq, nqoff+noff, nblk, &
           mpi_comm_rows, mpi_comm_self, l_col(noff+1), p_col_i(noff+1), &
           l_col(noff+1), p_col_o(noff+1), 0, 1, success)
      if (.not.(success)) return

    enddo

    n = 2*n

  enddo

  deallocate(limits, l_col, p_col_i, p_col_o, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi_col: error when deallocating l_col "//errorMessage
    stop
  endif

end subroutine solve_tridi_col_double

recursive subroutine solve_tridi_single_problem_double(nlen, d, e, q, ldq, success)

! Solves the symmetric, tridiagonal eigenvalue problem on a single processor.
! Takes precautions if DSTEDC fails or if the eigenvalues are not ordered correctly.

  use precision

  implicit none

  integer(kind=ik) :: nlen, ldq
  real(kind=rk8) :: d(nlen), e(nlen), q(ldq,nlen)

  real(kind=rk8), allocatable :: work(:), qtmp(:), ds(:), es(:)
  real(kind=rk8) :: dtmp

  integer(kind=ik) :: i, j, lwork, liwork, info
  integer(kind=ik), allocatable :: iwork(:)

  logical, intent(out) :: success
  integer(kind=ik) :: istat
  character(200) :: errorMessage

  success = .true.
  allocate(ds(nlen), es(nlen), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi_single: error when allocating ds "//errorMessage
    stop
  endif

! Save d and e for the case that dstedc fails

  ds(:) = d(:)
  es(:) = e(:)

! First try dstedc, this is normally faster but it may fail sometimes (why???)

  lwork = 1 + 4*nlen + nlen**2
  liwork = 3 + 5*nlen
  allocate(work(lwork), iwork(liwork), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi_single: error when allocating work "//errorMessage
    stop
  endif

  call DSTEDC('I', nlen, d, e, q, ldq, work, lwork, iwork, liwork, info)

  if (info /= 0) then

! DSTEDC failed, try DSTEQR. The workspace is enough for DSTEQR.

    write(error_unit,'(a,i8,a)') 'Warning: Lapack routine DSTEDC failed, info= ',info,', Trying DSTEQR!'

    d(:) = ds(:)
    e(:) = es(:)
    call DSTEQR('I', nlen, d, e, q, ldq, work, info)
! If DSTEQR fails also, we don't know what to do further ...

    if (info /= 0) then
      write(error_unit,'(a,i8,a)') 'ELPA1_solve_tridi_single: ERROR: Lapack routine DSTEQR failed, info= ',info,', Aborting!'
      success = .false.
      return
    endif
  end if

  deallocate(work,iwork,ds,es, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "solve_tridi_single: error when deallocating ds "//errorMessage
    stop
  endif

! Check if eigenvalues are monotonically increasing
! This seems to be not always the case (in the IBM implementation of dstedc ???)

  do i=1,nlen-1
    if (d(i+1)<d(i)) then

      if (abs(d(i+1) - d(i)) / abs(d(i+1) + d(i)) > 1e-14_rk8) then

        write(error_unit,'(a,i8,2g25.16)') '***WARNING: Monotony error dste**:',i+1,d(i),d(i+1)
      end if
      allocate(qtmp(nlen), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        write(error_unit,*) "solve_tridi_single: error when allocating qtmp "//errorMessage
        stop
      endif

      dtmp = d(i+1)
      qtmp(1:nlen) = q(1:nlen,i+1)
      do j=i,1,-1
        if (dtmp<d(j)) then
          d(j+1) = d(j)
          q(1:nlen,j+1) = q(1:nlen,j)
        else
          exit ! Loop
        endif
      enddo
      d(j+1) = dtmp
      q(1:nlen,j+1) = qtmp(1:nlen)
      deallocate(qtmp, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        write(error_unit,*) "solve_tridi_single: error when deallocating qtmp "//errorMessage
        stop
      endif
    endif
  enddo

end subroutine solve_tridi_single_problem_double

subroutine merge_systems_double( na, nm, d, e, q, ldq, nqoff, nblk, mpi_comm_rows, mpi_comm_cols, &
  l_col, p_col, l_col_out, p_col_out, npc_0, npc_n, success)

  use precision

  implicit none

  integer(kind=ik), intent(in) :: na, nm, ldq, nqoff, nblk, mpi_comm_rows, &
                                  mpi_comm_cols, npc_0, npc_n
  integer(kind=ik), intent(in) :: l_col(na), p_col(na), l_col_out(na), p_col_out(na)
  real(kind=rk8), intent(inout) :: d(na), e

  real(kind=rk8), intent(inout) :: q(ldq,*)

  logical, intent(out) :: success

  integer(kind=ik), parameter :: max_strip=128

  real(kind=rk8) :: DLAMCH, DLAPY2
  real(kind=rk8) :: beta, sig, s, c, t, tau, rho, eps, tol, &
                    qtrans(2,2), dmax, zmax, d1new, d2new
  real(kind=rk8) :: z(na), d1(na), d2(na), z1(na), delta(na), &
                    dbase(na), ddiff(na), ev_scale(na), tmp(na)
  real(kind=rk8) :: d1u(na), zu(na), d1l(na), zl(na)
  real(kind=rk8), allocatable :: qtmp1(:,:), qtmp2(:,:), ev(:,:)

  integer(kind=ik) :: i, j, na1, na2, l_rows, l_cols, l_rqs, l_rqe, &
                      l_rqm, ns, info
  integer(kind=ik) :: l_rnm, nnzu, nnzl, ndef, ncnt, max_local_cols, &
                      l_cols_qreorg, np, l_idx, nqcols1, nqcols2
  integer(kind=ik) :: my_proc, n_procs, my_prow, my_pcol, np_rows, &
                      np_cols, mpierr

  integer(kind=ik) :: np_next, np_prev, np_rem
  integer(kind=ik) :: idx(na), idx1(na), idx2(na)
  integer(kind=ik) :: coltyp(na), idxq1(na), idxq2(na)

  integer(kind=ik) :: istat
  character(200) :: errorMessage

  success = .true.
  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

! If my processor column isn't in the requested set, do nothing

  if (my_pcol<npc_0 .or. my_pcol>=npc_0+npc_n) then
    return
  endif
! Determine number of "next" and "prev" column for ring sends

  if (my_pcol == npc_0+npc_n-1) then
    np_next = npc_0
  else
    np_next = my_pcol + 1
  endif

  if (my_pcol == npc_0) then
    np_prev = npc_0+npc_n-1
  else
    np_prev = my_pcol - 1
  endif
  call check_monotony_double(nm,d,'Input1',success)
  if (.not.(success)) then
    return
  endif
  call check_monotony_double(na-nm,d(nm+1),'Input2',success)
  if (.not.(success)) then
    return
  endif
! Get global number of processors and my processor number.
! Please note that my_proc does not need to match any real processor number,
! it is just used for load balancing some loops.

  n_procs = np_rows*npc_n
  my_proc = my_prow*npc_n + (my_pcol-npc_0) ! Row major

! Local limits of the rows of Q

  l_rqs = local_index(nqoff+1 , my_prow, np_rows, nblk, +1) ! First row of Q
  l_rqm = local_index(nqoff+nm, my_prow, np_rows, nblk, -1) ! Last row <= nm
  l_rqe = local_index(nqoff+na, my_prow, np_rows, nblk, -1) ! Last row of Q

  l_rnm = l_rqm-l_rqs+1 ! Number of local rows <= nm
  l_rows = l_rqe-l_rqs+1 ! Total number of local rows

! My number of local columns

  l_cols = count(p_col(1:na)==my_pcol)

! Get max number of local columns

  max_local_cols = 0
  do np = npc_0, npc_0+npc_n-1
    max_local_cols = max(max_local_cols,count(p_col(1:na)==np))
  enddo

! Calculations start here

  beta = abs(e)
  sig = sign(1.0_rk8,e)

! Calculate rank-1 modifier z

  z(:) = 0

  if (mod((nqoff+nm-1)/nblk,np_rows)==my_prow) then
! nm is local on my row
    do i = 1, na
      if (p_col(i)==my_pcol) z(i) = q(l_rqm,l_col(i))
    enddo
  endif

  if (mod((nqoff+nm)/nblk,np_rows)==my_prow) then
! nm+1 is local on my row
    do i = 1, na
      if (p_col(i)==my_pcol) z(i) = z(i) + sig*q(l_rqm+1,l_col(i))
    enddo
  endif

  call global_gather_double(z, na)
! Normalize z so that norm(z) = 1. Since z is the concatenation of
! two normalized vectors, norm2(z) = sqrt(2).
  z = z/sqrt(2.0_rk8)
  rho = 2.0_rk8*beta
! Calculate index for merging both systems by ascending eigenvalues
  call DLAMRG( nm, na-nm, d, 1, 1, idx )

! Calculate the allowable deflation tolerance

  zmax = maxval(abs(z))
  dmax = maxval(abs(d))
  EPS = DLAMCH( 'Epsilon' )
  TOL = 8.0_rk8*EPS*max(dmax,zmax)

! If the rank-1 modifier is small enough, no more needs to be done
! except to reorganize D and Q

  if ( RHO*zmax <= TOL ) then

! Rearrange eigenvalues

    tmp = d
    do i=1,na
      d(i) = tmp(idx(i))
    enddo

! Rearrange eigenvectors
    call resort_ev_double(idx, na)

    return
  endif

! Merge and deflate system

  na1 = 0
  na2 = 0

! COLTYP:
! 1 : non-zero in the upper half only;
! 2 : dense;
! 3 : non-zero in the lower half only;
! 4 : deflated.

  coltyp(1:nm) = 1
  coltyp(nm+1:na) = 3

  do i=1,na

    if (rho*abs(z(idx(i))) <= tol) then

! Deflate due to small z component.

      na2 = na2+1
      d2(na2) = d(idx(i))
      idx2(na2) = idx(i)
      coltyp(idx(i)) = 4

    else if (na1>0) then

! Check if eigenvalues are close enough to allow deflation.

      S = Z(idx(i))
      C = Z1(na1)

! Find sqrt(a**2+b**2) without overflow or
! destructive underflow.
      TAU = DLAPY2( C, S )
      T = D1(na1) - D(idx(i))
      C = C / TAU
      S = -S / TAU
      if ( abs( T*C*S ) <= TOL ) then

! Deflation is possible.

        na2 = na2+1

        Z1(na1) = TAU

        d2new = D(idx(i))*C**2 + D1(na1)*S**2
        d1new = D(idx(i))*S**2 + D1(na1)*C**2

! D(idx(i)) >= D1(na1) and C**2 + S**2 == 1.0
! This means that after the above transformation it must be
!    D1(na1) <= d1new <= D(idx(i))
!    D1(na1) <= d2new <= D(idx(i))
!
! D1(na1) may get bigger but it is still smaller than the next D(idx(i+1))
! so there is no problem with sorting here.
! d2new <= D(idx(i)) which means that it might be smaller than D2(na2-1)
! which makes a check (and possibly a resort) necessary.
!
! The above relations may not hold exactly due to numeric differences
! so they have to be enforced in order not to get troubles with sorting.

        if (d1new<D1(na1)) d1new = D1(na1)
        if (d1new>D(idx(i))) d1new = D(idx(i))

        if (d2new<D1(na1)) d2new = D1(na1)
        if (d2new>D(idx(i))) d2new = D(idx(i))

        D1(na1) = d1new

        do j=na2-1,1,-1
          if (d2new<d2(j)) then
            d2(j+1) = d2(j)
            idx2(j+1) = idx2(j)
          else
            exit ! Loop
          endif
        enddo

        d2(j+1) = d2new
        idx2(j+1) = idx(i)

        qtrans(1,1) = C; qtrans(1,2) =-S
        qtrans(2,1) = S; qtrans(2,2) = C
        call transform_columns_double(idx(i), idx1(na1))
        if (coltyp(idx(i))==1 .and. coltyp(idx1(na1))/=1) coltyp(idx1(na1)) = 2
        if (coltyp(idx(i))==3 .and. coltyp(idx1(na1))/=3) coltyp(idx1(na1)) = 2

        coltyp(idx(i)) = 4

      else
        na1 = na1+1
        d1(na1) = d(idx(i))
        z1(na1) = z(idx(i))
        idx1(na1) = idx(i)
      endif
    else
      na1 = na1+1
      d1(na1) = d(idx(i))
      z1(na1) = z(idx(i))
      idx1(na1) = idx(i)
    endif

  enddo
  call check_monotony_double(na1,d1,'Sorted1',success)
  if (.not.(success)) then
    return
  endif
  call check_monotony_double(na2,d2,'Sorted2',success)
  if (.not.(success)) then
    return
  endif

  if (na1==1 .or. na1==2) then
    if (na1==1) then
      d(1) = d1(1) + rho*z1(1)**2 ! solve secular equation
    else ! na1==2
      call DLAED5(1, d1, z1, qtrans(1,1), rho, d(1))
      call DLAED5(2, d1, z1, qtrans(1,2), rho, d(2))

      call transform_columns_double(idx1(1), idx1(2))
    endif

! Add the deflated eigenvalues
    d(na1+1:na) = d2(1:na2)

! Calculate arrangement of all eigenvalues in output
    call DLAMRG( na1, na-na1, d, 1, 1, idx )
! Rearrange eigenvalues

    tmp = d
    do i=1,na
      d(i) = tmp(idx(i))
    enddo

! Rearrange eigenvectors

    do i=1,na
      if (idx(i)<=na1) then
        idxq1(i) = idx1(idx(i))
      else
        idxq1(i) = idx2(idx(i)-na1)
      endif
    enddo
    call resort_ev_double(idxq1, na)
  else if (na1>2) then

! Solve secular equation

    z(1:na1) = 1

    dbase(1:na1) = 0
    ddiff(1:na1) = 0

    info = 0

    do i = my_proc+1, na1, n_procs ! work distributed over all processors
      call DLAED4(na1, i, d1, z1, delta, rho, s, info) ! s is not used!
      if (info/=0) then
! If DLAED4 fails (may happen especially for LAPACK versions before 3.2)
! use the more stable bisection algorithm in solve_secular_equation
        call solve_secular_equation_double(na1, i, d1, z1, delta, rho, s)
      endif

! Compute updated z

      do j=1,na1
        if (i/=j) z(j) = z(j)*( delta(j) / (d1(j)-d1(i)) )
      enddo
      z(i) = z(i)*delta(i)

! store dbase/ddiff

      if (i<na1) then
        if (abs(delta(i+1)) < abs(delta(i))) then
          dbase(i) = d1(i+1)
          ddiff(i) = delta(i+1)
        else
          dbase(i) = d1(i)
          ddiff(i) = delta(i)
        endif
      else
        dbase(i) = d1(i)
        ddiff(i) = delta(i)
      endif
    enddo

    call global_product_double(z, na1)
    z(1:na1) = SIGN( SQRT( -z(1:na1) ), z1(1:na1) )

    call global_gather_double(dbase, na1)
    call global_gather_double(ddiff, na1)
    d(1:na1) = dbase(1:na1) - ddiff(1:na1)

! Calculate scale factors for eigenvectors
    ev_scale(:) = 0.0_rk8

    do i = my_proc+1, na1, n_procs ! work distributed over all processors

! tmp(1:na1) = z(1:na1) / delta(1:na1,i) ! original code
! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
! in exactly this order, but we want to prevent compiler optimization
!         ev_scale_val = ev_scale(i)
      call add_tmp_double(d1, dbase, ddiff, z, ev_scale(i), na1,i)
!         ev_scale(i) = ev_scale_val
    enddo

    call global_gather_double(ev_scale, na1)
! Add the deflated eigenvalues
    d(na1+1:na) = d2(1:na2)

! Calculate arrangement of all eigenvalues in output
    call DLAMRG( na1, na-na1, d, 1, 1, idx )

! Rearrange eigenvalues
    tmp = d
    do i=1,na
      d(i) = tmp(idx(i))
    enddo
    call check_monotony_double(na,d,'Output',success)

    if (.not.(success)) then
      return
    endif
! Eigenvector calculations

! Calculate the number of columns in the new local matrix Q
! which are updated from non-deflated/deflated eigenvectors.
! idxq1/2 stores the global column numbers.

    nqcols1 = 0 ! number of non-deflated eigenvectors
    nqcols2 = 0 ! number of deflated eigenvectors
    do i = 1, na
      if (p_col_out(i)==my_pcol) then
        if (idx(i)<=na1) then
          nqcols1 = nqcols1+1
          idxq1(nqcols1) = i
        else
          nqcols2 = nqcols2+1
          idxq2(nqcols2) = i
        endif
      endif
    enddo

    allocate(ev(max_local_cols,min(max_strip,max(1,nqcols1))), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "merge_systems: error when allocating ev "//errorMessage
      stop
    endif

    allocate(qtmp1(max(1,l_rows),max_local_cols), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "merge_systems: error when allocating qtmp1 "//errorMessage
      stop
    endif

    allocate(qtmp2(max(1,l_rows),min(max_strip,max(1,nqcols1))), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "merge_systems: error when allocating qtmp2 "//errorMessage
      stop
    endif

! Gather nonzero upper/lower components of old matrix Q
! which are needed for multiplication with new eigenvectors

    qtmp1 = 0 ! May contain empty (unset) parts
    qtmp2 = 0 ! Not really needed

    nnzu = 0
    nnzl = 0
    do i = 1, na1
      l_idx = l_col(idx1(i))
      if (p_col(idx1(i))==my_pcol) then
        if (coltyp(idx1(i))==1 .or. coltyp(idx1(i))==2) then
          nnzu = nnzu+1
          qtmp1(1:l_rnm,nnzu) = q(l_rqs:l_rqm,l_idx)
        endif
        if (coltyp(idx1(i))==3 .or. coltyp(idx1(i))==2) then
          nnzl = nnzl+1
          qtmp1(l_rnm+1:l_rows,nnzl) = q(l_rqm+1:l_rqe,l_idx)
        endif
      endif
    enddo

! Gather deflated eigenvalues behind nonzero components

    ndef = max(nnzu,nnzl)
    do i = 1, na2
      l_idx = l_col(idx2(i))
      if (p_col(idx2(i))==my_pcol) then
        ndef = ndef+1
        qtmp1(1:l_rows,ndef) = q(l_rqs:l_rqe,l_idx)
      endif
    enddo

    l_cols_qreorg = ndef ! Number of columns in reorganized matrix

! Set (output) Q to 0, it will sum up new Q

    do i = 1, na
      if(p_col_out(i)==my_pcol) q(l_rqs:l_rqe,l_col_out(i)) = 0
    enddo

    np_rem = my_pcol

    do np = 1, npc_n

! Do a ring send of qtmp1

      if (np>1) then

        if (np_rem==npc_0) then
          np_rem = npc_0+npc_n-1
        else
          np_rem = np_rem-1
        endif

        call MPI_Sendrecv_replace(qtmp1, l_rows*max_local_cols, MPI_REAL8, &
             np_next, 1111, np_prev, 1111, &
             mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

      endif

! Gather the parts in d1 and z which are fitting to qtmp1.
! This also delivers nnzu/nnzl for proc np_rem

      nnzu = 0
      nnzl = 0
      do i=1,na1
        if (p_col(idx1(i))==np_rem) then
          if (coltyp(idx1(i))==1 .or. coltyp(idx1(i))==2) then
            nnzu = nnzu+1
            d1u(nnzu) = d1(i)
            zu (nnzu) = z (i)
          endif
          if (coltyp(idx1(i))==3 .or. coltyp(idx1(i))==2) then
            nnzl = nnzl+1
            d1l(nnzl) = d1(i)
            zl (nnzl) = z (i)
          endif
        endif
      enddo

! Set the deflated eigenvectors in Q (comming from proc np_rem)

      ndef = max(nnzu,nnzl) ! Remote counter in input matrix
      do i = 1, na
        j = idx(i)
        if (j>na1) then
          if (p_col(idx2(j-na1))==np_rem) then
            ndef = ndef+1
            if (p_col_out(i)==my_pcol) &
              q(l_rqs:l_rqe,l_col_out(i)) = qtmp1(1:l_rows,ndef)
          endif
        endif
      enddo

      do ns = 0, nqcols1-1, max_strip ! strimining loop

        ncnt = min(max_strip,nqcols1-ns) ! number of columns in this strip

! Get partial result from (output) Q

        do i = 1, ncnt
          qtmp2(1:l_rows,i) = q(l_rqs:l_rqe,l_col_out(idxq1(i+ns)))
        enddo

! Compute eigenvectors of the rank-1 modified matrix.
! Parts for multiplying with upper half of Q:

        do i = 1, ncnt
          j = idx(idxq1(i+ns))
! Calculate the j-th eigenvector of the deflated system
! See above why we are doing it this way!
          tmp(1:nnzu) = d1u(1:nnzu)-dbase(j)
          call v_add_s_double(tmp,nnzu,ddiff(j))
          ev(1:nnzu,i) = zu(1:nnzu) / tmp(1:nnzu) * ev_scale(j)
        enddo

! Multiply old Q with eigenvectors (upper half)

        if (l_rnm>0 .and. ncnt>0 .and. nnzu>0) &
          call DGEMM('N', 'N', l_rnm, ncnt, nnzu, 1.0_rk8, qtmp1, ubound(qtmp1,dim=1), ev, ubound(ev,dim=1), &
               1.0_rk8, qtmp2(1,1), ubound(qtmp2,dim=1))
! Compute eigenvectors of the rank-1 modified matrix.
! Parts for multiplying with lower half of Q:

        do i = 1, ncnt
          j = idx(idxq1(i+ns))
! Calculate the j-th eigenvector of the deflated system
! See above why we are doing it this way!
          tmp(1:nnzl) = d1l(1:nnzl)-dbase(j)
          call v_add_s_double(tmp,nnzl,ddiff(j))
          ev(1:nnzl,i) = zl(1:nnzl) / tmp(1:nnzl) * ev_scale(j)
        enddo

! Multiply old Q with eigenvectors (lower half)

        if (l_rows-l_rnm>0 .and. ncnt>0 .and. nnzl>0) &
          call DGEMM('N', 'N', l_rows-l_rnm, ncnt, nnzl, 1.0_rk8, qtmp1(l_rnm+1,1), ubound(qtmp1,dim=1), ev, &
               ubound(ev,dim=1), 1.0_rk8, qtmp2(l_rnm+1,1), ubound(qtmp2,dim=1))
! Put partial result into (output) Q

        do i = 1, ncnt
          q(l_rqs:l_rqe,l_col_out(idxq1(i+ns))) = qtmp2(1:l_rows,i)
        enddo

      enddo
    enddo

    deallocate(ev, qtmp1, qtmp2, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "merge_systems: error when deallocating ev "//errorMessage
      stop
    endif
  endif

  return

  contains
  subroutine add_tmp_double(d1, dbase, ddiff, z, ev_scale_value, na1,i)

    use precision

    implicit none

    integer(kind=ik), intent(in) :: na1, i

    real(kind=rk8), intent(in) :: d1(:), dbase(:), ddiff(:), z(:)
    real(kind=rk8), intent(inout) :: ev_scale_value
    real(kind=rk8) :: tmp(1:na1)

! tmp(1:na1) = z(1:na1) / delta(1:na1,i) ! original code
! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
! in exactly this order, but we want to prevent compiler optimization

    tmp(1:na1) = d1(1:na1) -dbase(i)
    call v_add_s_double(tmp(1:na1),na1,ddiff(i))
    tmp(1:na1) = z(1:na1) / tmp(1:na1)
    ev_scale_value = 1.0_rk8/sqrt(dot_product(tmp(1:na1),tmp(1:na1)))

  end subroutine add_tmp_double

  subroutine resort_ev_double(idx_ev, nLength)

    use precision

    implicit none

    integer(kind=ik), intent(in) :: nLength
    integer(kind=ik) :: idx_ev(nLength)
    integer(kind=ik) :: i, nc, pc1, pc2, lc1, lc2, l_cols_out

    real(kind=rk8), allocatable :: qtmp(:,:)
    integer(kind=ik) :: istat
    character(200) :: errorMessage

    if (l_rows==0) return ! My processor column has no work to do

! Resorts eigenvectors so that q_new(:,i) = q_old(:,idx_ev(i))

    l_cols_out = count(p_col_out(1:na)==my_pcol)
    allocate(qtmp(l_rows,l_cols_out), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "resort_ev: error when allocating qtmp "//errorMessage
      stop
    endif

    nc = 0

    do i=1,na

      pc1 = p_col(idx_ev(i))
      lc1 = l_col(idx_ev(i))
      pc2 = p_col_out(i)

      if (pc2<0) cycle ! This column is not needed in output

      if (pc2==my_pcol) nc = nc+1 ! Counter for output columns

      if (pc1==my_pcol) then
        if (pc2==my_pcol) then
! send and recieve column are local
          qtmp(1:l_rows,nc) = q(l_rqs:l_rqe,lc1)
        else

          call MPI_Send(q(l_rqs,lc1), l_rows, MPI_REAL8, pc2, mod(i,4096), mpi_comm_cols, mpierr)

        endif
      else if (pc2==my_pcol) then

        call MPI_Recv(qtmp(1,nc), l_rows, MPI_REAL8, pc1, mod(i,4096), mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

      endif
    enddo

! Insert qtmp into (output) q

    nc = 0

    do i=1,na

      pc2 = p_col_out(i)
      lc2 = l_col_out(i)

      if (pc2==my_pcol) then
        nc = nc+1
        q(l_rqs:l_rqe,lc2) = qtmp(1:l_rows,nc)
      endif
    enddo

    deallocate(qtmp, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "resort_ev: error when deallocating qtmp "//errorMessage
      stop
    endif

  end subroutine resort_ev_double

  subroutine transform_columns_double(col1, col2)

    use precision

    implicit none

    integer(kind=ik) :: col1, col2
    integer(kind=ik) :: pc1, pc2, lc1, lc2

    if (l_rows==0) return ! My processor column has no work to do

    pc1 = p_col(col1)
    lc1 = l_col(col1)
    pc2 = p_col(col2)
    lc2 = l_col(col2)

    if (pc1==my_pcol) then
      if (pc2==my_pcol) then
! both columns are local
        tmp(1:l_rows) = q(l_rqs:l_rqe,lc1)*qtrans(1,1) + q(l_rqs:l_rqe,lc2)*qtrans(2,1)
        q(l_rqs:l_rqe,lc2) = q(l_rqs:l_rqe,lc1)*qtrans(1,2) + q(l_rqs:l_rqe,lc2)*qtrans(2,2)
        q(l_rqs:l_rqe,lc1) = tmp(1:l_rows)
      else

        call MPI_Sendrecv(q(l_rqs,lc1), l_rows, MPI_REAL8, pc2, 1, &
             tmp, l_rows, MPI_REAL8, pc2, 1, &
             mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

        q(l_rqs:l_rqe,lc1) = q(l_rqs:l_rqe,lc1)*qtrans(1,1) + tmp(1:l_rows)*qtrans(2,1)
      endif
    else if (pc2==my_pcol) then

      call MPI_Sendrecv(q(l_rqs,lc2), l_rows, MPI_REAL8, pc1, 1, &
           tmp, l_rows, MPI_REAL8, pc1, 1, &
           mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

      q(l_rqs:l_rqe,lc2) = tmp(1:l_rows)*qtrans(1,2) + q(l_rqs:l_rqe,lc2)*qtrans(2,2)
    endif

  end subroutine transform_columns_double

  subroutine global_gather_double(z, n)
! This routine sums up z over all processors.
! It should only be used for gathering distributed results,
! i.e. z(i) should be nonzero on exactly 1 processor column,
! otherways the results may be numerically different on different columns
    use precision

    implicit none

    integer(kind=ik) :: n
    real(kind=rk8) :: z(n)
    real(kind=rk8) :: tmp(n)

    if (npc_n==1 .and. np_rows==1) return ! nothing to do

! Do an mpi_allreduce over processor rows

    call MPI_Allreduce(z, tmp, n, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)

! If only 1 processor column, we are done
    if (npc_n==1) then
      z(:) = tmp(:)
      return
    endif

! If all processor columns are involved, we can use mpi_allreduce
    if (npc_n==np_cols) then

      call MPI_Allreduce(tmp, z, n, MPI_REAL8, MPI_SUM, mpi_comm_cols, mpierr)

      return
    endif

! Do a ring send over processor columns
    z(:) = 0
    do np = 1, npc_n
      z(:) = z(:) + tmp(:)

      call MPI_Sendrecv_replace(z, n, MPI_REAL8, np_next, 1111, np_prev, 1111, &
           mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

    enddo

  end subroutine global_gather_double

  subroutine global_product_double(z, n)
! This routine calculates the global product of z.
    use precision

    implicit none

    integer(kind=ik) :: n
    real(kind=rk8) :: z(n)

    real(kind=rk8) :: tmp(n)

    if (npc_n==1 .and. np_rows==1) return ! nothing to do

! Do an mpi_allreduce over processor rows

    call MPI_Allreduce(z, tmp, n, MPI_REAL8, MPI_PROD, mpi_comm_rows, mpierr)

! If only 1 processor column, we are done
    if (npc_n==1) then
      z(:) = tmp(:)
      return
    endif

! If all processor columns are involved, we can use mpi_allreduce
    if (npc_n==np_cols) then

      call MPI_Allreduce(tmp, z, n, MPI_REAL8, MPI_PROD, mpi_comm_cols, mpierr)

      return
    endif

! We send all vectors to the first proc, do the product there
! and redistribute the result.

    if (my_pcol == npc_0) then
      z(1:n) = tmp(1:n)
      do np = npc_0+1, npc_0+npc_n-1

        call MPI_Recv(tmp, n, MPI_REAL8, np, 1111, mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

        z(1:n) = z(1:n)*tmp(1:n)
      enddo
      do np = npc_0+1, npc_0+npc_n-1

        call MPI_Send(z, n, MPI_REAL8, np, 1111, mpi_comm_cols, mpierr)

      enddo
    else

      call MPI_Send(tmp, n, MPI_REAL8, npc_0, 1111, mpi_comm_cols, mpierr)
      call MPI_Recv(z ,n, MPI_REAL8, npc_0, 1111, mpi_comm_cols, MPI_STATUS_IGNORE, mpierr)

    endif

  end subroutine global_product_double

  subroutine check_monotony_double(n,d,text,success)
! This is a test routine for checking if the eigenvalues are monotonically increasing.
! It is for debug purposes only, an error should never be triggered!
    use precision

    implicit none

    integer(kind=ik) :: n
    real(kind=rk8) :: d(n)
    character*(*) :: text

    integer(kind=ik) :: i
    logical, intent(out) :: success

    success = .true.
    do i=1,n-1
      if (d(i+1)<d(i)) then
        write(error_unit,'(a,a,i8,2g25.17)') 'ELPA1_check_monotony: Monotony error on ',text,i,d(i),d(i+1)
        success = .false.
        return
      endif
    enddo

  end subroutine check_monotony_double

end subroutine merge_systems_double

subroutine v_add_s_double(v,n,s)

  use precision

  implicit none

  integer(kind=ik) :: n
  real(kind=rk8) :: v(n),s

  v(:) = v(:) + s
end subroutine v_add_s_double

subroutine distribute_global_column_double(g_col, l_col, noff, nlen, my_prow, np_rows, nblk)

  use precision

  implicit none

  integer(kind=ik) :: noff, nlen, my_prow, np_rows, nblk

  real(kind=rk8) :: g_col(nlen), l_col(*) ! chnage this to proper 2d 1d matching ! remove assumed size
  integer(kind=ik) :: nbs, nbe, jb, g_off, l_off, js, je

  nbs = noff/(nblk*np_rows)
  nbe = (noff+nlen-1)/(nblk*np_rows)

  do jb = nbs, nbe

    g_off = jb*nblk*np_rows + nblk*my_prow
    l_off = jb*nblk

    js = max(noff+1-g_off,1)
    je = min(noff+nlen-g_off,nblk)

    if (je<js) cycle

    l_col(l_off+js:l_off+je) = g_col(g_off+js-noff:g_off+je-noff)

  enddo
end subroutine distribute_global_column_double

subroutine solve_secular_equation_double(n, i, d, z, delta, rho, dlam)
!-------------------------------------------------------------------------------
! This routine solves the secular equation of a symmetric rank 1 modified
! diagonal matrix:
!
!    1. + rho*SUM(z(:)**2/(d(:)-x)) = 0
!
! It does the same as the LAPACK routine DLAED4 but it uses a bisection technique
! which is more robust (it always yields a solution) but also slower
! than the algorithm used in DLAED4.
!
! The same restictions than in DLAED4 hold, namely:
!
!   rho > 0   and   d(i+1) > d(i)
!
! but this routine will not terminate with error if these are not satisfied
! (it will normally converge to a pole in this case).
!
! The output in DELTA(j) is always (D(j) - lambda_I), even for the cases
! N=1 and N=2 which is not compatible with DLAED4.
! Thus this routine shouldn't be used for these cases as a simple replacement
! of DLAED4.
!
! The arguments are the same as in DLAED4 (with the exception of the INFO argument):
!
!
!  N      (input) INTEGER
!         The length of all arrays.
!
!  I      (input) INTEGER
!         The index of the eigenvalue to be computed.  1 <= I <= N.
!
!  D      (input) DOUBLE PRECISION array, dimension (N)
!         The original eigenvalues.  It is assumed that they are in
!         order, D(I) < D(J)  for I < J.
!
!  Z      (input) DOUBLE PRECISION array, dimension (N)
!         The components of the updating vector.
!
!  DELTA  (output) DOUBLE PRECISION array, dimension (N)
!         DELTA contains (D(j) - lambda_I) in its  j-th component.
!         See remark above about DLAED4 compatibility!
!
!  RHO    (input) DOUBLE PRECISION
!         The scalar in the symmetric updating formula.
!
!  DLAM   (output) DOUBLE PRECISION
!         The computed lambda_I, the I-th updated eigenvalue.
!-------------------------------------------------------------------------------

  use precision

  implicit none

  integer(kind=ik) :: n, i
  real(kind=rk8) :: d(n), z(n), delta(n), rho, dlam

  integer(kind=ik) :: iter
  real(kind=rk8) :: a, b, x, y, dshift

! In order to obtain sufficient numerical accuracy we have to shift the problem
! either by d(i) or d(i+1), whichever is closer to the solution

! Upper and lower bound of the shifted solution interval are a and b

  if (i==n) then

! Special case: Last eigenvalue
! We shift always by d(n), lower bound is d(n),
! upper bound is determined by a guess:

    dshift = d(n)
    delta(:) = d(:) - dshift

    a = 0.0_rk8 ! delta(n)
    b = rho*SUM(z(:)**2) + 1.0_rk8 ! rho*SUM(z(:)**2) is the lower bound for the guess
  else

! Other eigenvalues: lower bound is d(i), upper bound is d(i+1)
! We check the sign of the function in the midpoint of the interval
! in order to determine if eigenvalue is more close to d(i) or d(i+1)
    x = 0.5_rk8*(d(i)+d(i+1))
    y = 1.0_rk8 + rho*SUM(z(:)**2/(d(:)-x))
    if (y>0) then
! solution is next to d(i)
      dshift = d(i)
    else
! solution is next to d(i+1)
      dshift = d(i+1)
    endif

    delta(:) = d(:) - dshift
    a = delta(i)
    b = delta(i+1)

  endif

! Bisection:

  do iter=1,200

! Interval subdivision
    x = 0.5_rk8*(a+b)
    if (x==a .or. x==b) exit ! No further interval subdivisions possible

    if (abs(x) < 1.e-200_rk8) exit ! x next to pole

! evaluate value at x

    y = 1. + rho*SUM(z(:)**2/(delta(:)-x))

    if (y==0) then
! found exact solution
      exit
    elseif (y>0) then
      b = x
    else
      a = x
    endif

  enddo

! Solution:

  dlam = x + dshift
  delta(:) = delta(:) - x

end subroutine solve_secular_equation_double
!-------------------------------------------------------------------------------

subroutine hh_transform_real_double(alpha, xnorm_sq, xf, tau)
! Similar to LAPACK routine DLARFP, but uses ||x||**2 instead of x(:)
! and returns the factor xf by which x has to be scaled.
! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
! since this would be expensive for the parallel implementation.

  use precision

  implicit none

  real(kind=rk8), intent(inout) :: alpha
  real(kind=rk8), intent(in) :: xnorm_sq
  real(kind=rk8), intent(out) :: xf, tau

  real(kind=rk8) :: BETA

  if ( XNORM_SQ==0. ) then

    if ( ALPHA>=0. ) then
      TAU = 0.
    else
      TAU = 2.
      ALPHA = -ALPHA
    endif
    XF = 0.

  else

    BETA = SIGN( SQRT( ALPHA**2 + XNORM_SQ ), ALPHA )
    ALPHA = ALPHA + BETA
    if ( BETA<0 ) then
      BETA = -BETA
      TAU = -ALPHA / BETA
    else
      ALPHA = XNORM_SQ / ALPHA
      TAU = ALPHA / BETA
      ALPHA = -ALPHA
    end if
      XF = 1./ALPHA
      ALPHA = BETA
  endif
end subroutine hh_transform_real_double

! complex double precision

!> \brief Reduces a distributed symmetric matrix to tridiagonal form (like Scalapack Routine PDSYTRD)
!>
!  Parameters
!
!> \param na          Order of matrix
!>
!> \param a_mat(lda,matrixCols)    Distributed matrix which should be reduced.
!>              Distribution is like in Scalapack.
!>              Opposed to PDSYTRD, a(:,:) must be set completely (upper and lower half)
!>              a(:,:) is overwritten on exit with the Householder vectors
!>
!> \param lda         Leading dimension of a
!>
!> \param nblk        blocksize of cyclic distribution, must be the same in both directions!
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
subroutine tridiag_complex_double(na, a_mat, lda, nblk, mpi_comm_rows, mpi_comm_cols, &
  d_vec, e_vec, tau)

  use, intrinsic :: iso_c_binding
  use precision

  implicit none

  integer(kind=ik), intent(in) :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols

  complex(kind=ck8), intent(out) :: tau(na)

  complex(kind=ck8), intent(inout) :: a_mat(lda,*)

  real(kind=rk8), intent(out) :: d_vec(na), e_vec(na)

  integer(kind=ik), parameter :: max_stored_uv = 32

  complex(kind=ck8), parameter :: CZERO = (0.0_rk8,0.0_rk8), CONE = (1.0_rk8,0.0_rk8)

  integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr

  integer(kind=ik) :: totalblocks, max_loc_block_rows, max_loc_block_cols, max_local_rows, max_local_cols
  integer(kind=ik) :: l_cols, l_rows, n_stored_vecs
  integer(kind=ik) :: istep, i, j, l_col_beg, l_col_end, l_row_beg, l_row_end
  integer(kind=ik) :: tile_size, l_rows_per_tile, l_cols_per_tile

  real(kind=rk8) :: vnorm2
  complex(kind=ck8) :: vav, xc, aux(2*max_stored_uv), aux1(2), aux2(2), aux3(1), vrl, xf

  complex(kind=ck8), allocatable :: tmp(:), v_row(:), v_col(:), u_row(:), u_col(:), vu_stored_rows(:,:), uv_stored_cols(:,:)

  real(kind=rk8), allocatable :: tmp_real(:)
  integer(kind=ik) :: istat
  character(200) :: errorMessage

  aux3 = 0

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

! Matrix is split into tiles; work is done only for tiles on the diagonal or above
  tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
  tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

  l_rows_per_tile = tile_size/np_rows ! local rows of a tile
  l_cols_per_tile = tile_size/np_cols ! local cols of a tile

  totalblocks = (na-1)/nblk + 1
  max_loc_block_rows = (totalblocks-1)/np_rows + 1
  max_loc_block_cols = (totalblocks-1)/np_cols + 1

  max_local_rows = max_loc_block_rows*nblk
  max_local_cols = max_loc_block_cols*nblk

  allocate(tmp(max(max_local_rows,max_local_cols)), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_complex", "tmp", istat, errorMessage)

  allocate(v_row(max_local_rows+1), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_complex", "v_row", istat, errorMessage)

  allocate(u_row(max_local_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_complex", "u_row", istat, errorMessage)

  allocate(v_col(max_local_cols), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_complex", "v_col", istat, errorMessage)

  allocate(u_col(max_local_cols), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_complex", "u_col", istat, errorMessage)

  tmp = 0
  v_row = 0
  u_row = 0
  v_col = 0
  u_col = 0

  allocate(vu_stored_rows(max_local_rows,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_complex", "vu_stored_rows", istat, errorMessage)

  allocate(uv_stored_cols(max_local_cols,2*max_stored_uv), stat=istat, errmsg=errorMessage)
  call check_alloc("tridiag_complex", "uv_stored_cols", istat, errorMessage)

  d_vec(:) = 0
  e_vec(:) = 0
  tau(:) = 0

  n_stored_vecs = 0

  l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a_mat
  l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a_mat
  if (my_prow==prow(na, nblk, np_rows) .and. my_pcol==pcol(na, nblk, np_cols)) &
    d_vec(na) = a_mat(l_rows,l_cols)

! main cycle of tridiagonalization
! in each step, 1 Householder vector is calculated
  do istep=na,3,-1

! Calculate number of local rows and columns of the still remaining matrix
! on the local processor

    l_rows = local_index(istep-1, my_prow, np_rows, nblk, -1)
    l_cols = local_index(istep-1, my_pcol, np_cols, nblk, -1)

! Calculate vector for Householder transformation on all procs
! owning column istep

    if (my_pcol==pcol(istep, nblk, np_cols)) then

! Get vector to be transformed; distribute last element and norm of
! remaining elements to all procs in current column

! copy l_cols + 1 column of A to v_row
      v_row(1:l_rows) = a_mat(1:l_rows,l_cols+1)

      if (n_stored_vecs>0 .and. l_rows>0) then
        aux(1:2*n_stored_vecs) = conjg(uv_stored_cols(l_cols+1,1:2*n_stored_vecs))
        call ZGEMV('N', l_rows, 2*n_stored_vecs, &
             CONE, vu_stored_rows, ubound(vu_stored_rows,dim=1), &
             aux, 1, CONE, v_row, 1)
      endif

      if (my_prow==prow(istep-1, nblk, np_rows)) then
        aux1(1) = dot_product(v_row(1:l_rows-1),v_row(1:l_rows-1))
        aux1(2) = v_row(l_rows)
      else
        aux1(1) = dot_product(v_row(1:l_rows),v_row(1:l_rows))
        aux1(2) = 0.
      endif

      call MPI_Allreduce(aux1, aux2, 2, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)

      vnorm2 = aux2(1)
      vrl = aux2(2)

! Householder transformation
      call hh_transform_complex_double(vrl, vnorm2, xf, tau(istep))
! Scale v_row and store Householder vector for back transformation

      v_row(1:l_rows) = v_row(1:l_rows) * xf
      if (my_prow==prow(istep-1, nblk, np_rows)) then
        v_row(l_rows) = 1.

! vrl is newly computed off-diagonal element of the final tridiagonal matrix
        e_vec(istep-1) = vrl
      endif

! store Householder vector for back transformation
      a_mat(1:l_rows,l_cols+1) = v_row(1:l_rows)

! add tau after the end of actuall v_row, to be broadcasted with it
      v_row(l_rows+1) = tau(istep)
    endif !(my_pcol==pcol(istep, nblk, np_cols))

! Broadcast the Householder vector (and tau) along columns
    call MPI_Bcast(v_row, l_rows+1, MPI_DOUBLE_COMPLEX, pcol(istep, nblk, np_cols), mpi_comm_cols, mpierr)

!recover tau, which has been broadcasted together with v_row
    tau(istep) = v_row(l_rows+1)

! Transpose Householder vector v_row -> v_col
    call elpa_transpose_vectors_complex_double(v_row, ubound(v_row,dim=1), mpi_comm_rows, &
         v_col, ubound(v_col,dim=1), mpi_comm_cols, 1, istep-1, 1, nblk)

! Calculate u = (A + VU**T + UV**T)*v

! For cache efficiency, we use only the upper half of the matrix tiles for this,
! thus the result is partly in u_col(:) and partly in u_row(:)

    u_col(1:l_cols) = 0
    u_row(1:l_rows) = 0
    if (l_rows>0 .and. l_cols>0) then
      do i=0,(istep-2)/tile_size
        l_col_beg = i*l_cols_per_tile+1
        l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
        if (l_col_end<l_col_beg) cycle
        do j=0,i
          l_row_beg = j*l_rows_per_tile+1
          l_row_end = min(l_rows,(j+1)*l_rows_per_tile)
          if (l_row_end<l_row_beg) cycle

          call ZGEMV('C', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
               CONE, a_mat(l_row_beg,l_col_beg), lda, &
               v_row(l_row_beg), 1, &
               CONE, u_col(l_col_beg), 1)

          if (i/=j) then
            call ZGEMV('N', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, &
                 CONE, a_mat(l_row_beg,l_col_beg), lda, &
                 v_col(l_col_beg), 1, CONE, u_row(l_row_beg), 1)
          endif

        enddo ! j=0,i
      enddo ! i=0,(istep-2)/tile_size

      if (n_stored_vecs>0) then
        call ZGEMV('C', l_rows, 2*n_stored_vecs, &
             CONE, vu_stored_rows, ubound(vu_stored_rows,dim=1), &
             v_row, 1, CZERO, aux, 1)
        call ZGEMV('N', l_cols, 2*n_stored_vecs, &
             CONE, uv_stored_cols, ubound(uv_stored_cols,dim=1), &
             aux, 1, CONE, u_col, 1)
      endif

    endif ! (l_rows>0 .and. l_cols>0)

! Sum up all u_row(:) parts along rows and add them to the u_col(:) parts
! on the processors containing the diagonal
! This is only necessary if u_row has been calculated, i.e. if the
! global tile size is smaller than the global remaining matrix

    if (tile_size < istep-1) then
      call elpa_reduce_add_vectors_complex_double (u_row, ubound(u_row,dim=1), mpi_comm_rows, &
           u_col, ubound(u_col,dim=1), mpi_comm_cols, istep-1, 1, nblk)
    endif

! Sum up all the u_col(:) parts, transpose u_col -> u_row

    if (l_cols>0) then
      tmp(1:l_cols) = u_col(1:l_cols)

      call MPI_Allreduce(tmp, u_col, l_cols, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)

    endif

    call elpa_transpose_vectors_complex_double(u_col, ubound(u_col,dim=1), mpi_comm_cols, &
         u_row, ubound(u_row,dim=1), mpi_comm_rows, &
         1, (istep-1), 1, nblk)

! calculate u**T * v (same as v**T * (A + VU**T + UV**T) * v )

    xc = 0
    if (l_cols>0) &
      xc = dot_product(v_col(1:l_cols),u_col(1:l_cols))

    call MPI_Allreduce(xc, vav, 1 , MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_cols, mpierr)

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

! We have calculated another Hauseholder vector, number of implicitly stored increased
    n_stored_vecs = n_stored_vecs+1

! If the limit of max_stored_uv is reached, calculate A + VU**T + UV**T
    if (n_stored_vecs==max_stored_uv .or. istep==3) then

      do i=0,(istep-2)/tile_size
        l_col_beg = i*l_cols_per_tile+1
        l_col_end = min(l_cols,(i+1)*l_cols_per_tile)
        l_row_beg = 1
        l_row_end = min(l_rows,(i+1)*l_rows_per_tile)
        if (l_col_end<l_col_beg .or. l_row_end<l_row_beg) &
          cycle

        call ZGEMM('N', 'C', l_row_end-l_row_beg+1, l_col_end-l_col_beg+1, 2*n_stored_vecs, &
             CONE, vu_stored_rows(l_row_beg,1), ubound(vu_stored_rows,dim=1), &
             uv_stored_cols(l_col_beg,1), ubound(uv_stored_cols,dim=1), &
             CONE, a_mat(l_row_beg,l_col_beg), lda)
      enddo

      n_stored_vecs = 0

    endif

    if (my_prow==prow(istep-1, nblk, np_rows) .and. my_pcol==pcol(istep-1, nblk, np_cols)) then
      if (n_stored_vecs>0) then
        a_mat(l_rows,l_cols) = a_mat(l_rows,l_cols) &
          + dot_product(vu_stored_rows(l_rows,1:2*n_stored_vecs),uv_stored_cols(l_cols,1:2*n_stored_vecs))
      end if
      d_vec(istep-1) = a_mat(l_rows,l_cols)

    endif

  enddo ! main cycle over istep=na,3,-1

! Store e_vec(1) and d_vec(1)

  if (my_pcol==pcol(2, nblk, np_cols)) then
    if (my_prow==prow(1, nblk, np_rows)) then
! We use last l_cols value of loop above
      vrl = a_mat(1,l_cols)
      call hh_transform_complex_double(vrl, 0.0_rk8, xf, tau(2))
      e_vec(1) = vrl

      a_mat(1,l_cols) = 1. ! for consistency only
    endif

    call MPI_Bcast(tau(2), 1, MPI_DOUBLE_COMPLEX, prow(1, nblk, np_rows), mpi_comm_rows, mpierr)

  endif

  call MPI_Bcast(tau(2), 1, MPI_DOUBLE_COMPLEX, pcol(2, nblk, np_cols), mpi_comm_cols, mpierr)

  if (my_prow==prow(1, nblk, np_rows) .and. my_pcol==pcol(1, nblk, np_cols)) then
    d_vec(1) = DREAL(a_mat(1,1))
  endif

  deallocate(tmp, v_row, u_row, v_col, u_col, vu_stored_rows, uv_stored_cols, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "tridiag_complex: error when deallocating tmp "//errorMessage
    stop
  endif

! distribute the arrays d_vec and e_vec to all processors

  allocate(tmp_real(na), stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "tridiag_complex: error when allocating tmp_real "//errorMessage
    stop
  endif

  tmp_real = d_vec
  call MPI_Allreduce(tmp_real, d_vec, na, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)
  tmp_real = d_vec
  call MPI_Allreduce(tmp_real, d_vec, na, MPI_REAL8 ,MPI_SUM, mpi_comm_cols, mpierr)
  tmp_real = e_vec
  call MPI_Allreduce(tmp_real, e_vec, na, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)
  tmp_real = e_vec
  call MPI_Allreduce(tmp_real, e_vec, na, MPI_REAL8, MPI_SUM, mpi_comm_cols, mpierr)

  deallocate(tmp_real, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "tridiag_complex: error when deallocating tmp_real "//errorMessage
    stop
  endif

end subroutine tridiag_complex_double

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
!> \param mpi_comm_rows        MPI-Communicator for rows
!>
!> \param mpi_comm_cols        MPI-Communicator for columns
!>
subroutine trans_ev_complex_double(na, nqc, a_mat, lda, tau, q_mat, ldq, nblk, mpi_comm_rows, mpi_comm_cols)

  use, intrinsic :: iso_c_binding
  use precision

  implicit none

  integer(kind=ik), intent(in) :: na, nqc, lda, ldq, nblk, mpi_comm_rows, mpi_comm_cols
  complex(kind=ck8), intent(in) :: tau(na)

  complex(kind=ck8), intent(inout) :: a_mat(lda,*), q_mat(ldq,*)

  integer(kind=ik) :: max_stored_rows

  complex(kind=ck8), parameter :: CZERO = (0.0_rk8,0.0_rk8), CONE = (1.0_rk8,0.0_rk8)

  integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
  integer(kind=ik) :: totalblocks, max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
  integer(kind=ik) :: l_cols, l_rows, l_colh, nstor
  integer(kind=ik) :: istep, n, nc, ic, ics, ice, nb, cur_pcol

  complex(kind=ck8), allocatable :: tmp1(:), tmp2(:), hvb(:), hvm(:,:)
  complex(kind=ck8), allocatable :: tmat(:,:), h1(:), h2(:)
  integer(kind=ik) :: istat
  character(200) :: errorMessage

  call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
  call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
  call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
  call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

  totalblocks = (na-1)/nblk + 1
  max_blocks_row = (totalblocks-1)/np_rows + 1
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1 ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  max_stored_rows = (63/nblk+1)*nblk

  allocate(tmat(max_stored_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_complex", "tmat", istat, errorMessage)

  allocate(h1(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_complex", "h1", istat, errorMessage)

  allocate(h2(max_stored_rows*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_complex", "h2", istat, errorMessage)

  allocate(tmp1(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_complex", "tmp1", istat, errorMessage)

  allocate(tmp2(max_local_cols*max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_complex", "tmp2", istat, errorMessage)

  allocate(hvb(max_local_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_complex", "hvb", istat, errorMessage)

  allocate(hvm(max_local_rows,max_stored_rows), stat=istat, errmsg=errorMessage)
  call check_alloc("trans_ev_complex", "hvm", istat, errorMessage)

  hvm = 0 ! Must be set to 0 !!!
  hvb = 0 ! Safety only

  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  nstor = 0

! In the complex case tau(2) /= 0
  if (my_prow == prow(1, nblk, np_rows)) then
    q_mat(1,1:l_cols) = q_mat(1,1:l_cols)*(CONE-tau(2))
  endif

  do istep=1,na,nblk

    ics = max(istep,3)
    ice = min(istep+nblk-1,na)
    if (ice<ics) cycle

    cur_pcol = pcol(istep, nblk, np_cols)

    nb = 0
    do ic=ics,ice

      l_colh = local_index(ic, my_pcol, np_cols, nblk, -1) ! Column of Householder vector
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector

      if (my_pcol==cur_pcol) then
        hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)
        if (my_prow==prow(ic-1, nblk, np_rows)) then
          hvb(nb+l_rows) = 1.
        endif
      endif

      nb = nb+l_rows
    enddo

    if (nb>0) &
      call MPI_Bcast(hvb, nb, MPI_DOUBLE_COMPLEX, cur_pcol, mpi_comm_cols, mpierr)

    nb = 0
    do ic=ics,ice
      l_rows = local_index(ic-1, my_prow, np_rows, nblk, -1) ! # rows of Householder vector
      hvm(1:l_rows,nstor+1) = hvb(nb+1:nb+l_rows)
      nstor = nstor+1
      nb = nb+l_rows
    enddo

! Please note: for smaller matix sizes (na/np_rows<=256), a value of 32 for nstor is enough!
    if (nstor+nblk>max_stored_rows .or. istep+nblk>na .or. (na/np_rows<=256 .and. nstor>=32)) then

! Calculate scalar products of stored vectors.
! This can be done in different ways, we use zherk

      tmat = 0
      if (l_rows>0) &
        call ZHERK('U', 'C', nstor, l_rows, CONE, hvm, ubound(hvm,dim=1), CZERO, tmat, max_stored_rows)

      nc = 0
      do n=1,nstor-1
        h1(nc+1:nc+n) = tmat(1:n,n+1)
        nc = nc+n
      enddo

      if (nc>0) call MPI_Allreduce(h1, h2, nc, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)

! Calculate triangular matrix T

      nc = 0
      tmat(1,1) = tau(ice-nstor+1)
      do n=1,nstor-1
        call ZTRMV('L', 'C', 'N', n, &
             tmat, max_stored_rows, &
             h2(nc+1),1)
        tmat(n+1,1:n) = -conjg(h2(nc+1:nc+n))*tau(ice-nstor+n+1)
        tmat(n+1,n+1) = tau(ice-nstor+n+1)
        nc = nc+n
      enddo

! Q = Q - V * T * V**T * Q

      if (l_rows>0) then
        call ZGEMM('C', 'N', nstor, l_cols, l_rows, &
             CONE, hvm, ubound(hvm,dim=1), &
             q_mat, ldq, &
             CZERO, tmp1 ,nstor)
      else !l_rows>0
        tmp1(1:l_cols*nstor) = 0
      endif !l_rows>0

      call MPI_Allreduce(tmp1, tmp2, nstor*l_cols, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)

! copy back tmp2 - after reduction...
      if (l_rows>0) then
! tmp2 = tmat * tmp2
        call ZTRMM('L', 'L', 'N', 'N', nstor, l_cols, &
             CONE, tmat, max_stored_rows, &
             tmp2, nstor)
!q_mat = q_mat - hvm*tmp2
         call ZGEMM('N', 'N', l_rows, l_cols, nstor, &
              -CONE, hvm, ubound(hvm,dim=1), &
              tmp2, nstor, &
              CONE, q_mat, ldq)

      endif ! l_rows>0

      nstor = 0
    endif

  enddo ! istep=1,na,nblk

  deallocate(tmat, h1, h2, tmp1, tmp2, hvb, hvm, stat=istat, errmsg=errorMessage)
  if (istat .ne. 0) then
    write(error_unit,*) "trans_ev_complex: error when deallocating hvb "//errorMessage
    stop
  endif

end subroutine trans_ev_complex_double

subroutine hh_transform_complex_double(alpha, xnorm_sq, xf, tau)

! Similar to LAPACK routine ZLARFP, but uses ||x||**2 instead of x(:)
! and returns the factor xf by which x has to be scaled.
! It also hasn't the special handling for numbers < 1.d-300 or > 1.d150
! since this would be expensive for the parallel implementation.

  use precision

  implicit none

  complex(kind=ck8), intent(inout) :: alpha
  real(kind=rk8), intent(in) :: xnorm_sq
  complex(kind=ck8), intent(out) :: xf, tau

  real(kind=rk8) :: ALPHR, ALPHI, BETA

  ALPHR = real( ALPHA, kind=rk8 )
  ALPHI = DIMAG( ALPHA )
  if ( XNORM_SQ==0. .AND. ALPHI==0. ) then

    if ( ALPHR>=0. ) then
      TAU = 0.
    else
      TAU = 2.
      ALPHA = -ALPHA
    endif
    XF = 0.

  else

    BETA = SIGN( SQRT( ALPHR**2 + ALPHI**2 + XNORM_SQ ), ALPHR )
    ALPHA = ALPHA + BETA
    if ( BETA<0 ) then
      BETA = -BETA
      TAU = -ALPHA / BETA
    else
      ALPHR = ALPHI * (ALPHI/real( ALPHA , kind=rk8))
      ALPHR = ALPHR + XNORM_SQ/real( ALPHA, kind=rk8 )

      TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA )
      ALPHA = DCMPLX( -ALPHR, ALPHI )
    end if
    XF = 1.0_rk8/ALPHA
    ALPHA = BETA
  endif

end subroutine hh_transform_complex_double

end module ELPA1_compute
