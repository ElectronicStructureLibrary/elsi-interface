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

module ELPA1_AUXILIARY

  use elpa_utilities

  implicit none

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
  public :: solve_tridi

  interface cholesky_real
    module procedure elpa_cholesky_real_double
  end interface

  interface invert_trm_real
    module procedure elpa_invert_trm_real_double
  end interface

  interface cholesky_complex
    module procedure elpa_cholesky_real_double
  end interface

  interface invert_trm_complex
    module procedure elpa_invert_trm_complex_double
  end interface

  interface mult_at_b_real
    module procedure elpa_mult_at_b_real_double
  end interface

  interface mult_ah_b_complex
    module procedure elpa_mult_ah_b_complex_double
  end interface

  interface solve_tridi
    module procedure elpa_solve_tridi_double
  end interface

  contains

!> \brief  cholesky_real_double: Cholesky factorization of a double-precision real symmetric matrix
!> \details
!>
!> \param  na                   Order of matrix
!> \param  a(lda,matrixCols)    Distributed matrix which should be factorized.
!>                              Distribution is like in Scalapack.
!>                              Only upper triangle needs to be set.
!>                              On return, the upper triangle contains the Cholesky factor
!>                              and the lower triangle is set to 0.
!> \param  lda                  Leading dimension of a
!> \param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \param succes                logical, reports success or failure

  function elpa_cholesky_real_double(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols) &
    result(success)

    use elpa1_compute
    use elpa_utilities
    use elpa_mpi
    use precision

    implicit none

    integer(kind=ik) :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols

    real(kind=rk8) :: a(lda,*)

    integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
    integer(kind=ik) :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
    integer(kind=ik) :: n, nc, i, info
    integer(kind=ik) :: lcs, lce, lrs, lre
    integer(kind=ik) :: tile_size, l_rows_tile, l_cols_tile

    real(kind=rk8), allocatable :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)

    logical :: success
    integer(kind=ik) :: istat
    character(200) :: errorMessage

    call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
    call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
    call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

    success = .true.

! Matrix is split into tiles; work is done only for tiles on the diagonal or above

    tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
    tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

    l_rows_tile = tile_size/np_rows ! local rows of a tile
    l_cols_tile = tile_size/np_cols ! local cols of a tile

    l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
    l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

    allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_real: error when allocating tmp1 "//errorMessage
      stop
    endif

    allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_real: error when allocating tmp2 "//errorMessage
      stop
    endif

    tmp1 = 0
    tmp2 = 0

    allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_real: error when allocating tmatr "//errorMessage
      stop
    endif

    allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_real: error when allocating tmatc "//errorMessage
      stop
    endif

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

          call DPOTRF('U', na-n+1, a(l_row1,l_col1), lda, info)

          if (info/=0) then
            write(error_unit,*) "elpa_cholesky_real: Error in dpotrf 1: ",info
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

          call DPOTRF('U', nblk, a(l_row1,l_col1), lda, info)

          if (info/=0) then
            write(error_unit,*) "elpa_cholesky_real: Error in dpotrf 2: ",info
            success = .false.
            return
          endif

          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif

        call MPI_Bcast(tmp1, nblk*(nblk+1)/2, MPI_REAL8, pcol(n, nblk, np_cols), mpi_comm_cols, mpierr)

        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo

        if (l_cols-l_colx+1>0) &
          call DTRSM('L', 'U', 'T', 'N', nblk, l_cols-l_colx+1, 1.0_rk8, tmp2, ubound(tmp2,dim=1), a(l_row1,l_colx), lda)
      endif

      do i=1,nblk

        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = a(l_row1+i-1,l_colx:l_cols)

        if (l_cols-l_colx+1>0) then
          call MPI_Bcast(tmatc(l_colx,i), l_cols-l_colx+1, MPI_REAL8, prow(n, nblk, np_rows), mpi_comm_rows, mpierr)
        endif

      enddo
! this has to be checked since it was changed substantially when doing type safe

      call elpa_transpose_vectors_real_double(tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
           tmatr, ubound(tmatr,dim=1), mpi_comm_rows, n, na, nblk, nblk)

      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle

        call DGEMM('N', 'T', lre-lrs+1, lce-lcs+1, nblk, -1.0_rk8, &
             tmatr(lrs,1), ubound(tmatr,dim=1), tmatc(lcs,1), ubound(tmatc,dim=1), &
             1.0_rk8, a(lrs,lcs), lda)

      enddo

    enddo

    deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_real: error when deallocating tmp1 "//errorMessage
      stop
    endif

! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

    do i=1,na
      if (my_pcol==pcol(i, nblk, np_cols)) then
! column i is on local processor
        l_col1 = local_index(i , my_pcol, np_cols, nblk, +1) ! local column number
        l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
        a(l_row1:l_rows,l_col1) = 0
      endif
    enddo

  end function elpa_cholesky_real_double

!> \brief  elpa_invert_trm_real_double: Inverts a double-precision real upper triangular matrix
!> \details
!> \param  na                   Order of matrix
!> \param  a(lda,matrixCols)    Distributed matrix which should be inverted
!>                              Distribution is like in Scalapack.
!>                              Only upper triangle needs to be set.
!>                              The lower triangle is not referenced.
!> \param  lda                  Leading dimension of a
!> \param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \result succes               logical, reports success or failure
  function elpa_invert_trm_real_double(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols) result(success)

    use precision
    use elpa1_compute
    use elpa_utilities
    use elpa_mpi

    implicit none

    integer(kind=ik) :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols

    real(kind=rk8) :: a(lda,*)

    integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
    integer(kind=ik) :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
    integer(kind=ik) :: n, nc, i, info, ns, nb

    real(kind=rk8), allocatable :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)

    logical :: success
    integer(kind=ik) :: istat
    character(200) :: errorMessage

    call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
    call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
    call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

    success = .true.

    l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
    l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

    allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_real: error when allocating tmp1 "//errorMessage
      stop
    endif

    allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_real: error when allocating tmp2 "//errorMessage
      stop
    endif

    tmp1 = 0
    tmp2 = 0

    allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_real: error when allocating tmat1 "//errorMessage
      stop
    endif

    allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_real: error when allocating tmat2 "//errorMessage
      stop
    endif

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

          call DTRTRI('U', 'N', nb, a(l_row1,l_col1), lda, info)

          if (info/=0) then
            write(error_unit,*) "elpa_invert_trm_real: Error in DTRTRI"
            success = .false.
            return
          endif

          nc = 0
          do i=1,nb
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif

        call MPI_Bcast(tmp1, nb*(nb+1)/2, MPI_REAL8, pcol(n, nblk, np_cols), mpi_comm_cols, mpierr)

        nc = 0
        do i=1,nb
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo

        if (l_cols-l_colx+1>0) &
           call DTRMM('L', 'U', 'N', 'N', nb, l_cols-l_colx+1, 1.0_rk8, tmp2, ubound(tmp2,dim=1), a(l_row1,l_colx), lda)

        if (l_colx<=l_cols) tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
        if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

      endif

      if (l_row1>1) then
        if (my_pcol==pcol(n, nblk, np_cols)) then
          tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
          a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
        endif

        do i=1,nb

          call MPI_Bcast(tmat1(1,i), l_row1-1, MPI_REAL8, pcol(n, nblk, np_cols), mpi_comm_cols, mpierr)

        enddo
      endif

      if (l_cols-l_col1+1>0) &
        call MPI_Bcast(tmat2(1,l_col1), (l_cols-l_col1+1)*nblk, MPI_REAL8, prow(n, nblk, np_rows), mpi_comm_rows, mpierr)

      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call DGEMM('N', 'N', l_row1-1, l_cols-l_col1+1, nb, -1.0_rk8, &
             tmat1, ubound(tmat1,dim=1), tmat2(1,l_col1), ubound(tmat2,dim=1), &
             1.0_rk8, a(1,l_col1), lda)

    enddo

    deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_real: error when deallocating tmp1 "//errorMessage
      stop
    endif

  end function elpa_invert_trm_real_double

!> \brief  elpa_cholesky_complex_double: Cholesky factorization of a double-precision complex hermitian matrix
!> \details
!> \param  na                   Order of matrix
!> \param  a(lda,matrixCols)    Distributed matrix which should be factorized.
!>                              Distribution is like in Scalapack.
!>                              Only upper triangle needs to be set.
!>                              On return, the upper triangle contains the Cholesky factor
!>                              and the lower triangle is set to 0.
!> \param  lda                  Leading dimension of a
!> \param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \result succes               logical, reports success or failure
  function elpa_cholesky_complex_double(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols) result(success)

    use elpa1_compute
    use elpa_utilities
    use elpa_mpi
    use precision

    implicit none

    integer(kind=ik) :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols

    complex(kind=ck8) :: a(lda,*)

    integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
    integer(kind=ik) :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
    integer(kind=ik) :: n, nc, i, info
    integer(kind=ik) :: lcs, lce, lrs, lre
    integer(kind=ik) :: tile_size, l_rows_tile, l_cols_tile

    complex(kind=ck8), allocatable :: tmp1(:), tmp2(:,:), tmatr(:,:), tmatc(:,:)

    logical :: success
    integer(kind=ik) :: istat
    character(200) :: errorMessage

    success = .true.

    call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
    call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
    call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

! Matrix is split into tiles; work is done only for tiles on the diagonal or above

    tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
    tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

    l_rows_tile = tile_size/np_rows ! local rows of a tile
    l_cols_tile = tile_size/np_cols ! local cols of a tile

    l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
    l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

    allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_complex: error when allocating tmp1 "//errorMessage
      stop
    endif

    allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_complex: error when allocating tmp2 "//errorMessage
      stop
    endif

    tmp1 = 0
    tmp2 = 0

    allocate(tmatr(l_rows,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_complex: error when allocating tmatr "//errorMessage
      stop
    endif

    allocate(tmatc(l_cols,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_complex: error when allocating tmatc "//errorMessage
      stop
    endif

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

          call ZPOTRF('U', na-n+1, a(l_row1,l_col1),lda, info)

          if (info/=0) then
            write(error_unit,*) "elpa_cholesky_complex: Error in zpotrf"
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

          call ZPOTRF('U', nblk, a(l_row1,l_col1),lda, info)

          if (info/=0) then
            write(error_unit,*) "elpa_cholesky_complex: Error in zpotrf"
            success = .false.
            return
          endif

          nc = 0
          do i=1,nblk
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif

        call MPI_Bcast(tmp1, nblk*(nblk+1)/2, MPI_DOUBLE_COMPLEX, pcol(n, nblk, np_cols), mpi_comm_cols, mpierr)

        nc = 0
        do i=1,nblk
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo

        if (l_cols-l_colx+1>0) &
           call ZTRSM('L', 'U', 'C', 'N', nblk, l_cols-l_colx+1, (1.0_rk8,0.0_rk8), tmp2, ubound(tmp2,dim=1), &
                a(l_row1,l_colx), lda)

      endif

      do i=1,nblk

        if (my_prow==prow(n, nblk, np_rows)) tmatc(l_colx:l_cols,i) = conjg(a(l_row1+i-1,l_colx:l_cols))

        if (l_cols-l_colx+1>0) &
          call MPI_Bcast(tmatc(l_colx,i), l_cols-l_colx+1, MPI_DOUBLE_COMPLEX, prow(n, nblk, np_rows), &
               mpi_comm_rows, mpierr)

      enddo
! this has to be checked since it was changed substantially when doing type safe

      call elpa_transpose_vectors_complex_double (tmatc, ubound(tmatc,dim=1), mpi_comm_cols, &
           tmatr, ubound(tmatr,dim=1), mpi_comm_rows, n, na, nblk, nblk)

      do i=0,(na-1)/tile_size
        lcs = max(l_colx,i*l_cols_tile+1)
        lce = min(l_cols,(i+1)*l_cols_tile)
        lrs = l_rowx
        lre = min(l_rows,(i+1)*l_rows_tile)
        if (lce<lcs .or. lre<lrs) cycle

        call ZGEMM('N', 'C', lre-lrs+1, lce-lcs+1, nblk, (-1.0_rk8,0.0_rk8), &
             tmatr(lrs,1), ubound(tmatr,dim=1), tmatc(lcs,1), ubound(tmatc,dim=1), &
             (1.0_rk8,0.0_rk8), a(lrs,lcs), lda)

      enddo

    enddo

    deallocate(tmp1, tmp2, tmatr, tmatc, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_cholesky_complex: error when deallocating tmatr "//errorMessage
      stop
    endif

! Set the lower triangle to 0, it contains garbage (form the above matrix multiplications)

    do i=1,na
      if (my_pcol==pcol(i, nblk, np_cols)) then
! column i is on local processor
        l_col1 = local_index(i , my_pcol, np_cols, nblk, +1) ! local column number
        l_row1 = local_index(i+1, my_prow, np_rows, nblk, +1) ! first row below diagonal
        a(l_row1:l_rows,l_col1) = 0
      endif
    enddo

  end function elpa_cholesky_complex_double

!> \brief  elpa_invert_trm_complex_double: Inverts a double-precision complex upper triangular matrix
!> \details
!> \param  na                   Order of matrix
!> \param  a(lda,matrixCols)    Distributed matrix which should be inverted
!>                              Distribution is like in Scalapack.
!>                              Only upper triangle needs to be set.
!>                              The lower triangle is not referenced.
!> \param  lda                  Leading dimension of a
!> \param  nblk                 blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \result succes               logical, reports success or failure

  function elpa_invert_trm_complex_double(na, a, lda, nblk, mpi_comm_rows, mpi_comm_cols) result(success)

    use elpa1_compute
    use elpa_utilities
    use elpa_mpi
    use precision

    implicit none

    integer(kind=ik) :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols

    complex(kind=ck8) :: a(lda,*)

    integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
    integer(kind=ik) :: l_cols, l_rows, l_col1, l_row1, l_colx, l_rowx
    integer(kind=ik) :: n, nc, i, info, ns, nb

    complex(kind=ck8), allocatable :: tmp1(:), tmp2(:,:), tmat1(:,:), tmat2(:,:)

    logical :: success
    integer(kind=ik) :: istat
    character(200) :: errorMessage

    call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
    call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
    call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

    success = .true.

    l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
    l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local cols of a

    allocate(tmp1(nblk*nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_complex: error when allocating tmp1 "//errorMessage
      stop
    endif

    allocate(tmp2(nblk,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_complex: error when allocating tmp2 "//errorMessage
      stop
    endif

    tmp1 = 0
    tmp2 = 0

    allocate(tmat1(l_rows,nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_complex: error when allocating tmat1 "//errorMessage
      stop
    endif

    allocate(tmat2(nblk,l_cols), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_complex: error when allocating tmat2 "//errorMessage
      stop
    endif

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

          call ZTRTRI('U', 'N', nb, a(l_row1,l_col1), lda, info)

          if (info/=0) then
            write(error_unit,*) "elpa_invert_trm_complex: Error in ZTRTRI"
            success = .false.
            return
          endif

          nc = 0
          do i=1,nb
            tmp1(nc+1:nc+i) = a(l_row1:l_row1+i-1,l_col1+i-1)
            nc = nc+i
          enddo
        endif

        call MPI_Bcast(tmp1, nb*(nb+1)/2, MPI_DOUBLE_COMPLEX, pcol(n, nblk, np_cols), mpi_comm_cols, mpierr)

        nc = 0
        do i=1,nb
          tmp2(1:i,i) = tmp1(nc+1:nc+i)
          nc = nc+i
        enddo

        if (l_cols-l_colx+1>0) &
          call ZTRMM('L', 'U', 'N', 'N', nb, l_cols-l_colx+1, (1.0_rk8,0.0_rk8), tmp2, ubound(tmp2,dim=1), a(l_row1,l_colx), lda)

        if (l_colx<=l_cols) tmat2(1:nb,l_colx:l_cols) = a(l_row1:l_row1+nb-1,l_colx:l_cols)
        if (my_pcol==pcol(n, nblk, np_cols)) tmat2(1:nb,l_col1:l_col1+nb-1) = tmp2(1:nb,1:nb) ! tmp2 has the lower left triangle 0

      endif

      if (l_row1>1) then
        if (my_pcol==pcol(n, nblk, np_cols)) then
          tmat1(1:l_row1-1,1:nb) = a(1:l_row1-1,l_col1:l_col1+nb-1)
          a(1:l_row1-1,l_col1:l_col1+nb-1) = 0
        endif

        do i=1,nb

          call MPI_Bcast(tmat1(1,i), l_row1-1, MPI_DOUBLE_COMPLEX, pcol(n, nblk, np_cols), mpi_comm_cols, mpierr)

        enddo
      endif

      if (l_cols-l_col1+1>0) &
        call MPI_Bcast(tmat2(1,l_col1), (l_cols-l_col1+1)*nblk, MPI_DOUBLE_COMPLEX, prow(n, nblk, np_rows), &
             mpi_comm_rows, mpierr)

      if (l_row1>1 .and. l_cols-l_col1+1>0) &
        call ZGEMM('N', 'N', l_row1-1, l_cols-l_col1+1, nb, (-1.0_rk8,0.0_rk8), &
             tmat1, ubound(tmat1,dim=1), tmat2(1,l_col1), ubound(tmat2,dim=1), &
             (1.0_rk8,0.0_rk8), a(1,l_col1), lda)

    enddo

    deallocate(tmp1, tmp2, tmat1, tmat2, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_invert_trm_complex: error when deallocating tmp1 "//errorMessage
      stop
    endif

  end function elpa_invert_trm_complex_double

!> \brief  mult_at_b_real_double: Performs C : = A**T * B
!>         where   A is a square matrix (na,na) which is optionally upper or lower triangular
!>                 B is a (na,ncb) matrix
!>                 C is a (na,ncb) matrix where optionally only the upper or lower
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
!>                                           written to a certain extent, i.e. one shouldn't rely on the content there!
!> \param na                    Number of rows/columns of A, number of rows of B and C
!> \param ncb                   Number of columns of B and C
!> \param a                     matrix a
!> \param lda                   leading dimension of matrix a
!> \param b                     matrix b
!> \param ldb                   leading dimension of matrix b
!> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \param c                     matrix c
!> \param ldc                   leading dimension of matrix c
!> \result success

  function elpa_mult_at_b_real_double(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, &
    mpi_comm_rows, mpi_comm_cols, c, ldc) result(success)

    use elpa1_compute
    use elpa_mpi
    use precision

    implicit none

    character*1 :: uplo_a, uplo_c

    integer(kind=ik), intent(in) :: na, lda, ldb, ldc, nblk
    integer(kind=ik) :: ncb, mpi_comm_rows, mpi_comm_cols
    real(kind=rk8) :: a(lda,*), b(ldb,*), c(ldc,*)

    integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
    integer(kind=ik) :: l_cols, l_rows, l_rows_np
    integer(kind=ik) :: np, n, nb, nblk_mult, lrs, lre, lcs, lce
    integer(kind=ik) :: gcol_min, gcol, goff
    integer(kind=ik) :: nstor, nr_done, noff, np_bc, n_aux_bc, nvals
    integer(kind=ik), allocatable :: lrs_save(:), lre_save(:)

    logical :: a_lower, a_upper, c_lower, c_upper

    real(kind=rk8), allocatable :: aux_mat(:,:), aux_bc(:), tmp1(:,:), tmp2(:,:)
    integer(kind=ik) :: istat
    character(200) :: errorMessage
    logical :: success

    success = .true.

    call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
    call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
    call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

    l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and b
    l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

! Block factor for matrix multiplications, must be a multiple of nblk

    if (na/np_rows<=256) then
      nblk_mult = (31/nblk+1)*nblk
    else
      nblk_mult = (63/nblk+1)*nblk
    endif

    allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_mult_at_b_real: error when allocating aux_mat "//errorMessage
      stop
    endif

    allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_mult_at_b_real: error when allocating aux_bc "//errorMessage
      stop
    endif

    allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_mult_at_b_real: error when allocating lrs_save "//errorMessage
      stop
    endif

    allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_mult_at_b_real: error when allocating lre_save "//errorMessage
      stop
    endif

    a_lower = .false.
    a_upper = .false.
    c_lower = .false.
    c_upper = .false.

    if (uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
    if (uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
    if (uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
    if (uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

! Build up the result matrix by processor rows

    do np = 0, np_rows-1

! In this turn, procs of row np assemble the result

      l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

      nr_done = 0 ! Number of rows done
      aux_mat = 0
      nstor = 0 ! Number of columns stored in aux_mat

! Loop over the blocks on row np

      do nb=0,(l_rows_np-1)/nblk

        goff = nb*np_rows + np ! Global offset in blocks corresponding to nb

! Get the processor column which owns this block (A is transposed, so we need the column)
! and the offset in blocks within this column.
! The corresponding block column in A is then broadcast to all for multiplication with B

        np_bc = mod(goff,np_cols)
        noff = goff/np_cols
        n_aux_bc = 0

! Gather up the complete block column of A on the owner

        do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

          gcol = goff*nblk + n ! global column corresponding to n
          if (nstor==0 .and. n==1) gcol_min = gcol

          lrs = 1 ! 1st local row number for broadcast
          lre = l_rows ! last local row number for broadcast
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

        call MPI_Bcast(aux_bc, n_aux_bc, MPI_REAL8, np_bc, mpi_comm_cols, mpierr)

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

          lrs = 1 ! 1st local row number for multiply
          lre = l_rows ! last local row number for multiply
          if (a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
          if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

          lcs = 1 ! 1st local col number for multiply
          lce = l_cols ! last local col number for multiply
          if (c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
          if (c_lower) lce = min(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

          if (lcs<=lce) then
            allocate(tmp1(nstor,lcs:lce),tmp2(nstor,lcs:lce), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              write(error_unit,*) "elpa_mult_at_b_real: error when allocating tmp1 "//errorMessage
              stop
            endif

            if (lrs<=lre) then

              call DGEMM('T', 'N', nstor, lce-lcs+1, lre-lrs+1, 1.0_rk8, aux_mat(lrs,1), ubound(aux_mat,dim=1), &
                   b(lrs,lcs), ldb, 0.0_rk8, tmp1, nstor)

            else
              tmp1 = 0
            endif

! Sum up the results and send to processor row np

            call MPI_Reduce(tmp1, tmp2, nstor*(lce-lcs+1), MPI_REAL8, MPI_SUM, np, mpi_comm_rows, mpierr)

! Put the result into C
            if (my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

            deallocate(tmp1,tmp2, stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              write(error_unit,*) "elpa_mult_at_b_real: error when deallocating tmp1 "//errorMessage
              stop
            endif

          endif

          nr_done = nr_done+nstor
          nstor=0
          aux_mat(:,:)=0
        endif
      enddo
    enddo

    deallocate(aux_mat, aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
       write(error_unit,*) "elpa_mult_at_b_real: error when deallocating aux_mat "//errorMessage
       stop
    endif

  end function elpa_mult_at_b_real_double

!> \brief  elpa_mult_ah_b_complex_double: Performs C : = A**H * B
!>         where   A is a square matrix (na,na) which is optionally upper or lower triangular
!>                 B is a (na,ncb) matrix
!>                 C is a (na,ncb) matrix where optionally only the upper or lower
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
!>                                           written to a certain extent, i.e. one shouldn't rely on the content there!
!> \param na                    Number of rows/columns of A, number of rows of B and C
!> \param ncb                   Number of columns of B and C
!> \param a                     matrix a
!> \param lda                   leading dimension of matrix a
!> \param b                     matrix b
!> \param ldb                   leading dimension of matrix b
!> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
!> \param  mpi_comm_rows        MPI communicator for rows
!> \param  mpi_comm_cols        MPI communicator for columns
!> \param c                     matrix c
!> \param ldc                   leading dimension of matrix c
!> \result success

  function elpa_mult_ah_b_complex_double(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, &
    mpi_comm_rows, mpi_comm_cols, c, ldc) result(success)

    use precision
    use elpa1_compute
    use elpa_mpi

    implicit none

    character*1 :: uplo_a, uplo_c
    integer(kind=ik), intent(in) :: lda, ldb, ldc
    integer(kind=ik) :: na, ncb, nblk, mpi_comm_rows, mpi_comm_cols

    complex(kind=ck8):: a(lda,*), b(ldb,*), c(ldc,*)

    integer(kind=ik) :: my_prow, my_pcol, np_rows, np_cols, mpierr
    integer(kind=ik) :: l_cols, l_rows, l_rows_np
    integer(kind=ik) :: np, n, nb, nblk_mult, lrs, lre, lcs, lce
    integer(kind=ik) :: gcol_min, gcol, goff
    integer(kind=ik) :: nstor, nr_done, noff, np_bc, n_aux_bc, nvals
    integer(kind=ik), allocatable :: lrs_save(:), lre_save(:)

    logical :: a_lower, a_upper, c_lower, c_upper

    complex(kind=ck8), allocatable :: aux_mat(:,:), aux_bc(:), tmp1(:,:), tmp2(:,:)
    integer(kind=ik) :: istat
    character(200) :: errorMessage
    logical :: success

    success = .true.

    call MPI_Comm_rank(mpi_comm_rows,my_prow,mpierr)
    call MPI_Comm_size(mpi_comm_rows,np_rows,mpierr)
    call MPI_Comm_rank(mpi_comm_cols,my_pcol,mpierr)
    call MPI_Comm_size(mpi_comm_cols,np_cols,mpierr)

    l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a and b
    l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

! Block factor for matrix multiplications, must be a multiple of nblk

    if (na/np_rows<=256) then
      nblk_mult = (31/nblk+1)*nblk
    else
      nblk_mult = (63/nblk+1)*nblk
    endif

    allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_mult_ah_b_complex: error when allocating aux_mat "//errorMessage
      stop
    endif

    allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_mult_ah_b_complex: error when allocating aux_bc "//errorMessage
      stop
    endif

    allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
       write(error_unit,*) "elpa_mult_ah_b_complex: error when allocating lrs_save "//errorMessage
       stop
    endif

    allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      write(error_unit,*) "elpa_mult_ah_b_complex: error when allocating lre_save "//errorMessage
      stop
    endif

    a_lower = .false.
    a_upper = .false.
    c_lower = .false.
    c_upper = .false.

    if (uplo_a=='u' .or. uplo_a=='U') a_upper = .true.
    if (uplo_a=='l' .or. uplo_a=='L') a_lower = .true.
    if (uplo_c=='u' .or. uplo_c=='U') c_upper = .true.
    if (uplo_c=='l' .or. uplo_c=='L') c_lower = .true.

! Build up the result matrix by processor rows

    do np = 0, np_rows-1

! In this turn, procs of row np assemble the result

      l_rows_np = local_index(na, np, np_rows, nblk, -1) ! local rows on receiving processors

      nr_done = 0 ! Number of rows done
      aux_mat = 0
      nstor = 0 ! Number of columns stored in aux_mat

! Loop over the blocks on row np

      do nb=0,(l_rows_np-1)/nblk

        goff = nb*np_rows + np ! Global offset in blocks corresponding to nb

! Get the processor column which owns this block (A is transposed, so we need the column)
! and the offset in blocks within this column.
! The corresponding block column in A is then broadcast to all for multiplication with B

        np_bc = mod(goff,np_cols)
        noff = goff/np_cols
        n_aux_bc = 0

! Gather up the complete block column of A on the owner

        do n = 1, min(l_rows_np-nb*nblk,nblk) ! Loop over columns to be broadcast

          gcol = goff*nblk + n ! global column corresponding to n
          if (nstor==0 .and. n==1) gcol_min = gcol

          lrs = 1 ! 1st local row number for broadcast
          lre = l_rows ! last local row number for broadcast
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

        call MPI_Bcast(aux_bc, n_aux_bc, MPI_DOUBLE_COMPLEX, np_bc, mpi_comm_cols, mpierr)

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

          lrs = 1 ! 1st local row number for multiply
          lre = l_rows ! last local row number for multiply
          if (a_lower) lrs = local_index(gcol_min, my_prow, np_rows, nblk, +1)
          if (a_upper) lre = local_index(gcol, my_prow, np_rows, nblk, -1)

          lcs = 1 ! 1st local col number for multiply
          lce = l_cols ! last local col number for multiply
          if (c_upper) lcs = local_index(gcol_min, my_pcol, np_cols, nblk, +1)
          if (c_lower) lce = min(local_index(gcol, my_pcol, np_cols, nblk, -1),l_cols)

          if (lcs<=lce) then
            allocate(tmp1(nstor,lcs:lce),tmp2(nstor,lcs:lce), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              write(error_unit,*) "elpa_mult_ah_b_complex: error when allocating tmp1 "//errorMessage
              stop
            endif

            if (lrs<=lre) then

              call ZGEMM('C', 'N', nstor, lce-lcs+1, lre-lrs+1, (1.0_rk8,0.0_rk8), aux_mat(lrs,1), ubound(aux_mat,dim=1), &
                   b(lrs,lcs), ldb, (0.0_rk8,0.0_rk8), tmp1, nstor)

            else
              tmp1 = 0
            endif

! Sum up the results and send to processor row np

            call MPI_Reduce(tmp1, tmp2, nstor*(lce-lcs+1), MPI_DOUBLE_COMPLEX, MPI_SUM, np, mpi_comm_rows, mpierr)

! Put the result into C
            if (my_prow==np) c(nr_done+1:nr_done+nstor,lcs:lce) = tmp2(1:nstor,lcs:lce)

            deallocate(tmp1,tmp2, stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              write(error_unit,*) "elpa_mult_ah_b_complex: error when deallocating tmp1 "//errorMessage
              stop
            endif

          endif

          nr_done = nr_done+nstor
          nstor=0
          aux_mat(:,:)=0
        endif
      enddo
    enddo

    deallocate(aux_mat, aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
        write(error_unit,*) "elpa_mult_ah_b_complex: error when deallocating aux_mat "//errorMessage
        stop
    endif

  end function elpa_mult_ah_b_complex_double

!> \brief  elpa_solve_tridi_double: Solve tridiagonal eigensystem for a double-precision matrix with divide and conquer method
!> \details
!>
!> \param na                    Matrix dimension
!> \param nev                   number of eigenvalues/vectors to be computed
!> \param d                     array d(na) on input diagonal elements of tridiagonal matrix, on
!>                              output the eigenvalues in ascending order
!> \param e                     array e(na) on input subdiagonal elements of matrix, on exit destroyed
!> \param q                     on exit : matrix q(ldq,matrixCols) contains the eigenvectors
!> \param ldq                   leading dimension of matrix q
!> \param nblk                  blocksize of cyclic distribution, must be the same in both directions!
!> \param mpi_comm_rows         MPI communicator for rows
!> \param mpi_comm_cols         MPI communicator for columns
!> \result success              logical, .true. on success, else .false.

  function elpa_solve_tridi_double(na, nev, d, e, q, ldq, nblk, mpi_comm_rows, mpi_comm_cols) &
    result(success)

    use elpa1_compute, solve_tridi_double_private => solve_tridi_double
    use precision

    implicit none

    integer(kind=ik) :: na, nev, ldq, nblk, mpi_comm_rows, mpi_comm_cols
    real(kind=rk8) :: d(na), e(na)

    real(kind=rk8) :: q(ldq,*)

    logical :: success

    success = .false.

    call solve_tridi_double_private(na, nev, d, e, q, ldq, nblk, &
         mpi_comm_rows, mpi_comm_cols, success)

  end function

end module elpa1_auxiliary
