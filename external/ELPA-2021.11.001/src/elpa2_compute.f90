










!    This file is part of ELPA.
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
! ELPA2 -- 2-stage solver for ELPA
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
! Author: Andreas Marek, MPCDF


module ELPA2_compute

! Version 1.1.2, 2011-02-21

  use ELPA_utilities
  USE ELPA1_compute
  use elpa_pdgeqrf
  use precision
  use elpa_mpi
  use aligned_mem

  implicit none

  PRIVATE ! By default, all routines contained are private

  public :: bandred_real_double
  public :: tridiag_band_real_double
  public :: trans_ev_tridi_to_band_real_double
  public :: trans_ev_band_to_full_real_double

  public :: bandred_real_single
  public :: tridiag_band_real_single
  public :: trans_ev_tridi_to_band_real_single
  public :: trans_ev_band_to_full_real_single

  public :: bandred_complex_double
  public :: tridiag_band_complex_double
  public :: trans_ev_tridi_to_band_complex_double
  public :: trans_ev_band_to_full_complex_double

  public :: bandred_complex_single
  public :: tridiag_band_complex_single
  public :: trans_ev_tridi_to_band_complex_single
  public :: trans_ev_band_to_full_complex_single
  public :: band_band_real_double
!  public :: divide_band

  integer(kind=ik), public :: which_qr_decomposition = 1     ! defines, which QR-decomposition algorithm will be used
                                                    ! 0 for unblocked
                                                    ! 1 for blocked (maxrank: nblk)
  contains

! real double precision





















!cannot use "../src/elpa2/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



subroutine bandred_&
&real&
&_&
&double &
(obj, na, a_mat, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, tmat, &
wantDebug, useGPU, success, &
useQR, &
max_threads, isSkewsymmetric)

!-------------------------------------------------------------------------------
!  bandred_real/complex: Reduces a distributed symmetric matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a_mat(lda,matrixCols)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a_mat(:,:) must be set completely (upper and lower half)
!              a_mat(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,numBlocks)    where numBlocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa1_compute
  use precision
  use elpa_blas_interfaces
  use elpa_scalapack_interfaces
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

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                            :: na, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols

  real(kind=rck)                     :: a_mat(lda,*)
  real(kind=rck)                     :: tmat(nbw,nbw,*)

  real(kind=rk)                               :: eps
  logical, intent(in)                         :: useGPU
  logical, intent(in)                         :: isSkewsymmetric
  character(20)                               :: gpuString

  integer(kind=ik)                            :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                      :: mpierr,  my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                            :: l_cols, l_rows
  integer(kind=ik)                            :: vmrCols
  integer(kind=ik)                            :: i, j, lcs, lce, lre, lc, lr, cur_pcol, n_cols, nrow
  integer(kind=ik)                            :: istep, ncol, lch, lcx, nlc
  integer(kind=ik)                            :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rk)                              :: vnorm2
  real(kind=rck)                    :: xf, aux1(nbw), aux2(nbw), vrl, tau
  real(kind=rck)                    :: vav(nbw,nbw)

  real(kind=rck), allocatable :: tmpGPU(:)
  real(kind=rck), pointer     :: vmrGPU(:), umcGPU(:)
  real(kind=rck), allocatable :: tmpCPU(:,:), vmrCPU(:,:), umcCPU(:,:)
  real(kind=rck), allocatable :: vr(:)

  ! needed for blocked QR decomposition
  integer(kind=ik)                            :: PQRPARAM(11), work_size
  real(kind=rk)                               :: dwork_size(1)
  real(kind=rk), allocatable                  :: work_blocked(:), tauvector(:), blockheuristic(:)
  integer(kind=C_intptr_T)                    :: a_dev, vmr_dev, umc_dev, tmat_dev, vav_dev
  type(c_ptr)                                 :: vmr_host, umc_host
  !integer(kind=ik), external                  :: numroc -> use elpa_scalapack
  integer(kind=ik)                            :: ierr
  integer(kind=ik)                            :: cur_l_rows, cur_l_cols, vmr_size, umc_size
  integer(kind=ik)                            :: l_rows2, vmr_size2, umc_size2
  integer(kind=c_intptr_t)                    :: lc_start, lc_end
  integer(kind=ik)                            :: lr_end
  integer(kind=ik)                            :: na_cols
  integer(kind=BLAS_KIND)                     :: na_colsBLAS

  logical, intent(in)                         :: wantDebug
  logical, intent(out)                        :: success
  logical                                     :: successGPU
  integer(kind=ik)                            :: istat
  character(200)                              :: errorMessage
  integer(kind=ik)                            :: min_tile_size, error

  logical, intent(in)                         :: useQR
  integer(kind=ik)                            :: mystart, myend, m_way, n_way, work_per_thread, m_id, n_id, n_threads, &
                                                ii, pp
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                    &double&
                                                                    &_&
                                                                    &real

  logical                                     :: useGPU_reduction_lower_block_to_tridiagonal
  integer(kind=ik), intent(in)                :: max_threads
  logical                                     :: do_memcpy
  integer(kind=ik)                            :: i_blk,blk_off, blk_end
  logical                                     :: useIntelGPU


  integer(kind=MPI_KIND)                      :: bcast_request, allreduce_request1, allreduce_request2, &
                                                 allreduce_request3, allreduce_request4, allreduce_request5, &
                                                 allreduce_request6

  logical                                     :: useNonBlockingCollectivesCols
  logical                                     :: useNonBlockingCollectivesRows
  integer(kind=c_int)                         :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("bandred_&
  &real&
  &" // &
  "_double" // &
  gpuString )

  call obj%get("nbc_row_elpa2_full_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_bandred. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_full_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_bandred. Aborting..."
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

  useGPU_reduction_lower_block_to_tridiagonal = .false.

  if (useGPU .and. .not.(useIntelGPU)) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
    if (useQR) then
      !in this case switch off GPU usage for step "reduce current block to lower triangular form"
      ! since this is done by QR decomposition
      useGPU_reduction_lower_block_to_tridiagonal = .false.
    endif
  endif

  if (useIntelGPU) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
    if (useQR) then
      !in this case switch off GPU usage for step "reduce current block to lower triangular form"
      ! since this is done by QR decomposition
      useGPU_reduction_lower_block_to_tridiagonal = .false.
    endif
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")
  success = .true.


  ! Semibandwith nbw must be a multiple of blocksize nblk
  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_bandred_&
                             &real&
                             &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_bandred_&
                             &real&
                             &: ELPA2 works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  ! na_rows in used nowhere; only na_cols
  if (useGPU .and. .not.(useIntelGPU)) then
    na_colsBLAS = numroc(int(na,kind=BLAS_KIND), int(nblk,kind=BLAS_KIND), &
                         int(my_pcol,kind=BLAS_KIND), 0_BLAS_KIND, int(np_cols,kind=BLAS_KIND))
    na_cols = int(na_colsBLAS,kind=c_int)

    ! Here we convert the regular host array into a pinned host array
    successGPU = gpu_malloc(a_dev, lda*na_cols* size_of_datatype)
    call check_alloc_GPU_f("bandred: a_dev", 335,  successGPU)

    successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), &
                  nbw * nbw * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: vav", 340,  successGPU)

    successGPU = gpu_malloc(vav_dev, nbw*nbw* size_of_datatype)
    call check_alloc_GPU_f("bandred: vav_dev", 343,  successGPU)
  endif ! useGPU

  !if (useIntelGPU) then
  !  ! needed later when explicit copy
  !endif ! useIntelGPU

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

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

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  if (useQR) then

    if (which_qr_decomposition == 1) then
      call qr_pqrparam_init(obj,pqrparam(1:11),    nblk,'M',0,   nblk,'M',0,   nblk,'M',1,'s')
      allocate(tauvector(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: tauvector", 378,  istat,  errorMessage)

      allocate(blockheuristic(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: blockheuristic", 381,  istat,  errorMessage)

      l_rows = local_index(na, my_prow, np_rows, nblk, -1)
      allocate(vmrCPU(max(l_rows,1),na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vmrCPU", 385,  istat,  errorMessage)

      vmrCols = na

      call qr_pdgeqrf_2dcomm_&
           &double&
           &(obj, a_mat(1:lda,1:matrixCols), matrixCols, lda, vmrCPU(1:max(l_rows,1),1:vmrCols), max(l_rows,1), &
                             vmrCols, tauvector(1:na), na, tmat(1:nbw,1:nbw,1), nbw, &
                             nbw, dwork_size(1:1), 1, -1, na, nbw, nblk, nblk, na, na, 1, 0, PQRPARAM(1:11), &
                             mpi_comm_rows, mpi_comm_cols, blockheuristic)

      work_size = int(dwork_size(1))
      allocate(work_blocked(work_size), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: work_blocked", 407,  istat,  errorMessage)
      work_blocked = 0.0_rk
      deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("bandred: vmrCPU", 410,  istat,  errorMessage)

    endif ! which_qr_decomposition

  endif ! useQr

  blk_end = (na-1)/nbw
  if (useGPU .and. .not.(useIntelGPU)) then

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: a_mat", 422,  successGPU)

    cur_l_rows = 0
    cur_l_cols = 0

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("bandred: a_dev", 429,  successGPU)

    successGPU = gpu_malloc(tmat_dev, nbw*nbw*size_of_datatype)
    call check_alloc_GPU_f("bandred: tmat_dev", 432,  successGPU)

    istep = (na-1)/nbw
    blk_end = (na-1)/nbw
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows,1)
    cur_l_cols = max(l_cols,1)
    vmr_size = cur_l_rows*2*n_cols
    umc_size = cur_l_cols*2*n_cols

    istep = (na-1)/nbw - 1
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows2 = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows2,1)
    cur_l_cols = max(l_cols,1)
    vmr_size2 = cur_l_rows*2*n_cols
    umc_size2 = cur_l_cols*2*n_cols

    l_rows = max(l_rows,l_rows2)
    vmr_size = max(vmr_size,vmr_size2)
    umc_size = max(umc_size,umc_size2)

    allocate(vr(l_rows + 1), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      print *,"bandred_&
              &real&
              &: error when allocating vr "//errorMessage
      stop 1
    endif

    successGPU = gpu_malloc_host(vmr_host,vmr_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: vmr_host", 466,  successGPU)
    call c_f_pointer(vmr_host, vmrGPU, (/vmr_size/))

    successGPU = gpu_malloc(vmr_dev, vmr_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: vmr_dev", 470,  successGPU)

    successGPU = gpu_malloc_host(umc_host,umc_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: umc_host", 473,  successGPU)
    call c_f_pointer(umc_host, umcGPU, (/umc_size/))

    successGPU = gpu_malloc(umc_dev, umc_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: umc_dev", 477,  successGPU)

  endif ! useGPU

  !if (useIntelGPU) then
     ! needed later when explict mem copy
  !endif ! useIntelGPU


  do istep = blk_end, 1, -1

    n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

    ! Number of local columns/rows of remaining matrix
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

    ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

    if (useGPU) then
      if (useIntelGPU) then
        ! unify the the name vmr and vmrCPU, as well as vmrGPU
        ! the same for umcCPU and umcGPU
        ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

        allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vmrCPU", 503,  istat,  errorMessage)

        allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: umcCPU", 506,  istat,  errorMessage)

        allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vr", 509,  istat,  errorMessage)
      else
        cur_l_rows = max(l_rows, 1)
        cur_l_cols = max(l_cols, 1)
        vmr_size = cur_l_rows * 2 * n_cols
        umc_size = cur_l_cols * 2 * n_cols
      endif
    else ! GPU not used

      ! unify the the name vmr and vmrCPU, as well as vmrGPU
      ! the same for umcCPU and umcGPU
      ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

      allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vmrCPU", 523,  istat,  errorMessage)

      allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: umcCPU", 526,  istat,  errorMessage)

      allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vr", 529,  istat,  errorMessage)

    endif ! use GPU

    if (useGPU) then
      if (useIntelGPU) then
        vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
      else
        vmrGPU(1 : cur_l_rows * n_cols) = 0.0_rck
        umcGPU(1 : umc_size) = 0.0_rck
      endif
    else
      vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
    endif ! useGPU

    vr(:) = 0.0_rck
    tmat(:,:,istep) = 0.0_rck
    if (useGPU .and. .not.(useIntelGPU)) then
      lc_start = local_index(istep*nbw+1, my_pcol, np_cols, nblk, -1)
      lc_end   = local_index(istep*nbw+n_cols, my_pcol, np_cols, nblk, -1)
      lr_end   = local_index((istep-1)*nbw + n_cols, my_prow, np_rows, nblk, -1)

      if (lc_start .le. 0) lc_start = 1

      do_memcpy = .false.

      ! Note: mod(nbw,nblk) == 0
      do i_blk = 1, nbw/nblk
        blk_off = (i_blk-1) * nblk
        cur_pcol = pcol(istep*nbw+1+blk_off, nblk, np_cols)

        if (my_pcol == cur_pcol) then
          do_memcpy = .true.
        endif
      enddo

      if (do_memcpy) then
        successGPU = gpu_memcpy2d(int(loc(a_mat(1, lc_start)),kind=c_intptr_t), &
                      int((lda*size_of_datatype),kind=c_intptr_t), &
                      (a_dev + int( ( (lc_start-1) * lda*size_of_datatype),kind=c_intptr_t )), &
                      int(lda*size_of_datatype,kind=c_intptr_t), &
                      int(lr_end*size_of_datatype,kind=c_intptr_t), &
                      int((lc_end - lc_start+1),kind=c_intptr_t),int(gpuMemcpyDeviceToHost,kind=c_int))

        call check_memcpy_GPU_f("bandred: a_dev -> a_mat", 573,  successGPU)
      endif
    endif ! useGPU

    !if (useIntelGPU) then
       ! needed later when explict mem copy
    !endif ! useIntelGPU


    ! Reduce current block to lower triangular form
    if (useQR) then
      if (which_qr_decomposition == 1) then
        vmrCols = 2*n_cols
        call qr_pdgeqrf_2dcomm_&
             &double&
             &(obj, a_mat(1:lda,1:matrixCols), lda, matrixCols, vmrCPU(1:max(l_rows,1),1:vmrCols) ,   &
                                max(l_rows,1), vmrCols, tauvector(1:na), na, &
                                 tmat(1:nbw,1:nbw,istep), nbw, nbw, work_blocked(1:work_size), work_size, &
                                 work_size, na, n_cols, nblk, nblk,        &
                                 istep*nbw+n_cols-nbw, istep*nbw+n_cols, 1,&
                                 0, PQRPARAM(1:11), mpi_comm_rows, mpi_comm_cols,&
                                 blockheuristic)
      endif

    else !useQR
      do lc = n_cols, 1, -1

        ncol = istep*nbw + lc ! absolute column number of householder Vector
        nrow = ncol - nbw ! Absolute number of pivot row

        lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
        lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

        tau = 0

        if (nrow == 1) exit ! Nothing to do

        cur_pcol = pcol(ncol, nblk, np_cols) ! Processor column owning current block

        if (my_pcol==cur_pcol) then

          ! Get Vector to be transformed; distribute last element and norm of
          ! remaining elements to all procs in current column

          vr(1:lr) = a_mat(1:lr,lch) ! Vector to be transformed

          if (my_prow==prow(nrow, nblk, np_rows)) then
            aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
            aux1(2) = vr(lr)
          else
            aux1(1) = dot_product(vr(1:lr),vr(1:lr))
            aux1(2) = 0.0_rck
          endif

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL8, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             allreduce_request1, mpierr)

            call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)

            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL8, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             mpierr)

            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif


          vnorm2 = aux2(1)
          vrl    = aux2(2)

          ! Householder transformation
          call hh_transform_&
             &real&
             &_&
             &double &
                         (obj, vrl, vnorm2, xf, tau, wantDebug)
          ! Scale vr and store Householder Vector for back transformation

          vr(1:lr) = vr(1:lr) * xf
          if (my_prow==prow(nrow, nblk, np_rows)) then
            a_mat(1:lr-1,lch) = vr(1:lr-1)
            a_mat(lr,lch) = vrl
            vr(lr) = 1.0_rck
          else
            a_mat(1:lr,lch) = vr(1:lr)
          endif

        endif

        ! Broadcast Householder Vector and tau along columns

        vr(lr+1) = tau
        if (useNonBlockingCollectivesCols) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(vr, int(lr+1,kind=MPI_KIND), MPI_REAL8, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        bcast_request, mpierr)


          call mpi_wait(bcast_request, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_bcast(vr, int(lr+1,kind=MPI_KIND), MPI_REAL8, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
          vmrGPU(cur_l_rows * (lc - 1) + 1 : cur_l_rows * (lc - 1) + lr) = vr(1:lr)
        else
          vmrCPU(1:lr,lc) = vr(1:lr)
        endif
        tau = vr(lr+1)

        tmat(lc,lc,istep) = tau ! Store tau in diagonal of tmat
        ! Transform remaining columns in current block with Householder Vector
        ! Local dot product

        aux1 = 0.0_rck


        nlc = 0 ! number of local columns
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            if (lr>0) aux1(nlc) = dot_product(vr(1:lr),a_mat(1:lr,lcx))
          endif
        enddo

        ! Get global dot products
        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          if (nlc>0) then
            call mpi_iallreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_REAL8, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      allreduce_request3, mpierr)
            call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          if (nlc>0) then
            call mpi_allreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_REAL8, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif
        ! Transform

        nlc = 0
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            a_mat(1:lr,lcx) = a_mat(1:lr,lcx) - tau*aux2(nlc)*vr(1:lr)
          endif
        enddo
      enddo ! lc

      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        ! store column tiles back to GPU
        if (do_memcpy) then
          successGPU = gpu_memcpy2d((a_dev+ &
                        int(((lc_start-1)*lda*size_of_datatype),kind=c_intptr_t)), &
                        int(lda*size_of_datatype,kind=c_intptr_t), int(loc(a_mat(1,lc_start)),kind=c_intptr_t), &
                        int(lda*size_of_datatype,kind=c_intptr_t), &
                        int(lr_end*size_of_datatype,kind=c_intptr_t), &
                        int((lc_end - lc_start+1),kind=c_intptr_t), &
                        int(gpuMemcpyHostToDevice,kind=c_int))
          call check_memcpy_GPU_f("bandred: a_mat -> a_dev", 937,  successGPU)
        endif
      endif

      if (useGPU_reduction_lower_block_to_tridiagonal .and. useIntelGPU) then
        ! store column tiles back to GPU
      endif

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use dsyrk

      vav = 0
      call obj%timer%start("blas")
      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        if (l_rows>0) &
        call DSYRK('U', 'T',            &
                           int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                           vmrGPU, int(cur_l_rows,kind=BLAS_KIND), &
                           ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      else ! useGPU_reduction_to_tridiagonal
        if (l_rows>0) &
        call DSYRK('U', 'T',           &
                            int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, vmrCPU, &
                            int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      call symm_matrix_allreduce_&
         &double &
                         (obj, n_cols,vav, nbw, nbw,mpi_comm_rows, .true.)
         ! Calculate triangular matrix T for block Householder Transformation
      call obj%timer%start("blas")
      do lc=n_cols,1,-1
        tau = tmat(lc,lc,istep)
        if (lc<n_cols) then
          call DTRMV('U', 'T', 'N',&
                              int(n_cols-lc,kind=BLAS_KIND), tmat(lc+1,lc+1,istep), &
                              int(ubound(tmat,dim=1),kind=BLAS_KIND), vav(lc+1,lc), 1_BLAS_KIND)

          tmat(lc,lc+1:n_cols,istep) = -tau * vav(lc+1:n_cols,lc)
        endif
      enddo
      call obj%timer%stop("blas")
    endif !useQR

    if (useGPU .and. useQR ) then
      if (useIntelGPU) then
        ! copy the data for furhter usage
        ! qr worked on *CPU arrarys
        !vmrCUDA(1:cur_l_rows * n_cols) = vmrCPU(1:cur_l_rows,1:n_cols)
      else     
        ! copy the data for furhter usage
        ! qr worked on *CPU arrarys
        !vmrGPU(1:cur_l_rows * n_cols) = vmrCPU(1:cur_l_rows,1:n_cols)
        if (do_memcpy) then
          successGPU = gpu_memcpy2d((a_dev+ &
                        int(((lc_start-1)*lda*size_of_datatype),kind=c_intptr_t)), &
                        int(lda*size_of_datatype,kind=c_intptr_t), int(loc(a_mat(1,lc_start)),kind=c_intptr_t), &
                        int(lda*size_of_datatype,kind=c_intptr_t), &
                        int(lr_end*size_of_datatype,kind=c_intptr_t), &
                        int((lc_end - lc_start+1),kind=c_intptr_t), &
                        int(gpuMemcpyHostToDevice,kind=c_int))
          call check_memcpy_GPU_f("bandred: a_mat -> a_dev", 1022,  successGPU)
        endif
      endif
    endif

    ! Transpose vmr -> vmc (stored in umc, second half)
    if (useGPU) then
      if (useIntelGPU) then
        call elpa_transpose_vectors_&
             &real&
             &_&
             &double &
                                          (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                           umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                           1, istep*nbw, n_cols, nblk, max_threads, .true.)

      else
        call elpa_transpose_vectors_&
             &real&
             &_&
             &double &
                          (obj, vmrGPU(:), cur_l_rows, mpi_comm_rows, &
                           umcGPU(cur_l_cols * n_cols + 1:), cur_l_cols, &
                           mpi_comm_cols, 1, istep*nbw, n_cols, nblk, max_threads, .true.)
      endif
    else ! useGPU
      call elpa_transpose_vectors_&
           &real&
           &_&
           &double &
                                        (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                         umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                         1, istep*nbw, n_cols, nblk, max_threads, .true.)
    endif

    ! Calculate umc = A**T * vmr
    ! Note that the distributed A has to be transposed
    ! Opposed to direct tridiagonalization there is no need to use the cache locality
    ! of the tiles, so we can use strips of the matrix


!#if 0
!    ! original complex implemetation check for performance
!    umcCPU(1:l_cols,1:n_cols) = 0.0_rck
!    vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
!
!    if (l_cols>0 .and. l_rows>0) then
!      do i=0,(istep*nbw-1)/tile_size
!
!        lcs = i*l_cols_tile+1
!        lce = min(l_cols,(i+1)*l_cols_tile)
!        if (lce<lcs) cycle
!
!        lre = min(l_rows,(i+1)*l_rows_tile)
!
!          call obj%timer%start("blas")
!          call DGEMM('C', 'N', lce-lcs+1, n_cols, lre, ONE, a_mat(1,lcs), ubound(a_mat,dim=1), &
!                     vmrCPU, ubound(vmrCPU,dim=1), ONE, umcCPU(lcs,1), ubound(umcCPU,dim=1))
!          call obj%timer%stop("blas")
!
!        if (i==0) cycle
!        lre = min(l_rows,i*l_rows_tile)
!          call obj%timer%start("blas")
!          call DGEMM('N', 'N', lre, n_cols, lce-lcs+1, ONE, a_mat(1,lcs), lda, &
!                     umcCPU(lcs,n_cols+1), ubound(umcCPU,dim=1), ONE, vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1))
!          call obj%timer%stop("blas")
!      enddo
!
!    endif ! (l_cols>0 .and. l_rows>0)
!#endif 

    !Code for Algorithm 4

    ! n_way is actually a branch for the number of OpenMP threads
    n_way = 1

      if (.not. useGPU .or. useIntelGPU) then
        umcCPU(1:l_cols,1:n_cols) = 0.0_rck
        vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
      endif ! useGPU

      if (l_cols>0 .and. l_rows>0) then

        if (useGPU .and. .not.(useIntelGPU)) then
          successGPU = gpu_memset(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                        0, cur_l_rows*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: vmr_dev", 1193,  successGPU)

          successGPU = gpu_memcpy(vmr_dev, int(loc(vmrGPU(1)),kind=c_intptr_t), &
                        cur_l_rows*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: vmrGPU -> vmr_dev", 1197,  successGPU)

          successGPU = gpu_memset(umc_dev, 0, l_cols*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: umc_dev", 1200,  successGPU)

          successGPU = gpu_memcpy(umc_dev+l_cols*n_cols*size_of_datatype, &
                        int(loc(umcGPU(1+l_cols*n_cols)),kind=c_intptr_t), &
                        (umc_size-l_cols*n_cols)*size_of_datatype, &
                        gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev", 1206,  successGPU)
        endif ! useGPU

        do i=0,(istep*nbw-1)/tile_size

          lcs = i*l_cols_tile+1
          lce = min(l_cols,(i+1)*l_cols_tile)
          if (lce<lcs) cycle
          lre = min(l_rows,(i+1)*l_rows_tile)

          if (useGPU) then
            if (useIntelGPU) then
              call obj%timer%start("mkl_offload")
              call obj%timer%stop("mkl_offload")
              if (i==0) cycle
              lre = min(l_rows,i*l_rows_tile)
              call obj%timer%start("mkl_offload")

              if (isSkewsymmetric) then

              else
              endif
              call obj%timer%stop("mkl_offload")

            else
              call obj%timer%start("gpublas")
              call gpublas_DGEMM('T', 'N',                   &
                                         lce-lcs+1, n_cols, lre,     &
                                         ONE, (a_dev + ((lcs-1)*lda* &
                                         size_of_datatype)),         &
                                         lda, vmr_dev,cur_l_rows,    &
                                         ONE, (umc_dev+ (lcs-1)*     &
                                             size_of_datatype),      &
                                         cur_l_cols)

              call obj%timer%stop("gpublas")

              if(i==0) cycle
              call obj%timer%start("gpublas")

              lre = min(l_rows,i*l_rows_tile)
              if (isSkewsymmetric) then
                call gpublas_DGEMM('N', 'N', lre,n_cols, lce-lcs+1, -ONE, &
                              (a_dev+ ((lcs-1)*lda*                 &
                                    size_of_datatype)),             &
                         lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                size_of_datatype),              &
                                cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                              size_of_datatype),              &
                                cur_l_rows)
              else
                call gpublas_DGEMM('N', 'N', lre,n_cols, lce-lcs+1, ONE, &
                                            (a_dev+ ((lcs-1)*lda*                 &
                                                  size_of_datatype)),             &
                                       lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                              size_of_datatype),              &
                                              cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                                            size_of_datatype),              &
                                              cur_l_rows)
              endif
              call obj%timer%stop("gpublas")
            endif
          else ! useGPU

            call obj%timer%start("blas")
            call DGEMM('T', 'N',       &
                                int(lce-lcs+1,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lre,kind=BLAS_KIND), &
                                ONE, a_mat(1,lcs), int(ubound(a_mat,dim=1),kind=BLAS_KIND), &
                                vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ONE, umcCPU(lcs,1), &
                                int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
            call obj%timer%stop("blas")
            if (i==0) cycle
            lre = min(l_rows,i*l_rows_tile)
            call obj%timer%start("blas")

            if (isSkewsymmetric) then
              call DGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  -ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )

            else
              call DGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )
            endif
            call obj%timer%stop("blas")
          endif ! useGPU
        enddo ! i=0,(istep*nbw-1)/tile_size

        if (useGPU .and. .not.(useIntelGPU)) then
          if (tile_size < istep*nbw .or. n_way > 1) then
            successGPU = gpu_memcpy(int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                          vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                          (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyDeviceToHost)
            call check_memcpy_GPU_f("bandred: vmr_dev -> vmrGPU", 1342,  successGPU)
          endif

          successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                        umc_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU", 1347,  successGPU)
        endif ! useGPU
      endif ! l_cols>0 .and. l_rows>0

    ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if ur has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix

    ! Or if we used the Algorithm 4
    if (tile_size < istep*nbw .or. n_way > 1) then

      if (useGPU) then
        if (useIntelGPU) then
          call elpa_reduce_add_vectors_&
          &real&
          &_&
          &double &
                                      (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                       umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                      istep*nbw, n_cols, nblk, max_threads)

        else

          call elpa_reduce_add_vectors_&
               &real&
               &_&
               &double &
                               (obj, vmrGPU(cur_l_rows * n_cols + 1:),cur_l_rows,  &
                                mpi_comm_rows, umcGPU,                            &
                                cur_l_cols, mpi_comm_cols, istep*nbw, n_cols, nblk, max_threads)
        endif
      else ! useGPU

        call elpa_reduce_add_vectors_&
        &real&
        &_&
        &double &
                                         (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                          umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          istep*nbw, n_cols, nblk, max_threads)
      endif ! useGPU
    endif ! tile_size < istep*nbw .or. n_way > 1

    if (l_cols>0) then

      if (useGPU) then
        if (useIntelGPU) then
          allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpCPU", 1399,  istat,  errorMessage)

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL8,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           allreduce_request4, mpierr)
            call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL8,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("bandred: tmpCPU", 1421,  istat,  errorMessage)
        else
          allocate(tmpGPU(l_cols * n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpGPU", 1425,  istat,  errorMessage)
          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")

            call mpi_iallreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL8, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request5, mpierr)
            call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")

            call mpi_allreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL8, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          if (allocated(tmpGPU)) then
            deallocate(tmpGPU, stat=istat, errmsg=errorMessage)
            call check_deallocate_f("bandred: tmpGPU", 1447,  istat,  errorMessage)
          endif
        endif

      else ! useGPU

        allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: tmpCPU", 1454,  istat,  errorMessage)

        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL8,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6, mpierr)
          call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL8,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: tmpCPU", 1474,  istat,  errorMessage)
      endif ! useGPU
    endif ! l_cols > 0

    ! U = U * Tmat**T

    if (useGPU) then
      if (useIntelGPU) then
        call obj%timer%start("mkl_offload")

        call DTRMM('Right', 'Upper', 'T', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

        call DGEMM('T', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

        call DTRMM('Right', 'Upper', 'T', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
        call obj%timer%stop("mkl_offload")

      else
        successGPU = gpu_memcpy(umc_dev, int(loc(umcGPU(1)),kind=c_intptr_t), &
                      l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev ", 1526,  successGPU)

        successGPU = gpu_memcpy(tmat_dev,int(loc(tmat(1,1,istep)),kind=c_intptr_t), &
                      nbw*nbw*size_of_datatype,gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: tmat -> tmat_dev ", 1530,  successGPU)

        call obj%timer%start("gpublas")
        call gpublas_DTRMM('Right', 'Upper', 'T', 'Nonunit',  &
                            l_cols, n_cols, ONE, tmat_dev, nbw, umc_dev, cur_l_cols)
        call obj%timer%stop("gpublas")

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T
        call obj%timer%start("gpublas")
        call gpublas_DGEMM('T', 'N',             &
                                 n_cols, n_cols, l_cols, ONE, umc_dev, cur_l_cols, &
                                 (umc_dev+(cur_l_cols * n_cols )*size_of_datatype),cur_l_cols, &
                                 ZERO, vav_dev, nbw)

        call gpublas_DTRMM('Right', 'Upper', 'T', 'Nonunit',    &
           n_cols, n_cols, ONE, tmat_dev, nbw, vav_dev, nbw)
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(vav),kind=c_intptr_t), &
                    vav_dev, nbw*nbw*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: vav_dev -> vav ", 1550,  successGPU)
      endif
    else ! useGPU

      call obj%timer%start("blas")

      call DTRMM('Right', 'Upper', 'T', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      call DGEMM('T', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      call DTRMM('Right', 'Upper', 'T', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
      call obj%timer%stop("blas")

    endif ! useGPU

!#ifdef 1
    if (isSkewsymmetric) then
      call ssymm_matrix_allreduce_&
      &double &
      (obj, n_cols,vav, nbw, nbw ,mpi_comm_cols, .false.)
    else
!#endif
      call symm_matrix_allreduce_&
      &double &
      (obj, n_cols,vav, nbw, nbw ,mpi_comm_cols, .false.)
!#ifdef 1
    endif
!#endif

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(vav_dev, int(loc(vav),kind=c_intptr_t), &
                       nbw*nbw*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("bandred: vav -> vav_dev ", 1599,  successGPU)
    endif
    !if (useIntelGPU) then
      ! needed later
    !endif


    ! U = U - 0.5 * V * VAV

    if (useGPU) then
      if (useIntelGPU) then
      call obj%timer%start("mkl_offload")
      if (isSkewsymmetric) then
        call DGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                        &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      else
        call DGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            -0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                       &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      endif

      call obj%timer%stop("mkl_offload")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &real&
        &_&
        &double &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &real&
       &_&
       &double &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif

      else
        call obj%timer%start("gpublas")
        if (isSkewsymmetric) then
          call gpublas_DGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                    0.5_rk,                      &
                                    (umc_dev+(cur_l_cols * n_cols )* &
                                    size_of_datatype),   &
                                    cur_l_cols, vav_dev,nbw,        &
                                    ONE, umc_dev, cur_l_cols)
        else
          call gpublas_DGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                   -0.5_rk,                      &
                                   (umc_dev+(cur_l_cols * n_cols )* &
                                   size_of_datatype),   &
                                   cur_l_cols, vav_dev,nbw,        &
                                   ONE, umc_dev, cur_l_cols)
        endif
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                    umc_dev, umc_size*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU ", 1681,  successGPU)

        ! Transpose umc -> umr (stored in vmr, second half)
        if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
             &real&
             &_&
             &double &
                         (obj, umcGPU(:), cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        else
          call elpa_transpose_vectors_&
             &real&
             &_&
             &double &
                         (obj, umcGPU, cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        endif

        successGPU = gpu_memcpy(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                    int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                    (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: vmr -> vmrGPU ", 1705,  successGPU)
      endif

    else ! useGPU
      call obj%timer%start("blas")
      if (isSkewsymmetric) then
        call DGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                        &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      else
        call DGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            -0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                       &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      endif

      call obj%timer%stop("blas")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &real&
        &_&
        &double &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &real&
       &_&
       &double &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif
    endif  ! useGPU

    ! A = A - V*U**T - U*V**T


    do i=0,(istep*nbw-1)/tile_size
      lcs = i*l_cols_tile+1
      lce = min(l_cols,(i+1)*l_cols_tile)
      lre = min(l_rows,(i+1)*l_rows_tile)
      if (lce<lcs .or. lre<1) cycle

      if (useGPU) then
        if (useIntelGPU) then
          call obj%timer%start("mkl_offload")
          call DGEMM('N', 'T', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                              int(2*n_cols,kind=BLAS_KIND), &
                              -ONE, &
                              vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                              int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                              ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
          call obj%timer%stop("mkl_offload")

        else
          call obj%timer%start("gpublas")

          call gpublas_DGEMM('N', 'T',     &
                                     lre, lce-lcs+1, 2*n_cols, -ONE, &
                                     vmr_dev, cur_l_rows, (umc_dev +(lcs-1)*  &
                                     size_of_datatype), &
                                     cur_l_cols, ONE, (a_dev+(lcs-1)*lda* &
                                     size_of_datatype), lda)
          call obj%timer%stop("gpublas")
        endif
      else ! useGPU

        call obj%timer%start("blas")
        call DGEMM('N', 'T', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                            int(2*n_cols,kind=BLAS_KIND), &
                            -ONE, &
                            vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                            int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                            ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    enddo ! i=0,(istep*nbw-1)/tile_size

    if (.not.(useGPU) .or. useIntelGPU) then
      if (allocated(vr)) then
        deallocate(vr, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vr", 1871,  istat,  errorMessage)
      endif

      if (allocated(umcCPU)) then
        deallocate(umcCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: umcCPU", 1876,  istat,  errorMessage)
      endif

      if (allocated(vmrCPU)) then
        deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vmrCPU", 1881,  istat,  errorMessage)
      endif
    endif !useGPU

  enddo ! istep - loop

  if (useGPU .and. .not.(useIntelGPU)) then
    ! copy a_dev to a_mat
    ! we do it here, since a is needed on the host in the following routine
    ! (band to tridi). Previously, a has been kept on the device and then
    ! copied in redist_band (called from tridiag_band). However, it seems to
    ! be easier to do it here.
    successGPU = gpu_memcpy(int(loc(a_mat),kind=c_intptr_t), &
                  int(a_dev,kind=c_intptr_t), &
                  int(lda*matrixCols* size_of_datatype, kind=c_intptr_t), &
                  gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("bandred: a_dev -> a_mat ", 1897,  successGPU)

    successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: a_mat ", 1900,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("bandred: a_dev ", 1903,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("bandred: vav_dev ", 1906,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("bandred: tmat_dev ", 1909,  successGPU)

    successGPU = gpu_host_unregister(int(loc(vav),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: vav", 1912,  successGPU)

    if (associated(umcGPU)) then
      nullify(umcGPU)

      successGPU = gpu_free_host(umc_host)
      call check_host_dealloc_GPU_f("bandred: umc_host ", 1918,  successGPU)

      successGPU = gpu_free(umc_dev)
      call check_dealloc_GPU_f("bandred: umc_dev ", 1921,  successGPU)
    endif

    if (associated(vmrGPU)) then
      nullify(vmrGPU)

      successGPU = gpu_free_host(vmr_host)
      call check_host_dealloc_GPU_f("bandred: vmr_host ", 1928,  successGPU)

      successGPU = gpu_free(vmr_dev)
      call check_dealloc_GPU_f("bandred: vmr_dev ", 1931,  successGPU)
    endif
  endif ! useGPU

  !if (useIntelGPU) then
  !   ! needed later
  !endif

  if (allocated(vr)) then
    deallocate(vr, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vr", 1941,  istat,  errorMessage)
  endif

  if (allocated(umcCPU)) then
    deallocate(umcCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: umcCPU", 1946,  istat,  errorMessage)
  endif

  if (allocated(vmrCPU)) then
    deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vmrCPU", 1951,  istat,  errorMessage)
  endif

  if (useQR) then
    if (which_qr_decomposition == 1) then
      deallocate(work_blocked, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("bandred: work_blocked", 1958,  istat,  errorMessage)

      deallocate(tauvector, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("bandred: tauvector", 1961,  istat,  errorMessage)
    endif
  endif

  call obj%timer%stop("bandred_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

end subroutine bandred_&
&real&
&_&
&double

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







subroutine symm_matrix_allreduce_&
&double &
                    (obj, n, a, lda, ldb, comm, isRows)
!-------------------------------------------------------------------------------
!  symm_matrix_allreduce: Does an mpi_allreduce for a symmetric matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use precision
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)             :: n, lda, ldb, comm
  real(kind=rk8)     :: a(lda,*)
  integer(kind=ik)             :: i, nc
  integer(kind=MPI_KIND)       :: mpierr
  real(kind=rk8)     :: h1(n*n), h2(n*n)
  integer(kind=MPI_KIND)       :: allreduce_request1
  logical                      :: useNonBlockingCollectives
  logical                      :: useNonBlockingCollectivesCols
  logical                      :: useNonBlockingCollectivesRows
  logical, intent(in)          :: isRows
  integer(kind=c_int)          :: non_blocking_collectives_rows, error, &
                                  non_blocking_collectives_cols

  call obj%timer%start("&
          &symm_matrix_allreduce&
          &" // &
          &"_double"&
          )

  call obj%get("nbc_row_sym_allreduce", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa_sym_allreduce. Aborting..."
    stop
  endif

  call obj%get("nbc_col_sym_allreduce", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa_sym_allreduce. Aborting..."
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

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  nc = 0
  do i=1,n
    h1(nc+1:nc+i) = a(1:i,i)
    nc = nc+i
  enddo

  if (useNonBlockingCollectives) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                     int(comm,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                     int(comm,kind=MPI_KIND), mpierr)
    call obj%timer%stop("mpi_communication")
  endif
  nc = 0
  do i=1,n
    a(1:i,i) = h2(nc+1:nc+i)
    a(i,1:i-1) = a(1:i-1,i)
    nc = nc+i
  enddo

! nc = 0
! do i=1,n
!   a(1:i,i) = h2(nc+1:nc+i)
!   a(i,1:i-1) = a(1:i-1,i)
!   nc = nc+i
! enddo

  call obj%timer%stop("&
          &symm_matrix_allreduce&
          &" // &
          &"_double"&
          )

end subroutine symm_matrix_allreduce_&
&double



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







subroutine ssymm_matrix_allreduce_&
&double &
                    (obj, n, a, lda, ldb, comm, isRows)
!-------------------------------------------------------------------------------
!  symm_matrix_allreduce: Does an mpi_allreduce for a symmetric matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use precision
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)             :: n, lda, ldb, comm
  real(kind=rk8)     :: a(lda,*)
  integer(kind=ik)             :: i, nc
  integer(kind=MPI_KIND)       :: mpierr
  real(kind=rk8)     :: h1(n*n), h2(n*n)
  integer(kind=MPI_KIND)       :: allreduce_request1
  logical                      :: useNonBlockingCollectives
  logical                      :: useNonBlockingCollectivesCols
  logical                      :: useNonBlockingCollectivesRows
  logical, intent(in)          :: isRows
  integer(kind=c_int)          :: non_blocking_collectives_rows, error, &
                                  non_blocking_collectives_cols

  call obj%timer%start("&
          &ssymm_matrix_allreduce&
          &" // &
          &"_double"&
          )

  call obj%get("nbc_row_sym_allreduce", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa_sym_allreduce. Aborting..."
    stop
  endif

  call obj%get("nbc_col_sym_allreduce", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa_sym_allreduce. Aborting..."
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

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  nc = 0
  do i=1,n
    h1(nc+1:nc+i) = a(1:i,i)
    nc = nc+i
  enddo

  if (useNonBlockingCollectives) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                     int(comm,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                     int(comm,kind=MPI_KIND), mpierr)
    call obj%timer%stop("mpi_communication")
  endif
  nc = 0
  do i=1,n
    a(1:i,i) = h2(nc+1:nc+i)
    a(i,1:i-1) = - a(1:i-1,i)
    nc = nc+i
  enddo

! nc = 0
! do i=1,n
!   a(1:i,i) = h2(nc+1:nc+i)
!   a(i,1:i-1) = a(1:i-1,i)
!   nc = nc+i
! enddo

  call obj%timer%stop("&
          &ssymm_matrix_allreduce&
          &" // &
          &"_double"&
          )

end subroutine ssymm_matrix_allreduce_&
&double








subroutine trans_ev_band_to_full_&
    &real&
    &_&
    &double &
    (obj, na, nqc, nblk, nbw, a_mat, lda, tmat, q_mat, &
     ldq, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, useGPU &
     ,useQr)

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_real/complex:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a_mat, number of rows of matrix q_mat
!
!  nqc         Number of columns of matrix q_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a_mat(lda,matrixCols)    Matrix containing the Householder vectors (i.e. matrix a_mat after bandred_real/complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat and q_mat
!
!  tmat(nbw,nbw,numBlocks) Factors returned by bandred_real/complex
!
!  q_mat           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q_mat
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa_abstract_impl
  use elpa_blas_interfaces

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
  logical, intent(in)                    :: useGPU
  logical, intent(in)                     :: useQR
  integer(kind=ik)                       :: na, nqc, lda, ldq, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols
  real(kind=rck)                :: a_mat(lda,*)
  real(kind=rck)                :: q_mat(ldq,*), tmat(nbw,nbw,*)

  integer(kind=ik)                       :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                 :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, mpierr
  integer(kind=ik)                       :: max_blocks_row, max_blocks_col, max_local_rows, &
                                            max_local_cols
  integer(kind=ik)                       :: l_cols, l_rows, l_colh, n_cols
  integer(kind=ik)                       :: istep, lc, ncol, nrow, nb, ns

  real(kind=rck), allocatable   :: hvb(:)
  real(kind=rck), pointer       :: hvm(:,:), tmp1(:), tmp2(:)
  ! hvm_dev is fist used and set in this routine
  ! q_mat is changed in trans_ev_tridi on the host, copied to device and passed here. this can be adapted
  ! tmp_dev is first used in this routine
  ! tmat_dev is not passed along from bandred_real
  integer(kind=C_intptr_T)               :: hvm_dev, q_dev, tmp_dev, tmat_dev, t_tmp_dev, dev_offset
  type(c_ptr)                            :: hvm_host, tmp1_host, tmp2_host
  type(c_ptr)                            :: tmat_mpi_dev, t_tmp_mpi_dev
  real(kind=rck), pointer       :: tmat_mpi_deviceptr(:,:)
  real(kind=rck), pointer       :: t_tmp_mpi_deviceptr(:)

  type(c_ptr)                            :: hvm_mpi_dev
  real(kind=rck), pointer       :: hvm_mpi_deviceptr(:,:)


  integer(kind=ik)                       :: i

  real(kind=rck), allocatable, target   :: tmat_complete(:,:), t_tmp(:,:), t_tmp2(:,:)
  integer(kind=ik)                       :: t_cols, t_rows, ii, jj
  integer(kind=ik)                       :: cwy_blocking

  integer(kind=ik)                       :: istat
  character(200)                         :: errorMessage
  character(20)                          :: gpuString
  logical                                :: successGPU
  integer(kind=c_intptr_t), parameter    :: size_of_datatype = size_of_&
                                                               &double&
                                                               &_&
                                                               &real
  integer(kind=ik)                       :: blocking_factor, error, blk_end
  logical                                :: useIntelGPU
  integer(kind=MPI_KIND)                 :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                :: useNonBlockingCollectivesCols
  logical                                :: useNonBlockingCollectivesRows
  integer(kind=c_int)                    :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_band_to_full_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa2_band_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_band_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_band_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_band_to_full. Aborting..."
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

  call obj%get("blocking_in_band_to_full",blocking_factor,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem getting option for blocking_in_band_to_full. Aborting..."
    stop
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  max_blocks_row = ((na -1)/nblk)/np_rows + 1 ! Rows of a_mat
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1 ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  cwy_blocking = blocking_factor * nbw

  if (useGPU) then
    if (useIntelGPU) then
      allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp1", 249,  istat,  errorMessage)
  
      allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp2", 252,  istat,  errorMessage)
  
      allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: hvm", 255,  istat,  errorMessage)

    else
      ! copy q_mat to q_dev
      successGPU = gpu_malloc(q_dev,ldq*matrixCols*size_of_datatype)
      call check_alloc_GPU_f("trans_ev_band_to_full: q_dev", 260,  successGPU)

      successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuHostRegisterDefault)
      call check_host_register_GPU_f("trans_ev_band_to_full: q_mat", 264,  successGPU)

      successGPU = gpu_memcpy(q_dev,int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_mat -> q_dev", 268,  successGPU)

      successGPU = gpu_malloc_host(tmp1_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp1_host", 271,  successGPU)
      call c_f_pointer(tmp1_host, tmp1, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(tmp2_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp2_host", 275,  successGPU)
      call c_f_pointer(tmp2_host, tmp2, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(hvm_host,max_local_rows*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: hvm_host", 279,  successGPU)
      call c_f_pointer(hvm_host, hvm, (/max_local_rows,cwy_blocking/))
    endif
  else ! useGPU
    allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp1", 284,  istat,  errorMessage)

    allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp2", 287,  istat,  errorMessage)

    allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: hvm", 290,  istat,  errorMessage)
  endif !useGPU

  allocate(hvb(max_local_rows*cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: hvb", 294,  istat,  errorMessage)

  allocate(tmat_complete(cwy_blocking,cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: tmat_complete", 297,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_host_register(int(loc(tmat_complete),kind=c_intptr_t), &
                  cwy_blocking * cwy_blocking * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_band_to_full: tmat_complete", 303,  successGPU)
  endif
  !if (useIntelGPU) then
  !  ! needed later
  !endif



  if (blocking_factor > 1) then
    allocate(t_tmp(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp", 322,  istat,  errorMessage)

    allocate(t_tmp2(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp2", 325,  istat,  errorMessage)
  endif

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_malloc(hvm_dev,max_local_rows*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: hvm_dev", 330,  successGPU)

    successGPU = gpu_malloc(tmp_dev,max_local_cols*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmp_dev", 333,  successGPU)

    successGPU = gpu_malloc(tmat_dev,cwy_blocking*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmat_dev", 336,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  hvm = 0.0_rck ! Must be set to 0 !!!
  hvb = 0.0_rck ! Safety only
  tmp1 = 0.0_rck
  tmp2 = 0.0_rck
  tmat_complete = 0.0_rck
  if (blocking_factor > 1) then
     t_tmp = 0.0_rck ! Must be set to 0 !!!
     t_tmp2 = 0.0_rck

  endif
  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  blk_end = ((na-1)/nbw-1)/blocking_factor + 1
  do istep=1, blk_end

    ! This the call when using na >= ((blocking_factor+1)*nbw)
    ! n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw
    ! Number of columns in current step
    ! As an alternative we add some special case handling if na < cwy_blocking
    if (na < cwy_blocking) then
      n_cols = MAX(0, na-nbw)
      if ( n_cols .eq. 0 ) then
        exit
      end if
    else
      n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw ! Number of columns in current step
    end if

    ! Broadcast all Householder vectors for current step compressed in hvb

    nb = 0
    ns = 0

    do lc = 1, n_cols
      ncol = (istep-1)*cwy_blocking + nbw + lc ! absolute column number of householder Vector
      nrow = ncol - nbw ! absolute number of pivot row

      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
      l_colh = local_index(ncol , my_pcol, np_cols, nblk, -1) ! HV local column number

      if (my_pcol==pcol(ncol, nblk, np_cols)) hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)

      nb = nb+l_rows

      if (lc==n_cols .or. mod(ncol,nblk)==0) then
        if (useNonBlockingCollectivesCols) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_REAL8,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_REAL8,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         mpierr)
          call obj%timer%stop("mpi_communication")
        endif
        ns = nb
      endif
    enddo ! lc

    ! Expand compressed Householder vectors into matrix hvm

    nb = 0
    do lc = 1, n_cols
      nrow = (istep-1)*cwy_blocking + lc ! absolute number of pivot row
      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

      hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
      if (my_prow==prow(nrow, nblk, np_rows)) hvm(l_rows+1,lc) = 1.0_rck
      nb = nb+l_rows
    enddo

    l_rows = local_index(MIN(na,(istep+1)*cwy_blocking), my_prow, np_rows, nblk, -1)

    ! compute tmat2 out of tmat(:,:,)
    tmat_complete = 0
    do i = 1, blocking_factor
      t_cols = MIN(nbw, n_cols - (i-1)*nbw)
      if (t_cols <= 0) exit
      t_rows = (i - 1) * nbw
      tmat_complete(t_rows+1:t_rows+t_cols,t_rows+1:t_rows+t_cols) = tmat(1:t_cols,1:t_cols,(istep-1)*blocking_factor + i)

      if (i > 1) then
        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
          else 
          ! replace this with a cugemm if beneficial
          call obj%timer%start("blas")
          call DGEMM('T', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
          endif
        else ! useGPU
          call obj%timer%start("blas")
          call DGEMM('T', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU


        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_REAL8, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_REAL8, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif

        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
          else
            ! could use the cugemm here if beneficial
            call obj%timer%start("blas")
            call DTRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call DTRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                              tmat_complete(t_rows+1,t_rows+1), &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
            
          endif

        else ! useGPU
          call obj%timer%start("blas")
          call DTRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call DTRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                            tmat_complete(t_rows+1,t_rows+1), &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call obj%timer%stop("blas")
          tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
        endif ! useGPU


      endif
    enddo

    ! Q = Q - V * T**T * V**T * Q

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(hvm_dev, int(loc(hvm),kind=c_intptr_t), &
                          max_local_rows*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: hvm -> hvm_dev", 713,  successGPU)
          call obj%timer%start("gpublas")
          call gpublas_DGEMM('T', 'N', &
                                       n_cols, l_cols, l_rows, ONE, hvm_dev, max_local_rows, &
                                       q_dev, ldq , ZERO, tmp_dev, n_cols)
          call obj%timer%stop("gpublas")

          ! copy data from device to host for a later MPI_ALLREDUCE
          successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                        tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 726,  successGPU)
        endif
      else ! useGPU
        call obj%timer%start("blas")
        call DGEMM('T', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                            hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, &
                           int(n_cols,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    else ! l_rows>0
        tmp1(1:l_cols*n_cols) = 0.0_rck
    endif ! l_rows>0

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                   tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 764,  successGPU)

    endif

    if (useNonBlockingCollectivesRows) then
      call obj%timer%start("mpi_nbc_communication")
      call mpi_iallreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
      call obj%timer%stop("mpi_nbc_communication")
      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
    else
      call obj%timer%start("mpi_communication")
      call mpi_allreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_REAL8, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
    endif

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          call DTRMM('L', 'U', 'T', 'N', &
                              int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
          call DGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                              int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                              q_mat, int(ldq,kind=BLAS_KIND))
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2),kind=c_intptr_t), &
                        l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp2 -> tmp_dev", 819,  successGPU)

          ! needed: as long as not device to device copy
          successGPU = gpu_memcpy(tmat_dev, int(loc(tmat_complete),kind=c_intptr_t), &
                        cwy_blocking*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmat_complete -> tmat_dev", 825,  successGPU)

          call obj%timer%start("gpublas")
          call gpublas_DTRMM('L', 'U', 'T', 'N', &
                                   n_cols, l_cols, ONE, tmat_dev, cwy_blocking, tmp_dev, n_cols)
          call gpublas_DGEMM('N', 'N', l_rows, l_cols, n_cols, -ONE, hvm_dev, max_local_rows, tmp_dev, &
                                     n_cols, ONE, q_dev, ldq)
          call obj%timer%stop("gpublas")
        endif
      else
        call obj%timer%start("blas")
        call DTRMM('L', 'U', 'T', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
        call DGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                            int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                            int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                            q_mat, int(ldq,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU

    endif

  enddo ! istep

  deallocate(hvb, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: hvb", 902,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      deallocate(tmp1, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp1", 907,  istat,  errorMessage)

      deallocate(tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp2", 910,  istat,  errorMessage)

      deallocate(hvm, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: hvm", 913,  istat,  errorMessage)

    else
      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: hvm_dev", 917,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmp_dev", 920,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmat_dev", 923,  successGPU)

      ! final transfer of q_dev
      successGPU = gpu_memcpy(int(loc(q_mat),kind=c_intptr_t), q_dev, ldq*matrixCols*size_of_datatype, &
                    gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_dev -> q_mat", 928,  successGPU)

      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: q_dev", 931,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: q_mat", 934,  successGPU)
      nullify(tmp1)
      nullify(tmp2)
      nullify(hvm)

      ! take care of new pointers nullify them
      nullify(tmat_mpi_deviceptr)
      nullify(hvm_mpi_deviceptr)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp1_host", 948,  successGPU)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp2_host", 951,  successGPU)

      successGPU = gpu_free_host(hvm_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: hvm_host", 954,  successGPU)

      successGPU = gpu_host_unregister(int(loc(tmat_complete),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: tmat_complete", 957,  successGPU)
    endif
  else ! useGPU
    deallocate(tmp1, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp1", 961,  istat,  errorMessage)

    deallocate(tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp2", 964,  istat,  errorMessage)

    deallocate(hvm, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: hvm", 967,  istat,  errorMessage)
  endif ! useGPU

  deallocate(tmat_complete, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: tmat_complete", 971,  istat,  errorMessage)


  if (blocking_factor > 1) then
    deallocate(t_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp", 984,  istat,  errorMessage)

    deallocate(t_tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp2", 987,  istat,  errorMessage)
  endif

  call obj%timer%stop("trans_ev_band_to_full_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

end subroutine trans_ev_band_to_full_&
&real&
    &_&
    &double






subroutine tridiag_band_&
  &real&
  &_&
  &double &
  (obj, na, nb, nblk, a_mat, lda, d, e, matrixCols, &
  hh_trans, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, useGPU, wantDebug, nrThreads, isSkewsymmetric)
  !-------------------------------------------------------------------------------
  ! tridiag_band_real/complex:
  ! Reduces a real symmetric band matrix to tridiagonal form
  !
  !  na          Order of matrix a
  !
  !  nb          Semi bandwith
  !
  !  nblk        blocksize of cyclic distribution, must be the same in both directions!
  !
  !  a_mat(lda,matrixCols)    Distributed system matrix reduced to banded form in the upper diagonal
  !
  !  lda         Leading dimension of a
  !  matrixCols  local columns of matrix a
  !
  ! hh_trans : housholder vectors
  !
  !  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  mpi_comm_rows
  !  mpi_comm_cols
  !              MPI-Communicators for rows/columns
  !  mpi_comm_all
  !              MPI-Communicator for the total processor set
  !-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use precision
  use, intrinsic :: iso_c_binding
  use redist
  use elpa_blas_interfaces
  use elpa_skewsymmetric_blas
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

  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU, wantDebug
  logical, intent(in)                          :: isSkewsymmetric
  integer(kind=ik), intent(in)                 :: na, nb, nblk, lda, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  real(kind=rck), intent(in)         :: a_mat(lda,*)
  real(kind=rk), intent(out)        :: d(na), e(na) ! set only on PE 0
  real(kind=rck), intent(out), allocatable   :: hh_trans(:,:)

  real(kind=rk)                     :: vnorm2
  real(kind=rck)                     :: hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
  real(kind=rck)                     :: hd(nb), hs(nb)

  integer(kind=ik)                             :: i, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
  integer(kind=ik)                             :: my_pe, n_pes
  integer(kind=ik)                             :: my_prow, np_rows, my_pcol, np_cols
  integer(kind=MPI_KIND)                       :: my_peMPI, n_pesMPI, mpierr
  integer(kind=MPI_KIND)                       :: my_prowMPI, np_rowsMPI, my_pcolMPI, np_colsMPI
  integer(kind=MPI_KIND)                       :: ireq_ab, ireq_hv
  integer(kind=ik)                             :: na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off
  integer(kind=ik), intent(in)                 :: nrThreads
  integer(kind=ik), allocatable                :: global_id(:,:), hh_cnt(:), hh_dst(:)
  integer(kind=MPI_KIND), allocatable          :: ireq_hhr(:), ireq_hhs(:)
  integer(kind=ik), allocatable                :: limits(:), snd_limits(:,:)
  integer(kind=ik), allocatable                :: block_limits(:)
  real(kind=rck), allocatable         :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)
  integer                                      :: istat
  integer(kind=ik)                             :: nblockEnd
  character(200)                               :: errorMessage
  character(20)                                :: gpuString

  logical                                      :: useIntelGPU

   integer(kind=MPI_KIND)                      :: allreduce_request1, allreduce_request2
   logical                                     :: useNonBlockingCollectivesAll
   integer(kind=c_int)                         :: non_blocking_collectives, error

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_band_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_all_elpa2_band_to_tridi", non_blocking_collectives, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives in elpa2_band_to_tridi. Aborting..."
    stop
  endif

  if (non_blocking_collectives .eq. 1) then
    useNonBlockingCollectivesAll = .true.
  else
    useNonBlockingCollectivesAll = .false.
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND) ,my_peMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_all,kind=MPI_KIND) ,n_pesMPI ,mpierr)

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND),my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND),np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND),my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND),np_colsMPI ,mpierr)

  my_pe = int(my_peMPI,kind=MPI_KIND)
  n_pes = int(n_pesMPI,kind=MPI_KIND)
  my_prow = int(my_prowMPI,kind=MPI_KIND)
  np_rows = int(np_rowsMPI,kind=MPI_KIND)
  my_pcol = int(my_pcolMPI,kind=MPI_KIND)
  np_cols = int(np_colsMPI,kind=MPI_KIND)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Get global_id mapping 2D procssor coordinates to global id

  allocate(global_id(0:np_rows-1,0:np_cols-1), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: global_id", 202,  istat,  errorMessage)

  global_id(:,:) = 0
  global_id(my_prow, my_pcol) = my_pe


  if (wantDebug) call obj%timer%start("mpi_communication")
  if (useNonBlockingCollectivesAll) then
    call mpi_iallreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
  else
    call mpi_allreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), mpierr)
  endif
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Total number of blocks in the band:

  nblocks_total = (na-1)/nb + 1

  ! Set work distribution

  allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: block_limits", 246,  istat,  errorMessage)

  call divide_band(obj,nblocks_total, n_pes, block_limits)

  ! nblocks: the number of blocks for my task
  nblocks = block_limits(my_pe+1) - block_limits(my_pe)

  ! allocate the part of the band matrix which is needed by this PE
  ! The size is 1 block larger than needed to avoid extensive shifts
  allocate(ab(2*nb,(nblocks+1)*nb), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: ab", 256,  istat,  errorMessage)

  ab = 0.0_rck ! needed for lower half, the extra block should also be set to 0 for safety

  ! n_off: Offset of ab within band
  n_off = block_limits(my_pe)*nb

  ! Redistribute band in a to ab
  call redist_band_&
  &real&
  &_&
  &double&
  &(obj,a_mat, lda, na, nblk, nb, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, ab)

  ! Calculate the workload for each sweep in the back transformation
  ! and the space requirements to hold the HH vectors

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: limits", 274,  istat,  errorMessage)

  call determine_workload(obj,na, nb, np_rows, limits)
  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  do n = 1, nblocks_total
    call determine_workload(obj, nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    ! add to number of householder vectors
    ! please note: for nx==1 the one and only HH Vector is 0 and is neither calculated nor send below!
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_hh_vecs = num_hh_vecs + local_size
      num_chunks  = num_chunks+1
    endif
    nx = nx - nb
  enddo

  ! Allocate space for HH vectors

  allocate(hh_trans(nb,num_hh_vecs), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_trans", 297,  istat,  errorMessage)

  ! Allocate and init MPI requests

  allocate(ireq_hhr(num_chunks), stat=istat, errmsg=errorMessage) ! Recv requests
  call check_allocate_f("tridiag_band: ireq_hhr", 302,  istat,  errorMessage)
  allocate(ireq_hhs(nblocks), stat=istat, errmsg=errorMessage)    ! Send requests
  call check_allocate_f("tridiag_band: ireq_hhs", 304,  istat,  errorMessage)

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  nt = 0
  nBlockEnd=1
  do n = 1, nblocks_total
    call determine_workload(obj,nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_chunks  = num_chunks+1
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_irecv(hh_trans(1,num_hh_vecs+1), int(nb*local_size,kind=MPI_KIND),  MPI_REAL8,     &
                    int(nt,kind=MPI_KIND), int(10+n-block_limits(nt),kind=MPI_KIND), &
                    int(mpi_comm_all,kind=MPI_KIND), ireq_hhr(num_chunks), mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")

      num_hh_vecs = num_hh_vecs + local_size
    endif
    nx = nx - nb
    if (n == block_limits(nt+1)) then
      nt = nt + 1
    endif
  enddo
  ireq_hhs(:) = MPI_REQUEST_NULL
  ! Buffers for gathering/sending the HH vectors

  allocate(hh_gath(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! gathers HH vectors
  call check_allocate_f("tridiag_band: hh_gath", 341,  istat,  errorMessage)

  allocate(hh_send(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! send buffer for HH vectors
  call check_allocate_f("tridiag_band: hh_send", 344,  istat,  errorMessage)

  hh_gath(:,:,:) = 0.0_rck
  hh_send(:,:,:) = 0.0_rck

  ! Some counters

  allocate(hh_cnt(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_cnt", 352,  istat,  errorMessage)

  allocate(hh_dst(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_dst", 355,  istat,  errorMessage)

  hh_cnt(:) = 1 ! The first transfomation Vector is always 0 and not calculated at all
  hh_dst(:) = 0 ! PE number for receive
  ireq_ab = MPI_REQUEST_NULL
  ireq_hv = MPI_REQUEST_NULL
  ! Limits for sending

  allocate(snd_limits(0:np_rows,nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: snd_limits", 366,  istat,  errorMessage)

  do iblk=1,nblocks
    call determine_workload(obj, na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
  enddo


  ! ---------------------------------------------------------------------------
  ! Start of calculations

  na_s = block_limits(my_pe)*nb + 1

  if (my_pe>0 .and. na_s<=na) then
    ! send first column to previous PE
    ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
    ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
    if (wantDebug) call obj%timer%start("mpi_communication")
    call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_REAL8, &
                   int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), ireq_ab, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif


   do istep=1,na-nblockEnd

     if (my_pe==0) then
       n = MIN(na-na_s,nb) ! number of rows to be reduced
       hv(:) = 0.0_rck
       hd(:) = 0.0_rck
       tau = 0.0_rck

       ! Transform first column of remaining matrix
       ! The last step (istep=na-1) is only needed for sending the last HH vectors.
       ! We don't want the sign of the last element flipped (analogous to the other sweeps)

       if (istep < na-1) then
         ! Transform first column of remaining matrix
         vnorm2 = sum(ab(3:n+1,na_s-n_off)**2)

          call hh_transform_&
               &real&
               &_&
               &double &
                          (obj, ab(2,na_s-n_off), vnorm2, hf, tau, wantDebug)

          hv(1) = 1.0_rck
          hv(2:n) = ab(3:n+1,na_s-n_off)*hf
       endif

       if (isSkewsymmetric) then
         d(istep) = 0.0_rk
       else
         d(istep) = ab(1,na_s-n_off)
       endif
       e(istep) = ab(2,na_s-n_off)

       if (istep == na-1) then
         if (isSkewsymmetric) then
           d(na) = 0
         else
           d(na) = ab(1,na_s+1-n_off)
         endif

         e(na) = 0.0_rck
       endif
     else
       if (na>na_s) then
         ! Receive Householder Vector from previous task, from PE owning subdiagonal


         if (wantDebug) call obj%timer%start("mpi_communication")

         call mpi_recv(hv, int(nb,kind=MPI_KIND), MPI_REAL8, &
                       int(my_pe-1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")


         tau = hv(1)
         hv(1) = 1.0_rck
       endif
      endif

      na_s = na_s+1
      if (na_s-n_off > nb) then
        ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
        ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rck
        n_off = n_off + nb
      endif


      do iblk=1,nblocks
        ns = na_s + (iblk-1)*nb - n_off ! first column in block
        ne = ns+nb-1                    ! last column in block

        if (ns+n_off>na) exit

        ! Store Householder Vector for back transformation

        hh_cnt(iblk) = hh_cnt(iblk) + 1

        hh_gath(1   ,hh_cnt(iblk),iblk) = tau
        hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

        if (hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
          ! Wait for last transfer to finish
          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
          ! Copy vectors into send buffer
          hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
          ! Send to destination

          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_isend(hh_send(1,1,iblk), int(nb*hh_cnt(iblk),kind=MPI_KIND), MPI_REAL8, &
                         global_id(hh_dst(iblk), mod(iblk+block_limits(my_pe)-1,np_cols)), &
                         int(10+iblk,kind=MPI_KIND), int(mpi_comm_all,kind=MPI_KIND), ireq_hhs(iblk), mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! Reset counter and increase destination row
          hh_cnt(iblk) = 0
          hh_dst(iblk) = hh_dst(iblk)+1
        endif

        ! The following code is structured in a way to keep waiting times for
        ! other PEs at a minimum, especially if there is only one block.
        ! For this reason, it requests the last column as late as possible
        ! and sends the Householder Vector and the first column as early
        ! as possible.
        nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
        nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                      ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

        ! Multiply diagonal block and subdiagonal block with Householder Vector

        if (iblk==nblocks .and. nc==nb) then

          ! We need the last column from the next PE.
          ! First do the matrix multiplications without last column ...

          ! Diagonal block, the contribution of the last element is added below!
          ab(1,ne) = 0.0_rck

          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")

            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_dssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call DSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif
            ! Subdiagonal block
            if (nr>0) call DGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")

            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_dssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call DSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif
            ! Subdiagonal block
            if (nr>0) call DGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif

          ! ... then request last column ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_recv(ab(1,ne), int(nb+1,kind=MPI_KIND), MPI_REAL8,  &
                       int(my_pe+1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and complete the result
          hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
          hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

        else

          ! Normal matrix multiply
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_dssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call DSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                                hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif

            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_dssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call DSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                                hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif
            if (nr>0) call DGEMV('N', int(nr,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), tau, ab(nb+1,ns), &
                                        int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Calculate first column of subdiagonal block and calculate new
        ! Householder transformation for this column
        hv_new(:) = 0.0_rck ! Needed, last rows must be 0 for nr < nb
        tau_new = 0.0_rck
        if (nr>0) then

          ! complete (old) Householder transformation for first column

          ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

          ! calculate new Householder transformation ...
          if (nr>1) then
            vnorm2 = sum(ab(nb+2:nb+nr,ns)**2)

            call hh_transform_&
                &real&
                &_&
                &double &
               (obj, ab(nb+1,ns), vnorm2, hf, tau_new, wantDebug)
            hv_new(1) = 1.0_rck
            hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
            ab(nb+2:,ns) = 0.0_rck
          endif ! nr > 1

          ! ... and send it away immediatly if this is the last block

          if (iblk==nblocks) then
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            hv_s(1) = tau_new
            hv_s(2:) = hv_new(2:)

            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_isend(hv_s, int(nb,kind=MPI_KIND), MPI_REAL8, &
                           int(my_pe+1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                           ireq_hv, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

          endif

        endif

        ! Transform diagonal block
        if (.NOT. isSkewsymmetric) then
          x = dot_product(hv(1:nc),hd(1:nc))*tau
        endif

        if (.NOT. isSkewsymmetric) then
          hd(1:nc) = hd(1:nc) - 0.5_rk*x*hv(1:nc)
        endif
        if (my_pe>0 .and. iblk==1) then

          ! The first column of the diagonal block has to be send to the previous PE
          ! Calculate first column only ...
          if (isSkewsymmetric) then
            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*hv(1) + hv(1:nc)*hd(1)
          else
            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*hv(1) - hv(1:nc)*hd(1)
          endif
          ! ... send it away ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_wait(ireq_ab, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ab_s(1:nb+1) = ab(1:nb+1,ns)

          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_REAL8, &
                         int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                         ireq_ab, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and calculate remaining columns with rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (isSkewsymmetric) then
              if (nc>1) call elpa_dssr2(int(nc-1,kind=BLAS_KIND), hd(2), hv(2), ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            else
              if (nc>1) call DSYR2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                         hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            if (isSkewsymmetric) then 
              if (nc>1) call elpa_dssr2(int(nc-1,kind=BLAS_KIND), hd(2), hv(2), ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            else
              if (nc>1) call DSYR2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                         hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("blas")
          endif

        else
          ! No need to  send, just a rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (isSkewsymmetric) then
              call elpa_dssr2(int(nc,kind=BLAS_KIND), hd, hv, ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            else
              call DSYR2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND,  &
                                  hv, 1_BLAS_KIND, ab(1,ns), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")
            if (isSkewsymmetric) then 
              call elpa_dssr2(int(nc,kind=BLAS_KIND), hd, hv, ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            else
              call DSYR2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND,  &
                                hv, 1_BLAS_KIND, ab(1,ns), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

        if (nr>0) then
          if (nr>1) then
            if (useIntelGPU) then
              if (wantDebug) call obj%timer%start("mkl_offload")
              call DGEMV('T', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("mkl_offload")
            else
              if (wantDebug) call obj%timer%start("blas")
              call DGEMV('T', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("blas")
            endif
            x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
            h(2:nb) = h(2:nb) - x*hv(2:nb)
            ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update
            do i=2,nb
              ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*h(i) - hs(1:nr)*hv(i)
            enddo
          else
            ! No double Householder transformation for nr=1, just complete the row
            do i=2,nb
              ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*hv(i)
            enddo
          endif
        endif

        ! Use new HH Vector for the next block
        hv(:) = hv_new(:)
        tau = tau_new

      enddo


  enddo ! istep

  ! Finish the last outstanding requests


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
  call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

  call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
  call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_barrier(int(mpi_comm_all,kind=MPI_KIND),mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")
  deallocate(ab, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ab", 1363,  istat,  errorMessage)

  deallocate(ireq_hhr, ireq_hhs, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ireq_hhr", 1366,  istat,  errorMessage)

  deallocate(hh_cnt, hh_dst, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_dst", 1369,  istat,  errorMessage)

  deallocate(hh_gath, hh_send, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_gath", 1372,  istat,  errorMessage)

  deallocate(limits, snd_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: limits", 1375,  istat,  errorMessage)

  deallocate(block_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: block_limits", 1378,  istat,  errorMessage)

  deallocate(global_id, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: global_id", 1381,  istat,  errorMessage)

  call obj%timer%stop("tridiag_band_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

  ! intel compiler bug makes these ifdefs necessary
end subroutine tridiag_band_real_&
&double

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






subroutine trans_ev_tridi_to_band_&
&real&
&_&
&double &
(obj, na, nev, nblk, nbw, q, ldq, matrixCols,         &
 hh_trans, mpi_comm_rows, mpi_comm_cols, wantDebug, useGPU, max_threads, success, &
 kernel)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_real/complex:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!  matrixCols  local columns of matrix q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use pack_unpack_cpu
  use pack_unpack_gpu
  use compute_hh_trafo
  use elpa_gpu
  use precision
  use, intrinsic :: iso_c_binding
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

  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU

  integer(kind=ik), intent(in)                 :: kernel
  integer(kind=ik), intent(in)                 :: na, nev, nblk, nbw, ldq, matrixCols, mpi_comm_rows, mpi_comm_cols

  real(kind=rck), target              :: q(ldq,*)

  real(kind=rck), intent(in),target   :: hh_trans(:,:)
  integer(kind=c_intptr_t)                     :: hh_trans_dev
  type(c_ptr)                                  :: hh_trans_mpi_dev
  real(kind=rck), pointer             :: hh_trans_mpi_fortran_ptr(:,:)


  integer(kind=ik)                             :: np_rows, my_prow, np_cols, my_pcol
  integer(kind=MPI_KIND)                       :: np_rowsMPI, my_prowMPI, np_colsMPI, my_pcolMPI
  integer(kind=ik)                             :: i, j, ip, sweep, nbuf, l_nev, a_dim2
  integer(kind=ik)                             :: current_n, current_local_n, current_n_start, current_n_end
  integer(kind=ik)                             :: next_n, next_local_n, next_n_start, next_n_end
  integer(kind=ik)                             :: bottom_msg_length, top_msg_length, next_top_msg_length
  integer(kind=ik)                             :: stripe_width, last_stripe_width, stripe_count
  integer(kind=ik)                             :: num_result_blocks, num_result_buffers, num_bufs_recvd
  integer(kind=ik)                             :: a_off, current_tv_off, max_blk_size
  integer(kind=ik)                             :: src, src_offset, dst, offset, nfact, num_blk
  integer(kind=MPI_KIND)                       :: mpierr

  logical                                      :: flag
  real(kind=rck), pointer             :: aIntern(:,:,:)
  real(kind=rck)                      :: a_var

  type(c_ptr)                                  :: aIntern_ptr

  real(kind=rck), allocatable, target :: row(:)
  
  integer(kind=c_intptr_t)                     :: row_dev
  type(c_ptr)                                  :: row_mpi_dev
  real(kind=rck), pointer             :: row_mpi_fortran_ptr(:)

  real(kind=rck), pointer             :: row_group(:,:)

  real(kind=rck), allocatable, target :: top_border_send_buffer(:,:)
  real(kind=rck), allocatable, target :: top_border_recv_buffer(:,:)
  real(kind=rck), allocatable, target :: bottom_border_send_buffer(:,:)
  real(kind=rck), allocatable, target :: bottom_border_recv_buffer(:,:)

  integer(kind=c_intptr_t)                     :: top_border_recv_buffer_dev, top_border_send_buffer_dev
  type(c_ptr)                                  :: top_border_recv_buffer_mpi_dev, top_border_send_buffer_mpi_dev 
  real(kind=rck), pointer             :: top_border_recv_buffer_mpi_fortran_ptr(:,:), &
                                                  top_border_send_buffer_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                     :: bottom_border_send_buffer_dev, bottom_border_recv_buffer_dev
  type(c_ptr)                                  :: bottom_border_send_buffer_mpi_dev, bottom_border_recv_buffer_mpi_dev
  real(kind=rck), pointer             :: bottom_border_send_buffer_mpi_fortran_ptr(:,:), &
                                                  bottom_border_recv_buffer_mpi_fortran_ptr(:,:)
  type(c_ptr)                                  :: aIntern_mpi_dev
  real(kind=rck), pointer             :: aIntern_mpi_fortran_ptr(:,:,:)

  integer(kind=c_intptr_t)                     :: aIntern_dev
  integer(kind=c_intptr_t)                     :: bcast_buffer_dev

  type(c_ptr)                                  :: bcast_buffer_mpi_dev
  real(kind=rck), pointer             :: bcast_buffer_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: num
  integer(kind=c_intptr_t)                     :: dev_offset, dev_offset_1
  integer(kind=c_intptr_t)                     :: row_group_dev

  type(c_ptr)                                  :: row_group_mpi_dev
  real(kind=rck), pointer             :: row_group_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: q_dev
  type(c_ptr)                                  :: q_mpi_dev
  real(kind=rck), pointer             :: q_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: hh_tau_dev
  integer(kind=ik)                             :: row_group_size, unpack_idx

  type(c_ptr)                                  :: row_group_host, bcast_buffer_host

  integer(kind=ik)                             :: n_times
  integer(kind=ik)                             :: chunk, this_chunk

  real(kind=rck), allocatable,target  :: result_buffer(:,:,:)
  integer(kind=c_intptr_t)                     :: result_buffer_dev

  type(c_ptr)                                  :: result_buffer_mpi_dev
  real(kind=rck), pointer             :: result_buffer_mpi_fortran_ptr(:,:,:)

  real(kind=rck), pointer             :: bcast_buffer(:,:)

  integer(kind=ik)                             :: n_off

  integer(kind=MPI_KIND), allocatable          :: result_send_request(:), result_recv_request(:)
  integer(kind=ik), allocatable                :: limits(:)
  integer(kind=MPI_KIND), allocatable          :: top_send_request(:), bottom_send_request(:)
  integer(kind=MPI_KIND), allocatable          :: top_recv_request(:), bottom_recv_request(:)

  ! MPI send/recv tags, arbitrary

  integer(kind=ik), parameter                  :: bottom_recv_tag = 111
  integer(kind=ik), parameter                  :: top_recv_tag    = 222
  integer(kind=ik), parameter                  :: result_recv_tag = 333

  integer(kind=ik), intent(in)                 :: max_threads



  ! Just for measuring the kernel performance
  real(kind=c_double)                          :: kernel_time, kernel_time_recv ! MPI_WTIME always needs double
  ! long integer
  integer(kind=lik)                            :: kernel_flops, kernel_flops_recv

  logical, intent(in)                          :: wantDebug
  logical                                      :: success
  integer(kind=ik)                             :: istat, print_flops
  character(200)                               :: errorMessage
  character(20)                                :: gpuString
  logical                                      :: successGPU
  integer(kind=ik)                             :: error
  integer(kind=c_intptr_t), parameter          :: size_of_datatype = size_of_&
                                                                 &double&
                                                                 &_&
                                                                 &real
  logical, parameter                           :: allComputeOnGPU = .false.

  integer(kind=MPI_KIND)                     :: bcast_request1, allreduce_request1, allreduce_request2, &
                                                allreduce_request3, allreduce_request4
  logical                                    :: useNonBlockingCollectivesCols
  logical                                    :: useNonBlockingCollectivesRows
  integer(kind=c_int)                        :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_tridi_to_band_&
  &real&
  &" // &
  &"_double" //&
  gpuString)

  call obj%get("nbc_row_elpa2_tridi_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_tridi_to_band. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_tridi_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_tridi_to_band. Aborting..."
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

  n_times = 0
  if (useGPU) then
    unpack_idx = 0
    row_group_size = 0
  endif

  success = .true.
  kernel_time = 0.0
  kernel_flops = 0

  if (wantDebug) call obj%timer%start("mpi_communication")
  call MPI_Comm_rank(int(mpi_comm_rows,kind=MPI_KIND) , my_prowMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_rows,kind=MPI_KIND) , np_rowsMPI , mpierr)
  call MPI_Comm_rank(int(mpi_comm_cols,kind=MPI_KIND) , my_pcolMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_cols,kind=MPI_KIND) , np_colsMPI , mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")

  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &real&
                            &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &real&
                            &: band backtransform works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  nfact = nbw / nblk


  ! local number of eigenvectors
  l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

  if (l_nev==0) then
    stripe_width = 0
    stripe_count = 0
    last_stripe_width = 0

  else ! l_nev


    ! Suggested stripe width is 48 since 48*64 real*8 numbers should fit into
    ! every primary cache
    ! Suggested stripe width is 48 - should this be reduced for the complex case ???

    if (useGPU) then
      stripe_width = 1024 ! Must be a multiple of 4
      stripe_count = (l_nev - 1) / stripe_width + 1

    else ! useGPU
      call obj%get("stripewidth_real",stripe_width, error)

      !stripe_width = 48 ! Must be a multiple of 4


      stripe_count = (l_nev-1)/stripe_width + 1

      ! Adapt stripe width so that last one doesn't get too small

      stripe_width = (l_nev-1)/stripe_count + 1

      if (kernel .eq. ELPA_2STAGE_REAL_AVX512_BLOCK2 .or. &
          kernel .eq. ELPA_2STAGE_REAL_AVX512_BLOCK4 .or. &
          kernel .eq. ELPA_2STAGE_REAL_AVX512_BLOCK6 .or. &
          kernel .eq. ELPA_2STAGE_REAL_SVE512_BLOCK2 .or. &
          kernel .eq. ELPA_2STAGE_REAL_SVE512_BLOCK4 .or. &
          kernel .eq. ELPA_2STAGE_REAL_SVE512_BLOCK6  &
          ) then

        stripe_width = ((stripe_width+7)/8)*8 ! Must be a multiple of 8 because of AVX-512 memory alignment of 64 bytes
                                              ! (8 * sizeof(double) == 64)

      else
        stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 4 because of AVX/SSE memory alignment of 32 bytes
                                            ! (4 * sizeof(double) == 32)
      endif

   endif ! useGPU

   last_stripe_width = l_nev - (stripe_count-1)*stripe_width

  endif ! l_nev

  ! Determine the matrix distribution at the beginning

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: limits", 604,  istat,  errorMessage)
  call determine_workload(obj,na, nbw, np_rows, limits)

  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  a_dim2 = max_blk_size + nbw

  if (useGPU) then

    if (allComputeOnGPU) then
      if (wantDebug) call obj%timer%start("cuda_memcpy")

      successGPU = gpu_malloc(q_dev, ldq*matrixCols* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: q_dev", 617,  successGPU)
      successGPU =  gpu_memcpy(q_dev, int(loc(q(1,1)),kind=c_intptr_t),  &
                               ldq*matrixCols * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q -> q_dev", 621,  successGPU)

      ! associate with c_ptr
      q_mpi_dev = transfer(q_dev, q_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(q_mpi_dev, q_mpi_fortran_ptr, &
                       [ldq,matrixCols])
      if (wantDebug) call obj%timer%stop("cuda_memcpy")

      successGPU = gpu_malloc(hh_trans_dev, size(hh_trans,dim=1)*size(hh_trans,dim=2)* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 631,  successGPU)
      ! associate with c_ptr
      hh_trans_mpi_dev = transfer(hh_trans_dev, hh_trans_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(hh_trans_mpi_dev, hh_trans_mpi_fortran_ptr, &
                       [size(hh_trans,dim=1),size(hh_trans,dim=2)])
      successGPU =  gpu_memcpy(c_loc(hh_trans_mpi_fortran_ptr(1,1)),  &
                               c_loc(hh_trans(1,1)), &
                               size(hh_trans,dim=1)*size(hh_trans,dim=2) * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band: hh_trans -> hh_trans_dev", 641,  successGPU)

    endif ! allComputeOnGPU

    num =  (stripe_width*a_dim2*stripe_count)* size_of_datatype
    successGPU = gpu_malloc(aIntern_dev, stripe_width*a_dim2*stripe_count* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 647,  successGPU)

    successGPU = gpu_memset(aIntern_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 650,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      aIntern_mpi_dev = transfer(aIntern_dev, aIntern_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(aIntern_mpi_dev, aIntern_mpi_fortran_ptr, &
                       [stripe_width,a_dim2,stripe_count])
    endif ! allComputeOnGPU

    ! "row_group" and "row_group_dev" are needed for GPU optimizations
    successGPU = gpu_malloc_host(row_group_host,l_nev*nblk*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 667,  successGPU)
    call c_f_pointer(row_group_host, row_group, (/l_nev,nblk/))

    row_group(:, :) = 0.0_rck
    num =  (l_nev*nblk)* size_of_datatype
    successGPU = gpu_malloc(row_group_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 673,  successGPU)

    successGPU = gpu_memset(row_group_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_group_dev", 676,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      row_group_mpi_dev = transfer(row_group_dev, row_group_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(row_group_mpi_dev, row_group_mpi_fortran_ptr, &
                       [l_nev,nblk])
    endif ! allComputeOnGPU

  else ! GPUs are not used



    if (posix_memalign(aIntern_ptr, 64_c_intptr_t, stripe_width*a_dim2*stripe_count*  &
        C_SIZEOF(a_var)) /= 0) then
      print *,"trans_ev_tridi_to_band_real: error when allocating aIntern"//errorMessage
      stop 1
    endif

    call c_f_pointer(aIntern_ptr, aIntern,[stripe_width,a_dim2,stripe_count] )
    !allocate(aIntern(stripe_width,a_dim2,stripe_count), stat=istat, errmsg=errorMessage)

    aIntern(:,:,:) = 0.0_rck
  endif !useGPU

  allocate(row(l_nev), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: row", 723,  istat,  errorMessage)

  row(:) = 0.0_rck

  if (useGPU .and. allComputeOnGPU) then
    num =  (l_nev)* size_of_datatype
    successGPU = gpu_malloc(row_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_dev", 730,  successGPU)

    successGPU = gpu_memset(row_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_dev", 733,  successGPU)

    ! associate with c_ptr
    row_mpi_dev = transfer(row_dev, row_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(row_mpi_dev, row_mpi_fortran_ptr, &
                     [l_nev])
  endif

  ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
  ! and transpose the matrix using stripes of given stripe_width for cache blocking.

  ! The peculiar way it is done below is due to the fact that the last row should be
  ! ready first since it is the first one to start below



  if (wantDebug) call obj%timer%start("ip_loop")
  do ip = np_rows-1, 0, -1
    if (my_prow == ip) then
      ! Receive my rows which have not yet been received
      src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
      do i=limits(ip)+1,limits(ip+1)
        src = mod((i-1)/nblk, np_rows)

        if (src < my_prow) then

          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
            &real&
            &_gpu_&
            &double &
                          ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                                      stripe_width, last_stripe_width, a_dim2, l_nev,&
                                      row_group_size, nblk, unpack_idx, &
                                       i - limits(ip), .false., wantDebug, allComputeOnGPU)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                         int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            call unpack_row_&
                &real&
                &_cpu_&
                &double &
                (obj,aIntern, row,i-limits(ip), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        elseif (src == my_prow) then

          src_offset = src_offset+1

          if (useGPU) then

            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &real&
                 &_gpu_&
                 &double &
             ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                          stripe_width, last_stripe_width, a_dim2, l_nev,&
                          row_group_size, nblk, unpack_idx, &
                          i - limits(ip), .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_DCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_group_mpi_fortran_ptr(1,row_group_size)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            else ! allComputeOnGPU
              row_group(:, row_group_size) = q(src_offset, 1:l_nev)
            endif ! allComputeOnGPU
          else ! useGPU
            row(:) = q(src_offset, 1:l_nev)
          endif ! useGPU


          if (useGPU) then

          else
            call unpack_row_&
                 &real&
                 &_cpu_&
                 &double &
                            (obj,aIntern, row,i-limits(ip),  stripe_count, stripe_width, last_stripe_width)
          endif


        endif ! src == my_prow
      enddo ! i=limits(ip)+1,limits(ip+1)


      ! Send all rows which have not yet been send


      if (allComputeOnGPU .and. useGPU) then
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1

              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_DCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1011,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

              if (wantDebug) call obj%timer%start("cuda_mpi_communication")
              call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

            endif
          enddo
        enddo
      else !  allComputeOnGPU
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1
              row(:) = q(src_offset, 1:l_nev)

              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            endif
          enddo
        enddo
      endif ! allComputeOnGPU

    else if (my_prow < ip) then

      ! Send all rows going to PE ip
      if (allComputeOnGPU .and. useGPU) then
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1

            if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
            call gpublas_DCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
            if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

            ! is there a way to avoid this device_synchronize ?
            if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
            successGPU = gpu_devicesynchronize()
            call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1061,  successGPU)
            if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

            if (wantDebug) call obj%timer%start("cuda_mpi_communication")
            call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                          int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

          endif
        enddo
      else ! allComputeOnGPU 
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1
            row(:) = q(src_offset, 1:l_nev)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                        int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          endif
        enddo
      endif  ! allComputeOnGPU

      ! Receive all rows from PE ip
      do i=limits(my_prow)+1,limits(my_prow+1)
        src = mod((i-1)/nblk, np_rows)
        if (src == ip) then
          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &real&
                 &_gpu_&
                 &double&
                 &( obj, &
                 row_group, row_group_dev, aIntern_dev, stripe_count,  &
                 stripe_width, last_stripe_width, a_dim2, l_nev,       &
                 row_group_size, nblk, unpack_idx,                     &
                 i - limits(my_prow), .false., wantDebug, allComputeOnGPU)

            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_REAL8, &
                             int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
            call unpack_row_&
                 &real&
                 &_cpu_&
                 &double &
                 (obj,aIntern, row,i-limits(my_prow), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        endif
      enddo ! i=limits(my_prow)+1,limits(my_prow+1)
    endif ! (my_prow < ip)
  enddo ! ip = np_rows-1, 0, -1

  if (wantDebug) call obj%timer%stop("ip_loop")
 
  if (wantDebug) call obj%timer%start("allocate")

  if (useGPU) then
    ! Force an unpacking of all remaining rows that haven't been unpacked yet

    call unpack_and_prepare_row_group_&
         &real&
         &_gpu_&
         &double&
         &( obj, &
         row_group, row_group_dev, aIntern_dev, stripe_count, &
         stripe_width, last_stripe_width, &
         a_dim2, l_nev, row_group_size, nblk, unpack_idx,     &
         -1, .true., wantDebug, allComputeOnGPU)

  endif

  ! Set up result buffer queue

  num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

  num_result_buffers = 4*nfact
  allocate(result_buffer(l_nev,nblk,num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_buffer", 1266,  istat,  errorMessage)

  allocate(result_send_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_send_request", 1269,  istat,  errorMessage)

  allocate(result_recv_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_recv_request", 1272,  istat,  errorMessage)

  if (useGPU .and. allComputeOnGPU) then
    num_result_buffers = 4*nfact
    num =  (l_nev*nblk*num_result_buffers)* size_of_datatype
    successGPU = gpu_malloc(result_buffer_dev, num* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 1278,  successGPU)

    ! associate with c_ptr
    result_buffer_mpi_dev = transfer(result_buffer_dev, result_buffer_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(result_buffer_mpi_dev, result_buffer_mpi_fortran_ptr, &
                     [l_nev,nblk,num_result_buffers])
  endif

  result_send_request(:) = MPI_REQUEST_NULL
  result_recv_request(:) = MPI_REQUEST_NULL

  ! Queue up buffers
  if (wantDebug) call obj%timer%start("host_mpi_communication")
  if (my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
    do j = 1, min(num_result_buffers, num_result_blocks)
      if (useGPU) then
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL8,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      else ! useGPU
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL8,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      endif ! useGPU
    enddo
  endif
  if (wantDebug) call obj%timer%stop("host_mpi_communication")

  num_bufs_recvd = 0 ! No buffers received yet

  ! Initialize top/bottom requests

  allocate(top_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_send_request", 1335,  istat,  errorMessage)

  allocate(top_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_recv_request", 1338,  istat,  errorMessage)

  allocate(bottom_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_send_request", 1341,  istat,  errorMessage)

  allocate(bottom_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_recv_request", 1344,  istat,  errorMessage)

  top_send_request(:) = MPI_REQUEST_NULL
  top_recv_request(:) = MPI_REQUEST_NULL
  bottom_send_request(:) = MPI_REQUEST_NULL
  bottom_recv_request(:) = MPI_REQUEST_NULL


  allocate(top_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 1374,  istat,  errorMessage)

  allocate(top_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1377,  istat,  errorMessage)

  allocate(bottom_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1380,  istat,  errorMessage)

  allocate(bottom_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1383,  istat,  errorMessage)

  top_border_send_buffer(:,:) = 0.0_rck
  top_border_recv_buffer(:,:) = 0.0_rck
  bottom_border_send_buffer(:,:) = 0.0_rck
  bottom_border_recv_buffer(:,:) = 0.0_rck


  if (useGPU) then
    if (allComputeOnGPU) then
      ! top_border_recv_buffer and top_border_send_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(top_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1401,  successGPU)

      successGPU = gpu_malloc(top_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1404,  successGPU)

      successGPU = gpu_memset(top_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1407,  successGPU)

      successGPU = gpu_memset(top_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1410,  successGPU)

      ! associate with c_ptr
      top_border_recv_buffer_mpi_dev = transfer(top_border_recv_buffer_dev, top_border_recv_buffer_mpi_dev)
      top_border_send_buffer_mpi_dev = transfer(top_border_send_buffer_dev, top_border_send_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(top_border_recv_buffer_mpi_dev, top_border_recv_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])
      call c_f_pointer(top_border_send_buffer_mpi_dev, top_border_send_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])

      ! bottom_border_send_buffer and bottom_border_recv_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(bottom_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1433,  successGPU)
      successGPU = gpu_malloc(bottom_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1435,  successGPU)

      successGPU = gpu_memset(bottom_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1438,  successGPU)
      successGPU = gpu_memset(bottom_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1440,  successGPU)

      ! associate with c_ptr
      bottom_border_send_buffer_mpi_dev = transfer(bottom_border_send_buffer_dev, bottom_border_send_buffer_mpi_dev)
      bottom_border_recv_buffer_mpi_dev = transfer(bottom_border_recv_buffer_dev, bottom_border_recv_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(bottom_border_send_buffer_mpi_dev, bottom_border_send_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
      call c_f_pointer(bottom_border_recv_buffer_mpi_dev, bottom_border_recv_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
    endif ! allComputeOnGPU

    successGPU = gpu_host_register(int(loc(top_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 1462,  successGPU)

    successGPU = gpu_host_register(int(loc(top_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1467,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1472,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1477,  successGPU)
  endif ! useGPU


  ! Initialize broadcast buffer

  if (useGPU) then
    successGPU = gpu_malloc_host(bcast_buffer_host,nbw*max_blk_size*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 1485,  successGPU)
    call c_f_pointer(bcast_buffer_host, bcast_buffer, (/nbw,max_blk_size/))
  else
    allocate(bcast_buffer(nbw, max_blk_size), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_tridi_to_band: bcast_buffer", 1489,  istat,  errorMessage)
  endif

  bcast_buffer = 0.0_rck

  if (useGPU) then
    num =  ( nbw * max_blk_size) * size_of_datatype
    successGPU = gpu_malloc(bcast_buffer_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1497,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      bcast_buffer_mpi_dev = transfer(bcast_buffer_dev, bcast_buffer_mpi_dev)
      call c_f_pointer(bcast_buffer_mpi_dev, bcast_buffer_mpi_fortran_ptr, &
                      [nbw, max_blk_size])
    endif ! allComputeOnGPU

    successGPU = gpu_memset( bcast_buffer_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1507,  successGPU)

    num =  (max_blk_size)* size_of_datatype
    successGPU = gpu_malloc( hh_tau_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1511,  successGPU)

    successGPU = gpu_memset( hh_tau_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1514,  successGPU)
  endif ! useGPU

  current_tv_off = 0 ! Offset of next row to be broadcast

  ! ------------------- start of work loop -------------------

  a_off = 0 ! offset in aIntern (to avoid unnecessary shifts)

  top_msg_length = 0
  bottom_msg_length = 0
  if (wantDebug) call obj%timer%stop("allocate")

  if (wantDebug) call obj%timer%start("sweep_loop")
  do sweep = 0, (na-1)/nbw

    current_n = na - sweep*nbw
    call determine_workload(obj,current_n, nbw, np_rows, limits)
    current_n_start = limits(my_prow)
    current_n_end   = limits(my_prow+1)
    current_local_n = current_n_end - current_n_start

    next_n = max(current_n - nbw, 0)
    call determine_workload(obj,next_n, nbw, np_rows, limits)
    next_n_start = limits(my_prow)
    next_n_end   = limits(my_prow+1)
    next_local_n = next_n_end - next_n_start

    if (next_n_end < next_n) then
      bottom_msg_length = current_n_end - next_n_end
    else
      bottom_msg_length = 0
    endif

    if (next_local_n > 0) then
      next_top_msg_length = current_n_start - next_n_start
    else
      next_top_msg_length = 0
    endif

    if (sweep==0 .and. current_n_end < current_n .and. l_nev > 0) then
      if (wantDebug) call obj%timer%start("host_mpi_communication")
      do i = 1, stripe_count


        if (useGPU) then
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_REAL8, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        else !useGPU
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_REAL8, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        endif !useGPU

      enddo
      if (wantDebug) call obj%timer%stop("host_mpi_communication")
    endif

    if (current_local_n > 1) then
      if (useGPU .and. allComputeOnGPU) then
        if (my_pcol == mod(sweep,np_cols)) then
          if (wantDebug) call obj%timer%start("cuda_memcpy")
          successGPU =  gpu_memcpy(c_loc(bcast_buffer_mpi_fortran_ptr(1,1)), &
                                   c_loc(hh_trans_mpi_fortran_ptr(1,current_tv_off+1)),  &
                                     size(hh_trans,dim=1) * (current_tv_off+current_local_n-(current_tv_off+1)+1) * &   
                                     size_of_datatype, &
                                     gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1646,  successGPU)
          if (wantDebug) call obj%timer%stop("cuda_memcpy")
          current_tv_off = current_tv_off + current_local_n
        endif
      else !useGPU
        if (my_pcol == mod(sweep,np_cols)) then
          bcast_buffer(:,1:current_local_n) =    &
          hh_trans(:,current_tv_off+1:current_tv_off+current_local_n)
          current_tv_off = current_tv_off + current_local_n
        endif
      endif ! useGPU

      if (useNonBlockingCollectivesCols) then
        if (wantDebug) call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_REAL8, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_bcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_REAL8, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif
 
      if (useGPU .and. .not.(allComputeOnGPU)) then
        if (wantDebug) call obj%timer%start("memcpy")
        successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
                                 nbw * current_local_n *    &
                                 size_of_datatype, &
                                 gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1702,  successGPU)
        if (wantDebug) call obj%timer%stop("memcpy")
      endif ! useGPU

      if (useGPU) then
        !if (.not.(allComputeOnGPU)) then
        !  if (wantDebug) call obj%timer%start("memcpy")
        !  successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
        !                           nbw * current_local_n *    &
        !                           size_of_datatype, &
        !                           gpuMemcpyHostToDevice)
        !  call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1725,  successGPU)
        !  if (wantDebug) call obj%timer%stop("memcpy")
        !endif ! allComputeOnGPU

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &real&
             &_gpu_&
             &double&
             (bcast_buffer_dev, hh_tau_dev, nbw, &
             current_local_n, .false.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU

    else ! (current_local_n > 1) then

      ! for current_local_n == 1 the one and only HH Vector is 0 and not stored in hh_trans_real/complex
      bcast_buffer(:,1) = 0.0_rck
      if (useGPU) then
        successGPU = gpu_memset(bcast_buffer_dev, 0, nbw * size_of_datatype)
        call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1745,  successGPU)

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &real&
             &_gpu_&
             &double&
             &( &
             bcast_buffer_dev, hh_tau_dev, &
             nbw, 1, .true.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU
    endif ! (current_local_n > 1) then

    if (l_nev == 0) cycle

    if (current_local_n > 0) then

      do i = 1, stripe_count

        !wait_b
        if (current_n_end < current_n) then



          if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_recv")
          call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)

          if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_recv")
          n_off = current_local_n+a_off

          if (useGPU) then
            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_memcpy")
              successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                       c_loc(bottom_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                       stripe_width*nbw* size_of_datatype,      &
                                       gpuMemcpyDeviceToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1859,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_memcpy")
            else ! allComputeOnGPU
              if (wantDebug) call obj%timer%start("memcpy")
              dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * size_of_datatype
              successGPU =  gpu_memcpy( aIntern_dev + dev_offset , &
                                      int(loc(bottom_border_recv_buffer(1,i)),kind=c_intptr_t), &
                                       stripe_width*nbw*  size_of_datatype,    &
                                       gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1868,  successGPU)
              if (wantDebug) call obj%timer%stop("memcpy")
            endif ! allComputeOnGPU
          else ! useGPU
            aIntern(:,n_off+1:n_off+nbw,i) = reshape( &
            bottom_border_recv_buffer(1:stripe_width*nbw,i),(/stripe_width,nbw/))
          endif ! useGPU


          if (next_n_end < next_n) then

            if (useGPU) then
              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_REAL8, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            else ! useGPU
              if (wantDebug) call obj%timer%start("mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_REAL8, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("mpi_communication")
            endif ! useGPU

            endif ! (next_n_end < next_n)
          endif ! (current_n_end < current_n)

          if (current_local_n <= bottom_msg_length + top_msg_length) then

            !wait_t
            if (top_msg_length>0) then


              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! the MPI_IRECV will be done CUDA_AWARE we thus do not need a host to device copy
                  ! However, we have to copy from top_border_recv_buffer_mpi_fortran_ptr to aIntern_mpi_fortran_ptr
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  !Fortran pointer for indexing
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                          c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                          stripe_width*top_msg_length* size_of_datatype,      &
                                          gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2045,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  !             host_offset= (0 + (0 * stripe_width) + ( (i-1) * stripe_width * nbw ) ) * 8
                  successGPU =  gpu_memcpy( aIntern_dev+dev_offset , int(loc(top_border_recv_buffer(1,i)),kind=c_intptr_t),  &
                                             stripe_width*top_msg_length* size_of_datatype,      &
                                             gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2054,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
            endif ! top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
                &real&
                &_&
                &double&
                &(obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
                a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
                hh_tau_dev, kernel_flops, kernel_time, n_times, 0, current_local_n, i, &
                last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b        1
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")
            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")

            if (bottom_msg_length>0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer

                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                           c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2248,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")



              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 2255,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy( int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2263,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                      MPI_REAL8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                      int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                     MPI_REAL8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                     int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU
            endif !(bottom_msg_length>0)

          else ! current_local_n <= bottom_msg_length + top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &real&
             &_&
             &double&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, &
             current_local_n - bottom_msg_length, bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")

            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")
            if (bottom_msg_length > 0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                            c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                             stripe_width * bottom_msg_length * size_of_datatype,      &
                                             gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2516,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset,  &
                                           stripe_width*bottom_msg_length* size_of_datatype,  &
                                           gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2524,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_REAL8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
          
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_REAL8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU


            endif ! (bottom_msg_length > 0)


            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &real&
             &_&
             &double&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, top_msg_length, &
             current_local_n-top_msg_length-bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")


            !wait_t
            if (top_msg_length>0) then

              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
              if (useGPU) then
                if (allComputeOnGPU) then
                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                           c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                           stripe_width* top_msg_length* size_of_datatype,      &
                                           gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2739,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  ! copy top_border_recv_buffer to aIntern_dev, maybe not necessary if CUDA_AWARE IRECV
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(aIntern_dev + dev_offset ,int(loc( top_border_recv_buffer(:,i)),kind=c_intptr_t),  &
                                        stripe_width * top_msg_length * size_of_datatype,   &
                                        gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2748,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
           endif

           !compute

           if (wantDebug) call obj%timer%start("compute_hh_trafo")
           call compute_hh_trafo_&
             &real&
             &_&
             &double&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off, nbw, max_blk_size,  bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, 0, top_msg_length, i, &
             last_stripe_width, kernel)
           if (wantDebug) call obj%timer%stop("compute_hh_trafo")

         endif

         if (next_top_msg_length > 0) then
           !request top_border data

           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_REAL8, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else !useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_REAL8, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif !useGPU


         endif

         !send_t
         if (my_prow > 0) then

           if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
           call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")

           if (useGPU) then
             if (allComputeOnGPU) then
               if (wantDebug) call obj%timer%start("cuda_memcpy")
               successGPU =  gpu_memcpy(c_loc(top_border_send_buffer_mpi_fortran_ptr(1,i)),  &
                                        c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                        stripe_width* nbw* size_of_datatype,      &
                                        gpuMemcpyDeviceToDevice)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 3030,  successGPU)
               if (wantDebug) call obj%timer%stop("cuda_memcpy")
             else ! allComputeOnGPU
               if (wantDebug) call obj%timer%start("memcpy")
               dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
               successGPU =  gpu_memcpy(int(loc(top_border_send_buffer(:,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                         stripe_width*nbw * size_of_datatype, &
                                         gpuMemcpyDeviceToHost)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> top_border_send_buffer", 3038,  successGPU)
               if (wantDebug) call obj%timer%stop("memcpy")
             endif ! allComputeOnGPU
           else ! useGPU
             top_border_send_buffer(:,i) = reshape(aIntern(:,a_off+1:a_off+nbw,i),(/stripe_width*nbw/))
           endif ! useGPU
           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_REAL8, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else ! useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_REAL8, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif ! useGPU

         endif

         ! Care that there are not too many outstanding top_recv_request's
         if (stripe_count > 1) then
           if (i>1) then

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
             call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
           else

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
            call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

          endif
        endif
      enddo ! i = 1, stripe_count

      top_msg_length = next_top_msg_length

    else
      ! wait for last top_send_request

      do i = 1, stripe_count
        if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
        call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")
      enddo
    endif

    ! Care about the result

    if (my_prow == 0) then

      ! topmost process sends nbw rows to destination processes

      do j=0, nfact-1
        num_blk = sweep*nfact+j ! global number of destination block, 0 based
        if (num_blk*nblk >= na) exit

        nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block

        if (wantDebug) call obj%timer%start("host_mpi_wait_result_send")
        call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_result_send")
        dst = mod(num_blk, np_rows)

        if (dst == 0) then
          if (useGPU) then
            row_group_size = min(na - num_blk*nblk, nblk)

            call pack_row_group_&
                 &real&
                 &_gpu_&
                 &double&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
                         row_group(:, :), j * nblk + a_off, row_group_size, &
                         result_buffer_dev, nblk, num_result_buffers, nbuf, .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              ! memcpy DeviceToDevice row_group_dev -> q_dev
              do i = 1, row_group_size  
               if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
               call gpublas_DCOPY(l_nev, c_loc(row_group_mpi_fortran_ptr(1,i)), 1, &
                                           c_loc(q_mpi_fortran_ptr((num_blk / np_rows) * nblk + i,1)), ldq)
               if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
              enddo
            else ! allComputeOnGPU
              do i = 1, row_group_size
                q((num_blk / np_rows) * nblk + i, 1 : l_nev) = row_group(:, i)
              enddo
            endif ! allComputeOnGPU
          else ! useGPU

            do i = 1, min(na - num_blk*nblk, nblk)

              call pack_row_&
                   &real&
                   &_cpu_&
                   &double&
                   &(obj,aIntern, row, j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
              q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
            enddo
          endif ! useGPU

        else ! (dst == 0)

          if (useGPU) then

            call pack_row_group_&
                 &real&
                 &_gpu_&
                 &double&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, &
                   last_stripe_width, a_dim2, l_nev, &
                   result_buffer(:, :, nbuf), j * nblk + a_off, nblk, &
                   result_buffer_dev, nblk, num_result_buffers, nbuf, .true., wantDebug, allComputeOnGPU)

          else  ! useGPU
            do i = 1, nblk
              call pack_row_&
                   &real&
                   &_cpu_&
                   &double&
                   &(obj, aIntern, result_buffer(:,i,nbuf),j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
            enddo
          endif ! useGPU
         
          if (useGPU) then
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL8, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL8, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif ! useGPU
        endif ! (dst == 0)
      enddo  !j=0, nfact-1

    else ! (my_prow == 0)

      ! receive and store final result

      do j = num_bufs_recvd, num_result_blocks-1

        nbuf = mod(j, num_result_buffers) + 1 ! buffer number to get this block

        ! If there is still work to do, just test for the next result request
        ! and leave the loop if it is not ready, otherwise wait for all
        ! outstanding requests

        if (next_local_n > 0) then
            ! needed
            !successGPU = gpu_devicesynchronize()
          if (wantDebug) call obj%timer%start("host_mpi_test")
          call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_test")

          if (.not.flag) exit

        else ! (next_local_n > 0)
          if (wantDebug) call obj%timer%start("host_mpi_wait_result_recv")

          call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_wait_result_recv")
        endif ! (next_local_n > 0)

        ! Fill result buffer into q
        num_blk = j*np_rows + my_prow ! global number of current block, 0 based
        if (useGPU) then
          if (allComputeOnGPU) then
            do i = 1, min(na - num_blk*nblk, nblk)
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_DCOPY(l_nev, c_loc(result_buffer_mpi_fortran_ptr(1,i,nbuf)), 1, &
                                          c_loc(q_mpi_fortran_ptr(j*nblk + i,1)), ldq)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            enddo
          else ! allComputeOnGPU
            do i = 1, min(na - num_blk*nblk, nblk)
              q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
            enddo
          endif ! allComputeOnGPU
        else ! useGPU
          do i = 1, min(na - num_blk*nblk, nblk)
            q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
          enddo
        endif ! useGPU
        ! Queue result buffer again if there are outstanding blocks left

        if (useGPU) then
          if (wantDebug) call obj%timer%start("host_mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL8, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("host_mpi_communication")
        else ! useGPU
          if (wantDebug) call obj%timer%start("mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL8, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif !useGPU

      enddo ! j = num_bufs_recvd, num_result_blocks-1
      num_bufs_recvd = j

    endif ! (my_prow == 0)

    ! Shift the remaining rows to the front of aIntern (if necessary)

    offset = nbw - top_msg_length
    if (offset<0) then
      if (wantDebug) write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                                         &real&
                                         &: internal error, offset for shifting = ',offset
      success = .false.
      return
    endif

    a_off = a_off + offset
    if (a_off + next_local_n + nbw >= a_dim2) then
      do i = 1, stripe_count
        if (useGPU) then
          chunk = min(next_local_n,a_off)

          if (chunk < 1) exit

          if (wantDebug) call obj%timer%start("normal_memcpy")
          do j = top_msg_length+1, top_msg_length+next_local_n, chunk
            this_chunk = min(j+chunk-1,top_msg_length+next_local_n)-j+1
            dev_offset = ((j-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            dev_offset_1 = ((j+a_off-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            num = stripe_width*this_chunk*size_of_datatype
            successGPU = gpu_memcpy(aIntern_dev+dev_offset, aIntern_dev+dev_offset_1, num, gpuMemcpyDeviceToDevice)

            call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> aIntern_dev", 3528,  successGPU)
          end do
          if (wantDebug) call obj%timer%stop("normal_memcpy")
        else ! not useGPU
          do j = top_msg_length+1, top_msg_length+next_local_n
            aIntern(:,j,i) = aIntern(:,j+a_off,i)
          end do
        end if
      end do ! stripe_count

      a_off = 0
    end if
  end do ! sweep
  if (wantDebug) call obj%timer%stop("sweep_loop")

  ! Just for safety:
  if (ANY(top_send_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_send_request ***',my_prow,my_pcol
  if (ANY(bottom_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_send_request ***',my_prow,my_pcol
  if (ANY(top_recv_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_recv_request ***',my_prow,my_pcol
  if (ANY(bottom_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_recv_request ***',my_prow,my_pcol

  if (my_prow == 0) then

    if (wantDebug) call obj%timer%start("host_mpi_waitall")
    call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("host_mpi_waitall")
  endif

  if (ANY(result_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_send_request ***',my_prow,my_pcol
  if (ANY(result_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_recv_request ***',my_prow,my_pcol

  call obj%get("print_flops",print_flops,error)

  if (print_flops == 1) then

    if (useNonBlockingCollectivesRows) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          allreduce_request1, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          allreduce_request3, mpierr)

      call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif

    if (useNonBlockingCollectivesCols) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                          allreduce_request2, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          allreduce_request4, mpierr)

      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                         mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif
  endif

  if (useGPU .and. allComputeOnGPU) then
    ! finally copy q_dev to q
    if (wantDebug) call obj%timer%start("cuda_memcpy")
    successGPU =  gpu_memcpy(int(loc(q(1,1)),kind=c_intptr_t),  &
                             q_dev, &
                             ldq*matrixCols * size_of_datatype, &
                             gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q_dev -> q", 3634,  successGPU)
    if (wantDebug) call obj%timer%stop("cuda_memcpy")

  endif

  if (my_prow==0 .and. my_pcol==0 .and.print_flops == 1) &
      write(error_unit,'(" Kernel time:",f10.3," MFlops: ",es12.5)')  kernel_time, kernel_flops/kernel_time*1.d-6

  ! deallocate all working space

  if (.not.(useGPU)) then
    nullify(aIntern)
    call free(aIntern_ptr)
  endif

  deallocate(row, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: row", 3650,  istat,  errorMessage)

  deallocate(limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: limits", 3653,  istat,  errorMessage)

  deallocate(result_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_send_request", 3656,  istat,  errorMessage)

  deallocate(result_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_recv_request", 3659,  istat,  errorMessage)

  deallocate(result_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_buffer", 3662,  istat,  errorMessage)

  if (useGPU) then
    if (allComputeOnGPU) then
      successGPU = gpu_free(result_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 3667,  successGPU)
      nullify(result_buffer_mpi_fortran_ptr)
    endif

    nullify(bcast_buffer)

    successGPU = gpu_free_host(bcast_buffer_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 3674,  successGPU)
  else
    deallocate(bcast_buffer, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_tridi_to_band: bcast_buffer", 3677,  istat,  errorMessage)
  endif


  if (useGPU) then
    successGPU = gpu_free(aIntern_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 3683,  successGPU)

    if (allComputeOnGPU) then
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: q_dev", 3687,  successGPU)
      nullify(q_mpi_fortran_ptr)

      successGPU = gpu_free(hh_trans_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 3691,  successGPU)
      nullify(hh_trans_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 3695,  successGPU)
      nullify(top_border_recv_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 3699,  successGPU)
      nullify(top_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 3703,  successGPU)
      nullify(bottom_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 3707,  successGPU)
      nullify(bottom_border_recv_buffer_mpi_fortran_ptr)

      nullify(aIntern_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU = gpu_free(hh_tau_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 3714,  successGPU)

    nullify(row_group)

    successGPU = gpu_free_host(row_group_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 3719,  successGPU)

    successGPU = gpu_free(row_group_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 3722,  successGPU)

    if (allComputeOnGPU) then
      nullify(row_group_mpi_fortran_ptr)

      successGPU = gpu_free(row_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_dev", 3728,  successGPU)
      nullify(row_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU =  gpu_free(bcast_buffer_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 3733,  successGPU)

    if (allComputeOnGPU) then
      nullify(bcast_buffer_mpi_fortran_ptr)
    endif

    successGPU = gpu_host_unregister(int(loc(top_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 3740,  successGPU)

    successGPU = gpu_host_unregister(int(loc(top_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3743,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3746,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3749,  successGPU)
  endif ! useGPU

  deallocate(top_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 3753,  istat,  errorMessage)

  deallocate(top_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3756,  istat,  errorMessage)

  deallocate(bottom_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3759,  istat,  errorMessage)

  deallocate(bottom_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3762,  istat,  errorMessage)

  deallocate(top_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_send_request", 3765,  istat,  errorMessage)

  deallocate(top_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_recv_request", 3768,  istat,  errorMessage)

  deallocate(bottom_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_send_request", 3771,  istat,  errorMessage)

  deallocate(bottom_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_recv_request", 3774,  istat,  errorMessage)

  call obj%timer%stop("trans_ev_tridi_to_band_&
                      &real&
                      &" // &
                      &"_double" //&
                      gpuString)
  return

end subroutine

! vim: syntax=fortran



    subroutine band_band_real_&
&double &
                  (obj, na, nb, nbCol, nb2, nb2Col, ab, ab2, d, e, communicator)
    !-------------------------------------------------------------------------------
    ! band_band_real:
    ! Reduces a real symmetric banded matrix to a real symmetric matrix with smaller bandwidth. Householder transformations are not stored.
    ! Matrix size na and original bandwidth nb have to be a multiple of the target bandwidth nb2. (Hint: expand your matrix with
    ! zero entries, if this
    ! requirement doesn't hold)
    !
    !  na          Order of matrix
    !
    !  nb          Semi bandwidth of original matrix
    !
    !  nb2         Semi bandwidth of target matrix
    !
    !  ab          Input matrix with bandwidth nb. The leading dimension of the banded matrix has to be 2*nb. The parallel data layout
    !              has to be accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb+1 to min(na, block_limits(n+1)*nb)
    !              are located on rank n.
    !
    !  ab2         Output matrix with bandwidth nb2. The leading dimension of the banded matrix is 2*nb2. The parallel data layout is
    !              accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb2+1 to min(na, block_limits(n+1)*nb2) are located
    !              on rank n.
    !
    !  d(na)       Diagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
    !
    !  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
    !
    !  communicator
    !              MPI-Communicator for the total processor set
    !-------------------------------------------------------------------------------
      use elpa_abstract_impl
      use elpa2_workload
      use elpa_blas_interfaces

      use precision
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
      integer(kind=ik), intent(in)               :: na, nb, nbCol, nb2, nb2Col, communicator
      real(kind=rk), intent(inout)               :: ab(2*nb,nbCol) ! removed assumed size
      real(kind=rk), intent(inout)               :: ab2(2*nb2,nb2Col) ! removed assumed size
      real(kind=rk), intent(out)                 :: d(na), e(na) ! set only on PE 0

      real(kind=rk)                              :: hv(nb,nb2), w(nb,nb2), w_new(nb,nb2), tau(nb2), hv_new(nb,nb2), &
                                                  tau_new(nb2), ab_s(1+nb,nb2), ab_r(1+nb,nb2), ab_s2(2*nb2,nb2), hv_s(nb,nb2)

      real(kind=rk)                              :: work(nb*nb2), work2(nb2*nb2)
      integer(kind=ik)                         :: lwork, info
      integer(kind=BLAS_KIND)                  :: infoBLAS

      integer(kind=ik)                         :: istep, i, n, dest
      integer(kind=ik)                         :: n_off, na_s
      integer(kind=ik)                         :: my_pe, n_pes
      integer(kind=MPI_KIND)                   :: my_peMPI, n_pesMPI, mpierr
      integer(kind=ik)                         :: nblocks_total, nblocks
      integer(kind=ik)                         :: nblocks_total2, nblocks2
      integer(kind=MPI_KIND)                   :: ireq_ab, ireq_hv
!      integer(kind=ik)                         :: MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
!      integer(kind=ik), allocatable            :: mpi_statuses(:,:)
      integer(kind=ik), allocatable            :: block_limits(:), block_limits2(:)
      integer(kind=MPI_KIND), allocatable      :: ireq_ab2(:)

      integer(kind=ik)                         :: j, nc, nr, ns, ne, iblk
      integer(kind=ik)                         :: istat
      character(200)                           :: errorMessage

      call obj%timer%start("band_band_real" // "_double")

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(communicator,kind=MPI_KIND) ,my_peMPI ,mpierr)
      call mpi_comm_size(int(communicator,kind=MPI_KIND) ,n_pesMPI ,mpierr)

      my_pe = int(my_peMPI,kind=c_int)
      n_pes = int(n_pesMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")

      ! Total number of blocks in the band:
      nblocks_total = (na-1)/nb + 1
      nblocks_total2 = (na-1)/nb2 + 1

      ! Set work distribution
      allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating block_limits "//errorMessage
        stop 1
      endif
      call divide_band(obj, nblocks_total, n_pes, block_limits)

      allocate(block_limits2(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating block_limits2 "//errorMessage
        stop 1
      endif

      call divide_band(obj, nblocks_total2, n_pes, block_limits2)

      ! nblocks: the number of blocks for my task
      nblocks = block_limits(my_pe+1) - block_limits(my_pe)
      nblocks2 = block_limits2(my_pe+1) - block_limits2(my_pe)

      allocate(ireq_ab2(1:nblocks2), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating ireq_ab2 "//errorMessage
        stop 1
      endif

      call obj%timer%start("mpi_communication")

      ireq_ab2 = MPI_REQUEST_NULL

      if (nb2>1) then
        do i=0,nblocks2-1

          call mpi_irecv(ab2(1,i*nb2+1), int(2*nb2*nb2,kind=MPI_KIND), MPI_REAL8, &
                         0_MPI_KIND, 3_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_ab2(i+1), mpierr)
        enddo
      endif
      call obj%timer%stop("mpi_communication")

      ! n_off: Offset of ab within band
      n_off = block_limits(my_pe)*nb
      lwork = nb*nb2
      dest = 0
      ireq_ab = MPI_REQUEST_NULL
      ireq_hv = MPI_REQUEST_NULL
      ! ---------------------------------------------------------------------------
      ! Start of calculations

      na_s = block_limits(my_pe)*nb + 1

      if (my_pe>0 .and. na_s<=na) then
        ! send first nb2 columns to previous PE
        ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
        do i=1,nb2
          ab_s(1:nb+1,i) = ab(1:nb+1,na_s-n_off+i-1)
        enddo
        call obj%timer%start("mpi_communication")

        call mpi_isend(ab_s, int((nb+1)*nb2,kind=MPI_KIND), MPI_REAL8, int(my_pe-1,kind=MPI_KIND), &
                       1_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_ab, mpierr)
        call obj%timer%stop("mpi_communication")
      endif

      do istep=1,na/nb2

        if (my_pe==0) then

          n = MIN(na-na_s-nb2+1,nb) ! number of rows to be reduced
          hv(:,:) = 0.0_rk
          tau(:) = 0.0_rk

          ! The last step (istep=na-1) is only needed for sending the last HH vectors.
          ! We don't want the sign of the last element flipped (analogous to the other sweeps)
          if (istep < na/nb2) then

            ! Transform first block column of remaining matrix
      call obj%timer%start("blas")
            call DGEQRF(int(n,kind=BLAS_KIND), int(nb2,kind=BLAS_KIND), ab(1+nb2,na_s-n_off), &
                                 int(2*nb-1,kind=BLAs_KIND), tau, work, int(lwork,kind=BLAS_KIND), &
                                 infoBLAS)
      info = int(infoBLAS,kind=ik)
      call obj%timer%stop("blas")

            do i=1,nb2
              hv(i,i) = 1.0_rk
              hv(i+1:n,i) = ab(1+nb2+1:1+nb2+n-i,na_s-n_off+i-1)
              ab(1+nb2+1:2*nb,na_s-n_off+i-1) = 0.0_rk
            enddo

          endif

          if (nb2==1) then
            d(istep) = ab(1,na_s-n_off)
            e(istep) = ab(2,na_s-n_off)
            if (istep == na) then
              e(na) = 0.0_rk
            endif
          else
            ab_s2 = 0.0_rk
            ab_s2(:,:) = ab(1:nb2+1,na_s-n_off:na_s-n_off+nb2-1)
            if (block_limits2(dest+1)<istep) then
              dest = dest+1
            endif
            call obj%timer%start("mpi_communication")
            call mpi_send(ab_s2, int(2*nb2*nb2,kind=MPI_KIND), MPI_REAL8, int(dest,kind=MPI_KIND), &
                          3_MPI_KIND, int(communicator,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")


          endif

        else
          if (na>na_s+nb2-1) then
            ! Receive Householder vectors from previous task, from PE owning subdiagonal
            call obj%timer%start("mpi_communication")
            call mpi_recv(hv, int(nb*nb2,kind=MPI_KIND), MPI_REAL8, int(my_pe-1,kind=MPI_KIND), &
                          2_MPI_KIND, int(communicator,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            call obj%timer%stop("mpi_communication")


            do i=1,nb2
              tau(i) = hv(i,i)
              hv(i,i) = 1.0_rk
            enddo
          endif
        endif

        na_s = na_s+nb2
        if (na_s-n_off > nb) then
          ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
          ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rk
          n_off = n_off + nb
        endif

        do iblk=1,nblocks
          ns = na_s + (iblk-1)*nb - n_off ! first column in block
          ne = ns+nb-nb2                    ! last column in block

          if (ns+n_off>na) exit

            nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
            nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                          ! Note that nr>=0 implies that diagonal block is full (nc==nb)!
            call wy_gen_&
      &double&
      &(obj,nc,nb2,w,hv,tau,work,nb)

            if (iblk==nblocks .and. nc==nb) then
              !request last nb2 columns
              call obj%timer%start("mpi_communication")
              call mpi_recv(ab_r, int((nb+1)*nb2,kind=MPI_KIND), MPI_REAL8, int(my_pe+1,kind=MPI_KIND), &
                            1_MPI_KIND, int(communicator,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
              call obj%timer%stop("mpi_communication")

              do i=1,nb2
                ab(1:nb+1,ne+i-1) = ab_r(:,i)
              enddo
            endif
            hv_new(:,:) = 0.0_rk ! Needed, last rows must be 0 for nr < nb
            tau_new(:) = 0.0_rk

            if (nr>0) then
              call wy_right_&
        &double&
        &(obj,nr,nb,nb2,ab(nb+1,ns),2*nb-1,w,hv,work,nb)
        call obj%timer%start("blas")
              call DGEQRF(int(nr,kind=BLAS_KIND), int(nb2,kind=BLAS_KIND), ab(nb+1,ns), &
                                   int(2*nb-1,kind=BLAS_KIND), tau_new, work, int(lwork,kind=BLAS_KIND), &
                                   infoBLAS)
        info = int(infoBLAS,kind=ik)
        call obj%timer%stop("blas")
              do i=1,nb2
                hv_new(i,i) = 1.0_rk
                hv_new(i+1:,i) = ab(nb+2:2*nb-i+1,ns+i-1)
                ab(nb+2:,ns+i-1) = 0.0_rk
              enddo

              !send hh-Vector
              if (iblk==nblocks) then
                call obj%timer%start("mpi_communication")

                call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
                call obj%timer%stop("mpi_communication")

                hv_s = hv_new
                do i=1,nb2
                  hv_s(i,i) = tau_new(i)
                enddo
                call obj%timer%start("mpi_communication")
                call mpi_isend(hv_s, int(nb*nb2,kind=MPI_KIND), MPI_REAL8, int(my_pe+1,kind=MPI_KIND), &
                               2_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_hv, mpierr)
                call obj%timer%stop("mpi_communication")

              endif
            endif

            call wy_symm_&
      &double&
      &(obj,nc,nb2,ab(1,ns),2*nb-1,w,hv,work,work2,nb)

            if (my_pe>0 .and. iblk==1) then
              !send first nb2 columns to previous PE
              call obj%timer%start("mpi_communication")

              call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
              call obj%timer%stop("mpi_communication")

              do i=1,nb2
                ab_s(1:nb+1,i) = ab(1:nb+1,ns+i-1)
              enddo
              call obj%timer%start("mpi_communication")
              call mpi_isend(ab_s, int((nb+1)*nb2,kind=MPI_KIND), MPI_REAL8, int(my_pe-1,kind=MPI_KIND), &
                             1_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_ab, mpierr)
              call obj%timer%stop("mpi_communication")

            endif

            if (nr>0) then
              call wy_gen_&
        &double&
        &(obj,nr,nb2,w_new,hv_new,tau_new,work,nb)
              call wy_left_&
        &double&
        &(obj,nb-nb2,nr,nb2,ab(nb+1-nb2,ns+nb2),2*nb-1,w_new,hv_new,work,nb)
            endif

            ! Use new HH Vector for the next block
            hv(:,:) = hv_new(:,:)
            tau = tau_new
          enddo
        enddo

        ! Finish the last outstanding requests
         call obj%timer%start("mpi_communication")

        call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
        call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
!        allocate(mpi_statuses(MPI_STATUS_SIZE,nblocks2), stat=istat, errmsg=errorMessage)
!        if (istat .ne. 0) then
!          print *,"error allocating mpi_statuses "//errorMessage
!          stop 1
!        endif

        call mpi_waitall(nblocks2,ireq_ab2,MPI_STATUSES_IGNORE,mpierr)
!        deallocate(mpi_statuses, stat=istat, errmsg=errorMessage)
!        if (istat .ne. 0) then
!          print *,"error deallocating mpi_statuses "//errorMessage
!          stop 1
!        endif

        call mpi_barrier(int(communicator,kind=MPI_KIND) ,mpierr)
        call obj%timer%stop("mpi_communication")


        deallocate(block_limits, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating block_limits "//errorMessage
          stop 1
        endif

        deallocate(block_limits2, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating block_limits2 "//errorMessage
          stop 1
        endif

        deallocate(ireq_ab2, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating ireq_ab2 "//errorMessage
          stop 1
        endif

        call obj%timer%stop("band_band_real" // "_double")

    end subroutine

    subroutine wy_gen_&
    &double&
    &(obj, n, nb, W, Y, tau, mem, lda)

      use elpa_abstract_impl
      use elpa_blas_interfaces

      use precision
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
      integer(kind=ik), intent(in)            :: n      !length of householder-vectors
      integer(kind=ik), intent(in)            :: nb     !number of householder-vectors
      integer(kind=ik), intent(in)            :: lda        !leading dimension of Y and W
      real(kind=rk), intent(in)               :: Y(lda,nb)  !matrix containing nb householder-vectors of length b
      real(kind=rk), intent(in)               :: tau(nb)    !tau values
      real(kind=rk), intent(out)              :: W(lda,nb)  !output matrix W
      real(kind=rk), intent(in)               :: mem(nb)    !memory for a temporary matrix of size nb

      integer(kind=ik)                        :: i

   call obj%timer%start("wy_gen" // "_double")

   W(1:n,1) = tau(1)*Y(1:n,1)
   do i=2,nb
     W(1:n,i) = tau(i)*Y(1:n,i)
     call obj%timer%start("blas")
     call DGEMV('T', int(n,kind=BLAS_KIND), int(i-1,kind=BLAS_KIND),  1.0_rk, Y, int(lda,kind=BLAS_KIND), &
                         W(1,i), 1_BLAS_KIND, 0.0_rk, mem, 1_BLAS_KIND)
     call DGEMV('N', int(n,kind=BLAS_KIND), int(i-1,kind=BLAS_KIND), -1.0_rk, W, int(lda,kind=BLAS_KIND), &
                         mem, 1_BLAS_KIND, 1.0_rk, W(1,i), 1_BLAS_KIND)
     call obj%timer%stop("blas")
   enddo
   call obj%timer%stop("wy_gen" // "_double")
    end subroutine

    subroutine wy_left_&
    &double&
    &(obj, n, m, nb, A, lda, W, Y, mem, lda2)

      use precision
      use elpa_abstract_impl
      use elpa_blas_interfaces
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
      integer(kind=ik), intent(in)            :: n      !width of the matrix A
      integer(kind=ik), intent(in)            :: m      !length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk), intent(inout)            :: A(lda,*)   !matrix to be transformed   ! remove assumed size
      real(kind=rk), intent(in)               :: W(m,nb)    !blocked transformation matrix W
      real(kind=rk), intent(in)               :: Y(m,nb)    !blocked transformation matrix Y
      real(kind=rk), intent(inout)            :: mem(n,nb)  !memory for a temporary matrix of size n x nb

   call obj%timer%start("wy_left" // "_double")
   call obj%timer%start("blas")
   call DGEMM('T', 'N', int(nb,kind=BLAS_KIND), int(n,kind=BLAS_KIND), int(m,kind=BLAS_KIND), &
                       1.0_rk, W, int(lda2,kind=BLAS_KIND), A, int(lda,kind=BLAS_KIND), 0.0_rk, mem, &
                       int(nb,kind=BLAS_KIND))
   call DGEMM('N', 'N', int(m,kind=BLAS_KIND), int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), &
                       -1.0_rk, Y, int(lda2,kind=BLAS_KIND), mem, int(nb,kind=BLAS_KIND), 1.0_rk, A, int(lda,kind=BLAS_KIND))
   call obj%timer%stop("blas")
   call obj%timer%stop("wy_left" // "_double")
    end subroutine

    subroutine wy_right_&
    &double&
    &(obj, n, m, nb, A, lda, W, Y, mem, lda2)

      use precision
      use elpa_abstract_impl
      use elpa_blas_interfaces
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
      integer(kind=ik), intent(in)            :: n      !height of the matrix A
      integer(kind=ik), intent(in)            :: m      !length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk), intent(inout)            :: A(lda,*)   !matrix to be transformed  ! remove assumed size
      real(kind=rk), intent(in)               :: W(m,nb)    !blocked transformation matrix W
      real(kind=rk), intent(in)               :: Y(m,nb)    !blocked transformation matrix Y
      real(kind=rk), intent(inout)            :: mem(n,nb)  !memory for a temporary matrix of size n x nb


      call obj%timer%start("wy_right" // "_double")
      call obj%timer%start("blas")
      call DGEMM('N', 'N', int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), int(m,kind=BLAS_KIND), &
                          1.0_rk, A, int(lda,kind=BLAS_KIND), W, int(lda2,kind=BLAS_KIND), 0.0_rk, mem, int(n,kind=BLAS_KIND))
      call DGEMM('N', 'T', int(n,kind=BLAS_KIND), int(m,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), &
                          -1.0_rk, mem, int(n,kind=BLAS_KIND), Y, int(lda2,kind=BLAS_KIND), 1.0_rk, A, int(lda,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      call obj%timer%stop("wy_right" // "_double")

    end subroutine

    subroutine wy_symm_&
    &double&
    &(obj, n, nb, A, lda, W, Y, mem, mem2, lda2)

      use elpa_abstract_impl
      use elpa_blas_interfaces

      use precision
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
      integer(kind=ik), intent(in)            :: n      !width/heigth of the matrix A; length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk), intent(inout)            :: A(lda,*)   !matrix to be transformed  ! remove assumed size
      real(kind=rk), intent(in)               :: W(n,nb)    !blocked transformation matrix W
      real(kind=rk), intent(in)               :: Y(n,nb)    !blocked transformation matrix Y
      real(kind=rk)                           :: mem(n,nb)  !memory for a temporary matrix of size n x nb
      real(kind=rk)                           :: mem2(nb,nb)    !memory for a temporary matrix of size nb x nb

      call obj%timer%start("wy_symm" // "_double")
      call obj%timer%start("blas")
      call DSYMM('L', 'L', int(n, kind=BLAS_KIND), int(nb,kind=BLAS_KIND), 1.0_rk, A, &
                          int(lda,kind=BLAS_KIND), W, int(lda2,kind=BLAS_KIND), 0.0_rk, mem, int(n,kind=BLAS_KIND))
      call DGEMM('T', 'N', int(nb,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), int(n,kind=BLAS_KIND), &
                          1.0_rk, mem, int(n,kind=BLAS_KIND), W, int(lda2,kind=BLAS_KIND), 0.0_rk, mem2, &
                          int(nb,kind=BLAS_KIND))
      call DGEMM('N', 'N', int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), &
                          -0.5_rk, Y, int(lda2,kind=BLAS_KIND), mem2, int(nb,kind=BLAS_KIND), 1.0_rk, mem, int(n,kind=BLAS_KIND))
      call DSYR2K('L', 'N',int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -1.0_rk, Y, int(lda2,kind=BLAS_KIND), &
                           mem, int(n,kind=BLAS_KIND), 1.0_rk, A, int(lda,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      call obj%timer%stop("wy_symm" // "_double")

    end subroutine


! real single precision
























!cannot use "../src/elpa2/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)



subroutine bandred_&
&real&
&_&
&single &
(obj, na, a_mat, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, tmat, &
wantDebug, useGPU, success, &
useQR, &
max_threads, isSkewsymmetric)

!-------------------------------------------------------------------------------
!  bandred_real/complex: Reduces a distributed symmetric matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a_mat(lda,matrixCols)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a_mat(:,:) must be set completely (upper and lower half)
!              a_mat(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,numBlocks)    where numBlocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa1_compute
  use precision
  use elpa_blas_interfaces
  use elpa_scalapack_interfaces
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

  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                            :: na, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols

  real(kind=rck)                     :: a_mat(lda,*)
  real(kind=rck)                     :: tmat(nbw,nbw,*)

  real(kind=rk)                               :: eps
  logical, intent(in)                         :: useGPU
  logical, intent(in)                         :: isSkewsymmetric
  character(20)                               :: gpuString

  integer(kind=ik)                            :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                      :: mpierr,  my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                            :: l_cols, l_rows
  integer(kind=ik)                            :: vmrCols
  integer(kind=ik)                            :: i, j, lcs, lce, lre, lc, lr, cur_pcol, n_cols, nrow
  integer(kind=ik)                            :: istep, ncol, lch, lcx, nlc
  integer(kind=ik)                            :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rk)                              :: vnorm2
  real(kind=rck)                    :: xf, aux1(nbw), aux2(nbw), vrl, tau
  real(kind=rck)                    :: vav(nbw,nbw)

  real(kind=rck), allocatable :: tmpGPU(:)
  real(kind=rck), pointer     :: vmrGPU(:), umcGPU(:)
  real(kind=rck), allocatable :: tmpCPU(:,:), vmrCPU(:,:), umcCPU(:,:)
  real(kind=rck), allocatable :: vr(:)

  ! needed for blocked QR decomposition
  integer(kind=ik)                            :: PQRPARAM(11), work_size
  real(kind=rk)                               :: dwork_size(1)
  real(kind=rk), allocatable                  :: work_blocked(:), tauvector(:), blockheuristic(:)
  integer(kind=C_intptr_T)                    :: a_dev, vmr_dev, umc_dev, tmat_dev, vav_dev
  type(c_ptr)                                 :: vmr_host, umc_host
  !integer(kind=ik), external                  :: numroc -> use elpa_scalapack
  integer(kind=ik)                            :: ierr
  integer(kind=ik)                            :: cur_l_rows, cur_l_cols, vmr_size, umc_size
  integer(kind=ik)                            :: l_rows2, vmr_size2, umc_size2
  integer(kind=c_intptr_t)                    :: lc_start, lc_end
  integer(kind=ik)                            :: lr_end
  integer(kind=ik)                            :: na_cols
  integer(kind=BLAS_KIND)                     :: na_colsBLAS

  logical, intent(in)                         :: wantDebug
  logical, intent(out)                        :: success
  logical                                     :: successGPU
  integer(kind=ik)                            :: istat
  character(200)                              :: errorMessage
  integer(kind=ik)                            :: min_tile_size, error

  logical, intent(in)                         :: useQR
  integer(kind=ik)                            :: mystart, myend, m_way, n_way, work_per_thread, m_id, n_id, n_threads, &
                                                ii, pp
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                    &single&
                                                                    &_&
                                                                    &real

  logical                                     :: useGPU_reduction_lower_block_to_tridiagonal
  integer(kind=ik), intent(in)                :: max_threads
  logical                                     :: do_memcpy
  integer(kind=ik)                            :: i_blk,blk_off, blk_end
  logical                                     :: useIntelGPU


  integer(kind=MPI_KIND)                      :: bcast_request, allreduce_request1, allreduce_request2, &
                                                 allreduce_request3, allreduce_request4, allreduce_request5, &
                                                 allreduce_request6

  logical                                     :: useNonBlockingCollectivesCols
  logical                                     :: useNonBlockingCollectivesRows
  integer(kind=c_int)                         :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("bandred_&
  &real&
  &" // &
  "_single" // &
  gpuString )

  call obj%get("nbc_row_elpa2_full_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_bandred. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_full_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_bandred. Aborting..."
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

  useGPU_reduction_lower_block_to_tridiagonal = .false.

  if (useGPU .and. .not.(useIntelGPU)) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
    if (useQR) then
      !in this case switch off GPU usage for step "reduce current block to lower triangular form"
      ! since this is done by QR decomposition
      useGPU_reduction_lower_block_to_tridiagonal = .false.
    endif
  endif

  if (useIntelGPU) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
    if (useQR) then
      !in this case switch off GPU usage for step "reduce current block to lower triangular form"
      ! since this is done by QR decomposition
      useGPU_reduction_lower_block_to_tridiagonal = .false.
    endif
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")
  success = .true.


  ! Semibandwith nbw must be a multiple of blocksize nblk
  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_bandred_&
                             &real&
                             &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_bandred_&
                             &real&
                             &: ELPA2 works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  ! na_rows in used nowhere; only na_cols
  if (useGPU .and. .not.(useIntelGPU)) then
    na_colsBLAS = numroc(int(na,kind=BLAS_KIND), int(nblk,kind=BLAS_KIND), &
                         int(my_pcol,kind=BLAS_KIND), 0_BLAS_KIND, int(np_cols,kind=BLAS_KIND))
    na_cols = int(na_colsBLAS,kind=c_int)

    ! Here we convert the regular host array into a pinned host array
    successGPU = gpu_malloc(a_dev, lda*na_cols* size_of_datatype)
    call check_alloc_GPU_f("bandred: a_dev", 335,  successGPU)

    successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), &
                  nbw * nbw * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: vav", 340,  successGPU)

    successGPU = gpu_malloc(vav_dev, nbw*nbw* size_of_datatype)
    call check_alloc_GPU_f("bandred: vav_dev", 343,  successGPU)
  endif ! useGPU

  !if (useIntelGPU) then
  !  ! needed later when explicit copy
  !endif ! useIntelGPU

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

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

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile

  if (useQR) then

    if (which_qr_decomposition == 1) then
      call qr_pqrparam_init(obj,pqrparam(1:11),    nblk,'M',0,   nblk,'M',0,   nblk,'M',1,'s')
      allocate(tauvector(na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: tauvector", 378,  istat,  errorMessage)

      allocate(blockheuristic(nblk), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: blockheuristic", 381,  istat,  errorMessage)

      l_rows = local_index(na, my_prow, np_rows, nblk, -1)
      allocate(vmrCPU(max(l_rows,1),na), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vmrCPU", 385,  istat,  errorMessage)

      vmrCols = na

      call qr_pdgeqrf_2dcomm_&
           &single&
           &(obj, a_mat(1:lda,1:matrixCols), matrixCols, lda, vmrCPU(1:max(l_rows,1),1:vmrCols), max(l_rows,1), &
                             vmrCols, tauvector(1:na), na, tmat(1:nbw,1:nbw,1), nbw, &
                             nbw, dwork_size(1:1), 1, -1, na, nbw, nblk, nblk, na, na, 1, 0, PQRPARAM(1:11), &
                             mpi_comm_rows, mpi_comm_cols, blockheuristic)

      work_size = int(dwork_size(1))
      allocate(work_blocked(work_size), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: work_blocked", 407,  istat,  errorMessage)
      work_blocked = 0.0_rk
      deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("bandred: vmrCPU", 410,  istat,  errorMessage)

    endif ! which_qr_decomposition

  endif ! useQr

  blk_end = (na-1)/nbw
  if (useGPU .and. .not.(useIntelGPU)) then

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: a_mat", 422,  successGPU)

    cur_l_rows = 0
    cur_l_cols = 0

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("bandred: a_dev", 429,  successGPU)

    successGPU = gpu_malloc(tmat_dev, nbw*nbw*size_of_datatype)
    call check_alloc_GPU_f("bandred: tmat_dev", 432,  successGPU)

    istep = (na-1)/nbw
    blk_end = (na-1)/nbw
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows,1)
    cur_l_cols = max(l_cols,1)
    vmr_size = cur_l_rows*2*n_cols
    umc_size = cur_l_cols*2*n_cols

    istep = (na-1)/nbw - 1
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows2 = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows2,1)
    cur_l_cols = max(l_cols,1)
    vmr_size2 = cur_l_rows*2*n_cols
    umc_size2 = cur_l_cols*2*n_cols

    l_rows = max(l_rows,l_rows2)
    vmr_size = max(vmr_size,vmr_size2)
    umc_size = max(umc_size,umc_size2)

    allocate(vr(l_rows + 1), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      print *,"bandred_&
              &real&
              &: error when allocating vr "//errorMessage
      stop 1
    endif

    successGPU = gpu_malloc_host(vmr_host,vmr_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: vmr_host", 466,  successGPU)
    call c_f_pointer(vmr_host, vmrGPU, (/vmr_size/))

    successGPU = gpu_malloc(vmr_dev, vmr_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: vmr_dev", 470,  successGPU)

    successGPU = gpu_malloc_host(umc_host,umc_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: umc_host", 473,  successGPU)
    call c_f_pointer(umc_host, umcGPU, (/umc_size/))

    successGPU = gpu_malloc(umc_dev, umc_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: umc_dev", 477,  successGPU)

  endif ! useGPU

  !if (useIntelGPU) then
     ! needed later when explict mem copy
  !endif ! useIntelGPU


  do istep = blk_end, 1, -1

    n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

    ! Number of local columns/rows of remaining matrix
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

    ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

    if (useGPU) then
      if (useIntelGPU) then
        ! unify the the name vmr and vmrCPU, as well as vmrGPU
        ! the same for umcCPU and umcGPU
        ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

        allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vmrCPU", 503,  istat,  errorMessage)

        allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: umcCPU", 506,  istat,  errorMessage)

        allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vr", 509,  istat,  errorMessage)
      else
        cur_l_rows = max(l_rows, 1)
        cur_l_cols = max(l_cols, 1)
        vmr_size = cur_l_rows * 2 * n_cols
        umc_size = cur_l_cols * 2 * n_cols
      endif
    else ! GPU not used

      ! unify the the name vmr and vmrCPU, as well as vmrGPU
      ! the same for umcCPU and umcGPU
      ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

      allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vmrCPU", 523,  istat,  errorMessage)

      allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: umcCPU", 526,  istat,  errorMessage)

      allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vr", 529,  istat,  errorMessage)

    endif ! use GPU

    if (useGPU) then
      if (useIntelGPU) then
        vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
      else
        vmrGPU(1 : cur_l_rows * n_cols) = 0.0_rck
        umcGPU(1 : umc_size) = 0.0_rck
      endif
    else
      vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
    endif ! useGPU

    vr(:) = 0.0_rck
    tmat(:,:,istep) = 0.0_rck
    if (useGPU .and. .not.(useIntelGPU)) then
      lc_start = local_index(istep*nbw+1, my_pcol, np_cols, nblk, -1)
      lc_end   = local_index(istep*nbw+n_cols, my_pcol, np_cols, nblk, -1)
      lr_end   = local_index((istep-1)*nbw + n_cols, my_prow, np_rows, nblk, -1)

      if (lc_start .le. 0) lc_start = 1

      do_memcpy = .false.

      ! Note: mod(nbw,nblk) == 0
      do i_blk = 1, nbw/nblk
        blk_off = (i_blk-1) * nblk
        cur_pcol = pcol(istep*nbw+1+blk_off, nblk, np_cols)

        if (my_pcol == cur_pcol) then
          do_memcpy = .true.
        endif
      enddo

      if (do_memcpy) then
        successGPU = gpu_memcpy2d(int(loc(a_mat(1, lc_start)),kind=c_intptr_t), &
                      int((lda*size_of_datatype),kind=c_intptr_t), &
                      (a_dev + int( ( (lc_start-1) * lda*size_of_datatype),kind=c_intptr_t )), &
                      int(lda*size_of_datatype,kind=c_intptr_t), &
                      int(lr_end*size_of_datatype,kind=c_intptr_t), &
                      int((lc_end - lc_start+1),kind=c_intptr_t),int(gpuMemcpyDeviceToHost,kind=c_int))

        call check_memcpy_GPU_f("bandred: a_dev -> a_mat", 573,  successGPU)
      endif
    endif ! useGPU

    !if (useIntelGPU) then
       ! needed later when explict mem copy
    !endif ! useIntelGPU


    ! Reduce current block to lower triangular form
    if (useQR) then
      if (which_qr_decomposition == 1) then
        vmrCols = 2*n_cols
        call qr_pdgeqrf_2dcomm_&
             &single&
             &(obj, a_mat(1:lda,1:matrixCols), lda, matrixCols, vmrCPU(1:max(l_rows,1),1:vmrCols) ,   &
                                max(l_rows,1), vmrCols, tauvector(1:na), na, &
                                 tmat(1:nbw,1:nbw,istep), nbw, nbw, work_blocked(1:work_size), work_size, &
                                 work_size, na, n_cols, nblk, nblk,        &
                                 istep*nbw+n_cols-nbw, istep*nbw+n_cols, 1,&
                                 0, PQRPARAM(1:11), mpi_comm_rows, mpi_comm_cols,&
                                 blockheuristic)
      endif

    else !useQR
      do lc = n_cols, 1, -1

        ncol = istep*nbw + lc ! absolute column number of householder Vector
        nrow = ncol - nbw ! Absolute number of pivot row

        lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
        lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

        tau = 0

        if (nrow == 1) exit ! Nothing to do

        cur_pcol = pcol(ncol, nblk, np_cols) ! Processor column owning current block

        if (my_pcol==cur_pcol) then

          ! Get Vector to be transformed; distribute last element and norm of
          ! remaining elements to all procs in current column

          vr(1:lr) = a_mat(1:lr,lch) ! Vector to be transformed

          if (my_prow==prow(nrow, nblk, np_rows)) then
            aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
            aux1(2) = vr(lr)
          else
            aux1(1) = dot_product(vr(1:lr),vr(1:lr))
            aux1(2) = 0.0_rck
          endif

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL4, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             allreduce_request1, mpierr)

            call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)

            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_REAL4, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             mpierr)

            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif


          vnorm2 = aux2(1)
          vrl    = aux2(2)

          ! Householder transformation
          call hh_transform_&
             &real&
             &_&
             &single &
                         (obj, vrl, vnorm2, xf, tau, wantDebug)
          ! Scale vr and store Householder Vector for back transformation

          vr(1:lr) = vr(1:lr) * xf
          if (my_prow==prow(nrow, nblk, np_rows)) then
            a_mat(1:lr-1,lch) = vr(1:lr-1)
            a_mat(lr,lch) = vrl
            vr(lr) = 1.0_rck
          else
            a_mat(1:lr,lch) = vr(1:lr)
          endif

        endif

        ! Broadcast Householder Vector and tau along columns

        vr(lr+1) = tau
        if (useNonBlockingCollectivesCols) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(vr, int(lr+1,kind=MPI_KIND), MPI_REAL4, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        bcast_request, mpierr)


          call mpi_wait(bcast_request, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_bcast(vr, int(lr+1,kind=MPI_KIND), MPI_REAL4, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
          vmrGPU(cur_l_rows * (lc - 1) + 1 : cur_l_rows * (lc - 1) + lr) = vr(1:lr)
        else
          vmrCPU(1:lr,lc) = vr(1:lr)
        endif
        tau = vr(lr+1)

        tmat(lc,lc,istep) = tau ! Store tau in diagonal of tmat
        ! Transform remaining columns in current block with Householder Vector
        ! Local dot product

        aux1 = 0.0_rck


        nlc = 0 ! number of local columns
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            if (lr>0) aux1(nlc) = dot_product(vr(1:lr),a_mat(1:lr,lcx))
          endif
        enddo

        ! Get global dot products
        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          if (nlc>0) then
            call mpi_iallreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_REAL4, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      allreduce_request3, mpierr)
            call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          if (nlc>0) then
            call mpi_allreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_REAL4, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif
        ! Transform

        nlc = 0
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            a_mat(1:lr,lcx) = a_mat(1:lr,lcx) - tau*aux2(nlc)*vr(1:lr)
          endif
        enddo
      enddo ! lc

      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        ! store column tiles back to GPU
        if (do_memcpy) then
          successGPU = gpu_memcpy2d((a_dev+ &
                        int(((lc_start-1)*lda*size_of_datatype),kind=c_intptr_t)), &
                        int(lda*size_of_datatype,kind=c_intptr_t), int(loc(a_mat(1,lc_start)),kind=c_intptr_t), &
                        int(lda*size_of_datatype,kind=c_intptr_t), &
                        int(lr_end*size_of_datatype,kind=c_intptr_t), &
                        int((lc_end - lc_start+1),kind=c_intptr_t), &
                        int(gpuMemcpyHostToDevice,kind=c_int))
          call check_memcpy_GPU_f("bandred: a_mat -> a_dev", 937,  successGPU)
        endif
      endif

      if (useGPU_reduction_lower_block_to_tridiagonal .and. useIntelGPU) then
        ! store column tiles back to GPU
      endif

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use dsyrk

      vav = 0
      call obj%timer%start("blas")
      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        if (l_rows>0) &
        call SSYRK('U', 'T',            &
                           int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                           vmrGPU, int(cur_l_rows,kind=BLAS_KIND), &
                           ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      else ! useGPU_reduction_to_tridiagonal
        if (l_rows>0) &
        call SSYRK('U', 'T',           &
                            int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, vmrCPU, &
                            int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      call symm_matrix_allreduce_&
         &single &
                         (obj, n_cols,vav, nbw, nbw,mpi_comm_rows, .true.)
         ! Calculate triangular matrix T for block Householder Transformation
      call obj%timer%start("blas")
      do lc=n_cols,1,-1
        tau = tmat(lc,lc,istep)
        if (lc<n_cols) then
          call STRMV('U', 'T', 'N',&
                              int(n_cols-lc,kind=BLAS_KIND), tmat(lc+1,lc+1,istep), &
                              int(ubound(tmat,dim=1),kind=BLAS_KIND), vav(lc+1,lc), 1_BLAS_KIND)

          tmat(lc,lc+1:n_cols,istep) = -tau * vav(lc+1:n_cols,lc)
        endif
      enddo
      call obj%timer%stop("blas")
    endif !useQR

    if (useGPU .and. useQR ) then
      if (useIntelGPU) then
        ! copy the data for furhter usage
        ! qr worked on *CPU arrarys
        !vmrCUDA(1:cur_l_rows * n_cols) = vmrCPU(1:cur_l_rows,1:n_cols)
      else     
        ! copy the data for furhter usage
        ! qr worked on *CPU arrarys
        !vmrGPU(1:cur_l_rows * n_cols) = vmrCPU(1:cur_l_rows,1:n_cols)
        if (do_memcpy) then
          successGPU = gpu_memcpy2d((a_dev+ &
                        int(((lc_start-1)*lda*size_of_datatype),kind=c_intptr_t)), &
                        int(lda*size_of_datatype,kind=c_intptr_t), int(loc(a_mat(1,lc_start)),kind=c_intptr_t), &
                        int(lda*size_of_datatype,kind=c_intptr_t), &
                        int(lr_end*size_of_datatype,kind=c_intptr_t), &
                        int((lc_end - lc_start+1),kind=c_intptr_t), &
                        int(gpuMemcpyHostToDevice,kind=c_int))
          call check_memcpy_GPU_f("bandred: a_mat -> a_dev", 1022,  successGPU)
        endif
      endif
    endif

    ! Transpose vmr -> vmc (stored in umc, second half)
    if (useGPU) then
      if (useIntelGPU) then
        call elpa_transpose_vectors_&
             &real&
             &_&
             &single &
                                          (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                           umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                           1, istep*nbw, n_cols, nblk, max_threads, .true.)

      else
        call elpa_transpose_vectors_&
             &real&
             &_&
             &single &
                          (obj, vmrGPU(:), cur_l_rows, mpi_comm_rows, &
                           umcGPU(cur_l_cols * n_cols + 1:), cur_l_cols, &
                           mpi_comm_cols, 1, istep*nbw, n_cols, nblk, max_threads, .true.)
      endif
    else ! useGPU
      call elpa_transpose_vectors_&
           &real&
           &_&
           &single &
                                        (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                         umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                         1, istep*nbw, n_cols, nblk, max_threads, .true.)
    endif

    ! Calculate umc = A**T * vmr
    ! Note that the distributed A has to be transposed
    ! Opposed to direct tridiagonalization there is no need to use the cache locality
    ! of the tiles, so we can use strips of the matrix


!#if 0
!    ! original complex implemetation check for performance
!    umcCPU(1:l_cols,1:n_cols) = 0.0_rck
!    vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
!
!    if (l_cols>0 .and. l_rows>0) then
!      do i=0,(istep*nbw-1)/tile_size
!
!        lcs = i*l_cols_tile+1
!        lce = min(l_cols,(i+1)*l_cols_tile)
!        if (lce<lcs) cycle
!
!        lre = min(l_rows,(i+1)*l_rows_tile)
!
!          call obj%timer%start("blas")
!          call SGEMM('C', 'N', lce-lcs+1, n_cols, lre, ONE, a_mat(1,lcs), ubound(a_mat,dim=1), &
!                     vmrCPU, ubound(vmrCPU,dim=1), ONE, umcCPU(lcs,1), ubound(umcCPU,dim=1))
!          call obj%timer%stop("blas")
!
!        if (i==0) cycle
!        lre = min(l_rows,i*l_rows_tile)
!          call obj%timer%start("blas")
!          call SGEMM('N', 'N', lre, n_cols, lce-lcs+1, ONE, a_mat(1,lcs), lda, &
!                     umcCPU(lcs,n_cols+1), ubound(umcCPU,dim=1), ONE, vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1))
!          call obj%timer%stop("blas")
!      enddo
!
!    endif ! (l_cols>0 .and. l_rows>0)
!#endif 

    !Code for Algorithm 4

    ! n_way is actually a branch for the number of OpenMP threads
    n_way = 1

      if (.not. useGPU .or. useIntelGPU) then
        umcCPU(1:l_cols,1:n_cols) = 0.0_rck
        vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
      endif ! useGPU

      if (l_cols>0 .and. l_rows>0) then

        if (useGPU .and. .not.(useIntelGPU)) then
          successGPU = gpu_memset(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                        0, cur_l_rows*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: vmr_dev", 1193,  successGPU)

          successGPU = gpu_memcpy(vmr_dev, int(loc(vmrGPU(1)),kind=c_intptr_t), &
                        cur_l_rows*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: vmrGPU -> vmr_dev", 1197,  successGPU)

          successGPU = gpu_memset(umc_dev, 0, l_cols*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: umc_dev", 1200,  successGPU)

          successGPU = gpu_memcpy(umc_dev+l_cols*n_cols*size_of_datatype, &
                        int(loc(umcGPU(1+l_cols*n_cols)),kind=c_intptr_t), &
                        (umc_size-l_cols*n_cols)*size_of_datatype, &
                        gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev", 1206,  successGPU)
        endif ! useGPU

        do i=0,(istep*nbw-1)/tile_size

          lcs = i*l_cols_tile+1
          lce = min(l_cols,(i+1)*l_cols_tile)
          if (lce<lcs) cycle
          lre = min(l_rows,(i+1)*l_rows_tile)

          if (useGPU) then
            if (useIntelGPU) then
              call obj%timer%start("mkl_offload")
              call obj%timer%stop("mkl_offload")
              if (i==0) cycle
              lre = min(l_rows,i*l_rows_tile)
              call obj%timer%start("mkl_offload")

              if (isSkewsymmetric) then

              else
              endif
              call obj%timer%stop("mkl_offload")

            else
              call obj%timer%start("gpublas")
              call gpublas_SGEMM('T', 'N',                   &
                                         lce-lcs+1, n_cols, lre,     &
                                         ONE, (a_dev + ((lcs-1)*lda* &
                                         size_of_datatype)),         &
                                         lda, vmr_dev,cur_l_rows,    &
                                         ONE, (umc_dev+ (lcs-1)*     &
                                             size_of_datatype),      &
                                         cur_l_cols)

              call obj%timer%stop("gpublas")

              if(i==0) cycle
              call obj%timer%start("gpublas")

              lre = min(l_rows,i*l_rows_tile)
              if (isSkewsymmetric) then
                call gpublas_SGEMM('N', 'N', lre,n_cols, lce-lcs+1, -ONE, &
                              (a_dev+ ((lcs-1)*lda*                 &
                                    size_of_datatype)),             &
                         lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                size_of_datatype),              &
                                cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                              size_of_datatype),              &
                                cur_l_rows)
              else
                call gpublas_SGEMM('N', 'N', lre,n_cols, lce-lcs+1, ONE, &
                                            (a_dev+ ((lcs-1)*lda*                 &
                                                  size_of_datatype)),             &
                                       lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                              size_of_datatype),              &
                                              cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                                            size_of_datatype),              &
                                              cur_l_rows)
              endif
              call obj%timer%stop("gpublas")
            endif
          else ! useGPU

            call obj%timer%start("blas")
            call SGEMM('T', 'N',       &
                                int(lce-lcs+1,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lre,kind=BLAS_KIND), &
                                ONE, a_mat(1,lcs), int(ubound(a_mat,dim=1),kind=BLAS_KIND), &
                                vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ONE, umcCPU(lcs,1), &
                                int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
            call obj%timer%stop("blas")
            if (i==0) cycle
            lre = min(l_rows,i*l_rows_tile)
            call obj%timer%start("blas")

            if (isSkewsymmetric) then
              call SGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  -ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )

            else
              call SGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )
            endif
            call obj%timer%stop("blas")
          endif ! useGPU
        enddo ! i=0,(istep*nbw-1)/tile_size

        if (useGPU .and. .not.(useIntelGPU)) then
          if (tile_size < istep*nbw .or. n_way > 1) then
            successGPU = gpu_memcpy(int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                          vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                          (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyDeviceToHost)
            call check_memcpy_GPU_f("bandred: vmr_dev -> vmrGPU", 1342,  successGPU)
          endif

          successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                        umc_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU", 1347,  successGPU)
        endif ! useGPU
      endif ! l_cols>0 .and. l_rows>0

    ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if ur has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix

    ! Or if we used the Algorithm 4
    if (tile_size < istep*nbw .or. n_way > 1) then

      if (useGPU) then
        if (useIntelGPU) then
          call elpa_reduce_add_vectors_&
          &real&
          &_&
          &single &
                                      (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                       umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                      istep*nbw, n_cols, nblk, max_threads)

        else

          call elpa_reduce_add_vectors_&
               &real&
               &_&
               &single &
                               (obj, vmrGPU(cur_l_rows * n_cols + 1:),cur_l_rows,  &
                                mpi_comm_rows, umcGPU,                            &
                                cur_l_cols, mpi_comm_cols, istep*nbw, n_cols, nblk, max_threads)
        endif
      else ! useGPU

        call elpa_reduce_add_vectors_&
        &real&
        &_&
        &single &
                                         (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                          umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          istep*nbw, n_cols, nblk, max_threads)
      endif ! useGPU
    endif ! tile_size < istep*nbw .or. n_way > 1

    if (l_cols>0) then

      if (useGPU) then
        if (useIntelGPU) then
          allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpCPU", 1399,  istat,  errorMessage)

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL4,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           allreduce_request4, mpierr)
            call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL4,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("bandred: tmpCPU", 1421,  istat,  errorMessage)
        else
          allocate(tmpGPU(l_cols * n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpGPU", 1425,  istat,  errorMessage)
          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")

            call mpi_iallreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL4, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request5, mpierr)
            call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")

            call mpi_allreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL4, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          if (allocated(tmpGPU)) then
            deallocate(tmpGPU, stat=istat, errmsg=errorMessage)
            call check_deallocate_f("bandred: tmpGPU", 1447,  istat,  errorMessage)
          endif
        endif

      else ! useGPU

        allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: tmpCPU", 1454,  istat,  errorMessage)

        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL4,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6, mpierr)
          call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_REAL4,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: tmpCPU", 1474,  istat,  errorMessage)
      endif ! useGPU
    endif ! l_cols > 0

    ! U = U * Tmat**T

    if (useGPU) then
      if (useIntelGPU) then
        call obj%timer%start("mkl_offload")

        call STRMM('Right', 'Upper', 'T', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

        call SGEMM('T', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

        call STRMM('Right', 'Upper', 'T', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
        call obj%timer%stop("mkl_offload")

      else
        successGPU = gpu_memcpy(umc_dev, int(loc(umcGPU(1)),kind=c_intptr_t), &
                      l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev ", 1526,  successGPU)

        successGPU = gpu_memcpy(tmat_dev,int(loc(tmat(1,1,istep)),kind=c_intptr_t), &
                      nbw*nbw*size_of_datatype,gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: tmat -> tmat_dev ", 1530,  successGPU)

        call obj%timer%start("gpublas")
        call gpublas_STRMM('Right', 'Upper', 'T', 'Nonunit',  &
                            l_cols, n_cols, ONE, tmat_dev, nbw, umc_dev, cur_l_cols)
        call obj%timer%stop("gpublas")

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T
        call obj%timer%start("gpublas")
        call gpublas_SGEMM('T', 'N',             &
                                 n_cols, n_cols, l_cols, ONE, umc_dev, cur_l_cols, &
                                 (umc_dev+(cur_l_cols * n_cols )*size_of_datatype),cur_l_cols, &
                                 ZERO, vav_dev, nbw)

        call gpublas_STRMM('Right', 'Upper', 'T', 'Nonunit',    &
           n_cols, n_cols, ONE, tmat_dev, nbw, vav_dev, nbw)
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(vav),kind=c_intptr_t), &
                    vav_dev, nbw*nbw*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: vav_dev -> vav ", 1550,  successGPU)
      endif
    else ! useGPU

      call obj%timer%start("blas")

      call STRMM('Right', 'Upper', 'T', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      call SGEMM('T', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      call STRMM('Right', 'Upper', 'T', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
      call obj%timer%stop("blas")

    endif ! useGPU

!#ifdef 1
    if (isSkewsymmetric) then
      call ssymm_matrix_allreduce_&
      &single &
      (obj, n_cols,vav, nbw, nbw ,mpi_comm_cols, .false.)
    else
!#endif
      call symm_matrix_allreduce_&
      &single &
      (obj, n_cols,vav, nbw, nbw ,mpi_comm_cols, .false.)
!#ifdef 1
    endif
!#endif

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(vav_dev, int(loc(vav),kind=c_intptr_t), &
                       nbw*nbw*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("bandred: vav -> vav_dev ", 1599,  successGPU)
    endif
    !if (useIntelGPU) then
      ! needed later
    !endif


    ! U = U - 0.5 * V * VAV

    if (useGPU) then
      if (useIntelGPU) then
      call obj%timer%start("mkl_offload")
      if (isSkewsymmetric) then
        call SGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                        &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      else
        call SGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            -0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                       &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      endif

      call obj%timer%stop("mkl_offload")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &real&
        &_&
        &single &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &real&
       &_&
       &single &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif

      else
        call obj%timer%start("gpublas")
        if (isSkewsymmetric) then
          call gpublas_SGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                    0.5_rk,                      &
                                    (umc_dev+(cur_l_cols * n_cols )* &
                                    size_of_datatype),   &
                                    cur_l_cols, vav_dev,nbw,        &
                                    ONE, umc_dev, cur_l_cols)
        else
          call gpublas_SGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                   -0.5_rk,                      &
                                   (umc_dev+(cur_l_cols * n_cols )* &
                                   size_of_datatype),   &
                                   cur_l_cols, vav_dev,nbw,        &
                                   ONE, umc_dev, cur_l_cols)
        endif
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                    umc_dev, umc_size*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU ", 1681,  successGPU)

        ! Transpose umc -> umr (stored in vmr, second half)
        if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
             &real&
             &_&
             &single &
                         (obj, umcGPU(:), cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        else
          call elpa_transpose_vectors_&
             &real&
             &_&
             &single &
                         (obj, umcGPU, cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        endif

        successGPU = gpu_memcpy(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                    int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                    (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: vmr -> vmrGPU ", 1705,  successGPU)
      endif

    else ! useGPU
      call obj%timer%start("blas")
      if (isSkewsymmetric) then
        call SGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                        &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      else
        call SGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                            -0.5_rk, umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav,                       &
                            int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
      endif

      call obj%timer%stop("blas")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &real&
        &_&
        &single &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &real&
       &_&
       &single &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif
    endif  ! useGPU

    ! A = A - V*U**T - U*V**T


    do i=0,(istep*nbw-1)/tile_size
      lcs = i*l_cols_tile+1
      lce = min(l_cols,(i+1)*l_cols_tile)
      lre = min(l_rows,(i+1)*l_rows_tile)
      if (lce<lcs .or. lre<1) cycle

      if (useGPU) then
        if (useIntelGPU) then
          call obj%timer%start("mkl_offload")
          call SGEMM('N', 'T', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                              int(2*n_cols,kind=BLAS_KIND), &
                              -ONE, &
                              vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                              int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                              ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
          call obj%timer%stop("mkl_offload")

        else
          call obj%timer%start("gpublas")

          call gpublas_SGEMM('N', 'T',     &
                                     lre, lce-lcs+1, 2*n_cols, -ONE, &
                                     vmr_dev, cur_l_rows, (umc_dev +(lcs-1)*  &
                                     size_of_datatype), &
                                     cur_l_cols, ONE, (a_dev+(lcs-1)*lda* &
                                     size_of_datatype), lda)
          call obj%timer%stop("gpublas")
        endif
      else ! useGPU

        call obj%timer%start("blas")
        call SGEMM('N', 'T', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                            int(2*n_cols,kind=BLAS_KIND), &
                            -ONE, &
                            vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                            int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                            ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    enddo ! i=0,(istep*nbw-1)/tile_size

    if (.not.(useGPU) .or. useIntelGPU) then
      if (allocated(vr)) then
        deallocate(vr, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vr", 1871,  istat,  errorMessage)
      endif

      if (allocated(umcCPU)) then
        deallocate(umcCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: umcCPU", 1876,  istat,  errorMessage)
      endif

      if (allocated(vmrCPU)) then
        deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vmrCPU", 1881,  istat,  errorMessage)
      endif
    endif !useGPU

  enddo ! istep - loop

  if (useGPU .and. .not.(useIntelGPU)) then
    ! copy a_dev to a_mat
    ! we do it here, since a is needed on the host in the following routine
    ! (band to tridi). Previously, a has been kept on the device and then
    ! copied in redist_band (called from tridiag_band). However, it seems to
    ! be easier to do it here.
    successGPU = gpu_memcpy(int(loc(a_mat),kind=c_intptr_t), &
                  int(a_dev,kind=c_intptr_t), &
                  int(lda*matrixCols* size_of_datatype, kind=c_intptr_t), &
                  gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("bandred: a_dev -> a_mat ", 1897,  successGPU)

    successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: a_mat ", 1900,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("bandred: a_dev ", 1903,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("bandred: vav_dev ", 1906,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("bandred: tmat_dev ", 1909,  successGPU)

    successGPU = gpu_host_unregister(int(loc(vav),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: vav", 1912,  successGPU)

    if (associated(umcGPU)) then
      nullify(umcGPU)

      successGPU = gpu_free_host(umc_host)
      call check_host_dealloc_GPU_f("bandred: umc_host ", 1918,  successGPU)

      successGPU = gpu_free(umc_dev)
      call check_dealloc_GPU_f("bandred: umc_dev ", 1921,  successGPU)
    endif

    if (associated(vmrGPU)) then
      nullify(vmrGPU)

      successGPU = gpu_free_host(vmr_host)
      call check_host_dealloc_GPU_f("bandred: vmr_host ", 1928,  successGPU)

      successGPU = gpu_free(vmr_dev)
      call check_dealloc_GPU_f("bandred: vmr_dev ", 1931,  successGPU)
    endif
  endif ! useGPU

  !if (useIntelGPU) then
  !   ! needed later
  !endif

  if (allocated(vr)) then
    deallocate(vr, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vr", 1941,  istat,  errorMessage)
  endif

  if (allocated(umcCPU)) then
    deallocate(umcCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: umcCPU", 1946,  istat,  errorMessage)
  endif

  if (allocated(vmrCPU)) then
    deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vmrCPU", 1951,  istat,  errorMessage)
  endif

  if (useQR) then
    if (which_qr_decomposition == 1) then
      deallocate(work_blocked, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("bandred: work_blocked", 1958,  istat,  errorMessage)

      deallocate(tauvector, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("bandred: tauvector", 1961,  istat,  errorMessage)
    endif
  endif

  call obj%timer%stop("bandred_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

end subroutine bandred_&
&real&
&_&
&single

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







subroutine symm_matrix_allreduce_&
&single &
                    (obj, n, a, lda, ldb, comm, isRows)
!-------------------------------------------------------------------------------
!  symm_matrix_allreduce: Does an mpi_allreduce for a symmetric matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use precision
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)             :: n, lda, ldb, comm
  real(kind=rk4)     :: a(lda,*)
  integer(kind=ik)             :: i, nc
  integer(kind=MPI_KIND)       :: mpierr
  real(kind=rk4)     :: h1(n*n), h2(n*n)
  integer(kind=MPI_KIND)       :: allreduce_request1
  logical                      :: useNonBlockingCollectives
  logical                      :: useNonBlockingCollectivesCols
  logical                      :: useNonBlockingCollectivesRows
  logical, intent(in)          :: isRows
  integer(kind=c_int)          :: non_blocking_collectives_rows, error, &
                                  non_blocking_collectives_cols

  call obj%timer%start("&
          &symm_matrix_allreduce&
          &" // &
          &"_single"&
          )

  call obj%get("nbc_row_sym_allreduce", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa_sym_allreduce. Aborting..."
    stop
  endif

  call obj%get("nbc_col_sym_allreduce", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa_sym_allreduce. Aborting..."
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

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  nc = 0
  do i=1,n
    h1(nc+1:nc+i) = a(1:i,i)
    nc = nc+i
  enddo

  if (useNonBlockingCollectives) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                     int(comm,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                     int(comm,kind=MPI_KIND), mpierr)
    call obj%timer%stop("mpi_communication")
  endif
  nc = 0
  do i=1,n
    a(1:i,i) = h2(nc+1:nc+i)
    a(i,1:i-1) = a(1:i-1,i)
    nc = nc+i
  enddo

! nc = 0
! do i=1,n
!   a(1:i,i) = h2(nc+1:nc+i)
!   a(i,1:i-1) = a(1:i-1,i)
!   nc = nc+i
! enddo

  call obj%timer%stop("&
          &symm_matrix_allreduce&
          &" // &
          &"_single"&
          )

end subroutine symm_matrix_allreduce_&
&single



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







subroutine ssymm_matrix_allreduce_&
&single &
                    (obj, n, a, lda, ldb, comm, isRows)
!-------------------------------------------------------------------------------
!  symm_matrix_allreduce: Does an mpi_allreduce for a symmetric matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use precision
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)             :: n, lda, ldb, comm
  real(kind=rk4)     :: a(lda,*)
  integer(kind=ik)             :: i, nc
  integer(kind=MPI_KIND)       :: mpierr
  real(kind=rk4)     :: h1(n*n), h2(n*n)
  integer(kind=MPI_KIND)       :: allreduce_request1
  logical                      :: useNonBlockingCollectives
  logical                      :: useNonBlockingCollectivesCols
  logical                      :: useNonBlockingCollectivesRows
  logical, intent(in)          :: isRows
  integer(kind=c_int)          :: non_blocking_collectives_rows, error, &
                                  non_blocking_collectives_cols

  call obj%timer%start("&
          &ssymm_matrix_allreduce&
          &" // &
          &"_single"&
          )

  call obj%get("nbc_row_sym_allreduce", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa_sym_allreduce. Aborting..."
    stop
  endif

  call obj%get("nbc_col_sym_allreduce", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa_sym_allreduce. Aborting..."
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

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  nc = 0
  do i=1,n
    h1(nc+1:nc+i) = a(1:i,i)
    nc = nc+i
  enddo

  if (useNonBlockingCollectives) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                     int(comm,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                     int(comm,kind=MPI_KIND), mpierr)
    call obj%timer%stop("mpi_communication")
  endif
  nc = 0
  do i=1,n
    a(1:i,i) = h2(nc+1:nc+i)
    a(i,1:i-1) = - a(1:i-1,i)
    nc = nc+i
  enddo

! nc = 0
! do i=1,n
!   a(1:i,i) = h2(nc+1:nc+i)
!   a(i,1:i-1) = a(1:i-1,i)
!   nc = nc+i
! enddo

  call obj%timer%stop("&
          &ssymm_matrix_allreduce&
          &" // &
          &"_single"&
          )

end subroutine ssymm_matrix_allreduce_&
&single








subroutine trans_ev_band_to_full_&
    &real&
    &_&
    &single &
    (obj, na, nqc, nblk, nbw, a_mat, lda, tmat, q_mat, &
     ldq, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, useGPU &
     ,useQr)

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_real/complex:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a_mat, number of rows of matrix q_mat
!
!  nqc         Number of columns of matrix q_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a_mat(lda,matrixCols)    Matrix containing the Householder vectors (i.e. matrix a_mat after bandred_real/complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat and q_mat
!
!  tmat(nbw,nbw,numBlocks) Factors returned by bandred_real/complex
!
!  q_mat           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q_mat
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa_abstract_impl
  use elpa_blas_interfaces

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
  logical, intent(in)                    :: useGPU
  logical, intent(in)                     :: useQR
  integer(kind=ik)                       :: na, nqc, lda, ldq, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols
  real(kind=rck)                :: a_mat(lda,*)
  real(kind=rck)                :: q_mat(ldq,*), tmat(nbw,nbw,*)

  integer(kind=ik)                       :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                 :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, mpierr
  integer(kind=ik)                       :: max_blocks_row, max_blocks_col, max_local_rows, &
                                            max_local_cols
  integer(kind=ik)                       :: l_cols, l_rows, l_colh, n_cols
  integer(kind=ik)                       :: istep, lc, ncol, nrow, nb, ns

  real(kind=rck), allocatable   :: hvb(:)
  real(kind=rck), pointer       :: hvm(:,:), tmp1(:), tmp2(:)
  ! hvm_dev is fist used and set in this routine
  ! q_mat is changed in trans_ev_tridi on the host, copied to device and passed here. this can be adapted
  ! tmp_dev is first used in this routine
  ! tmat_dev is not passed along from bandred_real
  integer(kind=C_intptr_T)               :: hvm_dev, q_dev, tmp_dev, tmat_dev, t_tmp_dev, dev_offset
  type(c_ptr)                            :: hvm_host, tmp1_host, tmp2_host
  type(c_ptr)                            :: tmat_mpi_dev, t_tmp_mpi_dev
  real(kind=rck), pointer       :: tmat_mpi_deviceptr(:,:)
  real(kind=rck), pointer       :: t_tmp_mpi_deviceptr(:)

  type(c_ptr)                            :: hvm_mpi_dev
  real(kind=rck), pointer       :: hvm_mpi_deviceptr(:,:)


  integer(kind=ik)                       :: i

  real(kind=rck), allocatable, target   :: tmat_complete(:,:), t_tmp(:,:), t_tmp2(:,:)
  integer(kind=ik)                       :: t_cols, t_rows, ii, jj
  integer(kind=ik)                       :: cwy_blocking

  integer(kind=ik)                       :: istat
  character(200)                         :: errorMessage
  character(20)                          :: gpuString
  logical                                :: successGPU
  integer(kind=c_intptr_t), parameter    :: size_of_datatype = size_of_&
                                                               &single&
                                                               &_&
                                                               &real
  integer(kind=ik)                       :: blocking_factor, error, blk_end
  logical                                :: useIntelGPU
  integer(kind=MPI_KIND)                 :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                :: useNonBlockingCollectivesCols
  logical                                :: useNonBlockingCollectivesRows
  integer(kind=c_int)                    :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_band_to_full_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa2_band_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_band_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_band_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_band_to_full. Aborting..."
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

  call obj%get("blocking_in_band_to_full",blocking_factor,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem getting option for blocking_in_band_to_full. Aborting..."
    stop
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  max_blocks_row = ((na -1)/nblk)/np_rows + 1 ! Rows of a_mat
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1 ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  cwy_blocking = blocking_factor * nbw

  if (useGPU) then
    if (useIntelGPU) then
      allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp1", 249,  istat,  errorMessage)
  
      allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp2", 252,  istat,  errorMessage)
  
      allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: hvm", 255,  istat,  errorMessage)

    else
      ! copy q_mat to q_dev
      successGPU = gpu_malloc(q_dev,ldq*matrixCols*size_of_datatype)
      call check_alloc_GPU_f("trans_ev_band_to_full: q_dev", 260,  successGPU)

      successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuHostRegisterDefault)
      call check_host_register_GPU_f("trans_ev_band_to_full: q_mat", 264,  successGPU)

      successGPU = gpu_memcpy(q_dev,int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_mat -> q_dev", 268,  successGPU)

      successGPU = gpu_malloc_host(tmp1_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp1_host", 271,  successGPU)
      call c_f_pointer(tmp1_host, tmp1, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(tmp2_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp2_host", 275,  successGPU)
      call c_f_pointer(tmp2_host, tmp2, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(hvm_host,max_local_rows*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: hvm_host", 279,  successGPU)
      call c_f_pointer(hvm_host, hvm, (/max_local_rows,cwy_blocking/))
    endif
  else ! useGPU
    allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp1", 284,  istat,  errorMessage)

    allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp2", 287,  istat,  errorMessage)

    allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: hvm", 290,  istat,  errorMessage)
  endif !useGPU

  allocate(hvb(max_local_rows*cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: hvb", 294,  istat,  errorMessage)

  allocate(tmat_complete(cwy_blocking,cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: tmat_complete", 297,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_host_register(int(loc(tmat_complete),kind=c_intptr_t), &
                  cwy_blocking * cwy_blocking * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_band_to_full: tmat_complete", 303,  successGPU)
  endif
  !if (useIntelGPU) then
  !  ! needed later
  !endif



  if (blocking_factor > 1) then
    allocate(t_tmp(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp", 322,  istat,  errorMessage)

    allocate(t_tmp2(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp2", 325,  istat,  errorMessage)
  endif

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_malloc(hvm_dev,max_local_rows*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: hvm_dev", 330,  successGPU)

    successGPU = gpu_malloc(tmp_dev,max_local_cols*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmp_dev", 333,  successGPU)

    successGPU = gpu_malloc(tmat_dev,cwy_blocking*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmat_dev", 336,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  hvm = 0.0_rck ! Must be set to 0 !!!
  hvb = 0.0_rck ! Safety only
  tmp1 = 0.0_rck
  tmp2 = 0.0_rck
  tmat_complete = 0.0_rck
  if (blocking_factor > 1) then
     t_tmp = 0.0_rck ! Must be set to 0 !!!
     t_tmp2 = 0.0_rck

  endif
  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  blk_end = ((na-1)/nbw-1)/blocking_factor + 1
  do istep=1, blk_end

    ! This the call when using na >= ((blocking_factor+1)*nbw)
    ! n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw
    ! Number of columns in current step
    ! As an alternative we add some special case handling if na < cwy_blocking
    if (na < cwy_blocking) then
      n_cols = MAX(0, na-nbw)
      if ( n_cols .eq. 0 ) then
        exit
      end if
    else
      n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw ! Number of columns in current step
    end if

    ! Broadcast all Householder vectors for current step compressed in hvb

    nb = 0
    ns = 0

    do lc = 1, n_cols
      ncol = (istep-1)*cwy_blocking + nbw + lc ! absolute column number of householder Vector
      nrow = ncol - nbw ! absolute number of pivot row

      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
      l_colh = local_index(ncol , my_pcol, np_cols, nblk, -1) ! HV local column number

      if (my_pcol==pcol(ncol, nblk, np_cols)) hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)

      nb = nb+l_rows

      if (lc==n_cols .or. mod(ncol,nblk)==0) then
        if (useNonBlockingCollectivesCols) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_REAL4,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_REAL4,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         mpierr)
          call obj%timer%stop("mpi_communication")
        endif
        ns = nb
      endif
    enddo ! lc

    ! Expand compressed Householder vectors into matrix hvm

    nb = 0
    do lc = 1, n_cols
      nrow = (istep-1)*cwy_blocking + lc ! absolute number of pivot row
      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

      hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
      if (my_prow==prow(nrow, nblk, np_rows)) hvm(l_rows+1,lc) = 1.0_rck
      nb = nb+l_rows
    enddo

    l_rows = local_index(MIN(na,(istep+1)*cwy_blocking), my_prow, np_rows, nblk, -1)

    ! compute tmat2 out of tmat(:,:,)
    tmat_complete = 0
    do i = 1, blocking_factor
      t_cols = MIN(nbw, n_cols - (i-1)*nbw)
      if (t_cols <= 0) exit
      t_rows = (i - 1) * nbw
      tmat_complete(t_rows+1:t_rows+t_cols,t_rows+1:t_rows+t_cols) = tmat(1:t_cols,1:t_cols,(istep-1)*blocking_factor + i)

      if (i > 1) then
        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
          else 
          ! replace this with a cugemm if beneficial
          call obj%timer%start("blas")
          call SGEMM('T', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
          endif
        else ! useGPU
          call obj%timer%start("blas")
          call SGEMM('T', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU


        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_REAL4, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_REAL4, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif

        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
          else
            ! could use the cugemm here if beneficial
            call obj%timer%start("blas")
            call STRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call STRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                              tmat_complete(t_rows+1,t_rows+1), &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
            
          endif

        else ! useGPU
          call obj%timer%start("blas")
          call STRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call STRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                            tmat_complete(t_rows+1,t_rows+1), &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call obj%timer%stop("blas")
          tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
        endif ! useGPU


      endif
    enddo

    ! Q = Q - V * T**T * V**T * Q

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(hvm_dev, int(loc(hvm),kind=c_intptr_t), &
                          max_local_rows*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: hvm -> hvm_dev", 713,  successGPU)
          call obj%timer%start("gpublas")
          call gpublas_SGEMM('T', 'N', &
                                       n_cols, l_cols, l_rows, ONE, hvm_dev, max_local_rows, &
                                       q_dev, ldq , ZERO, tmp_dev, n_cols)
          call obj%timer%stop("gpublas")

          ! copy data from device to host for a later MPI_ALLREDUCE
          successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                        tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 726,  successGPU)
        endif
      else ! useGPU
        call obj%timer%start("blas")
        call SGEMM('T', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                            hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, &
                           int(n_cols,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    else ! l_rows>0
        tmp1(1:l_cols*n_cols) = 0.0_rck
    endif ! l_rows>0

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                   tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 764,  successGPU)

    endif

    if (useNonBlockingCollectivesRows) then
      call obj%timer%start("mpi_nbc_communication")
      call mpi_iallreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
      call obj%timer%stop("mpi_nbc_communication")
      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
    else
      call obj%timer%start("mpi_communication")
      call mpi_allreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_REAL4, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
    endif

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          call STRMM('L', 'U', 'T', 'N', &
                              int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
          call SGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                              int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                              q_mat, int(ldq,kind=BLAS_KIND))
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2),kind=c_intptr_t), &
                        l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp2 -> tmp_dev", 819,  successGPU)

          ! needed: as long as not device to device copy
          successGPU = gpu_memcpy(tmat_dev, int(loc(tmat_complete),kind=c_intptr_t), &
                        cwy_blocking*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmat_complete -> tmat_dev", 825,  successGPU)

          call obj%timer%start("gpublas")
          call gpublas_STRMM('L', 'U', 'T', 'N', &
                                   n_cols, l_cols, ONE, tmat_dev, cwy_blocking, tmp_dev, n_cols)
          call gpublas_SGEMM('N', 'N', l_rows, l_cols, n_cols, -ONE, hvm_dev, max_local_rows, tmp_dev, &
                                     n_cols, ONE, q_dev, ldq)
          call obj%timer%stop("gpublas")
        endif
      else
        call obj%timer%start("blas")
        call STRMM('L', 'U', 'T', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
        call SGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                            int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                            int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                            q_mat, int(ldq,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU

    endif

  enddo ! istep

  deallocate(hvb, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: hvb", 902,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      deallocate(tmp1, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp1", 907,  istat,  errorMessage)

      deallocate(tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp2", 910,  istat,  errorMessage)

      deallocate(hvm, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: hvm", 913,  istat,  errorMessage)

    else
      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: hvm_dev", 917,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmp_dev", 920,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmat_dev", 923,  successGPU)

      ! final transfer of q_dev
      successGPU = gpu_memcpy(int(loc(q_mat),kind=c_intptr_t), q_dev, ldq*matrixCols*size_of_datatype, &
                    gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_dev -> q_mat", 928,  successGPU)

      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: q_dev", 931,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: q_mat", 934,  successGPU)
      nullify(tmp1)
      nullify(tmp2)
      nullify(hvm)

      ! take care of new pointers nullify them
      nullify(tmat_mpi_deviceptr)
      nullify(hvm_mpi_deviceptr)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp1_host", 948,  successGPU)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp2_host", 951,  successGPU)

      successGPU = gpu_free_host(hvm_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: hvm_host", 954,  successGPU)

      successGPU = gpu_host_unregister(int(loc(tmat_complete),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: tmat_complete", 957,  successGPU)
    endif
  else ! useGPU
    deallocate(tmp1, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp1", 961,  istat,  errorMessage)

    deallocate(tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp2", 964,  istat,  errorMessage)

    deallocate(hvm, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: hvm", 967,  istat,  errorMessage)
  endif ! useGPU

  deallocate(tmat_complete, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: tmat_complete", 971,  istat,  errorMessage)


  if (blocking_factor > 1) then
    deallocate(t_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp", 984,  istat,  errorMessage)

    deallocate(t_tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp2", 987,  istat,  errorMessage)
  endif

  call obj%timer%stop("trans_ev_band_to_full_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

end subroutine trans_ev_band_to_full_&
&real&
    &_&
    &single






subroutine tridiag_band_&
  &real&
  &_&
  &single &
  (obj, na, nb, nblk, a_mat, lda, d, e, matrixCols, &
  hh_trans, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, useGPU, wantDebug, nrThreads, isSkewsymmetric)
  !-------------------------------------------------------------------------------
  ! tridiag_band_real/complex:
  ! Reduces a real symmetric band matrix to tridiagonal form
  !
  !  na          Order of matrix a
  !
  !  nb          Semi bandwith
  !
  !  nblk        blocksize of cyclic distribution, must be the same in both directions!
  !
  !  a_mat(lda,matrixCols)    Distributed system matrix reduced to banded form in the upper diagonal
  !
  !  lda         Leading dimension of a
  !  matrixCols  local columns of matrix a
  !
  ! hh_trans : housholder vectors
  !
  !  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  mpi_comm_rows
  !  mpi_comm_cols
  !              MPI-Communicators for rows/columns
  !  mpi_comm_all
  !              MPI-Communicator for the total processor set
  !-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use precision
  use, intrinsic :: iso_c_binding
  use redist
  use elpa_blas_interfaces
  use elpa_skewsymmetric_blas
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

  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU, wantDebug
  logical, intent(in)                          :: isSkewsymmetric
  integer(kind=ik), intent(in)                 :: na, nb, nblk, lda, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  real(kind=rck), intent(in)         :: a_mat(lda,*)
  real(kind=rk), intent(out)        :: d(na), e(na) ! set only on PE 0
  real(kind=rck), intent(out), allocatable   :: hh_trans(:,:)

  real(kind=rk)                     :: vnorm2
  real(kind=rck)                     :: hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
  real(kind=rck)                     :: hd(nb), hs(nb)

  integer(kind=ik)                             :: i, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
  integer(kind=ik)                             :: my_pe, n_pes
  integer(kind=ik)                             :: my_prow, np_rows, my_pcol, np_cols
  integer(kind=MPI_KIND)                       :: my_peMPI, n_pesMPI, mpierr
  integer(kind=MPI_KIND)                       :: my_prowMPI, np_rowsMPI, my_pcolMPI, np_colsMPI
  integer(kind=MPI_KIND)                       :: ireq_ab, ireq_hv
  integer(kind=ik)                             :: na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off
  integer(kind=ik), intent(in)                 :: nrThreads
  integer(kind=ik), allocatable                :: global_id(:,:), hh_cnt(:), hh_dst(:)
  integer(kind=MPI_KIND), allocatable          :: ireq_hhr(:), ireq_hhs(:)
  integer(kind=ik), allocatable                :: limits(:), snd_limits(:,:)
  integer(kind=ik), allocatable                :: block_limits(:)
  real(kind=rck), allocatable         :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)
  integer                                      :: istat
  integer(kind=ik)                             :: nblockEnd
  character(200)                               :: errorMessage
  character(20)                                :: gpuString

  logical                                      :: useIntelGPU

   integer(kind=MPI_KIND)                      :: allreduce_request1, allreduce_request2
   logical                                     :: useNonBlockingCollectivesAll
   integer(kind=c_int)                         :: non_blocking_collectives, error

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_band_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_all_elpa2_band_to_tridi", non_blocking_collectives, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives in elpa2_band_to_tridi. Aborting..."
    stop
  endif

  if (non_blocking_collectives .eq. 1) then
    useNonBlockingCollectivesAll = .true.
  else
    useNonBlockingCollectivesAll = .false.
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND) ,my_peMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_all,kind=MPI_KIND) ,n_pesMPI ,mpierr)

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND),my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND),np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND),my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND),np_colsMPI ,mpierr)

  my_pe = int(my_peMPI,kind=MPI_KIND)
  n_pes = int(n_pesMPI,kind=MPI_KIND)
  my_prow = int(my_prowMPI,kind=MPI_KIND)
  np_rows = int(np_rowsMPI,kind=MPI_KIND)
  my_pcol = int(my_pcolMPI,kind=MPI_KIND)
  np_cols = int(np_colsMPI,kind=MPI_KIND)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Get global_id mapping 2D procssor coordinates to global id

  allocate(global_id(0:np_rows-1,0:np_cols-1), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: global_id", 202,  istat,  errorMessage)

  global_id(:,:) = 0
  global_id(my_prow, my_pcol) = my_pe


  if (wantDebug) call obj%timer%start("mpi_communication")
  if (useNonBlockingCollectivesAll) then
    call mpi_iallreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
  else
    call mpi_allreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), mpierr)
  endif
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Total number of blocks in the band:

  nblocks_total = (na-1)/nb + 1

  ! Set work distribution

  allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: block_limits", 246,  istat,  errorMessage)

  call divide_band(obj,nblocks_total, n_pes, block_limits)

  ! nblocks: the number of blocks for my task
  nblocks = block_limits(my_pe+1) - block_limits(my_pe)

  ! allocate the part of the band matrix which is needed by this PE
  ! The size is 1 block larger than needed to avoid extensive shifts
  allocate(ab(2*nb,(nblocks+1)*nb), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: ab", 256,  istat,  errorMessage)

  ab = 0.0_rck ! needed for lower half, the extra block should also be set to 0 for safety

  ! n_off: Offset of ab within band
  n_off = block_limits(my_pe)*nb

  ! Redistribute band in a to ab
  call redist_band_&
  &real&
  &_&
  &single&
  &(obj,a_mat, lda, na, nblk, nb, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, ab)

  ! Calculate the workload for each sweep in the back transformation
  ! and the space requirements to hold the HH vectors

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: limits", 274,  istat,  errorMessage)

  call determine_workload(obj,na, nb, np_rows, limits)
  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  do n = 1, nblocks_total
    call determine_workload(obj, nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    ! add to number of householder vectors
    ! please note: for nx==1 the one and only HH Vector is 0 and is neither calculated nor send below!
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_hh_vecs = num_hh_vecs + local_size
      num_chunks  = num_chunks+1
    endif
    nx = nx - nb
  enddo

  ! Allocate space for HH vectors

  allocate(hh_trans(nb,num_hh_vecs), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_trans", 297,  istat,  errorMessage)

  ! Allocate and init MPI requests

  allocate(ireq_hhr(num_chunks), stat=istat, errmsg=errorMessage) ! Recv requests
  call check_allocate_f("tridiag_band: ireq_hhr", 302,  istat,  errorMessage)
  allocate(ireq_hhs(nblocks), stat=istat, errmsg=errorMessage)    ! Send requests
  call check_allocate_f("tridiag_band: ireq_hhs", 304,  istat,  errorMessage)

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  nt = 0
  nBlockEnd=1
  do n = 1, nblocks_total
    call determine_workload(obj,nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_chunks  = num_chunks+1
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_irecv(hh_trans(1,num_hh_vecs+1), int(nb*local_size,kind=MPI_KIND),  MPI_REAL4,     &
                    int(nt,kind=MPI_KIND), int(10+n-block_limits(nt),kind=MPI_KIND), &
                    int(mpi_comm_all,kind=MPI_KIND), ireq_hhr(num_chunks), mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")

      num_hh_vecs = num_hh_vecs + local_size
    endif
    nx = nx - nb
    if (n == block_limits(nt+1)) then
      nt = nt + 1
    endif
  enddo
  ireq_hhs(:) = MPI_REQUEST_NULL
  ! Buffers for gathering/sending the HH vectors

  allocate(hh_gath(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! gathers HH vectors
  call check_allocate_f("tridiag_band: hh_gath", 341,  istat,  errorMessage)

  allocate(hh_send(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! send buffer for HH vectors
  call check_allocate_f("tridiag_band: hh_send", 344,  istat,  errorMessage)

  hh_gath(:,:,:) = 0.0_rck
  hh_send(:,:,:) = 0.0_rck

  ! Some counters

  allocate(hh_cnt(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_cnt", 352,  istat,  errorMessage)

  allocate(hh_dst(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_dst", 355,  istat,  errorMessage)

  hh_cnt(:) = 1 ! The first transfomation Vector is always 0 and not calculated at all
  hh_dst(:) = 0 ! PE number for receive
  ireq_ab = MPI_REQUEST_NULL
  ireq_hv = MPI_REQUEST_NULL
  ! Limits for sending

  allocate(snd_limits(0:np_rows,nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: snd_limits", 366,  istat,  errorMessage)

  do iblk=1,nblocks
    call determine_workload(obj, na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
  enddo


  ! ---------------------------------------------------------------------------
  ! Start of calculations

  na_s = block_limits(my_pe)*nb + 1

  if (my_pe>0 .and. na_s<=na) then
    ! send first column to previous PE
    ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
    ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
    if (wantDebug) call obj%timer%start("mpi_communication")
    call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_REAL4, &
                   int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), ireq_ab, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif


   do istep=1,na-nblockEnd

     if (my_pe==0) then
       n = MIN(na-na_s,nb) ! number of rows to be reduced
       hv(:) = 0.0_rck
       hd(:) = 0.0_rck
       tau = 0.0_rck

       ! Transform first column of remaining matrix
       ! The last step (istep=na-1) is only needed for sending the last HH vectors.
       ! We don't want the sign of the last element flipped (analogous to the other sweeps)

       if (istep < na-1) then
         ! Transform first column of remaining matrix
         vnorm2 = sum(ab(3:n+1,na_s-n_off)**2)

          call hh_transform_&
               &real&
               &_&
               &single &
                          (obj, ab(2,na_s-n_off), vnorm2, hf, tau, wantDebug)

          hv(1) = 1.0_rck
          hv(2:n) = ab(3:n+1,na_s-n_off)*hf
       endif

       if (isSkewsymmetric) then
         d(istep) = 0.0_rk
       else
         d(istep) = ab(1,na_s-n_off)
       endif
       e(istep) = ab(2,na_s-n_off)

       if (istep == na-1) then
         if (isSkewsymmetric) then
           d(na) = 0
         else
           d(na) = ab(1,na_s+1-n_off)
         endif

         e(na) = 0.0_rck
       endif
     else
       if (na>na_s) then
         ! Receive Householder Vector from previous task, from PE owning subdiagonal


         if (wantDebug) call obj%timer%start("mpi_communication")

         call mpi_recv(hv, int(nb,kind=MPI_KIND), MPI_REAL4, &
                       int(my_pe-1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")


         tau = hv(1)
         hv(1) = 1.0_rck
       endif
      endif

      na_s = na_s+1
      if (na_s-n_off > nb) then
        ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
        ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rck
        n_off = n_off + nb
      endif


      do iblk=1,nblocks
        ns = na_s + (iblk-1)*nb - n_off ! first column in block
        ne = ns+nb-1                    ! last column in block

        if (ns+n_off>na) exit

        ! Store Householder Vector for back transformation

        hh_cnt(iblk) = hh_cnt(iblk) + 1

        hh_gath(1   ,hh_cnt(iblk),iblk) = tau
        hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

        if (hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
          ! Wait for last transfer to finish
          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
          ! Copy vectors into send buffer
          hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
          ! Send to destination

          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_isend(hh_send(1,1,iblk), int(nb*hh_cnt(iblk),kind=MPI_KIND), MPI_REAL4, &
                         global_id(hh_dst(iblk), mod(iblk+block_limits(my_pe)-1,np_cols)), &
                         int(10+iblk,kind=MPI_KIND), int(mpi_comm_all,kind=MPI_KIND), ireq_hhs(iblk), mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! Reset counter and increase destination row
          hh_cnt(iblk) = 0
          hh_dst(iblk) = hh_dst(iblk)+1
        endif

        ! The following code is structured in a way to keep waiting times for
        ! other PEs at a minimum, especially if there is only one block.
        ! For this reason, it requests the last column as late as possible
        ! and sends the Householder Vector and the first column as early
        ! as possible.
        nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
        nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                      ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

        ! Multiply diagonal block and subdiagonal block with Householder Vector

        if (iblk==nblocks .and. nc==nb) then

          ! We need the last column from the next PE.
          ! First do the matrix multiplications without last column ...

          ! Diagonal block, the contribution of the last element is added below!
          ab(1,ne) = 0.0_rck

          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")

            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_sssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call SSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif
            ! Subdiagonal block
            if (nr>0) call SGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")

            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_sssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call SSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif
            ! Subdiagonal block
            if (nr>0) call SGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif

          ! ... then request last column ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_recv(ab(1,ne), int(nb+1,kind=MPI_KIND), MPI_REAL4,  &
                       int(my_pe+1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and complete the result
          hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
          hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

        else

          ! Normal matrix multiply
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_sssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call SSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                                hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif

            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            if (isSkewsymmetric) then
              hd(:) = 0.0_rk
              call elpa_sssmv(int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), hv, hd)
            else
              call SSYMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                                hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            endif
            if (nr>0) call SGEMV('N', int(nr,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), tau, ab(nb+1,ns), &
                                        int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Calculate first column of subdiagonal block and calculate new
        ! Householder transformation for this column
        hv_new(:) = 0.0_rck ! Needed, last rows must be 0 for nr < nb
        tau_new = 0.0_rck
        if (nr>0) then

          ! complete (old) Householder transformation for first column

          ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

          ! calculate new Householder transformation ...
          if (nr>1) then
            vnorm2 = sum(ab(nb+2:nb+nr,ns)**2)

            call hh_transform_&
                &real&
                &_&
                &single &
               (obj, ab(nb+1,ns), vnorm2, hf, tau_new, wantDebug)
            hv_new(1) = 1.0_rck
            hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
            ab(nb+2:,ns) = 0.0_rck
          endif ! nr > 1

          ! ... and send it away immediatly if this is the last block

          if (iblk==nblocks) then
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            hv_s(1) = tau_new
            hv_s(2:) = hv_new(2:)

            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_isend(hv_s, int(nb,kind=MPI_KIND), MPI_REAL4, &
                           int(my_pe+1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                           ireq_hv, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

          endif

        endif

        ! Transform diagonal block
        if (.NOT. isSkewsymmetric) then
          x = dot_product(hv(1:nc),hd(1:nc))*tau
        endif

        if (.NOT. isSkewsymmetric) then
          hd(1:nc) = hd(1:nc) - 0.5_rk*x*hv(1:nc)
        endif
        if (my_pe>0 .and. iblk==1) then

          ! The first column of the diagonal block has to be send to the previous PE
          ! Calculate first column only ...
          if (isSkewsymmetric) then
            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*hv(1) + hv(1:nc)*hd(1)
          else
            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*hv(1) - hv(1:nc)*hd(1)
          endif
          ! ... send it away ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_wait(ireq_ab, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ab_s(1:nb+1) = ab(1:nb+1,ns)

          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_REAL4, &
                         int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                         ireq_ab, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and calculate remaining columns with rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (isSkewsymmetric) then
              if (nc>1) call elpa_sssr2(int(nc-1,kind=BLAS_KIND), hd(2), hv(2), ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            else
              if (nc>1) call SSYR2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                         hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            if (isSkewsymmetric) then 
              if (nc>1) call elpa_sssr2(int(nc-1,kind=BLAS_KIND), hd(2), hv(2), ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            else
              if (nc>1) call SSYR2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                         hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("blas")
          endif

        else
          ! No need to  send, just a rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (isSkewsymmetric) then
              call elpa_sssr2(int(nc,kind=BLAS_KIND), hd, hv, ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            else
              call SSYR2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND,  &
                                  hv, 1_BLAS_KIND, ab(1,ns), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")
            if (isSkewsymmetric) then 
              call elpa_sssr2(int(nc,kind=BLAS_KIND), hd, hv, ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            else
              call SSYR2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND,  &
                                hv, 1_BLAS_KIND, ab(1,ns), int(2*nb-1,kind=BLAS_KIND) )
            endif
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

        if (nr>0) then
          if (nr>1) then
            if (useIntelGPU) then
              if (wantDebug) call obj%timer%start("mkl_offload")
              call SGEMV('T', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("mkl_offload")
            else
              if (wantDebug) call obj%timer%start("blas")
              call SGEMV('T', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("blas")
            endif
            x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
            h(2:nb) = h(2:nb) - x*hv(2:nb)
            ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update
            do i=2,nb
              ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*h(i) - hs(1:nr)*hv(i)
            enddo
          else
            ! No double Householder transformation for nr=1, just complete the row
            do i=2,nb
              ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*hv(i)
            enddo
          endif
        endif

        ! Use new HH Vector for the next block
        hv(:) = hv_new(:)
        tau = tau_new

      enddo


  enddo ! istep

  ! Finish the last outstanding requests


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
  call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

  call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
  call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_barrier(int(mpi_comm_all,kind=MPI_KIND),mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")
  deallocate(ab, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ab", 1363,  istat,  errorMessage)

  deallocate(ireq_hhr, ireq_hhs, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ireq_hhr", 1366,  istat,  errorMessage)

  deallocate(hh_cnt, hh_dst, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_dst", 1369,  istat,  errorMessage)

  deallocate(hh_gath, hh_send, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_gath", 1372,  istat,  errorMessage)

  deallocate(limits, snd_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: limits", 1375,  istat,  errorMessage)

  deallocate(block_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: block_limits", 1378,  istat,  errorMessage)

  deallocate(global_id, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: global_id", 1381,  istat,  errorMessage)

  call obj%timer%stop("tridiag_band_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

  ! intel compiler bug makes these ifdefs necessary
end subroutine tridiag_band_real_&
&single

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






subroutine trans_ev_tridi_to_band_&
&real&
&_&
&single &
(obj, na, nev, nblk, nbw, q, ldq, matrixCols,         &
 hh_trans, mpi_comm_rows, mpi_comm_cols, wantDebug, useGPU, max_threads, success, &
 kernel)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_real/complex:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!  matrixCols  local columns of matrix q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use pack_unpack_cpu
  use pack_unpack_gpu
  use compute_hh_trafo
  use elpa_gpu
  use precision
  use, intrinsic :: iso_c_binding
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

  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU

  integer(kind=ik), intent(in)                 :: kernel
  integer(kind=ik), intent(in)                 :: na, nev, nblk, nbw, ldq, matrixCols, mpi_comm_rows, mpi_comm_cols

  real(kind=rck), target              :: q(ldq,*)

  real(kind=rck), intent(in),target   :: hh_trans(:,:)
  integer(kind=c_intptr_t)                     :: hh_trans_dev
  type(c_ptr)                                  :: hh_trans_mpi_dev
  real(kind=rck), pointer             :: hh_trans_mpi_fortran_ptr(:,:)


  integer(kind=ik)                             :: np_rows, my_prow, np_cols, my_pcol
  integer(kind=MPI_KIND)                       :: np_rowsMPI, my_prowMPI, np_colsMPI, my_pcolMPI
  integer(kind=ik)                             :: i, j, ip, sweep, nbuf, l_nev, a_dim2
  integer(kind=ik)                             :: current_n, current_local_n, current_n_start, current_n_end
  integer(kind=ik)                             :: next_n, next_local_n, next_n_start, next_n_end
  integer(kind=ik)                             :: bottom_msg_length, top_msg_length, next_top_msg_length
  integer(kind=ik)                             :: stripe_width, last_stripe_width, stripe_count
  integer(kind=ik)                             :: num_result_blocks, num_result_buffers, num_bufs_recvd
  integer(kind=ik)                             :: a_off, current_tv_off, max_blk_size
  integer(kind=ik)                             :: src, src_offset, dst, offset, nfact, num_blk
  integer(kind=MPI_KIND)                       :: mpierr

  logical                                      :: flag
  real(kind=rck), pointer             :: aIntern(:,:,:)
  real(kind=rck)                      :: a_var

  type(c_ptr)                                  :: aIntern_ptr

  real(kind=rck), allocatable, target :: row(:)
  
  integer(kind=c_intptr_t)                     :: row_dev
  type(c_ptr)                                  :: row_mpi_dev
  real(kind=rck), pointer             :: row_mpi_fortran_ptr(:)

  real(kind=rck), pointer             :: row_group(:,:)

  real(kind=rck), allocatable, target :: top_border_send_buffer(:,:)
  real(kind=rck), allocatable, target :: top_border_recv_buffer(:,:)
  real(kind=rck), allocatable, target :: bottom_border_send_buffer(:,:)
  real(kind=rck), allocatable, target :: bottom_border_recv_buffer(:,:)

  integer(kind=c_intptr_t)                     :: top_border_recv_buffer_dev, top_border_send_buffer_dev
  type(c_ptr)                                  :: top_border_recv_buffer_mpi_dev, top_border_send_buffer_mpi_dev 
  real(kind=rck), pointer             :: top_border_recv_buffer_mpi_fortran_ptr(:,:), &
                                                  top_border_send_buffer_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                     :: bottom_border_send_buffer_dev, bottom_border_recv_buffer_dev
  type(c_ptr)                                  :: bottom_border_send_buffer_mpi_dev, bottom_border_recv_buffer_mpi_dev
  real(kind=rck), pointer             :: bottom_border_send_buffer_mpi_fortran_ptr(:,:), &
                                                  bottom_border_recv_buffer_mpi_fortran_ptr(:,:)
  type(c_ptr)                                  :: aIntern_mpi_dev
  real(kind=rck), pointer             :: aIntern_mpi_fortran_ptr(:,:,:)

  integer(kind=c_intptr_t)                     :: aIntern_dev
  integer(kind=c_intptr_t)                     :: bcast_buffer_dev

  type(c_ptr)                                  :: bcast_buffer_mpi_dev
  real(kind=rck), pointer             :: bcast_buffer_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: num
  integer(kind=c_intptr_t)                     :: dev_offset, dev_offset_1
  integer(kind=c_intptr_t)                     :: row_group_dev

  type(c_ptr)                                  :: row_group_mpi_dev
  real(kind=rck), pointer             :: row_group_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: q_dev
  type(c_ptr)                                  :: q_mpi_dev
  real(kind=rck), pointer             :: q_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: hh_tau_dev
  integer(kind=ik)                             :: row_group_size, unpack_idx

  type(c_ptr)                                  :: row_group_host, bcast_buffer_host

  integer(kind=ik)                             :: n_times
  integer(kind=ik)                             :: chunk, this_chunk

  real(kind=rck), allocatable,target  :: result_buffer(:,:,:)
  integer(kind=c_intptr_t)                     :: result_buffer_dev

  type(c_ptr)                                  :: result_buffer_mpi_dev
  real(kind=rck), pointer             :: result_buffer_mpi_fortran_ptr(:,:,:)

  real(kind=rck), pointer             :: bcast_buffer(:,:)

  integer(kind=ik)                             :: n_off

  integer(kind=MPI_KIND), allocatable          :: result_send_request(:), result_recv_request(:)
  integer(kind=ik), allocatable                :: limits(:)
  integer(kind=MPI_KIND), allocatable          :: top_send_request(:), bottom_send_request(:)
  integer(kind=MPI_KIND), allocatable          :: top_recv_request(:), bottom_recv_request(:)

  ! MPI send/recv tags, arbitrary

  integer(kind=ik), parameter                  :: bottom_recv_tag = 111
  integer(kind=ik), parameter                  :: top_recv_tag    = 222
  integer(kind=ik), parameter                  :: result_recv_tag = 333

  integer(kind=ik), intent(in)                 :: max_threads



  ! Just for measuring the kernel performance
  real(kind=c_double)                          :: kernel_time, kernel_time_recv ! MPI_WTIME always needs double
  ! long integer
  integer(kind=lik)                            :: kernel_flops, kernel_flops_recv

  logical, intent(in)                          :: wantDebug
  logical                                      :: success
  integer(kind=ik)                             :: istat, print_flops
  character(200)                               :: errorMessage
  character(20)                                :: gpuString
  logical                                      :: successGPU
  integer(kind=ik)                             :: error
  integer(kind=c_intptr_t), parameter          :: size_of_datatype = size_of_&
                                                                 &single&
                                                                 &_&
                                                                 &real
  logical, parameter                           :: allComputeOnGPU = .false.

  integer(kind=MPI_KIND)                     :: bcast_request1, allreduce_request1, allreduce_request2, &
                                                allreduce_request3, allreduce_request4
  logical                                    :: useNonBlockingCollectivesCols
  logical                                    :: useNonBlockingCollectivesRows
  integer(kind=c_int)                        :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_tridi_to_band_&
  &real&
  &" // &
  &"_single" //&
  gpuString)

  call obj%get("nbc_row_elpa2_tridi_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_tridi_to_band. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_tridi_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_tridi_to_band. Aborting..."
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

  n_times = 0
  if (useGPU) then
    unpack_idx = 0
    row_group_size = 0
  endif

  success = .true.
  kernel_time = 0.0
  kernel_flops = 0

  if (wantDebug) call obj%timer%start("mpi_communication")
  call MPI_Comm_rank(int(mpi_comm_rows,kind=MPI_KIND) , my_prowMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_rows,kind=MPI_KIND) , np_rowsMPI , mpierr)
  call MPI_Comm_rank(int(mpi_comm_cols,kind=MPI_KIND) , my_pcolMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_cols,kind=MPI_KIND) , np_colsMPI , mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")

  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &real&
                            &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &real&
                            &: band backtransform works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  nfact = nbw / nblk


  ! local number of eigenvectors
  l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

  if (l_nev==0) then
    stripe_width = 0
    stripe_count = 0
    last_stripe_width = 0

  else ! l_nev


    ! Suggested stripe width is 48 since 48*64 real*8 numbers should fit into
    ! every primary cache
    ! Suggested stripe width is 48 - should this be reduced for the complex case ???

    if (useGPU) then
      stripe_width = 1024 ! Must be a multiple of 4
      stripe_count = (l_nev - 1) / stripe_width + 1

    else ! useGPU
      call obj%get("stripewidth_real",stripe_width, error)

      !stripe_width = 96 ! Must be a multiple of 8
      stripe_width = 2 * stripe_width


      stripe_count = (l_nev-1)/stripe_width + 1

      ! Adapt stripe width so that last one doesn't get too small

      stripe_width = (l_nev-1)/stripe_count + 1

      if (kernel .eq. ELPA_2STAGE_REAL_AVX512_BLOCK2 .or. &
          kernel .eq. ELPA_2STAGE_REAL_AVX512_BLOCK4 .or. &
          kernel .eq. ELPA_2STAGE_REAL_AVX512_BLOCK6 .or. &
          kernel .eq. ELPA_2STAGE_REAL_SVE512_BLOCK2 .or. &
          kernel .eq. ELPA_2STAGE_REAL_SVE512_BLOCK4 .or. &
          kernel .eq. ELPA_2STAGE_REAL_SVE512_BLOCK6  &
          ) then


       stripe_width = ((stripe_width+15)/16)*16 ! Must be a multiple of 16 because of AVX-512 memory alignment of 64 bytes
                                               ! (16 * sizeof(float) == 64)

     else
       stripe_width = ((stripe_width+7)/8)*8 ! Must be a multiple of 8 because of AVX/SSE memory alignment of 32 bytes
                                            ! (8 * sizeof(float) == 32)
     endif

   endif ! useGPU

   last_stripe_width = l_nev - (stripe_count-1)*stripe_width

  endif ! l_nev

  ! Determine the matrix distribution at the beginning

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: limits", 604,  istat,  errorMessage)
  call determine_workload(obj,na, nbw, np_rows, limits)

  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  a_dim2 = max_blk_size + nbw

  if (useGPU) then

    if (allComputeOnGPU) then
      if (wantDebug) call obj%timer%start("cuda_memcpy")

      successGPU = gpu_malloc(q_dev, ldq*matrixCols* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: q_dev", 617,  successGPU)
      successGPU =  gpu_memcpy(q_dev, int(loc(q(1,1)),kind=c_intptr_t),  &
                               ldq*matrixCols * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q -> q_dev", 621,  successGPU)

      ! associate with c_ptr
      q_mpi_dev = transfer(q_dev, q_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(q_mpi_dev, q_mpi_fortran_ptr, &
                       [ldq,matrixCols])
      if (wantDebug) call obj%timer%stop("cuda_memcpy")

      successGPU = gpu_malloc(hh_trans_dev, size(hh_trans,dim=1)*size(hh_trans,dim=2)* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 631,  successGPU)
      ! associate with c_ptr
      hh_trans_mpi_dev = transfer(hh_trans_dev, hh_trans_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(hh_trans_mpi_dev, hh_trans_mpi_fortran_ptr, &
                       [size(hh_trans,dim=1),size(hh_trans,dim=2)])
      successGPU =  gpu_memcpy(c_loc(hh_trans_mpi_fortran_ptr(1,1)),  &
                               c_loc(hh_trans(1,1)), &
                               size(hh_trans,dim=1)*size(hh_trans,dim=2) * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band: hh_trans -> hh_trans_dev", 641,  successGPU)

    endif ! allComputeOnGPU

    num =  (stripe_width*a_dim2*stripe_count)* size_of_datatype
    successGPU = gpu_malloc(aIntern_dev, stripe_width*a_dim2*stripe_count* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 647,  successGPU)

    successGPU = gpu_memset(aIntern_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 650,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      aIntern_mpi_dev = transfer(aIntern_dev, aIntern_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(aIntern_mpi_dev, aIntern_mpi_fortran_ptr, &
                       [stripe_width,a_dim2,stripe_count])
    endif ! allComputeOnGPU

    ! "row_group" and "row_group_dev" are needed for GPU optimizations
    successGPU = gpu_malloc_host(row_group_host,l_nev*nblk*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 667,  successGPU)
    call c_f_pointer(row_group_host, row_group, (/l_nev,nblk/))

    row_group(:, :) = 0.0_rck
    num =  (l_nev*nblk)* size_of_datatype
    successGPU = gpu_malloc(row_group_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 673,  successGPU)

    successGPU = gpu_memset(row_group_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_group_dev", 676,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      row_group_mpi_dev = transfer(row_group_dev, row_group_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(row_group_mpi_dev, row_group_mpi_fortran_ptr, &
                       [l_nev,nblk])
    endif ! allComputeOnGPU

  else ! GPUs are not used



    if (posix_memalign(aIntern_ptr, 64_c_intptr_t, stripe_width*a_dim2*stripe_count*  &
        C_SIZEOF(a_var)) /= 0) then
      print *,"trans_ev_tridi_to_band_real: error when allocating aIntern"//errorMessage
      stop 1
    endif

    call c_f_pointer(aIntern_ptr, aIntern,[stripe_width,a_dim2,stripe_count] )
    !allocate(aIntern(stripe_width,a_dim2,stripe_count), stat=istat, errmsg=errorMessage)

    aIntern(:,:,:) = 0.0_rck
  endif !useGPU

  allocate(row(l_nev), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: row", 723,  istat,  errorMessage)

  row(:) = 0.0_rck

  if (useGPU .and. allComputeOnGPU) then
    num =  (l_nev)* size_of_datatype
    successGPU = gpu_malloc(row_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_dev", 730,  successGPU)

    successGPU = gpu_memset(row_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_dev", 733,  successGPU)

    ! associate with c_ptr
    row_mpi_dev = transfer(row_dev, row_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(row_mpi_dev, row_mpi_fortran_ptr, &
                     [l_nev])
  endif

  ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
  ! and transpose the matrix using stripes of given stripe_width for cache blocking.

  ! The peculiar way it is done below is due to the fact that the last row should be
  ! ready first since it is the first one to start below



  if (wantDebug) call obj%timer%start("ip_loop")
  do ip = np_rows-1, 0, -1
    if (my_prow == ip) then
      ! Receive my rows which have not yet been received
      src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
      do i=limits(ip)+1,limits(ip+1)
        src = mod((i-1)/nblk, np_rows)

        if (src < my_prow) then

          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
            &real&
            &_gpu_&
            &single &
                          ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                                      stripe_width, last_stripe_width, a_dim2, l_nev,&
                                      row_group_size, nblk, unpack_idx, &
                                       i - limits(ip), .false., wantDebug, allComputeOnGPU)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                         int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            call unpack_row_&
                &real&
                &_cpu_&
                &single &
                (obj,aIntern, row,i-limits(ip), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        elseif (src == my_prow) then

          src_offset = src_offset+1

          if (useGPU) then

            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &real&
                 &_gpu_&
                 &single &
             ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                          stripe_width, last_stripe_width, a_dim2, l_nev,&
                          row_group_size, nblk, unpack_idx, &
                          i - limits(ip), .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_SCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_group_mpi_fortran_ptr(1,row_group_size)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            else ! allComputeOnGPU
              row_group(:, row_group_size) = q(src_offset, 1:l_nev)
            endif ! allComputeOnGPU
          else ! useGPU
            row(:) = q(src_offset, 1:l_nev)
          endif ! useGPU


          if (useGPU) then

          else
            call unpack_row_&
                 &real&
                 &_cpu_&
                 &single &
                            (obj,aIntern, row,i-limits(ip),  stripe_count, stripe_width, last_stripe_width)
          endif


        endif ! src == my_prow
      enddo ! i=limits(ip)+1,limits(ip+1)


      ! Send all rows which have not yet been send


      if (allComputeOnGPU .and. useGPU) then
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1

              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_SCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1011,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

              if (wantDebug) call obj%timer%start("cuda_mpi_communication")
              call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

            endif
          enddo
        enddo
      else !  allComputeOnGPU
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1
              row(:) = q(src_offset, 1:l_nev)

              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            endif
          enddo
        enddo
      endif ! allComputeOnGPU

    else if (my_prow < ip) then

      ! Send all rows going to PE ip
      if (allComputeOnGPU .and. useGPU) then
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1

            if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
            call gpublas_SCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
            if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

            ! is there a way to avoid this device_synchronize ?
            if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
            successGPU = gpu_devicesynchronize()
            call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1061,  successGPU)
            if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

            if (wantDebug) call obj%timer%start("cuda_mpi_communication")
            call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                          int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

          endif
        enddo
      else ! allComputeOnGPU 
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1
            row(:) = q(src_offset, 1:l_nev)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                        int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          endif
        enddo
      endif  ! allComputeOnGPU

      ! Receive all rows from PE ip
      do i=limits(my_prow)+1,limits(my_prow+1)
        src = mod((i-1)/nblk, np_rows)
        if (src == ip) then
          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &real&
                 &_gpu_&
                 &single&
                 &( obj, &
                 row_group, row_group_dev, aIntern_dev, stripe_count,  &
                 stripe_width, last_stripe_width, a_dim2, l_nev,       &
                 row_group_size, nblk, unpack_idx,                     &
                 i - limits(my_prow), .false., wantDebug, allComputeOnGPU)

            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_REAL4, &
                             int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
            call unpack_row_&
                 &real&
                 &_cpu_&
                 &single &
                 (obj,aIntern, row,i-limits(my_prow), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        endif
      enddo ! i=limits(my_prow)+1,limits(my_prow+1)
    endif ! (my_prow < ip)
  enddo ! ip = np_rows-1, 0, -1

  if (wantDebug) call obj%timer%stop("ip_loop")
 
  if (wantDebug) call obj%timer%start("allocate")

  if (useGPU) then
    ! Force an unpacking of all remaining rows that haven't been unpacked yet

    call unpack_and_prepare_row_group_&
         &real&
         &_gpu_&
         &single&
         &( obj, &
         row_group, row_group_dev, aIntern_dev, stripe_count, &
         stripe_width, last_stripe_width, &
         a_dim2, l_nev, row_group_size, nblk, unpack_idx,     &
         -1, .true., wantDebug, allComputeOnGPU)

  endif

  ! Set up result buffer queue

  num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

  num_result_buffers = 4*nfact
  allocate(result_buffer(l_nev,nblk,num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_buffer", 1266,  istat,  errorMessage)

  allocate(result_send_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_send_request", 1269,  istat,  errorMessage)

  allocate(result_recv_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_recv_request", 1272,  istat,  errorMessage)

  if (useGPU .and. allComputeOnGPU) then
    num_result_buffers = 4*nfact
    num =  (l_nev*nblk*num_result_buffers)* size_of_datatype
    successGPU = gpu_malloc(result_buffer_dev, num* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 1278,  successGPU)

    ! associate with c_ptr
    result_buffer_mpi_dev = transfer(result_buffer_dev, result_buffer_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(result_buffer_mpi_dev, result_buffer_mpi_fortran_ptr, &
                     [l_nev,nblk,num_result_buffers])
  endif

  result_send_request(:) = MPI_REQUEST_NULL
  result_recv_request(:) = MPI_REQUEST_NULL

  ! Queue up buffers
  if (wantDebug) call obj%timer%start("host_mpi_communication")
  if (my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
    do j = 1, min(num_result_buffers, num_result_blocks)
      if (useGPU) then
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL4,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      else ! useGPU
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL4,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      endif ! useGPU
    enddo
  endif
  if (wantDebug) call obj%timer%stop("host_mpi_communication")

  num_bufs_recvd = 0 ! No buffers received yet

  ! Initialize top/bottom requests

  allocate(top_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_send_request", 1335,  istat,  errorMessage)

  allocate(top_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_recv_request", 1338,  istat,  errorMessage)

  allocate(bottom_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_send_request", 1341,  istat,  errorMessage)

  allocate(bottom_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_recv_request", 1344,  istat,  errorMessage)

  top_send_request(:) = MPI_REQUEST_NULL
  top_recv_request(:) = MPI_REQUEST_NULL
  bottom_send_request(:) = MPI_REQUEST_NULL
  bottom_recv_request(:) = MPI_REQUEST_NULL


  allocate(top_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 1374,  istat,  errorMessage)

  allocate(top_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1377,  istat,  errorMessage)

  allocate(bottom_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1380,  istat,  errorMessage)

  allocate(bottom_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1383,  istat,  errorMessage)

  top_border_send_buffer(:,:) = 0.0_rck
  top_border_recv_buffer(:,:) = 0.0_rck
  bottom_border_send_buffer(:,:) = 0.0_rck
  bottom_border_recv_buffer(:,:) = 0.0_rck


  if (useGPU) then
    if (allComputeOnGPU) then
      ! top_border_recv_buffer and top_border_send_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(top_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1401,  successGPU)

      successGPU = gpu_malloc(top_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1404,  successGPU)

      successGPU = gpu_memset(top_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1407,  successGPU)

      successGPU = gpu_memset(top_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1410,  successGPU)

      ! associate with c_ptr
      top_border_recv_buffer_mpi_dev = transfer(top_border_recv_buffer_dev, top_border_recv_buffer_mpi_dev)
      top_border_send_buffer_mpi_dev = transfer(top_border_send_buffer_dev, top_border_send_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(top_border_recv_buffer_mpi_dev, top_border_recv_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])
      call c_f_pointer(top_border_send_buffer_mpi_dev, top_border_send_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])

      ! bottom_border_send_buffer and bottom_border_recv_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(bottom_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1433,  successGPU)
      successGPU = gpu_malloc(bottom_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1435,  successGPU)

      successGPU = gpu_memset(bottom_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1438,  successGPU)
      successGPU = gpu_memset(bottom_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1440,  successGPU)

      ! associate with c_ptr
      bottom_border_send_buffer_mpi_dev = transfer(bottom_border_send_buffer_dev, bottom_border_send_buffer_mpi_dev)
      bottom_border_recv_buffer_mpi_dev = transfer(bottom_border_recv_buffer_dev, bottom_border_recv_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(bottom_border_send_buffer_mpi_dev, bottom_border_send_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
      call c_f_pointer(bottom_border_recv_buffer_mpi_dev, bottom_border_recv_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
    endif ! allComputeOnGPU

    successGPU = gpu_host_register(int(loc(top_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 1462,  successGPU)

    successGPU = gpu_host_register(int(loc(top_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1467,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1472,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1477,  successGPU)
  endif ! useGPU


  ! Initialize broadcast buffer

  if (useGPU) then
    successGPU = gpu_malloc_host(bcast_buffer_host,nbw*max_blk_size*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 1485,  successGPU)
    call c_f_pointer(bcast_buffer_host, bcast_buffer, (/nbw,max_blk_size/))
  else
    allocate(bcast_buffer(nbw, max_blk_size), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_tridi_to_band: bcast_buffer", 1489,  istat,  errorMessage)
  endif

  bcast_buffer = 0.0_rck

  if (useGPU) then
    num =  ( nbw * max_blk_size) * size_of_datatype
    successGPU = gpu_malloc(bcast_buffer_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1497,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      bcast_buffer_mpi_dev = transfer(bcast_buffer_dev, bcast_buffer_mpi_dev)
      call c_f_pointer(bcast_buffer_mpi_dev, bcast_buffer_mpi_fortran_ptr, &
                      [nbw, max_blk_size])
    endif ! allComputeOnGPU

    successGPU = gpu_memset( bcast_buffer_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1507,  successGPU)

    num =  (max_blk_size)* size_of_datatype
    successGPU = gpu_malloc( hh_tau_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1511,  successGPU)

    successGPU = gpu_memset( hh_tau_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1514,  successGPU)
  endif ! useGPU

  current_tv_off = 0 ! Offset of next row to be broadcast

  ! ------------------- start of work loop -------------------

  a_off = 0 ! offset in aIntern (to avoid unnecessary shifts)

  top_msg_length = 0
  bottom_msg_length = 0
  if (wantDebug) call obj%timer%stop("allocate")

  if (wantDebug) call obj%timer%start("sweep_loop")
  do sweep = 0, (na-1)/nbw

    current_n = na - sweep*nbw
    call determine_workload(obj,current_n, nbw, np_rows, limits)
    current_n_start = limits(my_prow)
    current_n_end   = limits(my_prow+1)
    current_local_n = current_n_end - current_n_start

    next_n = max(current_n - nbw, 0)
    call determine_workload(obj,next_n, nbw, np_rows, limits)
    next_n_start = limits(my_prow)
    next_n_end   = limits(my_prow+1)
    next_local_n = next_n_end - next_n_start

    if (next_n_end < next_n) then
      bottom_msg_length = current_n_end - next_n_end
    else
      bottom_msg_length = 0
    endif

    if (next_local_n > 0) then
      next_top_msg_length = current_n_start - next_n_start
    else
      next_top_msg_length = 0
    endif

    if (sweep==0 .and. current_n_end < current_n .and. l_nev > 0) then
      if (wantDebug) call obj%timer%start("host_mpi_communication")
      do i = 1, stripe_count


        if (useGPU) then
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_REAL4, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        else !useGPU
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_REAL4, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        endif !useGPU

      enddo
      if (wantDebug) call obj%timer%stop("host_mpi_communication")
    endif

    if (current_local_n > 1) then
      if (useGPU .and. allComputeOnGPU) then
        if (my_pcol == mod(sweep,np_cols)) then
          if (wantDebug) call obj%timer%start("cuda_memcpy")
          successGPU =  gpu_memcpy(c_loc(bcast_buffer_mpi_fortran_ptr(1,1)), &
                                   c_loc(hh_trans_mpi_fortran_ptr(1,current_tv_off+1)),  &
                                     size(hh_trans,dim=1) * (current_tv_off+current_local_n-(current_tv_off+1)+1) * &   
                                     size_of_datatype, &
                                     gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1646,  successGPU)
          if (wantDebug) call obj%timer%stop("cuda_memcpy")
          current_tv_off = current_tv_off + current_local_n
        endif
      else !useGPU
        if (my_pcol == mod(sweep,np_cols)) then
          bcast_buffer(:,1:current_local_n) =    &
          hh_trans(:,current_tv_off+1:current_tv_off+current_local_n)
          current_tv_off = current_tv_off + current_local_n
        endif
      endif ! useGPU

      if (useNonBlockingCollectivesCols) then
        if (wantDebug) call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_REAL4, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_bcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_REAL4, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif
 
      if (useGPU .and. .not.(allComputeOnGPU)) then
        if (wantDebug) call obj%timer%start("memcpy")
        successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
                                 nbw * current_local_n *    &
                                 size_of_datatype, &
                                 gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1702,  successGPU)
        if (wantDebug) call obj%timer%stop("memcpy")
      endif ! useGPU

      if (useGPU) then
        !if (.not.(allComputeOnGPU)) then
        !  if (wantDebug) call obj%timer%start("memcpy")
        !  successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
        !                           nbw * current_local_n *    &
        !                           size_of_datatype, &
        !                           gpuMemcpyHostToDevice)
        !  call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1725,  successGPU)
        !  if (wantDebug) call obj%timer%stop("memcpy")
        !endif ! allComputeOnGPU

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &real&
             &_gpu_&
             &single&
             (bcast_buffer_dev, hh_tau_dev, nbw, &
             current_local_n, .false.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU

    else ! (current_local_n > 1) then

      ! for current_local_n == 1 the one and only HH Vector is 0 and not stored in hh_trans_real/complex
      bcast_buffer(:,1) = 0.0_rck
      if (useGPU) then
        successGPU = gpu_memset(bcast_buffer_dev, 0, nbw * size_of_datatype)
        call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1745,  successGPU)

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &real&
             &_gpu_&
             &single&
             &( &
             bcast_buffer_dev, hh_tau_dev, &
             nbw, 1, .true.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU
    endif ! (current_local_n > 1) then

    if (l_nev == 0) cycle

    if (current_local_n > 0) then

      do i = 1, stripe_count

        !wait_b
        if (current_n_end < current_n) then



          if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_recv")
          call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)

          if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_recv")
          n_off = current_local_n+a_off

          if (useGPU) then
            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_memcpy")
              successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                       c_loc(bottom_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                       stripe_width*nbw* size_of_datatype,      &
                                       gpuMemcpyDeviceToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1859,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_memcpy")
            else ! allComputeOnGPU
              if (wantDebug) call obj%timer%start("memcpy")
              dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * size_of_datatype
              successGPU =  gpu_memcpy( aIntern_dev + dev_offset , &
                                      int(loc(bottom_border_recv_buffer(1,i)),kind=c_intptr_t), &
                                       stripe_width*nbw*  size_of_datatype,    &
                                       gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1868,  successGPU)
              if (wantDebug) call obj%timer%stop("memcpy")
            endif ! allComputeOnGPU
          else ! useGPU
            aIntern(:,n_off+1:n_off+nbw,i) = reshape( &
            bottom_border_recv_buffer(1:stripe_width*nbw,i),(/stripe_width,nbw/))
          endif ! useGPU


          if (next_n_end < next_n) then

            if (useGPU) then
              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_REAL4, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            else ! useGPU
              if (wantDebug) call obj%timer%start("mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_REAL4, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("mpi_communication")
            endif ! useGPU

            endif ! (next_n_end < next_n)
          endif ! (current_n_end < current_n)

          if (current_local_n <= bottom_msg_length + top_msg_length) then

            !wait_t
            if (top_msg_length>0) then


              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! the MPI_IRECV will be done CUDA_AWARE we thus do not need a host to device copy
                  ! However, we have to copy from top_border_recv_buffer_mpi_fortran_ptr to aIntern_mpi_fortran_ptr
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  !Fortran pointer for indexing
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                          c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                          stripe_width*top_msg_length* size_of_datatype,      &
                                          gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2045,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  !             host_offset= (0 + (0 * stripe_width) + ( (i-1) * stripe_width * nbw ) ) * 8
                  successGPU =  gpu_memcpy( aIntern_dev+dev_offset , int(loc(top_border_recv_buffer(1,i)),kind=c_intptr_t),  &
                                             stripe_width*top_msg_length* size_of_datatype,      &
                                             gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2054,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
            endif ! top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
                &real&
                &_&
                &single&
                &(obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
                a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
                hh_tau_dev, kernel_flops, kernel_time, n_times, 0, current_local_n, i, &
                last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b        1
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")
            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")

            if (bottom_msg_length>0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer

                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                           c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2248,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")



              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 2255,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy( int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2263,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                      MPI_REAL4, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                      int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                     MPI_REAL4, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                     int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU
            endif !(bottom_msg_length>0)

          else ! current_local_n <= bottom_msg_length + top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &real&
             &_&
             &single&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, &
             current_local_n - bottom_msg_length, bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")

            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")
            if (bottom_msg_length > 0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                            c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                             stripe_width * bottom_msg_length * size_of_datatype,      &
                                             gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2516,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset,  &
                                           stripe_width*bottom_msg_length* size_of_datatype,  &
                                           gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2524,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_REAL4, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
          
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_REAL4, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU


            endif ! (bottom_msg_length > 0)


            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &real&
             &_&
             &single&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, top_msg_length, &
             current_local_n-top_msg_length-bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")


            !wait_t
            if (top_msg_length>0) then

              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
              if (useGPU) then
                if (allComputeOnGPU) then
                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                           c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                           stripe_width* top_msg_length* size_of_datatype,      &
                                           gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2739,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  ! copy top_border_recv_buffer to aIntern_dev, maybe not necessary if CUDA_AWARE IRECV
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(aIntern_dev + dev_offset ,int(loc( top_border_recv_buffer(:,i)),kind=c_intptr_t),  &
                                        stripe_width * top_msg_length * size_of_datatype,   &
                                        gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2748,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
           endif

           !compute

           if (wantDebug) call obj%timer%start("compute_hh_trafo")
           call compute_hh_trafo_&
             &real&
             &_&
             &single&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off, nbw, max_blk_size,  bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, 0, top_msg_length, i, &
             last_stripe_width, kernel)
           if (wantDebug) call obj%timer%stop("compute_hh_trafo")

         endif

         if (next_top_msg_length > 0) then
           !request top_border data

           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_REAL4, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else !useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_REAL4, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif !useGPU


         endif

         !send_t
         if (my_prow > 0) then

           if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
           call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")

           if (useGPU) then
             if (allComputeOnGPU) then
               if (wantDebug) call obj%timer%start("cuda_memcpy")
               successGPU =  gpu_memcpy(c_loc(top_border_send_buffer_mpi_fortran_ptr(1,i)),  &
                                        c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                        stripe_width* nbw* size_of_datatype,      &
                                        gpuMemcpyDeviceToDevice)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 3030,  successGPU)
               if (wantDebug) call obj%timer%stop("cuda_memcpy")
             else ! allComputeOnGPU
               if (wantDebug) call obj%timer%start("memcpy")
               dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
               successGPU =  gpu_memcpy(int(loc(top_border_send_buffer(:,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                         stripe_width*nbw * size_of_datatype, &
                                         gpuMemcpyDeviceToHost)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> top_border_send_buffer", 3038,  successGPU)
               if (wantDebug) call obj%timer%stop("memcpy")
             endif ! allComputeOnGPU
           else ! useGPU
             top_border_send_buffer(:,i) = reshape(aIntern(:,a_off+1:a_off+nbw,i),(/stripe_width*nbw/))
           endif ! useGPU
           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_REAL4, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else ! useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_REAL4, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif ! useGPU

         endif

         ! Care that there are not too many outstanding top_recv_request's
         if (stripe_count > 1) then
           if (i>1) then

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
             call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
           else

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
            call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

          endif
        endif
      enddo ! i = 1, stripe_count

      top_msg_length = next_top_msg_length

    else
      ! wait for last top_send_request

      do i = 1, stripe_count
        if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
        call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")
      enddo
    endif

    ! Care about the result

    if (my_prow == 0) then

      ! topmost process sends nbw rows to destination processes

      do j=0, nfact-1
        num_blk = sweep*nfact+j ! global number of destination block, 0 based
        if (num_blk*nblk >= na) exit

        nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block

        if (wantDebug) call obj%timer%start("host_mpi_wait_result_send")
        call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_result_send")
        dst = mod(num_blk, np_rows)

        if (dst == 0) then
          if (useGPU) then
            row_group_size = min(na - num_blk*nblk, nblk)

            call pack_row_group_&
                 &real&
                 &_gpu_&
                 &single&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
                         row_group(:, :), j * nblk + a_off, row_group_size, &
                         result_buffer_dev, nblk, num_result_buffers, nbuf, .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              ! memcpy DeviceToDevice row_group_dev -> q_dev
              do i = 1, row_group_size  
               if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
               call gpublas_SCOPY(l_nev, c_loc(row_group_mpi_fortran_ptr(1,i)), 1, &
                                           c_loc(q_mpi_fortran_ptr((num_blk / np_rows) * nblk + i,1)), ldq)
               if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
              enddo
            else ! allComputeOnGPU
              do i = 1, row_group_size
                q((num_blk / np_rows) * nblk + i, 1 : l_nev) = row_group(:, i)
              enddo
            endif ! allComputeOnGPU
          else ! useGPU

            do i = 1, min(na - num_blk*nblk, nblk)

              call pack_row_&
                   &real&
                   &_cpu_&
                   &single&
                   &(obj,aIntern, row, j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
              q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
            enddo
          endif ! useGPU

        else ! (dst == 0)

          if (useGPU) then

            call pack_row_group_&
                 &real&
                 &_gpu_&
                 &single&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, &
                   last_stripe_width, a_dim2, l_nev, &
                   result_buffer(:, :, nbuf), j * nblk + a_off, nblk, &
                   result_buffer_dev, nblk, num_result_buffers, nbuf, .true., wantDebug, allComputeOnGPU)

          else  ! useGPU
            do i = 1, nblk
              call pack_row_&
                   &real&
                   &_cpu_&
                   &single&
                   &(obj, aIntern, result_buffer(:,i,nbuf),j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
            enddo
          endif ! useGPU
         
          if (useGPU) then
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL4, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL4, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif ! useGPU
        endif ! (dst == 0)
      enddo  !j=0, nfact-1

    else ! (my_prow == 0)

      ! receive and store final result

      do j = num_bufs_recvd, num_result_blocks-1

        nbuf = mod(j, num_result_buffers) + 1 ! buffer number to get this block

        ! If there is still work to do, just test for the next result request
        ! and leave the loop if it is not ready, otherwise wait for all
        ! outstanding requests

        if (next_local_n > 0) then
            ! needed
            !successGPU = gpu_devicesynchronize()
          if (wantDebug) call obj%timer%start("host_mpi_test")
          call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_test")

          if (.not.flag) exit

        else ! (next_local_n > 0)
          if (wantDebug) call obj%timer%start("host_mpi_wait_result_recv")

          call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_wait_result_recv")
        endif ! (next_local_n > 0)

        ! Fill result buffer into q
        num_blk = j*np_rows + my_prow ! global number of current block, 0 based
        if (useGPU) then
          if (allComputeOnGPU) then
            do i = 1, min(na - num_blk*nblk, nblk)
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_SCOPY(l_nev, c_loc(result_buffer_mpi_fortran_ptr(1,i,nbuf)), 1, &
                                          c_loc(q_mpi_fortran_ptr(j*nblk + i,1)), ldq)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            enddo
          else ! allComputeOnGPU
            do i = 1, min(na - num_blk*nblk, nblk)
              q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
            enddo
          endif ! allComputeOnGPU
        else ! useGPU
          do i = 1, min(na - num_blk*nblk, nblk)
            q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
          enddo
        endif ! useGPU
        ! Queue result buffer again if there are outstanding blocks left

        if (useGPU) then
          if (wantDebug) call obj%timer%start("host_mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL4, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("host_mpi_communication")
        else ! useGPU
          if (wantDebug) call obj%timer%start("mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_REAL4, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif !useGPU

      enddo ! j = num_bufs_recvd, num_result_blocks-1
      num_bufs_recvd = j

    endif ! (my_prow == 0)

    ! Shift the remaining rows to the front of aIntern (if necessary)

    offset = nbw - top_msg_length
    if (offset<0) then
      if (wantDebug) write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                                         &real&
                                         &: internal error, offset for shifting = ',offset
      success = .false.
      return
    endif

    a_off = a_off + offset
    if (a_off + next_local_n + nbw >= a_dim2) then
      do i = 1, stripe_count
        if (useGPU) then
          chunk = min(next_local_n,a_off)

          if (chunk < 1) exit

          if (wantDebug) call obj%timer%start("normal_memcpy")
          do j = top_msg_length+1, top_msg_length+next_local_n, chunk
            this_chunk = min(j+chunk-1,top_msg_length+next_local_n)-j+1
            dev_offset = ((j-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            dev_offset_1 = ((j+a_off-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            num = stripe_width*this_chunk*size_of_datatype
            successGPU = gpu_memcpy(aIntern_dev+dev_offset, aIntern_dev+dev_offset_1, num, gpuMemcpyDeviceToDevice)

            call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> aIntern_dev", 3528,  successGPU)
          end do
          if (wantDebug) call obj%timer%stop("normal_memcpy")
        else ! not useGPU
          do j = top_msg_length+1, top_msg_length+next_local_n
            aIntern(:,j,i) = aIntern(:,j+a_off,i)
          end do
        end if
      end do ! stripe_count

      a_off = 0
    end if
  end do ! sweep
  if (wantDebug) call obj%timer%stop("sweep_loop")

  ! Just for safety:
  if (ANY(top_send_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_send_request ***',my_prow,my_pcol
  if (ANY(bottom_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_send_request ***',my_prow,my_pcol
  if (ANY(top_recv_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_recv_request ***',my_prow,my_pcol
  if (ANY(bottom_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_recv_request ***',my_prow,my_pcol

  if (my_prow == 0) then

    if (wantDebug) call obj%timer%start("host_mpi_waitall")
    call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("host_mpi_waitall")
  endif

  if (ANY(result_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_send_request ***',my_prow,my_pcol
  if (ANY(result_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_recv_request ***',my_prow,my_pcol

  call obj%get("print_flops",print_flops,error)

  if (print_flops == 1) then

    if (useNonBlockingCollectivesRows) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          allreduce_request1, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          allreduce_request3, mpierr)

      call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif

    if (useNonBlockingCollectivesCols) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                          allreduce_request2, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          allreduce_request4, mpierr)

      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                         mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif
  endif

  if (useGPU .and. allComputeOnGPU) then
    ! finally copy q_dev to q
    if (wantDebug) call obj%timer%start("cuda_memcpy")
    successGPU =  gpu_memcpy(int(loc(q(1,1)),kind=c_intptr_t),  &
                             q_dev, &
                             ldq*matrixCols * size_of_datatype, &
                             gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q_dev -> q", 3634,  successGPU)
    if (wantDebug) call obj%timer%stop("cuda_memcpy")

  endif

  if (my_prow==0 .and. my_pcol==0 .and.print_flops == 1) &
      write(error_unit,'(" Kernel time:",f10.3," MFlops: ",es12.5)')  kernel_time, kernel_flops/kernel_time*1.d-6

  ! deallocate all working space

  if (.not.(useGPU)) then
    nullify(aIntern)
    call free(aIntern_ptr)
  endif

  deallocate(row, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: row", 3650,  istat,  errorMessage)

  deallocate(limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: limits", 3653,  istat,  errorMessage)

  deallocate(result_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_send_request", 3656,  istat,  errorMessage)

  deallocate(result_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_recv_request", 3659,  istat,  errorMessage)

  deallocate(result_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_buffer", 3662,  istat,  errorMessage)

  if (useGPU) then
    if (allComputeOnGPU) then
      successGPU = gpu_free(result_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 3667,  successGPU)
      nullify(result_buffer_mpi_fortran_ptr)
    endif

    nullify(bcast_buffer)

    successGPU = gpu_free_host(bcast_buffer_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 3674,  successGPU)
  else
    deallocate(bcast_buffer, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_tridi_to_band: bcast_buffer", 3677,  istat,  errorMessage)
  endif


  if (useGPU) then
    successGPU = gpu_free(aIntern_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 3683,  successGPU)

    if (allComputeOnGPU) then
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: q_dev", 3687,  successGPU)
      nullify(q_mpi_fortran_ptr)

      successGPU = gpu_free(hh_trans_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 3691,  successGPU)
      nullify(hh_trans_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 3695,  successGPU)
      nullify(top_border_recv_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 3699,  successGPU)
      nullify(top_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 3703,  successGPU)
      nullify(bottom_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 3707,  successGPU)
      nullify(bottom_border_recv_buffer_mpi_fortran_ptr)

      nullify(aIntern_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU = gpu_free(hh_tau_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 3714,  successGPU)

    nullify(row_group)

    successGPU = gpu_free_host(row_group_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 3719,  successGPU)

    successGPU = gpu_free(row_group_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 3722,  successGPU)

    if (allComputeOnGPU) then
      nullify(row_group_mpi_fortran_ptr)

      successGPU = gpu_free(row_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_dev", 3728,  successGPU)
      nullify(row_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU =  gpu_free(bcast_buffer_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 3733,  successGPU)

    if (allComputeOnGPU) then
      nullify(bcast_buffer_mpi_fortran_ptr)
    endif

    successGPU = gpu_host_unregister(int(loc(top_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 3740,  successGPU)

    successGPU = gpu_host_unregister(int(loc(top_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3743,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3746,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3749,  successGPU)
  endif ! useGPU

  deallocate(top_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 3753,  istat,  errorMessage)

  deallocate(top_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3756,  istat,  errorMessage)

  deallocate(bottom_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3759,  istat,  errorMessage)

  deallocate(bottom_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3762,  istat,  errorMessage)

  deallocate(top_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_send_request", 3765,  istat,  errorMessage)

  deallocate(top_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_recv_request", 3768,  istat,  errorMessage)

  deallocate(bottom_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_send_request", 3771,  istat,  errorMessage)

  deallocate(bottom_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_recv_request", 3774,  istat,  errorMessage)

  call obj%timer%stop("trans_ev_tridi_to_band_&
                      &real&
                      &" // &
                      &"_single" //&
                      gpuString)
  return

end subroutine

! vim: syntax=fortran



    subroutine band_band_real_&
&single &
                  (obj, na, nb, nbCol, nb2, nb2Col, ab, ab2, d, e, communicator)
    !-------------------------------------------------------------------------------
    ! band_band_real:
    ! Reduces a real symmetric banded matrix to a real symmetric matrix with smaller bandwidth. Householder transformations are not stored.
    ! Matrix size na and original bandwidth nb have to be a multiple of the target bandwidth nb2. (Hint: expand your matrix with
    ! zero entries, if this
    ! requirement doesn't hold)
    !
    !  na          Order of matrix
    !
    !  nb          Semi bandwidth of original matrix
    !
    !  nb2         Semi bandwidth of target matrix
    !
    !  ab          Input matrix with bandwidth nb. The leading dimension of the banded matrix has to be 2*nb. The parallel data layout
    !              has to be accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb+1 to min(na, block_limits(n+1)*nb)
    !              are located on rank n.
    !
    !  ab2         Output matrix with bandwidth nb2. The leading dimension of the banded matrix is 2*nb2. The parallel data layout is
    !              accordant to divide_band(), i.e. the matrix columns block_limits(n)*nb2+1 to min(na, block_limits(n+1)*nb2) are located
    !              on rank n.
    !
    !  d(na)       Diagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
    !
    !  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0, set only if ab2 = 1 (output)
    !
    !  communicator
    !              MPI-Communicator for the total processor set
    !-------------------------------------------------------------------------------
      use elpa_abstract_impl
      use elpa2_workload
      use elpa_blas_interfaces

      use precision
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
      integer(kind=ik), intent(in)               :: na, nb, nbCol, nb2, nb2Col, communicator
      real(kind=rk), intent(inout)               :: ab(2*nb,nbCol) ! removed assumed size
      real(kind=rk), intent(inout)               :: ab2(2*nb2,nb2Col) ! removed assumed size
      real(kind=rk), intent(out)                 :: d(na), e(na) ! set only on PE 0

      real(kind=rk)                              :: hv(nb,nb2), w(nb,nb2), w_new(nb,nb2), tau(nb2), hv_new(nb,nb2), &
                                                  tau_new(nb2), ab_s(1+nb,nb2), ab_r(1+nb,nb2), ab_s2(2*nb2,nb2), hv_s(nb,nb2)

      real(kind=rk)                              :: work(nb*nb2), work2(nb2*nb2)
      integer(kind=ik)                         :: lwork, info
      integer(kind=BLAS_KIND)                  :: infoBLAS

      integer(kind=ik)                         :: istep, i, n, dest
      integer(kind=ik)                         :: n_off, na_s
      integer(kind=ik)                         :: my_pe, n_pes
      integer(kind=MPI_KIND)                   :: my_peMPI, n_pesMPI, mpierr
      integer(kind=ik)                         :: nblocks_total, nblocks
      integer(kind=ik)                         :: nblocks_total2, nblocks2
      integer(kind=MPI_KIND)                   :: ireq_ab, ireq_hv
!      integer(kind=ik)                         :: MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
!      integer(kind=ik), allocatable            :: mpi_statuses(:,:)
      integer(kind=ik), allocatable            :: block_limits(:), block_limits2(:)
      integer(kind=MPI_KIND), allocatable      :: ireq_ab2(:)

      integer(kind=ik)                         :: j, nc, nr, ns, ne, iblk
      integer(kind=ik)                         :: istat
      character(200)                           :: errorMessage

      call obj%timer%start("band_band_real" // "_single")

      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(communicator,kind=MPI_KIND) ,my_peMPI ,mpierr)
      call mpi_comm_size(int(communicator,kind=MPI_KIND) ,n_pesMPI ,mpierr)

      my_pe = int(my_peMPI,kind=c_int)
      n_pes = int(n_pesMPI,kind=c_int)
      call obj%timer%stop("mpi_communication")

      ! Total number of blocks in the band:
      nblocks_total = (na-1)/nb + 1
      nblocks_total2 = (na-1)/nb2 + 1

      ! Set work distribution
      allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating block_limits "//errorMessage
        stop 1
      endif
      call divide_band(obj, nblocks_total, n_pes, block_limits)

      allocate(block_limits2(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating block_limits2 "//errorMessage
        stop 1
      endif

      call divide_band(obj, nblocks_total2, n_pes, block_limits2)

      ! nblocks: the number of blocks for my task
      nblocks = block_limits(my_pe+1) - block_limits(my_pe)
      nblocks2 = block_limits2(my_pe+1) - block_limits2(my_pe)

      allocate(ireq_ab2(1:nblocks2), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating ireq_ab2 "//errorMessage
        stop 1
      endif

      call obj%timer%start("mpi_communication")

      ireq_ab2 = MPI_REQUEST_NULL

      if (nb2>1) then
        do i=0,nblocks2-1

          call mpi_irecv(ab2(1,i*nb2+1), int(2*nb2*nb2,kind=MPI_KIND), MPI_REAL4, &
                         0_MPI_KIND, 3_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_ab2(i+1), mpierr)
        enddo
      endif
      call obj%timer%stop("mpi_communication")

      ! n_off: Offset of ab within band
      n_off = block_limits(my_pe)*nb
      lwork = nb*nb2
      dest = 0
      ireq_ab = MPI_REQUEST_NULL
      ireq_hv = MPI_REQUEST_NULL
      ! ---------------------------------------------------------------------------
      ! Start of calculations

      na_s = block_limits(my_pe)*nb + 1

      if (my_pe>0 .and. na_s<=na) then
        ! send first nb2 columns to previous PE
        ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
        do i=1,nb2
          ab_s(1:nb+1,i) = ab(1:nb+1,na_s-n_off+i-1)
        enddo
        call obj%timer%start("mpi_communication")

        call mpi_isend(ab_s, int((nb+1)*nb2,kind=MPI_KIND), MPI_REAL4, int(my_pe-1,kind=MPI_KIND), &
                       1_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_ab, mpierr)
        call obj%timer%stop("mpi_communication")
      endif

      do istep=1,na/nb2

        if (my_pe==0) then

          n = MIN(na-na_s-nb2+1,nb) ! number of rows to be reduced
          hv(:,:) = 0.0_rk
          tau(:) = 0.0_rk

          ! The last step (istep=na-1) is only needed for sending the last HH vectors.
          ! We don't want the sign of the last element flipped (analogous to the other sweeps)
          if (istep < na/nb2) then

            ! Transform first block column of remaining matrix
      call obj%timer%start("blas")
            call SGEQRF(int(n,kind=BLAS_KIND), int(nb2,kind=BLAS_KIND), ab(1+nb2,na_s-n_off), &
                                 int(2*nb-1,kind=BLAs_KIND), tau, work, int(lwork,kind=BLAS_KIND), &
                                 infoBLAS)
      info = int(infoBLAS,kind=ik)
      call obj%timer%stop("blas")

            do i=1,nb2
              hv(i,i) = 1.0_rk
              hv(i+1:n,i) = ab(1+nb2+1:1+nb2+n-i,na_s-n_off+i-1)
              ab(1+nb2+1:2*nb,na_s-n_off+i-1) = 0.0_rk
            enddo

          endif

          if (nb2==1) then
            d(istep) = ab(1,na_s-n_off)
            e(istep) = ab(2,na_s-n_off)
            if (istep == na) then
              e(na) = 0.0_rk
            endif
          else
            ab_s2 = 0.0_rk
            ab_s2(:,:) = ab(1:nb2+1,na_s-n_off:na_s-n_off+nb2-1)
            if (block_limits2(dest+1)<istep) then
              dest = dest+1
            endif
            call obj%timer%start("mpi_communication")
            call mpi_send(ab_s2, int(2*nb2*nb2,kind=MPI_KIND), MPI_REAL4, int(dest,kind=MPI_KIND), &
                          3_MPI_KIND, int(communicator,kind=MPI_KIND), mpierr)
            call obj%timer%stop("mpi_communication")


          endif

        else
          if (na>na_s+nb2-1) then
            ! Receive Householder vectors from previous task, from PE owning subdiagonal
            call obj%timer%start("mpi_communication")
            call mpi_recv(hv, int(nb*nb2,kind=MPI_KIND), MPI_REAL4, int(my_pe-1,kind=MPI_KIND), &
                          2_MPI_KIND, int(communicator,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            call obj%timer%stop("mpi_communication")


            do i=1,nb2
              tau(i) = hv(i,i)
              hv(i,i) = 1.0_rk
            enddo
          endif
        endif

        na_s = na_s+nb2
        if (na_s-n_off > nb) then
          ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
          ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rk
          n_off = n_off + nb
        endif

        do iblk=1,nblocks
          ns = na_s + (iblk-1)*nb - n_off ! first column in block
          ne = ns+nb-nb2                    ! last column in block

          if (ns+n_off>na) exit

            nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
            nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                          ! Note that nr>=0 implies that diagonal block is full (nc==nb)!
            call wy_gen_&
      &single&
      &(obj,nc,nb2,w,hv,tau,work,nb)

            if (iblk==nblocks .and. nc==nb) then
              !request last nb2 columns
              call obj%timer%start("mpi_communication")
              call mpi_recv(ab_r, int((nb+1)*nb2,kind=MPI_KIND), MPI_REAL4, int(my_pe+1,kind=MPI_KIND), &
                            1_MPI_KIND, int(communicator,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
              call obj%timer%stop("mpi_communication")

              do i=1,nb2
                ab(1:nb+1,ne+i-1) = ab_r(:,i)
              enddo
            endif
            hv_new(:,:) = 0.0_rk ! Needed, last rows must be 0 for nr < nb
            tau_new(:) = 0.0_rk

            if (nr>0) then
              call wy_right_&
        &single&
        &(obj,nr,nb,nb2,ab(nb+1,ns),2*nb-1,w,hv,work,nb)
        call obj%timer%start("blas")
              call SGEQRF(int(nr,kind=BLAS_KIND), int(nb2,kind=BLAS_KIND), ab(nb+1,ns), &
                                   int(2*nb-1,kind=BLAS_KIND), tau_new, work, int(lwork,kind=BLAS_KIND), &
                                   infoBLAS)
        info = int(infoBLAS,kind=ik)
        call obj%timer%stop("blas")
              do i=1,nb2
                hv_new(i,i) = 1.0_rk
                hv_new(i+1:,i) = ab(nb+2:2*nb-i+1,ns+i-1)
                ab(nb+2:,ns+i-1) = 0.0_rk
              enddo

              !send hh-Vector
              if (iblk==nblocks) then
                call obj%timer%start("mpi_communication")

                call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
                call obj%timer%stop("mpi_communication")

                hv_s = hv_new
                do i=1,nb2
                  hv_s(i,i) = tau_new(i)
                enddo
                call obj%timer%start("mpi_communication")
                call mpi_isend(hv_s, int(nb*nb2,kind=MPI_KIND), MPI_REAL4, int(my_pe+1,kind=MPI_KIND), &
                               2_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_hv, mpierr)
                call obj%timer%stop("mpi_communication")

              endif
            endif

            call wy_symm_&
      &single&
      &(obj,nc,nb2,ab(1,ns),2*nb-1,w,hv,work,work2,nb)

            if (my_pe>0 .and. iblk==1) then
              !send first nb2 columns to previous PE
              call obj%timer%start("mpi_communication")

              call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
              call obj%timer%stop("mpi_communication")

              do i=1,nb2
                ab_s(1:nb+1,i) = ab(1:nb+1,ns+i-1)
              enddo
              call obj%timer%start("mpi_communication")
              call mpi_isend(ab_s, int((nb+1)*nb2,kind=MPI_KIND), MPI_REAL4, int(my_pe-1,kind=MPI_KIND), &
                             1_MPI_KIND, int(communicator,kind=MPI_KIND), ireq_ab, mpierr)
              call obj%timer%stop("mpi_communication")

            endif

            if (nr>0) then
              call wy_gen_&
        &single&
        &(obj,nr,nb2,w_new,hv_new,tau_new,work,nb)
              call wy_left_&
        &single&
        &(obj,nb-nb2,nr,nb2,ab(nb+1-nb2,ns+nb2),2*nb-1,w_new,hv_new,work,nb)
            endif

            ! Use new HH Vector for the next block
            hv(:,:) = hv_new(:,:)
            tau = tau_new
          enddo
        enddo

        ! Finish the last outstanding requests
         call obj%timer%start("mpi_communication")

        call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
        call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
!        allocate(mpi_statuses(MPI_STATUS_SIZE,nblocks2), stat=istat, errmsg=errorMessage)
!        if (istat .ne. 0) then
!          print *,"error allocating mpi_statuses "//errorMessage
!          stop 1
!        endif

        call mpi_waitall(nblocks2,ireq_ab2,MPI_STATUSES_IGNORE,mpierr)
!        deallocate(mpi_statuses, stat=istat, errmsg=errorMessage)
!        if (istat .ne. 0) then
!          print *,"error deallocating mpi_statuses "//errorMessage
!          stop 1
!        endif

        call mpi_barrier(int(communicator,kind=MPI_KIND) ,mpierr)
        call obj%timer%stop("mpi_communication")


        deallocate(block_limits, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating block_limits "//errorMessage
          stop 1
        endif

        deallocate(block_limits2, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating block_limits2 "//errorMessage
          stop 1
        endif

        deallocate(ireq_ab2, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating ireq_ab2 "//errorMessage
          stop 1
        endif

        call obj%timer%stop("band_band_real" // "_single")

    end subroutine

    subroutine wy_gen_&
    &single&
    &(obj, n, nb, W, Y, tau, mem, lda)

      use elpa_abstract_impl
      use elpa_blas_interfaces

      use precision
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
      integer(kind=ik), intent(in)            :: n      !length of householder-vectors
      integer(kind=ik), intent(in)            :: nb     !number of householder-vectors
      integer(kind=ik), intent(in)            :: lda        !leading dimension of Y and W
      real(kind=rk), intent(in)               :: Y(lda,nb)  !matrix containing nb householder-vectors of length b
      real(kind=rk), intent(in)               :: tau(nb)    !tau values
      real(kind=rk), intent(out)              :: W(lda,nb)  !output matrix W
      real(kind=rk), intent(in)               :: mem(nb)    !memory for a temporary matrix of size nb

      integer(kind=ik)                        :: i

   call obj%timer%start("wy_gen" // "_single")

   W(1:n,1) = tau(1)*Y(1:n,1)
   do i=2,nb
     W(1:n,i) = tau(i)*Y(1:n,i)
     call obj%timer%start("blas")
     call SGEMV('T', int(n,kind=BLAS_KIND), int(i-1,kind=BLAS_KIND),  1.0_rk, Y, int(lda,kind=BLAS_KIND), &
                         W(1,i), 1_BLAS_KIND, 0.0_rk, mem, 1_BLAS_KIND)
     call SGEMV('N', int(n,kind=BLAS_KIND), int(i-1,kind=BLAS_KIND), -1.0_rk, W, int(lda,kind=BLAS_KIND), &
                         mem, 1_BLAS_KIND, 1.0_rk, W(1,i), 1_BLAS_KIND)
     call obj%timer%stop("blas")
   enddo
   call obj%timer%stop("wy_gen" // "_single")
    end subroutine

    subroutine wy_left_&
    &single&
    &(obj, n, m, nb, A, lda, W, Y, mem, lda2)

      use precision
      use elpa_abstract_impl
      use elpa_blas_interfaces
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
      integer(kind=ik), intent(in)            :: n      !width of the matrix A
      integer(kind=ik), intent(in)            :: m      !length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk), intent(inout)            :: A(lda,*)   !matrix to be transformed   ! remove assumed size
      real(kind=rk), intent(in)               :: W(m,nb)    !blocked transformation matrix W
      real(kind=rk), intent(in)               :: Y(m,nb)    !blocked transformation matrix Y
      real(kind=rk), intent(inout)            :: mem(n,nb)  !memory for a temporary matrix of size n x nb

   call obj%timer%start("wy_left" // "_single")
   call obj%timer%start("blas")
   call SGEMM('T', 'N', int(nb,kind=BLAS_KIND), int(n,kind=BLAS_KIND), int(m,kind=BLAS_KIND), &
                       1.0_rk, W, int(lda2,kind=BLAS_KIND), A, int(lda,kind=BLAS_KIND), 0.0_rk, mem, &
                       int(nb,kind=BLAS_KIND))
   call SGEMM('N', 'N', int(m,kind=BLAS_KIND), int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), &
                       -1.0_rk, Y, int(lda2,kind=BLAS_KIND), mem, int(nb,kind=BLAS_KIND), 1.0_rk, A, int(lda,kind=BLAS_KIND))
   call obj%timer%stop("blas")
   call obj%timer%stop("wy_left" // "_single")
    end subroutine

    subroutine wy_right_&
    &single&
    &(obj, n, m, nb, A, lda, W, Y, mem, lda2)

      use precision
      use elpa_abstract_impl
      use elpa_blas_interfaces
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
      integer(kind=ik), intent(in)            :: n      !height of the matrix A
      integer(kind=ik), intent(in)            :: m      !length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk), intent(inout)            :: A(lda,*)   !matrix to be transformed  ! remove assumed size
      real(kind=rk), intent(in)               :: W(m,nb)    !blocked transformation matrix W
      real(kind=rk), intent(in)               :: Y(m,nb)    !blocked transformation matrix Y
      real(kind=rk), intent(inout)            :: mem(n,nb)  !memory for a temporary matrix of size n x nb


      call obj%timer%start("wy_right" // "_single")
      call obj%timer%start("blas")
      call SGEMM('N', 'N', int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), int(m,kind=BLAS_KIND), &
                          1.0_rk, A, int(lda,kind=BLAS_KIND), W, int(lda2,kind=BLAS_KIND), 0.0_rk, mem, int(n,kind=BLAS_KIND))
      call SGEMM('N', 'T', int(n,kind=BLAS_KIND), int(m,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), &
                          -1.0_rk, mem, int(n,kind=BLAS_KIND), Y, int(lda2,kind=BLAS_KIND), 1.0_rk, A, int(lda,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      call obj%timer%stop("wy_right" // "_single")

    end subroutine

    subroutine wy_symm_&
    &single&
    &(obj, n, nb, A, lda, W, Y, mem, mem2, lda2)

      use elpa_abstract_impl
      use elpa_blas_interfaces

      use precision
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
      integer(kind=ik), intent(in)            :: n      !width/heigth of the matrix A; length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk), intent(inout)            :: A(lda,*)   !matrix to be transformed  ! remove assumed size
      real(kind=rk), intent(in)               :: W(n,nb)    !blocked transformation matrix W
      real(kind=rk), intent(in)               :: Y(n,nb)    !blocked transformation matrix Y
      real(kind=rk)                           :: mem(n,nb)  !memory for a temporary matrix of size n x nb
      real(kind=rk)                           :: mem2(nb,nb)    !memory for a temporary matrix of size nb x nb

      call obj%timer%start("wy_symm" // "_single")
      call obj%timer%start("blas")
      call SSYMM('L', 'L', int(n, kind=BLAS_KIND), int(nb,kind=BLAS_KIND), 1.0_rk, A, &
                          int(lda,kind=BLAS_KIND), W, int(lda2,kind=BLAS_KIND), 0.0_rk, mem, int(n,kind=BLAS_KIND))
      call SGEMM('T', 'N', int(nb,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), int(n,kind=BLAS_KIND), &
                          1.0_rk, mem, int(n,kind=BLAS_KIND), W, int(lda2,kind=BLAS_KIND), 0.0_rk, mem2, &
                          int(nb,kind=BLAS_KIND))
      call SGEMM('N', 'N', int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), &
                          -0.5_rk, Y, int(lda2,kind=BLAS_KIND), mem2, int(nb,kind=BLAS_KIND), 1.0_rk, mem, int(n,kind=BLAS_KIND))
      call SSYR2K('L', 'N',int(n,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -1.0_rk, Y, int(lda2,kind=BLAS_KIND), &
                           mem, int(n,kind=BLAS_KIND), 1.0_rk, A, int(lda,kind=BLAS_KIND))
      call obj%timer%stop("blas")
      call obj%timer%stop("wy_symm" // "_single")

    end subroutine



! complex double precision























subroutine bandred_&
&complex&
&_&
&double &
(obj, na, a_mat, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, tmat, &
wantDebug, useGPU, success, &
max_threads, isSkewsymmetric)

!-------------------------------------------------------------------------------
!  bandred_real/complex: Reduces a distributed symmetric matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a_mat(lda,matrixCols)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a_mat(:,:) must be set completely (upper and lower half)
!              a_mat(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,numBlocks)    where numBlocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa1_compute
  use precision
  use elpa_blas_interfaces
  use elpa_scalapack_interfaces
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
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                            :: na, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols

  complex(kind=rck)                     :: a_mat(lda,*)
  complex(kind=rck)                     :: tmat(nbw,nbw,*)

  logical, intent(in)                         :: useGPU
  logical, intent(in)                         :: isSkewsymmetric
  character(20)                               :: gpuString

  integer(kind=ik)                            :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                      :: mpierr,  my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                            :: l_cols, l_rows
  integer(kind=ik)                            :: i, j, lcs, lce, lre, lc, lr, cur_pcol, n_cols, nrow
  integer(kind=ik)                            :: istep, ncol, lch, lcx, nlc
  integer(kind=ik)                            :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rk)                              :: vnorm2
  complex(kind=rck)                    :: xf, aux1(nbw), aux2(nbw), vrl, tau
  complex(kind=rck)                    :: vav(nbw,nbw)

  complex(kind=rck), allocatable :: tmpGPU(:)
  complex(kind=rck), pointer     :: vmrGPU(:), umcGPU(:)
  complex(kind=rck), allocatable :: tmpCPU(:,:), vmrCPU(:,:), umcCPU(:,:)
  complex(kind=rck), allocatable :: vr(:)

  integer(kind=C_intptr_T)                    :: a_dev, vmr_dev, umc_dev, tmat_dev, vav_dev
  type(c_ptr)                                 :: vmr_host, umc_host
  !integer(kind=ik), external                  :: numroc -> use elpa_scalapack
  integer(kind=ik)                            :: ierr
  integer(kind=ik)                            :: cur_l_rows, cur_l_cols, vmr_size, umc_size
  integer(kind=ik)                            :: l_rows2, vmr_size2, umc_size2
  integer(kind=c_intptr_t)                    :: lc_start, lc_end
  integer(kind=c_intptr_t)                    :: lce_1, lcs_1, lre_1
  integer(kind=ik)                            :: lr_end
  integer(kind=ik)                            :: na_cols
  integer(kind=BLAS_KIND)                     :: na_colsBLAS
  integer(kind=ik)                            :: na_rows
  integer(kind=BLAS_KIND)                     :: na_rowsBLAS

  logical, intent(in)                         :: wantDebug
  logical, intent(out)                        :: success
  logical                                     :: successGPU
  integer(kind=ik)                            :: istat
  character(200)                              :: errorMessage
  integer(kind=ik)                            :: min_tile_size, error

  integer(kind=ik)                            :: mystart, myend, m_way, n_way, work_per_thread, m_id, n_id, n_threads, &
                                                ii, pp
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                    &double&
                                                                    &_&
                                                                    &complex

  logical                                     :: useGPU_reduction_lower_block_to_tridiagonal
  integer(kind=ik), intent(in)                :: max_threads
  logical                                     :: do_memcpy
  integer(kind=ik)                            :: i_blk,blk_off, blk_end
  logical                                     :: useIntelGPU


  integer(kind=MPI_KIND)                      :: bcast_request, allreduce_request1, allreduce_request2, &
                                                 allreduce_request3, allreduce_request4, allreduce_request5, &
                                                 allreduce_request6

  logical                                     :: useNonBlockingCollectivesCols
  logical                                     :: useNonBlockingCollectivesRows
  integer(kind=c_int)                         :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("bandred_&
  &complex&
  &" // &
  "_double" // &
  gpuString )

  call obj%get("nbc_row_elpa2_full_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_bandred. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_full_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_bandred. Aborting..."
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

  useGPU_reduction_lower_block_to_tridiagonal = .false.

  if (useGPU .and. .not.(useIntelGPU)) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
  endif

  if (useIntelGPU) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")
  success = .true.


  ! Semibandwith nbw must be a multiple of blocksize nblk
  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_bandred_&
                             &complex&
                             &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_bandred_&
                             &complex&
                             &: ELPA2 works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  ! na_rows in used nowhere; only na_cols
  if (useGPU .and. .not.(useIntelGPU)) then
    na_rowsBLAS = numroc(int(na,kind=BLAS_KIND), int(nblk,kind=BLAS_KIND), &
                         int(my_prow,kind=BLAS_KIND), 0_BLAS_KIND, int(np_rows,kind=BLAS_KIND))
    na_rows = int(na_rowsBLAS,kind=c_int)
    na_colsBLAS = numroc(int(na,kind=BLAS_KIND), int(nblk,kind=BLAS_KIND), &
                         int(my_pcol,kind=BLAS_KIND), 0_BLAS_KIND, int(np_cols,kind=BLAS_KIND))
    na_cols = int(na_colsBLAS,kind=c_int)

    ! Here we convert the regular host array into a pinned host array
    successGPU = gpu_malloc(a_dev, lda*na_cols* size_of_datatype)
    call check_alloc_GPU_f("bandred: a_dev", 335,  successGPU)

    successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), &
                  nbw * nbw * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: vav", 340,  successGPU)

    successGPU = gpu_malloc(vav_dev, nbw*nbw* size_of_datatype)
    call check_alloc_GPU_f("bandred: vav_dev", 343,  successGPU)
  endif ! useGPU

  !if (useIntelGPU) then
  !  ! needed later when explicit copy
  !endif ! useIntelGPU

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

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

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile


  blk_end = (na-1)/nbw
  if (useGPU .and. .not.(useIntelGPU)) then

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: a_mat", 422,  successGPU)

    cur_l_rows = 0
    cur_l_cols = 0

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("bandred: a_dev", 429,  successGPU)

    successGPU = gpu_malloc(tmat_dev, nbw*nbw*size_of_datatype)
    call check_alloc_GPU_f("bandred: tmat_dev", 432,  successGPU)

    istep = (na-1)/nbw
    blk_end = (na-1)/nbw
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows,1)
    cur_l_cols = max(l_cols,1)
    vmr_size = cur_l_rows*2*n_cols
    umc_size = cur_l_cols*2*n_cols

    istep = (na-1)/nbw - 1
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows2 = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows2,1)
    cur_l_cols = max(l_cols,1)
    vmr_size2 = cur_l_rows*2*n_cols
    umc_size2 = cur_l_cols*2*n_cols

    l_rows = max(l_rows,l_rows2)
    vmr_size = max(vmr_size,vmr_size2)
    umc_size = max(umc_size,umc_size2)

    allocate(vr(l_rows + 1), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      print *,"bandred_&
              &complex&
              &: error when allocating vr "//errorMessage
      stop 1
    endif

    successGPU = gpu_malloc_host(vmr_host,vmr_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: vmr_host", 466,  successGPU)
    call c_f_pointer(vmr_host, vmrGPU, (/vmr_size/))

    successGPU = gpu_malloc(vmr_dev, vmr_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: vmr_dev", 470,  successGPU)

    successGPU = gpu_malloc_host(umc_host,umc_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: umc_host", 473,  successGPU)
    call c_f_pointer(umc_host, umcGPU, (/umc_size/))

    successGPU = gpu_malloc(umc_dev, umc_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: umc_dev", 477,  successGPU)

  endif ! useGPU

  !if (useIntelGPU) then
     ! needed later when explict mem copy
  !endif ! useIntelGPU


  do istep = blk_end, 1, -1

    n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

    ! Number of local columns/rows of remaining matrix
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

    ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

    if (useGPU) then
      if (useIntelGPU) then
        ! unify the the name vmr and vmrCPU, as well as vmrGPU
        ! the same for umcCPU and umcGPU
        ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

        allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vmrCPU", 503,  istat,  errorMessage)

        allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: umcCPU", 506,  istat,  errorMessage)

        allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vr", 509,  istat,  errorMessage)
      else
        cur_l_rows = max(l_rows, 1)
        cur_l_cols = max(l_cols, 1)
        vmr_size = cur_l_rows * 2 * n_cols
        umc_size = cur_l_cols * 2 * n_cols
      endif
    else ! GPU not used

      ! unify the the name vmr and vmrCPU, as well as vmrGPU
      ! the same for umcCPU and umcGPU
      ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

      allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vmrCPU", 523,  istat,  errorMessage)

      allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: umcCPU", 526,  istat,  errorMessage)

      allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vr", 529,  istat,  errorMessage)

    endif ! use GPU

    if (useGPU) then
      if (useIntelGPU) then
        vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
      else
        vmrGPU(1 : cur_l_rows * n_cols) = 0.0_rck
        umcGPU(1 : umc_size) = 0.0_rck
      endif
    else
      vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
    endif ! useGPU

    vr(:) = 0.0_rck
    tmat(:,:,istep) = 0.0_rck
    if (useGPU .and. .not.(useIntelGPU)) then
      lc_start = local_index(istep*nbw+1, my_pcol, np_cols, nblk, -1)
      lc_end   = local_index(istep*nbw+n_cols, my_pcol, np_cols, nblk, -1)
      lr_end   = local_index((istep-1)*nbw + n_cols, my_prow, np_rows, nblk, -1)

      if (lc_start .le. 0) lc_start = 1

      do_memcpy = .false.

      ! Note: mod(nbw,nblk) == 0
      do i_blk = 1, nbw/nblk
        blk_off = (i_blk-1) * nblk
        cur_pcol = pcol(istep*nbw+1+blk_off, nblk, np_cols)

        if (my_pcol == cur_pcol) then
          do_memcpy = .true.
        endif
      enddo

      if (do_memcpy) then
        successGPU = gpu_memcpy2d(int(loc(a_mat(1, lc_start)),kind=c_intptr_t), &
                      int((lda*size_of_datatype),kind=c_intptr_t), &
                      (a_dev + int( ( (lc_start-1) * lda*size_of_datatype),kind=c_intptr_t )), &
                      int(lda*size_of_datatype,kind=c_intptr_t), &
                      int(lr_end*size_of_datatype,kind=c_intptr_t), &
                      int((lc_end - lc_start+1),kind=c_intptr_t),int(gpuMemcpyDeviceToHost,kind=c_int))

        call check_memcpy_GPU_f("bandred: a_dev -> a_mat", 573,  successGPU)
      endif
    endif ! useGPU

    !if (useIntelGPU) then
       ! needed later when explict mem copy
    !endif ! useIntelGPU


    ! Reduce current block to lower triangular form
      do lc = n_cols, 1, -1

        ncol = istep*nbw + lc ! absolute column number of householder Vector
        nrow = ncol - nbw ! Absolute number of pivot row

        lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
        lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

        tau = 0

        if (nrow == 1) exit ! Nothing to do

        cur_pcol = pcol(ncol, nblk, np_cols) ! Processor column owning current block

        if (my_pcol==cur_pcol) then

          ! Get Vector to be transformed; distribute last element and norm of
          ! remaining elements to all procs in current column

          vr(1:lr) = a_mat(1:lr,lch) ! Vector to be transformed

          if (my_prow==prow(nrow, nblk, np_rows)) then
            aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
            aux1(2) = vr(lr)
          else
            aux1(1) = dot_product(vr(1:lr),vr(1:lr))
            aux1(2) = 0.0_rck
          endif

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_DOUBLE_COMPLEX, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             allreduce_request1, mpierr)

            call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)

            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_DOUBLE_COMPLEX, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             mpierr)

            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif


          vnorm2 = real(aux2(1),kind=rk)
          vrl    = aux2(2)

          ! Householder transformation
          call hh_transform_&
             &complex&
             &_&
             &double &
                         (obj, vrl, vnorm2, xf, tau, wantDebug)
          ! Scale vr and store Householder Vector for back transformation

          vr(1:lr) = vr(1:lr) * xf
          if (my_prow==prow(nrow, nblk, np_rows)) then
            a_mat(1:lr-1,lch) = vr(1:lr-1)
            a_mat(lr,lch) = vrl
            vr(lr) = 1.0_rck
          else
            a_mat(1:lr,lch) = vr(1:lr)
          endif

        endif

        ! Broadcast Householder Vector and tau along columns

        vr(lr+1) = tau
        if (useNonBlockingCollectivesCols) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(vr, int(lr+1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        bcast_request, mpierr)


          call mpi_wait(bcast_request, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_bcast(vr, int(lr+1,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
          vmrGPU(cur_l_rows * (lc - 1) + 1 : cur_l_rows * (lc - 1) + lr) = vr(1:lr)
        else
          vmrCPU(1:lr,lc) = vr(1:lr)
        endif
        tau = vr(lr+1)

        tmat(lc,lc,istep) = conjg(tau) ! Store tau in diagonal of tmat
        ! Transform remaining columns in current block with Householder Vector
        ! Local dot product

        aux1 = 0.0_rck


        nlc = 0 ! number of local columns
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            if (lr>0) aux1(nlc) = dot_product(vr(1:lr),a_mat(1:lr,lcx))
          endif
        enddo

        ! Get global dot products
        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          if (nlc>0) then
            call mpi_iallreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      allreduce_request3, mpierr)
            call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          if (nlc>0) then
            call mpi_allreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif
        ! Transform

        nlc = 0
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            a_mat(1:lr,lcx) = a_mat(1:lr,lcx) - conjg(tau)*aux2(nlc)*vr(1:lr)
          endif
        enddo
      enddo ! lc

      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        ! store column tiles back to GPU
        if (do_memcpy) then
          successGPU = gpu_memcpy2d((a_dev+ &
                        int(((lc_start-1)*lda*size_of_datatype),kind=c_intptr_t)), &
                        int(lda*size_of_datatype,kind=c_intptr_t), int(loc(a_mat(1,lc_start)),kind=c_intptr_t), &
                        int(lda*size_of_datatype,kind=c_intptr_t), &
                        int(lr_end*size_of_datatype,kind=c_intptr_t), &
                        int((lc_end - lc_start+1),kind=c_intptr_t), &
                        int(gpuMemcpyHostToDevice,kind=c_int))
          call check_memcpy_GPU_f("bandred: a_mat -> a_dev", 937,  successGPU)
        endif
      endif

      if (useGPU_reduction_lower_block_to_tridiagonal .and. useIntelGPU) then
        ! store column tiles back to GPU
      endif

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use dsyrk

      vav = 0
      call obj%timer%start("blas")
      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        if (l_rows>0) &
        call ZHERK('U', 'C',            &
                           int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                           vmrGPU, int(cur_l_rows,kind=BLAS_KIND), &
                           ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      else ! useGPU_reduction_to_tridiagonal
        if (l_rows>0) &
        call ZHERK('U', 'C',           &
                            int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, vmrCPU, &
                            int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      call herm_matrix_allreduce_&
         &double &
                         (obj, n_cols,vav, nbw, nbw,mpi_comm_rows, .true.)
         ! Calculate triangular matrix T for block Householder Transformation
      call obj%timer%start("blas")
      do lc=n_cols,1,-1
        tau = tmat(lc,lc,istep)
        if (lc<n_cols) then
          call ZTRMV('U', 'C', 'N',&
                              int(n_cols-lc,kind=BLAS_KIND), tmat(lc+1,lc+1,istep), &
                              int(ubound(tmat,dim=1),kind=BLAS_KIND), vav(lc+1,lc), 1_BLAS_KIND)

          tmat(lc,lc+1:n_cols,istep) = -tau * conjg(vav(lc+1:n_cols,lc))
        endif
      enddo
      call obj%timer%stop("blas")


    ! Transpose vmr -> vmc (stored in umc, second half)
    if (useGPU) then
      if (useIntelGPU) then
        call elpa_transpose_vectors_&
             &complex&
             &_&
             &double &
                                          (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                           umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                           1, istep*nbw, n_cols, nblk, max_threads, .true.)

      else
        call elpa_transpose_vectors_&
             &complex&
             &_&
             &double &
                          (obj, vmrGPU(:), cur_l_rows, mpi_comm_rows, &
                           umcGPU(cur_l_cols * n_cols + 1:), cur_l_cols, &
                           mpi_comm_cols, 1, istep*nbw, n_cols, nblk, max_threads, .true.)
      endif
    else ! useGPU
      call elpa_transpose_vectors_&
           &complex&
           &_&
           &double &
                                        (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                         umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                         1, istep*nbw, n_cols, nblk, max_threads, .true.)
    endif

    ! Calculate umc = A**T * vmr
    ! Note that the distributed A has to be transposed
    ! Opposed to direct tridiagonalization there is no need to use the cache locality
    ! of the tiles, so we can use strips of the matrix


!#if 0
!    ! original complex implemetation check for performance
!    umcCPU(1:l_cols,1:n_cols) = 0.0_rck
!    vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
!
!    if (l_cols>0 .and. l_rows>0) then
!      do i=0,(istep*nbw-1)/tile_size
!
!        lcs = i*l_cols_tile+1
!        lce = min(l_cols,(i+1)*l_cols_tile)
!        if (lce<lcs) cycle
!
!        lre = min(l_rows,(i+1)*l_rows_tile)
!
!          call obj%timer%start("blas")
!          call ZGEMM('C', 'N', lce-lcs+1, n_cols, lre, ONE, a_mat(1,lcs), ubound(a_mat,dim=1), &
!                     vmrCPU, ubound(vmrCPU,dim=1), ONE, umcCPU(lcs,1), ubound(umcCPU,dim=1))
!          call obj%timer%stop("blas")
!
!        if (i==0) cycle
!        lre = min(l_rows,i*l_rows_tile)
!          call obj%timer%start("blas")
!          call ZGEMM('N', 'N', lre, n_cols, lce-lcs+1, ONE, a_mat(1,lcs), lda, &
!                     umcCPU(lcs,n_cols+1), ubound(umcCPU,dim=1), ONE, vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1))
!          call obj%timer%stop("blas")
!      enddo
!
!    endif ! (l_cols>0 .and. l_rows>0)
!#endif 

    !Code for Algorithm 4

    ! n_way is actually a branch for the number of OpenMP threads
    n_way = 1

      if (.not. useGPU .or. useIntelGPU) then
        umcCPU(1:l_cols,1:n_cols) = 0.0_rck
        vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
      endif ! useGPU

      if (l_cols>0 .and. l_rows>0) then

        if (useGPU .and. .not.(useIntelGPU)) then
          successGPU = gpu_memset(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                        0, cur_l_rows*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: vmr_dev", 1193,  successGPU)

          successGPU = gpu_memcpy(vmr_dev, int(loc(vmrGPU(1)),kind=c_intptr_t), &
                        cur_l_rows*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: vmrGPU -> vmr_dev", 1197,  successGPU)

          successGPU = gpu_memset(umc_dev, 0, l_cols*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: umc_dev", 1200,  successGPU)

          successGPU = gpu_memcpy(umc_dev+l_cols*n_cols*size_of_datatype, &
                        int(loc(umcGPU(1+l_cols*n_cols)),kind=c_intptr_t), &
                        (umc_size-l_cols*n_cols)*size_of_datatype, &
                        gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev", 1206,  successGPU)
        endif ! useGPU

        do i=0,(istep*nbw-1)/tile_size

          lcs = i*l_cols_tile+1
          lce = min(l_cols,(i+1)*l_cols_tile)
          if (lce<lcs) cycle
          lre = min(l_rows,(i+1)*l_rows_tile)

          if (useGPU) then
            if (useIntelGPU) then
              call obj%timer%start("mkl_offload")
              call obj%timer%stop("mkl_offload")
              if (i==0) cycle
              lre = min(l_rows,i*l_rows_tile)
              call obj%timer%start("mkl_offload")

              if (isSkewsymmetric) then

              else
              endif
              call obj%timer%stop("mkl_offload")

            else
              call obj%timer%start("gpublas")
              call gpublas_ZGEMM('C', 'N',                   &
                                         lce-lcs+1, n_cols, lre,     &
                                         ONE, (a_dev + ((lcs-1)*lda* &
                                         size_of_datatype)),         &
                                         lda, vmr_dev,cur_l_rows,    &
                                         ONE, (umc_dev+ (lcs-1)*     &
                                             size_of_datatype),      &
                                         cur_l_cols)

              call obj%timer%stop("gpublas")

              if(i==0) cycle
              call obj%timer%start("gpublas")

              lre = min(l_rows,i*l_rows_tile)
              if (isSkewsymmetric) then
                call gpublas_ZGEMM('N', 'N', lre,n_cols, lce-lcs+1, -ONE, &
                              (a_dev+ ((lcs-1)*lda*                 &
                                    size_of_datatype)),             &
                         lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                size_of_datatype),              &
                                cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                              size_of_datatype),              &
                                cur_l_rows)
              else
                call gpublas_ZGEMM('N', 'N', lre,n_cols, lce-lcs+1, ONE, &
                                            (a_dev+ ((lcs-1)*lda*                 &
                                                  size_of_datatype)),             &
                                       lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                              size_of_datatype),              &
                                              cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                                            size_of_datatype),              &
                                              cur_l_rows)
              endif
              call obj%timer%stop("gpublas")
            endif
          else ! useGPU

            call obj%timer%start("blas")
            call ZGEMM('C', 'N',       &
                                int(lce-lcs+1,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lre,kind=BLAS_KIND), &
                                ONE, a_mat(1,lcs), int(ubound(a_mat,dim=1),kind=BLAS_KIND), &
                                vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ONE, umcCPU(lcs,1), &
                                int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
            call obj%timer%stop("blas")
            if (i==0) cycle
            lre = min(l_rows,i*l_rows_tile)
            call obj%timer%start("blas")

            if (isSkewsymmetric) then
              call ZGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  -ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )

            else
              call ZGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )
            endif
            call obj%timer%stop("blas")
          endif ! useGPU
        enddo ! i=0,(istep*nbw-1)/tile_size

        if (useGPU .and. .not.(useIntelGPU)) then
          if (tile_size < istep*nbw .or. n_way > 1) then
            successGPU = gpu_memcpy(int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                          vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                          (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyDeviceToHost)
            call check_memcpy_GPU_f("bandred: vmr_dev -> vmrGPU", 1342,  successGPU)
          endif

          successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                        umc_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU", 1347,  successGPU)
        endif ! useGPU
      endif ! l_cols>0 .and. l_rows>0

    ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if ur has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix

    ! Or if we used the Algorithm 4
    if (tile_size < istep*nbw .or. n_way > 1) then

      if (useGPU) then
        if (useIntelGPU) then
          call elpa_reduce_add_vectors_&
          &complex&
          &_&
          &double &
                                      (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                       umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                      istep*nbw, n_cols, nblk, max_threads)

        else

          call elpa_reduce_add_vectors_&
               &complex&
               &_&
               &double &
                               (obj, vmrGPU(cur_l_rows * n_cols + 1:),cur_l_rows,  &
                                mpi_comm_rows, umcGPU,                            &
                                cur_l_cols, mpi_comm_cols, istep*nbw, n_cols, nblk, max_threads)
        endif
      else ! useGPU

        call elpa_reduce_add_vectors_&
        &complex&
        &_&
        &double &
                                         (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                          umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          istep*nbw, n_cols, nblk, max_threads)
      endif ! useGPU
    endif ! tile_size < istep*nbw .or. n_way > 1

    if (l_cols>0) then

      if (useGPU) then
        if (useIntelGPU) then
          allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpCPU", 1399,  istat,  errorMessage)

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           allreduce_request4, mpierr)
            call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("bandred: tmpCPU", 1421,  istat,  errorMessage)
        else
          allocate(tmpGPU(l_cols * n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpGPU", 1425,  istat,  errorMessage)
          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")

            call mpi_iallreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request5, mpierr)
            call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")

            call mpi_allreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          if (allocated(tmpGPU)) then
            deallocate(tmpGPU, stat=istat, errmsg=errorMessage)
            call check_deallocate_f("bandred: tmpGPU", 1447,  istat,  errorMessage)
          endif
        endif

      else ! useGPU

        allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: tmpCPU", 1454,  istat,  errorMessage)

        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6, mpierr)
          call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: tmpCPU", 1474,  istat,  errorMessage)
      endif ! useGPU
    endif ! l_cols > 0

    ! U = U * Tmat**T

    if (useGPU) then
      if (useIntelGPU) then
        call obj%timer%start("mkl_offload")

        call ZTRMM('Right', 'Upper', 'C', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

        call ZGEMM('C', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

        call ZTRMM('Right', 'Upper', 'C', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
        call obj%timer%stop("mkl_offload")

      else
        successGPU = gpu_memcpy(umc_dev, int(loc(umcGPU(1)),kind=c_intptr_t), &
                      l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev ", 1526,  successGPU)

        successGPU = gpu_memcpy(tmat_dev,int(loc(tmat(1,1,istep)),kind=c_intptr_t), &
                      nbw*nbw*size_of_datatype,gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: tmat -> tmat_dev ", 1530,  successGPU)

        call obj%timer%start("gpublas")
        call gpublas_ZTRMM('Right', 'Upper', 'C', 'Nonunit',  &
                            l_cols, n_cols, ONE, tmat_dev, nbw, umc_dev, cur_l_cols)
        call obj%timer%stop("gpublas")

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T
        call obj%timer%start("gpublas")
        call gpublas_ZGEMM('C', 'N',             &
                                 n_cols, n_cols, l_cols, ONE, umc_dev, cur_l_cols, &
                                 (umc_dev+(cur_l_cols * n_cols )*size_of_datatype),cur_l_cols, &
                                 ZERO, vav_dev, nbw)

        call gpublas_ZTRMM('Right', 'Upper', 'C', 'Nonunit',    &
           n_cols, n_cols, ONE, tmat_dev, nbw, vav_dev, nbw)
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(vav),kind=c_intptr_t), &
                    vav_dev, nbw*nbw*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: vav_dev -> vav ", 1550,  successGPU)
      endif
    else ! useGPU

      call obj%timer%start("blas")

      call ZTRMM('Right', 'Upper', 'C', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      call ZGEMM('C', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      call ZTRMM('Right', 'Upper', 'C', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
      call obj%timer%stop("blas")

    endif ! useGPU

    call herm_matrix_allreduce_&
         &double &
         (obj, n_cols,vav, nbw, nbw ,mpi_comm_cols, .false.)

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(vav_dev, int(loc(vav),kind=c_intptr_t), &
                       nbw*nbw*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("bandred: vav -> vav_dev ", 1599,  successGPU)
    endif
    !if (useIntelGPU) then
      ! needed later
    !endif


    ! U = U - 0.5 * V * VAV

    if (useGPU) then
      if (useIntelGPU) then
      call obj%timer%start("mkl_offload")
      call ZGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                         (-0.5_rk, 0.0_rk),     &
                         umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav, &
                         int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      call obj%timer%stop("mkl_offload")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &complex&
        &_&
        &double &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &complex&
       &_&
       &double &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif

      else
        call obj%timer%start("gpublas")
        if (isSkewsymmetric) then
          call gpublas_ZGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                    (0.5_rk, 0.0_rk), &
                                    (umc_dev+(cur_l_cols * n_cols )* &
                                    size_of_datatype),   &
                                    cur_l_cols, vav_dev,nbw,        &
                                    ONE, umc_dev, cur_l_cols)
        else
          call gpublas_ZGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                   (-0.5_rk, 0.0_rk), &
                                   (umc_dev+(cur_l_cols * n_cols )* &
                                   size_of_datatype),   &
                                   cur_l_cols, vav_dev,nbw,        &
                                   ONE, umc_dev, cur_l_cols)
        endif
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                    umc_dev, umc_size*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU ", 1681,  successGPU)

        ! Transpose umc -> umr (stored in vmr, second half)
        if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
             &complex&
             &_&
             &double &
                         (obj, umcGPU(:), cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        else
          call elpa_transpose_vectors_&
             &complex&
             &_&
             &double &
                         (obj, umcGPU, cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        endif

        successGPU = gpu_memcpy(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                    int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                    (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: vmr -> vmrGPU ", 1705,  successGPU)
      endif

    else ! useGPU
      call obj%timer%start("blas")
      call ZGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                         (-0.5_rk, 0.0_rk),     &
                         umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav, &
                         int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      call obj%timer%stop("blas")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &complex&
        &_&
        &double &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &complex&
       &_&
       &double &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif
    endif  ! useGPU

    ! A = A - V*U**T - U*V**T


    do i=0,(istep*nbw-1)/tile_size
      lcs = i*l_cols_tile+1
      lce = min(l_cols,(i+1)*l_cols_tile)
      lre = min(l_rows,(i+1)*l_rows_tile)
      if (lce<lcs .or. lre<1) cycle

      if (useGPU) then
        if (useIntelGPU) then
          call obj%timer%start("mkl_offload")
          call ZGEMM('N', 'C', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                              int(2*n_cols,kind=BLAS_KIND), &
                              -ONE, &
                              vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                              int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                              ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
          call obj%timer%stop("mkl_offload")

        else
          call obj%timer%start("gpublas")

          call gpublas_ZGEMM('N', 'C',     &
                                     lre, lce-lcs+1, 2*n_cols, -ONE, &
                                     vmr_dev, cur_l_rows, (umc_dev +(lcs-1)*  &
                                     size_of_datatype), &
                                     cur_l_cols, ONE, (a_dev+(lcs-1)*lda* &
                                     size_of_datatype), lda)
          call obj%timer%stop("gpublas")
        endif
      else ! useGPU

        call obj%timer%start("blas")
        call ZGEMM('N', 'C', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                            int(2*n_cols,kind=BLAS_KIND), &
                            -ONE, &
                            vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                            int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                            ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    enddo ! i=0,(istep*nbw-1)/tile_size

    if (.not.(useGPU) .or. useIntelGPU) then
      if (allocated(vr)) then
        deallocate(vr, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vr", 1871,  istat,  errorMessage)
      endif

      if (allocated(umcCPU)) then
        deallocate(umcCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: umcCPU", 1876,  istat,  errorMessage)
      endif

      if (allocated(vmrCPU)) then
        deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vmrCPU", 1881,  istat,  errorMessage)
      endif
    endif !useGPU

  enddo ! istep - loop

  if (useGPU .and. .not.(useIntelGPU)) then
    ! copy a_dev to a_mat
    ! we do it here, since a is needed on the host in the following routine
    ! (band to tridi). Previously, a has been kept on the device and then
    ! copied in redist_band (called from tridiag_band). However, it seems to
    ! be easier to do it here.
    successGPU = gpu_memcpy(int(loc(a_mat),kind=c_intptr_t), &
                  int(a_dev,kind=c_intptr_t), &
                  int(lda*matrixCols* size_of_datatype, kind=c_intptr_t), &
                  gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("bandred: a_dev -> a_mat ", 1897,  successGPU)

    successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: a_mat ", 1900,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("bandred: a_dev ", 1903,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("bandred: vav_dev ", 1906,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("bandred: tmat_dev ", 1909,  successGPU)

    successGPU = gpu_host_unregister(int(loc(vav),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: vav", 1912,  successGPU)

    if (associated(umcGPU)) then
      nullify(umcGPU)

      successGPU = gpu_free_host(umc_host)
      call check_host_dealloc_GPU_f("bandred: umc_host ", 1918,  successGPU)

      successGPU = gpu_free(umc_dev)
      call check_dealloc_GPU_f("bandred: umc_dev ", 1921,  successGPU)
    endif

    if (associated(vmrGPU)) then
      nullify(vmrGPU)

      successGPU = gpu_free_host(vmr_host)
      call check_host_dealloc_GPU_f("bandred: vmr_host ", 1928,  successGPU)

      successGPU = gpu_free(vmr_dev)
      call check_dealloc_GPU_f("bandred: vmr_dev ", 1931,  successGPU)
    endif
  endif ! useGPU

  !if (useIntelGPU) then
  !   ! needed later
  !endif

  if (allocated(vr)) then
    deallocate(vr, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vr", 1941,  istat,  errorMessage)
  endif

  if (allocated(umcCPU)) then
    deallocate(umcCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: umcCPU", 1946,  istat,  errorMessage)
  endif

  if (allocated(vmrCPU)) then
    deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vmrCPU", 1951,  istat,  errorMessage)
  endif


  call obj%timer%stop("bandred_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

end subroutine bandred_&
&complex&
&_&
&double

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





subroutine herm_matrix_allreduce_&
&double &
                               (obj, n, a, lda, ldb, comm, isRows)
!-------------------------------------------------------------------------------
!  herm_matrix_allreduce: Does an mpi_allreduce for a hermitian matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
  use elpa_abstract_impl
  use precision
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)               :: n, lda, ldb, comm
  complex(kind=CK8) :: a(lda,ldb)

  integer(kind=ik)               :: i, nc
  integer(kind=MPI_KIND)         :: mpierr
  complex(kind=CK8) :: h1(n*n), h2(n*n)
  integer(kind=MPI_KIND)         :: allreduce_request1
  logical                        :: useNonBlockingCollectivesCols
  logical                        :: useNonBlockingCollectivesRows
  logical                        :: useNonBlockingCollectives
  logical, intent(in)            :: isRows
  integer(kind=c_int)            :: non_blocking_collectives_rows, error, &
                                    non_blocking_collectives_cols

  call obj%timer%start("herm_matrix_allreduce" // "_double")

  call obj%get("nbc_row_herm_allreduce", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa_herm_allreduce. Aborting..."
    stop
  endif
  call obj%get("nbc_col_herm_allreduce", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa_herm_allreduce. Aborting..."
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

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  nc = 0
  do i=1,n
    h1(nc+1:nc+i) = a(1:i,i)
    nc = nc+i
  enddo
  if (useNonBlockingCollectives) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                     int(comm,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                     int(comm,kind=MPI_KIND), mpierr)
    call obj%timer%stop("mpi_communication")
  endif

  nc = 0
  do i=1,n
    a(1:i,i) = h2(nc+1:nc+i)
    a(i,1:i-1) = conjg(a(1:i-1,i))
    nc = nc+i
  enddo


! nc = 0
! do i=1,n
!   a(1:i,i) = h2(nc+1:nc+i)
!   a(i,1:i-1) = conjg(a(1:i-1,i))
!   nc = nc+i
! enddo

  call obj%timer%stop("herm_matrix_allreduce" // "_double")

end subroutine herm_matrix_allreduce_&
&double







subroutine trans_ev_band_to_full_&
    &complex&
    &_&
    &double &
    (obj, na, nqc, nblk, nbw, a_mat, lda, tmat, q_mat, &
     ldq, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, useGPU &
     )

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_real/complex:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a_mat, number of rows of matrix q_mat
!
!  nqc         Number of columns of matrix q_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a_mat(lda,matrixCols)    Matrix containing the Householder vectors (i.e. matrix a_mat after bandred_real/complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat and q_mat
!
!  tmat(nbw,nbw,numBlocks) Factors returned by bandred_real/complex
!
!  q_mat           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q_mat
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa_abstract_impl
  use elpa_blas_interfaces

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
  class(elpa_abstract_impl_t), intent(inout) :: obj
  logical, intent(in)                    :: useGPU
  integer(kind=ik)                       :: na, nqc, lda, ldq, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols
  complex(kind=rck)                :: a_mat(lda,*)
  complex(kind=rck)                :: q_mat(ldq,*), tmat(nbw,nbw,*)

  integer(kind=ik)                       :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                 :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, mpierr
  integer(kind=ik)                       :: max_blocks_row, max_blocks_col, max_local_rows, &
                                            max_local_cols
  integer(kind=ik)                       :: l_cols, l_rows, l_colh, n_cols
  integer(kind=ik)                       :: istep, lc, ncol, nrow, nb, ns

  complex(kind=rck), allocatable   :: hvb(:)
  complex(kind=rck), pointer       :: hvm(:,:), tmp1(:), tmp2(:)
  ! hvm_dev is fist used and set in this routine
  ! q_mat is changed in trans_ev_tridi on the host, copied to device and passed here. this can be adapted
  ! tmp_dev is first used in this routine
  ! tmat_dev is not passed along from bandred_real
  integer(kind=C_intptr_T)               :: hvm_dev, q_dev, tmp_dev, tmat_dev, t_tmp_dev, dev_offset
  type(c_ptr)                            :: hvm_host, tmp1_host, tmp2_host
  type(c_ptr)                            :: tmat_mpi_dev, t_tmp_mpi_dev
  complex(kind=rck), pointer       :: tmat_mpi_deviceptr(:,:)
  complex(kind=rck), pointer       :: t_tmp_mpi_deviceptr(:)

  type(c_ptr)                            :: hvm_mpi_dev
  complex(kind=rck), pointer       :: hvm_mpi_deviceptr(:,:)


  integer(kind=ik)                       :: i

  complex(kind=rck), allocatable, target   :: tmat_complete(:,:), t_tmp(:,:), t_tmp2(:,:)
  integer(kind=ik)                       :: t_cols, t_rows, ii, jj
  integer(kind=ik)                       :: cwy_blocking

  integer(kind=ik)                       :: istat
  character(200)                         :: errorMessage
  character(20)                          :: gpuString
  logical                                :: successGPU
  integer(kind=c_intptr_t), parameter    :: size_of_datatype = size_of_&
                                                               &double&
                                                               &_&
                                                               &complex
  integer(kind=ik)                       :: blocking_factor, error, blk_end
  logical                                :: useIntelGPU
  integer(kind=MPI_KIND)                 :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                :: useNonBlockingCollectivesCols
  logical                                :: useNonBlockingCollectivesRows
  integer(kind=c_int)                    :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_band_to_full_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa2_band_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_band_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_band_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_band_to_full. Aborting..."
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

  call obj%get("blocking_in_band_to_full",blocking_factor,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem getting option for blocking_in_band_to_full. Aborting..."
    stop
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  max_blocks_row = ((na -1)/nblk)/np_rows + 1 ! Rows of a_mat
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1 ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  cwy_blocking = blocking_factor * nbw

  if (useGPU) then
    if (useIntelGPU) then
      allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp1", 249,  istat,  errorMessage)
  
      allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp2", 252,  istat,  errorMessage)
  
      allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: hvm", 255,  istat,  errorMessage)

    else
      ! copy q_mat to q_dev
      successGPU = gpu_malloc(q_dev,ldq*matrixCols*size_of_datatype)
      call check_alloc_GPU_f("trans_ev_band_to_full: q_dev", 260,  successGPU)

      successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuHostRegisterDefault)
      call check_host_register_GPU_f("trans_ev_band_to_full: q_mat", 264,  successGPU)

      successGPU = gpu_memcpy(q_dev,int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_mat -> q_dev", 268,  successGPU)

      successGPU = gpu_malloc_host(tmp1_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp1_host", 271,  successGPU)
      call c_f_pointer(tmp1_host, tmp1, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(tmp2_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp2_host", 275,  successGPU)
      call c_f_pointer(tmp2_host, tmp2, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(hvm_host,max_local_rows*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: hvm_host", 279,  successGPU)
      call c_f_pointer(hvm_host, hvm, (/max_local_rows,cwy_blocking/))
    endif
  else ! useGPU
    allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp1", 284,  istat,  errorMessage)

    allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp2", 287,  istat,  errorMessage)

    allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: hvm", 290,  istat,  errorMessage)
  endif !useGPU

  allocate(hvb(max_local_rows*cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: hvb", 294,  istat,  errorMessage)

  allocate(tmat_complete(cwy_blocking,cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: tmat_complete", 297,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_host_register(int(loc(tmat_complete),kind=c_intptr_t), &
                  cwy_blocking * cwy_blocking * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_band_to_full: tmat_complete", 303,  successGPU)
  endif
  !if (useIntelGPU) then
  !  ! needed later
  !endif



  if (blocking_factor > 1) then
    allocate(t_tmp(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp", 322,  istat,  errorMessage)

    allocate(t_tmp2(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp2", 325,  istat,  errorMessage)
  endif

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_malloc(hvm_dev,max_local_rows*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: hvm_dev", 330,  successGPU)

    successGPU = gpu_malloc(tmp_dev,max_local_cols*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmp_dev", 333,  successGPU)

    successGPU = gpu_malloc(tmat_dev,cwy_blocking*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmat_dev", 336,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  hvm = 0.0_rck ! Must be set to 0 !!!
  hvb = 0.0_rck ! Safety only
  tmp1 = 0.0_rck
  tmp2 = 0.0_rck
  tmat_complete = 0.0_rck
  if (blocking_factor > 1) then
     t_tmp = 0.0_rck ! Must be set to 0 !!!
     t_tmp2 = 0.0_rck

  endif
  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  blk_end = ((na-1)/nbw-1)/blocking_factor + 1
  do istep=1, blk_end

    ! This the call when using na >= ((blocking_factor+1)*nbw)
    ! n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw
    ! Number of columns in current step
    ! As an alternative we add some special case handling if na < cwy_blocking
    if (na < cwy_blocking) then
      n_cols = MAX(0, na-nbw)
      if ( n_cols .eq. 0 ) then
        exit
      end if
    else
      n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw ! Number of columns in current step
    end if

    ! Broadcast all Householder vectors for current step compressed in hvb

    nb = 0
    ns = 0

    do lc = 1, n_cols
      ncol = (istep-1)*cwy_blocking + nbw + lc ! absolute column number of householder Vector
      nrow = ncol - nbw ! absolute number of pivot row

      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
      l_colh = local_index(ncol , my_pcol, np_cols, nblk, -1) ! HV local column number

      if (my_pcol==pcol(ncol, nblk, np_cols)) hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)

      nb = nb+l_rows

      if (lc==n_cols .or. mod(ncol,nblk)==0) then
        if (useNonBlockingCollectivesCols) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_DOUBLE_COMPLEX,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         mpierr)
          call obj%timer%stop("mpi_communication")
        endif
        ns = nb
      endif
    enddo ! lc

    ! Expand compressed Householder vectors into matrix hvm

    nb = 0
    do lc = 1, n_cols
      nrow = (istep-1)*cwy_blocking + lc ! absolute number of pivot row
      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

      hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
      if (my_prow==prow(nrow, nblk, np_rows)) hvm(l_rows+1,lc) = 1.0_rck
      nb = nb+l_rows
    enddo

    l_rows = local_index(MIN(na,(istep+1)*cwy_blocking), my_prow, np_rows, nblk, -1)

    ! compute tmat2 out of tmat(:,:,)
    tmat_complete = 0
    do i = 1, blocking_factor
      t_cols = MIN(nbw, n_cols - (i-1)*nbw)
      if (t_cols <= 0) exit
      t_rows = (i - 1) * nbw
      tmat_complete(t_rows+1:t_rows+t_cols,t_rows+1:t_rows+t_cols) = tmat(1:t_cols,1:t_cols,(istep-1)*blocking_factor + i)

      if (i > 1) then
        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
          else 
          ! replace this with a cugemm if beneficial
          call obj%timer%start("blas")
          call ZGEMM('C', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
          endif
        else ! useGPU
          call obj%timer%start("blas")
          call ZGEMM('C', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU


        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif

        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
          else
            ! could use the cugemm here if beneficial
            call obj%timer%start("blas")
            call ZTRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call ZTRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                              tmat_complete(t_rows+1,t_rows+1), &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
            
          endif

        else ! useGPU
          call obj%timer%start("blas")
          call ZTRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call ZTRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                            tmat_complete(t_rows+1,t_rows+1), &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call obj%timer%stop("blas")
          tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
        endif ! useGPU


      endif
    enddo

    ! Q = Q - V * T**T * V**T * Q

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(hvm_dev, int(loc(hvm),kind=c_intptr_t), &
                          max_local_rows*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: hvm -> hvm_dev", 713,  successGPU)
          call obj%timer%start("gpublas")
          call gpublas_ZGEMM('C', 'N', &
                                       n_cols, l_cols, l_rows, ONE, hvm_dev, max_local_rows, &
                                       q_dev, ldq , ZERO, tmp_dev, n_cols)
          call obj%timer%stop("gpublas")

          ! copy data from device to host for a later MPI_ALLREDUCE
          successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                        tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 726,  successGPU)
        endif
      else ! useGPU
        call obj%timer%start("blas")
        call ZGEMM('C', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                            hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, &
                           int(n_cols,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    else ! l_rows>0
        tmp1(1:l_cols*n_cols) = 0.0_rck
    endif ! l_rows>0

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                   tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 764,  successGPU)

    endif

    if (useNonBlockingCollectivesRows) then
      call obj%timer%start("mpi_nbc_communication")
      call mpi_iallreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
      call obj%timer%stop("mpi_nbc_communication")
      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
    else
      call obj%timer%start("mpi_communication")
      call mpi_allreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
    endif

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          call ZTRMM('L', 'U', 'C', 'N', &
                              int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
          call ZGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                              int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                              q_mat, int(ldq,kind=BLAS_KIND))
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2),kind=c_intptr_t), &
                        l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp2 -> tmp_dev", 819,  successGPU)

          ! needed: as long as not device to device copy
          successGPU = gpu_memcpy(tmat_dev, int(loc(tmat_complete),kind=c_intptr_t), &
                        cwy_blocking*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmat_complete -> tmat_dev", 825,  successGPU)

          call obj%timer%start("gpublas")
          call gpublas_ZTRMM('L', 'U', 'C', 'N', &
                                   n_cols, l_cols, ONE, tmat_dev, cwy_blocking, tmp_dev, n_cols)
          call gpublas_ZGEMM('N', 'N', l_rows, l_cols, n_cols, -ONE, hvm_dev, max_local_rows, tmp_dev, &
                                     n_cols, ONE, q_dev, ldq)
          call obj%timer%stop("gpublas")
        endif
      else
        call obj%timer%start("blas")
        call ZTRMM('L', 'U', 'C', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
        call ZGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                            int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                            int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                            q_mat, int(ldq,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU

    endif

  enddo ! istep

  deallocate(hvb, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: hvb", 902,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      deallocate(tmp1, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp1", 907,  istat,  errorMessage)

      deallocate(tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp2", 910,  istat,  errorMessage)

      deallocate(hvm, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: hvm", 913,  istat,  errorMessage)

    else
      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: hvm_dev", 917,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmp_dev", 920,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmat_dev", 923,  successGPU)

      ! final transfer of q_dev
      successGPU = gpu_memcpy(int(loc(q_mat),kind=c_intptr_t), q_dev, ldq*matrixCols*size_of_datatype, &
                    gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_dev -> q_mat", 928,  successGPU)

      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: q_dev", 931,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: q_mat", 934,  successGPU)
      nullify(tmp1)
      nullify(tmp2)
      nullify(hvm)

      ! take care of new pointers nullify them
      nullify(tmat_mpi_deviceptr)
      nullify(hvm_mpi_deviceptr)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp1_host", 948,  successGPU)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp2_host", 951,  successGPU)

      successGPU = gpu_free_host(hvm_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: hvm_host", 954,  successGPU)

      successGPU = gpu_host_unregister(int(loc(tmat_complete),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: tmat_complete", 957,  successGPU)
    endif
  else ! useGPU
    deallocate(tmp1, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp1", 961,  istat,  errorMessage)

    deallocate(tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp2", 964,  istat,  errorMessage)

    deallocate(hvm, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: hvm", 967,  istat,  errorMessage)
  endif ! useGPU

  deallocate(tmat_complete, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: tmat_complete", 971,  istat,  errorMessage)


  if (blocking_factor > 1) then
    deallocate(t_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp", 984,  istat,  errorMessage)

    deallocate(t_tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp2", 987,  istat,  errorMessage)
  endif

  call obj%timer%stop("trans_ev_band_to_full_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

end subroutine trans_ev_band_to_full_&
&complex&
    &_&
    &double






subroutine tridiag_band_&
  &complex&
  &_&
  &double &
  (obj, na, nb, nblk, a_mat, lda, d, e, matrixCols, &
  hh_trans, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, useGPU, wantDebug, nrThreads, isSkewsymmetric)
  !-------------------------------------------------------------------------------
  ! tridiag_band_real/complex:
  ! Reduces a real symmetric band matrix to tridiagonal form
  !
  !  na          Order of matrix a
  !
  !  nb          Semi bandwith
  !
  !  nblk        blocksize of cyclic distribution, must be the same in both directions!
  !
  !  a_mat(lda,matrixCols)    Distributed system matrix reduced to banded form in the upper diagonal
  !
  !  lda         Leading dimension of a
  !  matrixCols  local columns of matrix a
  !
  ! hh_trans : housholder vectors
  !
  !  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  mpi_comm_rows
  !  mpi_comm_cols
  !              MPI-Communicators for rows/columns
  !  mpi_comm_all
  !              MPI-Communicator for the total processor set
  !-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use precision
  use, intrinsic :: iso_c_binding
  use redist
  use elpa_blas_interfaces
  use elpa_skewsymmetric_blas
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
  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU, wantDebug
  logical, intent(in)                          :: isSkewsymmetric
  integer(kind=ik), intent(in)                 :: na, nb, nblk, lda, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  complex(kind=rck), intent(in)         :: a_mat(lda,*)
  real(kind=rk), intent(out)        :: d(na), e(na) ! set only on PE 0
  complex(kind=rck), intent(out), allocatable   :: hh_trans(:,:)

  real(kind=rk)                     :: vnorm2
  complex(kind=rck)                     :: hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
  complex(kind=rck)                     :: hd(nb), hs(nb)

  integer(kind=ik)                             :: i, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
  integer(kind=ik)                             :: my_pe, n_pes
  integer(kind=ik)                             :: my_prow, np_rows, my_pcol, np_cols
  integer(kind=MPI_KIND)                       :: my_peMPI, n_pesMPI, mpierr
  integer(kind=MPI_KIND)                       :: my_prowMPI, np_rowsMPI, my_pcolMPI, np_colsMPI
  integer(kind=MPI_KIND)                       :: ireq_ab, ireq_hv
  integer(kind=ik)                             :: na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off
  integer(kind=ik), intent(in)                 :: nrThreads
  integer(kind=ik), allocatable                :: global_id(:,:), hh_cnt(:), hh_dst(:)
  integer(kind=MPI_KIND), allocatable          :: ireq_hhr(:), ireq_hhs(:)
  integer(kind=ik), allocatable                :: limits(:), snd_limits(:,:)
  integer(kind=ik), allocatable                :: block_limits(:)
  complex(kind=rck), allocatable         :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)
  integer                                      :: istat
  integer(kind=ik)                             :: nblockEnd
  character(200)                               :: errorMessage
  character(20)                                :: gpuString

  logical                                      :: useIntelGPU

   integer(kind=MPI_KIND)                      :: allreduce_request1, allreduce_request2
   logical                                     :: useNonBlockingCollectivesAll
   integer(kind=c_int)                         :: non_blocking_collectives, error

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_band_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_all_elpa2_band_to_tridi", non_blocking_collectives, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives in elpa2_band_to_tridi. Aborting..."
    stop
  endif

  if (non_blocking_collectives .eq. 1) then
    useNonBlockingCollectivesAll = .true.
  else
    useNonBlockingCollectivesAll = .false.
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND) ,my_peMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_all,kind=MPI_KIND) ,n_pesMPI ,mpierr)

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND),my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND),np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND),my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND),np_colsMPI ,mpierr)

  my_pe = int(my_peMPI,kind=MPI_KIND)
  n_pes = int(n_pesMPI,kind=MPI_KIND)
  my_prow = int(my_prowMPI,kind=MPI_KIND)
  np_rows = int(np_rowsMPI,kind=MPI_KIND)
  my_pcol = int(my_pcolMPI,kind=MPI_KIND)
  np_cols = int(np_colsMPI,kind=MPI_KIND)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Get global_id mapping 2D procssor coordinates to global id

  allocate(global_id(0:np_rows-1,0:np_cols-1), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: global_id", 202,  istat,  errorMessage)

  global_id(:,:) = 0
  global_id(my_prow, my_pcol) = my_pe


  if (wantDebug) call obj%timer%start("mpi_communication")
  if (useNonBlockingCollectivesAll) then
    call mpi_iallreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
  else
    call mpi_allreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), mpierr)
  endif
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Total number of blocks in the band:

  nblocks_total = (na-1)/nb + 1

  ! Set work distribution

  allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: block_limits", 246,  istat,  errorMessage)

  call divide_band(obj,nblocks_total, n_pes, block_limits)

  ! nblocks: the number of blocks for my task
  nblocks = block_limits(my_pe+1) - block_limits(my_pe)

  ! allocate the part of the band matrix which is needed by this PE
  ! The size is 1 block larger than needed to avoid extensive shifts
  allocate(ab(2*nb,(nblocks+1)*nb), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: ab", 256,  istat,  errorMessage)

  ab = 0.0_rck ! needed for lower half, the extra block should also be set to 0 for safety

  ! n_off: Offset of ab within band
  n_off = block_limits(my_pe)*nb

  ! Redistribute band in a to ab
  call redist_band_&
  &complex&
  &_&
  &double&
  &(obj,a_mat, lda, na, nblk, nb, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, ab)

  ! Calculate the workload for each sweep in the back transformation
  ! and the space requirements to hold the HH vectors

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: limits", 274,  istat,  errorMessage)

  call determine_workload(obj,na, nb, np_rows, limits)
  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  do n = 1, nblocks_total
    call determine_workload(obj, nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    ! add to number of householder vectors
    ! please note: for nx==1 the one and only HH Vector is 0 and is neither calculated nor send below!
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_hh_vecs = num_hh_vecs + local_size
      num_chunks  = num_chunks+1
    endif
    nx = nx - nb
  enddo

  ! Allocate space for HH vectors

  allocate(hh_trans(nb,num_hh_vecs), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_trans", 297,  istat,  errorMessage)

  ! Allocate and init MPI requests

  allocate(ireq_hhr(num_chunks), stat=istat, errmsg=errorMessage) ! Recv requests
  call check_allocate_f("tridiag_band: ireq_hhr", 302,  istat,  errorMessage)
  allocate(ireq_hhs(nblocks), stat=istat, errmsg=errorMessage)    ! Send requests
  call check_allocate_f("tridiag_band: ireq_hhs", 304,  istat,  errorMessage)

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  nt = 0
  nBlockEnd=1
  do n = 1, nblocks_total
    call determine_workload(obj,nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_chunks  = num_chunks+1
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_irecv(hh_trans(1,num_hh_vecs+1), int(nb*local_size,kind=MPI_KIND),  MPI_COMPLEX16,     &
                    int(nt,kind=MPI_KIND), int(10+n-block_limits(nt),kind=MPI_KIND), &
                    int(mpi_comm_all,kind=MPI_KIND), ireq_hhr(num_chunks), mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")

      num_hh_vecs = num_hh_vecs + local_size
    endif
    nx = nx - nb
    if (n == block_limits(nt+1)) then
      nt = nt + 1
    endif
  enddo
  ireq_hhs(:) = MPI_REQUEST_NULL
  ! Buffers for gathering/sending the HH vectors

  allocate(hh_gath(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! gathers HH vectors
  call check_allocate_f("tridiag_band: hh_gath", 341,  istat,  errorMessage)

  allocate(hh_send(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! send buffer for HH vectors
  call check_allocate_f("tridiag_band: hh_send", 344,  istat,  errorMessage)

  hh_gath(:,:,:) = 0.0_rck
  hh_send(:,:,:) = 0.0_rck

  ! Some counters

  allocate(hh_cnt(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_cnt", 352,  istat,  errorMessage)

  allocate(hh_dst(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_dst", 355,  istat,  errorMessage)

  hh_cnt(:) = 1 ! The first transfomation Vector is always 0 and not calculated at all
  hh_dst(:) = 0 ! PE number for receive
  ireq_ab = MPI_REQUEST_NULL
  ireq_hv = MPI_REQUEST_NULL
  ! Limits for sending

  allocate(snd_limits(0:np_rows,nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: snd_limits", 366,  istat,  errorMessage)

  do iblk=1,nblocks
    call determine_workload(obj, na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
  enddo


  ! ---------------------------------------------------------------------------
  ! Start of calculations

  na_s = block_limits(my_pe)*nb + 1

  if (my_pe>0 .and. na_s<=na) then
    ! send first column to previous PE
    ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
    ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
    if (wantDebug) call obj%timer%start("mpi_communication")
    call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_COMPLEX16, &
                   int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), ireq_ab, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif


   do istep=1,na-nblockEnd

     if (my_pe==0) then
       n = MIN(na-na_s,nb) ! number of rows to be reduced
       hv(:) = 0.0_rck
       hd(:) = 0.0_rck
       tau = 0.0_rck

       ! Transform first column of remaining matrix
       ! Opposed to the real case, the last step (istep=na-1) is needed here for making
       ! the last subdiagonal element a real number

         vnorm2 = sum(real(ab(3:n+1,na_s-n_off),kind=rk8)**2+dimag(ab(3:n+1,na_s-n_off))**2)
         if (n<2) vnorm2 = 0.0_rk ! Safety only

          call hh_transform_&
               &complex&
               &_&
               &double &
                          (obj, ab(2,na_s-n_off), vnorm2, hf, tau, wantDebug)

          hv(1) = 1.0_rck
          hv(2:n) = ab(3:n+1,na_s-n_off)*hf

       d(istep) = real(ab(1,na_s-n_off), kind=rk)
       e(istep) = real(ab(2,na_s-n_off), kind=rk)

       if (istep == na-1) then

         d(na) = real(ab(1,na_s+1-n_off),kind=rk)
         e(na) = 0.0_rck
       endif
     else
       if (na>na_s) then
         ! Receive Householder Vector from previous task, from PE owning subdiagonal


         if (wantDebug) call obj%timer%start("mpi_communication")

         call mpi_recv(hv, int(nb,kind=MPI_KIND), MPI_COMPLEX16, &
                       int(my_pe-1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")


         tau = hv(1)
         hv(1) = 1.0_rck
       endif
      endif

      na_s = na_s+1
      if (na_s-n_off > nb) then
        ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
        ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rck
        n_off = n_off + nb
      endif


      do iblk=1,nblocks
        ns = na_s + (iblk-1)*nb - n_off ! first column in block
        ne = ns+nb-1                    ! last column in block

        if (ns+n_off>na) exit

        ! Store Householder Vector for back transformation

        hh_cnt(iblk) = hh_cnt(iblk) + 1

        hh_gath(1   ,hh_cnt(iblk),iblk) = tau
        hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

        if (hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
          ! Wait for last transfer to finish
          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
          ! Copy vectors into send buffer
          hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
          ! Send to destination

          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_isend(hh_send(1,1,iblk), int(nb*hh_cnt(iblk),kind=MPI_KIND), MPI_COMPLEX16, &
                         global_id(hh_dst(iblk), mod(iblk+block_limits(my_pe)-1,np_cols)), &
                         int(10+iblk,kind=MPI_KIND), int(mpi_comm_all,kind=MPI_KIND), ireq_hhs(iblk), mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! Reset counter and increase destination row
          hh_cnt(iblk) = 0
          hh_dst(iblk) = hh_dst(iblk)+1
        endif

        ! The following code is structured in a way to keep waiting times for
        ! other PEs at a minimum, especially if there is only one block.
        ! For this reason, it requests the last column as late as possible
        ! and sends the Householder Vector and the first column as early
        ! as possible.
        nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
        nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                      ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

        ! Multiply diagonal block and subdiagonal block with Householder Vector

        if (iblk==nblocks .and. nc==nb) then

          ! We need the last column from the next PE.
          ! First do the matrix multiplications without last column ...

          ! Diagonal block, the contribution of the last element is added below!
          ab(1,ne) = 0.0_rck

          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")

            call ZHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            ! Subdiagonal block
            if (nr>0) call ZGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")

            call ZHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                             hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            ! Subdiagonal block
            if (nr>0) call ZGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif

          ! ... then request last column ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_recv(ab(1,ne), int(nb+1,kind=MPI_KIND), MPI_COMPLEX16,  &
                       int(my_pe+1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and complete the result
          hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
          hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

        else

          ! Normal matrix multiply
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            call ZHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                                hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)

            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            call ZHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            if (nr>0) call ZGEMV('N', int(nr,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), tau, ab(nb+1,ns), &
                                        int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Calculate first column of subdiagonal block and calculate new
        ! Householder transformation for this column
        hv_new(:) = 0.0_rck ! Needed, last rows must be 0 for nr < nb
        tau_new = 0.0_rck
        if (nr>0) then

          ! complete (old) Householder transformation for first column

          ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

          ! calculate new Householder transformation ...
          if (nr>1) then
            vnorm2 = sum(real(ab(nb+2:nb+nr,ns),kind=rk8)**2+dimag(ab(nb+2:nb+nr,ns))**2)

            call hh_transform_&
                &complex&
                &_&
                &double &
               (obj, ab(nb+1,ns), vnorm2, hf, tau_new, wantDebug)
            hv_new(1) = 1.0_rck
            hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
            ab(nb+2:,ns) = 0.0_rck
          endif ! nr > 1

          ! ... and send it away immediatly if this is the last block

          if (iblk==nblocks) then
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            hv_s(1) = tau_new
            hv_s(2:) = hv_new(2:)

            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_isend(hv_s, int(nb,kind=MPI_KIND), MPI_COMPLEX16, &
                           int(my_pe+1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                           ireq_hv, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

          endif

        endif

        ! Transform diagonal block
        x = dot_product(hv(1:nc),hd(1:nc))*conjg(tau)

          hd(1:nc) = hd(1:nc) - 0.5_rk*x*hv(1:nc)
        if (my_pe>0 .and. iblk==1) then

          ! The first column of the diagonal block has to be send to the previous PE
          ! Calculate first column only ...
          ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*conjg(hv(1)) - hv(1:nc)*conjg(hd(1))
          ! ... send it away ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_wait(ireq_ab, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ab_s(1:nb+1) = ab(1:nb+1,ns)

          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_COMPLEX16, &
                         int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                         ireq_ab, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and calculate remaining columns with rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (nc>1) call ZHER2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                         hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            if (nc>1) call ZHER2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                       hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            if (wantDebug) call obj%timer%stop("blas")
          endif

        else
          ! No need to  send, just a rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            call ZHER2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND, hv, 1_BLAS_KIND, &
                                ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")
            call ZHER2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND, hv, 1_BLAS_KIND, &
                                ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

        if (nr>0) then
          if (nr>1) then
            if (useIntelGPU) then
              if (wantDebug) call obj%timer%start("mkl_offload")
              call ZGEMV('C', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("mkl_offload")
            else
              if (wantDebug) call obj%timer%start("blas")
              call ZGEMV('C', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("blas")
            endif
            x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
            h(2:nb) = h(2:nb) - x*hv(2:nb)
            ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update
            do i=2,nb
              ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*conjg(h(i)) - hs(1:nr)*conjg(hv(i))
            enddo
          else
            ! No double Householder transformation for nr=1, just complete the row
            do i=2,nb
              ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*conjg(hv(i))
            enddo
          endif
        endif

        ! Use new HH Vector for the next block
        hv(:) = hv_new(:)
        tau = tau_new

      enddo


  enddo ! istep

  ! Finish the last outstanding requests


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
  call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

  call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
  call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_barrier(int(mpi_comm_all,kind=MPI_KIND),mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")
  deallocate(ab, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ab", 1363,  istat,  errorMessage)

  deallocate(ireq_hhr, ireq_hhs, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ireq_hhr", 1366,  istat,  errorMessage)

  deallocate(hh_cnt, hh_dst, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_dst", 1369,  istat,  errorMessage)

  deallocate(hh_gath, hh_send, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_gath", 1372,  istat,  errorMessage)

  deallocate(limits, snd_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: limits", 1375,  istat,  errorMessage)

  deallocate(block_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: block_limits", 1378,  istat,  errorMessage)

  deallocate(global_id, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: global_id", 1381,  istat,  errorMessage)

  call obj%timer%stop("tridiag_band_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

  ! intel compiler bug makes these ifdefs necessary
end subroutine tridiag_band_complex_&
&double

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






subroutine trans_ev_tridi_to_band_&
&complex&
&_&
&double &
(obj, na, nev, nblk, nbw, q, ldq, matrixCols,         &
 hh_trans, mpi_comm_rows, mpi_comm_cols, wantDebug, useGPU, max_threads, success, &
 kernel)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_real/complex:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!  matrixCols  local columns of matrix q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use pack_unpack_cpu
  use pack_unpack_gpu
  use compute_hh_trafo
  use elpa_gpu
  use precision
  use, intrinsic :: iso_c_binding
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
  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU

  integer(kind=ik), intent(in)                 :: kernel
  integer(kind=ik), intent(in)                 :: na, nev, nblk, nbw, ldq, matrixCols, mpi_comm_rows, mpi_comm_cols

  complex(kind=rck), target              :: q(ldq,*)

  complex(kind=rck), intent(in),target   :: hh_trans(:,:)
  integer(kind=c_intptr_t)                     :: hh_trans_dev
  type(c_ptr)                                  :: hh_trans_mpi_dev
  complex(kind=rck), pointer             :: hh_trans_mpi_fortran_ptr(:,:)


  integer(kind=ik)                             :: np_rows, my_prow, np_cols, my_pcol
  integer(kind=MPI_KIND)                       :: np_rowsMPI, my_prowMPI, np_colsMPI, my_pcolMPI
  integer(kind=ik)                             :: i, j, ip, sweep, nbuf, l_nev, a_dim2
  integer(kind=ik)                             :: current_n, current_local_n, current_n_start, current_n_end
  integer(kind=ik)                             :: next_n, next_local_n, next_n_start, next_n_end
  integer(kind=ik)                             :: bottom_msg_length, top_msg_length, next_top_msg_length
  integer(kind=ik)                             :: stripe_width, last_stripe_width, stripe_count
  integer(kind=ik)                             :: num_result_blocks, num_result_buffers, num_bufs_recvd
  integer(kind=ik)                             :: a_off, current_tv_off, max_blk_size
  integer(kind=ik)                             :: src, src_offset, dst, offset, nfact, num_blk
  integer(kind=MPI_KIND)                       :: mpierr

  logical                                      :: flag
  complex(kind=rck), pointer             :: aIntern(:,:,:)
  complex(kind=rck)                      :: a_var

  type(c_ptr)                                  :: aIntern_ptr

  complex(kind=rck), allocatable, target :: row(:)
  
  integer(kind=c_intptr_t)                     :: row_dev
  type(c_ptr)                                  :: row_mpi_dev
  complex(kind=rck), pointer             :: row_mpi_fortran_ptr(:)

  complex(kind=rck), pointer             :: row_group(:,:)

  complex(kind=rck), allocatable, target :: top_border_send_buffer(:,:)
  complex(kind=rck), allocatable, target :: top_border_recv_buffer(:,:)
  complex(kind=rck), allocatable, target :: bottom_border_send_buffer(:,:)
  complex(kind=rck), allocatable, target :: bottom_border_recv_buffer(:,:)

  integer(kind=c_intptr_t)                     :: top_border_recv_buffer_dev, top_border_send_buffer_dev
  type(c_ptr)                                  :: top_border_recv_buffer_mpi_dev, top_border_send_buffer_mpi_dev 
  complex(kind=rck), pointer             :: top_border_recv_buffer_mpi_fortran_ptr(:,:), &
                                                  top_border_send_buffer_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                     :: bottom_border_send_buffer_dev, bottom_border_recv_buffer_dev
  type(c_ptr)                                  :: bottom_border_send_buffer_mpi_dev, bottom_border_recv_buffer_mpi_dev
  complex(kind=rck), pointer             :: bottom_border_send_buffer_mpi_fortran_ptr(:,:), &
                                                  bottom_border_recv_buffer_mpi_fortran_ptr(:,:)
  type(c_ptr)                                  :: aIntern_mpi_dev
  complex(kind=rck), pointer             :: aIntern_mpi_fortran_ptr(:,:,:)

  integer(kind=c_intptr_t)                     :: aIntern_dev
  integer(kind=c_intptr_t)                     :: bcast_buffer_dev

  type(c_ptr)                                  :: bcast_buffer_mpi_dev
  complex(kind=rck), pointer             :: bcast_buffer_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: num
  integer(kind=c_intptr_t)                     :: dev_offset, dev_offset_1
  integer(kind=c_intptr_t)                     :: row_group_dev

  type(c_ptr)                                  :: row_group_mpi_dev
  complex(kind=rck), pointer             :: row_group_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: q_dev
  type(c_ptr)                                  :: q_mpi_dev
  complex(kind=rck), pointer             :: q_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: hh_tau_dev
  integer(kind=ik)                             :: row_group_size, unpack_idx

  type(c_ptr)                                  :: row_group_host, bcast_buffer_host

  integer(kind=ik)                             :: n_times
  integer(kind=ik)                             :: chunk, this_chunk

  complex(kind=rck), allocatable,target  :: result_buffer(:,:,:)
  integer(kind=c_intptr_t)                     :: result_buffer_dev

  type(c_ptr)                                  :: result_buffer_mpi_dev
  complex(kind=rck), pointer             :: result_buffer_mpi_fortran_ptr(:,:,:)

  complex(kind=rck), pointer             :: bcast_buffer(:,:)

  integer(kind=ik)                             :: n_off

  integer(kind=MPI_KIND), allocatable          :: result_send_request(:), result_recv_request(:)
  integer(kind=ik), allocatable                :: limits(:)
  integer(kind=MPI_KIND), allocatable          :: top_send_request(:), bottom_send_request(:)
  integer(kind=MPI_KIND), allocatable          :: top_recv_request(:), bottom_recv_request(:)

  ! MPI send/recv tags, arbitrary

  integer(kind=ik), parameter                  :: bottom_recv_tag = 111
  integer(kind=ik), parameter                  :: top_recv_tag    = 222
  integer(kind=ik), parameter                  :: result_recv_tag = 333

  integer(kind=ik), intent(in)                 :: max_threads



  ! Just for measuring the kernel performance
  real(kind=c_double)                          :: kernel_time, kernel_time_recv ! MPI_WTIME always needs double
  ! long integer
  integer(kind=lik)                            :: kernel_flops, kernel_flops_recv

  logical, intent(in)                          :: wantDebug
  logical                                      :: success
  integer(kind=ik)                             :: istat, print_flops
  character(200)                               :: errorMessage
  character(20)                                :: gpuString
  logical                                      :: successGPU
  integer(kind=ik)                             :: error
  integer(kind=c_intptr_t), parameter          :: size_of_datatype = size_of_&
                                                                 &double&
                                                                 &_&
                                                                 &complex
  logical, parameter                           :: allComputeOnGPU = .false.

  integer(kind=MPI_KIND)                     :: bcast_request1, allreduce_request1, allreduce_request2, &
                                                allreduce_request3, allreduce_request4
  logical                                    :: useNonBlockingCollectivesCols
  logical                                    :: useNonBlockingCollectivesRows
  integer(kind=c_int)                        :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_tridi_to_band_&
  &complex&
  &" // &
  &"_double" //&
  gpuString)

  call obj%get("nbc_row_elpa2_tridi_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_tridi_to_band. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_tridi_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_tridi_to_band. Aborting..."
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

  n_times = 0
  if (useGPU) then
    unpack_idx = 0
    row_group_size = 0
  endif

  success = .true.
  kernel_time = 0.0
  kernel_flops = 0

  if (wantDebug) call obj%timer%start("mpi_communication")
  call MPI_Comm_rank(int(mpi_comm_rows,kind=MPI_KIND) , my_prowMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_rows,kind=MPI_KIND) , np_rowsMPI , mpierr)
  call MPI_Comm_rank(int(mpi_comm_cols,kind=MPI_KIND) , my_pcolMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_cols,kind=MPI_KIND) , np_colsMPI , mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")

  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &complex&
                            &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &complex&
                            &: band backtransform works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  nfact = nbw / nblk


  ! local number of eigenvectors
  l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

  if (l_nev==0) then
    stripe_width = 0
    stripe_count = 0
    last_stripe_width = 0

  else ! l_nev


    ! Suggested stripe width is 48 since 48*64 real*8 numbers should fit into
    ! every primary cache
    ! Suggested stripe width is 48 - should this be reduced for the complex case ???

    if (useGPU) then
      stripe_width = 1024 ! Must be a multiple of 4
      stripe_count = (l_nev - 1) / stripe_width + 1

    else ! useGPU

      call obj%get("stripewidth_complex",stripe_width, error)

      !stripe_width = 48 ! Must be a multiple of 2

      stripe_count = (l_nev-1)/stripe_width + 1

      ! Adapt stripe width so that last one doesn't get too small

      stripe_width = (l_nev-1)/stripe_count + 1



     if (kernel .eq. ELPA_2STAGE_COMPLEX_AVX512_BLOCK1 .or. &
         kernel .eq. ELPA_2STAGE_COMPLEX_AVX512_BLOCK2 .or. &
         kernel .eq. ELPA_2STAGE_COMPLEX_SVE512_BLOCK1 .or. &
         kernel .eq. ELPA_2STAGE_COMPLEX_SVE512_BLOCK2  &
     ) then

       stripe_width = ((stripe_width+7)/8)*8 ! Must be a multiple of 4 because of AVX-512 memory alignment of 64 bytes
                                       ! (4 * sizeof(double complex) == 64)

     else

       stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 2 because of AVX/SSE memory alignment of 32 bytes
                                       ! (2 * sizeof(double complex) == 32)
     endif
   endif ! useGPU

   last_stripe_width = l_nev - (stripe_count-1)*stripe_width

  endif ! l_nev

  ! Determine the matrix distribution at the beginning

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: limits", 604,  istat,  errorMessage)
  call determine_workload(obj,na, nbw, np_rows, limits)

  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  a_dim2 = max_blk_size + nbw

  if (useGPU) then

    if (allComputeOnGPU) then
      if (wantDebug) call obj%timer%start("cuda_memcpy")

      successGPU = gpu_malloc(q_dev, ldq*matrixCols* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: q_dev", 617,  successGPU)
      successGPU =  gpu_memcpy(q_dev, int(loc(q(1,1)),kind=c_intptr_t),  &
                               ldq*matrixCols * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q -> q_dev", 621,  successGPU)

      ! associate with c_ptr
      q_mpi_dev = transfer(q_dev, q_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(q_mpi_dev, q_mpi_fortran_ptr, &
                       [ldq,matrixCols])
      if (wantDebug) call obj%timer%stop("cuda_memcpy")

      successGPU = gpu_malloc(hh_trans_dev, size(hh_trans,dim=1)*size(hh_trans,dim=2)* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 631,  successGPU)
      ! associate with c_ptr
      hh_trans_mpi_dev = transfer(hh_trans_dev, hh_trans_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(hh_trans_mpi_dev, hh_trans_mpi_fortran_ptr, &
                       [size(hh_trans,dim=1),size(hh_trans,dim=2)])
      successGPU =  gpu_memcpy(c_loc(hh_trans_mpi_fortran_ptr(1,1)),  &
                               c_loc(hh_trans(1,1)), &
                               size(hh_trans,dim=1)*size(hh_trans,dim=2) * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band: hh_trans -> hh_trans_dev", 641,  successGPU)

    endif ! allComputeOnGPU

    num =  (stripe_width*a_dim2*stripe_count)* size_of_datatype
    successGPU = gpu_malloc(aIntern_dev, stripe_width*a_dim2*stripe_count* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 647,  successGPU)

    successGPU = gpu_memset(aIntern_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 650,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      aIntern_mpi_dev = transfer(aIntern_dev, aIntern_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(aIntern_mpi_dev, aIntern_mpi_fortran_ptr, &
                       [stripe_width,a_dim2,stripe_count])
    endif ! allComputeOnGPU

    ! "row_group" and "row_group_dev" are needed for GPU optimizations
    successGPU = gpu_malloc_host(row_group_host,l_nev*nblk*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 667,  successGPU)
    call c_f_pointer(row_group_host, row_group, (/l_nev,nblk/))

    row_group(:, :) = 0.0_rck
    num =  (l_nev*nblk)* size_of_datatype
    successGPU = gpu_malloc(row_group_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 673,  successGPU)

    successGPU = gpu_memset(row_group_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_group_dev", 676,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      row_group_mpi_dev = transfer(row_group_dev, row_group_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(row_group_mpi_dev, row_group_mpi_fortran_ptr, &
                       [l_nev,nblk])
    endif ! allComputeOnGPU

  else ! GPUs are not used



    if (posix_memalign(aIntern_ptr, 64_c_intptr_t, stripe_width*a_dim2*stripe_count*  &
        C_SIZEOF(a_var)) /= 0) then
      print *,"trans_ev_tridi_to_band_real: error when allocating aIntern"//errorMessage
      stop 1
    endif

    call c_f_pointer(aIntern_ptr, aIntern,[stripe_width,a_dim2,stripe_count] )
    !allocate(aIntern(stripe_width,a_dim2,stripe_count), stat=istat, errmsg=errorMessage)

    aIntern(:,:,:) = 0.0_rck
  endif !useGPU

  allocate(row(l_nev), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: row", 723,  istat,  errorMessage)

  row(:) = 0.0_rck

  if (useGPU .and. allComputeOnGPU) then
    num =  (l_nev)* size_of_datatype
    successGPU = gpu_malloc(row_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_dev", 730,  successGPU)

    successGPU = gpu_memset(row_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_dev", 733,  successGPU)

    ! associate with c_ptr
    row_mpi_dev = transfer(row_dev, row_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(row_mpi_dev, row_mpi_fortran_ptr, &
                     [l_nev])
  endif

  ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
  ! and transpose the matrix using stripes of given stripe_width for cache blocking.

  ! The peculiar way it is done below is due to the fact that the last row should be
  ! ready first since it is the first one to start below



  if (wantDebug) call obj%timer%start("ip_loop")
  do ip = np_rows-1, 0, -1
    if (my_prow == ip) then
      ! Receive my rows which have not yet been received
      src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
      do i=limits(ip)+1,limits(ip+1)
        src = mod((i-1)/nblk, np_rows)

        if (src < my_prow) then

          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
            &complex&
            &_gpu_&
            &double &
                          ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                                      stripe_width, last_stripe_width, a_dim2, l_nev,&
                                      row_group_size, nblk, unpack_idx, &
                                       i - limits(ip), .false., wantDebug, allComputeOnGPU)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                         int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            call unpack_row_&
                &complex&
                &_cpu_&
                &double &
                (obj,aIntern, row,i-limits(ip), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        elseif (src == my_prow) then

          src_offset = src_offset+1

          if (useGPU) then

            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &complex&
                 &_gpu_&
                 &double &
             ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                          stripe_width, last_stripe_width, a_dim2, l_nev,&
                          row_group_size, nblk, unpack_idx, &
                          i - limits(ip), .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_ZCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_group_mpi_fortran_ptr(1,row_group_size)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            else ! allComputeOnGPU
              row_group(:, row_group_size) = q(src_offset, 1:l_nev)
            endif ! allComputeOnGPU
          else ! useGPU
            row(:) = q(src_offset, 1:l_nev)
          endif ! useGPU


          if (useGPU) then

          else
            call unpack_row_&
                 &complex&
                 &_cpu_&
                 &double &
                            (obj,aIntern, row,i-limits(ip),  stripe_count, stripe_width, last_stripe_width)
          endif


        endif ! src == my_prow
      enddo ! i=limits(ip)+1,limits(ip+1)


      ! Send all rows which have not yet been send


      if (allComputeOnGPU .and. useGPU) then
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1

              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_ZCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1011,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

              if (wantDebug) call obj%timer%start("cuda_mpi_communication")
              call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

            endif
          enddo
        enddo
      else !  allComputeOnGPU
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1
              row(:) = q(src_offset, 1:l_nev)

              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            endif
          enddo
        enddo
      endif ! allComputeOnGPU

    else if (my_prow < ip) then

      ! Send all rows going to PE ip
      if (allComputeOnGPU .and. useGPU) then
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1

            if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
            call gpublas_ZCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
            if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

            ! is there a way to avoid this device_synchronize ?
            if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
            successGPU = gpu_devicesynchronize()
            call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1061,  successGPU)
            if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

            if (wantDebug) call obj%timer%start("cuda_mpi_communication")
            call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                          int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

          endif
        enddo
      else ! allComputeOnGPU 
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1
            row(:) = q(src_offset, 1:l_nev)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                        int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          endif
        enddo
      endif  ! allComputeOnGPU

      ! Receive all rows from PE ip
      do i=limits(my_prow)+1,limits(my_prow+1)
        src = mod((i-1)/nblk, np_rows)
        if (src == ip) then
          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &complex&
                 &_gpu_&
                 &double&
                 &( obj, &
                 row_group, row_group_dev, aIntern_dev, stripe_count,  &
                 stripe_width, last_stripe_width, a_dim2, l_nev,       &
                 row_group_size, nblk, unpack_idx,                     &
                 i - limits(my_prow), .false., wantDebug, allComputeOnGPU)

            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX16, &
                             int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
            call unpack_row_&
                 &complex&
                 &_cpu_&
                 &double &
                 (obj,aIntern, row,i-limits(my_prow), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        endif
      enddo ! i=limits(my_prow)+1,limits(my_prow+1)
    endif ! (my_prow < ip)
  enddo ! ip = np_rows-1, 0, -1

  if (wantDebug) call obj%timer%stop("ip_loop")
 
  if (wantDebug) call obj%timer%start("allocate")

  if (useGPU) then
    ! Force an unpacking of all remaining rows that haven't been unpacked yet

    call unpack_and_prepare_row_group_&
         &complex&
         &_gpu_&
         &double&
         &( obj, &
         row_group, row_group_dev, aIntern_dev, stripe_count, &
         stripe_width, last_stripe_width, &
         a_dim2, l_nev, row_group_size, nblk, unpack_idx,     &
         -1, .true., wantDebug, allComputeOnGPU)

  endif

  ! Set up result buffer queue

  num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

  num_result_buffers = 4*nfact
  allocate(result_buffer(l_nev,nblk,num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_buffer", 1266,  istat,  errorMessage)

  allocate(result_send_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_send_request", 1269,  istat,  errorMessage)

  allocate(result_recv_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_recv_request", 1272,  istat,  errorMessage)

  if (useGPU .and. allComputeOnGPU) then
    num_result_buffers = 4*nfact
    num =  (l_nev*nblk*num_result_buffers)* size_of_datatype
    successGPU = gpu_malloc(result_buffer_dev, num* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 1278,  successGPU)

    ! associate with c_ptr
    result_buffer_mpi_dev = transfer(result_buffer_dev, result_buffer_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(result_buffer_mpi_dev, result_buffer_mpi_fortran_ptr, &
                     [l_nev,nblk,num_result_buffers])
  endif

  result_send_request(:) = MPI_REQUEST_NULL
  result_recv_request(:) = MPI_REQUEST_NULL

  ! Queue up buffers
  if (wantDebug) call obj%timer%start("host_mpi_communication")
  if (my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
    do j = 1, min(num_result_buffers, num_result_blocks)
      if (useGPU) then
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX16,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      else ! useGPU
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX16,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      endif ! useGPU
    enddo
  endif
  if (wantDebug) call obj%timer%stop("host_mpi_communication")

  num_bufs_recvd = 0 ! No buffers received yet

  ! Initialize top/bottom requests

  allocate(top_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_send_request", 1335,  istat,  errorMessage)

  allocate(top_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_recv_request", 1338,  istat,  errorMessage)

  allocate(bottom_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_send_request", 1341,  istat,  errorMessage)

  allocate(bottom_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_recv_request", 1344,  istat,  errorMessage)

  top_send_request(:) = MPI_REQUEST_NULL
  top_recv_request(:) = MPI_REQUEST_NULL
  bottom_send_request(:) = MPI_REQUEST_NULL
  bottom_recv_request(:) = MPI_REQUEST_NULL


  allocate(top_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 1374,  istat,  errorMessage)

  allocate(top_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1377,  istat,  errorMessage)

  allocate(bottom_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1380,  istat,  errorMessage)

  allocate(bottom_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1383,  istat,  errorMessage)

  top_border_send_buffer(:,:) = 0.0_rck
  top_border_recv_buffer(:,:) = 0.0_rck
  bottom_border_send_buffer(:,:) = 0.0_rck
  bottom_border_recv_buffer(:,:) = 0.0_rck


  if (useGPU) then
    if (allComputeOnGPU) then
      ! top_border_recv_buffer and top_border_send_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(top_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1401,  successGPU)

      successGPU = gpu_malloc(top_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1404,  successGPU)

      successGPU = gpu_memset(top_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1407,  successGPU)

      successGPU = gpu_memset(top_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1410,  successGPU)

      ! associate with c_ptr
      top_border_recv_buffer_mpi_dev = transfer(top_border_recv_buffer_dev, top_border_recv_buffer_mpi_dev)
      top_border_send_buffer_mpi_dev = transfer(top_border_send_buffer_dev, top_border_send_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(top_border_recv_buffer_mpi_dev, top_border_recv_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])
      call c_f_pointer(top_border_send_buffer_mpi_dev, top_border_send_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])

      ! bottom_border_send_buffer and bottom_border_recv_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(bottom_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1433,  successGPU)
      successGPU = gpu_malloc(bottom_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1435,  successGPU)

      successGPU = gpu_memset(bottom_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1438,  successGPU)
      successGPU = gpu_memset(bottom_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1440,  successGPU)

      ! associate with c_ptr
      bottom_border_send_buffer_mpi_dev = transfer(bottom_border_send_buffer_dev, bottom_border_send_buffer_mpi_dev)
      bottom_border_recv_buffer_mpi_dev = transfer(bottom_border_recv_buffer_dev, bottom_border_recv_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(bottom_border_send_buffer_mpi_dev, bottom_border_send_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
      call c_f_pointer(bottom_border_recv_buffer_mpi_dev, bottom_border_recv_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
    endif ! allComputeOnGPU

    successGPU = gpu_host_register(int(loc(top_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 1462,  successGPU)

    successGPU = gpu_host_register(int(loc(top_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1467,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1472,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1477,  successGPU)
  endif ! useGPU


  ! Initialize broadcast buffer

  if (useGPU) then
    successGPU = gpu_malloc_host(bcast_buffer_host,nbw*max_blk_size*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 1485,  successGPU)
    call c_f_pointer(bcast_buffer_host, bcast_buffer, (/nbw,max_blk_size/))
  else
    allocate(bcast_buffer(nbw, max_blk_size), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_tridi_to_band: bcast_buffer", 1489,  istat,  errorMessage)
  endif

  bcast_buffer = 0.0_rck

  if (useGPU) then
    num =  ( nbw * max_blk_size) * size_of_datatype
    successGPU = gpu_malloc(bcast_buffer_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1497,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      bcast_buffer_mpi_dev = transfer(bcast_buffer_dev, bcast_buffer_mpi_dev)
      call c_f_pointer(bcast_buffer_mpi_dev, bcast_buffer_mpi_fortran_ptr, &
                      [nbw, max_blk_size])
    endif ! allComputeOnGPU

    successGPU = gpu_memset( bcast_buffer_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1507,  successGPU)

    num =  (max_blk_size)* size_of_datatype
    successGPU = gpu_malloc( hh_tau_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1511,  successGPU)

    successGPU = gpu_memset( hh_tau_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1514,  successGPU)
  endif ! useGPU

  current_tv_off = 0 ! Offset of next row to be broadcast

  ! ------------------- start of work loop -------------------

  a_off = 0 ! offset in aIntern (to avoid unnecessary shifts)

  top_msg_length = 0
  bottom_msg_length = 0
  if (wantDebug) call obj%timer%stop("allocate")

  if (wantDebug) call obj%timer%start("sweep_loop")
  do sweep = 0, (na-1)/nbw

    current_n = na - sweep*nbw
    call determine_workload(obj,current_n, nbw, np_rows, limits)
    current_n_start = limits(my_prow)
    current_n_end   = limits(my_prow+1)
    current_local_n = current_n_end - current_n_start

    next_n = max(current_n - nbw, 0)
    call determine_workload(obj,next_n, nbw, np_rows, limits)
    next_n_start = limits(my_prow)
    next_n_end   = limits(my_prow+1)
    next_local_n = next_n_end - next_n_start

    if (next_n_end < next_n) then
      bottom_msg_length = current_n_end - next_n_end
    else
      bottom_msg_length = 0
    endif

    if (next_local_n > 0) then
      next_top_msg_length = current_n_start - next_n_start
    else
      next_top_msg_length = 0
    endif

    if (sweep==0 .and. current_n_end < current_n .and. l_nev > 0) then
      if (wantDebug) call obj%timer%start("host_mpi_communication")
      do i = 1, stripe_count


        if (useGPU) then
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        else !useGPU
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        endif !useGPU

      enddo
      if (wantDebug) call obj%timer%stop("host_mpi_communication")
    endif

    if (current_local_n > 1) then
      if (useGPU .and. allComputeOnGPU) then
        if (my_pcol == mod(sweep,np_cols)) then
          if (wantDebug) call obj%timer%start("cuda_memcpy")
          successGPU =  gpu_memcpy(c_loc(bcast_buffer_mpi_fortran_ptr(1,1)), &
                                   c_loc(hh_trans_mpi_fortran_ptr(1,current_tv_off+1)),  &
                                     size(hh_trans,dim=1) * (current_tv_off+current_local_n-(current_tv_off+1)+1) * &   
                                     size_of_datatype, &
                                     gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1646,  successGPU)
          if (wantDebug) call obj%timer%stop("cuda_memcpy")
          current_tv_off = current_tv_off + current_local_n
        endif
      else !useGPU
        if (my_pcol == mod(sweep,np_cols)) then
          bcast_buffer(:,1:current_local_n) =    &
          hh_trans(:,current_tv_off+1:current_tv_off+current_local_n)
          current_tv_off = current_tv_off + current_local_n
        endif
      endif ! useGPU

      if (useNonBlockingCollectivesCols) then
        if (wantDebug) call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_COMPLEX16, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_bcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_COMPLEX16, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif
 
      if (useGPU .and. .not.(allComputeOnGPU)) then
        if (wantDebug) call obj%timer%start("memcpy")
        successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
                                 nbw * current_local_n *    &
                                 size_of_datatype, &
                                 gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1702,  successGPU)
        if (wantDebug) call obj%timer%stop("memcpy")
      endif ! useGPU

      if (useGPU) then
        !if (.not.(allComputeOnGPU)) then
        !  if (wantDebug) call obj%timer%start("memcpy")
        !  successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
        !                           nbw * current_local_n *    &
        !                           size_of_datatype, &
        !                           gpuMemcpyHostToDevice)
        !  call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1725,  successGPU)
        !  if (wantDebug) call obj%timer%stop("memcpy")
        !endif ! allComputeOnGPU

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &complex&
             &_gpu_&
             &double&
             (bcast_buffer_dev, hh_tau_dev, nbw, &
             current_local_n, .false.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU

    else ! (current_local_n > 1) then

      ! for current_local_n == 1 the one and only HH Vector is 0 and not stored in hh_trans_real/complex
      bcast_buffer(:,1) = 0.0_rck
      if (useGPU) then
        successGPU = gpu_memset(bcast_buffer_dev, 0, nbw * size_of_datatype)
        call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1745,  successGPU)

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &complex&
             &_gpu_&
             &double&
             &( &
             bcast_buffer_dev, hh_tau_dev, &
             nbw, 1, .true.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU
    endif ! (current_local_n > 1) then

    if (l_nev == 0) cycle

    if (current_local_n > 0) then

      do i = 1, stripe_count

        !wait_b
        if (current_n_end < current_n) then



          if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_recv")
          call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)

          if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_recv")
          n_off = current_local_n+a_off

          if (useGPU) then
            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_memcpy")
              successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                       c_loc(bottom_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                       stripe_width*nbw* size_of_datatype,      &
                                       gpuMemcpyDeviceToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1859,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_memcpy")
            else ! allComputeOnGPU
              if (wantDebug) call obj%timer%start("memcpy")
              dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * size_of_datatype
              successGPU =  gpu_memcpy( aIntern_dev + dev_offset , &
                                      int(loc(bottom_border_recv_buffer(1,i)),kind=c_intptr_t), &
                                       stripe_width*nbw*  size_of_datatype,    &
                                       gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1868,  successGPU)
              if (wantDebug) call obj%timer%stop("memcpy")
            endif ! allComputeOnGPU
          else ! useGPU
            aIntern(:,n_off+1:n_off+nbw,i) = reshape( &
            bottom_border_recv_buffer(1:stripe_width*nbw,i),(/stripe_width,nbw/))
          endif ! useGPU


          if (next_n_end < next_n) then

            if (useGPU) then
              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            else ! useGPU
              if (wantDebug) call obj%timer%start("mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("mpi_communication")
            endif ! useGPU

            endif ! (next_n_end < next_n)
          endif ! (current_n_end < current_n)

          if (current_local_n <= bottom_msg_length + top_msg_length) then

            !wait_t
            if (top_msg_length>0) then


              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! the MPI_IRECV will be done CUDA_AWARE we thus do not need a host to device copy
                  ! However, we have to copy from top_border_recv_buffer_mpi_fortran_ptr to aIntern_mpi_fortran_ptr
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  !Fortran pointer for indexing
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                          c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                          stripe_width*top_msg_length* size_of_datatype,      &
                                          gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2045,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  !             host_offset= (0 + (0 * stripe_width) + ( (i-1) * stripe_width * nbw ) ) * 8
                  successGPU =  gpu_memcpy( aIntern_dev+dev_offset , int(loc(top_border_recv_buffer(1,i)),kind=c_intptr_t),  &
                                             stripe_width*top_msg_length* size_of_datatype,      &
                                             gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2054,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
            endif ! top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
                &complex&
                &_&
                &double&
                &(obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
                a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
                hh_tau_dev, kernel_flops, kernel_time, n_times, 0, current_local_n, i, &
                last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b        1
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")
            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")

            if (bottom_msg_length>0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer

                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                           c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2248,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")



              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 2255,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy( int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2263,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                      MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                      int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                     MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                     int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU
            endif !(bottom_msg_length>0)

          else ! current_local_n <= bottom_msg_length + top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &complex&
             &_&
             &double&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, &
             current_local_n - bottom_msg_length, bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")

            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")
            if (bottom_msg_length > 0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                            c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                             stripe_width * bottom_msg_length * size_of_datatype,      &
                                             gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2516,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset,  &
                                           stripe_width*bottom_msg_length* size_of_datatype,  &
                                           gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2524,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
          
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_COMPLEX16, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU




            endif ! (bottom_msg_length > 0)
            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &complex&
             &_&
             &double&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, top_msg_length, &
             current_local_n-top_msg_length-bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")


            !wait_t
            if (top_msg_length>0) then

              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
              if (useGPU) then
                if (allComputeOnGPU) then
                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                           c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                           stripe_width* top_msg_length* size_of_datatype,      &
                                           gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2739,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  ! copy top_border_recv_buffer to aIntern_dev, maybe not necessary if CUDA_AWARE IRECV
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(aIntern_dev + dev_offset ,int(loc( top_border_recv_buffer(:,i)),kind=c_intptr_t),  &
                                        stripe_width * top_msg_length * size_of_datatype,   &
                                        gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2748,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
           endif

           !compute

           if (wantDebug) call obj%timer%start("compute_hh_trafo")
           call compute_hh_trafo_&
             &complex&
             &_&
             &double&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off, nbw, max_blk_size,  bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, 0, top_msg_length, i, &
             last_stripe_width, kernel)
           if (wantDebug) call obj%timer%stop("compute_hh_trafo")

         endif

         if (next_top_msg_length > 0) then
           !request top_border data

           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX16, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else !useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX16, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif !useGPU


         endif

         !send_t
         if (my_prow > 0) then

           if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
           call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")

           if (useGPU) then
             if (allComputeOnGPU) then
               if (wantDebug) call obj%timer%start("cuda_memcpy")
               successGPU =  gpu_memcpy(c_loc(top_border_send_buffer_mpi_fortran_ptr(1,i)),  &
                                        c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                        stripe_width* nbw* size_of_datatype,      &
                                        gpuMemcpyDeviceToDevice)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 3030,  successGPU)
               if (wantDebug) call obj%timer%stop("cuda_memcpy")
             else ! allComputeOnGPU
               if (wantDebug) call obj%timer%start("memcpy")
               dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
               successGPU =  gpu_memcpy(int(loc(top_border_send_buffer(:,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                         stripe_width*nbw * size_of_datatype, &
                                         gpuMemcpyDeviceToHost)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> top_border_send_buffer", 3038,  successGPU)
               if (wantDebug) call obj%timer%stop("memcpy")
             endif ! allComputeOnGPU
           else ! useGPU
             top_border_send_buffer(:,i) = reshape(aIntern(:,a_off+1:a_off+nbw,i),(/stripe_width*nbw/))
           endif ! useGPU
           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_COMPLEX16, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else ! useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_COMPLEX16, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif ! useGPU

         endif

         ! Care that there are not too many outstanding top_recv_request's
         if (stripe_count > 1) then
           if (i>1) then

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
             call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
           else

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
            call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

          endif
        endif
      enddo ! i = 1, stripe_count

      top_msg_length = next_top_msg_length

    else
      ! wait for last top_send_request

      do i = 1, stripe_count
        if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
        call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")
      enddo
    endif

    ! Care about the result

    if (my_prow == 0) then

      ! topmost process sends nbw rows to destination processes

      do j=0, nfact-1
        num_blk = sweep*nfact+j ! global number of destination block, 0 based
        if (num_blk*nblk >= na) exit

        nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block

        if (wantDebug) call obj%timer%start("host_mpi_wait_result_send")
        call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_result_send")
        dst = mod(num_blk, np_rows)

        if (dst == 0) then
          if (useGPU) then
            row_group_size = min(na - num_blk*nblk, nblk)

            call pack_row_group_&
                 &complex&
                 &_gpu_&
                 &double&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
                         row_group(:, :), j * nblk + a_off, row_group_size, &
                         result_buffer_dev, nblk, num_result_buffers, nbuf, .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              ! memcpy DeviceToDevice row_group_dev -> q_dev
              do i = 1, row_group_size  
               if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
               call gpublas_ZCOPY(l_nev, c_loc(row_group_mpi_fortran_ptr(1,i)), 1, &
                                           c_loc(q_mpi_fortran_ptr((num_blk / np_rows) * nblk + i,1)), ldq)
               if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
              enddo
            else ! allComputeOnGPU
              do i = 1, row_group_size
                q((num_blk / np_rows) * nblk + i, 1 : l_nev) = row_group(:, i)
              enddo
            endif ! allComputeOnGPU
          else ! useGPU

            do i = 1, min(na - num_blk*nblk, nblk)

              call pack_row_&
                   &complex&
                   &_cpu_&
                   &double&
                   &(obj,aIntern, row, j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
              q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
            enddo
          endif ! useGPU

        else ! (dst == 0)

          if (useGPU) then

            call pack_row_group_&
                 &complex&
                 &_gpu_&
                 &double&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, &
                   last_stripe_width, a_dim2, l_nev, &
                   result_buffer(:, :, nbuf), j * nblk + a_off, nblk, &
                   result_buffer_dev, nblk, num_result_buffers, nbuf, .true., wantDebug, allComputeOnGPU)

          else  ! useGPU
            do i = 1, nblk
              call pack_row_&
                   &complex&
                   &_cpu_&
                   &double&
                   &(obj, aIntern, result_buffer(:,i,nbuf),j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
            enddo
          endif ! useGPU
         
          if (useGPU) then
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX16, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX16, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif ! useGPU
        endif ! (dst == 0)
      enddo  !j=0, nfact-1

    else ! (my_prow == 0)

      ! receive and store final result

      do j = num_bufs_recvd, num_result_blocks-1

        nbuf = mod(j, num_result_buffers) + 1 ! buffer number to get this block

        ! If there is still work to do, just test for the next result request
        ! and leave the loop if it is not ready, otherwise wait for all
        ! outstanding requests

        if (next_local_n > 0) then
            ! needed
            !successGPU = gpu_devicesynchronize()
          if (wantDebug) call obj%timer%start("host_mpi_test")
          call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_test")

          if (.not.flag) exit

        else ! (next_local_n > 0)
          if (wantDebug) call obj%timer%start("host_mpi_wait_result_recv")

          call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_wait_result_recv")
        endif ! (next_local_n > 0)

        ! Fill result buffer into q
        num_blk = j*np_rows + my_prow ! global number of current block, 0 based
        if (useGPU) then
          if (allComputeOnGPU) then
            do i = 1, min(na - num_blk*nblk, nblk)
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_ZCOPY(l_nev, c_loc(result_buffer_mpi_fortran_ptr(1,i,nbuf)), 1, &
                                          c_loc(q_mpi_fortran_ptr(j*nblk + i,1)), ldq)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            enddo
          else ! allComputeOnGPU
            do i = 1, min(na - num_blk*nblk, nblk)
              q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
            enddo
          endif ! allComputeOnGPU
        else ! useGPU
          do i = 1, min(na - num_blk*nblk, nblk)
            q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
          enddo
        endif ! useGPU
        ! Queue result buffer again if there are outstanding blocks left

        if (useGPU) then
          if (wantDebug) call obj%timer%start("host_mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX16, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("host_mpi_communication")
        else ! useGPU
          if (wantDebug) call obj%timer%start("mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX16, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif !useGPU

      enddo ! j = num_bufs_recvd, num_result_blocks-1
      num_bufs_recvd = j

    endif ! (my_prow == 0)

    ! Shift the remaining rows to the front of aIntern (if necessary)

    offset = nbw - top_msg_length
    if (offset<0) then
      if (wantDebug) write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                                         &complex&
                                         &: internal error, offset for shifting = ',offset
      success = .false.
      return
    endif

    a_off = a_off + offset
    if (a_off + next_local_n + nbw >= a_dim2) then
      do i = 1, stripe_count
        if (useGPU) then
          chunk = min(next_local_n,a_off)

          if (chunk < 1) exit

          if (wantDebug) call obj%timer%start("normal_memcpy")
          do j = top_msg_length+1, top_msg_length+next_local_n, chunk
            this_chunk = min(j+chunk-1,top_msg_length+next_local_n)-j+1
            dev_offset = ((j-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            dev_offset_1 = ((j+a_off-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            num = stripe_width*this_chunk*size_of_datatype
            successGPU = gpu_memcpy(aIntern_dev+dev_offset, aIntern_dev+dev_offset_1, num, gpuMemcpyDeviceToDevice)

            call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> aIntern_dev", 3528,  successGPU)
          end do
          if (wantDebug) call obj%timer%stop("normal_memcpy")
        else ! not useGPU
          do j = top_msg_length+1, top_msg_length+next_local_n
            aIntern(:,j,i) = aIntern(:,j+a_off,i)
          end do
        end if
      end do ! stripe_count

      a_off = 0
    end if
  end do ! sweep
  if (wantDebug) call obj%timer%stop("sweep_loop")

  ! Just for safety:
  if (ANY(top_send_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_send_request ***',my_prow,my_pcol
  if (ANY(bottom_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_send_request ***',my_prow,my_pcol
  if (ANY(top_recv_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_recv_request ***',my_prow,my_pcol
  if (ANY(bottom_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_recv_request ***',my_prow,my_pcol

  if (my_prow == 0) then

    if (wantDebug) call obj%timer%start("host_mpi_waitall")
    call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("host_mpi_waitall")
  endif

  if (ANY(result_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_send_request ***',my_prow,my_pcol
  if (ANY(result_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_recv_request ***',my_prow,my_pcol

  call obj%get("print_flops",print_flops,error)

  if (print_flops == 1) then

    if (useNonBlockingCollectivesRows) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          allreduce_request1, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          allreduce_request3, mpierr)

      call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif

    if (useNonBlockingCollectivesCols) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                          allreduce_request2, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          allreduce_request4, mpierr)

      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                         mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif
  endif

  if (useGPU .and. allComputeOnGPU) then
    ! finally copy q_dev to q
    if (wantDebug) call obj%timer%start("cuda_memcpy")
    successGPU =  gpu_memcpy(int(loc(q(1,1)),kind=c_intptr_t),  &
                             q_dev, &
                             ldq*matrixCols * size_of_datatype, &
                             gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q_dev -> q", 3634,  successGPU)
    if (wantDebug) call obj%timer%stop("cuda_memcpy")

  endif

  if (my_prow==0 .and. my_pcol==0 .and.print_flops == 1) &
      write(error_unit,'(" Kernel time:",f10.3," MFlops: ",es12.5)')  kernel_time, kernel_flops/kernel_time*1.d-6

  ! deallocate all working space

  if (.not.(useGPU)) then
    nullify(aIntern)
    call free(aIntern_ptr)
  endif

  deallocate(row, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: row", 3650,  istat,  errorMessage)

  deallocate(limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: limits", 3653,  istat,  errorMessage)

  deallocate(result_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_send_request", 3656,  istat,  errorMessage)

  deallocate(result_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_recv_request", 3659,  istat,  errorMessage)

  deallocate(result_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_buffer", 3662,  istat,  errorMessage)

  if (useGPU) then
    if (allComputeOnGPU) then
      successGPU = gpu_free(result_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 3667,  successGPU)
      nullify(result_buffer_mpi_fortran_ptr)
    endif

    nullify(bcast_buffer)

    successGPU = gpu_free_host(bcast_buffer_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 3674,  successGPU)
  else
    deallocate(bcast_buffer, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_tridi_to_band: bcast_buffer", 3677,  istat,  errorMessage)
  endif


  if (useGPU) then
    successGPU = gpu_free(aIntern_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 3683,  successGPU)

    if (allComputeOnGPU) then
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: q_dev", 3687,  successGPU)
      nullify(q_mpi_fortran_ptr)

      successGPU = gpu_free(hh_trans_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 3691,  successGPU)
      nullify(hh_trans_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 3695,  successGPU)
      nullify(top_border_recv_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 3699,  successGPU)
      nullify(top_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 3703,  successGPU)
      nullify(bottom_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 3707,  successGPU)
      nullify(bottom_border_recv_buffer_mpi_fortran_ptr)

      nullify(aIntern_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU = gpu_free(hh_tau_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 3714,  successGPU)

    nullify(row_group)

    successGPU = gpu_free_host(row_group_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 3719,  successGPU)

    successGPU = gpu_free(row_group_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 3722,  successGPU)

    if (allComputeOnGPU) then
      nullify(row_group_mpi_fortran_ptr)

      successGPU = gpu_free(row_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_dev", 3728,  successGPU)
      nullify(row_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU =  gpu_free(bcast_buffer_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 3733,  successGPU)

    if (allComputeOnGPU) then
      nullify(bcast_buffer_mpi_fortran_ptr)
    endif

    successGPU = gpu_host_unregister(int(loc(top_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 3740,  successGPU)

    successGPU = gpu_host_unregister(int(loc(top_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3743,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3746,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3749,  successGPU)
  endif ! useGPU

  deallocate(top_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 3753,  istat,  errorMessage)

  deallocate(top_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3756,  istat,  errorMessage)

  deallocate(bottom_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3759,  istat,  errorMessage)

  deallocate(bottom_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3762,  istat,  errorMessage)

  deallocate(top_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_send_request", 3765,  istat,  errorMessage)

  deallocate(top_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_recv_request", 3768,  istat,  errorMessage)

  deallocate(bottom_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_send_request", 3771,  istat,  errorMessage)

  deallocate(bottom_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_recv_request", 3774,  istat,  errorMessage)

  call obj%timer%stop("trans_ev_tridi_to_band_&
                      &complex&
                      &" // &
                      &"_double" //&
                      gpuString)
  return

end subroutine

! vim: syntax=fortran


! complex single precision























subroutine bandred_&
&complex&
&_&
&single &
(obj, na, a_mat, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, tmat, &
wantDebug, useGPU, success, &
max_threads, isSkewsymmetric)

!-------------------------------------------------------------------------------
!  bandred_real/complex: Reduces a distributed symmetric matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a_mat(lda,matrixCols)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a_mat(:,:) must be set completely (upper and lower half)
!              a_mat(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,numBlocks)    where numBlocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------

  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa1_compute
  use precision
  use elpa_blas_interfaces
  use elpa_scalapack_interfaces
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
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)                            :: na, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols

  complex(kind=rck)                     :: a_mat(lda,*)
  complex(kind=rck)                     :: tmat(nbw,nbw,*)

  logical, intent(in)                         :: useGPU
  logical, intent(in)                         :: isSkewsymmetric
  character(20)                               :: gpuString

  integer(kind=ik)                            :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                      :: mpierr,  my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI
  integer(kind=ik)                            :: l_cols, l_rows
  integer(kind=ik)                            :: i, j, lcs, lce, lre, lc, lr, cur_pcol, n_cols, nrow
  integer(kind=ik)                            :: istep, ncol, lch, lcx, nlc
  integer(kind=ik)                            :: tile_size, l_rows_tile, l_cols_tile

  real(kind=rk)                              :: vnorm2
  complex(kind=rck)                    :: xf, aux1(nbw), aux2(nbw), vrl, tau
  complex(kind=rck)                    :: vav(nbw,nbw)

  complex(kind=rck), allocatable :: tmpGPU(:)
  complex(kind=rck), pointer     :: vmrGPU(:), umcGPU(:)
  complex(kind=rck), allocatable :: tmpCPU(:,:), vmrCPU(:,:), umcCPU(:,:)
  complex(kind=rck), allocatable :: vr(:)

  integer(kind=C_intptr_T)                    :: a_dev, vmr_dev, umc_dev, tmat_dev, vav_dev
  type(c_ptr)                                 :: vmr_host, umc_host
  !integer(kind=ik), external                  :: numroc -> use elpa_scalapack
  integer(kind=ik)                            :: ierr
  integer(kind=ik)                            :: cur_l_rows, cur_l_cols, vmr_size, umc_size
  integer(kind=ik)                            :: l_rows2, vmr_size2, umc_size2
  integer(kind=c_intptr_t)                    :: lc_start, lc_end
  integer(kind=c_intptr_t)                    :: lce_1, lcs_1, lre_1
  integer(kind=ik)                            :: lr_end
  integer(kind=ik)                            :: na_cols
  integer(kind=BLAS_KIND)                     :: na_colsBLAS
  integer(kind=ik)                            :: na_rows
  integer(kind=BLAS_KIND)                     :: na_rowsBLAS

  logical, intent(in)                         :: wantDebug
  logical, intent(out)                        :: success
  logical                                     :: successGPU
  integer(kind=ik)                            :: istat
  character(200)                              :: errorMessage
  integer(kind=ik)                            :: min_tile_size, error

  integer(kind=ik)                            :: mystart, myend, m_way, n_way, work_per_thread, m_id, n_id, n_threads, &
                                                ii, pp
  integer(kind=c_intptr_t), parameter           :: size_of_datatype = size_of_&
                                                                    &single&
                                                                    &_&
                                                                    &complex

  logical                                     :: useGPU_reduction_lower_block_to_tridiagonal
  integer(kind=ik), intent(in)                :: max_threads
  logical                                     :: do_memcpy
  integer(kind=ik)                            :: i_blk,blk_off, blk_end
  logical                                     :: useIntelGPU


  integer(kind=MPI_KIND)                      :: bcast_request, allreduce_request1, allreduce_request2, &
                                                 allreduce_request3, allreduce_request4, allreduce_request5, &
                                                 allreduce_request6

  logical                                     :: useNonBlockingCollectivesCols
  logical                                     :: useNonBlockingCollectivesRows
  integer(kind=c_int)                         :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("bandred_&
  &complex&
  &" // &
  "_single" // &
  gpuString )

  call obj%get("nbc_row_elpa2_full_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_bandred. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_full_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_bandred. Aborting..."
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

  useGPU_reduction_lower_block_to_tridiagonal = .false.

  if (useGPU .and. .not.(useIntelGPU)) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
  endif

  if (useIntelGPU) then
    useGPU_reduction_lower_block_to_tridiagonal = .true.
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")
  success = .true.


  ! Semibandwith nbw must be a multiple of blocksize nblk
  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_bandred_&
                             &complex&
                             &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_bandred_&
                             &complex&
                             &: ELPA2 works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  ! na_rows in used nowhere; only na_cols
  if (useGPU .and. .not.(useIntelGPU)) then
    na_rowsBLAS = numroc(int(na,kind=BLAS_KIND), int(nblk,kind=BLAS_KIND), &
                         int(my_prow,kind=BLAS_KIND), 0_BLAS_KIND, int(np_rows,kind=BLAS_KIND))
    na_rows = int(na_rowsBLAS,kind=c_int)
    na_colsBLAS = numroc(int(na,kind=BLAS_KIND), int(nblk,kind=BLAS_KIND), &
                         int(my_pcol,kind=BLAS_KIND), 0_BLAS_KIND, int(np_cols,kind=BLAS_KIND))
    na_cols = int(na_colsBLAS,kind=c_int)

    ! Here we convert the regular host array into a pinned host array
    successGPU = gpu_malloc(a_dev, lda*na_cols* size_of_datatype)
    call check_alloc_GPU_f("bandred: a_dev", 335,  successGPU)

    successGPU = gpu_host_register(int(loc(vav),kind=c_intptr_t), &
                  nbw * nbw * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: vav", 340,  successGPU)

    successGPU = gpu_malloc(vav_dev, nbw*nbw* size_of_datatype)
    call check_alloc_GPU_f("bandred: vav_dev", 343,  successGPU)
  endif ! useGPU

  !if (useIntelGPU) then
  !  ! needed later when explicit copy
  !endif ! useIntelGPU

  ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

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

  l_rows_tile = tile_size/np_rows ! local rows of a tile
  l_cols_tile = tile_size/np_cols ! local cols of a tile


  blk_end = (na-1)/nbw
  if (useGPU .and. .not.(useIntelGPU)) then

    successGPU = gpu_host_register(int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuHostRegisterDefault)
    call check_host_register_GPU_f("bandred: a_mat", 422,  successGPU)

    cur_l_rows = 0
    cur_l_cols = 0

    successGPU = gpu_memcpy(a_dev, int(loc(a_mat),kind=c_intptr_t), &
                  lda*na_cols*size_of_datatype, gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("bandred: a_dev", 429,  successGPU)

    successGPU = gpu_malloc(tmat_dev, nbw*nbw*size_of_datatype)
    call check_alloc_GPU_f("bandred: tmat_dev", 432,  successGPU)

    istep = (na-1)/nbw
    blk_end = (na-1)/nbw
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows,1)
    cur_l_cols = max(l_cols,1)
    vmr_size = cur_l_rows*2*n_cols
    umc_size = cur_l_cols*2*n_cols

    istep = (na-1)/nbw - 1
    n_cols = min(na,(istep+1)*nbw)-istep*nbw
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows2 = local_index(istep*nbw, my_prow, np_rows, nblk, -1)
    cur_l_rows = max(l_rows2,1)
    cur_l_cols = max(l_cols,1)
    vmr_size2 = cur_l_rows*2*n_cols
    umc_size2 = cur_l_cols*2*n_cols

    l_rows = max(l_rows,l_rows2)
    vmr_size = max(vmr_size,vmr_size2)
    umc_size = max(umc_size,umc_size2)

    allocate(vr(l_rows + 1), stat=istat, errmsg=errorMessage)
    if (istat .ne. 0) then
      print *,"bandred_&
              &complex&
              &: error when allocating vr "//errorMessage
      stop 1
    endif

    successGPU = gpu_malloc_host(vmr_host,vmr_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: vmr_host", 466,  successGPU)
    call c_f_pointer(vmr_host, vmrGPU, (/vmr_size/))

    successGPU = gpu_malloc(vmr_dev, vmr_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: vmr_dev", 470,  successGPU)

    successGPU = gpu_malloc_host(umc_host,umc_size*size_of_datatype)
    call check_host_alloc_GPU_f("bandred: umc_host", 473,  successGPU)
    call c_f_pointer(umc_host, umcGPU, (/umc_size/))

    successGPU = gpu_malloc(umc_dev, umc_size*size_of_datatype)
    call check_alloc_GPU_f("bandred: umc_dev", 477,  successGPU)

  endif ! useGPU

  !if (useIntelGPU) then
     ! needed later when explict mem copy
  !endif ! useIntelGPU


  do istep = blk_end, 1, -1

    n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

    ! Number of local columns/rows of remaining matrix
    l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
    l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

    ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

    if (useGPU) then
      if (useIntelGPU) then
        ! unify the the name vmr and vmrCPU, as well as vmrGPU
        ! the same for umcCPU and umcGPU
        ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

        allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vmrCPU", 503,  istat,  errorMessage)

        allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: umcCPU", 506,  istat,  errorMessage)

        allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: vr", 509,  istat,  errorMessage)
      else
        cur_l_rows = max(l_rows, 1)
        cur_l_cols = max(l_cols, 1)
        vmr_size = cur_l_rows * 2 * n_cols
        umc_size = cur_l_cols * 2 * n_cols
      endif
    else ! GPU not used

      ! unify the the name vmr and vmrCPU, as well as vmrGPU
      ! the same for umcCPU and umcGPU
      ! Allocate vmr and umcCPU to their exact sizes so that they can be used in bcasts and reduces

      allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vmrCPU", 523,  istat,  errorMessage)

      allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: umcCPU", 526,  istat,  errorMessage)

      allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
      call check_allocate_f("bandred: vr", 529,  istat,  errorMessage)

    endif ! use GPU

    if (useGPU) then
      if (useIntelGPU) then
        vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
      else
        vmrGPU(1 : cur_l_rows * n_cols) = 0.0_rck
        umcGPU(1 : umc_size) = 0.0_rck
      endif
    else
      vmrCPU(1:l_rows,1:n_cols) = 0.0_rck
    endif ! useGPU

    vr(:) = 0.0_rck
    tmat(:,:,istep) = 0.0_rck
    if (useGPU .and. .not.(useIntelGPU)) then
      lc_start = local_index(istep*nbw+1, my_pcol, np_cols, nblk, -1)
      lc_end   = local_index(istep*nbw+n_cols, my_pcol, np_cols, nblk, -1)
      lr_end   = local_index((istep-1)*nbw + n_cols, my_prow, np_rows, nblk, -1)

      if (lc_start .le. 0) lc_start = 1

      do_memcpy = .false.

      ! Note: mod(nbw,nblk) == 0
      do i_blk = 1, nbw/nblk
        blk_off = (i_blk-1) * nblk
        cur_pcol = pcol(istep*nbw+1+blk_off, nblk, np_cols)

        if (my_pcol == cur_pcol) then
          do_memcpy = .true.
        endif
      enddo

      if (do_memcpy) then
        successGPU = gpu_memcpy2d(int(loc(a_mat(1, lc_start)),kind=c_intptr_t), &
                      int((lda*size_of_datatype),kind=c_intptr_t), &
                      (a_dev + int( ( (lc_start-1) * lda*size_of_datatype),kind=c_intptr_t )), &
                      int(lda*size_of_datatype,kind=c_intptr_t), &
                      int(lr_end*size_of_datatype,kind=c_intptr_t), &
                      int((lc_end - lc_start+1),kind=c_intptr_t),int(gpuMemcpyDeviceToHost,kind=c_int))

        call check_memcpy_GPU_f("bandred: a_dev -> a_mat", 573,  successGPU)
      endif
    endif ! useGPU

    !if (useIntelGPU) then
       ! needed later when explict mem copy
    !endif ! useIntelGPU


    ! Reduce current block to lower triangular form
      do lc = n_cols, 1, -1

        ncol = istep*nbw + lc ! absolute column number of householder Vector
        nrow = ncol - nbw ! Absolute number of pivot row

        lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
        lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

        tau = 0

        if (nrow == 1) exit ! Nothing to do

        cur_pcol = pcol(ncol, nblk, np_cols) ! Processor column owning current block

        if (my_pcol==cur_pcol) then

          ! Get Vector to be transformed; distribute last element and norm of
          ! remaining elements to all procs in current column

          vr(1:lr) = a_mat(1:lr,lch) ! Vector to be transformed

          if (my_prow==prow(nrow, nblk, np_rows)) then
            aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
            aux1(2) = vr(lr)
          else
            aux1(1) = dot_product(vr(1:lr),vr(1:lr))
            aux1(2) = 0.0_rck
          endif

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(aux1, aux2, 2_MPI_KIND, MPI_COMPLEX, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             allreduce_request1, mpierr)

            call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)

            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(aux1, aux2, 2_MPI_KIND, MPI_COMPLEX, &
                             MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                             mpierr)

            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif


          vnorm2 = real(aux2(1),kind=rk)
          vrl    = aux2(2)

          ! Householder transformation
          call hh_transform_&
             &complex&
             &_&
             &single &
                         (obj, vrl, vnorm2, xf, tau, wantDebug)
          ! Scale vr and store Householder Vector for back transformation

          vr(1:lr) = vr(1:lr) * xf
          if (my_prow==prow(nrow, nblk, np_rows)) then
            a_mat(1:lr-1,lch) = vr(1:lr-1)
            a_mat(lr,lch) = vrl
            vr(lr) = 1.0_rck
          else
            a_mat(1:lr,lch) = vr(1:lr)
          endif

        endif

        ! Broadcast Householder Vector and tau along columns

        vr(lr+1) = tau
        if (useNonBlockingCollectivesCols) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(vr, int(lr+1,kind=MPI_KIND), MPI_COMPLEX, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        bcast_request, mpierr)


          call mpi_wait(bcast_request, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_bcast(vr, int(lr+1,kind=MPI_KIND), MPI_COMPLEX, &
                        int(cur_pcol,kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                        mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
          vmrGPU(cur_l_rows * (lc - 1) + 1 : cur_l_rows * (lc - 1) + lr) = vr(1:lr)
        else
          vmrCPU(1:lr,lc) = vr(1:lr)
        endif
        tau = vr(lr+1)

        tmat(lc,lc,istep) = conjg(tau) ! Store tau in diagonal of tmat
        ! Transform remaining columns in current block with Householder Vector
        ! Local dot product

        aux1 = 0.0_rck


        nlc = 0 ! number of local columns
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            if (lr>0) aux1(nlc) = dot_product(vr(1:lr),a_mat(1:lr,lcx))
          endif
        enddo

        ! Get global dot products
        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          if (nlc>0) then
            call mpi_iallreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_COMPLEX, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      allreduce_request3, mpierr)
            call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          if (nlc>0) then
            call mpi_allreduce(aux1, aux2, int(nlc,kind=MPI_KIND), MPI_COMPLEX, &
                                      MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                                      mpierr)
          endif
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif
        ! Transform

        nlc = 0
        do j=1,lc-1
          lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
          if (lcx>0) then
            nlc = nlc+1
            a_mat(1:lr,lcx) = a_mat(1:lr,lcx) - conjg(tau)*aux2(nlc)*vr(1:lr)
          endif
        enddo
      enddo ! lc

      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        ! store column tiles back to GPU
        if (do_memcpy) then
          successGPU = gpu_memcpy2d((a_dev+ &
                        int(((lc_start-1)*lda*size_of_datatype),kind=c_intptr_t)), &
                        int(lda*size_of_datatype,kind=c_intptr_t), int(loc(a_mat(1,lc_start)),kind=c_intptr_t), &
                        int(lda*size_of_datatype,kind=c_intptr_t), &
                        int(lr_end*size_of_datatype,kind=c_intptr_t), &
                        int((lc_end - lc_start+1),kind=c_intptr_t), &
                        int(gpuMemcpyHostToDevice,kind=c_int))
          call check_memcpy_GPU_f("bandred: a_mat -> a_dev", 937,  successGPU)
        endif
      endif

      if (useGPU_reduction_lower_block_to_tridiagonal .and. useIntelGPU) then
        ! store column tiles back to GPU
      endif

      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use dsyrk

      vav = 0
      call obj%timer%start("blas")
      if (useGPU_reduction_lower_block_to_tridiagonal .and. .not.(useIntelGPU)) then
        if (l_rows>0) &
        call CHERK('U', 'C',            &
                           int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                           vmrGPU, int(cur_l_rows,kind=BLAS_KIND), &
                           ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      else ! useGPU_reduction_to_tridiagonal
        if (l_rows>0) &
        call CHERK('U', 'C',           &
                            int(n_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, vmrCPU, &
                            int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))
      endif
      call obj%timer%stop("blas")
      call herm_matrix_allreduce_&
         &single &
                         (obj, n_cols,vav, nbw, nbw,mpi_comm_rows, .true.)
         ! Calculate triangular matrix T for block Householder Transformation
      call obj%timer%start("blas")
      do lc=n_cols,1,-1
        tau = tmat(lc,lc,istep)
        if (lc<n_cols) then
          call CTRMV('U', 'C', 'N',&
                              int(n_cols-lc,kind=BLAS_KIND), tmat(lc+1,lc+1,istep), &
                              int(ubound(tmat,dim=1),kind=BLAS_KIND), vav(lc+1,lc), 1_BLAS_KIND)

          tmat(lc,lc+1:n_cols,istep) = -tau * conjg(vav(lc+1:n_cols,lc))
        endif
      enddo
      call obj%timer%stop("blas")


    ! Transpose vmr -> vmc (stored in umc, second half)
    if (useGPU) then
      if (useIntelGPU) then
        call elpa_transpose_vectors_&
             &complex&
             &_&
             &single &
                                          (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                           umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                           1, istep*nbw, n_cols, nblk, max_threads, .true.)

      else
        call elpa_transpose_vectors_&
             &complex&
             &_&
             &single &
                          (obj, vmrGPU(:), cur_l_rows, mpi_comm_rows, &
                           umcGPU(cur_l_cols * n_cols + 1:), cur_l_cols, &
                           mpi_comm_cols, 1, istep*nbw, n_cols, nblk, max_threads, .true.)
      endif
    else ! useGPU
      call elpa_transpose_vectors_&
           &complex&
           &_&
           &single &
                                        (obj, vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                         umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                         1, istep*nbw, n_cols, nblk, max_threads, .true.)
    endif

    ! Calculate umc = A**T * vmr
    ! Note that the distributed A has to be transposed
    ! Opposed to direct tridiagonalization there is no need to use the cache locality
    ! of the tiles, so we can use strips of the matrix


!#if 0
!    ! original complex implemetation check for performance
!    umcCPU(1:l_cols,1:n_cols) = 0.0_rck
!    vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
!
!    if (l_cols>0 .and. l_rows>0) then
!      do i=0,(istep*nbw-1)/tile_size
!
!        lcs = i*l_cols_tile+1
!        lce = min(l_cols,(i+1)*l_cols_tile)
!        if (lce<lcs) cycle
!
!        lre = min(l_rows,(i+1)*l_rows_tile)
!
!          call obj%timer%start("blas")
!          call CGEMM('C', 'N', lce-lcs+1, n_cols, lre, ONE, a_mat(1,lcs), ubound(a_mat,dim=1), &
!                     vmrCPU, ubound(vmrCPU,dim=1), ONE, umcCPU(lcs,1), ubound(umcCPU,dim=1))
!          call obj%timer%stop("blas")
!
!        if (i==0) cycle
!        lre = min(l_rows,i*l_rows_tile)
!          call obj%timer%start("blas")
!          call CGEMM('N', 'N', lre, n_cols, lce-lcs+1, ONE, a_mat(1,lcs), lda, &
!                     umcCPU(lcs,n_cols+1), ubound(umcCPU,dim=1), ONE, vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1))
!          call obj%timer%stop("blas")
!      enddo
!
!    endif ! (l_cols>0 .and. l_rows>0)
!#endif 

    !Code for Algorithm 4

    ! n_way is actually a branch for the number of OpenMP threads
    n_way = 1

      if (.not. useGPU .or. useIntelGPU) then
        umcCPU(1:l_cols,1:n_cols) = 0.0_rck
        vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rck
      endif ! useGPU

      if (l_cols>0 .and. l_rows>0) then

        if (useGPU .and. .not.(useIntelGPU)) then
          successGPU = gpu_memset(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                        0, cur_l_rows*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: vmr_dev", 1193,  successGPU)

          successGPU = gpu_memcpy(vmr_dev, int(loc(vmrGPU(1)),kind=c_intptr_t), &
                        cur_l_rows*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: vmrGPU -> vmr_dev", 1197,  successGPU)

          successGPU = gpu_memset(umc_dev, 0, l_cols*n_cols*size_of_datatype)
          call check_memset_GPU_f("bandred: umc_dev", 1200,  successGPU)

          successGPU = gpu_memcpy(umc_dev+l_cols*n_cols*size_of_datatype, &
                        int(loc(umcGPU(1+l_cols*n_cols)),kind=c_intptr_t), &
                        (umc_size-l_cols*n_cols)*size_of_datatype, &
                        gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev", 1206,  successGPU)
        endif ! useGPU

        do i=0,(istep*nbw-1)/tile_size

          lcs = i*l_cols_tile+1
          lce = min(l_cols,(i+1)*l_cols_tile)
          if (lce<lcs) cycle
          lre = min(l_rows,(i+1)*l_rows_tile)

          if (useGPU) then
            if (useIntelGPU) then
              call obj%timer%start("mkl_offload")
              call obj%timer%stop("mkl_offload")
              if (i==0) cycle
              lre = min(l_rows,i*l_rows_tile)
              call obj%timer%start("mkl_offload")

              if (isSkewsymmetric) then

              else
              endif
              call obj%timer%stop("mkl_offload")

            else
              call obj%timer%start("gpublas")
              call gpublas_CGEMM('C', 'N',                   &
                                         lce-lcs+1, n_cols, lre,     &
                                         ONE, (a_dev + ((lcs-1)*lda* &
                                         size_of_datatype)),         &
                                         lda, vmr_dev,cur_l_rows,    &
                                         ONE, (umc_dev+ (lcs-1)*     &
                                             size_of_datatype),      &
                                         cur_l_cols)

              call obj%timer%stop("gpublas")

              if(i==0) cycle
              call obj%timer%start("gpublas")

              lre = min(l_rows,i*l_rows_tile)
              if (isSkewsymmetric) then
                call gpublas_CGEMM('N', 'N', lre,n_cols, lce-lcs+1, -ONE, &
                              (a_dev+ ((lcs-1)*lda*                 &
                                    size_of_datatype)),             &
                         lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                size_of_datatype),              &
                                cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                              size_of_datatype),              &
                                cur_l_rows)
              else
                call gpublas_CGEMM('N', 'N', lre,n_cols, lce-lcs+1, ONE, &
                                            (a_dev+ ((lcs-1)*lda*                 &
                                                  size_of_datatype)),             &
                                       lda, (umc_dev+(cur_l_cols * n_cols+lcs-1)* &
                                              size_of_datatype),              &
                                              cur_l_cols, ONE, (vmr_dev+(cur_l_rows * n_cols)* &
                                            size_of_datatype),              &
                                              cur_l_rows)
              endif
              call obj%timer%stop("gpublas")
            endif
          else ! useGPU

            call obj%timer%start("blas")
            call CGEMM('C', 'N',       &
                                int(lce-lcs+1,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lre,kind=BLAS_KIND), &
                                ONE, a_mat(1,lcs), int(ubound(a_mat,dim=1),kind=BLAS_KIND), &
                                vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), ONE, umcCPU(lcs,1), &
                                int(ubound(umcCPU,dim=1),kind=BLAS_KIND) )
            call obj%timer%stop("blas")
            if (i==0) cycle
            lre = min(l_rows,i*l_rows_tile)
            call obj%timer%start("blas")

            if (isSkewsymmetric) then
              call CGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  -ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )

            else
              call CGEMM('N', 'N', int(lre,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(lce-lcs+1,kind=BLAS_KIND), &
                                  ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND),                                                   &
                                  umcCPU(lcs,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), ONE,                          &
                                  vmrCPU(1,n_cols+1), int(ubound(vmrCPU,dim=1), kind=BLAS_KIND) )
            endif
            call obj%timer%stop("blas")
          endif ! useGPU
        enddo ! i=0,(istep*nbw-1)/tile_size

        if (useGPU .and. .not.(useIntelGPU)) then
          if (tile_size < istep*nbw .or. n_way > 1) then
            successGPU = gpu_memcpy(int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                          vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                          (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyDeviceToHost)
            call check_memcpy_GPU_f("bandred: vmr_dev -> vmrGPU", 1342,  successGPU)
          endif

          successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                        umc_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU", 1347,  successGPU)
        endif ! useGPU
      endif ! l_cols>0 .and. l_rows>0

    ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
    ! on the processors containing the diagonal
    ! This is only necessary if ur has been calculated, i.e. if the
    ! global tile size is smaller than the global remaining matrix

    ! Or if we used the Algorithm 4
    if (tile_size < istep*nbw .or. n_way > 1) then

      if (useGPU) then
        if (useIntelGPU) then
          call elpa_reduce_add_vectors_&
          &complex&
          &_&
          &single &
                                      (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                       umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                      istep*nbw, n_cols, nblk, max_threads)

        else

          call elpa_reduce_add_vectors_&
               &complex&
               &_&
               &single &
                               (obj, vmrGPU(cur_l_rows * n_cols + 1:),cur_l_rows,  &
                                mpi_comm_rows, umcGPU,                            &
                                cur_l_cols, mpi_comm_cols, istep*nbw, n_cols, nblk, max_threads)
        endif
      else ! useGPU

        call elpa_reduce_add_vectors_&
        &complex&
        &_&
        &single &
                                         (obj, vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                          umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          istep*nbw, n_cols, nblk, max_threads)
      endif ! useGPU
    endif ! tile_size < istep*nbw .or. n_way > 1

    if (l_cols>0) then

      if (useGPU) then
        if (useIntelGPU) then
          allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpCPU", 1399,  istat,  errorMessage)

          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")
            call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           allreduce_request4, mpierr)
            call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), &
                           mpierr)
            umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
          call check_deallocate_f("bandred: tmpCPU", 1421,  istat,  errorMessage)
        else
          allocate(tmpGPU(l_cols * n_cols), stat=istat, errmsg=errorMessage)
          call check_allocate_f("bandred: tmpGPU", 1425,  istat,  errorMessage)
          if (useNonBlockingCollectivesRows) then
            if (wantDebug) call obj%timer%start("mpi_nbc_communication")

            call mpi_iallreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request5, mpierr)
            call mpi_wait(allreduce_request5, MPI_STATUS_IGNORE, mpierr)

            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
          else
            if (wantDebug) call obj%timer%start("mpi_communication")

            call mpi_allreduce(umcGPU, tmpGPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            umcGPU(1 : l_cols * n_cols) = tmpGPU(1 : l_cols * n_cols)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif

          if (allocated(tmpGPU)) then
            deallocate(tmpGPU, stat=istat, errmsg=errorMessage)
            call check_deallocate_f("bandred: tmpGPU", 1447,  istat,  errorMessage)
          endif
        endif

      else ! useGPU

        allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
        call check_allocate_f("bandred: tmpCPU", 1454,  istat,  errorMessage)

        if (useNonBlockingCollectivesRows) then
          if (wantDebug) call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request6, mpierr)
          call mpi_wait(allreduce_request6, MPI_STATUS_IGNORE, mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
        else
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_allreduce(umcCPU, tmpCPU, int(l_cols*n_cols,kind=MPI_KIND), MPI_COMPLEX,    &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif

        deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: tmpCPU", 1474,  istat,  errorMessage)
      endif ! useGPU
    endif ! l_cols > 0

    ! U = U * Tmat**T

    if (useGPU) then
      if (useIntelGPU) then
        call obj%timer%start("mkl_offload")

        call CTRMM('Right', 'Upper', 'C', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

        call CGEMM('C', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

        call CTRMM('Right', 'Upper', 'C', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
        call obj%timer%stop("mkl_offload")

      else
        successGPU = gpu_memcpy(umc_dev, int(loc(umcGPU(1)),kind=c_intptr_t), &
                      l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: umcGPU -> umc_dev ", 1526,  successGPU)

        successGPU = gpu_memcpy(tmat_dev,int(loc(tmat(1,1,istep)),kind=c_intptr_t), &
                      nbw*nbw*size_of_datatype,gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: tmat -> tmat_dev ", 1530,  successGPU)

        call obj%timer%start("gpublas")
        call gpublas_CTRMM('Right', 'Upper', 'C', 'Nonunit',  &
                            l_cols, n_cols, ONE, tmat_dev, nbw, umc_dev, cur_l_cols)
        call obj%timer%stop("gpublas")

        ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T
        call obj%timer%start("gpublas")
        call gpublas_CGEMM('C', 'N',             &
                                 n_cols, n_cols, l_cols, ONE, umc_dev, cur_l_cols, &
                                 (umc_dev+(cur_l_cols * n_cols )*size_of_datatype),cur_l_cols, &
                                 ZERO, vav_dev, nbw)

        call gpublas_CTRMM('Right', 'Upper', 'C', 'Nonunit',    &
           n_cols, n_cols, ONE, tmat_dev, nbw, vav_dev, nbw)
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(vav),kind=c_intptr_t), &
                    vav_dev, nbw*nbw*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: vav_dev -> vav ", 1550,  successGPU)
      endif
    else ! useGPU

      call obj%timer%start("blas")

      call CTRMM('Right', 'Upper', 'C', 'Nonunit',     &
                          int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep), &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), &
                          umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

      call CGEMM('C', 'N',              &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                          ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND), umcCPU(1,n_cols+1), &
                          int(ubound(umcCPU,dim=1),kind=BLAs_KIND), ZERO, vav, int(ubound(vav,dim=1),kind=BLAS_KIND))

      call CTRMM('Right', 'Upper', 'C', 'Nonunit',    &
                          int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), ONE, tmat(1,1,istep),    &
                          int(ubound(tmat,dim=1),kind=BLAS_KIND), vav, int(ubound(vav,dim=1),kind=BLAS_KIND) )
      call obj%timer%stop("blas")

    endif ! useGPU

    call herm_matrix_allreduce_&
         &single &
         (obj, n_cols,vav, nbw, nbw ,mpi_comm_cols, .false.)

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(vav_dev, int(loc(vav),kind=c_intptr_t), &
                       nbw*nbw*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("bandred: vav -> vav_dev ", 1599,  successGPU)
    endif
    !if (useIntelGPU) then
      ! needed later
    !endif


    ! U = U - 0.5 * V * VAV

    if (useGPU) then
      if (useIntelGPU) then
      call obj%timer%start("mkl_offload")
      call CGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                         (-0.5_rk, 0.0_rk),     &
                         umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav, &
                         int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      call obj%timer%stop("mkl_offload")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &complex&
        &_&
        &single &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &complex&
       &_&
       &single &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif

      else
        call obj%timer%start("gpublas")
        if (isSkewsymmetric) then
          call gpublas_CGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                    (0.5_rk, 0.0_rk), &
                                    (umc_dev+(cur_l_cols * n_cols )* &
                                    size_of_datatype),   &
                                    cur_l_cols, vav_dev,nbw,        &
                                    ONE, umc_dev, cur_l_cols)
        else
          call gpublas_CGEMM('N', 'N', l_cols, n_cols, n_cols,&
                                   (-0.5_rk, 0.0_rk), &
                                   (umc_dev+(cur_l_cols * n_cols )* &
                                   size_of_datatype),   &
                                   cur_l_cols, vav_dev,nbw,        &
                                   ONE, umc_dev, cur_l_cols)
        endif
        call obj%timer%stop("gpublas")

        successGPU = gpu_memcpy(int(loc(umcGPU(1)),kind=c_intptr_t), &
                    umc_dev, umc_size*size_of_datatype, gpuMemcpyDeviceToHost)
        call check_memcpy_GPU_f("bandred: umc_dev -> umcGPU ", 1681,  successGPU)

        ! Transpose umc -> umr (stored in vmr, second half)
        if (isSkewsymmetric) then
          call elpa_transpose_vectors_ss_&
             &complex&
             &_&
             &single &
                         (obj, umcGPU(:), cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        else
          call elpa_transpose_vectors_&
             &complex&
             &_&
             &single &
                         (obj, umcGPU, cur_l_cols, mpi_comm_cols, &
                          vmrGPU(cur_l_rows * n_cols + 1:), cur_l_rows, mpi_comm_rows, &
                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
        endif

        successGPU = gpu_memcpy(vmr_dev+cur_l_rows*n_cols*size_of_datatype, &
                    int(loc(vmrGPU(1+cur_l_rows*n_cols)),kind=c_intptr_t), &
                    (vmr_size-cur_l_rows*n_cols)*size_of_datatype, gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("bandred: vmr -> vmrGPU ", 1705,  successGPU)
      endif

    else ! useGPU
      call obj%timer%start("blas")
      call CGEMM('N', 'N', int(l_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND), int(n_cols,kind=BLAS_KIND),     &
                         (-0.5_rk, 0.0_rk),     &
                         umcCPU(1,n_cols+1), int(ubound(umcCPU,dim=1),kind=BLAS_KIND), vav, &
                         int(ubound(vav,dim=1),kind=BLAS_KIND), ONE, umcCPU, int(ubound(umcCPU,dim=1),kind=BLAS_KIND))

      call obj%timer%stop("blas")

      ! Transpose umc -> umr (stored in vmr, second half)
      if (isSkewsymmetric) then
        call elpa_transpose_vectors_ss_&
          &complex&
        &_&
        &single &
                                 (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                        vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                        1, istep*nbw, n_cols, nblk, max_threads, .false.)
      else
       call elpa_transpose_vectors_&
       &complex&
       &_&
       &single &
                                (obj, umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                          vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                          1, istep*nbw, n_cols, nblk, max_threads, .false.)
      endif
    endif  ! useGPU

    ! A = A - V*U**T - U*V**T


    do i=0,(istep*nbw-1)/tile_size
      lcs = i*l_cols_tile+1
      lce = min(l_cols,(i+1)*l_cols_tile)
      lre = min(l_rows,(i+1)*l_rows_tile)
      if (lce<lcs .or. lre<1) cycle

      if (useGPU) then
        if (useIntelGPU) then
          call obj%timer%start("mkl_offload")
          call CGEMM('N', 'C', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                              int(2*n_cols,kind=BLAS_KIND), &
                              -ONE, &
                              vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                              int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                              ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
          call obj%timer%stop("mkl_offload")

        else
          call obj%timer%start("gpublas")

          call gpublas_CGEMM('N', 'C',     &
                                     lre, lce-lcs+1, 2*n_cols, -ONE, &
                                     vmr_dev, cur_l_rows, (umc_dev +(lcs-1)*  &
                                     size_of_datatype), &
                                     cur_l_cols, ONE, (a_dev+(lcs-1)*lda* &
                                     size_of_datatype), lda)
          call obj%timer%stop("gpublas")
        endif
      else ! useGPU

        call obj%timer%start("blas")
        call CGEMM('N', 'C', int(lre,kind=BLAS_KIND),int(lce-lcs+1,kind=BLAS_KIND), &
                            int(2*n_cols,kind=BLAS_KIND), &
                            -ONE, &
                            vmrCPU, int(ubound(vmrCPU,dim=1),kind=BLAS_KIND), umcCPU(lcs,1), &
                            int(ubound(umcCPU,dim=1),kind=BLAS_KIND), &
                            ONE, a_mat(1,lcs), int(lda,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    enddo ! i=0,(istep*nbw-1)/tile_size

    if (.not.(useGPU) .or. useIntelGPU) then
      if (allocated(vr)) then
        deallocate(vr, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vr", 1871,  istat,  errorMessage)
      endif

      if (allocated(umcCPU)) then
        deallocate(umcCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: umcCPU", 1876,  istat,  errorMessage)
      endif

      if (allocated(vmrCPU)) then
        deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("bandred: vmrCPU", 1881,  istat,  errorMessage)
      endif
    endif !useGPU

  enddo ! istep - loop

  if (useGPU .and. .not.(useIntelGPU)) then
    ! copy a_dev to a_mat
    ! we do it here, since a is needed on the host in the following routine
    ! (band to tridi). Previously, a has been kept on the device and then
    ! copied in redist_band (called from tridiag_band). However, it seems to
    ! be easier to do it here.
    successGPU = gpu_memcpy(int(loc(a_mat),kind=c_intptr_t), &
                  int(a_dev,kind=c_intptr_t), &
                  int(lda*matrixCols* size_of_datatype, kind=c_intptr_t), &
                  gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("bandred: a_dev -> a_mat ", 1897,  successGPU)

    successGPU = gpu_host_unregister(int(loc(a_mat),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: a_mat ", 1900,  successGPU)

    successGPU = gpu_free(a_dev)
    call check_dealloc_GPU_f("bandred: a_dev ", 1903,  successGPU)

    successGPU = gpu_free(vav_dev)
    call check_dealloc_GPU_f("bandred: vav_dev ", 1906,  successGPU)

    successGPU = gpu_free(tmat_dev)
    call check_dealloc_GPU_f("bandred: tmat_dev ", 1909,  successGPU)

    successGPU = gpu_host_unregister(int(loc(vav),kind=c_intptr_t))
    call check_host_unregister_GPU_f("bandred: vav", 1912,  successGPU)

    if (associated(umcGPU)) then
      nullify(umcGPU)

      successGPU = gpu_free_host(umc_host)
      call check_host_dealloc_GPU_f("bandred: umc_host ", 1918,  successGPU)

      successGPU = gpu_free(umc_dev)
      call check_dealloc_GPU_f("bandred: umc_dev ", 1921,  successGPU)
    endif

    if (associated(vmrGPU)) then
      nullify(vmrGPU)

      successGPU = gpu_free_host(vmr_host)
      call check_host_dealloc_GPU_f("bandred: vmr_host ", 1928,  successGPU)

      successGPU = gpu_free(vmr_dev)
      call check_dealloc_GPU_f("bandred: vmr_dev ", 1931,  successGPU)
    endif
  endif ! useGPU

  !if (useIntelGPU) then
  !   ! needed later
  !endif

  if (allocated(vr)) then
    deallocate(vr, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vr", 1941,  istat,  errorMessage)
  endif

  if (allocated(umcCPU)) then
    deallocate(umcCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: umcCPU", 1946,  istat,  errorMessage)
  endif

  if (allocated(vmrCPU)) then
    deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("bandred: vmrCPU", 1951,  istat,  errorMessage)
  endif


  call obj%timer%stop("bandred_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

end subroutine bandred_&
&complex&
&_&
&single

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





subroutine herm_matrix_allreduce_&
&single &
                               (obj, n, a, lda, ldb, comm, isRows)
!-------------------------------------------------------------------------------
!  herm_matrix_allreduce: Does an mpi_allreduce for a hermitian matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
  use elpa_abstract_impl
  use precision
  implicit none
  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)               :: n, lda, ldb, comm
  complex(kind=CK4) :: a(lda,ldb)

  integer(kind=ik)               :: i, nc
  integer(kind=MPI_KIND)         :: mpierr
  complex(kind=CK4) :: h1(n*n), h2(n*n)
  integer(kind=MPI_KIND)         :: allreduce_request1
  logical                        :: useNonBlockingCollectivesCols
  logical                        :: useNonBlockingCollectivesRows
  logical                        :: useNonBlockingCollectives
  logical, intent(in)            :: isRows
  integer(kind=c_int)            :: non_blocking_collectives_rows, error, &
                                    non_blocking_collectives_cols

  call obj%timer%start("herm_matrix_allreduce" // "_single")

  call obj%get("nbc_row_herm_allreduce", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa_herm_allreduce. Aborting..."
    stop
  endif
  call obj%get("nbc_col_herm_allreduce", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa_herm_allreduce. Aborting..."
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

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  if (isRows) then
    useNonBlockingCollectives = useNonBlockingCollectivesRows
  else
    useNonBlockingCollectives = useNonBlockingCollectivesCols
  endif

  nc = 0
  do i=1,n
    h1(nc+1:nc+i) = a(1:i,i)
    nc = nc+i
  enddo
  if (useNonBlockingCollectives) then
    call obj%timer%start("mpi_nbc_communication")
    call mpi_iallreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                     int(comm,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
    call obj%timer%stop("mpi_nbc_communication")
  else
    call obj%timer%start("mpi_communication")
    call mpi_allreduce(h1, h2, int(nc,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                     int(comm,kind=MPI_KIND), mpierr)
    call obj%timer%stop("mpi_communication")
  endif

  nc = 0
  do i=1,n
    a(1:i,i) = h2(nc+1:nc+i)
    a(i,1:i-1) = conjg(a(1:i-1,i))
    nc = nc+i
  enddo


! nc = 0
! do i=1,n
!   a(1:i,i) = h2(nc+1:nc+i)
!   a(i,1:i-1) = conjg(a(1:i-1,i))
!   nc = nc+i
! enddo

  call obj%timer%stop("herm_matrix_allreduce" // "_single")

end subroutine herm_matrix_allreduce_&
&single







subroutine trans_ev_band_to_full_&
    &complex&
    &_&
    &single &
    (obj, na, nqc, nblk, nbw, a_mat, lda, tmat, q_mat, &
     ldq, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, useGPU &
     )

!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_real/complex:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a_mat, number of rows of matrix q_mat
!
!  nqc         Number of columns of matrix q_mat
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a_mat(lda,matrixCols)    Matrix containing the Householder vectors (i.e. matrix a_mat after bandred_real/complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a_mat
!  matrixCols  local columns of matrix a_mat and q_mat
!
!  tmat(nbw,nbw,numBlocks) Factors returned by bandred_real/complex
!
!  q_mat           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q_mat
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------
  use precision
  use elpa_gpu
  use, intrinsic :: iso_c_binding
  use elpa_abstract_impl
  use elpa_blas_interfaces

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
  class(elpa_abstract_impl_t), intent(inout) :: obj
  logical, intent(in)                    :: useGPU
  integer(kind=ik)                       :: na, nqc, lda, ldq, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols
  complex(kind=rck)                :: a_mat(lda,*)
  complex(kind=rck)                :: q_mat(ldq,*), tmat(nbw,nbw,*)

  integer(kind=ik)                       :: my_prow, my_pcol, np_rows, np_cols
  integer(kind=MPI_KIND)                 :: my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI, mpierr
  integer(kind=ik)                       :: max_blocks_row, max_blocks_col, max_local_rows, &
                                            max_local_cols
  integer(kind=ik)                       :: l_cols, l_rows, l_colh, n_cols
  integer(kind=ik)                       :: istep, lc, ncol, nrow, nb, ns

  complex(kind=rck), allocatable   :: hvb(:)
  complex(kind=rck), pointer       :: hvm(:,:), tmp1(:), tmp2(:)
  ! hvm_dev is fist used and set in this routine
  ! q_mat is changed in trans_ev_tridi on the host, copied to device and passed here. this can be adapted
  ! tmp_dev is first used in this routine
  ! tmat_dev is not passed along from bandred_real
  integer(kind=C_intptr_T)               :: hvm_dev, q_dev, tmp_dev, tmat_dev, t_tmp_dev, dev_offset
  type(c_ptr)                            :: hvm_host, tmp1_host, tmp2_host
  type(c_ptr)                            :: tmat_mpi_dev, t_tmp_mpi_dev
  complex(kind=rck), pointer       :: tmat_mpi_deviceptr(:,:)
  complex(kind=rck), pointer       :: t_tmp_mpi_deviceptr(:)

  type(c_ptr)                            :: hvm_mpi_dev
  complex(kind=rck), pointer       :: hvm_mpi_deviceptr(:,:)


  integer(kind=ik)                       :: i

  complex(kind=rck), allocatable, target   :: tmat_complete(:,:), t_tmp(:,:), t_tmp2(:,:)
  integer(kind=ik)                       :: t_cols, t_rows, ii, jj
  integer(kind=ik)                       :: cwy_blocking

  integer(kind=ik)                       :: istat
  character(200)                         :: errorMessage
  character(20)                          :: gpuString
  logical                                :: successGPU
  integer(kind=c_intptr_t), parameter    :: size_of_datatype = size_of_&
                                                               &single&
                                                               &_&
                                                               &complex
  integer(kind=ik)                       :: blocking_factor, error, blk_end
  logical                                :: useIntelGPU
  integer(kind=MPI_KIND)                 :: bcast_request1, allreduce_request1, allreduce_request2
  logical                                :: useNonBlockingCollectivesCols
  logical                                :: useNonBlockingCollectivesRows
  integer(kind=c_int)                    :: non_blocking_collectives_rows, non_blocking_collectives_cols


  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_band_to_full_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_row_elpa2_band_to_full", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_band_to_full. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_band_to_full", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_band_to_full. Aborting..."
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

  call obj%get("blocking_in_band_to_full",blocking_factor,error)
  if (error .ne. ELPA_OK) then
    print *,"Problem getting option for blocking_in_band_to_full. Aborting..."
    stop
  endif


  call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI ,mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)
  call obj%timer%stop("mpi_communication")

  max_blocks_row = ((na -1)/nblk)/np_rows + 1 ! Rows of a_mat
  max_blocks_col = ((nqc-1)/nblk)/np_cols + 1 ! Columns of q_mat!

  max_local_rows = max_blocks_row*nblk
  max_local_cols = max_blocks_col*nblk

  cwy_blocking = blocking_factor * nbw

  if (useGPU) then
    if (useIntelGPU) then
      allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp1", 249,  istat,  errorMessage)
  
      allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: tmp2", 252,  istat,  errorMessage)
  
      allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
      call check_allocate_f("trans_ev_band_to_full: hvm", 255,  istat,  errorMessage)

    else
      ! copy q_mat to q_dev
      successGPU = gpu_malloc(q_dev,ldq*matrixCols*size_of_datatype)
      call check_alloc_GPU_f("trans_ev_band_to_full: q_dev", 260,  successGPU)

      successGPU = gpu_host_register(int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuHostRegisterDefault)
      call check_host_register_GPU_f("trans_ev_band_to_full: q_mat", 264,  successGPU)

      successGPU = gpu_memcpy(q_dev,int(loc(q_mat),kind=c_intptr_t),&
                    ldq*matrixCols*size_of_datatype, gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_mat -> q_dev", 268,  successGPU)

      successGPU = gpu_malloc_host(tmp1_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp1_host", 271,  successGPU)
      call c_f_pointer(tmp1_host, tmp1, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(tmp2_host,max_local_cols*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: tmp2_host", 275,  successGPU)
      call c_f_pointer(tmp2_host, tmp2, (/max_local_cols*cwy_blocking/))

      successGPU = gpu_malloc_host(hvm_host,max_local_rows*cwy_blocking*size_of_datatype)
      call check_host_alloc_GPU_f("trans_ev_band_to_full: hvm_host", 279,  successGPU)
      call c_f_pointer(hvm_host, hvm, (/max_local_rows,cwy_blocking/))
    endif
  else ! useGPU
    allocate(tmp1(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp1", 284,  istat,  errorMessage)

    allocate(tmp2(max_local_cols*cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: tmp2", 287,  istat,  errorMessage)

    allocate(hvm(max_local_rows,cwy_blocking), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: hvm", 290,  istat,  errorMessage)
  endif !useGPU

  allocate(hvb(max_local_rows*cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: hvb", 294,  istat,  errorMessage)

  allocate(tmat_complete(cwy_blocking,cwy_blocking), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_band_to_full: tmat_complete", 297,  istat,  errorMessage)

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_host_register(int(loc(tmat_complete),kind=c_intptr_t), &
                  cwy_blocking * cwy_blocking * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_band_to_full: tmat_complete", 303,  successGPU)
  endif
  !if (useIntelGPU) then
  !  ! needed later
  !endif



  if (blocking_factor > 1) then
    allocate(t_tmp(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp", 322,  istat,  errorMessage)

    allocate(t_tmp2(cwy_blocking,nbw), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_band_to_full: t_tmp2", 325,  istat,  errorMessage)
  endif

  if (useGPU .and. .not.(useIntelGPU)) then
    successGPU = gpu_malloc(hvm_dev,max_local_rows*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: hvm_dev", 330,  successGPU)

    successGPU = gpu_malloc(tmp_dev,max_local_cols*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmp_dev", 333,  successGPU)

    successGPU = gpu_malloc(tmat_dev,cwy_blocking*cwy_blocking*size_of_datatype)
    call check_alloc_GPU_f("trans_ev_band_to_full: tmat_dev", 336,  successGPU)
  endif

  !if (useIntelGPU) then
  !  ! needed later
  !endif


  hvm = 0.0_rck ! Must be set to 0 !!!
  hvb = 0.0_rck ! Safety only
  tmp1 = 0.0_rck
  tmp2 = 0.0_rck
  tmat_complete = 0.0_rck
  if (blocking_factor > 1) then
     t_tmp = 0.0_rck ! Must be set to 0 !!!
     t_tmp2 = 0.0_rck

  endif
  l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q_mat

  blk_end = ((na-1)/nbw-1)/blocking_factor + 1
  do istep=1, blk_end

    ! This the call when using na >= ((blocking_factor+1)*nbw)
    ! n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw
    ! Number of columns in current step
    ! As an alternative we add some special case handling if na < cwy_blocking
    if (na < cwy_blocking) then
      n_cols = MAX(0, na-nbw)
      if ( n_cols .eq. 0 ) then
        exit
      end if
    else
      n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw ! Number of columns in current step
    end if

    ! Broadcast all Householder vectors for current step compressed in hvb

    nb = 0
    ns = 0

    do lc = 1, n_cols
      ncol = (istep-1)*cwy_blocking + nbw + lc ! absolute column number of householder Vector
      nrow = ncol - nbw ! absolute number of pivot row

      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
      l_colh = local_index(ncol , my_pcol, np_cols, nblk, -1) ! HV local column number

      if (my_pcol==pcol(ncol, nblk, np_cols)) hvb(nb+1:nb+l_rows) = a_mat(1:l_rows,l_colh)

      nb = nb+l_rows

      if (lc==n_cols .or. mod(ncol,nblk)==0) then
        if (useNonBlockingCollectivesCols) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_ibcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_COMPLEX,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         bcast_request1, mpierr)
          call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
        else
          call obj%timer%start("mpi_communication")
          call mpi_bcast(hvb(ns+1), int(nb-ns,kind=MPI_KIND), MPI_COMPLEX,&
                         int(pcol(ncol, nblk, np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), &
                         mpierr)
          call obj%timer%stop("mpi_communication")
        endif
        ns = nb
      endif
    enddo ! lc

    ! Expand compressed Householder vectors into matrix hvm

    nb = 0
    do lc = 1, n_cols
      nrow = (istep-1)*cwy_blocking + lc ! absolute number of pivot row
      l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

      hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
      if (my_prow==prow(nrow, nblk, np_rows)) hvm(l_rows+1,lc) = 1.0_rck
      nb = nb+l_rows
    enddo

    l_rows = local_index(MIN(na,(istep+1)*cwy_blocking), my_prow, np_rows, nblk, -1)

    ! compute tmat2 out of tmat(:,:,)
    tmat_complete = 0
    do i = 1, blocking_factor
      t_cols = MIN(nbw, n_cols - (i-1)*nbw)
      if (t_cols <= 0) exit
      t_rows = (i - 1) * nbw
      tmat_complete(t_rows+1:t_rows+t_cols,t_rows+1:t_rows+t_cols) = tmat(1:t_cols,1:t_cols,(istep-1)*blocking_factor + i)

      if (i > 1) then
        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
          else 
          ! replace this with a cugemm if beneficial
          call obj%timer%start("blas")
          call CGEMM('C', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
          endif
        else ! useGPU
          call obj%timer%start("blas")
          call CGEMM('C', 'N', &
                            int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, hvm, &
                            int(max_local_rows,kind=BLAS_KIND), hvm(:,(i-1)*nbw+1:), &
                            int(max_local_rows,kind=BLAS_KIND), ZERO, t_tmp, int(cwy_blocking, kind=BLAS_KIND))
          call obj%timer%stop("blas")
        endif ! useGPU


        if (useNonBlockingCollectivesRows) then
          call obj%timer%start("mpi_nbc_communication")
          call mpi_iallreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), allreduce_request1, mpierr)
          call obj%timer%stop("mpi_nbc_communication")
          call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
        else
          call obj%timer%start("mpi_communication")
          call mpi_allreduce(t_tmp, t_tmp2, int(cwy_blocking*nbw,kind=MPI_KIND), MPI_COMPLEX, &
                           MPI_SUM, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
          call obj%timer%stop("mpi_communication")
        endif

        if (useGPU) then
          if (useIntelGPU) then
            !call obj%timer%start("mkl_offload")
            !call obj%timer%stop("mkl_offload")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
          else
            ! could use the cugemm here if beneficial
            call obj%timer%start("blas")
            call CTRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call CTRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                              tmat_complete(t_rows+1,t_rows+1), &
                              int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
            call obj%timer%stop("blas")
            tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
            
          endif

        else ! useGPU
          call obj%timer%start("blas")
          call CTRMM('L', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call CTRMM('R', 'U', 'N', 'N', int(t_rows,kind=BLAS_KIND), int(t_cols,kind=BLAS_KIND), -ONE, &
                            tmat_complete(t_rows+1,t_rows+1), &
                            int(cwy_blocking,kind=BLAS_KIND), t_tmp2, int(cwy_blocking,kind=BLAS_KIND))
          call obj%timer%stop("blas")
          tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
        endif ! useGPU


      endif
    enddo

    ! Q = Q - V * T**T * V**T * Q

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(hvm_dev, int(loc(hvm),kind=c_intptr_t), &
                          max_local_rows*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: hvm -> hvm_dev", 713,  successGPU)
          call obj%timer%start("gpublas")
          call gpublas_CGEMM('C', 'N', &
                                       n_cols, l_cols, l_rows, ONE, hvm_dev, max_local_rows, &
                                       q_dev, ldq , ZERO, tmp_dev, n_cols)
          call obj%timer%stop("gpublas")

          ! copy data from device to host for a later MPI_ALLREDUCE
          successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                        tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 726,  successGPU)
        endif
      else ! useGPU
        call obj%timer%start("blas")
        call CGEMM('C', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), int(l_rows,kind=BLAS_KIND), ONE, &
                            hvm, int(ubound(hvm,dim=1),kind=BLAS_KIND), q_mat, int(ldq,kind=BLAS_KIND), ZERO, tmp1, &
                           int(n_cols,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU
    else ! l_rows>0
        tmp1(1:l_cols*n_cols) = 0.0_rck
    endif ! l_rows>0

    if (useGPU .and. .not.(useIntelGPU)) then
      successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                   tmp_dev, l_cols*n_cols*size_of_datatype, gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: tmp_dev -> tmp1", 764,  successGPU)

    endif

    if (useNonBlockingCollectivesRows) then
      call obj%timer%start("mpi_nbc_communication")
      call mpi_iallreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), allreduce_request2, mpierr)
      call obj%timer%stop("mpi_nbc_communication")
      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
    else
      call obj%timer%start("mpi_communication")
      call mpi_allreduce(tmp1, tmp2, int(n_cols*l_cols,kind=MPI_KIND), MPI_COMPLEX, MPI_SUM, &
                       int(mpi_comm_rows,kind=MPI_KIND), mpierr)
      call obj%timer%stop("mpi_communication")
    endif

    if (l_rows>0) then
      if (useGPU) then
        if (useIntelGPU) then
          !call obj%timer%start("mkl_offload")
          call CTRMM('L', 'U', 'C', 'N', &
                              int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                              int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
          call CGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                              int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                              int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                              q_mat, int(ldq,kind=BLAS_KIND))
          !call obj%timer%stop("mkl_offload")

        else
          successGPU = gpu_memcpy(tmp_dev, int(loc(tmp2),kind=c_intptr_t), &
                        l_cols*n_cols*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmp2 -> tmp_dev", 819,  successGPU)

          ! needed: as long as not device to device copy
          successGPU = gpu_memcpy(tmat_dev, int(loc(tmat_complete),kind=c_intptr_t), &
                        cwy_blocking*cwy_blocking*size_of_datatype, gpuMemcpyHostToDevice)
          call check_memcpy_GPU_f("trans_ev_band_to_full: tmat_complete -> tmat_dev", 825,  successGPU)

          call obj%timer%start("gpublas")
          call gpublas_CTRMM('L', 'U', 'C', 'N', &
                                   n_cols, l_cols, ONE, tmat_dev, cwy_blocking, tmp_dev, n_cols)
          call gpublas_CGEMM('N', 'N', l_rows, l_cols, n_cols, -ONE, hvm_dev, max_local_rows, tmp_dev, &
                                     n_cols, ONE, q_dev, ldq)
          call obj%timer%stop("gpublas")
        endif
      else
        call obj%timer%start("blas")
        call CTRMM('L', 'U', 'C', 'N', &
                            int(n_cols,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), ONE, tmat_complete, &
                            int(cwy_blocking,kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND))
        call CGEMM('N', 'N', int(l_rows,kind=BLAS_KIND), int(l_cols,kind=BLAS_KIND), &
                            int(n_cols,kind=BLAS_KIND), -ONE, hvm, &
                            int(ubound(hvm,dim=1),kind=BLAS_KIND), tmp2, int(n_cols,kind=BLAS_KIND), ONE, &
                            q_mat, int(ldq,kind=BLAS_KIND))
        call obj%timer%stop("blas")
      endif ! useGPU

    endif

  enddo ! istep

  deallocate(hvb, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: hvb", 902,  istat,  errorMessage)

  if (useGPU) then
    if (useIntelGPU) then
      deallocate(tmp1, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp1", 907,  istat,  errorMessage)

      deallocate(tmp2, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: tmp2", 910,  istat,  errorMessage)

      deallocate(hvm, stat=istat, errmsg=errorMessage)
      call check_deallocate_f("trans_ev_band_to_full: hvm", 913,  istat,  errorMessage)

    else
      successGPU = gpu_free(hvm_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: hvm_dev", 917,  successGPU)

      successGPU = gpu_free(tmp_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmp_dev", 920,  successGPU)

      successGPU = gpu_free(tmat_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: tmat_dev", 923,  successGPU)

      ! final transfer of q_dev
      successGPU = gpu_memcpy(int(loc(q_mat),kind=c_intptr_t), q_dev, ldq*matrixCols*size_of_datatype, &
                    gpuMemcpyDeviceToHost)
      call check_memcpy_GPU_f("trans_ev_band_to_full: q_dev -> q_mat", 928,  successGPU)

      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_band_to_full: q_dev", 931,  successGPU)

      successGPU = gpu_host_unregister(int(loc(q_mat),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: q_mat", 934,  successGPU)
      nullify(tmp1)
      nullify(tmp2)
      nullify(hvm)

      ! take care of new pointers nullify them
      nullify(tmat_mpi_deviceptr)
      nullify(hvm_mpi_deviceptr)

      successGPU = gpu_free_host(tmp1_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp1_host", 948,  successGPU)

      successGPU = gpu_free_host(tmp2_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: tmp2_host", 951,  successGPU)

      successGPU = gpu_free_host(hvm_host)
      call check_host_dealloc_GPU_f("trans_ev_band_to_full: hvm_host", 954,  successGPU)

      successGPU = gpu_host_unregister(int(loc(tmat_complete),kind=c_intptr_t))
      call check_host_unregister_GPU_f("trans_ev_band_to_full: tmat_complete", 957,  successGPU)
    endif
  else ! useGPU
    deallocate(tmp1, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp1", 961,  istat,  errorMessage)

    deallocate(tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: tmp2", 964,  istat,  errorMessage)

    deallocate(hvm, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: hvm", 967,  istat,  errorMessage)
  endif ! useGPU

  deallocate(tmat_complete, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_band_to_full: tmat_complete", 971,  istat,  errorMessage)


  if (blocking_factor > 1) then
    deallocate(t_tmp, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp", 984,  istat,  errorMessage)

    deallocate(t_tmp2, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_band_to_full: t_tmp2", 987,  istat,  errorMessage)
  endif

  call obj%timer%stop("trans_ev_band_to_full_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

end subroutine trans_ev_band_to_full_&
&complex&
    &_&
    &single






subroutine tridiag_band_&
  &complex&
  &_&
  &single &
  (obj, na, nb, nblk, a_mat, lda, d, e, matrixCols, &
  hh_trans, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, useGPU, wantDebug, nrThreads, isSkewsymmetric)
  !-------------------------------------------------------------------------------
  ! tridiag_band_real/complex:
  ! Reduces a real symmetric band matrix to tridiagonal form
  !
  !  na          Order of matrix a
  !
  !  nb          Semi bandwith
  !
  !  nblk        blocksize of cyclic distribution, must be the same in both directions!
  !
  !  a_mat(lda,matrixCols)    Distributed system matrix reduced to banded form in the upper diagonal
  !
  !  lda         Leading dimension of a
  !  matrixCols  local columns of matrix a
  !
  ! hh_trans : housholder vectors
  !
  !  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
  !
  !  mpi_comm_rows
  !  mpi_comm_cols
  !              MPI-Communicators for rows/columns
  !  mpi_comm_all
  !              MPI-Communicator for the total processor set
  !-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use precision
  use, intrinsic :: iso_c_binding
  use redist
  use elpa_blas_interfaces
  use elpa_skewsymmetric_blas
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
  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU, wantDebug
  logical, intent(in)                          :: isSkewsymmetric
  integer(kind=ik), intent(in)                 :: na, nb, nblk, lda, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  complex(kind=rck), intent(in)         :: a_mat(lda,*)
  real(kind=rk), intent(out)        :: d(na), e(na) ! set only on PE 0
  complex(kind=rck), intent(out), allocatable   :: hh_trans(:,:)

  real(kind=rk)                     :: vnorm2
  complex(kind=rck)                     :: hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
  complex(kind=rck)                     :: hd(nb), hs(nb)

  integer(kind=ik)                             :: i, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
  integer(kind=ik)                             :: my_pe, n_pes
  integer(kind=ik)                             :: my_prow, np_rows, my_pcol, np_cols
  integer(kind=MPI_KIND)                       :: my_peMPI, n_pesMPI, mpierr
  integer(kind=MPI_KIND)                       :: my_prowMPI, np_rowsMPI, my_pcolMPI, np_colsMPI
  integer(kind=MPI_KIND)                       :: ireq_ab, ireq_hv
  integer(kind=ik)                             :: na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off
  integer(kind=ik), intent(in)                 :: nrThreads
  integer(kind=ik), allocatable                :: global_id(:,:), hh_cnt(:), hh_dst(:)
  integer(kind=MPI_KIND), allocatable          :: ireq_hhr(:), ireq_hhs(:)
  integer(kind=ik), allocatable                :: limits(:), snd_limits(:,:)
  integer(kind=ik), allocatable                :: block_limits(:)
  complex(kind=rck), allocatable         :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)
  integer                                      :: istat
  integer(kind=ik)                             :: nblockEnd
  character(200)                               :: errorMessage
  character(20)                                :: gpuString

  logical                                      :: useIntelGPU

   integer(kind=MPI_KIND)                      :: allreduce_request1, allreduce_request2
   logical                                     :: useNonBlockingCollectivesAll
   integer(kind=c_int)                         :: non_blocking_collectives, error

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("tridiag_band_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

  useIntelGPU = .false.
  if (useGPU) then
    if (gpu_vendor() == INTEL_GPU) then
      useIntelGPU = .true.
    endif
  endif

  call obj%get("nbc_all_elpa2_band_to_tridi", non_blocking_collectives, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives in elpa2_band_to_tridi. Aborting..."
    stop
  endif

  if (non_blocking_collectives .eq. 1) then
    useNonBlockingCollectivesAll = .true.
  else
    useNonBlockingCollectivesAll = .false.
  endif

  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_comm_rank(int(mpi_comm_all,kind=MPI_KIND) ,my_peMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_all,kind=MPI_KIND) ,n_pesMPI ,mpierr)

  call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND),my_prowMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND),np_rowsMPI ,mpierr)
  call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND),my_pcolMPI ,mpierr)
  call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND),np_colsMPI ,mpierr)

  my_pe = int(my_peMPI,kind=MPI_KIND)
  n_pes = int(n_pesMPI,kind=MPI_KIND)
  my_prow = int(my_prowMPI,kind=MPI_KIND)
  np_rows = int(np_rowsMPI,kind=MPI_KIND)
  my_pcol = int(my_pcolMPI,kind=MPI_KIND)
  np_cols = int(np_colsMPI,kind=MPI_KIND)
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Get global_id mapping 2D procssor coordinates to global id

  allocate(global_id(0:np_rows-1,0:np_cols-1), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: global_id", 202,  istat,  errorMessage)

  global_id(:,:) = 0
  global_id(my_prow, my_pcol) = my_pe


  if (wantDebug) call obj%timer%start("mpi_communication")
  if (useNonBlockingCollectivesAll) then
    call mpi_iallreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), allreduce_request1, mpierr)
    call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
  else
    call mpi_allreduce(mpi_in_place, global_id, int(np_rows*np_cols,kind=MPI_KIND), mpi_integer, &
                         mpi_sum, int(mpi_comm_all,kind=MPI_KIND), mpierr)
  endif
  if (wantDebug) call obj%timer%stop("mpi_communication")

  ! Total number of blocks in the band:

  nblocks_total = (na-1)/nb + 1

  ! Set work distribution

  allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: block_limits", 246,  istat,  errorMessage)

  call divide_band(obj,nblocks_total, n_pes, block_limits)

  ! nblocks: the number of blocks for my task
  nblocks = block_limits(my_pe+1) - block_limits(my_pe)

  ! allocate the part of the band matrix which is needed by this PE
  ! The size is 1 block larger than needed to avoid extensive shifts
  allocate(ab(2*nb,(nblocks+1)*nb), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: ab", 256,  istat,  errorMessage)

  ab = 0.0_rck ! needed for lower half, the extra block should also be set to 0 for safety

  ! n_off: Offset of ab within band
  n_off = block_limits(my_pe)*nb

  ! Redistribute band in a to ab
  call redist_band_&
  &complex&
  &_&
  &single&
  &(obj,a_mat, lda, na, nblk, nb, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm_all, ab)

  ! Calculate the workload for each sweep in the back transformation
  ! and the space requirements to hold the HH vectors

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: limits", 274,  istat,  errorMessage)

  call determine_workload(obj,na, nb, np_rows, limits)
  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  do n = 1, nblocks_total
    call determine_workload(obj, nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    ! add to number of householder vectors
    ! please note: for nx==1 the one and only HH Vector is 0 and is neither calculated nor send below!
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_hh_vecs = num_hh_vecs + local_size
      num_chunks  = num_chunks+1
    endif
    nx = nx - nb
  enddo

  ! Allocate space for HH vectors

  allocate(hh_trans(nb,num_hh_vecs), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_trans", 297,  istat,  errorMessage)

  ! Allocate and init MPI requests

  allocate(ireq_hhr(num_chunks), stat=istat, errmsg=errorMessage) ! Recv requests
  call check_allocate_f("tridiag_band: ireq_hhr", 302,  istat,  errorMessage)
  allocate(ireq_hhs(nblocks), stat=istat, errmsg=errorMessage)    ! Send requests
  call check_allocate_f("tridiag_band: ireq_hhs", 304,  istat,  errorMessage)

  num_hh_vecs = 0
  num_chunks  = 0
  nx = na
  nt = 0
  nBlockEnd=1
  do n = 1, nblocks_total
    call determine_workload(obj,nx, nb, np_rows, limits)
    local_size = limits(my_prow+1) - limits(my_prow)
    if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
      num_chunks  = num_chunks+1
      if (wantDebug) call obj%timer%start("mpi_communication")
      call mpi_irecv(hh_trans(1,num_hh_vecs+1), int(nb*local_size,kind=MPI_KIND),  MPI_COMPLEX8,     &
                    int(nt,kind=MPI_KIND), int(10+n-block_limits(nt),kind=MPI_KIND), &
                    int(mpi_comm_all,kind=MPI_KIND), ireq_hhr(num_chunks), mpierr)
      if (wantDebug) call obj%timer%stop("mpi_communication")

      num_hh_vecs = num_hh_vecs + local_size
    endif
    nx = nx - nb
    if (n == block_limits(nt+1)) then
      nt = nt + 1
    endif
  enddo
  ireq_hhs(:) = MPI_REQUEST_NULL
  ! Buffers for gathering/sending the HH vectors

  allocate(hh_gath(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! gathers HH vectors
  call check_allocate_f("tridiag_band: hh_gath", 341,  istat,  errorMessage)

  allocate(hh_send(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! send buffer for HH vectors
  call check_allocate_f("tridiag_band: hh_send", 344,  istat,  errorMessage)

  hh_gath(:,:,:) = 0.0_rck
  hh_send(:,:,:) = 0.0_rck

  ! Some counters

  allocate(hh_cnt(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_cnt", 352,  istat,  errorMessage)

  allocate(hh_dst(nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: hh_dst", 355,  istat,  errorMessage)

  hh_cnt(:) = 1 ! The first transfomation Vector is always 0 and not calculated at all
  hh_dst(:) = 0 ! PE number for receive
  ireq_ab = MPI_REQUEST_NULL
  ireq_hv = MPI_REQUEST_NULL
  ! Limits for sending

  allocate(snd_limits(0:np_rows,nblocks), stat=istat, errmsg=errorMessage)
  call check_allocate_f("tridiag_band: snd_limits", 366,  istat,  errorMessage)

  do iblk=1,nblocks
    call determine_workload(obj, na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
  enddo


  ! ---------------------------------------------------------------------------
  ! Start of calculations

  na_s = block_limits(my_pe)*nb + 1

  if (my_pe>0 .and. na_s<=na) then
    ! send first column to previous PE
    ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
    ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
    if (wantDebug) call obj%timer%start("mpi_communication")
    call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_COMPLEX8, &
                   int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), ireq_ab, mpierr)
    if (wantDebug) call obj%timer%stop("mpi_communication")
  endif


   do istep=1,na-nblockEnd

     if (my_pe==0) then
       n = MIN(na-na_s,nb) ! number of rows to be reduced
       hv(:) = 0.0_rck
       hd(:) = 0.0_rck
       tau = 0.0_rck

       ! Transform first column of remaining matrix
       ! Opposed to the real case, the last step (istep=na-1) is needed here for making
       ! the last subdiagonal element a real number

         vnorm2 = sum(real(ab(3:n+1,na_s-n_off),kind=rk4)**2+aimag(ab(3:n+1,na_s-n_off))**2)
         if (n<2) vnorm2 = 0.0_rk ! Safety only

          call hh_transform_&
               &complex&
               &_&
               &single &
                          (obj, ab(2,na_s-n_off), vnorm2, hf, tau, wantDebug)

          hv(1) = 1.0_rck
          hv(2:n) = ab(3:n+1,na_s-n_off)*hf

       d(istep) = real(ab(1,na_s-n_off), kind=rk)
       e(istep) = real(ab(2,na_s-n_off), kind=rk)

       if (istep == na-1) then

         d(na) = real(ab(1,na_s+1-n_off),kind=rk)
         e(na) = 0.0_rck
       endif
     else
       if (na>na_s) then
         ! Receive Householder Vector from previous task, from PE owning subdiagonal


         if (wantDebug) call obj%timer%start("mpi_communication")

         call mpi_recv(hv, int(nb,kind=MPI_KIND), MPI_COMPLEX8, &
                       int(my_pe-1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
         if (wantDebug) call obj%timer%stop("mpi_communication")


         tau = hv(1)
         hv(1) = 1.0_rck
       endif
      endif

      na_s = na_s+1
      if (na_s-n_off > nb) then
        ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
        ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rck
        n_off = n_off + nb
      endif


      do iblk=1,nblocks
        ns = na_s + (iblk-1)*nb - n_off ! first column in block
        ne = ns+nb-1                    ! last column in block

        if (ns+n_off>na) exit

        ! Store Householder Vector for back transformation

        hh_cnt(iblk) = hh_cnt(iblk) + 1

        hh_gath(1   ,hh_cnt(iblk),iblk) = tau
        hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

        if (hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
          ! Wait for last transfer to finish
          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")
          ! Copy vectors into send buffer
          hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
          ! Send to destination

          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_isend(hh_send(1,1,iblk), int(nb*hh_cnt(iblk),kind=MPI_KIND), MPI_COMPLEX8, &
                         global_id(hh_dst(iblk), mod(iblk+block_limits(my_pe)-1,np_cols)), &
                         int(10+iblk,kind=MPI_KIND), int(mpi_comm_all,kind=MPI_KIND), ireq_hhs(iblk), mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! Reset counter and increase destination row
          hh_cnt(iblk) = 0
          hh_dst(iblk) = hh_dst(iblk)+1
        endif

        ! The following code is structured in a way to keep waiting times for
        ! other PEs at a minimum, especially if there is only one block.
        ! For this reason, it requests the last column as late as possible
        ! and sends the Householder Vector and the first column as early
        ! as possible.
        nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
        nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
                                      ! Note that nr>=0 implies that diagonal block is full (nc==nb)!

        ! Multiply diagonal block and subdiagonal block with Householder Vector

        if (iblk==nblocks .and. nc==nb) then

          ! We need the last column from the next PE.
          ! First do the matrix multiplications without last column ...

          ! Diagonal block, the contribution of the last element is added below!
          ab(1,ne) = 0.0_rck

          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")

            call CHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            ! Subdiagonal block
            if (nr>0) call CGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")

            call CHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                             hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            ! Subdiagonal block
            if (nr>0) call CGEMV('N', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                          tau, ab(nb+1,ns), int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, &
                                          ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif

          ! ... then request last column ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_recv(ab(1,ne), int(nb+1,kind=MPI_KIND), MPI_COMPLEX8,  &
                       int(my_pe+1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                       MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and complete the result
          hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
          hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

        else

          ! Normal matrix multiply
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            call CHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                                hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)

            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            call CHEMV('L', int(nc,kind=BLAS_KIND), tau, ab(1,ns), int(2*nb-1,kind=BLAS_KIND), &
                               hv, 1_BLAS_KIND, ZERO, hd, 1_BLAS_KIND)
            if (nr>0) call CGEMV('N', int(nr,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), tau, ab(nb+1,ns), &
                                        int(2*nb-1,kind=BLAS_KIND), hv, 1_BLAS_KIND, ZERO, hs, 1_BLAS_KIND)
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Calculate first column of subdiagonal block and calculate new
        ! Householder transformation for this column
        hv_new(:) = 0.0_rck ! Needed, last rows must be 0 for nr < nb
        tau_new = 0.0_rck
        if (nr>0) then

          ! complete (old) Householder transformation for first column

          ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

          ! calculate new Householder transformation ...
          if (nr>1) then
            vnorm2 = sum(real(ab(nb+2:nb+nr,ns),kind=rk4)**2+aimag(ab(nb+2:nb+nr,ns))**2)

            call hh_transform_&
                &complex&
                &_&
                &single &
               (obj, ab(nb+1,ns), vnorm2, hf, tau_new, wantDebug)
            hv_new(1) = 1.0_rck
            hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
            ab(nb+2:,ns) = 0.0_rck
          endif ! nr > 1

          ! ... and send it away immediatly if this is the last block

          if (iblk==nblocks) then
            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            hv_s(1) = tau_new
            hv_s(2:) = hv_new(2:)

            if (wantDebug) call obj%timer%start("mpi_communication")
            call mpi_isend(hv_s, int(nb,kind=MPI_KIND), MPI_COMPLEX8, &
                           int(my_pe+1,kind=MPI_KIND), 2_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                           ireq_hv, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

          endif

        endif

        ! Transform diagonal block
        x = dot_product(hv(1:nc),hd(1:nc))*conjg(tau)

          hd(1:nc) = hd(1:nc) - 0.5_rk*x*hv(1:nc)
        if (my_pe>0 .and. iblk==1) then

          ! The first column of the diagonal block has to be send to the previous PE
          ! Calculate first column only ...
          ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*conjg(hv(1)) - hv(1:nc)*conjg(hd(1))
          ! ... send it away ...
          if (wantDebug) call obj%timer%start("mpi_communication")
          call mpi_wait(ireq_ab, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ab_s(1:nb+1) = ab(1:nb+1,ns)

          if (wantDebug) call obj%timer%start("mpi_communication")

          call mpi_isend(ab_s, int(nb+1,kind=MPI_KIND), MPI_COMPLEX8, &
                         int(my_pe-1,kind=MPI_KIND), 1_MPI_KIND, int(mpi_comm_all,kind=MPI_KIND), &
                         ireq_ab, mpierr)
          if (wantDebug) call obj%timer%stop("mpi_communication")

          ! ... and calculate remaining columns with rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            if (nc>1) call CHER2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                         hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            if (wantDebug) call obj%timer%stop("mkl_offload")
          else
            if (wantDebug) call obj%timer%start("blas")
            if (nc>1) call CHER2('L', int(nc-1,kind=BLAS_KIND), -ONE, hd(2), 1_BLAS_KIND, &
                                       hv(2), 1_BLAS_KIND, ab(1,ns+1), int(2*nb-1,kind=BLAS_KIND) )
            if (wantDebug) call obj%timer%stop("blas")
          endif

        else
          ! No need to  send, just a rank-2 update
          if (useIntelGPU) then
            if (wantDebug) call obj%timer%start("mkl_offload")
            call CHER2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND, hv, 1_BLAS_KIND, &
                                ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("mkl_offload")

          else
            if (wantDebug) call obj%timer%start("blas")
            call CHER2('L', int(nc,kind=BLAS_KIND), -ONE, hd, 1_BLAS_KIND, hv, 1_BLAS_KIND, &
                                ab(1,ns), int(2*nb-1,kind=BLAS_KIND))
            if (wantDebug) call obj%timer%stop("blas")
          endif
        endif

        ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

        if (nr>0) then
          if (nr>1) then
            if (useIntelGPU) then
              if (wantDebug) call obj%timer%start("mkl_offload")
              call CGEMV('C', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("mkl_offload")
            else
              if (wantDebug) call obj%timer%start("blas")
              call CGEMV('C', int(nr,kind=BLAS_KIND), int(nb-1,kind=BLAS_KIND), &
                                  tau_new, ab(nb,ns+1), int(2*nb-1,kind=BLAS_KIND), &
                                  hv_new, 1_BLAS_KIND, ZERO, h(2), 1_BLAS_KIND)
              if (wantDebug) call obj%timer%stop("blas")
            endif
            x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
            h(2:nb) = h(2:nb) - x*hv(2:nb)
            ! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update
            do i=2,nb
              ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*conjg(h(i)) - hs(1:nr)*conjg(hv(i))
            enddo
          else
            ! No double Householder transformation for nr=1, just complete the row
            do i=2,nb
              ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*conjg(hv(i))
            enddo
          endif
        endif

        ! Use new HH Vector for the next block
        hv(:) = hv_new(:)
        tau = tau_new

      enddo


  enddo ! istep

  ! Finish the last outstanding requests


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
  call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

  call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
  call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")


  if (wantDebug) call obj%timer%start("mpi_communication")
  call mpi_barrier(int(mpi_comm_all,kind=MPI_KIND),mpierr)
  if (wantDebug) call obj%timer%stop("mpi_communication")
  deallocate(ab, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ab", 1363,  istat,  errorMessage)

  deallocate(ireq_hhr, ireq_hhs, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: ireq_hhr", 1366,  istat,  errorMessage)

  deallocate(hh_cnt, hh_dst, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_dst", 1369,  istat,  errorMessage)

  deallocate(hh_gath, hh_send, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: hh_gath", 1372,  istat,  errorMessage)

  deallocate(limits, snd_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: limits", 1375,  istat,  errorMessage)

  deallocate(block_limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: block_limits", 1378,  istat,  errorMessage)

  deallocate(global_id, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("tridiag_band: global_id", 1381,  istat,  errorMessage)

  call obj%timer%stop("tridiag_band_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

  ! intel compiler bug makes these ifdefs necessary
end subroutine tridiag_band_complex_&
&single

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






subroutine trans_ev_tridi_to_band_&
&complex&
&_&
&single &
(obj, na, nev, nblk, nbw, q, ldq, matrixCols,         &
 hh_trans, mpi_comm_rows, mpi_comm_cols, wantDebug, useGPU, max_threads, success, &
 kernel)

!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_real/complex:
!  Transforms the eigenvectors of a tridiagonal matrix back to the eigenvectors of the band matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nev         Number eigenvectors to compute (= columns of matrix q)
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nb          semi bandwith
!
!  q           On input: Eigenvectors of tridiagonal matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!  matrixCols  local columns of matrix q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------
  use elpa_abstract_impl
  use elpa2_workload
  use pack_unpack_cpu
  use pack_unpack_gpu
  use compute_hh_trafo
  use elpa_gpu
  use precision
  use, intrinsic :: iso_c_binding
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
  class(elpa_abstract_impl_t), intent(inout)   :: obj
  logical, intent(in)                          :: useGPU

  integer(kind=ik), intent(in)                 :: kernel
  integer(kind=ik), intent(in)                 :: na, nev, nblk, nbw, ldq, matrixCols, mpi_comm_rows, mpi_comm_cols

  complex(kind=rck), target              :: q(ldq,*)

  complex(kind=rck), intent(in),target   :: hh_trans(:,:)
  integer(kind=c_intptr_t)                     :: hh_trans_dev
  type(c_ptr)                                  :: hh_trans_mpi_dev
  complex(kind=rck), pointer             :: hh_trans_mpi_fortran_ptr(:,:)


  integer(kind=ik)                             :: np_rows, my_prow, np_cols, my_pcol
  integer(kind=MPI_KIND)                       :: np_rowsMPI, my_prowMPI, np_colsMPI, my_pcolMPI
  integer(kind=ik)                             :: i, j, ip, sweep, nbuf, l_nev, a_dim2
  integer(kind=ik)                             :: current_n, current_local_n, current_n_start, current_n_end
  integer(kind=ik)                             :: next_n, next_local_n, next_n_start, next_n_end
  integer(kind=ik)                             :: bottom_msg_length, top_msg_length, next_top_msg_length
  integer(kind=ik)                             :: stripe_width, last_stripe_width, stripe_count
  integer(kind=ik)                             :: num_result_blocks, num_result_buffers, num_bufs_recvd
  integer(kind=ik)                             :: a_off, current_tv_off, max_blk_size
  integer(kind=ik)                             :: src, src_offset, dst, offset, nfact, num_blk
  integer(kind=MPI_KIND)                       :: mpierr

  logical                                      :: flag
  complex(kind=rck), pointer             :: aIntern(:,:,:)
  complex(kind=rck)                      :: a_var

  type(c_ptr)                                  :: aIntern_ptr

  complex(kind=rck), allocatable, target :: row(:)
  
  integer(kind=c_intptr_t)                     :: row_dev
  type(c_ptr)                                  :: row_mpi_dev
  complex(kind=rck), pointer             :: row_mpi_fortran_ptr(:)

  complex(kind=rck), pointer             :: row_group(:,:)

  complex(kind=rck), allocatable, target :: top_border_send_buffer(:,:)
  complex(kind=rck), allocatable, target :: top_border_recv_buffer(:,:)
  complex(kind=rck), allocatable, target :: bottom_border_send_buffer(:,:)
  complex(kind=rck), allocatable, target :: bottom_border_recv_buffer(:,:)

  integer(kind=c_intptr_t)                     :: top_border_recv_buffer_dev, top_border_send_buffer_dev
  type(c_ptr)                                  :: top_border_recv_buffer_mpi_dev, top_border_send_buffer_mpi_dev 
  complex(kind=rck), pointer             :: top_border_recv_buffer_mpi_fortran_ptr(:,:), &
                                                  top_border_send_buffer_mpi_fortran_ptr(:,:)
  integer(kind=c_intptr_t)                     :: bottom_border_send_buffer_dev, bottom_border_recv_buffer_dev
  type(c_ptr)                                  :: bottom_border_send_buffer_mpi_dev, bottom_border_recv_buffer_mpi_dev
  complex(kind=rck), pointer             :: bottom_border_send_buffer_mpi_fortran_ptr(:,:), &
                                                  bottom_border_recv_buffer_mpi_fortran_ptr(:,:)
  type(c_ptr)                                  :: aIntern_mpi_dev
  complex(kind=rck), pointer             :: aIntern_mpi_fortran_ptr(:,:,:)

  integer(kind=c_intptr_t)                     :: aIntern_dev
  integer(kind=c_intptr_t)                     :: bcast_buffer_dev

  type(c_ptr)                                  :: bcast_buffer_mpi_dev
  complex(kind=rck), pointer             :: bcast_buffer_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: num
  integer(kind=c_intptr_t)                     :: dev_offset, dev_offset_1
  integer(kind=c_intptr_t)                     :: row_group_dev

  type(c_ptr)                                  :: row_group_mpi_dev
  complex(kind=rck), pointer             :: row_group_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: q_dev
  type(c_ptr)                                  :: q_mpi_dev
  complex(kind=rck), pointer             :: q_mpi_fortran_ptr(:,:)

  integer(kind=c_intptr_t)                     :: hh_tau_dev
  integer(kind=ik)                             :: row_group_size, unpack_idx

  type(c_ptr)                                  :: row_group_host, bcast_buffer_host

  integer(kind=ik)                             :: n_times
  integer(kind=ik)                             :: chunk, this_chunk

  complex(kind=rck), allocatable,target  :: result_buffer(:,:,:)
  integer(kind=c_intptr_t)                     :: result_buffer_dev

  type(c_ptr)                                  :: result_buffer_mpi_dev
  complex(kind=rck), pointer             :: result_buffer_mpi_fortran_ptr(:,:,:)

  complex(kind=rck), pointer             :: bcast_buffer(:,:)

  integer(kind=ik)                             :: n_off

  integer(kind=MPI_KIND), allocatable          :: result_send_request(:), result_recv_request(:)
  integer(kind=ik), allocatable                :: limits(:)
  integer(kind=MPI_KIND), allocatable          :: top_send_request(:), bottom_send_request(:)
  integer(kind=MPI_KIND), allocatable          :: top_recv_request(:), bottom_recv_request(:)

  ! MPI send/recv tags, arbitrary

  integer(kind=ik), parameter                  :: bottom_recv_tag = 111
  integer(kind=ik), parameter                  :: top_recv_tag    = 222
  integer(kind=ik), parameter                  :: result_recv_tag = 333

  integer(kind=ik), intent(in)                 :: max_threads



  ! Just for measuring the kernel performance
  real(kind=c_double)                          :: kernel_time, kernel_time_recv ! MPI_WTIME always needs double
  ! long integer
  integer(kind=lik)                            :: kernel_flops, kernel_flops_recv

  logical, intent(in)                          :: wantDebug
  logical                                      :: success
  integer(kind=ik)                             :: istat, print_flops
  character(200)                               :: errorMessage
  character(20)                                :: gpuString
  logical                                      :: successGPU
  integer(kind=ik)                             :: error
  integer(kind=c_intptr_t), parameter          :: size_of_datatype = size_of_&
                                                                 &single&
                                                                 &_&
                                                                 &complex
  logical, parameter                           :: allComputeOnGPU = .false.

  integer(kind=MPI_KIND)                     :: bcast_request1, allreduce_request1, allreduce_request2, &
                                                allreduce_request3, allreduce_request4
  logical                                    :: useNonBlockingCollectivesCols
  logical                                    :: useNonBlockingCollectivesRows
  integer(kind=c_int)                        :: non_blocking_collectives_rows, non_blocking_collectives_cols

  if(useGPU) then
    gpuString = "_gpu"
  else
    gpuString = ""
  endif

  call obj%timer%start("trans_ev_tridi_to_band_&
  &complex&
  &" // &
  &"_single" //&
  gpuString)

  call obj%get("nbc_row_elpa2_tridi_to_band", non_blocking_collectives_rows, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for rows in elpa2_tridi_to_band. Aborting..."
    stop
  endif

  call obj%get("nbc_col_elpa2_tridi_to_band", non_blocking_collectives_cols, error)
  if (error .ne. ELPA_OK) then
    print *,"Problem setting option for non blocking collectives for cols in elpa2_tridi_to_band. Aborting..."
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

  n_times = 0
  if (useGPU) then
    unpack_idx = 0
    row_group_size = 0
  endif

  success = .true.
  kernel_time = 0.0
  kernel_flops = 0

  if (wantDebug) call obj%timer%start("mpi_communication")
  call MPI_Comm_rank(int(mpi_comm_rows,kind=MPI_KIND) , my_prowMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_rows,kind=MPI_KIND) , np_rowsMPI , mpierr)
  call MPI_Comm_rank(int(mpi_comm_cols,kind=MPI_KIND) , my_pcolMPI , mpierr)
  call MPI_Comm_size(int(mpi_comm_cols,kind=MPI_KIND) , np_colsMPI , mpierr)

  my_prow = int(my_prowMPI,kind=c_int)
  my_pcol = int(my_pcolMPI,kind=c_int)
  np_rows = int(np_rowsMPI,kind=c_int)
  np_cols = int(np_colsMPI,kind=c_int)

  if (wantDebug) call obj%timer%stop("mpi_communication")

  if (mod(nbw,nblk)/=0) then
    if (my_prow==0 .and. my_pcol==0) then
      if (wantDebug) then
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &complex&
                            &: ERROR: nbw=',nbw,', nblk=',nblk
        write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                            &complex&
                            &: band backtransform works only for nbw==n*nblk'
      endif
      success = .false.
      return
    endif
  endif

  nfact = nbw / nblk


  ! local number of eigenvectors
  l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

  if (l_nev==0) then
    stripe_width = 0
    stripe_count = 0
    last_stripe_width = 0

  else ! l_nev


    ! Suggested stripe width is 48 since 48*64 real*8 numbers should fit into
    ! every primary cache
    ! Suggested stripe width is 48 - should this be reduced for the complex case ???

    if (useGPU) then
      stripe_width = 1024 ! Must be a multiple of 4
      stripe_count = (l_nev - 1) / stripe_width + 1

    else ! useGPU

      call obj%get("stripewidth_complex",stripe_width, error)

      !stripe_width = 48 ! Must be a multiple of 4

      stripe_count = (l_nev-1)/stripe_width + 1

      ! Adapt stripe width so that last one doesn't get too small

      stripe_width = (l_nev-1)/stripe_count + 1



     if (kernel .eq. ELPA_2STAGE_COMPLEX_AVX512_BLOCK1 .or. &
         kernel .eq. ELPA_2STAGE_COMPLEX_AVX512_BLOCK2 .or. &
         kernel .eq. ELPA_2STAGE_COMPLEX_SVE512_BLOCK1 .or. &
         kernel .eq. ELPA_2STAGE_COMPLEX_SVE512_BLOCK2  &
         ) then

       stripe_width = ((stripe_width+15)/16)*16 ! Must be a multiple of 8 because of AVX-512 memory alignment of 64 bytes
                                       ! (8 * sizeof(float complex) == 64)

     else
       stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 4 because of AVX/SSE memory alignment of 32 bytes
                                       ! (4 * sizeof(float complex) == 32)
     endif
   endif ! useGPU

   last_stripe_width = l_nev - (stripe_count-1)*stripe_width

  endif ! l_nev

  ! Determine the matrix distribution at the beginning

  allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: limits", 604,  istat,  errorMessage)
  call determine_workload(obj,na, nbw, np_rows, limits)

  max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

  a_dim2 = max_blk_size + nbw

  if (useGPU) then

    if (allComputeOnGPU) then
      if (wantDebug) call obj%timer%start("cuda_memcpy")

      successGPU = gpu_malloc(q_dev, ldq*matrixCols* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: q_dev", 617,  successGPU)
      successGPU =  gpu_memcpy(q_dev, int(loc(q(1,1)),kind=c_intptr_t),  &
                               ldq*matrixCols * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q -> q_dev", 621,  successGPU)

      ! associate with c_ptr
      q_mpi_dev = transfer(q_dev, q_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(q_mpi_dev, q_mpi_fortran_ptr, &
                       [ldq,matrixCols])
      if (wantDebug) call obj%timer%stop("cuda_memcpy")

      successGPU = gpu_malloc(hh_trans_dev, size(hh_trans,dim=1)*size(hh_trans,dim=2)* size_of_datatype)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 631,  successGPU)
      ! associate with c_ptr
      hh_trans_mpi_dev = transfer(hh_trans_dev, hh_trans_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(hh_trans_mpi_dev, hh_trans_mpi_fortran_ptr, &
                       [size(hh_trans,dim=1),size(hh_trans,dim=2)])
      successGPU =  gpu_memcpy(c_loc(hh_trans_mpi_fortran_ptr(1,1)),  &
                               c_loc(hh_trans(1,1)), &
                               size(hh_trans,dim=1)*size(hh_trans,dim=2) * size_of_datatype, &
                               gpuMemcpyHostToDevice)
      call check_memcpy_GPU_f("trans_ev_tridi_to_band: hh_trans -> hh_trans_dev", 641,  successGPU)

    endif ! allComputeOnGPU

    num =  (stripe_width*a_dim2*stripe_count)* size_of_datatype
    successGPU = gpu_malloc(aIntern_dev, stripe_width*a_dim2*stripe_count* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 647,  successGPU)

    successGPU = gpu_memset(aIntern_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 650,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      aIntern_mpi_dev = transfer(aIntern_dev, aIntern_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(aIntern_mpi_dev, aIntern_mpi_fortran_ptr, &
                       [stripe_width,a_dim2,stripe_count])
    endif ! allComputeOnGPU

    ! "row_group" and "row_group_dev" are needed for GPU optimizations
    successGPU = gpu_malloc_host(row_group_host,l_nev*nblk*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 667,  successGPU)
    call c_f_pointer(row_group_host, row_group, (/l_nev,nblk/))

    row_group(:, :) = 0.0_rck
    num =  (l_nev*nblk)* size_of_datatype
    successGPU = gpu_malloc(row_group_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 673,  successGPU)

    successGPU = gpu_memset(row_group_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_group_dev", 676,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      row_group_mpi_dev = transfer(row_group_dev, row_group_mpi_dev)
      ! and associate a fortran pointer
      call c_f_pointer(row_group_mpi_dev, row_group_mpi_fortran_ptr, &
                       [l_nev,nblk])
    endif ! allComputeOnGPU

  else ! GPUs are not used



    if (posix_memalign(aIntern_ptr, 64_c_intptr_t, stripe_width*a_dim2*stripe_count*  &
        C_SIZEOF(a_var)) /= 0) then
      print *,"trans_ev_tridi_to_band_real: error when allocating aIntern"//errorMessage
      stop 1
    endif

    call c_f_pointer(aIntern_ptr, aIntern,[stripe_width,a_dim2,stripe_count] )
    !allocate(aIntern(stripe_width,a_dim2,stripe_count), stat=istat, errmsg=errorMessage)

    aIntern(:,:,:) = 0.0_rck
  endif !useGPU

  allocate(row(l_nev), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: row", 723,  istat,  errorMessage)

  row(:) = 0.0_rck

  if (useGPU .and. allComputeOnGPU) then
    num =  (l_nev)* size_of_datatype
    successGPU = gpu_malloc(row_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: row_dev", 730,  successGPU)

    successGPU = gpu_memset(row_dev , 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: row_dev", 733,  successGPU)

    ! associate with c_ptr
    row_mpi_dev = transfer(row_dev, row_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(row_mpi_dev, row_mpi_fortran_ptr, &
                     [l_nev])
  endif

  ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
  ! and transpose the matrix using stripes of given stripe_width for cache blocking.

  ! The peculiar way it is done below is due to the fact that the last row should be
  ! ready first since it is the first one to start below



  if (wantDebug) call obj%timer%start("ip_loop")
  do ip = np_rows-1, 0, -1
    if (my_prow == ip) then
      ! Receive my rows which have not yet been received
      src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
      do i=limits(ip)+1,limits(ip+1)
        src = mod((i-1)/nblk, np_rows)

        if (src < my_prow) then

          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
            &complex&
            &_gpu_&
            &single &
                          ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                                      stripe_width, last_stripe_width, a_dim2, l_nev,&
                                      row_group_size, nblk, unpack_idx, &
                                       i - limits(ip), .false., wantDebug, allComputeOnGPU)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                         int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")

            call unpack_row_&
                &complex&
                &_cpu_&
                &single &
                (obj,aIntern, row,i-limits(ip), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        elseif (src == my_prow) then

          src_offset = src_offset+1

          if (useGPU) then

            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &complex&
                 &_gpu_&
                 &single &
             ( obj, &
                          row_group, row_group_dev, aIntern_dev, stripe_count, &
                          stripe_width, last_stripe_width, a_dim2, l_nev,&
                          row_group_size, nblk, unpack_idx, &
                          i - limits(ip), .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_CCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_group_mpi_fortran_ptr(1,row_group_size)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            else ! allComputeOnGPU
              row_group(:, row_group_size) = q(src_offset, 1:l_nev)
            endif ! allComputeOnGPU
          else ! useGPU
            row(:) = q(src_offset, 1:l_nev)
          endif ! useGPU


          if (useGPU) then

          else
            call unpack_row_&
                 &complex&
                 &_cpu_&
                 &single &
                            (obj,aIntern, row,i-limits(ip),  stripe_count, stripe_width, last_stripe_width)
          endif


        endif ! src == my_prow
      enddo ! i=limits(ip)+1,limits(ip+1)


      ! Send all rows which have not yet been send


      if (allComputeOnGPU .and. useGPU) then
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1

              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_CCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1011,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

              if (wantDebug) call obj%timer%start("cuda_mpi_communication")
              call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

            endif
          enddo
        enddo
      else !  allComputeOnGPU
        src_offset = 0
        do dst = 0, ip-1
          do i=limits(dst)+1,limits(dst+1)
            if (mod((i-1)/nblk, np_rows) == my_prow) then
              src_offset = src_offset+1
              row(:) = q(src_offset, 1:l_nev)

              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                                int(dst,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            endif
          enddo
        enddo
      endif ! allComputeOnGPU

    else if (my_prow < ip) then

      ! Send all rows going to PE ip
      if (allComputeOnGPU .and. useGPU) then
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1

            if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
            call gpublas_CCOPY(l_nev, c_loc(q_mpi_fortran_ptr(src_offset,1)), ldq, &
                                          c_loc(row_mpi_fortran_ptr(1)), 1)
            if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")

            ! is there a way to avoid this device_synchronize ?
            if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
            successGPU = gpu_devicesynchronize()
            call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 1061,  successGPU)
            if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")

            if (wantDebug) call obj%timer%start("cuda_mpi_communication")
            call MPI_Send(row_mpi_fortran_ptr, int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                          int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("cuda_mpi_communication")

          endif
        enddo
      else ! allComputeOnGPU 
        src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
        do i=limits(ip)+1,limits(ip+1)
          src = mod((i-1)/nblk, np_rows)
          if (src == my_prow) then
            src_offset = src_offset+1
            row(:) = q(src_offset, 1:l_nev)
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Send(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                        int(ip,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          endif
        enddo
      endif  ! allComputeOnGPU

      ! Receive all rows from PE ip
      do i=limits(my_prow)+1,limits(my_prow+1)
        src = mod((i-1)/nblk, np_rows)
        if (src == ip) then
          if (useGPU) then
            ! An unpacking of the current row group may occur before queuing the next row

            call unpack_and_prepare_row_group_&
                 &complex&
                 &_gpu_&
                 &single&
                 &( obj, &
                 row_group, row_group_dev, aIntern_dev, stripe_count,  &
                 stripe_width, last_stripe_width, a_dim2, l_nev,       &
                 row_group_size, nblk, unpack_idx,                     &
                 i - limits(my_prow), .false., wantDebug, allComputeOnGPU)

            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Recv(row_group(:, row_group_size), int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                          int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")

          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Recv(row, int(l_nev,kind=MPI_KIND), MPI_COMPLEX8, &
                             int(src,kind=MPI_KIND), 0_MPI_KIND, int(mpi_comm_rows,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
            call unpack_row_&
                 &complex&
                 &_cpu_&
                 &single &
                 (obj,aIntern, row,i-limits(my_prow), stripe_count, stripe_width, last_stripe_width)
          endif ! useGPU


        endif
      enddo ! i=limits(my_prow)+1,limits(my_prow+1)
    endif ! (my_prow < ip)
  enddo ! ip = np_rows-1, 0, -1

  if (wantDebug) call obj%timer%stop("ip_loop")
 
  if (wantDebug) call obj%timer%start("allocate")

  if (useGPU) then
    ! Force an unpacking of all remaining rows that haven't been unpacked yet

    call unpack_and_prepare_row_group_&
         &complex&
         &_gpu_&
         &single&
         &( obj, &
         row_group, row_group_dev, aIntern_dev, stripe_count, &
         stripe_width, last_stripe_width, &
         a_dim2, l_nev, row_group_size, nblk, unpack_idx,     &
         -1, .true., wantDebug, allComputeOnGPU)

  endif

  ! Set up result buffer queue

  num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

  num_result_buffers = 4*nfact
  allocate(result_buffer(l_nev,nblk,num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_buffer", 1266,  istat,  errorMessage)

  allocate(result_send_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_send_request", 1269,  istat,  errorMessage)

  allocate(result_recv_request(num_result_buffers), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: result_recv_request", 1272,  istat,  errorMessage)

  if (useGPU .and. allComputeOnGPU) then
    num_result_buffers = 4*nfact
    num =  (l_nev*nblk*num_result_buffers)* size_of_datatype
    successGPU = gpu_malloc(result_buffer_dev, num* size_of_datatype)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 1278,  successGPU)

    ! associate with c_ptr
    result_buffer_mpi_dev = transfer(result_buffer_dev, result_buffer_mpi_dev)
    ! and associate a fortran pointer
    call c_f_pointer(result_buffer_mpi_dev, result_buffer_mpi_fortran_ptr, &
                     [l_nev,nblk,num_result_buffers])
  endif

  result_send_request(:) = MPI_REQUEST_NULL
  result_recv_request(:) = MPI_REQUEST_NULL

  ! Queue up buffers
  if (wantDebug) call obj%timer%start("host_mpi_communication")
  if (my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
    do j = 1, min(num_result_buffers, num_result_blocks)
      if (useGPU) then
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX8,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      else ! useGPU
        call MPI_Irecv(result_buffer(1,1,j), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX8,     &
                       0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),          &
                      result_recv_request(j), mpierr)
      endif ! useGPU
    enddo
  endif
  if (wantDebug) call obj%timer%stop("host_mpi_communication")

  num_bufs_recvd = 0 ! No buffers received yet

  ! Initialize top/bottom requests

  allocate(top_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_send_request", 1335,  istat,  errorMessage)

  allocate(top_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_recv_request", 1338,  istat,  errorMessage)

  allocate(bottom_send_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_send_request", 1341,  istat,  errorMessage)

  allocate(bottom_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_recv_request", 1344,  istat,  errorMessage)

  top_send_request(:) = MPI_REQUEST_NULL
  top_recv_request(:) = MPI_REQUEST_NULL
  bottom_send_request(:) = MPI_REQUEST_NULL
  bottom_recv_request(:) = MPI_REQUEST_NULL


  allocate(top_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 1374,  istat,  errorMessage)

  allocate(top_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1377,  istat,  errorMessage)

  allocate(bottom_border_send_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1380,  istat,  errorMessage)

  allocate(bottom_border_recv_buffer(stripe_width*nbw, stripe_count), stat=istat, errmsg=errorMessage)
  call check_allocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1383,  istat,  errorMessage)

  top_border_send_buffer(:,:) = 0.0_rck
  top_border_recv_buffer(:,:) = 0.0_rck
  bottom_border_send_buffer(:,:) = 0.0_rck
  bottom_border_recv_buffer(:,:) = 0.0_rck


  if (useGPU) then
    if (allComputeOnGPU) then
      ! top_border_recv_buffer and top_border_send_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(top_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1401,  successGPU)

      successGPU = gpu_malloc(top_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1404,  successGPU)

      successGPU = gpu_memset(top_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 1407,  successGPU)

      successGPU = gpu_memset(top_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 1410,  successGPU)

      ! associate with c_ptr
      top_border_recv_buffer_mpi_dev = transfer(top_border_recv_buffer_dev, top_border_recv_buffer_mpi_dev)
      top_border_send_buffer_mpi_dev = transfer(top_border_send_buffer_dev, top_border_send_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(top_border_recv_buffer_mpi_dev, top_border_recv_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])
      call c_f_pointer(top_border_send_buffer_mpi_dev, top_border_send_buffer_mpi_fortran_ptr, [stripe_width*nbw, stripe_count])

      ! bottom_border_send_buffer and bottom_border_recv_buffer
      num =  ( stripe_width*nbw*stripe_count) * size_of_datatype
      successGPU = gpu_malloc(bottom_border_send_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1433,  successGPU)
      successGPU = gpu_malloc(bottom_border_recv_buffer_dev, num)
      call check_alloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1435,  successGPU)

      successGPU = gpu_memset(bottom_border_send_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 1438,  successGPU)
      successGPU = gpu_memset(bottom_border_recv_buffer_dev, 0, num)
      call check_memset_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 1440,  successGPU)

      ! associate with c_ptr
      bottom_border_send_buffer_mpi_dev = transfer(bottom_border_send_buffer_dev, bottom_border_send_buffer_mpi_dev)
      bottom_border_recv_buffer_mpi_dev = transfer(bottom_border_recv_buffer_dev, bottom_border_recv_buffer_mpi_dev)
      ! and create a fortran pointer
      call c_f_pointer(bottom_border_send_buffer_mpi_dev, bottom_border_send_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
      call c_f_pointer(bottom_border_recv_buffer_mpi_dev, bottom_border_recv_buffer_mpi_fortran_ptr, &
                       [stripe_width*nbw, stripe_count])
    endif ! allComputeOnGPU

    successGPU = gpu_host_register(int(loc(top_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 1462,  successGPU)

    successGPU = gpu_host_register(int(loc(top_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 1467,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_send_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 1472,  successGPU)

    successGPU = gpu_host_register(int(loc(bottom_border_recv_buffer),kind=c_intptr_t), &
                  stripe_width*nbw* stripe_count * size_of_datatype,&
                  gpuHostRegisterDefault)
    call check_host_register_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 1477,  successGPU)
  endif ! useGPU


  ! Initialize broadcast buffer

  if (useGPU) then
    successGPU = gpu_malloc_host(bcast_buffer_host,nbw*max_blk_size*size_of_datatype)
    call check_host_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 1485,  successGPU)
    call c_f_pointer(bcast_buffer_host, bcast_buffer, (/nbw,max_blk_size/))
  else
    allocate(bcast_buffer(nbw, max_blk_size), stat=istat, errmsg=errorMessage)
    call check_allocate_f("trans_ev_tridi_to_band: bcast_buffer", 1489,  istat,  errorMessage)
  endif

  bcast_buffer = 0.0_rck

  if (useGPU) then
    num =  ( nbw * max_blk_size) * size_of_datatype
    successGPU = gpu_malloc(bcast_buffer_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1497,  successGPU)

    if (allComputeOnGPU) then
      ! associate with c_ptr
      bcast_buffer_mpi_dev = transfer(bcast_buffer_dev, bcast_buffer_mpi_dev)
      call c_f_pointer(bcast_buffer_mpi_dev, bcast_buffer_mpi_fortran_ptr, &
                      [nbw, max_blk_size])
    endif ! allComputeOnGPU

    successGPU = gpu_memset( bcast_buffer_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1507,  successGPU)

    num =  (max_blk_size)* size_of_datatype
    successGPU = gpu_malloc( hh_tau_dev, num)
    call check_alloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1511,  successGPU)

    successGPU = gpu_memset( hh_tau_dev, 0, num)
    call check_memset_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 1514,  successGPU)
  endif ! useGPU

  current_tv_off = 0 ! Offset of next row to be broadcast

  ! ------------------- start of work loop -------------------

  a_off = 0 ! offset in aIntern (to avoid unnecessary shifts)

  top_msg_length = 0
  bottom_msg_length = 0
  if (wantDebug) call obj%timer%stop("allocate")

  if (wantDebug) call obj%timer%start("sweep_loop")
  do sweep = 0, (na-1)/nbw

    current_n = na - sweep*nbw
    call determine_workload(obj,current_n, nbw, np_rows, limits)
    current_n_start = limits(my_prow)
    current_n_end   = limits(my_prow+1)
    current_local_n = current_n_end - current_n_start

    next_n = max(current_n - nbw, 0)
    call determine_workload(obj,next_n, nbw, np_rows, limits)
    next_n_start = limits(my_prow)
    next_n_end   = limits(my_prow+1)
    next_local_n = next_n_end - next_n_start

    if (next_n_end < next_n) then
      bottom_msg_length = current_n_end - next_n_end
    else
      bottom_msg_length = 0
    endif

    if (next_local_n > 0) then
      next_top_msg_length = current_n_start - next_n_start
    else
      next_top_msg_length = 0
    endif

    if (sweep==0 .and. current_n_end < current_n .and. l_nev > 0) then
      if (wantDebug) call obj%timer%start("host_mpi_communication")
      do i = 1, stripe_count


        if (useGPU) then
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        else !useGPU
          call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), &
                         int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                         bottom_recv_request(i), mpierr)
        endif !useGPU

      enddo
      if (wantDebug) call obj%timer%stop("host_mpi_communication")
    endif

    if (current_local_n > 1) then
      if (useGPU .and. allComputeOnGPU) then
        if (my_pcol == mod(sweep,np_cols)) then
          if (wantDebug) call obj%timer%start("cuda_memcpy")
          successGPU =  gpu_memcpy(c_loc(bcast_buffer_mpi_fortran_ptr(1,1)), &
                                   c_loc(hh_trans_mpi_fortran_ptr(1,current_tv_off+1)),  &
                                     size(hh_trans,dim=1) * (current_tv_off+current_local_n-(current_tv_off+1)+1) * &   
                                     size_of_datatype, &
                                     gpuMemcpyDeviceToDevice)
          call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1646,  successGPU)
          if (wantDebug) call obj%timer%stop("cuda_memcpy")
          current_tv_off = current_tv_off + current_local_n
        endif
      else !useGPU
        if (my_pcol == mod(sweep,np_cols)) then
          bcast_buffer(:,1:current_local_n) =    &
          hh_trans(:,current_tv_off+1:current_tv_off+current_local_n)
          current_tv_off = current_tv_off + current_local_n
        endif
      endif ! useGPU

      if (useNonBlockingCollectivesCols) then
        if (wantDebug) call obj%timer%start("mpi_nbc_communication")
        call mpi_ibcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_COMPLEX8, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), bcast_request1, mpierr)
        call mpi_wait(bcast_request1, MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("mpi_nbc_communication")
      else
        if (wantDebug) call obj%timer%start("mpi_communication")
        call mpi_bcast(bcast_buffer, int(nbw*current_local_n,kind=MPI_KIND), MPI_COMPLEX8, &
                     int(mod(sweep,np_cols),kind=MPI_KIND), int(mpi_comm_cols,kind=MPI_KIND), mpierr)
        if (wantDebug) call obj%timer%stop("mpi_communication")
      endif
 
      if (useGPU .and. .not.(allComputeOnGPU)) then
        if (wantDebug) call obj%timer%start("memcpy")
        successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
                                 nbw * current_local_n *    &
                                 size_of_datatype, &
                                 gpuMemcpyHostToDevice)
        call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1702,  successGPU)
        if (wantDebug) call obj%timer%stop("memcpy")
      endif ! useGPU

      if (useGPU) then
        !if (.not.(allComputeOnGPU)) then
        !  if (wantDebug) call obj%timer%start("memcpy")
        !  successGPU =  gpu_memcpy(bcast_buffer_dev, int(loc(bcast_buffer(1,1)),kind=c_intptr_t),  &
        !                           nbw * current_local_n *    &
        !                           size_of_datatype, &
        !                           gpuMemcpyHostToDevice)
        !  call check_memcpy_GPU_f("trans_ev_tridi_to_band: bcast_buffer -> bcast_buffer_dev", 1725,  successGPU)
        !  if (wantDebug) call obj%timer%stop("memcpy")
        !endif ! allComputeOnGPU

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &complex&
             &_gpu_&
             &single&
             (bcast_buffer_dev, hh_tau_dev, nbw, &
             current_local_n, .false.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU

    else ! (current_local_n > 1) then

      ! for current_local_n == 1 the one and only HH Vector is 0 and not stored in hh_trans_real/complex
      bcast_buffer(:,1) = 0.0_rck
      if (useGPU) then
        successGPU = gpu_memset(bcast_buffer_dev, 0, nbw * size_of_datatype)
        call check_memset_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 1745,  successGPU)

        if (wantDebug) call obj%timer%start("extract_hh")
        call extract_hh_tau_&
             &complex&
             &_gpu_&
             &single&
             &( &
             bcast_buffer_dev, hh_tau_dev, &
             nbw, 1, .true.)
        if (wantDebug) call obj%timer%stop("extract_hh")
      endif ! useGPU
    endif ! (current_local_n > 1) then

    if (l_nev == 0) cycle

    if (current_local_n > 0) then

      do i = 1, stripe_count

        !wait_b
        if (current_n_end < current_n) then



          if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_recv")
          call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)

          if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_recv")
          n_off = current_local_n+a_off

          if (useGPU) then
            if (allComputeOnGPU) then
              if (wantDebug) call obj%timer%start("cuda_memcpy")
              successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                       c_loc(bottom_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                       stripe_width*nbw* size_of_datatype,      &
                                       gpuMemcpyDeviceToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1859,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_memcpy")
            else ! allComputeOnGPU
              if (wantDebug) call obj%timer%start("memcpy")
              dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * size_of_datatype
              successGPU =  gpu_memcpy( aIntern_dev + dev_offset , &
                                      int(loc(bottom_border_recv_buffer(1,i)),kind=c_intptr_t), &
                                       stripe_width*nbw*  size_of_datatype,    &
                                       gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer -> aIntern_dev", 1868,  successGPU)
              if (wantDebug) call obj%timer%stop("memcpy")
            endif ! allComputeOnGPU
          else ! useGPU
            aIntern(:,n_off+1:n_off+nbw,i) = reshape( &
            bottom_border_recv_buffer(1:stripe_width*nbw,i),(/stripe_width,nbw/))
          endif ! useGPU


          if (next_n_end < next_n) then

            if (useGPU) then
              if (wantDebug) call obj%timer%start("host_mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_communication")
            else ! useGPU
              if (wantDebug) call obj%timer%start("mpi_communication")
              call MPI_Irecv(bottom_border_recv_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), &
                            MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), &
                            int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),      &
                            bottom_recv_request(i), mpierr)
              if (wantDebug) call obj%timer%stop("mpi_communication")
            endif ! useGPU

            endif ! (next_n_end < next_n)
          endif ! (current_n_end < current_n)

          if (current_local_n <= bottom_msg_length + top_msg_length) then

            !wait_t
            if (top_msg_length>0) then


              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! the MPI_IRECV will be done CUDA_AWARE we thus do not need a host to device copy
                  ! However, we have to copy from top_border_recv_buffer_mpi_fortran_ptr to aIntern_mpi_fortran_ptr
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  !Fortran pointer for indexing
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                          c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                          stripe_width*top_msg_length* size_of_datatype,      &
                                          gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2045,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  !             host_offset= (0 + (0 * stripe_width) + ( (i-1) * stripe_width * nbw ) ) * 8
                  successGPU =  gpu_memcpy( aIntern_dev+dev_offset , int(loc(top_border_recv_buffer(1,i)),kind=c_intptr_t),  &
                                             stripe_width*top_msg_length* size_of_datatype,      &
                                             gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2054,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
            endif ! top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
                &complex&
                &_&
                &single&
                &(obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
                a_off, nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
                hh_tau_dev, kernel_flops, kernel_time, n_times, 0, current_local_n, i, &
                last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b        1
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")
            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")

            if (bottom_msg_length>0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer

                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                           c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2248,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")



              if (wantDebug) call obj%timer%start("cuda_aware_device_synchronize")
              successGPU = gpu_devicesynchronize()
              call check_memcpy_GPU_f("trans_ev_tridi_to_band: device_synchronize", 2255,  successGPU)
              if (wantDebug) call obj%timer%stop("cuda_aware_device_synchronize")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy( int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                            stripe_width * bottom_msg_length * size_of_datatype,      &
                                            gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2263,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                      MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                      int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND),  &
                     MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                     int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU
            endif !(bottom_msg_length>0)

          else ! current_local_n <= bottom_msg_length + top_msg_length

            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &complex&
             &_&
             &single&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, &
             current_local_n - bottom_msg_length, bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")

            !send_b
            if (wantDebug) call obj%timer%start("host_mpi_wait_bottom_send")

            call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_wait_bottom_send")
            if (bottom_msg_length > 0) then
              n_off = current_local_n+nbw-bottom_msg_length+a_off

              if (useGPU) then
                if (allComputeOnGPU) then
                  ! send should be done on GPU, send_buffer must be created and filled first
                  ! memcpy from aIntern_dev to bottom_border_send_buffer_dev
                  ! either with two offsets or with indexed pointer
                  if (wantDebug) call obj%timer%start("cuda_memcpy")

                  successGPU =  gpu_memcpy( c_loc(bottom_border_send_buffer_mpi_fortran_ptr(1,i)), &
                                            c_loc(aIntern_mpi_fortran_ptr(1,n_off+1,i)), &
                                             stripe_width * bottom_msg_length * size_of_datatype,      &
                                             gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer_dev", 2516,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(int(loc(bottom_border_send_buffer(1,i)),kind=c_intptr_t), aIntern_dev + dev_offset,  &
                                           stripe_width*bottom_msg_length* size_of_datatype,  &
                                           gpuMemcpyDeviceToHost)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> bottom_border_send_buffer", 2524,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU

                if (wantDebug) call obj%timer%start("host_mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("host_mpi_communication")
          
              else !useGPU
                bottom_border_send_buffer(1:stripe_width*bottom_msg_length,i) = reshape(&
                aIntern(:,n_off+1:n_off+bottom_msg_length,i),(/stripe_width*bottom_msg_length/))
                if (wantDebug) call obj%timer%start("mpi_communication")
                call MPI_Isend(bottom_border_send_buffer(1,i), int(bottom_msg_length*stripe_width,kind=MPI_KIND), &
                           MPI_COMPLEX8, int(my_prow+1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                           int(mpi_comm_rows,kind=MPI_KIND), bottom_send_request(i), mpierr)
                if (wantDebug) call obj%timer%stop("mpi_communication")
              endif !useGPU




            endif ! (bottom_msg_length > 0)
            !compute

            if (wantDebug) call obj%timer%start("compute_hh_trafo")
            call compute_hh_trafo_&
             &complex&
             &_&
             &single&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, top_msg_length, &
             current_local_n-top_msg_length-bottom_msg_length, i, &
             last_stripe_width, kernel)
            if (wantDebug) call obj%timer%stop("compute_hh_trafo")


            !wait_t
            if (top_msg_length>0) then

              if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
              call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)
              if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
              if (useGPU) then
                if (allComputeOnGPU) then
                  if (wantDebug) call obj%timer%start("cuda_memcpy")
                  successGPU =  gpu_memcpy(c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                           c_loc(top_border_recv_buffer_mpi_fortran_ptr(1,i)),  &
                                           stripe_width* top_msg_length* size_of_datatype,      &
                                           gpuMemcpyDeviceToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2739,  successGPU)
                  if (wantDebug) call obj%timer%stop("cuda_memcpy")
                else ! allComputeOnGPU
                  if (wantDebug) call obj%timer%start("memcpy")
                  ! copy top_border_recv_buffer to aIntern_dev, maybe not necessary if CUDA_AWARE IRECV
                  dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
                  successGPU =  gpu_memcpy(aIntern_dev + dev_offset ,int(loc( top_border_recv_buffer(:,i)),kind=c_intptr_t),  &
                                        stripe_width * top_msg_length * size_of_datatype,   &
                                        gpuMemcpyHostToDevice)
                  call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 2748,  successGPU)
                  if (wantDebug) call obj%timer%stop("memcpy")
                endif ! allComputeOnGPU
              else ! useGPU
                aIntern(:,a_off+1:a_off+top_msg_length,i) = &
                reshape(top_border_recv_buffer(1:stripe_width*top_msg_length,i),(/stripe_width,top_msg_length/))
              endif ! useGPU
           endif

           !compute

           if (wantDebug) call obj%timer%start("compute_hh_trafo")
           call compute_hh_trafo_&
             &complex&
             &_&
             &single&
             (obj, useGPU, wantDebug, aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count, max_threads, &
             a_off, nbw, max_blk_size,  bcast_buffer, bcast_buffer_dev, &
             hh_tau_dev, kernel_flops, kernel_time, n_times, 0, top_msg_length, i, &
             last_stripe_width, kernel)
           if (wantDebug) call obj%timer%stop("compute_hh_trafo")

         endif

         if (next_top_msg_length > 0) then
           !request top_border data

           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX8, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else !useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Irecv(top_border_recv_buffer(1,i), int(next_top_msg_length*stripe_width,kind=MPI_KIND), &
                         MPI_COMPLEX8, int(my_prow-1,kind=MPI_KIND), int(top_recv_tag,kind=MPI_KIND), &
                         int(mpi_comm_rows,kind=MPI_KIND), top_recv_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif !useGPU


         endif

         !send_t
         if (my_prow > 0) then

           if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
           call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
           if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")

           if (useGPU) then
             if (allComputeOnGPU) then
               if (wantDebug) call obj%timer%start("cuda_memcpy")
               successGPU =  gpu_memcpy(c_loc(top_border_send_buffer_mpi_fortran_ptr(1,i)),  &
                                        c_loc(aIntern_mpi_fortran_ptr(1,a_off+1,i)), &
                                        stripe_width* nbw* size_of_datatype,      &
                                        gpuMemcpyDeviceToDevice)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer -> aIntern_dev", 3030,  successGPU)
               if (wantDebug) call obj%timer%stop("cuda_memcpy")
             else ! allComputeOnGPU
               if (wantDebug) call obj%timer%start("memcpy")
               dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_datatype
               successGPU =  gpu_memcpy(int(loc(top_border_send_buffer(:,i)),kind=c_intptr_t), aIntern_dev + dev_offset, &
                                         stripe_width*nbw * size_of_datatype, &
                                         gpuMemcpyDeviceToHost)
               call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> top_border_send_buffer", 3038,  successGPU)
               if (wantDebug) call obj%timer%stop("memcpy")
             endif ! allComputeOnGPU
           else ! useGPU
             top_border_send_buffer(:,i) = reshape(aIntern(:,a_off+1:a_off+nbw,i),(/stripe_width*nbw/))
           endif ! useGPU
           if (useGPU) then
             if (wantDebug) call obj%timer%start("host_mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_COMPLEX8, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_communication")
           else ! useGPU
             if (wantDebug) call obj%timer%start("mpi_communication")
             call MPI_Isend(top_border_send_buffer(1,i), int(nbw*stripe_width,kind=MPI_KIND), MPI_COMPLEX8, &
                            int(my_prow-1,kind=MPI_KIND), int(bottom_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND),   &
                            top_send_request(i), mpierr)
             if (wantDebug) call obj%timer%stop("mpi_communication")
           endif ! useGPU

         endif

         ! Care that there are not too many outstanding top_recv_request's
         if (stripe_count > 1) then
           if (i>1) then

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
             call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")
           else

             if (wantDebug) call obj%timer%start("host_mpi_wait_top_recv")
            call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)
             if (wantDebug) call obj%timer%stop("host_mpi_wait_top_recv")

          endif
        endif
      enddo ! i = 1, stripe_count

      top_msg_length = next_top_msg_length

    else
      ! wait for last top_send_request

      do i = 1, stripe_count
        if (wantDebug) call obj%timer%start("host_mpi_wait_top_send")
        call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_top_send")
      enddo
    endif

    ! Care about the result

    if (my_prow == 0) then

      ! topmost process sends nbw rows to destination processes

      do j=0, nfact-1
        num_blk = sweep*nfact+j ! global number of destination block, 0 based
        if (num_blk*nblk >= na) exit

        nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block

        if (wantDebug) call obj%timer%start("host_mpi_wait_result_send")
        call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)
        if (wantDebug) call obj%timer%stop("host_mpi_wait_result_send")
        dst = mod(num_blk, np_rows)

        if (dst == 0) then
          if (useGPU) then
            row_group_size = min(na - num_blk*nblk, nblk)

            call pack_row_group_&
                 &complex&
                 &_gpu_&
                 &single&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, last_stripe_width, a_dim2, l_nev, &
                         row_group(:, :), j * nblk + a_off, row_group_size, &
                         result_buffer_dev, nblk, num_result_buffers, nbuf, .false., wantDebug, allComputeOnGPU)

            if (allComputeOnGPU) then
              ! memcpy DeviceToDevice row_group_dev -> q_dev
              do i = 1, row_group_size  
               if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
               call gpublas_CCOPY(l_nev, c_loc(row_group_mpi_fortran_ptr(1,i)), 1, &
                                           c_loc(q_mpi_fortran_ptr((num_blk / np_rows) * nblk + i,1)), ldq)
               if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
              enddo
            else ! allComputeOnGPU
              do i = 1, row_group_size
                q((num_blk / np_rows) * nblk + i, 1 : l_nev) = row_group(:, i)
              enddo
            endif ! allComputeOnGPU
          else ! useGPU

            do i = 1, min(na - num_blk*nblk, nblk)

              call pack_row_&
                   &complex&
                   &_cpu_&
                   &single&
                   &(obj,aIntern, row, j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
              q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
            enddo
          endif ! useGPU

        else ! (dst == 0)

          if (useGPU) then

            call pack_row_group_&
                 &complex&
                 &_gpu_&
                 &single&
                 &(obj, row_group_dev, aIntern_dev, stripe_count, stripe_width, &
                   last_stripe_width, a_dim2, l_nev, &
                   result_buffer(:, :, nbuf), j * nblk + a_off, nblk, &
                   result_buffer_dev, nblk, num_result_buffers, nbuf, .true., wantDebug, allComputeOnGPU)

          else  ! useGPU
            do i = 1, nblk
              call pack_row_&
                   &complex&
                   &_cpu_&
                   &single&
                   &(obj, aIntern, result_buffer(:,i,nbuf),j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)
            enddo
          endif ! useGPU
         
          if (useGPU) then
            if (wantDebug) call obj%timer%start("host_mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX8, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("host_mpi_communication")
          else ! useGPU
            if (wantDebug) call obj%timer%start("mpi_communication")
            call MPI_Isend(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX8, &
                           int(dst,kind=MPI_KIND), int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_send_request(nbuf), mpierr)
            if (wantDebug) call obj%timer%stop("mpi_communication")
          endif ! useGPU
        endif ! (dst == 0)
      enddo  !j=0, nfact-1

    else ! (my_prow == 0)

      ! receive and store final result

      do j = num_bufs_recvd, num_result_blocks-1

        nbuf = mod(j, num_result_buffers) + 1 ! buffer number to get this block

        ! If there is still work to do, just test for the next result request
        ! and leave the loop if it is not ready, otherwise wait for all
        ! outstanding requests

        if (next_local_n > 0) then
            ! needed
            !successGPU = gpu_devicesynchronize()
          if (wantDebug) call obj%timer%start("host_mpi_test")
          call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_test")

          if (.not.flag) exit

        else ! (next_local_n > 0)
          if (wantDebug) call obj%timer%start("host_mpi_wait_result_recv")

          call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)
          if (wantDebug) call obj%timer%stop("host_mpi_wait_result_recv")
        endif ! (next_local_n > 0)

        ! Fill result buffer into q
        num_blk = j*np_rows + my_prow ! global number of current block, 0 based
        if (useGPU) then
          if (allComputeOnGPU) then
            do i = 1, min(na - num_blk*nblk, nblk)
              if (wantDebug) call obj%timer%start("cuda_aware_gpublas")
              call gpublas_CCOPY(l_nev, c_loc(result_buffer_mpi_fortran_ptr(1,i,nbuf)), 1, &
                                          c_loc(q_mpi_fortran_ptr(j*nblk + i,1)), ldq)
              if (wantDebug) call obj%timer%stop("cuda_aware_gpublas")
            enddo
          else ! allComputeOnGPU
            do i = 1, min(na - num_blk*nblk, nblk)
              q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
            enddo
          endif ! allComputeOnGPU
        else ! useGPU
          do i = 1, min(na - num_blk*nblk, nblk)
            q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
          enddo
        endif ! useGPU
        ! Queue result buffer again if there are outstanding blocks left

        if (useGPU) then
          if (wantDebug) call obj%timer%start("host_mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX8, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("host_mpi_communication")
        else ! useGPU
          if (wantDebug) call obj%timer%start("mpi_communication")

          if (j+num_result_buffers < num_result_blocks) &
          call MPI_Irecv(result_buffer(1,1,nbuf), int(l_nev*nblk,kind=MPI_KIND), MPI_COMPLEX8, &
                           0_MPI_KIND, int(result_recv_tag,kind=MPI_KIND), int(mpi_comm_rows,kind=MPI_KIND), &
                           result_recv_request(nbuf), mpierr)

          ! carefull the "recieve" has to be done at the corresponding wait or send
          !         if (j+num_result_buffers < num_result_blocks) &
          !                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)
          if (wantDebug) call obj%timer%stop("mpi_communication")
        endif !useGPU

      enddo ! j = num_bufs_recvd, num_result_blocks-1
      num_bufs_recvd = j

    endif ! (my_prow == 0)

    ! Shift the remaining rows to the front of aIntern (if necessary)

    offset = nbw - top_msg_length
    if (offset<0) then
      if (wantDebug) write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_&
                                         &complex&
                                         &: internal error, offset for shifting = ',offset
      success = .false.
      return
    endif

    a_off = a_off + offset
    if (a_off + next_local_n + nbw >= a_dim2) then
      do i = 1, stripe_count
        if (useGPU) then
          chunk = min(next_local_n,a_off)

          if (chunk < 1) exit

          if (wantDebug) call obj%timer%start("normal_memcpy")
          do j = top_msg_length+1, top_msg_length+next_local_n, chunk
            this_chunk = min(j+chunk-1,top_msg_length+next_local_n)-j+1
            dev_offset = ((j-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            dev_offset_1 = ((j+a_off-1)*stripe_width+(i-1)*stripe_width*a_dim2)*size_of_datatype
            num = stripe_width*this_chunk*size_of_datatype
            successGPU = gpu_memcpy(aIntern_dev+dev_offset, aIntern_dev+dev_offset_1, num, gpuMemcpyDeviceToDevice)

            call check_memcpy_GPU_f("trans_ev_tridi_to_band: aIntern_dev -> aIntern_dev", 3528,  successGPU)
          end do
          if (wantDebug) call obj%timer%stop("normal_memcpy")
        else ! not useGPU
          do j = top_msg_length+1, top_msg_length+next_local_n
            aIntern(:,j,i) = aIntern(:,j+a_off,i)
          end do
        end if
      end do ! stripe_count

      a_off = 0
    end if
  end do ! sweep
  if (wantDebug) call obj%timer%stop("sweep_loop")

  ! Just for safety:
  if (ANY(top_send_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_send_request ***',my_prow,my_pcol
  if (ANY(bottom_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_send_request ***',my_prow,my_pcol
  if (ANY(top_recv_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_recv_request ***',my_prow,my_pcol
  if (ANY(bottom_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_recv_request ***',my_prow,my_pcol

  if (my_prow == 0) then

    if (wantDebug) call obj%timer%start("host_mpi_waitall")
    call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)
    if (wantDebug) call obj%timer%stop("host_mpi_waitall")
  endif

  if (ANY(result_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_send_request ***',my_prow,my_pcol
  if (ANY(result_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_recv_request ***',my_prow,my_pcol

  call obj%get("print_flops",print_flops,error)

  if (print_flops == 1) then

    if (useNonBlockingCollectivesRows) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          allreduce_request1, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          allreduce_request3, mpierr)

      call mpi_wait(allreduce_request1, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request3, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_ROWS, &
                          mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_ROWS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif

    if (useNonBlockingCollectivesCols) then
      call mpi_iallreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                          allreduce_request2, mpierr)

      call mpi_iallreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          allreduce_request4, mpierr)

      call mpi_wait(allreduce_request2, MPI_STATUS_IGNORE, mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_wait(allreduce_request4, MPI_STATUS_IGNORE, mpierr)
      kernel_time_recv = kernel_time
    else
      call mpi_allreduce(kernel_flops, kernel_flops_recv, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_COLS, &
                         mpierr)
      kernel_flops = kernel_flops_recv

      call mpi_allreduce(kernel_time, kernel_time_recv, 1, MPI_REAL8, MPI_MAX, MPI_COMM_COLS, &
                          mpierr)
      kernel_time_recv = kernel_time
    endif
  endif

  if (useGPU .and. allComputeOnGPU) then
    ! finally copy q_dev to q
    if (wantDebug) call obj%timer%start("cuda_memcpy")
    successGPU =  gpu_memcpy(int(loc(q(1,1)),kind=c_intptr_t),  &
                             q_dev, &
                             ldq*matrixCols * size_of_datatype, &
                             gpuMemcpyDeviceToHost)
    call check_memcpy_GPU_f("trans_ev_tridi_to_band 1: q_dev -> q", 3634,  successGPU)
    if (wantDebug) call obj%timer%stop("cuda_memcpy")

  endif

  if (my_prow==0 .and. my_pcol==0 .and.print_flops == 1) &
      write(error_unit,'(" Kernel time:",f10.3," MFlops: ",es12.5)')  kernel_time, kernel_flops/kernel_time*1.d-6

  ! deallocate all working space

  if (.not.(useGPU)) then
    nullify(aIntern)
    call free(aIntern_ptr)
  endif

  deallocate(row, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: row", 3650,  istat,  errorMessage)

  deallocate(limits, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: limits", 3653,  istat,  errorMessage)

  deallocate(result_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_send_request", 3656,  istat,  errorMessage)

  deallocate(result_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_recv_request", 3659,  istat,  errorMessage)

  deallocate(result_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: result_buffer", 3662,  istat,  errorMessage)

  if (useGPU) then
    if (allComputeOnGPU) then
      successGPU = gpu_free(result_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: result_buffer_dev", 3667,  successGPU)
      nullify(result_buffer_mpi_fortran_ptr)
    endif

    nullify(bcast_buffer)

    successGPU = gpu_free_host(bcast_buffer_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_host", 3674,  successGPU)
  else
    deallocate(bcast_buffer, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("trans_ev_tridi_to_band: bcast_buffer", 3677,  istat,  errorMessage)
  endif


  if (useGPU) then
    successGPU = gpu_free(aIntern_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: aIntern_dev", 3683,  successGPU)

    if (allComputeOnGPU) then
      successGPU = gpu_free(q_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: q_dev", 3687,  successGPU)
      nullify(q_mpi_fortran_ptr)

      successGPU = gpu_free(hh_trans_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_trans_dev", 3691,  successGPU)
      nullify(hh_trans_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer_dev", 3695,  successGPU)
      nullify(top_border_recv_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(top_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer_dev", 3699,  successGPU)
      nullify(top_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_send_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer_dev", 3703,  successGPU)
      nullify(bottom_border_send_buffer_mpi_fortran_ptr)

      successGPU = gpu_free(bottom_border_recv_buffer_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer_dev", 3707,  successGPU)
      nullify(bottom_border_recv_buffer_mpi_fortran_ptr)

      nullify(aIntern_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU = gpu_free(hh_tau_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: hh_tau_dev", 3714,  successGPU)

    nullify(row_group)

    successGPU = gpu_free_host(row_group_host)
    call check_host_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_host", 3719,  successGPU)

    successGPU = gpu_free(row_group_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_group_dev", 3722,  successGPU)

    if (allComputeOnGPU) then
      nullify(row_group_mpi_fortran_ptr)

      successGPU = gpu_free(row_dev)
      call check_dealloc_GPU_f("trans_ev_tridi_to_band: row_dev", 3728,  successGPU)
      nullify(row_mpi_fortran_ptr)
    endif ! allComputeOnGPU

    successGPU =  gpu_free(bcast_buffer_dev)
    call check_dealloc_GPU_f("trans_ev_tridi_to_band: bcast_buffer_dev", 3733,  successGPU)

    if (allComputeOnGPU) then
      nullify(bcast_buffer_mpi_fortran_ptr)
    endif

    successGPU = gpu_host_unregister(int(loc(top_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_send_buffer", 3740,  successGPU)

    successGPU = gpu_host_unregister(int(loc(top_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3743,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_send_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3746,  successGPU)

    successGPU = gpu_host_unregister(int(loc(bottom_border_recv_buffer),kind=c_intptr_t))
    call check_host_unregister_GPU_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3749,  successGPU)
  endif ! useGPU

  deallocate(top_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_send_buffer", 3753,  istat,  errorMessage)

  deallocate(top_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_border_recv_buffer", 3756,  istat,  errorMessage)

  deallocate(bottom_border_send_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_send_buffer", 3759,  istat,  errorMessage)

  deallocate(bottom_border_recv_buffer, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_border_recv_buffer", 3762,  istat,  errorMessage)

  deallocate(top_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_send_request", 3765,  istat,  errorMessage)

  deallocate(top_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: top_recv_request", 3768,  istat,  errorMessage)

  deallocate(bottom_send_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_send_request", 3771,  istat,  errorMessage)

  deallocate(bottom_recv_request, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("trans_ev_tridi_to_band: bottom_recv_request", 3774,  istat,  errorMessage)

  call obj%timer%stop("trans_ev_tridi_to_band_&
                      &complex&
                      &" // &
                      &"_single" //&
                      gpuString)
  return

end subroutine

! vim: syntax=fortran



end module ELPA2_compute
