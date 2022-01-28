










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
  use elpa_cholesky
  use elpa_invert_trm

  implicit none

  public :: elpa_mult_at_b_real_double_impl      !< Multiply double-precision real matrices A**T * B

  public :: elpa_mult_ah_b_complex_double_impl   !< Multiply double-precision complex matrices A**H * B

  public :: elpa_solve_tridi_double_impl         !< Solve tridiagonal eigensystem for a double-precision matrix with divide and conquer method

  public :: elpa_mult_at_b_real_single_impl      !< Multiply single-precision real matrices A**T * B
  public :: elpa_solve_tridi_single_impl         !< Solve tridiagonal eigensystem for a single-precision matrix with divide and conquer method

  public :: elpa_mult_ah_b_complex_single_impl   !< Multiply single-precision complex matrices A**H * B

  contains

















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




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


  use elpa1_compute
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding
  use elpa_gpu
  use mod_check_for_gpu
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
  real(kind=rck), pointer, contiguous :: aux_mat(:,:), tmp1(:,:)
  real(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
  integer(kind=ik)              :: istat
  character(200)                :: errorMessage
  character(20)                 :: gpuString
  logical                       :: success
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
  integer(kind=c_intptr_t)      :: aux_dev, b_dev, tmp1_dev
  type(c_ptr)                   :: aux_host, tmp1_host
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: aux_off, b_off
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &real

  success = .true.

  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
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
  l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
  l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

  ! Block factor for matrix multiplications, must be a multiple of nblk

  if (na/np_rows <= 256) then
     nblk_mult = (31/nblk+1)*nblk
  else
     nblk_mult = (63/nblk+1)*nblk
  endif

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")

    ! copy b to b_dev
    num = ldb*ldbCols*size_of_datatype
    successGPU = gpu_malloc(b_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: b_dev", 218,  successGPU)

    successGPU = gpu_host_register(int(loc(b),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)

    call check_host_register_GPU_f("elpa_mult_at_b: b", 223,  successGPU)

    successGPU = gpu_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
                  gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_mult_at_b: b to b_dev", 227,  successGPU)

    num = l_rows*nblk_mult*size_of_datatype
    successGPU = gpu_malloc_host(aux_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: aux_host", 231,  successGPU)

    call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

    successGPU = gpu_malloc(aux_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: aux_dev", 236,  successGPU)

    num = nblk_mult*l_cols*size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: tmp1_host", 240,  successGPU)

    call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

    successGPU = gpu_malloc(tmp1_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: tmp1_dev", 245,  successGPU)
  else ! useGPU
    allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_mult_at_b: aux_mat", 248,  istat,  errorMessage)
  endif ! useGPU

  allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: aux_bc", 252,  istat,  errorMessage)

  allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lrs_save", 255,  istat,  errorMessage)

  allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lre_save", 258,  istat,  errorMessage)

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
              successGPU = gpu_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                            num, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("elpa_mult_at_b: aux_mat to aux_dev", 373,  successGPU)

              aux_off = (lrs-1)*size_of_datatype
              b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

              call obj%timer%start("gpublas")
              call gpublas_DGEMM('T', 'N', nstor, lce-lcs+1, &
                   lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                   tmp1_dev, nstor)
              call obj%timer%stop("gpublas")

              num = nstor*(lce-lcs+1)*size_of_datatype
              successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                            tmp1_dev, num, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("elpa_mult_at_b: tmp1_dev to tmp1", 387,  successGPU)
            else ! useGPU
              call obj%timer%start("blas")
              call DGEMM('T', 'N', int(nstor,kind=BLAS_KIND), &
                                int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                                ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                                b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                                int(nstor,kind=BLAS_KIND))
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
          call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 418,  istat,  errorMessage)
        endif

        nr_done = nr_done+nstor
        nstor=0
        aux_mat(:,:)=0
      endif
    enddo
  enddo

  if (useGPU) then
    successGPU = gpu_free(b_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: b_dev", 430,  successGPU)

    successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 433,  successGPU)

    nullify(aux_mat)
    nullify(tmp1)

    successGPU = gpu_free_host(aux_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: aux_host", 439,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: aux_dev", 442,  successGPU)

    successGPU = gpu_free_host(tmp1_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: tmp1_host", 445,  successGPU)

    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: tmp1_dev", 448,  successGPU)
  else ! useGPU
    deallocate(aux_mat, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_mult_at_b: aux_mat", 451,  istat,  errorMessage)
  endif ! useGPU

  deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 455,  istat,  errorMessage)


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




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


  use elpa1_compute
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding
  use elpa_gpu
  use mod_check_for_gpu
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
  real(kind=rck), pointer, contiguous :: aux_mat(:,:), tmp1(:,:)
  real(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
  integer(kind=ik)              :: istat
  character(200)                :: errorMessage
  character(20)                 :: gpuString
  logical                       :: success
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
  integer(kind=c_intptr_t)      :: aux_dev, b_dev, tmp1_dev
  type(c_ptr)                   :: aux_host, tmp1_host
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: aux_off, b_off
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &real

  success = .true.

  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
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
  l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
  l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

  ! Block factor for matrix multiplications, must be a multiple of nblk

  if (na/np_rows <= 256) then
     nblk_mult = (31/nblk+1)*nblk
  else
     nblk_mult = (63/nblk+1)*nblk
  endif

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")

    ! copy b to b_dev
    num = ldb*ldbCols*size_of_datatype
    successGPU = gpu_malloc(b_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: b_dev", 218,  successGPU)

    successGPU = gpu_host_register(int(loc(b),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)

    call check_host_register_GPU_f("elpa_mult_at_b: b", 223,  successGPU)

    successGPU = gpu_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
                  gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_mult_at_b: b to b_dev", 227,  successGPU)

    num = l_rows*nblk_mult*size_of_datatype
    successGPU = gpu_malloc_host(aux_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: aux_host", 231,  successGPU)

    call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

    successGPU = gpu_malloc(aux_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: aux_dev", 236,  successGPU)

    num = nblk_mult*l_cols*size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: tmp1_host", 240,  successGPU)

    call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

    successGPU = gpu_malloc(tmp1_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: tmp1_dev", 245,  successGPU)
  else ! useGPU
    allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_mult_at_b: aux_mat", 248,  istat,  errorMessage)
  endif ! useGPU

  allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: aux_bc", 252,  istat,  errorMessage)

  allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lrs_save", 255,  istat,  errorMessage)

  allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lre_save", 258,  istat,  errorMessage)

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
              successGPU = gpu_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                            num, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("elpa_mult_at_b: aux_mat to aux_dev", 373,  successGPU)

              aux_off = (lrs-1)*size_of_datatype
              b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

              call obj%timer%start("gpublas")
              call gpublas_SGEMM('T', 'N', nstor, lce-lcs+1, &
                   lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                   tmp1_dev, nstor)
              call obj%timer%stop("gpublas")

              num = nstor*(lce-lcs+1)*size_of_datatype
              successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                            tmp1_dev, num, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("elpa_mult_at_b: tmp1_dev to tmp1", 387,  successGPU)
            else ! useGPU
              call obj%timer%start("blas")
              call SGEMM('T', 'N', int(nstor,kind=BLAS_KIND), &
                                int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                                ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                                b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                                int(nstor,kind=BLAS_KIND))
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
          call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 418,  istat,  errorMessage)
        endif

        nr_done = nr_done+nstor
        nstor=0
        aux_mat(:,:)=0
      endif
    enddo
  enddo

  if (useGPU) then
    successGPU = gpu_free(b_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: b_dev", 430,  successGPU)

    successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 433,  successGPU)

    nullify(aux_mat)
    nullify(tmp1)

    successGPU = gpu_free_host(aux_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: aux_host", 439,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: aux_dev", 442,  successGPU)

    successGPU = gpu_free_host(tmp1_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: tmp1_host", 445,  successGPU)

    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: tmp1_dev", 448,  successGPU)
  else ! useGPU
    deallocate(aux_mat, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_mult_at_b: aux_mat", 451,  istat,  errorMessage)
  endif ! useGPU

  deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 455,  istat,  errorMessage)


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




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


  use elpa1_compute
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding
  use elpa_gpu
  use mod_check_for_gpu
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
  complex(kind=rck), pointer, contiguous :: aux_mat(:,:), tmp1(:,:)
  complex(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
  integer(kind=ik)              :: istat
  character(200)                :: errorMessage
  character(20)                 :: gpuString
  logical                       :: success
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
  integer(kind=c_intptr_t)      :: aux_dev, b_dev, tmp1_dev
  type(c_ptr)                   :: aux_host, tmp1_host
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: aux_off, b_off
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &double&
                                                            &_&
                                                            &complex

  success = .true.

  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
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
  l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
  l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

  ! Block factor for matrix multiplications, must be a multiple of nblk

  if (na/np_rows <= 256) then
     nblk_mult = (31/nblk+1)*nblk
  else
     nblk_mult = (63/nblk+1)*nblk
  endif

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")

    ! copy b to b_dev
    num = ldb*ldbCols*size_of_datatype
    successGPU = gpu_malloc(b_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: b_dev", 218,  successGPU)

    successGPU = gpu_host_register(int(loc(b),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)

    call check_host_register_GPU_f("elpa_mult_at_b: b", 223,  successGPU)

    successGPU = gpu_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
                  gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_mult_at_b: b to b_dev", 227,  successGPU)

    num = l_rows*nblk_mult*size_of_datatype
    successGPU = gpu_malloc_host(aux_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: aux_host", 231,  successGPU)

    call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

    successGPU = gpu_malloc(aux_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: aux_dev", 236,  successGPU)

    num = nblk_mult*l_cols*size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: tmp1_host", 240,  successGPU)

    call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

    successGPU = gpu_malloc(tmp1_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: tmp1_dev", 245,  successGPU)
  else ! useGPU
    allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_mult_at_b: aux_mat", 248,  istat,  errorMessage)
  endif ! useGPU

  allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: aux_bc", 252,  istat,  errorMessage)

  allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lrs_save", 255,  istat,  errorMessage)

  allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lre_save", 258,  istat,  errorMessage)

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
              successGPU = gpu_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                            num, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("elpa_mult_at_b: aux_mat to aux_dev", 373,  successGPU)

              aux_off = (lrs-1)*size_of_datatype
              b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

              call obj%timer%start("gpublas")
              call gpublas_ZGEMM('C', 'N', nstor, lce-lcs+1, &
                   lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                   tmp1_dev, nstor)
              call obj%timer%stop("gpublas")

              num = nstor*(lce-lcs+1)*size_of_datatype
              successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                            tmp1_dev, num, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("elpa_mult_at_b: tmp1_dev to tmp1", 387,  successGPU)
            else ! useGPU
              call obj%timer%start("blas")
              call ZGEMM('C', 'N', int(nstor,kind=BLAS_KIND), &
                                int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                                ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                                b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                                int(nstor,kind=BLAS_KIND))
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
          call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 418,  istat,  errorMessage)
        endif

        nr_done = nr_done+nstor
        nstor=0
        aux_mat(:,:)=0
      endif
    enddo
  enddo

  if (useGPU) then
    successGPU = gpu_free(b_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: b_dev", 430,  successGPU)

    successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 433,  successGPU)

    nullify(aux_mat)
    nullify(tmp1)

    successGPU = gpu_free_host(aux_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: aux_host", 439,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: aux_dev", 442,  successGPU)

    successGPU = gpu_free_host(tmp1_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: tmp1_host", 445,  successGPU)

    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: tmp1_dev", 448,  successGPU)
  else ! useGPU
    deallocate(aux_mat, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_mult_at_b: aux_mat", 451,  istat,  errorMessage)
  endif ! useGPU

  deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 455,  istat,  errorMessage)


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




!cannot use "../src/elpa1/../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


  use elpa1_compute
  use elpa_mpi
  use precision
  use elpa_abstract_impl
  use, intrinsic :: iso_c_binding
  use elpa_gpu
  use mod_check_for_gpu
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
  complex(kind=rck), pointer, contiguous :: aux_mat(:,:), tmp1(:,:)
  complex(kind=rck), allocatable :: aux_bc(:), tmp2(:,:)
  integer(kind=ik)              :: istat
  character(200)                :: errorMessage
  character(20)                 :: gpuString
  logical                       :: success
  logical                       :: successGPU
  logical                       :: useGPU
  integer(kind=c_int)           :: gpu, numGPU
  integer(kind=ik)              :: mpi_comm_rows, mpi_comm_cols, mpi_comm_all
  integer(kind=ik)              :: nblk, matrixRows, matrixCols, error
  integer(kind=c_intptr_t)      :: aux_dev, b_dev, tmp1_dev
  type(c_ptr)                   :: aux_host, tmp1_host
  integer(kind=c_intptr_t)      :: num
  integer(kind=c_intptr_t)      :: aux_off, b_off
  integer(kind=c_intptr_t), parameter :: size_of_datatype = size_of_&
                                                            &single&
                                                            &_&
                                                            &complex

  success = .true.

  ! GPU settings
  if (gpu_vendor() == NVIDIA_GPU) then
    call obj%get("gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for GPU. Aborting..."
      stop
    endif
    if (gpu .eq. 1) then
      print *,"You still use the deprecated option 'gpu', consider switching to 'nvidia-gpu'. Will set the new &
              & keyword 'nvidia-gpu'"
      call obj%set("nvidia-gpu",gpu,error)
      if (error .ne. ELPA_OK) then
        print *,"Problem setting option for NVIDIA GPU. Aborting..."
        stop
      endif
    endif

    call obj%get("nvidia-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for NVIDIA GPU. Aborting..."
      stop
    endif
  else if (gpu_vendor() == AMD_GPU) then
    call obj%get("amd-gpu",gpu,error)
    if (error .ne. ELPA_OK) then
      print *,"Problem getting option for AMD GPU. Aborting..."
      stop
    endif
  else
    gpu = 0
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
  l_rows = local_index(na,  my_prow, np_rows, nblk, -1) ! Local rows of a and b
  l_cols = local_index(ncb, my_pcol, np_cols, nblk, -1) ! Local cols of b

  ! Block factor for matrix multiplications, must be a multiple of nblk

  if (na/np_rows <= 256) then
     nblk_mult = (31/nblk+1)*nblk
  else
     nblk_mult = (63/nblk+1)*nblk
  endif

  if (useGPU) then
    call obj%timer%start("check_for_gpu")
    if (check_for_gpu(obj, myid, numGPU)) then
      ! set the neccessary parameters
      call set_gpu_parameters()
    else
      print *,"GPUs are requested but not detected! Aborting..."
      success = .false.
      return
    endif
    call obj%timer%stop("check_for_gpu")

    ! copy b to b_dev
    num = ldb*ldbCols*size_of_datatype
    successGPU = gpu_malloc(b_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: b_dev", 218,  successGPU)

    successGPU = gpu_host_register(int(loc(b),kind=c_intptr_t),num,&
                  gpuHostRegisterDefault)

    call check_host_register_GPU_f("elpa_mult_at_b: b", 223,  successGPU)

    successGPU = gpu_memcpy(b_dev,int(loc(b),kind=c_intptr_t),num,&
                  gpuMemcpyHostToDevice)
    call check_memcpy_GPU_f("elpa_mult_at_b: b to b_dev", 227,  successGPU)

    num = l_rows*nblk_mult*size_of_datatype
    successGPU = gpu_malloc_host(aux_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: aux_host", 231,  successGPU)

    call c_f_pointer(aux_host,aux_mat,(/l_rows,nblk_mult/))

    successGPU = gpu_malloc(aux_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: aux_dev", 236,  successGPU)

    num = nblk_mult*l_cols*size_of_datatype
    successGPU = gpu_malloc_host(tmp1_host,num)
    call check_host_alloc_GPU_f("elpa_mult_at_b: tmp1_host", 240,  successGPU)

    call c_f_pointer(tmp1_host,tmp1,(/nblk_mult,l_cols/))

    successGPU = gpu_malloc(tmp1_dev,num)
    call check_alloc_GPU_f("elpa_mult_at_b: tmp1_dev", 245,  successGPU)
  else ! useGPU
    allocate(aux_mat(l_rows,nblk_mult), stat=istat, errmsg=errorMessage)
    call check_allocate_f("elpa_mult_at_b: aux_mat", 248,  istat,  errorMessage)
  endif ! useGPU

  allocate(aux_bc(l_rows*nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: aux_bc", 252,  istat,  errorMessage)

  allocate(lrs_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lrs_save", 255,  istat,  errorMessage)

  allocate(lre_save(nblk), stat=istat, errmsg=errorMessage)
  call check_allocate_f("elpa_mult_at_b: lre_save", 258,  istat,  errorMessage)

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
              successGPU = gpu_memcpy(aux_dev, int(loc(aux_mat),kind=c_intptr_t), &
                            num, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("elpa_mult_at_b: aux_mat to aux_dev", 373,  successGPU)

              aux_off = (lrs-1)*size_of_datatype
              b_off = ((lcs-1)*ldb+lrs-1)*size_of_datatype

              call obj%timer%start("gpublas")
              call gpublas_CGEMM('C', 'N', nstor, lce-lcs+1, &
                   lre-lrs+1, ONE, aux_dev+aux_off, l_rows, b_dev+b_off, ldb, ZERO, &
                   tmp1_dev, nstor)
              call obj%timer%stop("gpublas")

              num = nstor*(lce-lcs+1)*size_of_datatype
              successGPU = gpu_memcpy(int(loc(tmp1),kind=c_intptr_t), &
                            tmp1_dev, num, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("elpa_mult_at_b: tmp1_dev to tmp1", 387,  successGPU)
            else ! useGPU
              call obj%timer%start("blas")
              call CGEMM('C', 'N', int(nstor,kind=BLAS_KIND), &
                                int(lce-lcs+1,kind=BLAS_KIND), int(lre-lrs+1,kind=BLAS_KIND), &
                                ONE, aux_mat(lrs:lre,1:nstor), int(lre-lrs+1,kind=BLAS_KIND), &
                                b(lrs,lcs), int(ldb,kind=BLAS_KIND), ZERO, tmp1, &
                                int(nstor,kind=BLAS_KIND))
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
          call check_deallocate_f("elpa_mult_at_b: tmp1, tmp2", 418,  istat,  errorMessage)
        endif

        nr_done = nr_done+nstor
        nstor=0
        aux_mat(:,:)=0
      endif
    enddo
  enddo

  if (useGPU) then
    successGPU = gpu_free(b_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: b_dev", 430,  successGPU)

    successGPU = gpu_host_unregister(int(loc(b),kind=c_intptr_t))
    call check_host_unregister_GPU_f("elpa_multiply_a_b: b", 433,  successGPU)

    nullify(aux_mat)
    nullify(tmp1)

    successGPU = gpu_free_host(aux_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: aux_host", 439,  successGPU)

    successGPU = gpu_free(aux_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: aux_dev", 442,  successGPU)

    successGPU = gpu_free_host(tmp1_host)
    call check_host_dealloc_GPU_f("elpa_multiply_a_b: tmp1_host", 445,  successGPU)

    successGPU = gpu_free(tmp1_dev)
    call check_dealloc_GPU_f("elpa_multiply_a_b: tmp1_dev", 448,  successGPU)
  else ! useGPU
    deallocate(aux_mat, stat=istat, errmsg=errorMessage)
    call check_deallocate_f("elpa_mult_at_b: aux_mat", 451,  istat,  errorMessage)
  endif ! useGPU

  deallocate(aux_bc, lrs_save, lre_save, stat=istat, errmsg=errorMessage)
  call check_deallocate_f("elpa_mult_at_b: aux_bc, lrs_save, lre_save", 455,  istat,  errorMessage)


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







      use elpa1_compute, solve_tridi_&
                         &double&
                         &_private_impl => solve_tridi_&
                         &double&
                         &_impl
      use precision
      use elpa_abstract_impl
      use elpa_omp
      use solve_tridi
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)         :: na, nev, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      integer(kind=ik)         :: mpi_comm_all
      real(kind=rk8) :: d(obj%na), e(obj%na)
      real(kind=rk8) :: q(obj%local_nrows,*)

      logical                  :: wantDebug
      logical                  :: success

      integer                  :: debug, error
      integer                  :: nrThreads, limitThreads

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
      call obj%get("mpi_comm_parent", mpi_comm_all,error)
      if (error .ne. ELPA_OK) then
        print *,"Problem getting option for mpi_comm_all. Aborting..."
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
               mpi_comm_all, mpi_comm_rows, mpi_comm_cols,.false., wantDebug, success, &
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







      use elpa1_compute, solve_tridi_&
                         &single&
                         &_private_impl => solve_tridi_&
                         &single&
                         &_impl
      use precision
      use elpa_abstract_impl
      use elpa_omp
      use solve_tridi
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)         :: na, nev, matrixRows, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols
      integer(kind=ik)         :: mpi_comm_all
      real(kind=rk4) :: d(obj%na), e(obj%na)
      real(kind=rk4) :: q(obj%local_nrows,*)

      logical                  :: wantDebug
      logical                  :: success

      integer                  :: debug, error
      integer                  :: nrThreads, limitThreads

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
      call obj%get("mpi_comm_parent", mpi_comm_all,error)
      if (error .ne. ELPA_OK) then
        print *,"Problem getting option for mpi_comm_all. Aborting..."
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
               mpi_comm_all, mpi_comm_rows, mpi_comm_cols,.false., wantDebug, success, &
               nrThreads)


      ! restore original OpenMP settings


      call obj%timer%stop("elpa_solve_tridi_public_&
      &real&
      &_&
      &single&
      &")




    end function




end module elpa1_auxiliary_impl

