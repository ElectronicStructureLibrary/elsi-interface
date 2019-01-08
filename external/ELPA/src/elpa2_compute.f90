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


module ELPA2_compute

! Version 1.1.2, 2011-02-21

  use ELPA_utilities
  USE ELPA1_compute
  use elpa1, only : elpa_print_times, time_evp_back, time_evp_fwd, time_evp_solve
  use elpa2_utilities
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



  public :: bandred_complex_double
  public :: tridiag_band_complex_double
  public :: trans_ev_tridi_to_band_complex_double
  public :: trans_ev_band_to_full_complex_double


  public :: band_band_real_double
  public :: divide_band

  integer(kind=ik), public :: which_qr_decomposition = 1     ! defines, which QR-decomposition algorithm will be used
! 0 for unblocked
! 1 for blocked (maxrank: nblk)
  contains

! real double precision first






! --------------------------------------------------------------------------------------------------
! redist_band: redistributes band from 2D block cyclic form to 1D band




subroutine redist_band_real_double(r_a, lda, na, nblk, nbw, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm, r_ab)









   use precision
   implicit none

   integer(kind=ik), intent(in)     :: lda, na, nblk, nbw, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm

   real(kind=rk8), intent(in)        :: r_a(lda, matrixCols)





   real(kind=rk8), intent(out)       :: r_ab(:,:)




   integer(kind=ik), allocatable    :: ncnt_s(:), nstart_s(:), ncnt_r(:), nstart_r(:), &
                                       global_id(:,:), global_id_tmp(:,:), block_limits(:)

   real(kind=rk8), allocatable       :: r_sbuf(:,:,:), r_rbuf(:,:,:), r_buf(:,:)



   integer(kind=ik)                 :: i, j, my_pe, n_pes, my_prow, np_rows, my_pcol, np_cols, &
                                       nfact, np, npr, npc, mpierr, is, js
   integer(kind=ik)                 :: nblocks_total, il, jl, l_rows, l_cols, n_off



   call mpi_comm_rank(mpi_comm,my_pe,mpierr)
   call mpi_comm_size(mpi_comm,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


! Get global_id mapping 2D procssor coordinates to global id

   allocate(global_id(0:np_rows-1,0:np_cols-1))

   global_id(:,:) = 0
   global_id(my_prow, my_pcol) = my_pe



   call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)



! Set work distribution

   nblocks_total = (na-1)/nbw + 1

   allocate(block_limits(0:n_pes))
   call divide_band(nblocks_total, n_pes, block_limits)


   allocate(ncnt_s(0:n_pes-1))
   allocate(nstart_s(0:n_pes-1))
   allocate(ncnt_r(0:n_pes-1))
   allocate(nstart_r(0:n_pes-1))


   nfact = nbw/nblk

! Count how many blocks go to which PE

   ncnt_s(:) = 0
   np = 0 ! receiver PE number
   do j=0,(na-1)/nblk ! loop over rows of blocks
     if (j/nfact==block_limits(np+1)) np = np+1
     if (mod(j,np_rows) == my_prow) then
       do i=0,nfact
         if (mod(i+j,np_cols) == my_pcol) then
           ncnt_s(np) = ncnt_s(np) + 1
         endif
       enddo
     endif
   enddo

! Allocate send buffer


   allocate(r_sbuf(nblk,nblk,sum(ncnt_s)))
   r_sbuf(:,:,:) = 0.



! Determine start offsets in send buffer

   nstart_s(0) = 0
   do i=1,n_pes-1
     nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   enddo

! Fill send buffer

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of a

   np = 0
   do j=0,(na-1)/nblk ! loop over rows of blocks
     if (j/nfact==block_limits(np+1)) np = np+1
     if (mod(j,np_rows) == my_prow) then
       do i=0,nfact
         if (mod(i+j,np_cols) == my_pcol) then
           nstart_s(np) = nstart_s(np) + 1
           js = (j/np_rows)*nblk
           is = ((i+j)/np_cols)*nblk
           jl = MIN(nblk,l_rows-js)
           il = MIN(nblk,l_cols-is)


           r_sbuf(1:jl,1:il,nstart_s(np)) = r_a(js+1:js+jl,is+1:is+il)


         endif
       enddo
      endif
   enddo

! Count how many blocks we get from which PE

   ncnt_r(:) = 0
   do j=block_limits(my_pe)*nfact,min(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
     npr = mod(j,np_rows)
     do i=0,nfact
       npc = mod(i+j,np_cols)
       np = global_id(npr,npc)
       ncnt_r(np) = ncnt_r(np) + 1
     enddo
   enddo

! Allocate receive buffer


   allocate(r_rbuf(nblk,nblk,sum(ncnt_r)))



! Set send counts/send offsets, receive counts/receive offsets
! now actually in variables, not in blocks

   ncnt_s(:) = ncnt_s(:)*nblk*nblk

   nstart_s(0) = 0
   do i=1,n_pes-1
     nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   enddo

   ncnt_r(:) = ncnt_r(:)*nblk*nblk

   nstart_r(0) = 0
   do i=1,n_pes-1
     nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   enddo

! Exchange all data with MPI_Alltoallv






    call MPI_Alltoallv(r_sbuf, ncnt_s, nstart_s, MPI_REAL8, r_rbuf, ncnt_r, nstart_r, MPI_REAL8, mpi_comm, mpierr)








! set band from receive buffer

   ncnt_r(:) = ncnt_r(:)/(nblk*nblk)

   nstart_r(0) = 0
   do i=1,n_pes-1
     nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   enddo


   allocate(r_buf((nfact+1)*nblk,nblk))



! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nbw

   do j=block_limits(my_pe)*nfact,min(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
     npr = mod(j,np_rows)
     do i=0,nfact
       npc = mod(i+j,np_cols)
       np = global_id(npr,npc)
       nstart_r(np) = nstart_r(np) + 1

       r_buf(i*nblk+1:i*nblk+nblk,:) = transpose(r_rbuf(:,:,nstart_r(np)))


     enddo
     do i=1,MIN(nblk,na-j*nblk)

       r_ab(1:nbw+1,i+j*nblk-n_off) = r_buf(i:i+nbw,i)


     enddo
   enddo

   deallocate(ncnt_s, nstart_s)
   deallocate(ncnt_r, nstart_r)
   deallocate(global_id)
   deallocate(block_limits)


   deallocate(r_sbuf, r_rbuf, r_buf)






end subroutine






! single precision


! double precision






! --------------------------------------------------------------------------------------------------
! redist_band: redistributes band from 2D block cyclic form to 1D band






subroutine redist_band_complex_double(c_a, lda, na, nblk, nbw, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm, c_ab)







   use precision
   implicit none

   integer(kind=ik), intent(in)     :: lda, na, nblk, nbw, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm


   complex(kind=ck8), intent(in)     :: c_a(lda, matrixCols)






   complex(kind=ck8), intent(out)    :: c_ab(:,:)


   integer(kind=ik), allocatable    :: ncnt_s(:), nstart_s(:), ncnt_r(:), nstart_r(:), &
                                       global_id(:,:), global_id_tmp(:,:), block_limits(:)



   complex(kind=ck8), allocatable    :: c_sbuf(:,:,:), c_rbuf(:,:,:), c_buf(:,:)

   integer(kind=ik)                 :: i, j, my_pe, n_pes, my_prow, np_rows, my_pcol, np_cols, &
                                       nfact, np, npr, npc, mpierr, is, js
   integer(kind=ik)                 :: nblocks_total, il, jl, l_rows, l_cols, n_off



   call mpi_comm_rank(mpi_comm,my_pe,mpierr)
   call mpi_comm_size(mpi_comm,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


! Get global_id mapping 2D procssor coordinates to global id

   allocate(global_id(0:np_rows-1,0:np_cols-1))

   global_id(:,:) = 0
   global_id(my_prow, my_pcol) = my_pe



   call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)



! Set work distribution

   nblocks_total = (na-1)/nbw + 1

   allocate(block_limits(0:n_pes))
   call divide_band(nblocks_total, n_pes, block_limits)


   allocate(ncnt_s(0:n_pes-1))
   allocate(nstart_s(0:n_pes-1))
   allocate(ncnt_r(0:n_pes-1))
   allocate(nstart_r(0:n_pes-1))


   nfact = nbw/nblk

! Count how many blocks go to which PE

   ncnt_s(:) = 0
   np = 0 ! receiver PE number
   do j=0,(na-1)/nblk ! loop over rows of blocks
     if (j/nfact==block_limits(np+1)) np = np+1
     if (mod(j,np_rows) == my_prow) then
       do i=0,nfact
         if (mod(i+j,np_cols) == my_pcol) then
           ncnt_s(np) = ncnt_s(np) + 1
         endif
       enddo
     endif
   enddo

! Allocate send buffer



   allocate(c_sbuf(nblk,nblk,sum(ncnt_s)))
   c_sbuf(:,:,:) = 0.


! Determine start offsets in send buffer

   nstart_s(0) = 0
   do i=1,n_pes-1
     nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   enddo

! Fill send buffer

   l_rows = local_index(na, my_prow, np_rows, nblk, -1) ! Local rows of a
   l_cols = local_index(na, my_pcol, np_cols, nblk, -1) ! Local columns of a

   np = 0
   do j=0,(na-1)/nblk ! loop over rows of blocks
     if (j/nfact==block_limits(np+1)) np = np+1
     if (mod(j,np_rows) == my_prow) then
       do i=0,nfact
         if (mod(i+j,np_cols) == my_pcol) then
           nstart_s(np) = nstart_s(np) + 1
           js = (j/np_rows)*nblk
           is = ((i+j)/np_cols)*nblk
           jl = MIN(nblk,l_rows-js)
           il = MIN(nblk,l_cols-is)



           c_sbuf(1:jl,1:il,nstart_s(np)) = c_a(js+1:js+jl,is+1:is+il)

         endif
       enddo
      endif
   enddo

! Count how many blocks we get from which PE

   ncnt_r(:) = 0
   do j=block_limits(my_pe)*nfact,min(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
     npr = mod(j,np_rows)
     do i=0,nfact
       npc = mod(i+j,np_cols)
       np = global_id(npr,npc)
       ncnt_r(np) = ncnt_r(np) + 1
     enddo
   enddo

! Allocate receive buffer



   allocate(c_rbuf(nblk,nblk,sum(ncnt_r)))


! Set send counts/send offsets, receive counts/receive offsets
! now actually in variables, not in blocks

   ncnt_s(:) = ncnt_s(:)*nblk*nblk

   nstart_s(0) = 0
   do i=1,n_pes-1
     nstart_s(i) = nstart_s(i-1) + ncnt_s(i-1)
   enddo

   ncnt_r(:) = ncnt_r(:)*nblk*nblk

   nstart_r(0) = 0
   do i=1,n_pes-1
     nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   enddo

! Exchange all data with MPI_Alltoallv








    call MPI_Alltoallv(c_sbuf, ncnt_s, nstart_s, MPI_COMPLEX16, c_rbuf, ncnt_r, nstart_r, MPI_COMPLEX16, mpi_comm, mpierr)







! set band from receive buffer

   ncnt_r(:) = ncnt_r(:)/(nblk*nblk)

   nstart_r(0) = 0
   do i=1,n_pes-1
     nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   enddo



   allocate(c_buf((nfact+1)*nblk,nblk))


! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nbw

   do j=block_limits(my_pe)*nfact,min(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
     npr = mod(j,np_rows)
     do i=0,nfact
       npc = mod(i+j,np_cols)
       np = global_id(npr,npc)
       nstart_r(np) = nstart_r(np) + 1


       c_buf(i*nblk+1:i*nblk+nblk,:) = conjg(transpose(c_rbuf(:,:,nstart_r(np))))

     enddo
     do i=1,MIN(nblk,na-j*nblk)


       c_ab(1:nbw+1,i+j*nblk-n_off) = c_buf(i:i+nbw,i)

     enddo
   enddo

   deallocate(ncnt_s, nstart_s)
   deallocate(ncnt_r, nstart_r)
   deallocate(global_id)
   deallocate(block_limits)



   deallocate(c_sbuf, c_rbuf, c_buf)





end subroutine










! real double precision










































































    subroutine bandred_real_double(na, a, a_dev, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, &
                            tmat, tmat_dev, wantDebug, useGPU, success, useQR)

!-------------------------------------------------------------------------------
!  bandred_real: Reduces a distributed symmetric matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,matrixCols)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a
!  matrixCols  local columns of matrix a
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

      use cuda_functions
      use, intrinsic :: iso_c_binding
      use elpa1_compute


      use precision
      implicit none

      integer(kind=ik)                      :: na, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols

      real(kind=rk8)              :: a(lda,*), tmat(nbw,nbw,*)

      real(kind=rk8)              :: eps
      logical, intent(in)                   :: useGPU

      integer(kind=ik)                      :: my_prow, my_pcol, np_rows, np_cols, mpierr
      integer(kind=ik)                      :: l_cols, l_rows, vmrCols
      integer(kind=ik)                      :: i, j, lcs, lce, lrs, lre, lc, lr, cur_pcol, n_cols, nrow
      integer(kind=ik)                      :: istep, ncol, lch, lcx, nlc, mynlc
      integer(kind=ik)                      :: tile_size, l_rows_tile, l_cols_tile

      real(kind=rk8)              :: vnorm2, xf, aux1(nbw), aux2(nbw), vrl, tau, vav(nbw,nbw)

      real(kind=rk8), allocatable :: tmpCUDA(:),  vmrCUDA(:),  umcCUDA(:)
      real(kind=rk8), allocatable :: tmpCPU(:,:), vmrCPU(:,:), umcCPU(:,:)
      real(kind=rk8), allocatable :: vr(:)
! needed for blocked QR decomposition
      integer(kind=ik)                      :: PQRPARAM(11), work_size
      real(kind=rk8)              :: dwork_size(1)
      real(kind=rk8), allocatable :: work_blocked(:), tauvector(:), blockheuristic(:)

! a_dev is passed from bandred_real to trans_ev_band
      integer(kind=C_intptr_T)              :: a_dev, vmr_dev, umc_dev, tmat_dev, vav_dev

      integer(kind=ik), external            :: numroc

      integer(kind=ik)                      :: ierr
      integer(kind=ik)                      :: cur_l_rows, cur_l_cols, vmr_size, umc_size
      integer(kind=c_size_t)                :: lc_start, lc_end
      integer(kind=ik)                      :: lr_end
      integer(kind=ik)                      :: na_cols !, na_rows

      logical, intent(in)                   :: wantDebug
      logical, intent(out)                  :: success
      logical                               :: successCUDA
      integer(kind=ik)                      :: istat
      character(200)                        :: errorMessage

      logical, intent(in)                   :: useQR

      integer(kind=ik)                      :: mystart, myend, m_way, n_way, work_per_thread, m_id, n_id, n_threads, &
                                               ii, pp, transformChunkSize



      call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
      call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
      call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
      call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

      success = .true.


      successCUDA = .true.
      umc_dev = 0
      vmr_dev = 0


! Semibandwith nbw must be a multiple of blocksize nblk
      if (mod(nbw,nblk)/=0) then
        if (my_prow==0 .and. my_pcol==0) then
          if (wantDebug) then
            write(error_unit,*) 'ELPA2_bandred_real: ERROR: nbw=',nbw,', nblk=',nblk
            write(error_unit,*) 'ELPA2_bandred_real: ELPA2 works only for nbw==n*nblk'
          endif
          success = .false.
          return
        endif
      endif

! na_rows in used nowhere; only na_cols
      if (useGPU) then

!        na_rows = numroc(na, nblk, my_prow, 0, np_rows)
        na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

      endif ! useGPU

! Matrix is split into tiles; work is done only for tiles on the diagonal or above

      tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
      tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

      l_rows_tile = tile_size/np_rows ! local rows of a tile
      l_cols_tile = tile_size/np_cols ! local cols of a tile

      if (useQR) then

        if (useGPU) then
          print *,"qr decomposition at the moment not supported with GPU"
          stop
        endif

        if (which_qr_decomposition == 1) then
          call qr_pqrparam_init(pqrparam(1:11),    nblk,'M',0,   nblk,'M',0,   nblk,'M',1,'s')
          allocate(tauvector(na), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when allocating tauvector "//errorMessage
            stop
          endif

          allocate(blockheuristic(nblk), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when allocating blockheuristic "//errorMessage
            stop
          endif

          l_rows = local_index(na, my_prow, np_rows, nblk, -1)
          allocate(vmrCPU(max(l_rows,1),na), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when allocating vmrCPU "//errorMessage
            stop
          endif

          vmrCols = na


          call qr_pdgeqrf_2dcomm_double(a(1:lda,1:matrixCols), matrixCols, lda, vmrCPU(1:max(l_rows,1),1:vmrCols), max(l_rows,1), &
                                 vmrCols, tauvector(1:na), na, tmat(1:nbw,1:nbw,1), nbw, &
                                 nbw, dwork_size(1:1), 1, -1, na, nbw, nblk, nblk, na, na, 1, 0, PQRPARAM(1:11), &
                                 mpi_comm_rows, mpi_comm_cols, blockheuristic)


          work_size = dwork_size(1)
          allocate(work_blocked(work_size), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when allocating work_blocked "//errorMessage
            stop
          endif
          work_blocked = 0.0_rk8
          deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when deallocating vmrCPU "//errorMessage
            stop
          endif

        endif ! which_qr_decomposition

      endif ! useQr

      if (useGPU) then
! Here we convert the regular host array into a pinned host array
        successCUDA = cuda_malloc(a_dev, lda*na_cols*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"bandred_real: error in cudaMalloc"
          stop
        endif
        successCUDA = cuda_malloc(tmat_dev, nbw*nbw*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"bandred_real: error in cudaMalloc"
          stop
        endif
        successCUDA = cuda_malloc(vav_dev, nbw*nbw*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"bandred_real: error in cudaMalloc"
          stop
        endif

        cur_l_rows = 0
        cur_l_cols = 0

        if (.not.(successCUDA)) then
          print *,"bandred_real: error in cudaMemcpy"
          stop
        endif
      endif ! useGPU


      do istep = (na-1)/nbw, 1, -1

        n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

! Number of local columns/rows of remaining matrix
        l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
        l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

        if (useGPU) then
          cur_l_rows = max(l_rows, 1)
          cur_l_cols = max(l_cols, 1)

          vmr_size = cur_l_rows * 2 * n_cols
          umc_size = cur_l_cols * 2 * n_cols

! Allocate vmr and umc only if the inew size exceeds their current capacity
! Added for FORTRAN CALLS
          if ((.not. allocated(vr)) .or. (l_rows + 1 .gt. ubound(vr, dim=1))) then
            if (allocated(vr)) then
              deallocate(vr, stat=istat, errmsg=errorMessage)
              if (istat .ne. 0) then
                print *,"bandred_real: error when deallocating vr "//errorMessage
                stop
              endif
            endif
            allocate(vr(l_rows + 1), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              print *,"bandred_real: error when allocating vr "//errorMessage
              stop
            endif

          endif

          if ((.not. allocated(vmrCUDA)) .or. (vmr_size .gt. ubound(vmrCUDA, dim=1))) then
            if (allocated(vmrCUDA)) then
              deallocate(vmrCUDA, stat=istat, errmsg=errorMessage)
              if (istat .ne. 0) then
                print *,"bandred_real: error when allocating vmrCUDA "//errorMessage
                stop
              endif

              successCUDA = cuda_free(vmr_dev)
              if (.not.(successCUDA)) then
                print *,"bandred_real: error in cuda_free"
                stop
              endif
            endif

            allocate(vmrCUDA(vmr_size), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              print *,"bandred_real: error when allocating vmrCUDA "//errorMessage
              stop
            endif
            successCUDA = cuda_malloc(vmr_dev, vmr_size*size_of_double_real_datatype)
            if (.not.(successCUDA)) then
              print *,"bandred_real: error in cudaMalloc"
              stop
            endif

          endif

          if ((.not. allocated(umcCUDA)) .or. (umc_size .gt. ubound(umcCUDA, dim=1))) then
            if (allocated(umcCUDA)) then
              deallocate(umcCUDA, stat=istat, errmsg=errorMessage)
              if (istat .ne. 0) then
                print *,"bandred_real: error when deallocating umcCUDA "//errorMessage
                stop
              endif

              successCUDA = cuda_free(umc_dev)
              if (.not.(successCUDA)) then
                 print *,"bandred_real: error in cudaFree"
                 stop
              endif

            endif

            allocate(umcCUDA(umc_size), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              print *,"bandred_real: error when deallocating umcCUDA "//errorMessage
              stop
            endif
            successCUDA = cuda_malloc(umc_dev, umc_size*size_of_double_real_datatype)
            if (.not.(successCUDA)) then
              print *,"bandred_real: error in cudaMalloc"
              stop
            endif

          endif

        else ! GPU not used

! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

          allocate(vmrCPU(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when allocating vmrCPU "//errorMessage
            stop
          endif

          allocate(umcCPU(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when allocating umcCPU "//errorMessage
            stop
          endif

          allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_real: error when allocating vr "//errorMessage
            stop
          endif

        endif ! use GPU

        if (useGPU) then
          vmrCUDA(1 : cur_l_rows * n_cols) = 0.0_rk8
        else
          vmrCPU(1:l_rows,1:n_cols) = 0.0_rk8
        endif
        vr(:) = 0.0_rk8
        tmat(:,:,istep) = 0.0_rk8

        if (useGPU) then
          umcCUDA(1 : umc_size) = 0.0_rk8
          lc_start = local_index(istep*nbw+1, my_pcol, np_cols, nblk, -1)
          lc_end   = local_index(istep*nbw+n_cols, my_pcol, np_cols, nblk, -1)
          lr_end   = local_index((istep-1)*nbw + n_cols, my_prow, np_rows, nblk, -1)

          if(lc_start .le. 0) lc_start = 1

! Here we assume that the processor grid and the block grid are aligned
          cur_pcol = pcol(istep*nbw+1, nblk, np_cols)

          if(my_pcol == cur_pcol) then

            if (.not.(successCUDA)) then
              print *,"bandred_real: error in cudaMemcpy2d"
              stop
            endif

          endif
        endif ! useGPU

! Reduce current block to lower triangular form

        if (useQR) then
          if (which_qr_decomposition == 1) then
            vmrCols = 2*n_cols

            call qr_pdgeqrf_2dcomm_double(a(1:lda,1:matrixCols), lda, matrixCols, vmrCPU(1:max(l_rows,1),1:vmrCols) ,   &
                                    max(l_rows,1), vmrCols, tauvector(1:na), na, &
                                     tmat(1:nbw,1:nbw,istep), nbw, nbw, work_blocked(1:work_size), work_size, &
                                     work_size, na, n_cols, nblk, nblk,        &
                                     istep*nbw+n_cols-nbw, istep*nbw+n_cols, 1,&
                                     0, PQRPARAM(1:11), mpi_comm_rows, mpi_comm_cols,&
                                     blockheuristic)

          endif

       else !useQR

         do lc = n_cols, 1, -1

           ncol = istep*nbw + lc ! absolute column number of householder vector
           nrow = ncol - nbw ! Absolute number of pivot row

           lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
           lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

           tau = 0

           if (nrow == 1) exit ! Nothing to do

           cur_pcol = pcol(ncol, nblk, np_cols) ! Processor column owning current block

           if (my_pcol==cur_pcol) then

! Get vector to be transformed; distribute last element and norm of
! remaining elements to all procs in current column

             vr(1:lr) = a(1:lr,lch) ! vector to be transformed

             if (my_prow==prow(nrow, nblk, np_rows)) then
               aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
               aux1(2) = vr(lr)
             else
               aux1(1) = dot_product(vr(1:lr),vr(1:lr))
               aux1(2) = 0.0_rk8
             endif



             call mpi_allreduce(aux1, aux2, 2, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)




             vnorm2 = aux2(1)
             vrl    = aux2(2)

! Householder transformation
             call hh_transform_real_double(vrl, vnorm2, xf, tau)

! Scale vr and store Householder vector for back transformation

             vr(1:lr) = vr(1:lr) * xf
             if (my_prow==prow(nrow, nblk, np_rows)) then
               a(1:lr-1,lch) = vr(1:lr-1)
               a(lr,lch) = vrl
               vr(lr) = 1.0_rk8
             else
               a(1:lr,lch) = vr(1:lr)
             endif

           endif

! Broadcast Householder vector and tau along columns

           vr(lr+1) = tau


           call MPI_Bcast(vr, lr+1, MPI_REAL8, cur_pcol, mpi_comm_cols, mpierr)



           if (useGPU) then
             vmrCUDA(cur_l_rows * (lc - 1) + 1 : cur_l_rows * (lc - 1) + lr) = vr(1:lr)
           else
             vmrCPU(1:lr,lc) = vr(1:lr)
           endif

           tau = vr(lr+1)
           tmat(lc,lc,istep) = tau ! Store tau in diagonal of tmat

! Transform remaining columns in current block with Householder vector
! Local dot product

           aux1 = 0



           nlc = 0 ! number of local columns
           do j=1,lc-1
             lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
             if (lcx>0) then
               nlc = nlc+1
               if (lr>0) aux1(nlc) = dot_product(vr(1:lr),a(1:lr,lcx))
             endif
           enddo

! Get global dot products


           if (nlc>0) call mpi_allreduce(aux1, aux2, nlc, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)


! Transform

           nlc = 0
           do j=1,lc-1
             lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
             if (lcx>0) then
               nlc = nlc+1
               a(1:lr,lcx) = a(1:lr,lcx) - tau*aux2(nlc)*vr(1:lr)
             endif
           enddo

         enddo ! lc

         if (useGPU) then
! store column tiles back to GPU
           cur_pcol = pcol(istep*nbw+1, nblk, np_cols)
           if (my_pcol == cur_pcol) then

             if (.not.(successCUDA)) then
               print *,"bandred_real: error in cudaMemcpy2d"
               stop
             endif

           endif
         endif

! Calculate scalar products of stored Householder vectors.
! This can be done in different ways, we use dsyrk

         vav = 0

         if (useGPU) then
           if (l_rows>0) &
             call DSYRK('U', 'T', n_cols, l_rows, 1.0_rk8, vmrCUDA, cur_l_rows, 0.0_rk8, vav, ubound(vav,dim=1))
         else
           if (l_rows>0) &
             call DSYRK('U', 'T', n_cols, l_rows, 1.0_rk8, vmrCPU, ubound(vmrCPU,dim=1), 0.0_rk8, vav, ubound(vav,dim=1))
         endif

         call symm_matrix_allreduce_double(n_cols,vav, nbw, nbw,mpi_comm_rows)

! Calculate triangular matrix T for block Householder Transformation
         do lc=n_cols,1,-1
           tau = tmat(lc,lc,istep)
           if (lc<n_cols) then
             call DTRMV('U', 'T', 'N', n_cols-lc, tmat(lc+1,lc+1,istep), ubound(tmat,dim=1), vav(lc+1,lc), 1)
             tmat(lc,lc+1:n_cols,istep) = -tau * vav(lc+1:n_cols,lc)
           endif
         enddo
       endif

! Transpose vmr -> vmc (stored in umc, second half)

       if (useGPU) then
         call elpa_transpose_vectors_real_double  (vmrCUDA, cur_l_rows, mpi_comm_rows, &
                                            umcCUDA(cur_l_cols * n_cols + 1), cur_l_cols, mpi_comm_cols, &
                                            1, istep*nbw, n_cols, nblk)
       else
         call elpa_transpose_vectors_real_double  (vmrCPU, ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                            umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), mpi_comm_cols, &
                                            1, istep*nbw, n_cols, nblk)
       endif

! Calculate umc = A**T * vmr
! Note that the distributed A has to be transposed
! Opposed to direct tridiagonalization there is no need to use the cache locality
! of the tiles, so we can use strips of the matrix

! here the GPU version and CPU version diverged substantially, due to the newest
! optimizations due to Intel. The GPU version has to be re-written
       if (useGPU) then
         umcCUDA(1 : l_cols * n_cols) = 0.0_rk8
         vmrCUDA(cur_l_rows * n_cols + 1 : cur_l_rows * n_cols * 2) = 0.0_rk8

         if (l_cols>0 .and. l_rows>0) then

           if (.not.(successCUDA)) then
             print *,"bandred_real: error in cudaMemcpy"
             stop
           endif

           if (.not.(successCUDA)) then
             print *,"bandred_real: error in cudaMemcpy"
             stop
           endif

           do i=0,(istep*nbw-1)/tile_size

             lcs = i*l_cols_tile+1
             lce = min(l_cols,(i+1)*l_cols_tile)
             if (lce<lcs) cycle

             lre = min(l_rows,(i+1)*l_rows_tile)
             call cublas_dgemm('T', 'N', lce-lcs+1, n_cols, lre, &
                               1.0_rk8, (a_dev + ((lcs-1)*lda*size_of_double_real_datatype)), lda, vmr_dev,cur_l_rows, &
                               1.0_rk8, (umc_dev+ (lcs-1)*size_of_double_real_datatype), cur_l_cols)
             if(i==0) cycle
             lre = min(l_rows,i*l_rows_tile)
             call cublas_dgemm('N', 'N', lre,n_cols, lce-lcs+1,&
                               1.0_rk8, (a_dev+ ((lcs-1)*lda*size_of_double_real_datatype)), lda,                  &
                               (umc_dev+(cur_l_cols * n_cols+lcs-1)*size_of_double_real_datatype), cur_l_cols, &
                               1.0_rk8, (vmr_dev+(cur_l_rows * n_cols)*size_of_double_real_datatype), cur_l_rows)
           enddo

           if (.not.(successCUDA)) then
             print *,"bandred_real: error in cudaMemcpy"
             stop
           endif

           if (.not.(successCUDA)) then
             print *,"bandred_real: error in cudaMemcpy"
             stop
           endif

         endif ! l_cols>0 .and. l_rows>0

       else ! do not useGPU version

!Code for Algorithm 4

         n_way = 1

!umc(1:l_cols,1:n_cols) = 0.d0
!vmr(1:l_rows,n_cols+1:2*n_cols) = 0

         if (n_way > 1) then
!$omp do
           do i=1,min(l_cols_tile, l_cols)
             umcCPU(i,1:n_cols) = 0.0_rk8
           enddo

!$omp do
           do i=1,l_rows
             vmrCPU(i,n_cols+1:2*n_cols) = 0.0_rk8
           enddo
           if (l_cols>0 .and. l_rows>0) then

!SYMM variant 4
!Partitioned Matrix Expression:
! Ct = Atl Bt + Atr Bb
! Cb = Atr' Bt + Abl Bb
!
!Loop invariant:
! Ct = Atl Bt + Atr Bb
!
!Update:
! C1 = A10'B0 + A11B1 + A21 B2
!
!This algorithm chosen because in this algoirhtm, the loop around the dgemm calls
!is easily parallelized, and regardless of choise of algorithm,
!the startup cost for parallelizing the dgemms inside the loop is too great

!$omp do schedule(static,1)
             do i=0,(istep*nbw-1)/tile_size
               lcs = i*l_cols_tile+1                   ! local column start
               lce = min(l_cols, (i+1)*l_cols_tile)    ! local column end

               lrs = i*l_rows_tile+1                   ! local row start
               lre = min(l_rows, (i+1)*l_rows_tile)    ! local row end

!C1 += [A11 A12] [B1
!                 B2]
               if ( lre > lrs .and. l_cols > lcs ) then
                 call DGEMM('N', 'N', lre-lrs+1, n_cols, l_cols-lcs+1,          &
                            1.0_rk8, a(lrs,lcs), ubound(a,dim=1),                 &
                                  umcCPU(lcs,n_cols+1), ubound(umcCPU,dim=1),  &
                            0.0_rk8, vmrCPU(lrs,n_cols+1), ubound(vmrCPU,dim=1))
               endif

! C1 += A10' B0
               if ( lce > lcs .and. i > 0 ) then
                 call DGEMM('T', 'N', lce-lcs+1, n_cols, lrs-1,           &
                            1.0_rk8, a(1,lcs),   ubound(a,dim=1),           &
                                  vmrCPU(1,1),   ubound(vmrCPU,dim=1),   &
                            0.0_rk8, umcCPU(lcs,1), ubound(umcCPU,dim=1))
               endif
             enddo
           endif ! l_cols>0 .and. l_rows>0
         else ! n_way > 1
           umcCPU(1:l_cols,1:n_cols) = 0.0_rk8
           vmrCPU(1:l_rows,n_cols+1:2*n_cols) = 0.0_rk8
           if (l_cols>0 .and. l_rows>0) then
             do i=0,(istep*nbw-1)/tile_size

               lcs = i*l_cols_tile+1
               lce = min(l_cols,(i+1)*l_cols_tile)
               if (lce<lcs) cycle

               lre = min(l_rows,(i+1)*l_rows_tile)
               call DGEMM('T', 'N', lce-lcs+1, n_cols, lre, 1.0_rk8, a(1,lcs), ubound(a,dim=1), &
                            vmrCPU, ubound(vmrCPU,dim=1), 1.0_rk8, umcCPU(lcs,1), ubound(umcCPU,dim=1))

               if (i==0) cycle
                 lre = min(l_rows,i*l_rows_tile)
                 call DGEMM('N', 'N', lre, n_cols, lce-lcs+1, 1.0_rk8, a(1,lcs), lda, &
                            umcCPU(lcs,n_cols+1), ubound(umcCPU,dim=1), 1.0_rk8, vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1))
             enddo
           endif
         endif ! n_way > 1

       endif ! do not useGPU version

! Sum up all ur(:) parts along rows and add them to the uc(:) parts
! on the processors containing the diagonal
! This is only necessary if ur has been calculated, i.e. if the
! global tile size is smaller than the global remaining matrix

       if (useGPU) then
! here the GPU version and CPU version divereged due to the same reasons as above

         if (tile_size < istep*nbw) then
           call elpa_reduce_add_vectors_real_double  (vmrCUDA(cur_l_rows * n_cols + 1),cur_l_rows,mpi_comm_rows, &
                                               umcCUDA, cur_l_cols, mpi_comm_cols, &
                                               istep*nbw, n_cols, nblk)
         endif

         if (l_cols>0) then
           allocate(tmpCUDA(l_cols * n_cols), stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_real: error when allocating tmpCUDA "//errorMessage
             stop
           endif




           call mpi_allreduce(umcCUDA, tmpCUDA, l_cols*n_cols, MPI_REAL8, MPI_SUM, mpi_comm_rows, ierr)
           umcCUDA(1 : l_cols * n_cols) = tmpCUDA(1 : l_cols * n_cols)



           if (allocated(tmpCUDA)) then
             deallocate(tmpCUDA, stat=istat, errmsg=errorMessage)
             if (istat .ne. 0) then
               print *,"bandred_real: error when deallocating tmpCUDA "//errorMessage
               stop
             endif
           endif
         endif ! l_cols

! U = U * Tmat**T

         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif

         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif
         call cublas_dtrmm('Right', 'Upper', 'Trans', 'Nonunit', l_cols, n_cols, &
                           1.0_rk8, tmat_dev, nbw, umc_dev, cur_l_cols)
! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif
         call cublas_dgemm('T', 'N', n_cols, n_cols, l_cols, &
                           1.0_rk8, umc_dev, cur_l_cols, (umc_dev+(cur_l_cols * n_cols )*size_of_double_real_datatype),cur_l_cols, &
                           0.0_rk8, vav_dev, nbw)

         call cublas_dtrmm('Right', 'Upper', 'Trans', 'Nonunit', n_cols, n_cols, &
                           1.0_rk8, tmat_dev, nbw, vav_dev, nbw)


         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif

         call symm_matrix_allreduce_double(n_cols,vav, nbw,nbw,mpi_comm_cols)


         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif

! U = U - 0.5 * V * VAV
         call cublas_dgemm('N', 'N', l_cols, n_cols, n_cols,&
                           -0.5_rk8, (umc_dev+(cur_l_cols * n_cols )*size_of_double_real_datatype),cur_l_cols, vav_dev,nbw,&
                           1.0_rk8, umc_dev, cur_l_cols)


         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif

! Transpose umc -> umr (stored in vmr, second half)
         call elpa_transpose_vectors_real_double  (umcCUDA, cur_l_cols, mpi_comm_cols, &
                                            vmrCUDA(cur_l_rows * n_cols + 1), cur_l_rows, mpi_comm_rows, &
                                            1, istep*nbw, n_cols, nblk)


         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif


         if (.not.(successCUDA)) then
           print *,"bandred_real: error in cudaMemcpy"
           stop
         endif

! A = A - V*U**T - U*V**T
         do i=0,(istep*nbw-1)/tile_size
           lcs = i*l_cols_tile+1
           lce = min(l_cols,(i+1)*l_cols_tile)
           lre = min(l_rows,(i+1)*l_rows_tile)
           if (lce<lcs .or. lre<1) cycle
           call cublas_dgemm('N', 'T', lre, lce-lcs+1, 2*n_cols, -1.0_rk8, &
                             vmr_dev, cur_l_rows, (umc_dev +(lcs-1)*size_of_double_real_datatype), cur_l_cols, &
                             1.0_rk8, (a_dev+(lcs-1)*lda*size_of_double_real_datatype), lda)
         enddo

       else ! do not useGPU

! Or if we used the Algorithm 4
         if (tile_size < istep*nbw .or. n_way > 1) then
           call elpa_reduce_add_vectors_real_double  (vmrCPU(1,n_cols+1),ubound(vmrCPU,dim=1),mpi_comm_rows, &
                                             umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                             istep*nbw, n_cols, nblk)
         endif

         if (l_cols>0) then
           allocate(tmpCPU(l_cols,n_cols), stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_real: error when allocating tmpCPU "//errorMessage
             stop
           endif



           call mpi_allreduce(umcCPU, tmpCPU, l_cols*n_cols, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)
           umcCPU(1:l_cols,1:n_cols) = tmpCPU(1:l_cols,1:n_cols)



           deallocate(tmpCPU, stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_real: error when deallocating tmpCPU "//errorMessage
             stop
           endif
         endif

! U = U * Tmat**T
         call DTRMM('Right', 'Upper', 'Trans', 'Nonunit', l_cols,n_cols, 1.0_rk8, tmat(1,1,istep), ubound(tmat,dim=1), &
                    umcCPU, ubound(umcCPU,dim=1))

! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T

         call DGEMM('T', 'N', n_cols, n_cols, l_cols, 1.0_rk8, umcCPU, ubound(umcCPU,dim=1), umcCPU(1,n_cols+1), &
                    ubound(umcCPU,dim=1), 0.0_rk8, vav, ubound(vav,dim=1))

         call DTRMM('Right', 'Upper', 'Trans', 'Nonunit', n_cols, n_cols, 1.0_rk8, tmat(1,1,istep),    &
                    ubound(tmat,dim=1), vav, ubound(vav,dim=1))

         call symm_matrix_allreduce_double(n_cols,vav, nbw, nbw ,mpi_comm_cols)

! U = U - 0.5 * V * VAV
         call DGEMM('N', 'N', l_cols, n_cols, n_cols, -0.5_rk8, umcCPU(1,n_cols+1), ubound(umcCPU,dim=1), vav, &
                     ubound(vav,dim=1), 1.0_rk8, umcCPU, ubound(umcCPU,dim=1))

! Transpose umc -> umr (stored in vmr, second half)
         call elpa_transpose_vectors_real_double(umcCPU, ubound(umcCPU,dim=1), mpi_comm_cols, &
                                         vmrCPU(1,n_cols+1), ubound(vmrCPU,dim=1), mpi_comm_rows, &
                                         1, istep*nbw, n_cols, nblk)

! A = A - V*U**T - U*V**T



         do i=0,(istep*nbw-1)/tile_size
           lcs = i*l_cols_tile+1
           lce = min(l_cols,(i+1)*l_cols_tile)
           lre = min(l_rows,(i+1)*l_rows_tile)
           if (lce<lcs .or. lre<1) cycle
           call DGEMM('N', 'T', lre,lce-lcs+1, 2*n_cols, -1.0_rk8, &
                       vmrCPU, ubound(vmrCPU,dim=1), umcCPU(lcs,1), ubound(umcCPU,dim=1), &
                       1.0_rk8, a(1,lcs), lda)
         enddo


       endif ! useGPU

       if (.not.(useGPU)) then
         if (allocated(vr)) then
           deallocate(vr, stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_real: error when deallocating vr "//errorMessage
             stop
           endif
         endif

         if (allocated(umcCPU)) then
           deallocate(umcCPU, stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_real: error when deallocating vmrCPU "//errorMessage
             stop
           endif
         endif

         if (allocated(vmrCPU)) then
           deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_real: error when deallocating vmrCPU "//errorMessage
             stop
           endif
         endif

       endif !useGPU

     enddo ! istep

     if (useGPU) then
! this is not needed since a_dev is passed along from one subroutine to the other

       if (.not.(successCUDA)) then
         print *,"bandred_real: error in cudaMemcpy"
         stop
       endif

       successCUDA = cuda_free(a_dev)
       if (.not.(successCUDA)) then
         print *,"bandred_real: error in cudaFree"
         stop
       endif

!#ifdef 1
! it should be possible to keep tmat dev on the device and not copy it arround
! this is not necessary tmat_dev is passed (unchanged) from one routine to the other
       successCUDA = cuda_free(tmat_dev)
       if (.not.(successCUDA)) then
         print *,"bandred_real: error in cudaFree"
         stop
       endif
!#endif

       successCUDA = cuda_free(vav_dev)
       if (.not.(successCUDA)) then
         print *,"bandred_real: error in cudaFree"
         stop
       endif
     endif ! useGPU

     if (allocated(vr)) then
       deallocate(vr, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"bandred_real: error when deallocating vr "//errorMessage
         stop
       endif
     endif

     if (allocated(umcCPU)) then
       deallocate(umcCPU, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"bandred_real: error when deallocating umcCPU "//errorMessage
         stop
       endif
     endif

     if (allocated(vmrCPU)) then
       deallocate(vmrCPU, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"bandred_real: error when deallocating vmrCPU "//errorMessage
         stop
       endif
     endif

     if (useGPU) then
       successCUDA = cuda_free(vmr_dev)
       if (.not.(successCUDA)) then
         print *,"bandred_real: error in cudaFree"
         stop
       endif

       successCUDA = cuda_free(umc_dev)
       if (.not.(successCUDA)) then
         print *,"bandred_real: error in cudaFree"
         stop
       endif
       if (allocated(umcCUDA)) then
         deallocate(umcCUDA, stat=istat, errmsg=errorMessage)
         if (istat .ne. 0) then
           print *,"bandred_real: error when deallocating umcCUDA "//errorMessage
           stop
         endif
       endif
       if (allocated(vmrCUDA)) then
         deallocate(vmrCUDA, stat=istat, errmsg=errorMessage)
         if (istat .ne. 0) then
           print *,"bandred_real: error when deallocating vmrCUDA "//errorMessage
           stop
         endif
       endif

     endif ! useGPU

     if (useQR) then
       if (which_qr_decomposition == 1) then
         deallocate(work_blocked, stat=istat, errmsg=errorMessage)
         if (istat .ne. 0) then
           print *,"bandred_real: error when deallocating work_blocked "//errorMessage
           stop
         endif

         deallocate(tauvector, stat=istat, errmsg=errorMessage)
         if (istat .ne. 0) then
           print *,"bandred_real: error when deallocating tauvector "//errorMessage
           stop
         endif
       endif
     endif



   end subroutine bandred_real_double ! slower for gpu on 10000 10000 ???

    subroutine symm_matrix_allreduce_double(n,a,lda,ldb,comm)
!-------------------------------------------------------------------------------
!  symm_matrix_allreduce: Does an mpi_allreduce for a symmetric matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------

      use precision
      implicit none
      integer(kind=ik)             :: n, lda, ldb, comm

      real(kind=rk8)     :: a(lda,*)

      integer(kind=ik)             :: i, nc, mpierr
      real(kind=rk8)     :: h1(n*n), h2(n*n)



      nc = 0
      do i=1,n
        h1(nc+1:nc+i) = a(1:i,i)
        nc = nc+i
      enddo



      call mpi_allreduce(h1, h2, nc, MPI_REAL8, MPI_SUM, comm, mpierr)

      nc = 0
      do i=1,n
        a(1:i,i) = h2(nc+1:nc+i)
        a(i,1:i-1) = a(1:i-1,i)
        nc = nc+i
      enddo


!      nc = 0
!      do i=1,n
!        a(1:i,i) = h2(nc+1:nc+i)
!        a(i,1:i-1) = a(1:i-1,i)
!        nc = nc+i
!      enddo



    end subroutine symm_matrix_allreduce_double

    subroutine trans_ev_band_to_full_real_double(na, nqc, nblk, nbw, a, a_dev, lda, tmat, tmat_dev, q, q_dev, ldq, matrixCols, &
                                                       numBlocks, mpi_comm_rows, mpi_comm_cols, useGPU, useQR)
!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_real:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a(lda,matrixCols)    Matrix containing the Householder vectors (i.e. matrix a after bandred_real)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!  matrixCols  local columns of matrix a and q
!
!  tmat(nbw,nbw,numBlocks) Factors returned by bandred_real
!
!  q           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

      use precision
      use cuda_functions
      use, intrinsic :: iso_c_binding

      implicit none

      integer(kind=ik)                       :: na, nqc, lda, ldq, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols

      real(kind=rk8)               :: a(lda,*), q(ldq,*), tmat(nbw,nbw,*)


      integer(kind=C_intptr_T)               :: a_dev ! passed from bandred_real at the moment not used since copied in bandred_real

      integer(kind=ik)                       :: my_prow, my_pcol, np_rows, np_cols, mpierr
      integer(kind=ik)                       :: max_blocks_row, max_blocks_col, max_local_rows, &
                                                max_local_cols
      integer(kind=ik)                       :: l_cols, l_rows, l_colh, n_cols
      integer(kind=ik)                       :: istep, lc, ncol, nrow, nb, ns

      real(kind=rk8), allocatable  :: tmp1(:), tmp2(:), hvb(:), hvm(:,:)

! hvm_dev is fist used and set in this routine
! q is changed in trans_ev_tridi on the host, copied to device and passed here. this can be adapted
! tmp_dev is first used in this routine
! tmat_dev is passed along from bandred_real
      integer(kind=C_intptr_T)               :: hvm_dev, q_dev, tmp_dev, tmat_dev

      integer(kind=ik)                       :: i

      real(kind=rk8), allocatable  :: tmat_complete(:,:), t_tmp(:,:), t_tmp2(:,:)
      integer(kind=ik)                       :: cwy_blocking, t_blocking, t_cols, t_rows
      logical, intent(in)                    :: useQR, useGPU
      integer(kind=ik)                       :: istat
      character(200)                         :: errorMessage
      logical                                :: successCUDA

      successCUDA = .true.



      call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
      call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
      call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
      call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

      max_blocks_row = ((na -1)/nblk)/np_rows + 1  ! Rows of A
      max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

      max_local_rows = max_blocks_row*nblk
      max_local_cols = max_blocks_col*nblk

      if (useGPU) then

! here the GPU and CPU version diverged: the CPU version now always uses the useQR path which
! is not implemented in the GPU version
        allocate(tmp1(max_local_cols*nbw), stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_band_to_full_real: error when allocating tmp1 "//errorMessage
          stop
        endif

        allocate(tmp2(max_local_cols*nbw), stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_band_to_full_real: error when allocating tmp2 "//errorMessage
          stop
        endif

        allocate(hvb(max_local_rows*nbw), stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_band_to_full_real: error when allocating hvb "//errorMessage
          stop
        endif

        allocate(hvm(max_local_rows,nbw), stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_band_to_full_real: error when allocating hvm "//errorMessage
          stop
        endif
        successCUDA = cuda_malloc(hvm_dev, (max_local_rows)*nbw*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaMalloc"
          stop
        endif
        successCUDA = cuda_malloc(tmp_dev, (max_local_cols)*nbw*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaMalloc"
          stop
        endif

!#ifdef 1
! it should be possible to keep tmat dev on the device and not copy it around
! already existent on GPU
        successCUDA = cuda_malloc(tmat_dev, nbw*nbw*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaMalloc"
          stop
        endif
!#endif

! q_dev already living on device
!        successCUDA = cuda_malloc(q_dev, ldq*matrixCols*size_of_double_real_datatype)
!        if (.not.(successCUDA)) then
!          print *,"trans_ev_band_to_full_real: error in cudaMalloc"
!          stop
!        endif

!      q_temp(:,:) = 0.0
!      q_temp(1:ldq,1:na_cols) = q(1:ldq,1:na_cols)

!        ! copy q_dev to device, maybe this can be avoided if q_dev can be kept on device in trans_ev_tridi_to_band
!        successCUDA = cuda_memcpy(q_dev, loc(q), (ldq)*(matrixCols)*size_of_double_real_datatype, cudaMemcpyHostToDevice)
!        if (.not.(successCUDA)) then
!          print *,"trans_ev_band_to_full_real: error in cudaMalloc"
!          stop
!        endif

! if MPI is NOT used the following steps could be done on the GPU and memory transfers could be avoided
        successCUDA = cuda_memset(hvm_dev, 0, (max_local_rows)*(nbw)*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaMalloc"
          stop
        endif
        hvm = 0.0_rk8   ! Must be set to 0 !!!
        hvb = 0.0_rk8   ! Safety only
        
        l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

        do istep=1,(na-1)/nbw

          n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

! Broadcast all Householder vectors for current step compressed in hvb

          nb = 0
          ns = 0

          do lc = 1, n_cols
            ncol = istep*nbw + lc ! absolute column number of householder vector
            nrow = ncol - nbw ! absolute number of pivot row

            l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
            l_colh = local_index(ncol  , my_pcol, np_cols, nblk, -1) ! HV local column number

            if (my_pcol==pcol(ncol, nblk, np_cols)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

            nb = nb+l_rows

            if (lc==n_cols .or. mod(ncol,nblk)==0) then


              call MPI_Bcast(hvb(ns+1), nb-ns, MPI_REAL8, pcol(ncol, nblk, np_cols), mpi_comm_cols, mpierr)



              ns = nb
            endif
          enddo

! Expand compressed Householder vectors into matrix hvm

          nb = 0
          do lc = 1, n_cols
            nrow = (istep-1)*nbw+lc ! absolute number of pivot row
            l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

            hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
            if (my_prow==prow(nrow, nblk, np_rows)) hvm(l_rows+1,lc) = 1.0_rk8

            nb = nb+l_rows
          enddo

          if (.not.(successCUDA)) then
            print *,"trans_ev_band_to_full_real: error in cudaMemcpy"
            stop

          endif

          l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)

! Q = Q - V * T**T * V**T * Q

          if (l_rows>0) then
            call cublas_dgemm('T', 'N', n_cols, l_cols, l_rows, 1.0_rk8, hvm_dev, max_local_rows, &
                              q_dev, ldq , 0.0_rk8, tmp_dev, n_cols)


! copy data from device to host for a later MPI_ALLREDUCE

! copy to host maybe this can be avoided this is needed if MPI is used (allreduce)

            if (.not.(successCUDA)) then
              print *,"trans_ev_band_to_full_real: error in cudaMemcpy"
              stop
            endif


          else ! l_rows>0


            tmp1(1:l_cols*n_cols) = 0


          endif ! l_rows>0

!#ifdef WITH_GPU_VERSION
!       istat = cuda_memcpy(loc(tmp1), tmp_dev, max_local_cols*nbw*size_of_real_datatype,cudaMemcpyDeviceToHost)
!       if (istat .ne. 0) then
!         print *,"error in cudaMemcpy"
!         stop
!       endif
!#endif


          call mpi_allreduce(tmp1, tmp2, n_cols*l_cols, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)


!#ifdef WITH_GPU_VERSION
!       istat = cuda_memcpy(tmp_dev, loc(tmp2), max_local_cols*nbw*size_of_real_datatype,cudaMemcpyHostToDevice)
!       if (istat .ne. 0) then
!         print *,"error in cudaMemcpy"
!         stop
!       endif
!#endif

          if (l_rows>0) then

! after the mpi_allreduce we have to copy back to the device
! copy back to device

            if (.not.(successCUDA)) then
              print *,"trans_ev_band_to_full_real: error in cudaMemcpy"
              stop
            endif


!#ifdef 1
! it should be possible to keep tmat on the device and not copy it aroud
!            ! copy to device, maybe this can be avoided tmat is input from bandred_real


            if (.not.(successCUDA)) then
              print *,"trans_ev_band_to_full_real: error in cudaMemcpy"
              stop
            endif
!#endif /* 1 */

            call cublas_dtrmm('L', 'U', 'T', 'N', n_cols, l_cols, 1.0_rk8, tmat_dev, nbw, tmp_dev, n_cols)
            call cublas_dgemm('N', 'N', l_rows, l_cols, n_cols, -1.0_rk8, hvm_dev, max_local_rows, &
                              tmp_dev, n_cols, 1.0_rk8, q_dev, ldq)

! copy to host maybe this can be avoided
! this is not necessary hvm is not used anymore

            if (.not.(successCUDA)) then
              print *,"trans_ev_band_to_full_real: error in cudaMemcpy"
              stop
            endif

          endif ! l_rows > 0
!#ifdef WITH_GPU_VERSION
!       istat = cuda_memcpy(loc(hvm), hvm_dev, ((max_local_rows)*nbw*size_of_real_datatype),cudaMemcpyDeviceToHost)
!       if (istat .ne. 0) then
!         print *,"error in cudaMemcpy"
!         stop
!       endif
!
!#endif
        enddo ! istep

      else ! do not useGPU

! t_blocking was formerly 2; 3 is a better choice
        t_blocking = 3 ! number of matrices T (tmat) which are aggregated into a new (larger) T matrix (tmat_complete) and applied at once

! we only use the t_blocking if we could call it fully, this is might be better but needs to benchmarked.
!       if ( na >= ((t_blocking+1)*nbw) ) then
        cwy_blocking = t_blocking * nbw

        allocate(tmp1(max_local_cols*cwy_blocking))
        allocate(tmp2(max_local_cols*cwy_blocking))
        allocate(hvb(max_local_rows*cwy_blocking))
        allocate(hvm(max_local_rows,cwy_blocking))
        allocate(tmat_complete(cwy_blocking,cwy_blocking))
        allocate(t_tmp(cwy_blocking,nbw))
        allocate(t_tmp2(cwy_blocking,nbw))
!        else
!          allocate(tmp1(max_local_cols*nbw))
!          allocate(tmp2(max_local_cols*nbw))
!          allocate(hvb(max_local_rows*nbw))
!          allocate(hvm(max_local_rows,nbw))
!        endif
        hvm = 0.0_rk8   ! Must be set to 0 !!!
        hvb = 0.0_rk8   ! Safety only
        l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

!       if ( na >= ((t_blocking+1)*nbw) ) then

        do istep=1,((na-1)/nbw-1)/t_blocking + 1
! This the call when using  na >= ((t_blocking+1)*nbw)
!      n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw ! Number of columns in current step
! As an alternative we add some special case handling if na < cwy_blocking
          IF (na < cwy_blocking) THEN
            n_cols = MAX(0, na-nbw)
            IF ( n_cols .eq. 0 ) THEN
              EXIT
            END IF
          ELSE
            n_cols = MIN(na,istep*cwy_blocking+nbw) - (istep-1)*cwy_blocking - nbw ! Number of columns in current step
          END IF
! Broadcast all Householder vectors for current step compressed in hvb

          nb = 0
          ns = 0

          do lc = 1, n_cols
            ncol = (istep-1)*cwy_blocking + nbw + lc ! absolute column number of householder vector
            nrow = ncol - nbw ! absolute number of pivot row

            l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
            l_colh = local_index(ncol  , my_pcol, np_cols, nblk, -1) ! HV local column number

            if (my_pcol==pcol(ncol, nblk, np_cols)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

            nb = nb+l_rows

            if (lc==n_cols .or. mod(ncol,nblk)==0) then


              call MPI_Bcast(hvb(ns+1), nb-ns, MPI_REAL8, pcol(ncol, nblk, np_cols), mpi_comm_cols, mpierr)



              ns = nb
            endif
          enddo

! Expand compressed Householder vectors into matrix hvm

          nb = 0
          do lc = 1, n_cols
            nrow = (istep-1)*cwy_blocking + lc ! absolute number of pivot row
            l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

            hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
            if (my_prow==prow(nrow, nblk, np_rows)) hvm(l_rows+1,lc) = 1.0_rk8

            nb = nb+l_rows
          enddo

          l_rows = local_index(MIN(na,(istep+1)*cwy_blocking), my_prow, np_rows, nblk, -1)

! compute tmat2 out of tmat(:,:,)
          tmat_complete = 0
          do i = 1, t_blocking
            t_cols = MIN(nbw, n_cols - (i-1)*nbw)
            if (t_cols <= 0) exit
            t_rows = (i - 1) * nbw
            tmat_complete(t_rows+1:t_rows+t_cols,t_rows+1:t_rows+t_cols) = tmat(1:t_cols,1:t_cols,(istep-1)*t_blocking + i)
            if (i > 1) then
              call DGEMM('T', 'N', t_rows, t_cols, l_rows, 1.0_rk8, hvm(1,1), max_local_rows, hvm(1,(i-1)*nbw+1), &
                        max_local_rows, 0.0_rk8, t_tmp, cwy_blocking)



              call mpi_allreduce(t_tmp, t_tmp2, cwy_blocking*nbw, MPI_REAL8, MPI_SUM, mpi_comm_rows, mpierr)


              call DTRMM('L', 'U', 'N', 'N', t_rows, t_cols, 1.0_rk8, tmat_complete, cwy_blocking, t_tmp2, cwy_blocking)
              call DTRMM('R', 'U', 'N', 'N', t_rows, t_cols, -1.0_rk8, tmat_complete(t_rows+1,t_rows+1), cwy_blocking, &
                         t_tmp2, cwy_blocking)
              tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)



!              call DTRMM('L', 'U', 'N', 'N', t_rows, t_cols, 1.0_rk8, tmat_complete, cwy_blocking, t_tmp2, cwy_blocking)
!              call DTRMM('R', 'U', 'N', 'N', t_rows, t_cols, -1.0_rk8, tmat_complete(t_rows+1,t_rows+1), cwy_blocking, &
!                         t_tmp2, cwy_blocking)

!              tmat_complete(1:t_rows,t_rows+1:t_rows+t_cols) = t_tmp2(1:t_rows,1:t_cols)
             endif
          enddo

! Q = Q - V * T**T * V**T * Q

          if (l_rows>0) then
            call DGEMM('T', 'N', n_cols, l_cols, l_rows, 1.0_rk8, hvm, ubound(hvm,dim=1), &
                       q, ldq, 0.0_rk8, tmp1, n_cols)

          else ! l_rows>0

            tmp1(1:l_cols*n_cols) = 0.0_rk8
          endif ! l_rows>0



          call mpi_allreduce(tmp1, tmp2, n_cols*l_cols, MPI_REAL8, MPI_SUM, mpi_comm_rows ,mpierr)



          if (l_rows>0) then
            call DTRMM('L', 'U', 'T', 'N', n_cols, l_cols, 1.0_rk8, tmat_complete, cwy_blocking, tmp2, n_cols)
            call DGEMM('N', 'N', l_rows, l_cols, n_cols, -1.0_rk8, hvm, ubound(hvm,dim=1), tmp2, n_cols, 1.0_rk8, q, ldq)

          endif


!          if (l_rows>0) then
!            call DTRMM('L', 'U', 'T', 'N', n_cols, l_cols, 1.0_rk8, tmat_complete, cwy_blocking, tmp2, n_cols)
!            call DGEMM('N', 'N', l_rows, l_cols, n_cols, -1.0_rk8, hvm, ubound(hvm,dim=1), tmp2, n_cols, 1.0_rk8, q, ldq)
!          endif
        enddo ! istep

      endif ! useGPU

      deallocate(tmp1, tmp2, hvb, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_band_to_full_real: error when deallocating tmp1 tmp2 hvb "//errorMessage
        stop
      endif

      if (useGPU) then
        successCUDA = cuda_free(hvm_dev)
        if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaFree"
          stop
        endif

        successCUDA = cuda_free(tmp_dev)
        if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaFree"
          stop
        endif

        successCUDA = cuda_free(tmat_dev)
        if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaFree"
          stop
        endif
! final transfer of q_dev

         if (.not.(successCUDA)) then
          print *,"trans_ev_band_to_full_real: error in cudaFree"
          stop
         endif

!   q(1:ldq,1:na_cols) = q_temp(1:ldq,1:na_cols)

         successCUDA = cuda_free(q_dev)
         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_real: error in cudaFree"
           stop
         endif

!   deallocate(q_temp, stat=istat, errmsg=errorMessage)
!   if (istat .ne. 0) then
!     print *,"error when deallocating q_temp "//errorMessage
!     stop
!   endif
!   deallocate(tmat_temp, stat=istat, errmsg=errorMessage)
!   if (istat .ne. 0) then
!     print *,"trans_ev_band_to_full_real: error when deallocating tmat_temp "//errorMessage
!     stop
!   endif

      endif ! useGPU

      deallocate(hvm, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_band_to_full_real: error when deallocating hvm "//errorMessage
        stop
      endif

      if (useQr) then
        deallocate(tmat_complete, t_tmp, t_tmp2, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_band_to_full_real: error when deallocating tmat_complete, t_tmp, t_tmp2 "//errorMessage
          stop
        endif

      endif



    end subroutine trans_ev_band_to_full_real_double

    subroutine tridiag_band_real_double(na, nb, nblk, a, lda, d, e, matrixCols, hh_trans_real, &
                                 mpi_comm_rows, mpi_comm_cols, mpi_comm)
!-------------------------------------------------------------------------------
! tridiag_band_real:
! Reduces a real symmetric band matrix to tridiagonal form
!
!  na          Order of matrix a
!
!  nb          Semi bandwith
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  a(lda,matrixCols)    Distributed system matrix reduced to banded form in the upper diagonal
!
!  lda         Leading dimension of a
!  matrixCols  local columns of matrix a
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------

      use precision
      implicit none

      integer(kind=ik), intent(in)             ::  na, nb, nblk, lda, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm

      real(kind=rk8), intent(in)     :: a(lda,*)

      real(kind=rk8), intent(out)    :: d(na), e(na) ! set only on PE 0
      real(kind=rk8), intent(out), &
          allocatable                          :: hh_trans_real(:,:)

      real(kind=rk8)                 :: vnorm2, hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
      real(kind=rk8)                 :: hd(nb), hs(nb)

      integer(kind=ik)                         :: i, j, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
      integer(kind=ik)                         :: my_pe, n_pes, mpierr
      integer(kind=ik)                         :: my_prow, np_rows, my_pcol, np_cols
      integer(kind=ik)                         :: ireq_ab, ireq_hv
      integer(kind=ik)                         :: na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off

      integer(kind=ik), allocatable            :: ireq_hhr(:), ireq_hhs(:), global_id(:,:), hh_cnt(:), hh_dst(:)
      integer(kind=ik), allocatable            :: limits(:), snd_limits(:,:)
      integer(kind=ik), allocatable            :: block_limits(:)
      real(kind=rk8), allocatable    :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)


      integer                                  :: istat
      character(200)                           :: errorMessage





      call mpi_comm_rank(mpi_comm,my_pe,mpierr)
      call mpi_comm_size(mpi_comm,n_pes,mpierr)

      call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
      call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
      call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
      call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

! Get global_id mapping 2D procssor coordinates to global id

      allocate(global_id(0:np_rows-1,0:np_cols-1), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating global_id "//errorMessage
        stop
      endif


      global_id(:,:) = 0
      global_id(my_prow, my_pcol) = my_pe






      call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)




! Total number of blocks in the band:

      nblocks_total = (na-1)/nb + 1

! Set work distribution

      allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating block_limits"//errorMessage
        stop
      endif

      call divide_band(nblocks_total, n_pes, block_limits)

! nblocks: the number of blocks for my task
      nblocks = block_limits(my_pe+1) - block_limits(my_pe)

! allocate the part of the band matrix which is needed by this PE
! The size is 1 block larger than needed to avoid extensive shifts
      allocate(ab(2*nb,(nblocks+1)*nb), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating ab"//errorMessage
        stop
      endif

      ab = 0 ! needed for lower half, the extra block should also be set to 0 for safety

! n_off: Offset of ab within band
      n_off = block_limits(my_pe)*nb

! Redistribute band in a to ab
      call redist_band_real_double(a, lda, na, nblk, nb, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm, ab)
! Calculate the workload for each sweep in the back transformation
! and the space requirements to hold the HH vectors

      allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating limits"//errorMessage
        stop
      endif

      call determine_workload(na, nb, np_rows, limits)
      max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

      num_hh_vecs = 0
      num_chunks  = 0
      nx = na
      do n = 1, nblocks_total
        call determine_workload(nx, nb, np_rows, limits)
        local_size = limits(my_prow+1) - limits(my_prow)
! add to number of householder vectors
! please note: for nx==1 the one and only HH vector is 0 and is neither calculated nor send below!
        if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
          num_hh_vecs = num_hh_vecs + local_size
          num_chunks  = num_chunks+1
        endif
        nx = nx - nb
      enddo

! Allocate space for HH vectors

      allocate(hh_trans_real(nb,num_hh_vecs), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating hh_trans_real"//errorMessage
        stop
      endif


! Allocate and init MPI requests

      allocate(ireq_hhr(num_chunks), stat=istat, errmsg=errorMessage) ! Recv requests
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating ireq_hhr"//errorMessage
        stop
      endif
      allocate(ireq_hhs(nblocks), stat=istat, errmsg=errorMessage)    ! Send requests
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating ireq_hhs"//errorMessage
        stop
      endif

      num_hh_vecs = 0
      num_chunks  = 0
      nx = na
      nt = 0
      do n = 1, nblocks_total
        call determine_workload(nx, nb, np_rows, limits)
        local_size = limits(my_prow+1) - limits(my_prow)
        if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
          num_chunks  = num_chunks+1


          call mpi_irecv(hh_trans_real(1,num_hh_vecs+1), nb*local_size, MPI_REAL8, nt, &
                           10+n-block_limits(nt), mpi_comm, ireq_hhr(num_chunks), mpierr)



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
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating hh_gath"//errorMessage
        stop
      endif

      allocate(hh_send(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! send buffer for HH vectors
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating hh_send"//errorMessage
        stop
      endif
      hh_gath(:,:,:) = 0.0_rk8
      hh_send(:,:,:) = 0.0_rk8

! Some counters

      allocate(hh_cnt(nblocks), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating hh_cnt"//errorMessage
        stop
      endif

      allocate(hh_dst(nblocks), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating hh_dst"//errorMessage
        stop
      endif

      hh_cnt(:) = 1 ! The first transfomation vector is always 0 and not calculated at all
      hh_dst(:) = 0 ! PE number for receive

      ireq_ab = MPI_REQUEST_NULL
      ireq_hv = MPI_REQUEST_NULL

! Limits for sending

      allocate(snd_limits(0:np_rows,nblocks), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when allocating snd_limits"//errorMessage
        stop
      endif
      do iblk=1,nblocks
        call determine_workload(na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
      enddo



! ---------------------------------------------------------------------------
! Start of calculations

      na_s = block_limits(my_pe)*nb + 1

      if (my_pe>0 .and. na_s<=na) then
! send first column to previous PE
! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
        ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)


        call mpi_isend(ab_s, nb+1, MPI_REAL8, my_pe-1, 1, mpi_comm, ireq_ab, mpierr)


      endif





      do istep=1,na-1


        if (my_pe==0) then
          n = MIN(na-na_s,nb) ! number of rows to be reduced
          hv(:) = 0.0_rk8
          tau = 0.0_rk8
! The last step (istep=na-1) is only needed for sending the last HH vectors.
! We don't want the sign of the last element flipped (analogous to the other sweeps)
          if (istep < na-1) then
! Transform first column of remaining matrix
            vnorm2 = sum(ab(3:n+1,na_s-n_off)**2)
            call hh_transform_real_double(ab(2,na_s-n_off),vnorm2,hf,tau)
            hv(1) = 1.0_rk8
            hv(2:n) = ab(3:n+1,na_s-n_off)*hf
          endif
          d(istep) = ab(1,na_s-n_off)
          e(istep) = ab(2,na_s-n_off)
          if (istep == na-1) then
            d(na) = ab(1,na_s+1-n_off)
            e(na) = 0.0_rk8
          endif
        else
          if (na>na_s) then
! Receive Householder vector from previous task, from PE owning subdiagonal






            call mpi_recv(hv, nb, MPI_REAL8, my_pe-1, 2, mpi_comm, MPI_STATUS_IGNORE, mpierr)





            tau = hv(1)
            hv(1) = 1.0_rk8
          endif
        endif

        na_s = na_s+1
        if (na_s-n_off > nb) then
          ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
          ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rk8
          n_off = n_off + nb
        endif



          do iblk=1,nblocks
            ns = na_s + (iblk-1)*nb - n_off ! first column in block
            ne = ns+nb-1                    ! last column in block

            if (ns+n_off>na) exit

! Store Householder vector for back transformation

            hh_cnt(iblk) = hh_cnt(iblk) + 1

            hh_gath(1   ,hh_cnt(iblk),iblk) = tau
            hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)


            if (hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
! Wait for last transfer to finish



              call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)


! Copy vectors into send buffer
              hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
! Send to destination



              call mpi_isend(hh_send(1,1,iblk), nb*hh_cnt(iblk), MPI_REAL8, &
                           global_id(hh_dst(iblk), mod(iblk+block_limits(my_pe)-1,np_cols)), &
                           10+iblk, mpi_comm, ireq_hhs(iblk), mpierr)



! Reset counter and increase destination row
              hh_cnt(iblk) = 0.0_rk8
              hh_dst(iblk) = hh_dst(iblk)+1
            endif

! The following code is structured in a way to keep waiting times for
! other PEs at a minimum, especially if there is only one block.
! For this reason, it requests the last column as late as possible
! and sends the Householder vector and the first column as early
! as possible.

            nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
            nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
! Note that nr>=0 implies that diagonal block is full (nc==nb)!

! Multiply diagonal block and subdiagonal block with Householder vector

            if (iblk==nblocks .and. nc==nb) then

! We need the last column from the next PE.
! First do the matrix multiplications without last column ...

! Diagonal block, the contribution of the last element is added below!
              ab(1,ne) = 0.0_rk8

              call DSYMV('L', nc, tau, ab(1,ns), 2*nb-1, hv, 1, 0.0_rk8, hd, 1)

! Subdiagonal block
              if (nr>0) call DGEMV('N', nr, nb-1, tau, ab(nb+1,ns), 2*nb-1, hv, 1, 0.0_rk8, hs, 1)

! ... then request last column ...



              call mpi_recv(ab(1,ne), nb+1, MPI_REAL8, my_pe+1, 1, mpi_comm, MPI_STATUS_IGNORE, mpierr)




! ... and complete the result
              hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
              hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

            else

! Normal matrix multiply
              call DSYMV('L', nc, tau, ab(1,ns), 2*nb-1, hv, 1, 0.0_rk8, hd, 1)
              if (nr>0) call DGEMV('N', nr, nb, tau, ab(nb+1,ns), 2*nb-1, hv, 1, 0.0_rk8, hs, 1)
            endif

! Calculate first column of subdiagonal block and calculate new
! Householder transformation for this column
            hv_new(:) = 0.0_rk8 ! Needed, last rows must be 0 for nr < nb
            tau_new = 0.0_rk8
            if (nr>0) then

! complete (old) Householder transformation for first column

              ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

! calculate new Householder transformation ...
              if (nr>1) then
                vnorm2 = sum(ab(nb+2:nb+nr,ns)**2)
                call hh_transform_real_double(ab(nb+1,ns),vnorm2,hf,tau_new)
                hv_new(1) = 1.0_rk8
                hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
                ab(nb+2:,ns) = 0.0_rk8
              endif

! ... and send it away immediatly if this is the last block

              if (iblk==nblocks) then



                call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)




                hv_s(1) = tau_new
                hv_s(2:) = hv_new(2:)



                call mpi_isend(hv_s, nb, MPI_REAL8, my_pe+1, 2, mpi_comm, ireq_hv, mpierr)



              endif

            endif

! Transform diagonal block
            x = dot_product(hv(1:nc),hd(1:nc))*tau
            hd(1:nc) = hd(1:nc) - 0.5_rk8*x*hv(1:nc)

            if (my_pe>0 .and. iblk==1) then

! The first column of the diagonal block has to be send to the previous PE
! Calculate first column only ...

              ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*hv(1) - hv(1:nc)*hd(1)

! ... send it away ...




              call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)




              ab_s(1:nb+1) = ab(1:nb+1,ns)




              call mpi_isend(ab_s, nb+1, MPI_REAL8, my_pe-1, 1, mpi_comm, ireq_ab, mpierr)



! ... and calculate remaining columns with rank-2 update
              if (nc>1) call DSYR2('L', nc-1, -1.0_rk8, hd(2), 1, hv(2), 1, ab(1,ns+1), 2*nb-1)
            else
! No need to  send, just a rank-2 update
              call DSYR2('L', nc, -1.0_rk8, hd, 1, hv, 1, ab(1,ns), 2*nb-1)
            endif

! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

            if (nr>0) then
              if (nr>1) then
                call DGEMV('T', nr, nb-1, tau_new, ab(nb,ns+1), 2*nb-1, hv_new, 1, 0.0_rk8, h(2), 1)
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

! Use new HH vector for the next block
            hv(:) = hv_new(:)
            tau = tau_new

          enddo


      enddo ! istep

! Finish the last outstanding requests





      call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
      call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

      call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
      call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)







      call mpi_barrier(mpi_comm,mpierr)


      deallocate(ab, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when deallocating ab"//errorMessage
        stop
      endif

      deallocate(ireq_hhr, ireq_hhs, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_real: error when deallocating ireq_hhr, ireq_hhs"//errorMessage
        stop
      endif

      deallocate(hh_cnt, hh_dst, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"tridiag_band_real: error when deallocating hh_cnt, hh_dst"//errorMessage
         stop
       endif

      deallocate(hh_gath, hh_send, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"tridiag_band_real: error when deallocating hh_gath, hh_send"//errorMessage
         stop
       endif

      deallocate(limits, snd_limits, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"tridiag_band_real: error when deallocating limits, send_limits"//errorMessage
         stop
       endif

      deallocate(block_limits, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"tridiag_band_real: error when deallocating block_limits"//errorMessage
         stop
       endif

      deallocate(global_id, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"tridiag_band_real: error when allocating global_id"//errorMessage
         stop
       endif


    end subroutine tridiag_band_real_double

    subroutine trans_ev_tridi_to_band_real_double(na, nev, nblk, nbw, q, q_dev, ldq, matrixCols, hh_trans_real, &
                                           mpi_comm_rows, mpi_comm_cols, wantDebug, useGPU, success, &
                                           THIS_REAL_ELPA_KERNEL)
!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_real:
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

      use cuda_functions
      use precision
      use pack_unpack_real
      use pack_unpack_real_gpu
      use compute_hh_trafo_real
      use, intrinsic :: iso_c_binding
      implicit none
      logical, intent(in) :: useGPU

      integer(kind=ik), intent(in)             :: THIS_REAL_ELPA_KERNEL
      integer(kind=ik), intent(in)             :: na, nev, nblk, nbw, ldq, matrixCols, mpi_comm_rows, mpi_comm_cols

      real(kind=rk8)                 :: q(ldq,*)

      integer(kind=c_intptr_t)                 :: q_dev


      real(kind=rk8), intent(in)     :: hh_trans_real(:,:)
      integer(kind=ik)                         :: np_rows, my_prow, np_cols, my_pcol

      integer(kind=ik)                         :: i, j, ip, sweep, nbuf, l_nev, a_dim2
      integer(kind=ik)                         :: current_n, current_local_n, current_n_start, current_n_end
      integer(kind=ik)                         :: next_n, next_local_n, next_n_start, next_n_end
      integer(kind=ik)                         :: bottom_msg_length, top_msg_length, next_top_msg_length
      integer(kind=ik)                         :: stripe_width, last_stripe_width, stripe_count

      integer(kind=ik)                         :: num_result_blocks, num_result_buffers, num_bufs_recvd
      integer(kind=ik)                         :: a_off, current_tv_off, max_blk_size
      integer(kind=ik)                         :: mpierr, src, src_offset, dst, offset, nfact, num_blk



      logical                                  :: flag


      real(kind=rk8), pointer        :: aIntern(:,:,:)

      real(kind=rk8)                 :: a_real
      type(c_ptr)                              :: aIntern_ptr
      real(kind=rk8)   , allocatable :: row(:)
      real(kind=rk8)   , allocatable :: row_group(:,:)


      real(kind=rk8), allocatable    :: top_border_send_buffer(:,:,:), top_border_recv_buffer(:,:,:)
      real(kind=rk8), allocatable    :: bottom_border_send_buffer(:,:,:), bottom_border_recv_buffer(:,:,:)

      real(kind=rk8), allocatable    :: result_buffer(:,:,:)
      real(kind=rk8), allocatable    :: bcast_buffer(:,:)
      integer(kind=ik)                         :: tmp

!      real*8, allocatable, device              :: a_dev(:,:,:)
!      real*8, allocatable, device              :: bcast_buffer_dev(:,:)
!      real*8, allocatable, device              :: row_dev(:)
!      real*8, allocatable, device              :: row_group_dev(:,:)
!      real*8, allocatable, device              :: hh_dot_dev(:)
!      real*8, allocatable, device              :: hh_tau_dev(:)

      integer(kind=c_intptr_t)                 :: aIntern_dev
      integer(kind=c_intptr_t)                 :: bcast_buffer_dev
      integer(kind=c_size_t)                   :: num
      integer(kind=c_size_t)                   :: dev_offset, dev_offset_1


      integer(kind=c_intptr_t)                 :: row_dev
      integer(kind=c_intptr_t)                 :: row_group_dev
      integer(kind=c_intptr_t)                 :: hh_dot_dev
      integer(kind=c_intptr_t)                 :: hh_tau_dev
      Integer(kind=ik)                         :: top, chunk, this_chunk
      integer(kind=ik)                         :: row_group_size, unpack_idx

      integer(kind=ik)                         :: n_off
      integer(kind=ik), allocatable            :: result_send_request(:), result_recv_request(:), limits(:)
      integer(kind=ik), allocatable            :: top_send_request(:), bottom_send_request(:)
      integer(kind=ik), allocatable            :: top_recv_request(:), bottom_recv_request(:)

! MPI send/recv tags, arbitrary

      integer(kind=ik), parameter              :: bottom_recv_tag = 111
      integer(kind=ik), parameter              :: top_recv_tag    = 222
      integer(kind=ik), parameter              :: result_recv_tag = 333

! Just for measuring the kernel performance
      real(kind=c_double)                      :: kernel_time, kernel_time_recv ! MPI_WTIME always needs double
! long integer
      integer(kind=lik)                        :: kernel_flops, kernel_flops_recv



      logical, intent(in)                      :: wantDebug
      logical                                  :: success
      integer(kind=ik)                         :: istat
      character(200)                           :: errorMessage
      logical                                  :: successCUDA

      integer(kind=C_SIZE_T)                   :: aux_int


      successCUDA = .true.




      success = .true.
      kernel_time = 0.0
      kernel_flops = 0



      call MPI_Comm_rank(mpi_comm_rows, my_prow, mpierr)
      call MPI_Comm_size(mpi_comm_rows, np_rows, mpierr)
      call MPI_Comm_rank(mpi_comm_cols, my_pcol, mpierr)
      call MPI_Comm_size(mpi_comm_cols, np_cols, mpierr)

      if (mod(nbw,nblk)/=0) then
        if (my_prow==0 .and. my_pcol==0) then
          if (wantDebug) then
            write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_real: ERROR: nbw=',nbw,', nblk=',nblk
            write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_real: band backtransform works only for nbw==n*nblk'
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
      else

! Suggested stripe width is 48 since 48*64 real*8 numbers should fit into
! every primary cache
        if (.not.(useGPU)) then



          stripe_width = 48 ! Must be a multiple of 4


          stripe_count = (l_nev-1)/stripe_width + 1

! Adapt stripe width so that last one doesn't get too small

          stripe_width = (l_nev-1)/stripe_count + 1


          stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 4 because of AVX/SSE memory alignment of 32 bytes
! (4 * sizeof(double) == 32)

        else ! GPUs are used
          stripe_width = 256 ! Must be a multiple of 4
          stripe_count = (l_nev - 1) / stripe_width + 1
        endif

        last_stripe_width = l_nev - (stripe_count-1)*stripe_width
      endif

! Determine the matrix distribution at the beginning

      allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_real: error when allocating limits"//errorMessage
        stop
      endif
      call determine_workload(na, nbw, np_rows, limits)

      max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

      a_dim2 = max_blk_size + nbw

      if (useGPU) then
        num =  (stripe_width*a_dim2*stripe_count)*size_of_double_real_datatype
        successCUDA = cuda_malloc(aIntern_dev, stripe_width*a_dim2*stripe_count*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc"//errorMessage
          stop
        endif

        successCUDA = cuda_memset(aIntern_dev , 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMemset"//errorMessage
          stop
        endif

      else ! GPUs are not used

        aux_int = int(stripe_width*a_dim2*stripe_count*8,kind=C_SIZE_T)

        if (posix_memalign(aIntern_ptr, 64_C_SIZE_T, aux_int) /= 0) then
          print *,"trans_ev_tridi_to_band_real: error when allocating aIntern"//errorMessage
          stop
        endif

        call c_f_pointer(aIntern_ptr, aIntern,[stripe_width,a_dim2,stripe_count] )
!allocate(aIntern(stripe_width,a_dim2,stripe_count), stat=istat, errmsg=errorMessage)

        aIntern(:,:,:) = 0.0_rk8


      endif !useGPU

      allocate(row(l_nev), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_real: error when allocating row"//errorMessage
        stop
      endif
      row(:) = 0.0_rk8
      if (useGPU) then
        num =  (l_nev)*size_of_double_real_datatype
        successCUDA = cuda_malloc( row_dev,num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc "//errorMessage
          stop
        endif

        successCUDA = cuda_memset(row_dev , 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMemset"//errorMessage
          stop
        endif

! "row_group" and "row_group_dev" are needed for GPU optimizations
        allocate(row_group(l_nev, nblk), stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_real: error when allocating row_group"//errorMessage
          stop
        endif

        row_group(:, :) = 0.0_rk8

        num =  (l_nev*nblk)*size_of_double_real_datatype
!    call cuda_malloc2d( row_group_dev,l_nev*size_of_double_real_datatype,nblk*size_of_real_datatype)
        successCUDA = cuda_malloc(row_group_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc"//errorMessage
          stop
        endif
        successCUDA = cuda_memset(row_group_dev , 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMemset"//errorMessage
          stop
        endif

      endif ! useGPU

! Copy q from a block cyclic distribution into a distribution with contiguous rows,
! and transpose the matrix using stripes of given stripe_width for cache blocking.

! The peculiar way it is done below is due to the fact that the last row should be
! ready first since it is the first one to start below



      do ip = np_rows-1, 0, -1
        if (my_prow == ip) then
! Receive my rows which have not yet been received
          src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
          do i=limits(ip)+1,limits(ip+1)
            src = mod((i-1)/nblk, np_rows)
            if (src < my_prow) then

              if (useGPU) then
! An unpacking of the current row group may occur before queuing the next row
                call unpack_and_prepare_row_group_real_gpu_double(row_group, row_group_dev, aIntern_dev, stripe_count, &
                                                           stripe_width, last_stripe_width, a_dim2, l_nev,&
                                                           row_group_size, nblk, unpack_idx, i - limits(ip), .false.)


                call MPI_Recv(row_group(:, row_group_size), l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)




              else


                call MPI_Recv(row, l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)



                call unpack_row_real_cpu_double(aIntern, row,i-limits(ip), stripe_count, stripe_width, last_stripe_width)
              endif


            elseif (src==my_prow) then
              src_offset = src_offset+1
              if (.not.(useGPU)) row(:) = q(src_offset, 1:l_nev)



              if (useGPU) then
! An unpacking of the current row group may occur before queuing the next row
                call unpack_and_prepare_row_group_real_gpu_double(row_group, row_group_dev, aIntern_dev, stripe_count, &
                                                           stripe_width, last_stripe_width, a_dim2, l_nev,&
                                                           row_group_size, nblk, unpack_idx, i - limits(ip), .false.)
                row_group(:, row_group_size) = q(src_offset, 1:l_nev)
              else
                call unpack_row_real_cpu_double(aIntern, row,i-limits(ip),  stripe_count, stripe_width, last_stripe_width)
              endif



            endif
          enddo

! Send all rows which have not yet been send
          src_offset = 0
          do dst = 0, ip-1
            do i=limits(dst)+1,limits(dst+1)
              if (mod((i-1)/nblk, np_rows) == my_prow) then
                src_offset = src_offset+1
                row(:) = q(src_offset, 1:l_nev)


                call MPI_Send(row, l_nev, MPI_REAL8, dst, 0, mpi_comm_rows, mpierr)


              endif
            enddo
          enddo

        else if (my_prow < ip) then
! Send all rows going to PE ip
          src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
          do i=limits(ip)+1,limits(ip+1)
            src = mod((i-1)/nblk, np_rows)
            if (src == my_prow) then
              src_offset = src_offset+1
              row(:) = q(src_offset, 1:l_nev)


              call MPI_Send(row, l_nev, MPI_REAL8, ip, 0, mpi_comm_rows, mpierr)


            endif
          enddo
! Receive all rows from PE ip
          do i=limits(my_prow)+1,limits(my_prow+1)
            src = mod((i-1)/nblk, np_rows)
            if (src == ip) then

              if (useGPU) then
! An unpacking of the current row group may occur before queuing the next row
                call unpack_and_prepare_row_group_real_gpu_double(row_group, row_group_dev, aIntern_dev, stripe_count,  &
                                                          stripe_width, last_stripe_width, a_dim2, l_nev,       &
                                  row_group_size, nblk, unpack_idx, i - limits(my_prow), .false.)


               call MPI_Recv(row_group(:, row_group_size), l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)


              else


                call MPI_Recv(row, l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)


                call unpack_row_real_cpu_double(aIntern, row,i-limits(my_prow), stripe_count, stripe_width, last_stripe_width)
              endif



            endif
          enddo
        endif
      enddo

      if (useGPU) then
! Force an unpacking of all remaining rows that haven't been unpacked yet
        call unpack_and_prepare_row_group_real_gpu_double(row_group, row_group_dev, aIntern_dev, stripe_count, &
                                                      stripe_width, last_stripe_width, &
                                                   a_dim2, l_nev, row_group_size, nblk, unpack_idx, -1, .true.)
        successCUDA = cuda_devicesynchronize()

         if (.not.(successCUDA)) then
           print *,"trans_ev_tridi_to_band_real: error in cudaDeviceSynchronize"//errorMessage
           stop
         endif
      endif

! Set up result buffer queue

      num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

      num_result_buffers = 4*nfact
      allocate(result_buffer(l_nev,nblk,num_result_buffers), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_real: error when allocating result_buffer"//errorMessage
        stop
      endif

      allocate(result_send_request(num_result_buffers), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_real: error when allocating result_send_request"//errorMessage
        stop
      endif

      allocate(result_recv_request(num_result_buffers), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_real: error when allocating result_recv_request"//errorMessage
        stop
      endif


      result_send_request(:) = MPI_REQUEST_NULL
      result_recv_request(:) = MPI_REQUEST_NULL

! Queue up buffers

      if (my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
        do j = 1, min(num_result_buffers, num_result_blocks)


          call MPI_Irecv(result_buffer(1,1,j), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                              mpi_comm_rows, result_recv_request(j), mpierr)


        enddo
      endif

      num_bufs_recvd = 0 ! No buffers received yet

! Initialize top/bottom requests

      allocate(top_send_request(stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating top_send_request"//errorMessage
         stop
       endif

      allocate(top_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating top_recv_request"//errorMessage
         stop
       endif

      allocate(bottom_send_request(stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating bottom_send_request"//errorMessage
         stop
       endif

      allocate(bottom_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating bottom_recv_request"//errorMessage
         stop
       endif


      top_send_request(:) = MPI_REQUEST_NULL
      top_recv_request(:) = MPI_REQUEST_NULL
      bottom_send_request(:) = MPI_REQUEST_NULL
      bottom_recv_request(:) = MPI_REQUEST_NULL




       allocate(top_border_send_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating top_border_send_bufer"//errorMessage
         stop
       endif

      allocate(top_border_recv_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating top_border_recv_buffer"//errorMessage
         stop
       endif

      allocate(bottom_border_send_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating bottom_border_send_buffer"//errorMessage
         stop
       endif

      allocate(bottom_border_recv_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating bottom_border_recv_buffer"//errorMessage
         stop
       endif

      top_border_send_buffer(:,:,:) = 0.0_rk8
      top_border_recv_buffer(:,:,:) = 0.0_rk8
      bottom_border_send_buffer(:,:,:) = 0.0_rk8
      bottom_border_recv_buffer(:,:,:) = 0.0_rk8



      allocate(bcast_buffer(nbw, max_blk_size), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when allocating bcast_buffer"//errorMessage
         stop
       endif

      bcast_buffer = 0.0_rk8

      if (useGPU) then
        num =  ( nbw * max_blk_size) * size_of_double_real_datatype
        successCUDA = cuda_malloc(bcast_buffer_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc"
          stop
        endif

        successCUDA = cuda_memset( bcast_buffer_dev, 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMemset"
          stop
        endif

        num =  ((max_blk_size-1))*size_of_double_real_datatype
        successCUDA = cuda_malloc( hh_dot_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc"
          stop
        endif

        successCUDA = cuda_memset( hh_dot_dev, 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMemset"
          stop
        endif

        num =  (max_blk_size)*size_of_double_real_datatype
        successCUDA = cuda_malloc( hh_tau_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc"
          stop
        endif

        successCUDA = cuda_memset( hh_tau_dev, 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMemset"
          stop
        endif
      endif ! useGPU

      current_tv_off = 0 ! Offset of next row to be broadcast

! ------------------- start of work loop -------------------

      a_off = 0 ! offset in aIntern (to avoid unnecessary shifts)

      top_msg_length = 0
      bottom_msg_length = 0

      do sweep = 0, (na-1)/nbw

        current_n = na - sweep*nbw
        call determine_workload(current_n, nbw, np_rows, limits)
        current_n_start = limits(my_prow)
        current_n_end   = limits(my_prow+1)
        current_local_n = current_n_end - current_n_start

        next_n = max(current_n - nbw, 0)
        call determine_workload(next_n, nbw, np_rows, limits)
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
          do i = 1, stripe_count





            call MPI_Irecv(bottom_border_recv_buffer(1,1,i), nbw*stripe_width, MPI_REAL8, my_prow+1, bottom_recv_tag, &
                        mpi_comm_rows, bottom_recv_request(i), mpierr)





          enddo
        endif

        if (current_local_n > 1) then
          if (my_pcol == mod(sweep,np_cols)) then
            bcast_buffer(:,1:current_local_n) = hh_trans_real(:,current_tv_off+1:current_tv_off+current_local_n)
            current_tv_off = current_tv_off + current_local_n
          endif


          call mpi_bcast(bcast_buffer, nbw*current_local_n, MPI_REAL8, mod(sweep,np_cols), mpi_comm_cols, mpierr)




          if (useGPU) then

            if (.not.(successCUDA)) then
              print *,"trans_ev_tridi_to_band_real: error in cudaMemcpy"
              stop
            endif

            call extract_hh_tau_real_gpu_double(bcast_buffer_dev, hh_tau_dev, nbw, current_local_n, .false.)
            call compute_hh_dot_products_real_gpu_double(bcast_buffer_dev, hh_dot_dev, nbw, current_local_n)
          endif

        else ! (current_local_n > 1) then

! for current_local_n == 1 the one and only HH vector is 0 and not stored in hh_trans_real
          bcast_buffer(:,1) = 0.0_rk8

          if (useGPU) then
            successCUDA = cuda_memset(bcast_buffer_dev, 0, nbw * size_of_double_real_datatype)
            if (.not.(successCUDA)) then
              print *,"trans_ev_tridi_to_band_real: error in cudaMemset"
              stop
            endif

            call extract_hh_tau_real_gpu_double(bcast_buffer_dev, hh_tau_dev, nbw, 1, .true.)
          endif
        endif ! (current_local_n > 1) then

        if (l_nev == 0) cycle

        if (current_local_n > 0) then

          do i = 1, stripe_count


!wait_b
            if (current_n_end < current_n) then





              call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)



              n_off = current_local_n+a_off

              if (useGPU) then
                dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) *size_of_double_real_datatype

                if (.not.(successCUDA)) then
                  print *,"trans_ev_tridi_to_band_real: error in cudaMemcpy"
                  stop
                endif

              else
                aIntern(:,n_off+1:n_off+nbw,i) = bottom_border_recv_buffer(:,1:nbw,i)
              endif



           if (next_n_end < next_n) then





             call MPI_Irecv(bottom_border_recv_buffer(1,1,i), nbw*stripe_width, MPI_REAL8, my_prow+1, bottom_recv_tag, &
                                   mpi_comm_rows, bottom_recv_request(i), mpierr)





           endif
         endif

         if (current_local_n <= bottom_msg_length + top_msg_length) then

!wait_t
           if (top_msg_length>0) then




             call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)



             if (useGPU) then
               dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *size_of_double_real_datatype
!             host_offset= (0 + (0 * stripe_width) + ( (i-1) * stripe_width * nbw ) ) * 8

               if (.not.(successCUDA)) then
                 print *,"trans_ev_tridi_to_band_real: error in cudaMemcpy"
                 stop
                endif
             else
               aIntern(:,a_off+1:a_off+top_msg_length,i) = top_border_recv_buffer(:,1:top_msg_length,i)
             endif ! useGPU

           endif ! top_msg_length

!compute

           call compute_hh_trafo_real_cpu_double(aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count,       &
                                          a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, hh_dot_dev, &
                                          hh_tau_dev, kernel_flops, kernel_time, 0, current_local_n, i,          &
                                          last_stripe_width, THIS_REAL_ELPA_KERNEL)


!send_b




           call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)


           if (bottom_msg_length>0) then
             n_off = current_local_n+nbw-bottom_msg_length+a_off

             if (useGPU) then
               dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *size_of_double_real_datatype

               if (.not.(successCUDA)) then
                 print *,"trans_ev_tridi_to_band_real: error in cudaMemcpy"
                 stop
               endif
             else
               bottom_border_send_buffer(:,1:bottom_msg_length,i) = aIntern(:,n_off+1:n_off+bottom_msg_length,i)
             endif


             call MPI_Isend(bottom_border_send_buffer(1,1,i), bottom_msg_length*stripe_width, MPI_REAL8, my_prow+1, &
                            top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)



           endif

         else ! current_local_n <= bottom_msg_length + top_msg_length

!compute

        call compute_hh_trafo_real_cpu_double(aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count,       &
                                       a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, hh_dot_dev, &
                                       hh_tau_dev, kernel_flops, kernel_time,                                 &
                                       current_local_n - bottom_msg_length, bottom_msg_length, i,             &
                                       last_stripe_width, THIS_REAL_ELPA_KERNEL)
!send_b



        call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)


        if (bottom_msg_length > 0) then
          n_off = current_local_n+nbw-bottom_msg_length+a_off

          if (useGPU) then
            dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *size_of_double_real_datatype

            if (.not.(successCUDA)) then
              print *,"trans_ev_tridi_to_band_real: error cudaMemcpy"
              stop
            endif
          else
            bottom_border_send_buffer(:,1:bottom_msg_length,i) = aIntern(:,n_off+1:n_off+bottom_msg_length,i)
          endif


          call MPI_Isend(bottom_border_send_buffer(1,1,i), bottom_msg_length*stripe_width, MPI_REAL8, my_prow+1, &
                         top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)


        endif


!compute

        call compute_hh_trafo_real_cpu_double(aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count,           &
                                       a_off,  nbw, max_blk_size, bcast_buffer, bcast_buffer_dev, hh_dot_dev,     &
                                       hh_tau_dev, kernel_flops, kernel_time,  top_msg_length,                    &
                                       current_local_n-top_msg_length-bottom_msg_length, i,                       &
                                       last_stripe_width, THIS_REAL_ELPA_KERNEL)


!wait_t
        if (top_msg_length>0) then




          call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)


          if (useGPU) then
            dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *size_of_double_real_datatype

            if (.not.(successCUDA)) then
              print *,"trans_ev_tridi_to_band_real: error in cudaMemcpy"
              stop
            endif
          else
            aIntern(:,a_off+1:a_off+top_msg_length,i) = top_border_recv_buffer(:,1:top_msg_length,i)
          endif

        endif

!compute

        call compute_hh_trafo_real_cpu_double(aIntern, aIntern_dev, stripe_width, a_dim2, stripe_count,           &
                                       a_off, nbw, max_blk_size,  bcast_buffer, bcast_buffer_dev, hh_dot_dev,     &
                                       hh_tau_dev, kernel_flops, kernel_time, 0, top_msg_length, i,               &
                                       last_stripe_width, THIS_REAL_ELPA_KERNEL)

      endif

      if (next_top_msg_length > 0) then
!request top_border data




        call MPI_Irecv(top_border_recv_buffer(1,1,i), next_top_msg_length*stripe_width, MPI_REAL8, my_prow-1, &
                       top_recv_tag, mpi_comm_rows, top_recv_request(i), mpierr)





      endif

!send_t
      if (my_prow > 0) then




        call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)


        if (useGPU) then
          dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * size_of_double_real_datatype

          if (.not.(successCUDA)) then
            print *,"trans_ev_tridi_to_band_real: error in cudaMemcpy"
            stop
          endif

        else
          top_border_send_buffer(:,1:nbw,i) = aIntern(:,a_off+1:a_off+nbw,i)
        endif


        call MPI_Isend(top_border_send_buffer(1,1,i), nbw*stripe_width, MPI_REAL8, my_prow-1, bottom_recv_tag, &
                       mpi_comm_rows, top_send_request(i), mpierr)




      endif

! Care that there are not too many outstanding top_recv_request's
      if (stripe_count > 1) then
        if (i>1) then





          call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)




        else





          call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)




        endif
      endif

    enddo

    top_msg_length = next_top_msg_length

  else
! wait for last top_send_request
    do i = 1, stripe_count




      call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)




      enddo
    endif

! Care about the result

    if (my_prow == 0) then

! topmost process sends nbw rows to destination processes

      do j=0, nfact-1
        num_blk = sweep*nfact+j ! global number of destination block, 0 based
        if (num_blk*nblk >= na) exit

        nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block






        call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)





        dst = mod(num_blk, np_rows)

        if (dst == 0) then
          if (useGPU) then
            row_group_size = min(na - num_blk*nblk, nblk)
            call pack_row_group_real_gpu_double(row_group_dev, aIntern_dev, stripe_count, stripe_width, &
                                     last_stripe_width, a_dim2, l_nev, &
                                         row_group(:, :), j * nblk + a_off, row_group_size)
            do i = 1, row_group_size
              q((num_blk / np_rows) * nblk + i, 1 : l_nev) = row_group(:, i)
            enddo
          else ! useGPU

            do i = 1, min(na - num_blk*nblk, nblk)

              call pack_row_real_cpu_double(aIntern, row, j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)

              q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
            enddo
          endif ! useGPU

        else ! (dst == 0)

          if (useGPU) then
            call pack_row_group_real_gpu_double(row_group_dev, aIntern_dev, stripe_count, stripe_width, &
                                       last_stripe_width, a_dim2, l_nev, &
                                           result_buffer(:, :, nbuf), j * nblk + a_off, nblk)
          else  ! useGPU
            do i = 1, nblk

              call pack_row_real_cpu_double(aIntern, result_buffer(:,i,nbuf),j*nblk+i+a_off, stripe_width, &
                                      last_stripe_width, stripe_count)

            enddo
          endif ! useGPU


          call MPI_Isend(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, dst, &
                                    result_recv_tag, mpi_comm_rows, result_send_request(nbuf), mpierr)



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





          call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)






          if (.not.flag) exit

        else ! (next_local_n > 0)






          call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)






        endif ! (next_local_n > 0)

! Fill result buffer into q
        num_blk = j*np_rows + my_prow ! global number of current block, 0 based
        do i = 1, min(na - num_blk*nblk, nblk)
          q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
        enddo

! Queue result buffer again if there are outstanding blocks left



        if (j+num_result_buffers < num_result_blocks) &
            call MPI_Irecv(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                         mpi_comm_rows, result_recv_request(nbuf), mpierr)
! carefull the "recieve" has to be done at the corresponding wait or send
!         if (j+num_result_buffers < num_result_blocks) &
!                result_buffer(1:l_nev*nblk,1,nbuf) =  result_buffer(1:l_nev*nblk,1,nbuf)




      enddo ! j = num_bufs_recvd, num_result_blocks-1
      num_bufs_recvd = j

    endif ! (my_prow == 0)

! Shift the remaining rows to the front of aIntern (if necessary)

    offset = nbw - top_msg_length
    if (offset<0) then
      if (wantDebug) write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_real: internal error, offset for shifting = ',offset
      success = .false.
      return
    endif

    a_off = a_off + offset
    if (a_off + next_local_n + nbw > a_dim2) then

         do i = 1, stripe_count
           if (useGPU) then
             chunk = min(next_local_n - 1, a_off)
             do j = top_msg_length + 1, top_msg_length + next_local_n, chunk
               top = min(j + chunk, top_msg_length + next_local_n)
               this_chunk = top - j + 1
               dev_offset = (0 + ( (j-1) * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *size_of_double_real_datatype
               dev_offset_1 = (0 + ( (j + a_off-1) * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * &
                             size_of_double_real_datatype
! it is not logical to set here always the value for the parameter
! "cudaMemcpyDeviceToDevice" do this ONCE at startup
!               tmp = cuda_d2d(1)

               if (.not.(successCUDA)) then
                 print *,"trans_ev_tridi_to_band_real: error cudaMemcpy"
                 stop
               endif
             enddo
           else ! not useGPU
             do j = top_msg_length+1, top_msg_length+next_local_n
               aIntern(:,j,i) = aIntern(:,j+a_off,i)
             enddo
           endif
         enddo ! stripe_count



         a_off = 0
       endif

     enddo

! Just for safety:

     if (ANY(top_send_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_send_request ***',my_prow,my_pcol
     if (ANY(bottom_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_send_request ***',my_prow,my_pcol
     if (ANY(top_recv_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_recv_request ***',my_prow,my_pcol
     if (ANY(bottom_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_recv_request ***',my_prow,my_pcol

     if (my_prow == 0) then






       call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)





     endif

     if (ANY(result_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_send_request ***',my_prow,my_pcol
     if (ANY(result_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_recv_request ***',my_prow,my_pcol





! copy q to q_dev needed in trans_ev_band_to_full
        successCUDA = cuda_malloc(q_dev, ldq*matrixCols*size_of_double_real_datatype)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc"
          stop
        endif

! copy q_dev to device, maybe this can be avoided if q_dev can be kept on device in trans_ev_tridi_to_band

        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_real: error in cudaMalloc"
          stop
        endif

! deallocate all working space

     if (.not.(useGPU)) then
       nullify(aIntern)
       call free(aIntern_ptr)
!       deallocate(aIntern, stat=istat, errmsg=errorMessage)
!       if (istat .ne. 0) then
!         print *,"trans_ev_tridi_to_band_real: error when deallocating aIntern "//errorMessage
!         stop
!       endif
     endif

     deallocate(row, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating row "//errorMessage
       stop
     endif

     deallocate(limits, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating limits"//errorMessage
       stop
     endif

     deallocate(result_send_request, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating result_send_request "//errorMessage
       stop
     endif

     deallocate(result_recv_request, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating result_recv_request "//errorMessage
       stop
     endif

     deallocate(top_border_send_buffer, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating top_border_send_buffer "//errorMessage
       stop
     endif

     deallocate(top_border_recv_buffer, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating top_border_recv_buffer "//errorMessage
       stop
     endif

     deallocate(bottom_border_send_buffer, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating bottom_border_send_buffer "//errorMessage
       stop
     endif

     deallocate(bottom_border_recv_buffer, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating bottom_border_recv_buffer "//errorMessage
       stop
     endif

     deallocate(result_buffer, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating result_buffer "//errorMessage
       stop
     endif

     deallocate(bcast_buffer, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating bcast_buffer "//errorMessage
       stop
     endif

     deallocate(top_send_request, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating top_send_request "//errorMessage
       stop
     endif

     deallocate(top_recv_request, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating top_recv_request "//errorMessage
       stop
     endif

     deallocate(bottom_send_request, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating bottom_send_request "//errorMessage
       stop
     endif

     deallocate(bottom_recv_request, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"trans_ev_tridi_to_band_real: error when deallocating bottom_recv_request "//errorMessage
       stop
     endif

     if (useGPU) then
       successCUDA = cuda_free(hh_dot_dev)
       if (.not.(successCUDA)) then
         print *,"trans_ev_tridi_to_band_real: error in cudaFree "//errorMessage
         stop
       endif

       successCUDA = cuda_free(hh_tau_dev)
       if (.not.(successCUDA)) then
         print *,"trans_ev_tridi_to_band_real: error in cudaFree "//errorMessage
         stop
       endif

       successCUDA = cuda_free(row_dev)
       if (.not.(successCUDA)) then
         print *,"trans_ev_tridi_to_band_real: error in cudaFree "//errorMessage
         stop
       endif

       deallocate(row_group, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_tridi_to_band_real: error when deallocating row_group "//errorMessage
         stop
       endif

       successCUDA = cuda_free(row_group_dev)
       if (.not.(successCUDA)) then
         print *,"trans_ev_tridi_to_band_real: error in cudaFree "//errorMessage
         stop
       endif

       successCUDA =  cuda_free(bcast_buffer_dev)
       if (.not.(successCUDA)) then
         print *,"trans_ev_tridi_to_band_real: error in cudaFree "//errorMessage
         stop
       endif
     endif ! useGPU



     return

    end subroutine trans_ev_tridi_to_band_real_double


    subroutine determine_workload(na, nb, nprocs, limits)

      use precision
      implicit none

      integer(kind=ik), intent(in)  :: na, nb, nprocs
      integer(kind=ik), intent(out) :: limits(0:nprocs)

      integer(kind=ik)              :: i



      if (na <= 0) then
        limits(:) = 0

        return
      endif

      if (nb*nprocs > na) then
! there is not enough work for all
        do i = 0, nprocs
          limits(i) = min(na, i*nb)
        enddo
      else
         do i = 0, nprocs
           limits(i) = (i*na)/nprocs
         enddo
      endif


    end subroutine

!---------------------------------------------------------------------------------------------------
! divide_band: sets the work distribution in band
! Proc n works on blocks block_limits(n)+1 .. block_limits(n+1)

    subroutine divide_band(nblocks_total, n_pes, block_limits)

      use precision
      implicit none
      integer(kind=ik), intent(in)  :: nblocks_total ! total number of blocks in band
      integer(kind=ik), intent(in)  :: n_pes         ! number of PEs for division
      integer(kind=ik), intent(out) :: block_limits(0:n_pes)

      integer(kind=ik)              :: n, nblocks, nblocks_left



      block_limits(0) = 0
      if (nblocks_total < n_pes) then
! Not enough work for all: The first tasks get exactly 1 block
        do n=1,n_pes
          block_limits(n) = min(nblocks_total,n)
        enddo
      else
! Enough work for all. If there is no exact loadbalance,
! the LAST tasks get more work since they are finishing earlier!
        nblocks = nblocks_total/n_pes
        nblocks_left = nblocks_total - n_pes*nblocks
        do n=1,n_pes
          if (n<=n_pes-nblocks_left) then
            block_limits(n) = block_limits(n-1) + nblocks
          else
            block_limits(n) = block_limits(n-1) + nblocks + 1
          endif
        enddo
      endif



    end subroutine



    subroutine band_band_real_double(na, nb, nbCol, nb2, nb2Col, ab, ab2, d, e, mpi_comm)
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
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------

      use precision
      implicit none

      integer(kind=ik), intent(in)             :: na, nb, nbCol, nb2, nb2Col, mpi_comm
      real(kind=rk8), intent(inout)  :: ab(2*nb,nbCol) ! removed assumed size
      real(kind=rk8), intent(inout)  :: ab2(2*nb2,nb2Col) ! removed assumed size
      real(kind=rk8), intent(out)    :: d(na), e(na) ! set only on PE 0

      real(kind=rk8)                 :: hv(nb,nb2), w(nb,nb2), w_new(nb,nb2), tau(nb2), hv_new(nb,nb2), &
                                                  tau_new(nb2), ab_s(1+nb,nb2), ab_r(1+nb,nb2), ab_s2(2*nb2,nb2), hv_s(nb,nb2)

      real(kind=rk8)                 :: work(nb*nb2), work2(nb2*nb2)
      integer(kind=ik)                         :: lwork, info

      integer(kind=ik)                         :: istep, i, n, dest
      integer(kind=ik)                         :: n_off, na_s
      integer(kind=ik)                         :: my_pe, n_pes, mpierr
      integer(kind=ik)                         :: nblocks_total, nblocks
      integer(kind=ik)                         :: nblocks_total2, nblocks2
      integer(kind=ik)                         :: ireq_ab, ireq_hv

!      integer(kind=ik)                         :: MPI_STATUS_IGNORE(MPI_STATUS_SIZE)

!      integer(kind=ik), allocatable            :: mpi_statuses(:,:)
      integer(kind=ik), allocatable            :: block_limits(:), block_limits2(:), ireq_ab2(:)

      integer(kind=ik)                         :: j, nc, nr, ns, ne, iblk
      integer(kind=ik)                         :: istat
      character(200)                           :: errorMessage


!      if (na .lt. 2*nb) then
!        print *,"na lt 2*nb ",na,2*nb
!        stop
!      endif
!      if (na .lt. 2*nb2) then
!        print *,"na lt 2*nb2 ",na,2*nb2
!        stop
!      endif
!      if (na .lt. nbCol) then
!        print *,"na lt nbCol ",na,nbCol
!        stop
!      endif
!      if (na .lt. nb2Col) then
!        print *,"na lt nb2Col ",na,nb2Col
!        stop
!      endif


      call mpi_comm_rank(mpi_comm,my_pe,mpierr)
      call mpi_comm_size(mpi_comm,n_pes,mpierr)


! Total number of blocks in the band:
      nblocks_total = (na-1)/nb + 1
      nblocks_total2 = (na-1)/nb2 + 1

! Set work distribution
      allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating block_limits "//errorMessage
        stop
      endif
      call divide_band(nblocks_total, n_pes, block_limits)

      allocate(block_limits2(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating block_limits2 "//errorMessage
        stop
      endif

      call divide_band(nblocks_total2, n_pes, block_limits2)

! nblocks: the number of blocks for my task
      nblocks = block_limits(my_pe+1) - block_limits(my_pe)
      nblocks2 = block_limits2(my_pe+1) - block_limits2(my_pe)

      allocate(ireq_ab2(1:nblocks2), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"error allocating ireq_ab2 "//errorMessage
        stop
      endif




      ireq_ab2 = MPI_REQUEST_NULL

      if (nb2>1) then
        do i=0,nblocks2-1

          call mpi_irecv(ab2(1,i*nb2+1), 2*nb2*nb2, MPI_REAL8, 0, 3, mpi_comm, ireq_ab2(i+1), mpierr)
        enddo
      endif



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



        call mpi_isend(ab_s, (nb+1)*nb2, MPI_REAL8, my_pe-1, 1, mpi_comm, ireq_ab, mpierr)


      endif

      do istep=1,na/nb2

        if (my_pe==0) then

          n = MIN(na-na_s-nb2+1,nb) ! number of rows to be reduced
          hv(:,:) = 0.0_rk8
          tau(:) = 0.0_rk8

! The last step (istep=na-1) is only needed for sending the last HH vectors.
! We don't want the sign of the last element flipped (analogous to the other sweeps)
          if (istep < na/nb2) then

! Transform first block column of remaining matrix
            call dgeqrf(n, nb2, ab(1+nb2,na_s-n_off), 2*nb-1, tau, work, lwork, info)

            do i=1,nb2
              hv(i,i) = 1.0_rk8
              hv(i+1:n,i) = ab(1+nb2+1:1+nb2+n-i,na_s-n_off+i-1)
              ab(1+nb2+1:2*nb,na_s-n_off+i-1) = 0.0_rk8
            enddo

          endif

          if (nb2==1) then
            d(istep) = ab(1,na_s-n_off)
            e(istep) = ab(2,na_s-n_off)
            if (istep == na) then
              e(na) = 0.0_rk8
            endif
          else
            ab_s2 = 0.0_rk8
            ab_s2(:,:) = ab(1:nb2+1,na_s-n_off:na_s-n_off+nb2-1)
            if (block_limits2(dest+1)<istep) then
              dest = dest+1
            endif


            call mpi_send(ab_s2, 2*nb2*nb2, MPI_REAL8, dest, 3, mpi_comm, mpierr)




          endif

        else
          if (na>na_s+nb2-1) then
! Receive Householder vectors from previous task, from PE owning subdiagonal


            call mpi_recv(hv, nb*nb2, MPI_REAL8, my_pe-1, 2, mpi_comm, MPI_STATUS_IGNORE, mpierr)




            do i=1,nb2
              tau(i) = hv(i,i)
              hv(i,i) = 1.0_rk8
            enddo
          endif
        endif

        na_s = na_s+nb2
        if (na_s-n_off > nb) then
          ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
          ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0.0_rk8
          n_off = n_off + nb
        endif

        do iblk=1,nblocks
          ns = na_s + (iblk-1)*nb - n_off ! first column in block
          ne = ns+nb-nb2                    ! last column in block

          if (ns+n_off>na) exit

            nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
            nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
! Note that nr>=0 implies that diagonal block is full (nc==nb)!
            call wy_gen_double(nc,nb2,w,hv,tau,work,nb)

            if (iblk==nblocks .and. nc==nb) then
!request last nb2 columns


              call mpi_recv(ab_r,(nb+1)*nb2, MPI_REAL8, my_pe+1, 1, mpi_comm, MPI_STATUS_IGNORE, mpierr)



              do i=1,nb2
                ab(1:nb+1,ne+i-1) = ab_r(:,i)
              enddo
            endif
            hv_new(:,:) = 0.0_rk8 ! Needed, last rows must be 0 for nr < nb
            tau_new(:) = 0.0_rk8

            if (nr>0) then
              call wy_right_double(nr,nb,nb2,ab(nb+1,ns),2*nb-1,w,hv,work,nb)
              call dgeqrf(nr, nb2, ab(nb+1,ns), 2*nb-1, tau_new, work, lwork, info)
              do i=1,nb2
                hv_new(i,i) = 1.0_rk8
                hv_new(i+1:,i) = ab(nb+2:2*nb-i+1,ns+i-1)
                ab(nb+2:,ns+i-1) = 0.0_rk8
              enddo

!send hh-vector
              if (iblk==nblocks) then



                call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)



                hv_s = hv_new
                do i=1,nb2
                  hv_s(i,i) = tau_new(i)
                enddo


                call mpi_isend(hv_s,nb*nb2, MPI_REAL8, my_pe+1, 2, mpi_comm, ireq_hv, mpierr)



              endif
            endif

            call wy_symm_double(nc,nb2,ab(1,ns),2*nb-1,w,hv,work,work2,nb)

            if (my_pe>0 .and. iblk==1) then
!send first nb2 columns to previous PE



              call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)



              do i=1,nb2
                ab_s(1:nb+1,i) = ab(1:nb+1,ns+i-1)
              enddo


              call mpi_isend(ab_s,(nb+1)*nb2, MPI_REAL8, my_pe-1, 1, mpi_comm, ireq_ab, mpierr)



            endif

            if (nr>0) then
              call wy_gen_double(nr,nb2,w_new,hv_new,tau_new,work,nb)
              call wy_left_double(nb-nb2,nr,nb2,ab(nb+1-nb2,ns+nb2),2*nb-1,w_new,hv_new,work,nb)
            endif

! Use new HH vector for the next block
            hv(:,:) = hv_new(:,:)
            tau = tau_new
          enddo
        enddo

! Finish the last outstanding requests



        call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
        call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
!        allocate(mpi_statuses(MPI_STATUS_SIZE,nblocks2), stat=istat, errmsg=errorMessage)
!        if (istat .ne. 0) then
!          print *,"error allocating mpi_statuses "//errorMessage
!          stop
!        endif

        call mpi_waitall(nblocks2,ireq_ab2,MPI_STATUSES_IGNORE,mpierr)
!        deallocate(mpi_statuses, stat=istat, errmsg=errorMessage)
!        if (istat .ne. 0) then
!          print *,"error deallocating mpi_statuses "//errorMessage
!          stop
!        endif

        call mpi_barrier(mpi_comm,mpierr)




        deallocate(block_limits, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating block_limits "//errorMessage
          stop
        endif

        deallocate(block_limits2, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating block_limits2 "//errorMessage
          stop
        endif

        deallocate(ireq_ab2, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"error deallocating ireq_ab2 "//errorMessage
          stop
        endif



    end subroutine

    subroutine wy_gen_double(n, nb, W, Y, tau, mem, lda)


      use precision
      implicit none
      integer(kind=ik), intent(in)            :: n      !length of householder-vectors
      integer(kind=ik), intent(in)            :: nb     !number of householder-vectors
      integer(kind=ik), intent(in)            :: lda        !leading dimension of Y and W
      real(kind=rk8), intent(in)    :: Y(lda,nb)  !matrix containing nb householder-vectors of length b
      real(kind=rk8), intent(in)    :: tau(nb)    !tau values
      real(kind=rk8), intent(out)   :: W(lda,nb)  !output matrix W
      real(kind=rk8), intent(in)    :: mem(nb)    !memory for a temporary matrix of size nb

      integer(kind=ik)             :: i



   W(1:n,1) = tau(1)*Y(1:n,1)
   do i=2,nb
     W(1:n,i) = tau(i)*Y(1:n,i)
     call DGEMV('T', n, i-1,  1.0_rk8, Y, lda, W(1,i), 1, 0.0_rk8, mem,1)
     call DGEMV('N', n, i-1, -1.0_rk8, W, lda, mem, 1, 1.0_rk8, W(1,i),1)
   enddo

    end subroutine

    subroutine wy_left_double(n, m, nb, A, lda, W, Y, mem, lda2)


      use precision
      implicit none
      integer(kind=ik), intent(in)            :: n      !width of the matrix A
      integer(kind=ik), intent(in)            :: m      !length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk8), intent(inout) :: A(lda,*)   !matrix to be transformed   ! remove assumed size
      real(kind=rk8), intent(in)    :: W(m,nb)    !blocked transformation matrix W
      real(kind=rk8), intent(in)    :: Y(m,nb)    !blocked transformation matrix Y
      real(kind=rk8), intent(inout) :: mem(n,nb)  !memory for a temporary matrix of size n x nb



   call DGEMM('T', 'N', nb, n, m, 1.0_rk8, W, lda2, A, lda, 0.0_rk8, mem, nb)
   call DGEMM('N', 'N', m, n, nb, -1.0_rk8, Y, lda2, mem, nb, 1.0_rk8, A, lda)


    end subroutine

    subroutine wy_right_double(n, m, nb, A, lda, W, Y, mem, lda2)


      use precision
      implicit none
      integer(kind=ik), intent(in)            :: n      !height of the matrix A
      integer(kind=ik), intent(in)            :: m      !length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk8), intent(inout) :: A(lda,*)   !matrix to be transformed  ! remove assumed size
      real(kind=rk8), intent(in)    :: W(m,nb)    !blocked transformation matrix W
      real(kind=rk8), intent(in)    :: Y(m,nb)    !blocked transformation matrix Y
      real(kind=rk8), intent(inout) :: mem(n,nb)  !memory for a temporary matrix of size n x nb




   call DGEMM('N', 'N', n, nb, m, 1.0_rk8, A, lda, W, lda2, 0.0_rk8, mem, n)
   call DGEMM('N', 'T', n, m, nb, -1.0_rk8, mem, n, Y, lda2, 1.0_rk8, A, lda)



    end subroutine

    subroutine wy_symm_double(n, nb, A, lda, W, Y, mem, mem2, lda2)


      use precision
      implicit none
      integer(kind=ik), intent(in)            :: n      !width/heigth of the matrix A; length of matrix W and Y
      integer(kind=ik), intent(in)            :: nb     !width of matrix W and Y
      integer(kind=ik), intent(in)            :: lda        !leading dimension of A
      integer(kind=ik), intent(in)            :: lda2       !leading dimension of W and Y
      real(kind=rk8), intent(inout) :: A(lda,*)   !matrix to be transformed  ! remove assumed size
      real(kind=rk8), intent(in)    :: W(n,nb)    !blocked transformation matrix W
      real(kind=rk8), intent(in)    :: Y(n,nb)    !blocked transformation matrix Y
      real(kind=rk8)                :: mem(n,nb)  !memory for a temporary matrix of size n x nb
      real(kind=rk8)                :: mem2(nb,nb)    !memory for a temporary matrix of size nb x nb



   call DSYMM('L', 'L', n, nb, 1.0_rk8, A, lda, W, lda2, 0.0_rk8, mem, n)
   call DGEMM('T', 'N', nb, nb, n, 1.0_rk8, mem, n, W, lda2, 0.0_rk8, mem2, nb)
   call DGEMM('N', 'N', n, nb, nb, -0.5_rk8, Y, lda2, mem2, nb, 1.0_rk8, mem, n)
   call DSYR2K('L', 'N', n, nb, -1.0_rk8, Y, lda2, mem, n, 1.0_rk8, A, lda)



    end subroutine





! real single precision


! complex double precision







    subroutine bandred_complex_double(na, a, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols, tmat, wantDebug, &
                             useGPU, success)

!-------------------------------------------------------------------------------
!  bandred_complex: Reduces a distributed hermitian matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,matrixCols)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a
!  matrixCols  local columns of matrix a
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

      use precision
      use cuda_functions
      use, intrinsic :: iso_c_binding

      implicit none

      logical, intent(in)                         :: useGPU

      integer(kind=ik)                            :: na, lda, nblk, nbw, matrixCols, numBlocks, mpi_comm_rows, mpi_comm_cols

      complex(kind=ck8)              :: a(lda,*), tmat(nbw,nbw,*)



      complex(kind=ck8), parameter   :: CZERO = (0.0_rk8, 0.0_rk8), CONE = (1.0_rk8, 0.0_rk8)


      integer(kind=ik)                            :: my_prow, my_pcol, np_rows, np_cols, mpierr
      integer(kind=ik)                            :: l_cols, l_rows
      integer(kind=ik)                            :: i, j, lcs, lce, lre, lc, lr, cur_pcol, n_cols, nrow
      integer(kind=ik)                            :: istep, ncol, lch, lcx, nlc
      integer(kind=ik)                            :: tile_size, l_rows_tile, l_cols_tile

      real(kind=rk8)                    :: vnorm2
      complex(kind=ck8)              :: xf, aux1(nbw), aux2(nbw), vrl, tau, vav(nbw,nbw)

      complex(kind=ck8), allocatable :: tmp(:,:), vr(:), vmr(:,:), umc(:,:)
      integer(kind=c_intptr_t)                    :: umc_dev, tmat_dev,vav_dev,vmr_dev,a_dev
      integer(kind=ik)                            :: cur_l_rows, cur_l_cols,vmr_size ,umc_size
      integer(kind=c_size_t)                      :: lc_start, lc_end, lr_end, lce_1, lcs_1,lre_1
      integer(kind=ik)                            :: na_rows, na_cols

      integer(kind=ik), external                  :: numroc


      logical, intent(in)                         :: wantDebug
      logical, intent(out)                        :: success
      character(200)                              :: errorMessage
      integer(kind=ik)                            :: istat
      logical                                     :: successCUDA


      istat = 0
      a_dev = 0
      umc_dev = 0
      vmr_dev = 0
      successCUDA = .true.





      call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
      call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
      call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
      call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


      success = .true.

! Semibandwith nbw must be a multiple of blocksize nblk

      if (mod(nbw,nblk)/=0) then
        if (my_prow==0 .and. my_pcol==0) then
          if (wantDebug) then
            write(error_unit,*) 'ELPA2_bandred_complex: ERROR: nbw=',nbw,', nblk=',nblk
            write(error_unit,*) 'ELPA2_bandred_complex: ELPA2 works only for nbw==n*nblk'
          endif
          success = .false.
          return
        endif
      endif
      if (useGPU) then

        na_rows = numroc(na, nblk, my_prow, 0, np_rows)
        na_cols = numroc(na, nblk, my_pcol, 0, np_cols)



        successCUDA = cuda_malloc(tmat_dev, nbw*nbw*size_of_double_complex_datatype)

        if (.not.(successCUDA)) then
          print *, " bandred_complex: cuda malloc failed tmat_dev ", istat
          stop
        endif


        successCUDA = cuda_malloc(vav_dev, nbw*nbw*size_of_double_complex_datatype)

        if (.not.(successCUDA)) then
          print *, "bandred_complex:  cuda malloc failed vav_dev ", istat
          stop
        endif


        successCUDA = cuda_malloc(a_dev, lda*na_cols*size_of_double_complex_datatype)

        if (.not.(successCUDA)) then
          print *, "bandred_complex:  cuda malloc failed a_dev ", istat
          stop
        endif
      endif ! useGPU

! Matrix is split into tiles; work is done only for tiles on the diagonal or above

      tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
      tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

      l_rows_tile = tile_size/np_rows ! local rows of a tile
      l_cols_tile = tile_size/np_cols ! local cols of a tile

      if (useGPU) then



        if (.not.(successCUDA)) then
          print *, "bandred_complex:  cuda memcpy faild a_dev ", istat
          stop
        endif
      endif

      do istep = (na-1)/nbw, 1, -1

        n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

! Number of local columns/rows of remaining matrix
        l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
        l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

        if (useGPU) then
          cur_l_rows = max(l_rows, 1)
          cur_l_cols = max(l_cols, 1)

          vmr_size = cur_l_rows * 2 * n_cols
          umc_size = cur_l_cols * 2 * n_cols

          if ((.not. allocated(umc)) .or. (umc_size .gt. ubound(umc, dim=1))) then
            if (allocated(umc)) then
              deallocate(umc, stat=istat, errmsg=errorMessage)
              if (istat .ne. 0) then
                print *,"bandred_complex: error when allocating umc "//errorMessage
                stop
              endif
              successCUDA = cuda_free(umc_dev)
              if (.not.(successCUDA))then
                print *,"bandred_complex: error in cudaFree"
                stop
              endif
            endif

            allocate(umc(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              print *,"bandred_complex: error when allocating umc "//errorMessage
              stop
            endif

            if (max(l_cols,1) * 2*n_cols .gt. umc_size) then
              print *,"bandred_complex: umc_size ",max(l_cols,1) * 2*n_cols,umc_size
            endif

            successCUDA = cuda_malloc(umc_dev, umc_size*size_of_double_complex_datatype)

            if (.not.(successCUDA)) then
              print *, "bandred_complex:  cuda malloc failed umc_dev ", istat
              stop
            endif
          endif

          if ((.not. allocated(vmr)) .or. (vmr_size .gt. ubound(vmr, dim=1))) then
            if (allocated(vmr)) then
              deallocate(vmr, stat=istat, errmsg=errorMessage)
              if (istat .ne. 0) then
                print *,"bandred_complex: error when deallocating vmr "//errorMessage
                stop
              endif
              successCUDA = cuda_free(vmr_dev)
              if (.not.(successCUDA))then
                print *,"bandred_complex: error in cudaFree"
                stop
              endif
            endif

            allocate(vmr(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              print *,"bandred_complex: error when allocating vmr "//errorMessage
              stop
            endif

            if (max(l_rows,1) * 2*n_cols .gt. vmr_size) then
              print *,"bandred_complex: vmc_size ",max(l_rows,1) * 2*n_cols,vmr_size
            endif


            successCUDA = cuda_malloc(vmr_dev, vmr_size*size_of_double_complex_datatype)

            if (.not.(successCUDA)) then
              print *, "bandred_complex:  cuda malloc failed vmr_dev ", istat
              stop
            endif

          endif

          if ((.not. allocated(vr)) .or. (l_rows + 1 .gt. ubound(vr, dim=1))) then
            if (allocated(vr)) then
              deallocate(vr, stat=istat, errmsg=errorMessage)
              if (istat .ne. 0) then
                print *,"bandred_complex: error when deallocating vr "//errorMessage
                stop
              endif
            endif

            allocate(vr(l_rows + 1), stat=istat, errmsg=errorMessage)
            if (istat .ne. 0) then
              print *,"bandred_complex: error when allocating vr "//errorMessage
              stop
            endif
          endif

        else ! GPU not used
          allocate(vmr(max(l_rows,1),2*n_cols), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_complex: error when allocating vmr "//errorMessage
            stop
          endif

          allocate(umc(max(l_cols,1),2*n_cols), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_complex: error when allocating umc "//errorMessage
            stop
          endif

          allocate(vr(l_rows+1), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_complex: error when allocating vr "//errorMessage
            stop
          endif
        endif ! useGPU


        vmr(1:l_rows,1:n_cols) = 0._ck4
        vr(:) = 0._ck4
        tmat(:,:,istep) = 0._ck4


        if (useGPU) then
          lc_start = local_index(istep*nbw+1, my_pcol, np_cols, nblk, -1)
          lc_end   = local_index(istep*nbw+n_cols, my_pcol, np_cols, nblk, -1)
          lr_end   = local_index((istep-1)*nbw + n_cols, my_prow, np_rows, nblk, -1)

          if (lc_start .le. 0) lc_start = 1
          cur_pcol = pcol(istep*nbw+1, nblk, np_cols)
          if (my_pcol == cur_pcol) then

            if (.not.(successCUDA)) then
              print *, "bandred_complex: error in cudaMemcpy2"
              stop
            endif
          endif
        endif

! Reduce current block to lower triangular form

        do lc = n_cols, 1, -1

          ncol = istep*nbw + lc ! absolute column number of householder vector
          nrow = ncol - nbw ! Absolute number of pivot row

          lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
          lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

          tau = 0

          if(nrow == 1) exit ! Nothing to do

          cur_pcol = pcol(ncol, nblk, np_cols) ! Processor column owning current block

          if (my_pcol==cur_pcol) then

! Get vector to be transformed; distribute last element and norm of
! remaining elements to all procs in current column

            vr(1:lr) = a(1:lr,lch) ! vector to be transformed

            if (my_prow==prow(nrow, nblk, np_rows)) then
              aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
              aux1(2) = vr(lr)
            else
              aux1(1) = dot_product(vr(1:lr),vr(1:lr))

              aux1(2) = 0._ck8

            endif




            call mpi_allreduce(aux1, aux2, 2, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)




            vnorm2 = aux2(1)
            vrl    = aux2(2)

! Householder transformation

            call hh_transform_complex_double(vrl, vnorm2, xf, tau)

! Scale vr and store Householder vector for back transformation

            vr(1:lr) = vr(1:lr) * xf
            if (my_prow==prow(nrow, nblk, np_rows)) then
              a(1:lr-1,lch) = vr(1:lr-1)
              a(lr,lch) = vrl

              vr(lr) = 1._ck8

            else
              a(1:lr,lch) = vr(1:lr)
            endif

          endif

! Broadcast Householder vector and tau along columns

          vr(lr+1) = tau




          call MPI_Bcast(vr, lr+1, MPI_DOUBLE_COMPLEX, cur_pcol, mpi_comm_cols, mpierr)




          vmr(1:lr,lc) = vr(1:lr)
          tau = vr(lr+1)
          tmat(lc,lc,istep) = conjg(tau) ! Store tau in diagonal of tmat

! Transform remaining columns in current block with Householder vector

! Local dot product

          aux1 = 0._ck8


          nlc = 0 ! number of local columns
          do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if (lcx>0) then
              nlc = nlc+1
              aux1(nlc) = dot_product(vr(1:lr),a(1:lr,lcx))
            endif
          enddo

! Get global dot products




          if (nlc>0) call mpi_allreduce(aux1, aux2, nlc, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)


! Transform

          nlc = 0
          do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if (lcx>0) then
              nlc = nlc+1
              a(1:lr,lcx) = a(1:lr,lcx) - conjg(tau)*aux2(nlc)*vr(1:lr)
            endif
          enddo




!
!          ! Transform
!
!          nlc = 0
!          do j=1,lc-1
!            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
!            if (lcx>0) then
!              nlc = nlc+1
!              a(1:lr,lcx) = a(1:lr,lcx) - conjg(tau)*aux2(nlc)*vr(1:lr)
!            endif
!          enddo

        enddo

! Calculate scalar products of stored Householder vectors.
! This can be done in different ways, we use zherk

        if (useGPU) then
          cur_pcol = pcol(istep*nbw+1, nblk, np_cols)
          if (my_pcol == cur_pcol) then

            if (.not.(successCUDA)) then
              print *, "bandred_complex: cuda memcpy a_dev  failed ", istat
              stop
            endif
          endif
        endif

        vav = 0
        if (l_rows>0) &

           call zherk('U', 'C', n_cols, l_rows, CONE, vmr, ubound(vmr,dim=1), CZERO, vav, ubound(vav,dim=1))
        call herm_matrix_allreduce_double(n_cols,vav, nbw,nbw,mpi_comm_rows)



! Calculate triangular matrix T for block Householder Transformation

        do lc=n_cols,1,-1
          tau = tmat(lc,lc,istep)
          if (lc<n_cols) then

            call ztrmv('U', 'C', 'N', n_cols-lc, tmat(lc+1,lc+1,istep), ubound(tmat,dim=1), vav(lc+1,lc), 1)

            tmat(lc,lc+1:n_cols,istep) = -tau * conjg(vav(lc+1:n_cols,lc))
          endif
        enddo

! Transpose vmr -> vmc (stored in umc, second half)

        call elpa_transpose_vectors_complex_double  (vmr, ubound(vmr,dim=1), mpi_comm_rows, &
                                      umc(1,n_cols+1), ubound(umc,dim=1), mpi_comm_cols, &
                                      1, istep*nbw, n_cols, nblk)


! Calculate umc = A**T * vmr
! Note that the distributed A has to be transposed
! Opposed to direct tridiagonalization there is no need to use the cache locality
! of the tiles, so we can use strips of the matrix

        umc(1:l_cols,1:n_cols) = 0.0_ck8
        vmr(1:l_rows,n_cols+1:2*n_cols) = 0._ck8

        if (l_cols>0 .and. l_rows>0) then
          if (useGPU) then
            if (size(vmr,dim=1)*size(vmr,dim=2) .gt. vmr_size) then
              print *,"bandred_complex: vmr size 2 :",size(vmr,dim=1)*size(vmr,dim=2),vmr_size
              stop
            endif

            if (.not.(successCUDA)) then
              print *, "bandred_complex:  cuda memcpy vmr_dev failed ", istat
              stop
            endif
            if (size(umc,dim=1)*size(umc,dim=2) .gt. umc_size) then
              print *,"bandred_complex: umc size 2 :",size(umc,dim=1)*size(umc,dim=2),umc_size
              stop
            endif

            if (.not.(successCUDA)) then
              print *, "bandred_complex:  cuda memcpy umc_dev failed  ", istat
              stop
            endif
          endif
          do i=0,(istep*nbw-1)/tile_size

            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            if (lce<lcs) cycle

            lre = min(l_rows,(i+1)*l_rows_tile)

            if (useGPU) then

              call cublas_ZGEMM('C', 'N', lce-lcs+1, n_cols, lre, CONE, (a_dev + ((lcs-1)*lda* &
	                        size_of_double_complex_datatype)), lda, &
                                vmr_dev, cur_l_rows, CONE, (umc_dev +(lcs-1)*size_of_double_complex_datatype), cur_l_cols)

            else

              call ZGEMM('C', 'N', lce-lcs+1, n_cols, lre, CONE, a(1,lcs), ubound(a,dim=1), &
                         vmr, ubound(vmr,dim=1), CONE, umc(lcs,1), ubound(umc,dim=1))

            endif

            if (i==0) cycle
            lre = min(l_rows,i*l_rows_tile)

            if (useGPU) then
              call cublas_ZGEMM('N', 'N', lre, n_cols, lce-lcs+1, CONE, (a_dev+((lcs-1)*lda* &
	                        size_of_double_complex_datatype)),lda,  &
                                (umc_dev+(cur_l_cols * n_cols+lcs-1)*size_of_double_complex_datatype), cur_l_cols,CONE,         &
                                (vmr_dev+(cur_l_rows * n_cols)*size_of_double_complex_datatype), cur_l_rows)
            else
              call ZGEMM('N', 'N', lre, n_cols, lce-lcs+1, CONE, a(1,lcs), lda, &
                         umc(lcs,n_cols+1), ubound(umc,dim=1), CONE, vmr(1,n_cols+1), ubound(vmr,dim=1))
            endif

          enddo

          if (useGPU) then
            if (size(vmr,dim=1)*size(vmr,dim=2) .gt. vmr_size) then
              print *,"bandred_complex: vmr size 3 :",size(vmr,dim=1)*size(vmr,dim=2),vmr_size
              stop
            endif

            if (.not.(successCUDA)) then
              print *, "bandred_complex:  cuad memcpy failed vmr ", istat
              stop
            endif
            if (size(umc,dim=1)*size(umc,dim=2) .gt. umc_size) then
              print *,"bandred_complex: umc size 3 :",size(umc,dim=1)*size(umc,dim=2),umc_size
              stop
            endif

            if (.not.(successCUDA)) then
              print *, "bandred_complex:  cuad memcpy failed umc ", istat
              stop
            endif
          endif ! useGPU
        endif

! Sum up all ur(:) parts along rows and add them to the uc(:) parts
! on the processors containing the diagonal
! This is only necessary if ur has been calculated, i.e. if the
! global tile size is smaller than the global remaining matrix

        if (tile_size < istep*nbw) then

          call elpa_reduce_add_vectors_complex_double  (vmr(1,n_cols+1),ubound(vmr,dim=1),mpi_comm_rows, &
                                          umc, ubound(umc,dim=1), mpi_comm_cols, &
                                          istep*nbw, n_cols, nblk)

        endif


        if (l_cols>0) then
          allocate(tmp(l_cols,n_cols), stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_complex: error when allocating tmp "//errorMessage
            stop
          endif



          call mpi_allreduce(umc, tmp, l_cols*n_cols, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)



          umc(1:l_cols,1:n_cols) = tmp(1:l_cols,1:n_cols)
          deallocate(tmp, stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"bandred_complex: error when deallocating tmp "//errorMessage
            stop
          endif
        endif



! U = U * Tmat**T
        if (useGPU) then
          if (size(umc,dim=1)*size(umc,dim=2) .gt. umc_size) then
            print *,"bandred_complex: umc size 4 :",size(umc,dim=1)*size(umc,dim=2),umc_size
            stop
          endif

          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuad memcpy failed umc_dev ", istat
            stop
          endif

          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuad memcpy failed tmat_dev ", istat
            stop
          endif

          call  cublas_ztrmm('Right', 'Upper', 'C', 'Nonunit', l_cols, n_cols, CONE, tmat_dev, nbw, umc_dev, cur_l_cols)

        else ! not useGPU

          call ztrmm('Right', 'Upper', 'C', 'Nonunit', l_cols, n_cols, CONE, tmat(1,1,istep), ubound(tmat,dim=1), &
                     umc, ubound(umc,dim=1))

        endif

! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T
        if (useGPU) then

          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuad memcpy failed vav_dev ", istat
            stop
          endif

          call cublas_zgemm('C', 'N', n_cols, n_cols, l_cols, CONE, umc_dev, cur_l_cols, (umc_dev +( cur_l_cols *n_cols) &
                            *size_of_double_complex_datatype ), cur_l_cols, CZERO, vav_dev, nbw)

          call cublas_ztrmm('Right', 'Upper', 'C', 'Nonunit', n_cols, n_cols, CONE, tmat_dev, nbw, vav_dev, nbw)



          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuad memcpy failed vav ", istat
            stop
          endif


          call herm_matrix_allreduce_double(n_cols,vav, nbw, nbw,mpi_comm_cols)



          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuad memcpy failed vav_dev ", istat
            stop
          endif
        else

          call zgemm('C', 'N', n_cols, n_cols, l_cols, CONE, umc, ubound(umc,dim=1), umc(1,n_cols+1), &
                     ubound(umc,dim=1), CZERO, vav, ubound(vav,dim=1))
          call ztrmm('Right', 'Upper', 'C', 'Nonunit', n_cols, n_cols, CONE, tmat(1,1,istep), &
                     ubound(tmat,dim=1), vav, ubound(vav,dim=1))



          call herm_matrix_allreduce_double(n_cols,vav,nbw,nbw,mpi_comm_cols)

        endif

! U = U - 0.5 * V * VAV

        if (useGPU) then

          call cublas_zgemm('N', 'N', l_cols, n_cols, n_cols, (-0.5_rk8, 0.0_rk8), (umc_dev +  &
                            (cur_l_cols * n_cols )*size_of_double_complex_datatype), &
                            cur_l_cols, vav_dev, nbw, CONE, umc_dev, cur_l_cols)

! Transpose umc -> umr (stored in vmr, second half)

          if (size(umc,dim=1)*size(umc,dim=2) .gt. umc_size) then
            print *,"bandred_complex: umc size 5 :",size(umc,dim=1)*size(umc,dim=2),umc_size
            stop
          endif

          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuad memcpy failed umc ", istat
            stop
          endif

          call elpa_transpose_vectors_complex_double  (umc, ubound(umc,dim=1), mpi_comm_cols, &
                                                vmr(1,n_cols+1), ubound(vmr,dim=1), mpi_comm_rows, &
                                                1, istep*nbw, n_cols, nblk)

          if (size(vmr,dim=1)*size(vmr,dim=2) .gt. vmr_size) then
            print *,"bandred_complex: vmr size 4 :",size(vmr,dim=1)*size(vmr,dim=2),vmr_size
            stop
          endif

          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuda memcpy failed vav_dev", istat
            stop
          endif

          if (size(umc,dim=1)*size(umc,dim=2) .gt. umc_size) then
            print *,"bandred_complex: umc size 6 :",size(umc,dim=1)*size(umc,dim=2),umc_size
            stop
          endif

          if (.not.(successCUDA)) then
            print *, "bandred_complex:  cuda memcpy failed umc_dev ", istat
            stop
          endif
        else ! not useGPU

          call zgemm('N', 'N', l_cols, n_cols, n_cols, (-0.5_rk8, 0.0_rk8), umc(1,n_cols+1), ubound(umc,dim=1), &
                     vav, ubound(vav,dim=1), CONE, umc, ubound(umc,dim=1))

! Transpose umc -> umr (stored in vmr, second half)


          call elpa_transpose_vectors_complex_double  (umc, ubound(umc,dim=1), mpi_comm_cols, &
                                                vmr(1,n_cols+1), ubound(vmr,dim=1), mpi_comm_rows, &
                                                1, istep*nbw, n_cols, nblk)

        endif
! A = A - V*U**T - U*V**T

        do i=0,(istep*nbw-1)/tile_size
          lcs = i*l_cols_tile+1
          lce = min(l_cols,(i+1)*l_cols_tile)
          lre = min(l_rows,(i+1)*l_rows_tile)
          if (lce<lcs .or. lre<1) cycle

            if (useGPU) then
              call cublas_zgemm('N', 'C', lre, lce-lcs+1, 2*n_cols, -CONE, &
                                vmr_dev ,cur_l_rows, (umc_dev +(lcs-1)*size_of_double_complex_datatype),cur_l_cols, &
                                CONE, (a_dev + (lcs-1)*lda*size_of_double_complex_datatype),lda)
            else
              call zgemm('N', 'C', lre,lce-lcs+1, 2*n_cols, -CONE, &
                         vmr, ubound(vmr,dim=1), umc(lcs,1), ubound(umc,dim=1), &
                         CONE, a(1,lcs), lda)
            endif

          enddo

         if (.not.(useGPU)) then

           if (allocated(vr)) then
             deallocate(vr, stat=istat, errmsg=errorMessage)
             if (istat .ne. 0) then
               print *,"bandred_complex: error when deallocating vr "//errorMessage
               stop
             endif
           endif
           if (allocated(vmr)) then
             deallocate(vmr, stat=istat, errmsg=errorMessage)
             if (istat .ne. 0) then
               print *,"bandred_complex: error when deallocating vmr "//errorMessage
               stop
             endif
           endif

           if (allocated(umc)) then
             deallocate(umc, stat=istat, errmsg=errorMessage)
             if (istat .ne. 0) then
               print *,"bandred_complex: error when deallocating umc "//errorMessage
               stop
             endif
           endif


         endif ! not useGPU

       enddo ! istep

       if (useGPU) then



         if (.not.(successCUDA)) then
           print *, "bandred_complex:  cuad memcpy failed a ", istat
           stop
         endif

         successCUDA = cuda_free(a_dev)
         if (.not.(successCUDA)) then
           print *,"bandred_complex: error in cudaFree"
           stop
         endif

         successCUDA = cuda_free(tmat_dev)
         if (.not.(successCUDA)) then
           print *,"bandred_complex: error in cudaFree"
           stop
         endif

         successCUDA = cuda_free(vav_dev)
         if (.not.(successCUDA)) then
           print *,"bandred_complex: error in cudaFree"
           stop
         endif

         if (allocated(vr)) then
           deallocate(vr, stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_complex: error when deallocating vr "//errorMessage
             stop
           endif
         endif
         if (allocated(vmr)) then
           deallocate(vmr, stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_complex: error when deallocating vmr "//errorMessage
             stop
           endif

           successCUDA = cuda_free(vmr_dev)
           if (.not.(successCUDA)) then
             print *,"bandred_complex: error in cudaFree"
             stop
           endif
         endif

         if (allocated(umc)) then
           deallocate(umc, stat=istat, errmsg=errorMessage)
           if (istat .ne. 0) then
             print *,"bandred_complex: error when deallocating umc "//errorMessage
             stop
           endif

           successCUDA = cuda_free(umc_dev)
           if (.not.(successCUDA)) then
             print *,"bandred_complex: error in cudaFree"
             stop
           endif
         endif
       endif ! use GPU




     end subroutine bandred_complex_double



     subroutine herm_matrix_allreduce_double(n,a,lda,ldb,comm)

!-------------------------------------------------------------------------------
!  herm_matrix_allreduce: Does an mpi_allreduce for a hermitian matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set


      use precision
      implicit none
      integer(kind=ik)               :: n, lda, ldb, comm
      complex(kind=ck8) :: a(lda,ldb)

      integer(kind=ik)               :: i, nc, mpierr
      complex(kind=ck8) :: h1(n*n), h2(n*n)



       nc = 0
       do i=1,n
         h1(nc+1:nc+i) = a(1:i,i)
         nc = nc+i
       enddo




       call mpi_allreduce(h1, h2, nc, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, mpierr)



       nc = 0
       do i=1,n
         a(1:i,i) = h2(nc+1:nc+i)
         a(i,1:i-1) = conjg(a(1:i-1,i))
         nc = nc+i
       enddo




!       nc = 0
!       do i=1,n
!         a(1:i,i) = h2(nc+1:nc+i)
!         a(i,1:i-1) = conjg(a(1:i-1,i))
!         nc = nc+i
!       enddo




     end subroutine herm_matrix_allreduce_double



     subroutine trans_ev_band_to_full_complex_double(na, nqc, nblk, nbw, a, lda, tmat, q, ldq, matrixCols, numBlocks, &
                                         mpi_comm_rows, mpi_comm_cols, useGPU)


!-------------------------------------------------------------------------------
!  trans_ev_band_to_full_complex:
!  Transforms the eigenvectors of a band matrix back to the eigenvectors of the original matrix
!
!  Parameters
!
!  na          Order of matrix a, number of rows of matrix q
!
!  nqc         Number of columns of matrix q
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith
!
!  a(lda,matrixCols)    Matrix containing the Householder vectors (i.e. matrix a after bandred_complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!  matrixCols  local columns of matrix a and q
!
!  tmat(nbw,nbw,numBlocks) Factors returned by bandred_complex
!
!  q           On input: Eigenvectors of band matrix
!              On output: Transformed eigenvectors
!              Distribution is like in Scalapack.
!
!  ldq         Leading dimension of q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!-------------------------------------------------------------------------------

       use cuda_functions
       use, intrinsic :: iso_c_binding
       use precision

       implicit none

       logical, intent(in)                         :: useGPU
       integer(kind=ik)                            :: na, nqc, lda, ldq, nblk, nbw, matrixCols, numBlocks, &
                                                      mpi_comm_rows, mpi_comm_cols

      complex(kind=ck8)               :: a(lda,*), q(ldq,*), tmat(nbw,nbw,*)



       complex(kind=ck8), parameter   :: CZERO = (0.0_rk8,0.0_rk8), CONE = (1.0_rk8,0.0_rk8)


       integer(kind=ik)                            :: my_prow, my_pcol, np_rows, np_cols, mpierr
       integer(kind=ik)                            :: max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
       integer(kind=ik)                            :: l_cols, l_rows, l_colh, n_cols
       integer(kind=ik)                            :: istep, lc, ncol, nrow, nb, ns

       complex(kind=ck8), allocatable :: tmp1(:), tmp2(:), hvb(:), hvm(:,:)

       integer(kind=ik)                            :: i
       integer(kind=C_intptr_T)                    :: hvm_dev, q_dev, tmat_dev, tmp_dev

       integer(kind=ik)                            :: istat
       character(200)                              :: errorMessage
       logical                                     :: successCUDA

       successCUDA = .true.




       call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
       call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
       call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
       call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


       max_blocks_row = ((na -1)/nblk)/np_rows + 1  ! Rows of A
       max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

       max_local_rows = max_blocks_row*nblk
       max_local_cols = max_blocks_col*nblk

       allocate(tmp1(max_local_cols*nbw), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_band_to_full_complex: error when allocating tmp1 "//errorMessage
         stop
       endif

       allocate(tmp2(max_local_cols*nbw), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_band_to_full_complex: error when allocating tmp2 "//errorMessage
         stop
       endif

       allocate(hvb(max_local_rows*nbw), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_band_to_full_complex: error when allocating hvb "//errorMessage
         stop
       endif

       allocate(hvm(max_local_rows,nbw), stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_band_to_full_complex: error when allocating hvm "//errorMessage
         stop
       endif

       if (useGPU) then
!   allocate(q_temp(ldq,max_local_cols), stat=istat, errmsg=errorMessage)
!   if (istat .ne. 0) then
!     print *,"trans_ev_band_to_full_complex: error when allocating q_temp "//errorMessage
!   endif

! allocate(tmat_temp(nbw,nbw), stat=istat, errmsg=errorMessage)
! if (istat .ne. 0) then
! print *,"trans_ev_band_to_full_complex: error when allocating tmat_temp "//errorMessage
! endif

         successCUDA = cuda_malloc(hvm_dev, max_local_rows*nbw*size_of_double_complex_datatype)

         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaMalloc"
           stop
         endif

         successCUDA = cuda_malloc(tmat_dev, nbw*nbw*size_of_double_complex_datatype)

         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaMalloc"
           stop
         endif

         successCUDA = cuda_malloc(q_dev, ldq*matrixCols*size_of_double_complex_datatype)

         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaMalloc"
           stop
         endif

         successCUDA = cuda_malloc(tmp_dev, max_local_cols*nbw*size_of_double_complex_datatype)

         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaMalloc"
           stop
         endif

!!e   istat = cuda_memset(tmp_dev, 0, (max_local_rows)*(nbw)*size_of_complex_datatype)
!   istat = cuda_memset(tmp_dev, 0, (max_local_cols)*(nbw)*size_of_complex_datatype)
!   if (istat .ne. 0) then
!     print *,"trans_ev_band_to_full_complex: error in cudaMalloc"
!     stop
!   endif
       endif

       hvm = 0._ck8   ! Must be set to 0 !!!
       hvb = 0._ck8   ! Safety only

       if (useGPU) then
!   q_temp(:,:) = 0.0
!   q_temp(1:ldq,1:na_cols) = q(1:ldq,1:na_cols)


         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaMemcpy"
           stop
         endif

         successCUDA = cuda_memset(hvm_dev, 0, (max_local_rows)*(nbw)*size_of_double_complex_datatype)

         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaMemset"
           stop
         endif
       endif

       l_cols = local_index(nqc, my_pcol, np_cols, nblk, -1) ! Local columns of q

       do istep=1,(na-1)/nbw

         n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

! Broadcast all Householder vectors for current step compressed in hvb

         nb = 0
         ns = 0

         do lc = 1, n_cols
           ncol = istep*nbw + lc ! absolute column number of householder vector
           nrow = ncol - nbw ! absolute number of pivot row

           l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast
           l_colh = local_index(ncol  , my_pcol, np_cols, nblk, -1) ! HV local column number

           if (my_pcol==pcol(ncol, nblk, np_cols)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

           nb = nb+l_rows

           if (lc==n_cols .or. mod(ncol,nblk)==0) then





             call MPI_Bcast(hvb(ns+1), nb-ns, MPI_DOUBLE_COMPLEX, pcol(ncol, nblk, np_cols), mpi_comm_cols, mpierr)




             ns = nb
           endif
         enddo

! Expand compressed Householder vectors into matrix hvm

         nb = 0
         do lc = 1, n_cols
           nrow = (istep-1)*nbw+lc ! absolute number of pivot row
           l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

           hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
           if (my_prow==prow(nrow, nblk, np_rows)) hvm(l_rows+1,lc) = 1.

           nb = nb+l_rows
         enddo

         if (useGPU) then

           if (.not.(successCUDA)) then
             print *,"trans_ev_band_to_full_complex: error in cudaMemcpy"
             stop
           endif
         endif

         l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)

! Q = Q - V * T**T * V**T * Q

         if (l_rows > 0) then
           if (useGPU) then

             call cublas_zgemm('C', 'N', n_cols, l_cols, l_rows, CONE, hvm_dev, max_local_rows, &
                               q_dev, ldq, CZERO, tmp_dev, n_cols)



             if (.not.(successCUDA)) then
               print *,"trans_ev_band_to_full_complex: error in cudaMemcpy"
               stop
             endif
           else

             call zgemm('C', 'N', n_cols, l_cols, l_rows, CONE, hvm, ubound(hvm,dim=1), &
                        q, ldq, CZERO, tmp1, n_cols)

           endif
         else ! l_rows > 0
           if (useGPU) then
             if (l_cols*n_cols .gt. (max_local_cols)*(nbw)) then
               print *,"trans_ev_band_to_full_complex: tmp_dev ",l_cols*n_cols,max_local_cols
               stop
             endif

!       istat = cuda_memset(tmp_dev, 0, l_cols*n_cols*size_of_complex_datatype)
!       if (istat .ne. 0) then
!         print *,"trans_ev_band_to_full_complex: error in cudaMemset"
!         stop
!       endif
           endif

           tmp1(1:l_cols*n_cols) = 0._ck8


         endif




         call mpi_allreduce(tmp1, tmp2, n_cols*l_cols, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm_rows, mpierr)



         if (l_rows>0) then

           if (useGPU) then




             if (.not.(successCUDA)) then
               print *,"trans_ev_band_to_full_complex: error in cudaMemcpy"
               stop
             endif

! tmat_temp(1:nbw,1:nbw) = tmat(1:nbw,1:nbw,istep)

             if (.not.(successCUDA)) then
               print *,"trans_ev_band_to_full_complex: error in cudaMemcpy"
               stop
             endif

             call cublas_ztrmm('L', 'U', 'C', 'N', n_cols, l_cols, CONE, tmat_dev, nbw, tmp_dev, n_cols)
             call cublas_zgemm('N', 'N', l_rows, l_cols, n_cols, -CONE, hvm_dev, max_local_rows, &
                              tmp_dev, n_cols, CONE, q_dev, ldq)

           else ! not useGPU




             call ztrmm('L', 'U', 'C', 'N', n_cols, l_cols, CONE, tmat(1,1,istep), ubound(tmat,dim=1), tmp2, n_cols)
             call zgemm('N', 'N', l_rows, l_cols, n_cols, -CONE, hvm, ubound(hvm,dim=1), &
                        tmp2, n_cols, CONE, q, ldq)



           endif
         endif

!#ifdef WITH_GPU_VERSION
!     istat =cuda_memcpy(loc(hvm(1,1)),hvm_dev,((max_local_rows)*nbw*size_of_complex_datatype),cudaMemcpyDeviceToHost)
!     if (istat .ne. 0) then
!       print *,"trans_ev_band_to_full_complex: error in cudaMemcpy"
!       stop
!     endif
!#endif

       enddo

       deallocate(tmp1, tmp2, hvb, hvm, stat=istat, errmsg=errorMessage)
       if (istat .ne. 0) then
         print *,"trans_ev_band_to_full_complex: error when deallocating tmp1, tmp2, hvb, hvm "//errorMessage
         stop
       endif

       if (useGPU) then

         successCUDA = cuda_free(hvm_dev)
         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaFree"
           stop
         endif

         successCUDA = cuda_free(tmp_dev)
         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaFree"
           stop
         endif

         successCUDA = cuda_free(tmat_dev)
         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaFree"
           stop
         endif

         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaMemcpy"
           stop
         endif
!   q(1:ldq,1:na_cols) = q_temp(1:ldq,1:na_cols)

         successCUDA = cuda_free(q_dev)
         if (.not.(successCUDA)) then
           print *,"trans_ev_band_to_full_complex: error in cudaFree"
           stop
         endif

!   deallocate(q_temp, stat=istat, errmsg=errorMessage)
!   if (istat .ne. 0) then
!     print *,"trans_ev_band_to_full_complex: error when deallocating q_temp "//errorMessage
!   endif

!deallocate(tmat_temp, stat=istat, errmsg=errorMessage)
!if (istat .ne. 0) then
!print *,"trans_ev_band_to_full_complex: error when deallocating tmat_temp "//errorMessage
!endif
       endif ! use GPU



     end subroutine trans_ev_band_to_full_complex_double



    subroutine tridiag_band_complex_double(na, nb, nblk, a, lda, d, e, matrixCols, hh_trans_complex, &
                                    mpi_comm_rows, mpi_comm_cols, mpi_comm)


!-------------------------------------------------------------------------------
! tridiag_band_complex:
! Reduces a complex hermitian symmetric band matrix to tridiagonal form
!
!  na          Order of matrix a
!
!  nb          Semi bandwith
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  a(lda,matrixCols)    Distributed system matrix reduced to banded form in the upper diagonal
!
!  lda         Leading dimension of a
!  matrixCols  local columns of matrix a
!
!  d(na)       Diagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  e(na)       Subdiagonal of tridiagonal matrix, set only on PE 0 (output)
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!  mpi_comm
!              MPI-Communicator for the total processor set
!-------------------------------------------------------------------------------

      use precision
      implicit none

!#ifdef WITH_GPU_VERSION
!   integer(C_SIZE_T)                        :: h_dev, hv_new_dev ,ab_dev,x_dev,hs_dev,tau_new_dev,hv_dev,hd_dev
!   complex*16, allocatable                  :: ab_temp(:,:)
!#endif

      integer(kind=ik), intent(in)                 ::  na, nb, nblk, lda, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm

      complex(kind=ck8),intent(in)    :: a(lda,*)

      real(kind=rk8), intent(out)        :: d(na), e(na) ! set only on PE 0
      complex(kind=ck8), intent(inout), &
          allocatable                              :: hh_trans_complex(:,:)

      real(kind=rk8)                     :: vnorm2
      complex(kind=ck8)               :: hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
      complex(kind=ck8)               :: hd(nb), hs(nb)

      integer(kind=ik)                             :: i, j, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
      integer(kind=ik)                             :: my_pe, n_pes, mpierr
      integer(kind=ik)                             :: my_prow, np_rows, my_pcol, np_cols
      integer(kind=ik)                             :: ireq_ab, ireq_hv
      integer(kind=ik)                             :: na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off

      integer(kind=ik), allocatable                :: ireq_hhr(:), ireq_hhs(:), global_id(:,:), hh_cnt(:), hh_dst(:)
      integer(kind=ik), allocatable                :: limits(:), snd_limits(:,:)
      integer(kind=ik), allocatable                :: block_limits(:)
      complex(kind=ck8), allocatable  :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:)
      integer(kind=ik)                             :: istat
      character(200)                               :: errorMessage


!   ! dummies for calling redist_band
!   real*8                   :: r_a(1,1), r_ab(1,1)



      call mpi_comm_rank(mpi_comm,my_pe,mpierr)
      call mpi_comm_size(mpi_comm,n_pes,mpierr)

      call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
      call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
      call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
      call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

!#ifdef WITH_GPU_VERSION
!   t_1 = 0
!   t_2 = 0
!#endif
! Get global_id mapping 2D procssor coordinates to global id

      allocate(global_id(0:np_rows-1,0:np_cols-1), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating global_id "//errorMessage
        stop
      endif
      global_id(:,:) = 0
      global_id(my_prow, my_pcol) = my_pe


      call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)



! Total number of blocks in the band:

      nblocks_total = (na-1)/nb + 1

! Set work distribution

      allocate(block_limits(0:n_pes), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating block_limits "//errorMessage
        stop
      endif

      call divide_band(nblocks_total, n_pes, block_limits)

! nblocks: the number of blocks for my task
      nblocks = block_limits(my_pe+1) - block_limits(my_pe)

! allocate the part of the band matrix which is needed by this PE
! The size is 1 block larger than needed to avoid extensive shifts
      allocate(ab(2*nb,(nblocks+1)*nb), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating ab "//errorMessage
        stop
      endif

!#ifdef WITH_GPU_VERSION
!   allocate(ab_temp(2*nb,nblocks*nb), stat=istat, errmsg=errorMessage)
!   if (istat .ne. 0) then
!     print *,"error when allocating ab_temp "//errorMessage
!     stop
!   endif
!#endif
      ab = 0 ! needed for lower half, the extra block should also be set to 0 for safety



!#ifdef WITH_GPU_VERSION
!
!   istat = cuda_malloc(ab_dev, 2*nb*(nblocks+1)*nb*size_of_complex_datatype)
!   if (istat .ne. 0) print *, " cuda malloc failed ab_dev", istat
!
!   istat = cuda_malloc(hv_new_dev, nb*size_of_complex_datatype )
!   if (istat .ne. 0) print *, " cuda malloc failed hv_new_dev", istat
!
!!   istat = cuda_malloc(temp_c_dev,  nb*nb*size_of_complex_datatype )
!!   if(istat .ne. 0) print *, " cuda malloc failed temp_c", istat
!
!   istat = cuda_malloc(h_dev , nb*size_of_complex_datatype)
!   if (istat .ne. 0) print *, " cuda malloc failed h_dev", istat
!
!   istat = cuda_malloc(hs_dev , nb*size_of_complex_datatype)
!   if (istat .ne. 0) print *, " cuda malloc failed hs_dev", istat
!
!   istat = cuda_malloc(x_dev , 1*size_of_complex_datatype)
!   if (istat .ne. 0) print *, " cuda malloc failed x_dev", istat
!
!   istat = cuda_malloc( tau_new_dev , 1*size_of_complex_datatype)
!   if (istat .ne. 0) print *, " cuda malloc failed tau_new_dev", istat
!
!   istat = cuda_malloc(hv_dev , nb*size_of_complex_datatype)
!   if (istat .ne. 0) print *, " cuda malloc failed hv_dev", istat
!
!   istat = cuda_malloc(hd_dev , nb*size_of_complex_datatype)
!   if (istat .ne. 0) print *, " cuda malloc failed hd_dev", istat
!#endif
! n_off: Offset of ab within band
      n_off = block_limits(my_pe)*nb

! Redistribute band in a to ab

      call redist_band_complex_double(a, lda, na, nblk, nb, matrixCols, mpi_comm_rows, mpi_comm_cols, mpi_comm, ab)

! Calculate the workload for each sweep in the back transformation
! and the space requirements to hold the HH vectors

      allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating limits "//errorMessage
        stop
      endif

      call determine_workload(na, nb, np_rows, limits)
      max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

      num_hh_vecs = 0
      num_chunks  = 0
      nx = na
      do n = 1, nblocks_total
        call determine_workload(nx, nb, np_rows, limits)
        local_size = limits(my_prow+1) - limits(my_prow)
! add to number of householder vectors
! please note: for nx==1 the one and only HH vector is 0 and is neither calculated nor send below!
        if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
          num_hh_vecs = num_hh_vecs + local_size
          num_chunks  = num_chunks+1
        endif
        nx = nx - nb
      enddo

! Allocate space for HH vectors

      allocate(hh_trans_complex(nb,num_hh_vecs), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating hh_trans_comples "//errorMessage
        stop
      endif
! Allocate and init MPI requests

      allocate(ireq_hhr(num_chunks), stat=istat, errmsg=errorMessage) ! Recv requests
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating ireq_hhr "//errorMessage
        stop
      endif

      allocate(ireq_hhs(nblocks), stat=istat, errmsg=errorMessage)    ! Send requests
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating ireq_hhs "//errorMessage
        stop
      endif

      num_hh_vecs = 0
      num_chunks  = 0
      nx = na
      nt = 0
      do n = 1, nblocks_total
        call determine_workload(nx, nb, np_rows, limits)
        local_size = limits(my_prow+1) - limits(my_prow)
        if (mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
          num_chunks  = num_chunks+1



          call mpi_irecv(hh_trans_complex(1,num_hh_vecs+1), nb*local_size, MPI_COMPLEX16, nt, &
                           10+n-block_limits(nt), mpi_comm, ireq_hhr(num_chunks), mpierr)



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
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating hh_gath "//errorMessage
        stop
      endif

      allocate(hh_send(nb,max_blk_size,nblocks), stat=istat, errmsg=errorMessage) ! send buffer for HH vectors
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating hh_sebd "//errorMessage
        stop
      endif


      hh_gath(:,:,:) = 0._ck8
      hh_send(:,:,:) = 0._ck8


! Some counters

      allocate(hh_cnt(nblocks), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating hh_cnt "//errorMessage
        stop
      endif
      allocate(hh_dst(nblocks), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating hh_dst "//errorMessage
        stop
      endif

      hh_cnt(:) = 1 ! The first transfomation vector is always 0 and not calculated at all
      hh_dst(:) = 0 ! PE number for receive

      ireq_ab = MPI_REQUEST_NULL
      ireq_hv = MPI_REQUEST_NULL

! Limits for sending

      allocate(snd_limits(0:np_rows,nblocks), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when allocating snd_limits "//errorMessage
        stop
      endif

      do iblk=1,nblocks
        call determine_workload(na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
      enddo



! ---------------------------------------------------------------------------
! Start of calculations

      na_s = block_limits(my_pe)*nb + 1

      if (my_pe>0 .and. na_s<=na) then
! send first column to previous PE
! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
        ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)



        call mpi_isend(ab_s, nb+1, MPI_COMPLEX16, my_pe-1, 1, mpi_comm, ireq_ab, mpierr)



      endif




      do istep=1,na-1

      if (my_pe==0) then
        n = MIN(na-na_s,nb) ! number of rows to be reduced

        hv(:) = 0._ck8
        tau = 0._ck8

! Transform first column of remaining matrix
! Opposed to the real case, the last step (istep=na-1) is needed here for making
! the last subdiagonal element a real number

        vnorm2 = sum(real(ab(3:n+1,na_s-n_off),kind=rk8)**2+dimag(ab(3:n+1,na_s-n_off))**2)

        if (n<2) vnorm2 = 0. ! Safety only

        call hh_transform_complex_double(ab(2,na_s-n_off),vnorm2,hf,tau)



        hv(1) = 1._ck8

        hv(2:n) = ab(3:n+1,na_s-n_off)*hf

        d(istep) = ab(1,na_s-n_off)
        e(istep) = ab(2,na_s-n_off)
        if (istep == na-1) then
          d(na) = ab(1,na_s+1-n_off)

          e(na) = 0._rk4

        endif
      else
        if (na>na_s) then
! Receive Householder vector from previous task, from PE owning subdiagonal





          call mpi_recv(hv, nb, MPI_COMPLEX16, my_pe-1, 2, mpi_comm, MPI_STATUS_IGNORE, mpierr)





          tau = hv(1)

          hv(1) = 1._ck8

        endif
      endif

      na_s = na_s+1
      if (na_s-n_off > nb) then
!#ifdef WITH_GPU_VERSION
!       ab_temp(:,1:nblocks*nb) =  ab(:,nb+1:(nblocks +1)*nb)
!       ab(:, 1:nblocks*nb) = ab_temp(:, 1:nblocks*nb)
!#else
        ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
!#endif
        ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
        n_off = n_off + nb
      endif


!#ifdef WITH_GPU_VERSION
!       call cpu_time(start)
!#endif
        do iblk=1,nblocks

          ns = na_s + (iblk-1)*nb - n_off ! first column in block
          ne = ns+nb-1                    ! last column in block

          if (ns+n_off>na) exit

! Store Householder vector for back transformation

          hh_cnt(iblk) = hh_cnt(iblk) + 1

          hh_gath(1   ,hh_cnt(iblk),iblk) = tau
          hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)



          if (hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
! Wait for last transfer to finish



            call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)


! Copy vectors into send buffer
            hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
! Send to destination



            call mpi_isend(hh_send(1,1,iblk), nb*hh_cnt(iblk), MPI_COMPLEX16, &
                            global_id(hh_dst(iblk), mod(iblk+block_limits(my_pe)-1,np_cols)), &
                            10+iblk, mpi_comm, ireq_hhs(iblk), mpierr)



! Reset counter and increase destination row
            hh_cnt(iblk) = 0
            hh_dst(iblk) = hh_dst(iblk)+1
          endif

! The following code is structured in a way to keep waiting times for
! other PEs at a minimum, especially if there is only one block.
! For this reason, it requests the last column as late as possible
! and sends the Householder vector and the first column as early
! as possible.


          nc = MIN(na-ns-n_off+1,nb) ! number of columns in diagonal block
          nr = MIN(na-nb-ns-n_off+1,nb) ! rows in subdiagonal block (may be < 0!!!)
! Note that nr>=0 implies that diagonal block is full (nc==nb)!


! Multiply diagonal block and subdiagonal block with Householder vector

          if (iblk==nblocks .and. nc==nb) then

! We need the last column from the next PE.
! First do the matrix multiplications without last column ...

! Diagonal block, the contribution of the last element is added below

            ab(1,ne) = 0._ck8
            call ZHEMV('L', nc, tau, ab(1,ns), 2*nb-1, hv, 1,(0.0_rk8,0.0_rk8),hd,1)

! Subdiagonal block
            if (nr>0) call ZGEMV('N', nr, nb-1, tau, ab(nb+1,ns), 2*nb-1, hv, 1,(0.0_rk8,0.0_rk8),hs,1)

! ... then request last column ...
! ... then request last column ...



            call mpi_recv(ab(1,ne), nb+1, MPI_COMPLEX16, my_pe+1, 1, mpi_comm, MPI_STATUS_IGNORE, mpierr)




! ... and complete the result
            hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
            hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

          else
! Normal matrix multiply
            call ZHEMV('L', nc, tau, ab(1,ns), 2*nb-1, hv, 1, (0.0_rk8,0.0_rk8), hd, 1)
            if (nr>0) call ZGEMV('N', nr, nb, tau, ab(nb+1,ns), 2*nb-1, hv, 1, (0.0_rk8,0.0_rk8), hs, 1)

          endif


! Calculate first column of subdiagonal block and calculate new
! Householder transformation for this column

            hv_new(:) = 0 ! Needed, last rows must be 0 for nr < nb
!#ifdef WITH_GPU_VERSION
!           istat = cuda_memset(hv_new_dev, 0,nb*size_of_complex_datatype )
!           if (istat .ne. 0) print *, " cuda memset failed hv_new_dev", istat
!#endif
            tau_new = 0

            if (nr>0) then

! complete (old) Householder transformation for first column

              ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

! calculate new Householder transformation ...
              if (nr>1) then

                vnorm2 = sum(real(ab(nb+2:nb+nr,ns),kind=rk8)**2+dimag(ab(nb+2:nb+nr,ns))**2)



                call hh_transform_complex_double(ab(nb+1,ns),vnorm2,hf,tau_new)



                hv_new(1) = 1._ck8
                hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
                ab(nb+2:,ns) = 0._ck8

              endif

! ... and send it away immediatly if this is the last block

              if (iblk==nblocks) then



                call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)




                hv_s(1) = tau_new
                hv_s(2:) = hv_new(2:)




                call mpi_isend(hv_s, nb, MPI_COMPLEX16, my_pe+1, 2 ,mpi_comm, ireq_hv, mpierr)




              endif

            endif


! Transform diagonal block
           x = dot_product(hv(1:nc),hd(1:nc))*conjg(tau)

           hd(1:nc) = hd(1:nc) - 0.5_rk8*x*hv(1:nc)


!#ifdef WITH_GPU_VERSION
!         istat = cuda_memcpy2d((ab_dev +  (ns-1)*2*nb*size_of_complex_datatype), 2*nb*size_of_complex_datatype,loc(a(1,ns)), 2*nb*size_of_complex_datatype, 2*size_of_complex_datatype , &
!                               2*nb*size_of_complex_datatype,cudaMemcpyHostToDevice)
!         if (istat .ne. 0) print *, "cuda memcpy a_dev H2D failed ", istat
!         istat =cuda_memcpy(hv_dev,loc(hv),nc*size_of_complex_datatype,cudaMemcpyHostToDevice)
!         if (istat .ne. 0) print *,"cuda memcpy failed hv_dev", istat
!         istat =cuda_memcpy(hd_dev,loc(hd), nb*size_of_complex_datatype,cudaMemcpyHostToDevice)
!         if (istat .ne. 0) print *,"cuda memcpy failed hd_dev", istat
!#endif

           if (my_pe>0 .and. iblk==1) then

! The first column of the diagonal block has to be send to the previous PE
! Calculate first column only ...

!#ifdef WITH_GPU_VERSION
!            call double_hh_transform_2( ns, nc, nb  )
!            istat=cuda_memcpy(loc(ab),ab_dev,(2*nb*(nblocks+1)*nb)*size_of_complex_datatype,cudaMemcpyDeviceToHost)
!            if (istat .ne. 0) print *, " cuda memcpy failed ab ", istat
!#else
             ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*conjg(hv(1)) - hv(1:nc)*conjg(hd(1))
!#endif

! ... send it away ...





             call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)



             ab_s(1:nb+1) = ab(1:nb+1,ns)



             call mpi_isend(ab_s, nb+1, MPI_COMPLEX16, my_pe-1, 1, mpi_comm, ireq_ab, mpierr)




! ... and calculate remaining columns with rank-2 update
             if (nc>1) then

!#ifdef WITH_GPU_VERSION
!              call cublas_ZHER2( 'L',nc -1,(-1.d0,0.d0), hd_dev + 1*16, 1, hv_dev +1*16, 1 , ab_dev + (ns*2*nb )*16, 2*nb-1)
!#else

               call ZHER2('L', nc-1, (-1.0_rk8,0.0_rk8), hd(2), 1, hv(2), 1, ab(1,ns+1), 2*nb-1)

!#endif
             endif
           else

! No need to  send, just a rank-2 update
!#ifdef WITH_GPU_VERSION
!            call cublas_ZHER2( 'L',nc ,(-1.d0,0.d0), hd_dev, 1, hv_dev, 1 , ab_dev + ((ns-1)*2*nb )*16, 2*nb-1)
!#else

             call ZHER2('L', nc, (-1.0_rk8,0.0_rk8), hd, 1, hv, 1, ab(1,ns), 2*nb-1)

!#endif
           endif

!#ifdef WITH_GPU_VERSION
!          istat=cuda_memcpy( loc(hd),hd_dev,nb*size_of_complex_datatype,cudaMemcpyDeviceToHost)
!          if (istat .ne. 0) print *,"cuda memcpy failed hd_dev", istat
!#endif

! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

!#ifdef WITH_GPU_VERSION
!         istat =cuda_memcpy(hs_dev,loc(hs),nb*size_of_complex_datatype,cudaMemcpyHostToDevice)
!         if (istat .ne. 0) print *,"cuda memcpy failed hs_dev", istat
!#endif

           if (nr>0) then
             if (nr>1) then
!#ifdef WITH_GPU_VERSION
!              istat = cuda_memcpy(hv_new_dev,loc(hv_new),nb*size_of_complex_datatype,cudaMemcpyHostToDevice)
!              if (istat .ne. 0) print *,"cuda memcpy failed hv_new_dev", istat
!
!              istat = cuda_memcpy(h_dev,loc(h),nb*size_of_complex_datatype,cudaMemcpyHostToDevice)
!              if (istat .ne. 0) print *,"cuda memcpy failed h_dev", istat
!
!              call cublas_ZGEMV('C',nr,nb-1,tau_new,ab_dev + (nb-1 + ns *2*nb)*16,2*nb-1,hv_new_dev,1,(0.d0,0.d0),h_dev + 1* 16,1)
!
!              istat = cuda_memcpy(tau_new_dev,loc(tau_new),1*size_of_complex_datatype,cudaMemcpyHostToDevice)
!              if (istat .ne. 0) print *,"cuda memcpy failed tau_new_dev", istat
!
!              call dot_product_kernel(nr , tau_new)
!              call dot_product_kernel_1( nb, nr , ns)
!
!              istat = cuda_memcpy(loc(x),x_dev,1*size_of_complex_datatype,cudaMemcpyDeviceToHost)
!              if (istat .ne. 0) print *, " cuda memcpy failed x_dev ", istat
!
!              istat =cuda_memcpy(loc(h),h_dev,nb*size_of_complex_datatype,cudaMemcpyDeviceToHost)
!              if (istat .ne. 0) print *, " cuda memcpy failed h ", istat
!#else

               call ZGEMV('C', nr, nb-1, tau_new, ab(nb,ns+1), 2*nb-1, hv_new, 1, (0.0_rk8, 0.0_rk8), h(2), 1)

               x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
               h(2:nb) = h(2:nb) - x*hv(2:nb)
! Unfortunately there is no BLAS routine like DSYR2 for a nonsymmetric rank 2 update
               do i=2,nb
                 ab(2+nb-i:1+nb+nr-i,i+ns-1) = ab(2+nb-i:1+nb+nr-i,i+ns-1) - hv_new(1:nr)*conjg(h(i)) - hs(1:nr)*conjg(hv(i))
               enddo
!#endif
             else
! No double Householder transformation for nr=1, just complete the row
!#ifdef WITH_GPU_VERSION
!               call double_hh_transform_1(nb, ns)
!#else
               do i=2,nb
                 ab(2+nb-i,i+ns-1) = ab(2+nb-i,i+ns-1) - hs(1)*conjg(hv(i))
               enddo
!#endif

             endif
           endif

! Use new HH vector for the next block
           hv(:) = hv_new(:)
           tau = tau_new

         enddo
!#ifdef WITH_GPU_VERSION
!      call cpu_time(finish)
!      tstep2 = finish-start
!      t_2 = t_2 + tstep2
!#endif



     enddo
!#ifdef WITH_GPU_VERSION
!     call cpu_time(finish_1)
!     tstep1 = finish_1-start_1
!     t_1 = t_1 + tstep1
!#endif

! Finish the last outstanding requests




     call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
     call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

     call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
     call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)







     call mpi_barrier(mpi_comm,mpierr)


     deallocate(ab, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"tridiag_band_complex: error when deallocating ab "//errorMessage
       stop
     endif

!#ifdef WITH_GPU_VERSION
!     deallocate(ab_temp, stat=istat, errmsg=errorMessage)
!     if (istat .ne. 0) then
!       print *,"error when deallocating ab_temp "//errorMessage
!       stop
!     endif
!
!#endif

     deallocate(ireq_hhr, ireq_hhs, stat=istat, errmsg=errorMessage)
     if (istat .ne. 0) then
       print *,"tridiag_band_complex: error when deallocating ireq_hhr, ireq_hhs "//errorMessage
       stop
     endif

      deallocate(hh_cnt, hh_dst, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when deallocating hh_cnt, hh_dst "//errorMessage
        stop
      endif

      deallocate(hh_gath, hh_send, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when deallocating hh_gath, hh_send,  "//errorMessage
        stop
      endif

      deallocate(limits, snd_limits, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when deallocating limits, snd_limits  "//errorMessage
        stop
      endif

      deallocate(block_limits, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when deallocating block_limits,  "//errorMessage
        stop
      endif

      deallocate(global_id, stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"tridiag_band_complex: error when deallocating global_id,  "//errorMessage
        stop
      endif

!#ifdef WITH_GPU_VERSION
!     istat = cuda_free(ab_dev)
!     if (istat .ne. 0) then
!       print *,"error in cudaFree"
!       stop
!     endif
!
!     istat = cuda_free(hv_new_dev)
!     if (istat .ne. 0) then
!       print *,"error in cudaFree"
!       stop
!     endif
!
!     istat = cuda_free(hs_dev)
!     if (istat .ne. 0) then
!       print *,"error in cudaFree"
!       stop
!     endif
!
!     istat = cuda_free(h_dev)
!     if (istat .ne. 0) then
!       print *,"error in cudaFree"
!       stop
!     endif
!
!     istat = cuda_free(tau_new_dev)
!     if (istat .ne. 0) then
!       print *,"error in cudaFree"
!       stop
!     endif
!
!     istat = cuda_free(x_dev)
!     if (istat .ne. 0) then
!       print *,"error in cudaFree"
!       stop
!     endif
!
!#endif


!#ifdef WITH_GPU_VERSION
!   contains
!
!     subroutine dot_product_kernel(nr,tau_new)
!       implicit none
!       integer, intent(in)    :: nr
!       complex*16, intent(in) :: tau_new
!
!       call launch_dot_product_kernel( hs_dev,hv_new_dev,tau_new,x_dev,h_dev,hv_dev, nr )
!     end subroutine
!
!     subroutine dot_product_kernel_1( nb , nr , ns)
!       implicit none
!       integer, intent(in) ::  nb, nr, ns
!
!       call launch_dot_product_kernel_1(ab_dev,hs_dev, hv_new_dev,x_dev,h_dev,hv_dev,nb , nr, ns)
!     end subroutine
!
!     subroutine double_hh_transform_1( nb , ns)
!       implicit none
!       integer, intent(in) ::  nb, ns
!
!       call launch_double_hh_transform_1(ab_dev,hs_dev,hv_dev,nb , ns)
!     end subroutine
!
!     subroutine double_hh_transform_2( ns,nc, nb)
!       implicit none
!       integer, intent(in) ::  nc, ns, nb
!
!       call launch_double_hh_transform_2(ab_dev,hd_dev,hv_dev,nc , ns, nb)
!     end subroutine
!#endif

    end subroutine tridiag_band_complex_double ! has to be checked for GPU



    subroutine trans_ev_tridi_to_band_complex_double(na, nev, nblk, nbw, q, ldq, matrixCols,  &
                                              hh_trans_complex, mpi_comm_rows, mpi_comm_cols, &
                                              wantDebug, useGPU, success, THIS_COMPLEX_ELPA_KERNEL)


!-------------------------------------------------------------------------------
!  trans_ev_tridi_to_band_complex:
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
! matrixCols   local columns of matrix q
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------

      use pack_unpack_complex
      use compute_hh_trafo_complex
      use precision
      use cuda_functions
      use, intrinsic :: iso_c_binding
      implicit none

      logical, intent(in)                         :: useGPU
      integer(kind=ik), intent(in)                :: THIS_COMPLEX_ELPA_KERNEL
      integer(kind=ik), intent(in)                :: na, nev, nblk, nbw, ldq, matrixCols, mpi_comm_rows, mpi_comm_cols

      complex(kind=ck8)              :: q(ldq,*)

      complex(kind=ck8)              :: hh_trans_complex(:,:)
      integer(kind=ik)                            :: np_rows, my_prow, np_cols, my_pcol
      integer(kind=ik)                            :: tmp

      integer(kind=ik)                            :: i, j, ip, sweep, nbuf, l_nev, a_dim2
      integer(kind=ik)                            :: current_n, current_local_n, current_n_start, current_n_end
      integer(kind=ik)                            :: next_n, next_local_n, next_n_start, next_n_end
      integer(kind=ik)                            :: bottom_msg_length, top_msg_length, next_top_msg_length
      integer(kind=ik)                            :: stripe_width, last_stripe_width, stripe_count

      integer(kind=ik)                            :: num_result_blocks, num_result_buffers, num_bufs_recvd
      integer(kind=ik)                            :: a_off, current_tv_off, max_blk_size
      integer(kind=ik)                            :: mpierr, src, src_offset, dst, offset, nfact, num_blk
      logical                                     :: flag


      complex(kind=ck8), pointer     :: aIntern(:,:,:)

      complex(kind=ck8)              :: a_complex
      type(c_ptr)                                 :: aIntern_ptr
      complex(kind=ck8), allocatable :: row(:)
      complex(kind=ck8), allocatable :: row_group(:,:)


      complex(kind=ck8), allocatable :: top_border_send_buffer(:,:,:), top_border_recv_buffer(:,:,:)
      complex(kind=ck8), allocatable :: bottom_border_send_buffer(:,:,:), bottom_border_recv_buffer(:,:,:)

      integer(kind=c_intptr_t)                    :: aIntern_dev
      integer(kind=c_intptr_t)                    :: bcast_buffer_dev
      integer(kind=c_size_t)                      :: num
      integer(kind=c_size_t)                      :: dev_offset, dev_offset_1, dev_offset_2


      integer(kind=c_intptr_t)                    :: row_dev
      integer(kind=c_intptr_t)                    :: row_group_dev
      integer(kind=c_intptr_t)                    :: hh_tau_dev
      integer(kind=c_intptr_t)                    :: hh_dot_dev
      integer(kind=ik)                            :: row_group_size, unpack_idx
      integer(kind=ik)                            :: n_times

      integer(kind=ik)                            :: top, chunk, this_chunk
      complex(kind=ck8), allocatable :: result_buffer(:,:,:)
      complex(kind=ck8), allocatable :: bcast_buffer(:,:)

      integer(kind=ik)                            :: n_off
      integer(kind=ik), allocatable               :: result_send_request(:), result_recv_request(:), limits(:)
      integer(kind=ik), allocatable               :: top_send_request(:), bottom_send_request(:)
      integer(kind=ik), allocatable               :: top_recv_request(:), bottom_recv_request(:)



      integer(kind=ik), external                  :: numroc

      integer(kind=ik)                            :: na_rows, na_cols
!    real*8                                       :: ttt0, ttt1, ttt2, t2_compute_kernel, t0_compute_kernel,t1_compute_kernel, &
!                                                    t0_mpi_time, t1_mpi_time,t2_mpi_time
!    real*8                                       :: t0_cpu_code,t1_cpu_code,t2_cpu_code,t0_block_time,t1_block_time,t2_block_time,t0_cuda_memcpy
!    real*8                                       :: t0_inner_do_time, t1_inner_do_time , t2_inner_do_time,t0_outer_do_time ,t1_outer_do_time , &
!                                                    t2_outer_do_time ,t0_result_time ,t1_result_time, t2_result_time,t0_mpi_recv_time,         &
!                                                    t1_mpi_recv_time,t2_mpi_recv_time
!   real*8                                        :: t1_mpi_wait_time,t0_mpi_wait_time,t2_mpi_wait_time,t1_memcpy_time,t0_memcpy_time,t2_memcpy_time, &
!                                                    t1_mpi_irecv_time,t0_mpi_irecv_time,t2_mpi_irecv_time,t0_mpi_outer_wait_time,t1_mpi_outer_wait_time,&
!                                                    t2_mpi_outer_wait_time, time0
!   real*4                                        :: time1

! MPI send/recv tags, arbitrary

      integer(kind=ik), parameter                 :: bottom_recv_tag = 111
      integer(kind=ik), parameter                 :: top_recv_tag    = 222
      integer(kind=ik), parameter                 :: result_recv_tag = 333



! Just for measuring the kernel performance
      real(kind=c_double)                         :: kernel_time, kernel_time_recv  ! MPI_WTIME always needs double
! long integer
      integer(kind=lik)                           :: kernel_flops, kernel_flops_recv

      logical, intent(in)                         :: wantDebug
      integer(kind=ik)                            :: istat
      character(200)                              :: errorMessage
      logical                                     :: successCUDA
      logical                                     :: success

      integer(kind=C_SIZE_T)                      :: aux_int

      successCUDA = .true.





      if (useGPU) then
        n_times =0
!    n_times_1 =0
        unpack_idx = 0
        row_group_size = 0
!    time0=0
!    t0_compute_kernel=0
      endif

      kernel_time = 0.0
      kernel_flops = 0



      call MPI_Comm_rank(mpi_comm_rows, my_prow, mpierr)
      call MPI_Comm_size(mpi_comm_rows, np_rows, mpierr)
      call MPI_Comm_rank(mpi_comm_cols, my_pcol, mpierr)
      call MPI_Comm_size(mpi_comm_cols, np_cols, mpierr)

      if (useGPU) then

        na_rows = numroc(na, nblk, my_prow, 0, np_rows)
        na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

      endif

      success = .true.
      if (mod(nbw,nblk)/=0) then
        if (my_prow==0 .and. my_pcol==0) then
          if (wantDebug) then
            write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_complex: ERROR: nbw=',nbw,', nblk=',nblk
            write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_complex: band backtransform works only for nbw==n*nblk'
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
      else
! Suggested stripe width is 48 - should this be reduced for the complex case ???

        if (useGPU) then
          stripe_width = 256
        else

          stripe_width = 48 ! Must be a multiple of 2

        endif



        stripe_count = (l_nev-1)/stripe_width + 1


! Adapt stripe width so that last one doesn't get too small


        if (.not.(useGPU)) then
          stripe_width = (l_nev-1)/stripe_count + 1
        endif


        if (.not.(useGPU)) then

          stripe_width = ((stripe_width+3)/4)*4 ! Must be a multiple of 2 because of AVX/SSE memory alignment of 32 bytes
! (2 * sizeof(double complex) == 32)

        endif

        last_stripe_width = l_nev - (stripe_count-1)*stripe_width

      endif

! Determine the matrix distribution at the beginning

      allocate(limits(0:np_rows), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error when allocating limits "//errorMessage
        stop
      endif

      call determine_workload(na, nbw, np_rows, limits)

      max_blk_size = maxval(limits(1:np_rows) - limits(0:np_rows-1))

      a_dim2 = max_blk_size + nbw




      if (.not.(useGPU)) then
        aux_int = int(stripe_width*a_dim2*stripe_count*16,kind=C_SIZE_T)

        if (posix_memalign(aIntern_ptr, 64_C_SIZE_T, aux_int) /= 0) then
          print *,"trans_ev_tridi_to_band_complex: error allocating a "//errorMessage
          stop
        endif

        call c_f_pointer(aIntern_ptr, aIntern, [stripe_width,a_dim2,stripe_count] )

!        allocate(aIntern(stripe_width,a_dim2,stripe_count), stat=istat, errmsg=errorMessage)

        aIntern(:,:,:) = 0
      endif



      allocate(row(l_nev), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating row "//errorMessage
        stop
      endif

      row(:) = 0

      if (useGPU) then

        num =  (stripe_width*a_dim2*stripe_count)*size_of_double_complex_datatype
        if (na_rows * na_cols .lt. stripe_width*a_dim2*stripe_count) then
          print *,"trans_ev_tridi_to_band_complex aIntern_dev ",na_rows * na_cols, stripe_width*a_dim2*stripe_count
!      stop
        endif

        successCUDA = cuda_malloc(aIntern_dev, stripe_width*a_dim2*stripe_count*size_of_double_complex_datatype)

        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc "
          stop
        endif

        if (num .gt. na_rows * na_cols) then
          print *,"trans_ev_tridi_to_band_complex aIntern_dev 1",num, na_rows * na_cols
!      stop
        endif
        successCUDA = cuda_memset(aIntern_dev , 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMemset "
          stop
        endif

        num =  (l_nev)*size_of_double_complex_datatype

        successCUDA = cuda_malloc( row_dev,num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc "
          stop
        endif

        successCUDA = cuda_memset(row_dev , 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMemset "
          stop
        endif

! "row_group" and "row_group_dev" are needed for GPU optimizations
        allocate(row_group(l_nev, nblk), stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error allocating row_group "//errorMessage
          stop
        endif

        row_group(:, :) = 0._ck8



        num =  (l_nev*nblk)*size_of_double_complex_datatype

        successCUDA = cuda_malloc(row_group_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc "
          stop
        endif

        successCUDA = cuda_memset(row_group_dev , 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMemset "
          stop
        endif
      endif ! useGPU

! Copy q from a block cyclic distribution into a distribution with contiguous rows,
! and transpose the matrix using stripes of given stripe_width for cache blocking.

! The peculiar way it is done below is due to the fact that the last row should be
! ready first since it is the first one to start below


      do ip = np_rows-1, 0, -1
        if (my_prow == ip) then
! Receive my rows which have not yet been received
          src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
          do i=limits(ip)+1,limits(ip+1)
            src = mod((i-1)/nblk, np_rows)
            if (src < my_prow) then




              if (useGPU) then
                call unpack_and_prepare_row_group_complex_gpu_double(i - limits(ip), .false.)



                call MPI_Recv(row_group(:, row_group_size), l_nev,MPI_COMPLEX16, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)



              else



                call MPI_Recv(row, l_nev, MPI_COMPLEX16, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)



              endif






              if (.not.(useGPU)) then

                call unpack_row_complex_cpu_double(aIntern, row,i-limits(ip), stripe_count, stripe_width, last_stripe_width)

              endif



            elseif (src==my_prow) then
              src_offset = src_offset+1
              if (useGPU) then

                call unpack_and_prepare_row_group_complex_gpu_double(i - limits(ip),.false.)

                row_group(:, row_group_size) = q(src_offset, 1:l_nev)
              else
                row(:) = q(src_offset, 1:l_nev)
              endif



              if (.not.(useGPU)) then

                call unpack_row_complex_cpu_double(aIntern, row,i-limits(ip), stripe_count, stripe_width, last_stripe_width)

              endif



            endif
          enddo
! Send all rows which have not yet been send
          src_offset = 0
          do dst = 0, ip-1
            do i=limits(dst)+1,limits(dst+1)
              if(mod((i-1)/nblk, np_rows) == my_prow) then
                  src_offset = src_offset+1
                  row(:) = q(src_offset, 1:l_nev)





                  call MPI_Send(row, l_nev, MPI_COMPLEX16, dst, 0, mpi_comm_rows, mpierr)




              endif
            enddo
          enddo
        else if(my_prow < ip) then
! Send all rows going to PE ip
          src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
          do i=limits(ip)+1,limits(ip+1)
            src = mod((i-1)/nblk, np_rows)
            if (src == my_prow) then
              src_offset = src_offset+1
              row(:) = q(src_offset, 1:l_nev)




              call MPI_Send(row, l_nev, MPI_COMPLEX16, ip, 0, mpi_comm_rows, mpierr)



            endif
          enddo
! Receive all rows from PE ip
          do i=limits(my_prow)+1,limits(my_prow+1)
            src = mod((i-1)/nblk, np_rows)
            if (src == ip) then




              if (useGPU) then
                call unpack_and_prepare_row_group_complex_gpu_double(i - limits(my_prow), .false.)


                call MPI_Recv(row_group(:, row_group_size), l_nev,MPI_COMPLEX16, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)


              else


                call MPI_Recv(row, l_nev, MPI_COMPLEX16, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)


              endif






              if (.not.(useGPU)) then

                call unpack_row_complex_cpu_double(aIntern, row,i-limits(my_prow), stripe_count, stripe_width, last_stripe_width)

              endif



            endif
          enddo
        endif
      enddo

      if (useGPU) then

        call unpack_and_prepare_row_group_complex_gpu_double(-1, .true.)

        successCUDA = cuda_devicesynchronize()
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaDeviceSynchronize"
          stop
        endif
      endif

! Set up result buffer queue

      num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

      num_result_buffers = 4*nfact
      allocate(result_buffer(l_nev,nblk,num_result_buffers), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating result_buffer "//errorMessage
        stop
      endif

      allocate(result_send_request(num_result_buffers), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating result_send_request "//errorMessage
        stop
      endif

      allocate(result_recv_request(num_result_buffers), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating result_recv_request "//errorMessage
        stop
      endif


      result_send_request(:) = MPI_REQUEST_NULL
      result_recv_request(:) = MPI_REQUEST_NULL


! Queue up buffers



      if (my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
        do j = 1, min(num_result_buffers, num_result_blocks)

          call MPI_Irecv(result_buffer(1,1,j), l_nev*nblk, MPI_COMPLEX16, 0, result_recv_tag, &
                             mpi_comm_rows, result_recv_request(j), mpierr)

        enddo
      endif


      num_bufs_recvd = 0 ! No buffers received yet

! Initialize top/bottom requests

      allocate(top_send_request(stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating top_send_request "//errorMessage
        stop
      endif

      allocate(top_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating top_recv_request "//errorMessage
        stop
      endif

      allocate(bottom_send_request(stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating bottom_send_request "//errorMessage
        stop
      endif

      allocate(bottom_recv_request(stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating bottom_recv_request "//errorMessage
        stop
      endif


      top_send_request(:) = MPI_REQUEST_NULL
      top_recv_request(:) = MPI_REQUEST_NULL
      bottom_send_request(:) = MPI_REQUEST_NULL
      bottom_recv_request(:) = MPI_REQUEST_NULL



      allocate(top_border_send_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating top_border_send_buffer "//errorMessage
        stop
      endif

      allocate(top_border_recv_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating top_border_recv_buffer "//errorMessage
        stop
      endif

      allocate(bottom_border_send_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating bottom_border_send_buffer "//errorMessage
        stop
      endif

      allocate(bottom_border_recv_buffer(stripe_width, nbw, stripe_count), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating bottom_border_recv_buffer "//errorMessage
        stop
      endif


      top_border_send_buffer(:,:,:) = 0._ck8
      top_border_recv_buffer(:,:,:) = 0._ck8
      bottom_border_send_buffer(:,:,:) = 0._ck8
      bottom_border_recv_buffer(:,:,:) = 0._ck8




! Initialize broadcast buffer

      allocate(bcast_buffer(nbw, max_blk_size), stat=istat, errmsg=errorMessage)
      if (istat .ne. 0) then
        print *,"trans_ev_tridi_to_band_complex: error allocating bcast_buffer "//errorMessage
        stop
      endif
      bcast_buffer = 0

      if (useGPU) then

        num =  ( nbw * max_blk_size) * size_of_double_complex_datatype

        successCUDA = cuda_malloc(bcast_buffer_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc"
          stop
        endif

        successCUDA = cuda_memset( bcast_buffer_dev, 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMemset"
          stop
        endif

        num =  ((max_blk_size-1))*size_of_double_complex_datatype

        successCUDA = cuda_malloc( hh_dot_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc"
          stop
        endif

        successCUDA = cuda_memset( hh_dot_dev, 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMemset"
          stop
        endif

        num =  (max_blk_size)*size_of_double_complex_datatype

        successCUDA = cuda_malloc( hh_tau_dev, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc"
          stop
        endif

        successCUDA = cuda_memset( hh_tau_dev, 0, num)
        if (.not.(successCUDA)) then
          print *,"trans_ev_tridi_to_band_complex: error in cudaMemset"
          stop
        endif
      endif ! useGPU

      current_tv_off = 0 ! Offset of next row to be broadcast


! ------------------- start of work loop -------------------

      a_off = 0 ! offset in A (to avoid unnecessary shifts)

      top_msg_length = 0
      bottom_msg_length = 0



      do sweep = 0, (na-1)/nbw



        current_n = na - sweep*nbw
        call determine_workload(current_n, nbw, np_rows, limits)
        current_n_start = limits(my_prow)
        current_n_end   = limits(my_prow+1)
        current_local_n = current_n_end - current_n_start

        next_n = max(current_n - nbw, 0)
        call determine_workload(next_n, nbw, np_rows, limits)
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
          do i = 1, stripe_count






            call MPI_Irecv(bottom_border_recv_buffer(1,1,i), nbw*stripe_width, MPI_COMPLEX16, my_prow+1, bottom_recv_tag, &
                           mpi_comm_rows, bottom_recv_request(i), mpierr)







          enddo
        endif


        if (current_local_n > 1) then
          if (my_pcol == mod(sweep,np_cols)) then
            bcast_buffer(:,1:current_local_n) = hh_trans_complex(:,current_tv_off+1:current_tv_off+current_local_n)
            current_tv_off = current_tv_off + current_local_n
          endif




          call mpi_bcast(bcast_buffer, nbw*current_local_n, MPI_COMPLEX16, mod(sweep,np_cols), mpi_comm_cols, mpierr)




          if (useGPU) then



            call extract_hh_tau_complex_gpu_double(nbw, current_local_n, .false.)
            call compute_hh_dot_products_complex_gpu_double(nbw, current_local_n)

          endif
        else
! for current_local_n == 1 the one and only HH vector is 0 and not stored in hh_trans_complex

          bcast_buffer(:,1) = 0._ck8


          if (useGPU) then

            successCUDA = cuda_memset(bcast_buffer_dev, 0, nbw * size_of_double_complex_datatype)

            if (.not.(successCUDA)) then
              print *,"trans_ev_tridi_to_band_complex: error in cudaMemset"
              stop
            endif

            call extract_hh_tau_complex_gpu_double(nbw, 1, .true.)

!NOTE(ca): I commented out the following line
!        istat =  cuda_memcpy(loc(bcast_buffer(1,1)),bcast_buffer_dev,nbw*current_local_n * size_of_complex_datatype ,
!        cudaMemcpyDeviceToHost)
!        if (istat .ne. 0) then
!          print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc"
!          stop
!        endif

          endif ! useGPU
        endif



        if (l_nev == 0) cycle

          if (current_local_n > 0) then


            do i = 1, stripe_count



!wait_b
              if (current_n_end < current_n) then






                call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)









                n_off = current_local_n+a_off
                if (useGPU) then
!                t1_memcpy_time =MPI_Wtime()

                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * size_of_double_complex_datatype


                  if (.not.(successCUDA)) then
                    print *,"trans_ev_tridi_to_band_complex: error in cudaMalloc"
                    stop
                  endif

!                t2_memcpy_time =MPI_Wtime()
!                t0_memcpy_time = t0_memcpy_time + ( t2_memcpy_time - t1_memcpy_time)
                else
                  aIntern(:,n_off+1:n_off+nbw,i) = bottom_border_recv_buffer(:,1:nbw,i)
                endif



                if (next_n_end < next_n) then






                  call MPI_Irecv(bottom_border_recv_buffer(1,1,i), nbw*stripe_width, MPI_COMPLEX16, my_prow+1, bottom_recv_tag, &
                                     mpi_comm_rows, bottom_recv_request(i), mpierr)





                endif
              endif

              if (current_local_n <= bottom_msg_length + top_msg_length) then

!wait_t
                if (top_msg_length>0) then






                  call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)


                  if (useGPU) then
!                t2_mpi_wait_time =MPI_Wtime()
!                t0_mpi_wait_time = t0_mpi_wait_time + ( t2_mpi_wait_time -t1_mpi_wait_time)
!                t1_memcpy_time =MPI_Wtime()
!

                    dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * &
                                  size_of_double_complex_datatype

!                   host_offset= (0 + (0 * stripe_width) + ( (i-1) * stripe_width * nbw ))* 16

                    if (.not.(successCUDA)) then
                      print *,"trans_ev_tridi_to_band_complex: error in cudaMemcpy"
                      stop
                    endif

!                   t2_memcpy_time =MPI_Wtime()
!                   t0_memcpy_time = t0_memcpy_time + ( t2_memcpy_time - t1_memcpy_time)
                  else
                    aIntern(:,a_off+1:a_off+top_msg_length,i) = top_border_recv_buffer(:,1:top_msg_length,i)
                  endif ! useGPU


                endif

!compute

                if (useGPU) then

                  call compute_hh_trafo_complex_gpu_double(0, current_local_n, i, a_off, dev_offset, dev_offset_1, dev_offset_2)
!                call compute_hh_trafo_complex_gpu_double(0, current_local_n, i)

                else

                  call compute_hh_trafo_complex_cpu_double(aIntern, stripe_width, a_dim2, stripe_count,            &
                                                    a_off, nbw, max_blk_size, bcast_buffer, kernel_flops, kernel_time, &
                                                    0, current_local_n, i, last_stripe_width,                          &
                                                    THIS_COMPLEX_ELPA_KERNEL)

                endif



!send_b






                call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)







                if (bottom_msg_length>0) then
                  n_off = current_local_n+nbw-bottom_msg_length+a_off



                  if (useGPU) then
!                  t1_memcpy_time =MPI_Wtime()

                    dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * &
                                 size_of_double_complex_datatype


                    if (.not.(successCUDA)) then
                      print *,"trans_ev_tridi_to_band_complex: error in cudaMemcpy"
                      stop
                    endif

!                   t2_memcpy_time =MPI_Wtime()
!                   t0_memcpy_time = t0_memcpy_time + ( t2_memcpy_time -t1_memcpy_time)
                  else
                    bottom_border_send_buffer(:,1:bottom_msg_length,i) = aIntern(:,n_off+1:n_off+bottom_msg_length,i)
                  endif



                  call MPI_Isend(bottom_border_send_buffer(1,1,i), bottom_msg_length*stripe_width, MPI_COMPLEX16, my_prow+1, &
                              top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)





                endif

              else

!compute



              if (useGPU) then
                call compute_hh_trafo_complex_gpu_double(current_local_n -bottom_msg_length, bottom_msg_length, i, a_off, &
                                                  dev_offset, dev_offset_1, dev_offset_2)
!              call compute_hh_trafo_complex_gpu_double(current_local_n -bottom_msg_length, bottom_msg_length, i)
              else
                call compute_hh_trafo_complex_cpu_double(aIntern, stripe_width, a_dim2, stripe_count,              &
                                                  a_off, nbw, max_blk_size, bcast_buffer, kernel_flops, kernel_time, &
                                                  current_local_n - bottom_msg_length, bottom_msg_length, i,         &
                                                  last_stripe_width, THIS_COMPLEX_ELPA_KERNEL)

              endif




!send_b







              call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)






              if (bottom_msg_length > 0) then
                n_off = current_local_n+nbw-bottom_msg_length+a_off


                if (useGPU) then
!                t1_memcpy_time =MPI_Wtime()

                  dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * &
                                size_of_double_complex_datatype


                  if (.not.(successCUDA)) then
                    print *,"trans_ev_tridi_to_band_complex: error in cudaMemcpy"
                    stop
                  endif

!                t2_memcpy_time =MPI_Wtime()
!                t0_memcpy_time = t0_memcpy_time + ( t2_memcpy_time -t1_memcpy_time)
                else
                  bottom_border_send_buffer(:,1:bottom_msg_length,i) = aIntern(:,n_off+1:n_off+bottom_msg_length,i)
                endif





                call MPI_Isend(bottom_border_send_buffer(1,1,i), bottom_msg_length*stripe_width, MPI_COMPLEX16, my_prow+1, &
                              top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)






              endif

!compute



              if (useGPU) then
!              call compute_hh_trafo_complex_gpu_double(top_msg_length,current_local_n-top_msg_length-bottom_msg_length, i)

                call compute_hh_trafo_complex_gpu_double(top_msg_length,current_local_n-top_msg_length-bottom_msg_length, &
                                                 i, a_off, &
                                                  dev_offset, dev_offset_1, dev_offset_2)
              else
               call compute_hh_trafo_complex_cpu_double(aIntern, stripe_width, a_dim2, stripe_count,                  &
                                                 a_off, nbw, max_blk_size, bcast_buffer, kernel_flops, kernel_time,   &
                                                 top_msg_length, current_local_n-top_msg_length-bottom_msg_length, i, &
                                                 last_stripe_width, THIS_COMPLEX_ELPA_KERNEL)

              endif


!wait_t
              if (top_msg_length>0) then






                call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)


                if (useGPU) then
!                t2_mpi_wait_time =MPI_Wtime()
!                t0_mpi_wait_time = t0_mpi_wait_time +(t2_mpi_wait_time-t1_mpi_wait_time)
!
!                t1_memcpy_time =MPI_Wtime()

                   dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * &
                                 size_of_double_complex_datatype


                   if (.not.(successCUDA)) then
                     print *,"trans_ev_tridi_to_band_complex: error in cudaMemcpy"
                     stop
                   endif

!
!                t2_memcpy_time =MPI_Wtime()
!                t0_memcpy_time = t0_memcpy_time + ( t2_memcpy_time-t1_memcpy_time)
                else
                  aIntern(:,a_off+1:a_off+top_msg_length,i) = top_border_recv_buffer(:,1:top_msg_length,i)
                endif


              endif

!compute



              if (useGPU) then
                call compute_hh_trafo_complex_gpu_double(0, top_msg_length, i, a_off, dev_offset, dev_offset_1, dev_offset_2)
!              call compute_hh_trafo_complex_gpu_double(0, top_msg_length, i)
              else
                call compute_hh_trafo_complex_cpu_double(aIntern, stripe_width, a_dim2, stripe_count,               &
                                                  a_off, nbw, max_blk_size, bcast_buffer, kernel_flops, kernel_time,  &
                                                  0, top_msg_length, i, last_stripe_width,                            &
                                                  THIS_COMPLEX_ELPA_KERNEL)
              endif




            endif

            if (next_top_msg_length > 0) then
!request top_border data






              call MPI_Irecv(top_border_recv_buffer(1,1,i), next_top_msg_length*stripe_width, MPI_COMPLEX16, my_prow-1, &
                               top_recv_tag, mpi_comm_rows, top_recv_request(i), mpierr)






            endif

!send_t
            if (my_prow > 0) then






              call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)









              if (useGPU) then
!              t1_memcpy_time =MPI_Wtime()


                dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) * &
                              size_of_double_complex_datatype


                if (.not.(successCUDA)) then
                  print *,"trans_ev_tridi_to_band_complex: error in cudaMemcpy"
                  stop
                endif

!              t2_memcpy_time =MPI_Wtime()
!              t0_memcpy_time = t0_memcpy_time + (t2_memcpy_time-t1_memcpy_time)
!
              else
                top_border_send_buffer(:,1:nbw,i) = aIntern(:,a_off+1:a_off+nbw,i)
              endif




              call MPI_Isend(top_border_send_buffer(1,1,i), nbw*stripe_width, MPI_COMPLEX16, my_prow-1, bottom_recv_tag, &
                               mpi_comm_rows, top_send_request(i), mpierr)






            endif

! Care that there are not too many outstanding top_recv_request's

            if (stripe_count > 1) then
              if (i>1) then




                call MPI_Wait(top_recv_request(i-1), MPI_STATUS_IGNORE, mpierr)




              else




                call MPI_Wait(top_recv_request(stripe_count), MPI_STATUS_IGNORE, mpierr)




              endif
            endif

          enddo



          top_msg_length = next_top_msg_length

        else
! wait for last top_send_request


          do i = 1, stripe_count




            call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)




          enddo

        endif


! Care about the result

        if (my_prow == 0) then

! topmost process sends nbw rows to destination processes

          do j=0,nfact-1

            num_blk = sweep*nfact+j ! global number of destination block, 0 based
            if (num_blk*nblk >= na) exit

            nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block





            call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)





            dst = mod(num_blk, np_rows)

            if (dst == 0) then
              if (useGPU) then
                row_group_size = min(na - num_blk*nblk, nblk)

                call pack_row_group_complex_gpu_double(row_group(:, :), j * nblk + a_off,row_group_size)

                do i = 1, row_group_size
                  q((num_blk / np_rows) * nblk + i, 1 : l_nev) = row_group(:, i)
                enddo
              else
                do i = 1, min(na - num_blk*nblk, nblk)



                  call pack_row_complex_cpu_double(aIntern, row, j*nblk+i+a_off, stripe_width, last_stripe_width, stripe_count)



                  q((num_blk/np_rows)*nblk+i,1:l_nev) = row(:)
                enddo
              endif
            else
              if (useGPU) then

                call pack_row_group_complex_gpu_double(result_buffer(:, :, nbuf), j * nblk + a_off, nblk)

              else
                do i = 1, nblk



                  call pack_row_complex_cpu_double(aIntern, result_buffer(:,i,nbuf),j*nblk+i+a_off, stripe_width, &
		                            last_stripe_width, stripe_count)



                enddo
              endif





              call MPI_Isend(result_buffer(1,1,nbuf), l_nev*nblk, MPI_COMPLEX16, dst, &
                                   result_recv_tag, mpi_comm_rows, result_send_request(nbuf), mpierr)



            endif
          enddo

        else

! receive and store final result

          do j = num_bufs_recvd, num_result_blocks-1

            nbuf = mod(j, num_result_buffers) + 1 ! buffer number to get this block

! If there is still work to do, just test for the next result request
! and leave the loop if it is not ready, otherwise wait for all
! outstanding requests

            if (next_local_n > 0) then




              call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)




              if (.not.flag) exit
            else




              call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)





            endif

! Fill result buffer into q
              num_blk = j*np_rows + my_prow ! global number of current block, 0 based
              do i = 1, min(na - num_blk*nblk, nblk)
                q(j*nblk+i, 1:l_nev) = result_buffer(1:l_nev, i, nbuf)
              enddo

! Queue result buffer again if there are outstanding blocks left


              if (j+num_result_buffers < num_result_blocks) &


               call MPI_Irecv(result_buffer(1,1,nbuf), l_nev*nblk, MPI_COMPLEX16, 0, result_recv_tag, &
                              mpi_comm_rows, result_recv_request(nbuf), mpierr)




            enddo
            num_bufs_recvd = j

          endif


! Shift the remaining rows to the front of A (if necessary)

          offset = nbw - top_msg_length

          if (offset<0) then
            if (wantDebug) then
              write(error_unit,*) 'ELPA2_trans_ev_tridi_to_band_complex: internal error, offset for shifting = ',offset
            endif
            success = .false.
            return
          endif

          a_off = a_off + offset
          if (a_off + next_local_n + nbw > a_dim2) then

            do i = 1, stripe_count
              if (useGPU) then
                chunk = min(next_local_n - 1, a_off)
                do j = top_msg_length + 1, top_msg_length + next_local_n, chunk
                  top = min(j + chunk, top_msg_length + next_local_n)
                  this_chunk = top - j + 1

                  dev_offset = (0 + ( (j-1) * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * &
                                size_of_double_complex_datatype
                  dev_offset_1 = (0 + ( (j + a_off-1) * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * &
                                 size_of_double_complex_datatype
! it is not logical to set here always the parameter "cudaMemcpyDeviceToDevice" do this ONCE at startup
!                  tmp = cuda_d2d(1)


                  if (.not.(successCUDA)) then
                    print *,"trans_ev_tridi_to_band_complex: error in cudaMemcpy"
                    stop
                  endif

                enddo
              else ! useGPU
                do j = top_msg_length+1, top_msg_length+next_local_n
                  aIntern(:,j,i) = aIntern(:,j+a_off,i)
                enddo
              endif
            enddo

            a_off = 0
          endif
        enddo



! Just for safety:

        if (ANY(top_send_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_send_request ***',my_prow,my_pcol
        if (ANY(bottom_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_send_request ***',my_prow,my_pcol
        if (ANY(top_recv_request    /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR top_recv_request ***',my_prow,my_pcol
        if (ANY(bottom_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR bottom_recv_request ***',my_prow,my_pcol


        if (my_prow == 0) then





          call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)




        endif


        if (ANY(result_send_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_send_request ***',my_prow,my_pcol
        if (ANY(result_recv_request /= MPI_REQUEST_NULL)) write(error_unit,*) '*** ERROR result_recv_request ***',my_prow,my_pcol




! deallocate all working space

        if (.not.(useGPU)) then
	  nullify(aIntern)
	  call free(aIntern_ptr)
!          deallocate(a, stat=istat, errmsg=errorMessage)
!          if (istat .ne. 0) then
!            print *,"trans_ev_tridi_to_band_complex: error deallocating a "//errorMessage
!            stop
!          endif
        endif

        deallocate(row, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating row "//errorMessage
          stop
        endif

        deallocate(limits, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating limits "//errorMessage
          stop
        endif

        deallocate(result_send_request, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating result_send_request "//errorMessage
          stop
        endif

        deallocate(result_recv_request, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating result_recv_request "//errorMessage
          stop
        endif

        deallocate(top_border_send_buffer, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating top_border_send_buffer "//errorMessage
          stop
        endif

        deallocate(top_border_recv_buffer, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating top_border_recv_buffer "//errorMessage
          stop
        endif

        deallocate(bottom_border_send_buffer, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating top_border_send_buffer "//errorMessage
          stop
        endif

        deallocate(bottom_border_recv_buffer, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating bottom_border_recv_buffer "//errorMessage
          stop
        endif

        deallocate(result_buffer, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating result_buffer "//errorMessage
          stop
        endif

        deallocate(bcast_buffer, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating bcast_buffer "//errorMessage
          stop
        endif

        deallocate(top_send_request, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating top_send_request "//errorMessage
          stop
        endif

        deallocate(top_recv_request, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating top_recv_request "//errorMessage
          stop
        endif

        deallocate(bottom_send_request, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating bottom_send_request "//errorMessage
          stop
        endif

        deallocate(bottom_recv_request, stat=istat, errmsg=errorMessage)
        if (istat .ne. 0) then
          print *,"trans_ev_tridi_to_band_complex: error deallocating bottom_recv_request "//errorMessage
          stop
        endif

        if (useGPU) then
          successCUDA = cuda_free(aIntern_dev)
          if (.not.(successCUDA)) then
            print *,"trans_ev_tridi_to_band_complex: error in cudaFree"
            stop
          endif

          successCUDA = cuda_free(hh_tau_dev)
          if (.not.(successCUDA)) then
            print *,"trans_ev_tridi_to_band_complex: error in cudaFree"
            stop
          endif

          successCUDA = cuda_free(hh_dot_dev)
          if (.not.(successCUDA)) then
            print *,"trans_ev_tridi_to_band_complex: error in cudaFree"
            stop
          endif

          successCUDA = cuda_free(row_dev)
          if (.not.(successCUDA)) then
            print *,"trans_ev_tridi_to_band_complex: error in cudaFree"
            stop
          endif

          deallocate(row_group, stat=istat, errmsg=errorMessage)
          if (istat .ne. 0) then
            print *,"trans_ev_tridi_to_band_complex: error deallocating row_group "//errorMessage
            stop
          endif

          successCUDA= cuda_free(row_group_dev)
          if (.not.(successCUDA)) then
            print *,"trans_ev_tridi_to_band_complex: error in cudaFree"
            stop
          endif

          successCUDA =  cuda_free(bcast_buffer_dev)
          if (.not.(successCUDA)) then
            print *,"trans_ev_tridi_to_band_complex: error in cudaFree"
            stop
          endif

        endif ! useGPU

       return
       contains

! The host wrapper for extracting "tau" from the HH reflectors (see the
! kernel below)
    subroutine extract_hh_tau_complex_gpu_double(nbw, n, is_zero)

      use cuda_c_kernel
      use precision
      implicit none
      integer(kind=ik), value :: nbw, n
      logical, value          :: is_zero
      integer(kind=ik)        :: val_is_zero

      if (is_zero) then
        val_is_zero = 1
      else
        val_is_zero = 0
      endif

      call launch_extract_hh_tau_c_kernel_complex_double(bcast_buffer_dev,hh_tau_dev, nbw, n,val_is_zero)

    end subroutine


    subroutine compute_hh_dot_products_complex_gpu_double(nbw, n)

      use cuda_c_kernel
      use precision
      implicit none
      integer(kind=ik), value :: nbw, n

      if (n .le. 1) return

      call launch_compute_hh_dotp_c_kernel_complex_double( bcast_buffer_dev, hh_dot_dev, nbw,n)

     end subroutine


     subroutine pack_row_group_complex_gpu_double(rows, n_offset, row_count)

       use cuda_c_kernel
       use precision
       implicit none
       integer(kind=ik), intent(in) :: n_offset, row_count
       complex(kind=ck8)             :: rows(:,:)
       integer(kind=ik)             :: max_idx
       logical                      :: successCUDA

       successCUDA = .true.

       max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

       call launch_my_pack_c_kernel_complex_double(row_count, n_offset, max_idx, stripe_width,a_dim2, stripe_count, &
                                            l_nev, aIntern_dev, row_group_dev)


       if (.not.(successCUDA)) then
         print *,"pack_row_group_complex_gpu: error in cudaMemcpy"
         stop
       endif

     end subroutine

     subroutine unpack_row_group_complex_gpu_double(rows, n_offset, row_count)

       use cuda_c_kernel
       use precision
       implicit none
       integer(kind=ik), intent(in)    :: n_offset, row_count
       complex(kind=ck8), intent(in)    :: rows(:, :)
       integer(kind=ik)                :: max_idx
       integer(kind=ik)                :: i
       logical                         :: successCUDA

       successCUDA = .true.

       max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

       if (.not.(successCUDA)) then
         print *,"unpack_row_group_complex_gpu: error in cudaMemcpy"
         stop
       endif

       call launch_my_unpack_c_kernel_complex_double( row_count, n_offset,max_idx,stripe_width,a_dim2, stripe_count, l_nev, &
                                              row_group_dev,aIntern_dev)

     end subroutine


     subroutine unpack_and_prepare_row_group_complex_gpu_double(next_unpack_idx, force)


       use precision
       implicit none
       integer(kind=ik), intent(in) :: next_unpack_idx
       logical, intent(in)          :: force

       if (row_group_size == 0) then
! Nothing to flush, just prepare for the upcoming row
         row_group_size = 1
       else
         if (force .or. (row_group_size == nblk) .or. (unpack_idx + 1 /=next_unpack_idx)) then
! A flush and a reset must  performed

           call unpack_row_group_complex_gpu_double(row_group(:, :), unpack_idx - row_group_size, row_group_size)

           row_group_size = 1
         else
! Just prepare for the upcoming row
           row_group_size = row_group_size + 1
         endif
       endif
! Always update the index for the upcoming row
       unpack_idx = next_unpack_idx

    end subroutine


    subroutine compute_hh_trafo_complex_gpu_double(off, ncols, istripe, a_off, dev_offset, dev_offset_1, dev_offset_2)


      use, intrinsic :: iso_c_binding
      use cuda_c_kernel
      use precision
      implicit none
      integer(kind=ik), intent(in) :: off, ncols, istripe
      integer(kind=ik)             :: nl
      real(kind=c_double)          :: ttt ! MPI_WTIME always needs double

      integer(kind=ik)             :: a_off
      integer(kind=c_size_t)       :: dev_offset, dev_offset_1, dev_offset_2

      if (ncols < 1) return
      ttt = mpi_wtime()
      nl = merge(stripe_width, last_stripe_width, istripe < stripe_count)


      dev_offset = (0 + ( (  a_off + off-1 )* stripe_width) + ( (istripe - 1)*stripe_width*a_dim2 )) * &
                    size_of_double_complex_datatype
      dev_offset_1 = (0 +  (  off-1 )* nbw) *size_of_double_complex_datatype
      dev_offset_2 =( off-1 )*size_of_double_complex_datatype

!      t1_compute_kernel =MPI_Wtime()
      call launch_compute_hh_trafo_c_kernel_complex_double(aIntern_dev + dev_offset,bcast_buffer_dev + dev_offset_1, &
                                                    hh_tau_dev + dev_offset_2, nl, nbw,stripe_width, off,ncols)


!      time0 = time0 + time1
!      t2_compute_kernel =MPI_Wtime()
!      t0_compute_kernel =  t0_compute_kernel + t2_compute_kernel-t1_compute_kernel

      kernel_flops = kernel_flops + 4 * int(nl, lik) * int(ncols, lik) * int(nbw, lik)
      kernel_time = kernel_time + mpi_wtime() - ttt
      n_times =n_times +1
    end subroutine
end subroutine







! complex single precision


end module ELPA2_compute
