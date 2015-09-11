!-------------------------------------------------------------------------------
subroutine bandred_complex(na, a, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols, tmat)

!-------------------------------------------------------------------------------
!  bandred_complex: Reduces a distributed hermitian matrix to band form
!
!  Parameters
!
!  na          Order of matrix
!
!  a(lda,*)    Distributed matrix which should be reduced.
!              Distribution is like in Scalapack.
!              Opposed to Scalapack, a(:,:) must be set completely (upper and lower half)
!              a(:,:) is overwritten on exit with the band and the Householder vectors
!              in the upper half.
!
!  lda         Leading dimension of a
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  nbw         semi bandwith of output matrix
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns
!
!  tmat(nbw,nbw,num_blocks)    where num_blocks = (na-1)/nbw + 1
!              Factors for the Householder vectors (returned), needed for back transformation
!
!-------------------------------------------------------------------------------
use cuda_routines
use iso_c_binding

   implicit none

   integer na, lda, nblk, nbw, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*), tmat(nbw,nbw,*)

   complex*16, parameter :: CZERO = (0.d0,0.d0), CONE = (1.d0,0.d0)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer l_cols, l_rows
   integer lcs, lce, lre
   integer i, j, lc, lr, cur_pcol, n_cols, nrow
   integer istep, ncol, lch, lcx, nlc
   integer tile_size, l_rows_tile, l_cols_tile
 
   real*8 vnorm2
   complex*16 xf, aux1(nbw), aux2(nbw), vrl, tau, vav(nbw,nbw)

   complex*16, allocatable:: tmp(:,:), vr(:), vmr(:,:), umc(:,:)

   integer(c_size_t) :: umc_dev, tmat_dev,vav_dev,vmr_dev,a_dev
   integer  cur_l_rows, cur_l_cols,vmr_size ,umc_size,istat 
   integer(c_size_t) ::  lc_start, lc_end, lr_end, lce_1, lcs_1,lre_1
   integer :: na_rows, na_cols
   integer, external :: numroc 
    
   integer pcol, prow
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   ! Semibandwith nbw must be a multiple of blocksize nblk

   if(mod(nbw,nblk)/=0) then
      if(my_prow==0 .and. my_pcol==0) then
         print *,'ERROR: nbw=',nbw,', nblk=',nblk
         print *,'ELPA2 works only for nbw==n*nblk'
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif
   endif

  na_rows = numroc(na, nblk, my_prow, 0, np_rows)
  na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

  istat = cuda_malloc(tmat_dev, nbw*nbw*16_8)
  if(istat .ne. 0) print *, " cuda malloc failed tmat_dev ", istat
  istat = cuda_malloc(vav_dev, nbw*nbw*16_8)
  if(istat .ne. 0) print *, " cuda malloc failed vav_dev ", istat
  istat = cuda_malloc(a_dev, lda*na_cols*16_8)
  if(istat .ne. 0) print *, " cuda malloc failed a_dev ", istat
 
   ! Matrix is split into tiles; work is done only for tiles on the diagonal or above

   tile_size = nblk*least_common_multiple(np_rows,np_cols) ! minimum global tile size
!  tile_size = ((128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide

   tile_size = (( 128*max(np_rows,np_cols)-1)/tile_size+1)*tile_size ! make local tiles at least 128 wide


   l_rows_tile = tile_size/np_rows ! local rows of a tile
   l_cols_tile = tile_size/np_cols ! local cols of a tile


   istat = cuda_memcpy(a_dev, loc(a(1,1)),(lda)*(na_cols)*16_8,cudaMemcpyHostToDevice)
   if(istat .ne. 0) print *, " cuda memcpy faild a_dev ", istat

   do istep = (na-1)/nbw, 1, -1

      n_cols = MIN(na,(istep+1)*nbw) - istep*nbw ! Number of columns in current step

      ! Number of local columns/rows of remaining matrix
      l_cols = local_index(istep*nbw, my_pcol, np_cols, nblk, -1)
      l_rows = local_index(istep*nbw, my_prow, np_rows, nblk, -1)

      ! Allocate vmr and umc to their exact sizes so that they can be used in bcasts and reduces

      cur_l_rows = max(l_rows, 1)
      cur_l_cols = max(l_cols, 1)

      vmr_size = cur_l_rows * 2 * n_cols
      umc_size = cur_l_cols * 2 * n_cols

      if ((.not. allocated(umc)) .or. (umc_size .gt. ubound(umc, 1))) then
       if (allocated(umc)) then
              deallocate(umc)
              istat = cuda_free(umc_dev)
          endif
           allocate(umc(max(l_cols,1),2*n_cols))
         istat = cuda_malloc(umc_dev, umc_size*16_8)
          if(istat .ne. 0) then 
                print *, " cuda malloc failed umc_dev ", istat
                exit
          endif
           

      endif

      if ((.not. allocated(vmr)) .or. (vmr_size .gt. ubound(vmr, 1))) then
           if (allocated(vmr)) then
         deallocate(vmr)
         istat = cuda_free(vmr_dev)
       endif
           allocate(vmr(max(l_rows,1),2*n_cols))
          istat = cuda_malloc(vmr_dev, vmr_size*16_8)
       if(istat .ne. 0) then 
        print *, " cuda malloc failed vmr_dev ", istat
        exit
       endif

      endif

        if ((.not. allocated(vr)) .or. (l_rows + 1 .gt. ubound(vr, 1))) then
          if (allocated(vr)) deallocate(vr)
          allocate(vr(l_rows + 1))
        endif

      vmr(1:l_rows,1:n_cols) = 0.
      vr(:) = 0
      tmat(:,:,istep) = 0
      
      lc_start = local_index(istep*nbw+1, my_pcol, np_cols, nblk, -1)
      lc_end   = local_index(istep*nbw+n_cols, my_pcol, np_cols, nblk, -1)
      lr_end   = local_index((istep-1)*nbw + n_cols, my_prow, np_rows, nblk, -1)

      if(lc_start .le. 0) lc_start = 1
      cur_pcol = pcol(istep*nbw+1)
       if(my_pcol == cur_pcol) then
                istat = cuda_memcpy2d(loc(a(1, lc_start)), lda*16_8, (a_dev + ((lc_start-1) * lda*16_8)), lda*16_8, &
                             lr_end*16_8, (lc_end - lc_start+1),cudaMemcpyDeviceToHost)
               if(istat .ne. 0) print *, "data transfer error 1 ", istat
      endif

      ! Reduce current block to lower triangular form

      do lc = n_cols, 1, -1

         ncol = istep*nbw + lc ! absolute column number of householder vector
         nrow = ncol - nbw ! Absolute number of pivot row

         lr  = local_index(nrow, my_prow, np_rows, nblk, -1) ! current row length
         lch = local_index(ncol, my_pcol, np_cols, nblk, -1) ! HV local column number

         tau = 0

         if(nrow == 1) exit ! Nothing to do

         cur_pcol = pcol(ncol) ! Processor column owning current block

         if(my_pcol==cur_pcol) then

            ! Get vector to be transformed; distribute last element and norm of
            ! remaining elements to all procs in current column

            vr(1:lr) = a(1:lr,lch) ! vector to be transformed

            if(my_prow==prow(nrow)) then
               aux1(1) = dot_product(vr(1:lr-1),vr(1:lr-1))
               aux1(2) = vr(lr)
            else
               aux1(1) = dot_product(vr(1:lr),vr(1:lr))
               aux1(2) = 0.
            endif

            call mpi_allreduce(aux1,aux2,2,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

            vnorm2 = aux2(1)
            vrl    = aux2(2)

            ! Householder transformation

            call hh_transform_complex(vrl, vnorm2, xf, tau)

            ! Scale vr and store Householder vector for back transformation

            vr(1:lr) = vr(1:lr) * xf
            if(my_prow==prow(nrow)) then
               a(1:lr-1,lch) = vr(1:lr-1)
               a(lr,lch) = vrl
               vr(lr) = 1.
            else
               a(1:lr,lch) = vr(1:lr)
            endif

         endif

         ! Broadcast Householder vector and tau along columns

         vr(lr+1) = tau
         call MPI_Bcast(vr,lr+1,MPI_DOUBLE_COMPLEX,cur_pcol,mpi_comm_cols,mpierr)
         vmr(1:lr,lc) = vr(1:lr)
         tau = vr(lr+1)
         tmat(lc,lc,istep) = conjg(tau) ! Store tau in diagonal of tmat

         ! Transform remaining columns in current block with Householder vector

         ! Local dot product

         aux1 = 0

         nlc = 0 ! number of local columns
         do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if(lcx>0) then
               nlc = nlc+1
               aux1(nlc) = dot_product(vr(1:lr),a(1:lr,lcx))
            endif
         enddo

         ! Get global dot products
         if(nlc>0) call mpi_allreduce(aux1,aux2,nlc,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)

         ! Transform

         nlc = 0
         do j=1,lc-1
            lcx = local_index(istep*nbw+j, my_pcol, np_cols, nblk, 0)
            if(lcx>0) then
               nlc = nlc+1
               a(1:lr,lcx) = a(1:lr,lcx) - conjg(tau)*aux2(nlc)*vr(1:lr)
            endif
         enddo

      enddo
      ! Calculate scalar products of stored Householder vectors.
      ! This can be done in different ways, we use zherk
     cur_pcol = pcol(istep*nbw+1)
     if(my_pcol == cur_pcol) then
        istat = cuda_memcpy2d((a_dev+((lc_start-1)*lda*16_8)), lda*16_8, loc(a(1,lc_start)), &
                lda*16_8,  lr_end*16_8, (lc_end - lc_start+1),cudaMemcpyHostToDevice)
        if(istat .ne. 0) print *, "cuda memcpy a_dev  failed ", istat
     endif

      vav = 0
      if(l_rows>0) &
         call zherk('U','C',n_cols,l_rows,CONE,vmr,ubound(vmr,1),CZERO,vav,ubound(vav,1))
      call herm_matrix_allreduce(n_cols,vav,ubound(vav,1),mpi_comm_rows)

      ! Calculate triangular matrix T for block Householder Transformation

      do lc=n_cols,1,-1
         tau = tmat(lc,lc,istep)
         if(lc<n_cols) then
            call ztrmv('U','C','N',n_cols-lc,tmat(lc+1,lc+1,istep),ubound(tmat,1),vav(lc+1,lc),1)
            tmat(lc,lc+1:n_cols,istep) = -tau * conjg(vav(lc+1:n_cols,lc))
         endif
      enddo


      ! Transpose vmr -> vmc (stored in umc, second half)

      call elpa_transpose_vectors  (vmr, 2*ubound(vmr,1), mpi_comm_rows, &
                                    umc(1,n_cols+1), 2*ubound(umc,1), mpi_comm_cols, &
                                    1, 2*istep*nbw, n_cols, 2*nblk)

      ! Calculate umc = A**T * vmr
      ! Note that the distributed A has to be transposed
      ! Opposed to direct tridiagonalization there is no need to use the cache locality
      ! of the tiles, so we can use strips of the matrix

      umc(1:l_cols,1:n_cols) = 0.d0
      vmr(1:l_rows,n_cols+1:2*n_cols) = 0

      if(l_cols>0 .and. l_rows>0) then

       ! print *,' the value of vmr_size and umc_size ', vmr_size,umc_size
        istat = cuda_memcpy(vmr_dev, loc(vmr(1,1)),vmr_size*16_8,cudaMemcpyHostToDevice)
         if(istat .ne. 0) print *, " cuda memcpy vmr_dev failed ", istat
        istat = cuda_memcpy(umc_dev, loc(umc(1,1)),umc_size*16_8,cudaMemcpyHostToDevice)
         if(istat .ne. 0) print *, " cuda memcpy umc_dev failed  ", istat


         do i=0,(istep*nbw-1)/tile_size

            lcs = i*l_cols_tile+1
            lce = min(l_cols,(i+1)*l_cols_tile)
            if(lce<lcs) cycle

            lre = min(l_rows,(i+1)*l_rows_tile)
            call cublas_ZGEMM('C','N',lce-lcs+1,n_cols,lre,CONE, (a_dev + ((lcs-1)*lda*16_8)), lda, &
                      vmr_dev,cur_l_rows,CONE,(umc_dev +(lcs-1)*16_8), cur_l_cols)

            if(i==0) cycle
            lre = min(l_rows,i*l_rows_tile)
            call cublas_ZGEMM('N','N',lre,n_cols,lce-lcs+1,CONE,  (a_dev+((lcs-1)*lda*16_8)),lda,  &
                       (umc_dev+(cur_l_cols * n_cols+lcs-1)*16_8), cur_l_cols,CONE, &
                        (vmr_dev+(cur_l_rows * n_cols)*16_8), cur_l_rows)
         enddo

          istat = cuda_memcpy(loc(vmr(1,1)),vmr_dev,vmr_size*16_8,cudaMemcpyDeviceToHost)
           if(istat .ne. 0) print *, " cuad memcpy failed vmr ", istat
         istat = cuda_memcpy(loc(umc(1,1)), umc_dev,umc_size*16_8,cudaMemcpyDeviceToHost)
         if(istat .ne. 0) print *, " cuad memcpy failed umc ", istat

      
    endif
   ! Sum up all ur(:) parts along rows and add them to the uc(:) parts
      ! on the processors containing the diagonal
      ! This is only necessary if ur has been calculated, i.e. if the
      ! global tile size is smaller than the global remaining matrix

      if(tile_size < istep*nbw) then
         call elpa_reduce_add_vectors  (vmr(1,n_cols+1),2*ubound(vmr,1),mpi_comm_rows, &
                                        umc, 2*ubound(umc,1), mpi_comm_cols, &
                                        2*istep*nbw, n_cols, 2*nblk)
      endif

      if(l_cols>0) then
         allocate(tmp(l_cols,n_cols))
         call mpi_allreduce(umc,tmp,l_cols*n_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
         umc(1:l_cols,1:n_cols) = tmp(1:l_cols,1:n_cols)
         deallocate(tmp)
      endif

   

      ! U = U * Tmat**T
      istat = cuda_memcpy(umc_dev, loc(umc(1,1)),umc_size*16_8,cudaMemcpyHostToDevice)
      if(istat .ne. 0) print *, " cuad memcpy failed umc_dev ", istat

      istat = cuda_memcpy(tmat_dev,loc(tmat(1,1,istep)),nbw*nbw*16_8,cudaMemcpyHostToDevice)
      if(istat .ne. 0) print *, " cuad memcpy failed tmat_dev ", istat

      call  cublas_ztrmm('Right','Upper','C','Nonunit',l_cols,n_cols,CONE,tmat_dev,nbw,umc_dev,cur_l_cols)

      ! VAV = Tmat * V**T * A * V * Tmat**T = (U*Tmat**T)**T * V * Tmat**T
      istat = cuda_memcpy(vav_dev,loc(vav(1,1)), nbw*nbw*16_8,cudaMemcpyHostToDevice)
      if(istat .ne. 0) print *, " cuad memcpy failed vav_dev ", istat

      call cublas_zgemm('C','N',n_cols,n_cols,l_cols,CONE,umc_dev,cur_l_cols,(umc_dev +( cur_l_cols *n_cols) *16_8 ), cur_l_cols,CZERO,vav_dev, nbw)

       call cublas_ztrmm('Right','Upper','C','Nonunit',n_cols,n_cols,CONE,tmat_dev,nbw,vav_dev,nbw)
       
       istat = cuda_memcpy(loc(vav(1,1)), vav_dev,nbw*nbw*16_8,cudaMemcpyDeviceToHost)
       if(istat .ne. 0) print *, " cuad memcpy failed vav ", istat

       call herm_matrix_allreduce(n_cols,vav,ubound(vav,1),mpi_comm_cols)

       istat = cuda_memcpy(vav_dev,loc(vav(1,1)),nbw*nbw*16_8,cudaMemcpyHostToDevice)
       if(istat .ne. 0) print *, " cuad memcpy failed vav_dev ", istat

      ! U = U - 0.5 * V * VAV
       call cublas_zgemm('N','N',l_cols,n_cols,n_cols,(-0.5d0,0.d0),(umc_dev + (cur_l_cols * n_cols )*16_8),cur_l_cols,vav_dev,nbw,CONE,umc_dev,cur_l_cols)
      ! Transpose umc -> umr (stored in vmr, second half)
       
       istat = cuda_memcpy(loc(umc(1,1)),umc_dev,umc_size*16_8,cudaMemcpyDeviceToHost)
       if(istat .ne. 0) print *, " cuad memcpy failed umc ", istat

       call elpa_transpose_vectors  (umc, 2*ubound(umc,1), mpi_comm_cols, &
                                     vmr(1,n_cols+1), 2*ubound(vmr,1), mpi_comm_rows, &
                                     1, 2*istep*nbw, n_cols, 2*nblk)



      istat = cuda_memcpy(vmr_dev,loc(vmr(1,1)),vmr_size*16_8,cudaMemcpyHostToDevice)
      if(istat .ne. 0) print *, " cuda memcpy failed vav_dev", istat
      
      istat = cuda_memcpy(umc_dev,loc(umc(1,1)),umc_size*16_8,cudaMemcpyHostToDevice)
      if(istat .ne. 0) print *, " cuda memcpy failed umc_dev ", istat
      ! A = A - V*U**T - U*V**T
      do i=0,(istep*nbw-1)/tile_size
         lcs = i*l_cols_tile+1
         lce = min(l_cols,(i+1)*l_cols_tile)
         lre = min(l_rows,(i+1)*l_rows_tile)
         if(lce<lcs .or. lre<1) cycle

         call cublas_zgemm('N','C',lre,lce-lcs+1,2*n_cols,-CONE, &
                    vmr_dev ,cur_l_rows,(umc_dev +(lcs-1)*16_8),cur_l_cols, &
                   CONE,(a_dev + (lcs-1)*lda*16_8),lda)

      enddo

   enddo

   istat = cuda_memcpy ( loc(a(1,1)), a_dev, lda*na_cols*16_8,cudaMemcpyDeviceToHost)
   if(istat .ne. 0) print *, " cuad memcpy failed a ", istat

   istat = cuda_free(a_dev)
   istat = cuda_free(tmat_dev)
   istat = cuda_free(vav_dev)

 ! Free used memory
   if (allocated(vr)) deallocate(vr)

   if (allocated(vmr)) then
        deallocate(vmr)
        istat = cuda_free(vmr_dev)
   endif


   if (allocated(umc)) then
        deallocate(umc)
        istat = cuda_free(umc_dev)
   endif

end subroutine bandred_complex

!-------------------------------------------------------------------------------

subroutine herm_matrix_allreduce(n,a,lda,comm)

!-------------------------------------------------------------------------------
!  herm_matrix_allreduce: Does an mpi_allreduce for a hermitian matrix A.
!  On entry, only the upper half of A needs to be set
!  On exit, the complete matrix is set
!-------------------------------------------------------------------------------

   implicit none
   integer n, lda, comm
   complex*16 a(lda,*)

   integer i, nc, mpierr
   complex*16 h1(n*n), h2(n*n)

   nc = 0
   do i=1,n
      h1(nc+1:nc+i) = a(1:i,i)
      nc = nc+i
   enddo

   call mpi_allreduce(h1,h2,nc,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,mpierr)

   nc = 0
   do i=1,n
      a(1:i,i) = h2(nc+1:nc+i)
      a(i,1:i-1) = conjg(a(1:i-1,i))
      nc = nc+i
   enddo

end subroutine herm_matrix_allreduce

!-------------------------------------------------------------------------------
