subroutine trans_ev_band_to_full_real(na, nqc, nblk, nbw, a, lda, tmat, q, ldq, sdq, mpi_comm_rows, mpi_comm_cols)

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
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after bandred_real)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tmat(nbw,nbw,.) Factors returned by bandred_real
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

!use cublas
use cuda_routines
use iso_c_binding

   implicit none

   integer :: istat
   integer na, nqc, lda, ldq, sdq, nblk, nbw, mpi_comm_rows, mpi_comm_cols
   real*8 a(lda,*), q(ldq,sdq), tmat(nbw, nbw, *)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, l_colh, n_cols
   integer istep, lc, ncol, nrow, nb, ns

   integer(C_SIZE_T) :: hvm_dev, q_dev, tmp_dev, tmat_dev
   real*8, allocatable:: tmp1(:), tmp2(:), hvb(:), hvm(:,:)

   integer pcol, prow, i
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number


   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)


   max_blocks_row = ((na -1)/nblk)/np_rows + 1  ! Rows of A
   max_blocks_col = ((nqc-1)/nblk)/np_cols + 1  ! Columns of q!

   max_local_rows = max_blocks_row*nblk
   max_local_cols = max_blocks_col*nblk

   allocate(tmp1(max_local_cols*nbw))
   allocate(tmp2(max_local_cols*nbw))
   allocate(hvb(max_local_rows*nbw))
   allocate(hvm(max_local_rows,nbw))
   istat = cuda_malloc(hvm_dev, (max_local_rows)*nbw*8_8)
   istat = cuda_malloc(tmp_dev, (max_local_cols)*nbw*8_8)
   istat = cuda_malloc(tmat_dev, nbw*nbw*8_8)  
   istat = cuda_malloc(q_dev, ldq*sdq*8_8)

!P   q_dev = q
   istat = cuda_memcpy(q_dev, loc(q), (ldq)*(sdq)*8_8, cudaMemcpyHostToDevice)
!P   hvm_dev = 0   ! Must be set to 0 !!!
   istat = cuda_memset(hvm_dev, 0, (max_local_rows)*(nbw)*8_8)

   hvb = 0   ! Safety only
   hvm = 0

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

         if(my_pcol==pcol(ncol)) hvb(nb+1:nb+l_rows) = a(1:l_rows,l_colh)

         nb = nb+l_rows

         if(lc==n_cols .or. mod(ncol,nblk)==0) then
            call MPI_Bcast(hvb(ns+1),nb-ns,MPI_REAL8,pcol(ncol),mpi_comm_cols,mpierr)
            ns = nb
         endif
      enddo

      ! Expand compressed Householder vectors into matrix hvm

      nb = 0
      do lc = 1, n_cols
         nrow = (istep-1)*nbw+lc ! absolute number of pivot row
         l_rows = local_index(nrow-1, my_prow, np_rows, nblk, -1) ! row length for bcast

         hvm(1:l_rows,lc) = hvb(nb+1:nb+l_rows)
         if(my_prow==prow(nrow)) hvm(l_rows+1,lc) = 1.

         nb = nb+l_rows
      enddo
!P      hvm_dev = hvm
      istat = cuda_memcpy(hvm_dev, loc(hvm), ((max_local_rows)*nbw*8_8),cudaMemcpyHostToDevice)

      l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)

      ! Q = Q - V * T**T * V**T * Q

      if(l_rows>0) then
         call cublas_dgemm('T','N',n_cols,l_cols,l_rows,1.d0,hvm_dev,max_local_rows, &
                    q_dev,ldq ,0.d0,tmp_dev,n_cols)
!P      tmp1 = tmp_dev
        istat = cuda_memcpy(loc(tmp1), tmp_dev, l_cols*n_cols*8_8,cudaMemcpyDeviceToHost)
      else
         tmp1(1:l_cols*n_cols) = 0
      endif
! re-introduced this
      call mpi_allreduce(tmp1,tmp2,n_cols*l_cols,MPI_REAL8,MPI_SUM,mpi_comm_rows,mpierr)
      if(l_rows>0) then
!P      tmp_dev = tmp2
         istat = cuda_memcpy(tmp_dev, loc(tmp2), n_cols*l_cols*8_8,cudaMemcpyHostToDevice)
!P         tmat_dev(:,:) = tmat(:,:,istep)
         istat = cuda_memcpy(tmat_dev, loc(tmat(1,1,istep)), nbw*nbw*8_8,cudaMemcpyHostToDevice)
         call cublas_dtrmm('L','U','T','N',n_cols,l_cols,1.0d0, tmat_dev, nbw, tmp_dev, n_cols)
         call cublas_dgemm('N','N',l_rows,l_cols,n_cols,-1.d0,hvm_dev,max_local_rows, &
                    tmp_dev,n_cols,1.d0,q_dev,ldq)
!P     hvm = hvm_dev
       istat = cuda_memcpy(loc(hvm), hvm_dev, ((max_local_rows)*nbw*8_8),cudaMemcpyDeviceToHost)
      endif

   enddo

   deallocate(tmp1, tmp2, hvb, hvm)

   istat = cuda_free(hvm_dev)
   istat = cuda_free(tmp_dev)
   istat = cuda_free(tmat_dev)
   istat = cuda_memcpy(loc(q), q_dev, ldq*sdq*8_8, cudaMemcpyDeviceToHost)
   istat = cuda_free(q_dev)

end subroutine trans_ev_band_to_full_real
