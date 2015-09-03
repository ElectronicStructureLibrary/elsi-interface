
!-------------------------------------------------------------------------------

subroutine trans_ev_band_to_full_complex(na, nqc, nblk, nbw, a, lda, tmat, q, ldq, sdq, mpi_comm_rows, mpi_comm_cols)

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
!  a(lda,*)    Matrix containing the Householder vectors (i.e. matrix a after bandred_complex)
!              Distribution is like in Scalapack.
!
!  lda         Leading dimension of a
!
!  tmat(nbw,nbw,.) Factors returned by bandred_complex
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
use cuda_routines
use iso_c_binding


   implicit none

   integer :: istat
   integer na, nqc, lda, ldq, sdq, nblk, nbw, mpi_comm_rows, mpi_comm_cols
   complex*16 a(lda,*), q(ldq,sdq), tmat(nbw, nbw, *)

   complex*16, parameter :: CZERO = (0.d0,0.d0), CONE = (1.d0,0.d0)

   integer my_prow, my_pcol, np_rows, np_cols, mpierr
   integer max_blocks_row, max_blocks_col, max_local_rows, max_local_cols
   integer l_cols, l_rows, l_colh, n_cols
   integer istep, lc, ncol, nrow, nb, ns

   integer(C_SIZE_T) :: hvm_dev, q_dev, tmat_dev, tmp_dev
   complex*16, allocatable:: tmp1(:), tmp2(:), hvb(:), hvm(:,:)

   integer pcol, prow, i
   pcol(i) = MOD((i-1)/nblk,np_cols) !Processor col for global col number
   prow(i) = MOD((i-1)/nblk,np_rows) !Processor row for global row number

!call mpi_barrier(mpi_comm_rows, mpierr)
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
   istat = cuda_malloc(hvm_dev, max_local_rows*nbw*16_8)
   istat = cuda_malloc(tmp_dev, max_local_cols*nbw*16_8)
   istat = cuda_malloc(tmat_dev, nbw*nbw*16_8)
   istat = cuda_malloc(q_dev, ldq*sdq*16_8)

!P   q_dev = q
   istat = cuda_memcpy(q_dev, loc(q),ldq*sdq*16_8, cudaMemcpyHostToDevice)

   hvm = 0   ! Must be set to 0 !!!
   hvb = 0   ! Safety only

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
            call MPI_Bcast(hvb(ns+1),nb-ns,MPI_DOUBLE_COMPLEX,pcol(ncol),mpi_comm_cols,mpierr)
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
      istat = cuda_memcpy(hvm_dev, loc(hvm), ((max_local_rows)*nbw*16_8),cudaMemcpyHostToDevice)

      l_rows = local_index(MIN(na,(istep+1)*nbw), my_prow, np_rows, nblk, -1)
      
      ! Q = Q - V * T**T * V**T * Q

      if(l_rows>0) then
          call cublas_zgemm('C','N',n_cols,l_cols,l_rows,CONE,hvm_dev,max_local_rows, &
                 q_dev,ldq,CZERO,tmp_dev,n_cols)
        istat = cuda_memcpy(loc(tmp1), tmp_dev, n_cols*l_cols*16_8, &
                   cudaMemcpyDeviceToHost)
      else
          tmp1(1:l_cols*n_cols) = 0
      endif
! re-introduced this
      call mpi_allreduce(tmp1,tmp2,n_cols*l_cols,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_rows,mpierr)
      if(l_rows>0) then
        istat = cuda_memcpy(tmat_dev, loc(tmat(1,1,istep)),nbw*nbw*16_8,cudaMemcpyHostToDevice)
        istat = cuda_memcpy(tmp_dev,loc(tmp2),n_cols*l_cols*16_8, &
                  cudaMemcpyHostToDevice)

        call cublas_ztrmm('L','U','C','N',n_cols,l_cols,CONE,tmat_dev,nbw,tmp_dev,n_cols)  
        call cublas_zgemm('N','N',l_rows,l_cols,n_cols,-CONE,hvm_dev, max_local_rows, &
                    tmp_dev,n_cols,CONE,q_dev,ldq)
      endif

   enddo

   deallocate(tmp1, tmp2, hvb, hvm)
   istat = cuda_free(hvm_dev)
   istat = cuda_free(tmp_dev)
   istat = cuda_free(tmat_dev)
   istat = cuda_memcpy(loc(q), q_dev,ldq*sdq*16_8, cudaMemcpyDeviceToHost)
   istat = cuda_free(q_dev)

end subroutine trans_ev_band_to_full_complex

