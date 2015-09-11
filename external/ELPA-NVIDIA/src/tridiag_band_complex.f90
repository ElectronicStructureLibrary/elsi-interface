!-------------------------------------------------------------------------------
subroutine tridiag_band_complex(na, nb, nblk, a, lda, d, e, mpi_comm_rows, mpi_comm_cols, mpi_comm)
!-------------------------------------------------------------------------------
! tridiag_band_complex:
! Reduces a real symmetric band matrix to tridiagonal form
!
!  na          Order of matrix a
!
!  nb          Semi bandwith
!
!  nblk        blocksize of cyclic distribution, must be the same in both directions!
!
!  a(lda,*)    Distributed system matrix reduced to banded form in the upper diagonal
!
!  lda         Leading dimension of a
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

   implicit none

   integer, intent(in) ::  na, nb, nblk, lda, mpi_comm_rows, mpi_comm_cols, mpi_comm
   complex*16, intent(in) :: a(lda,*)
   real*8, intent(out) :: d(na), e(na) ! set only on PE 0

   real*8 tstep,tstep1,tstep2,t_1,t_2,t_3, start,finish,start_1,finish_1,start_2,finish_2  
   real*8 vnorm2
   complex*16 hv(nb), tau, x, h(nb), ab_s(1+nb), hv_s(nb), hv_new(nb), tau_new, hf
   complex*16 hd(nb), hs(nb)

   integer i, j, n, nc, nr, ns, ne, istep, iblk, nblocks_total, nblocks, nt
   integer my_pe, n_pes, mpierr
   integer my_prow, np_rows, my_pcol, np_cols
   integer ireq_ab, ireq_hv
   integer na_s, nx, num_hh_vecs, num_chunks, local_size, max_blk_size, n_off
   integer, allocatable :: ireq_hhr(:), ireq_hhs(:), global_id(:,:), hh_cnt(:), hh_dst(:)
   integer, allocatable :: limits(:), snd_limits(:,:)
   integer, allocatable :: block_limits(:)
   complex*16, allocatable :: ab(:,:), hh_gath(:,:,:), hh_send(:,:,:),ab_temp(:,:)
   ! dummies for calling redist_band
   real*8 :: r_a(1,1), r_ab(1,1)

   call mpi_comm_rank(mpi_comm,my_pe,mpierr)
   call mpi_comm_size(mpi_comm,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   t_1 = 0
   t_2 = 0 
   ! Get global_id mapping 2D procssor coordinates to global id

   allocate(global_id(0:np_rows-1,0:np_cols-1))
   global_id(:,:) = 0
   global_id(my_prow, my_pcol) = my_pe

   call mpi_allreduce(mpi_in_place, global_id, np_rows*np_cols, mpi_integer, mpi_sum, mpi_comm, mpierr)

   ! Total number of blocks in the band:

   nblocks_total = (na-1)/nb + 1

   ! Set work distribution

   allocate(block_limits(0:n_pes))
   call divide_band(nblocks_total, n_pes, block_limits)

   ! nblocks: the number of blocks for my task
   nblocks = block_limits(my_pe+1) - block_limits(my_pe)
   ! allocate the part of the band matrix which is needed by this PE
   ! The size is 1 block larger than needed to avoid extensive shifts
   allocate(ab(2*nb,(nblocks+1)*nb))
   allocate(ab_temp(2*nb,nblocks*nb))
   ab = 0 ! needed for lower half, the extra block should also be set to 0 for safety

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nb

   ! Redistribute band in a to ab
   call redist_band(.false., r_a, a, lda, na, nblk, nb, mpi_comm_rows, mpi_comm_cols, mpi_comm, r_ab, ab)

   ! Calculate the workload for each sweep in the back transformation
   ! and the space requirements to hold the HH vectors
   allocate(limits(0:np_rows))
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
      if(mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
         num_hh_vecs = num_hh_vecs + local_size
         num_chunks  = num_chunks+1
      endif
      nx = nx - nb
   enddo

   ! Allocate space for HH vectors

   allocate(hh_trans_complex(nb,num_hh_vecs))

   ! Allocate and init MPI requests

   allocate(ireq_hhr(num_chunks)) ! Recv requests
   allocate(ireq_hhs(nblocks))    ! Send requests
   num_hh_vecs = 0
   num_chunks  = 0
   nx = na
   nt = 0
   do n = 1, nblocks_total
      call determine_workload(nx, nb, np_rows, limits)
      local_size = limits(my_prow+1) - limits(my_prow)
      if(mod(n-1,np_cols) == my_pcol .and. local_size>0 .and. nx>1) then
         num_chunks  = num_chunks+1
         call mpi_irecv(hh_trans_complex(1,num_hh_vecs+1), nb*local_size, MPI_COMPLEX16, nt, &
                        10+n-block_limits(nt), mpi_comm, ireq_hhr(num_chunks), mpierr)
         num_hh_vecs = num_hh_vecs + local_size
      endif
      nx = nx - nb
      if(n == block_limits(nt+1)) then
         nt = nt + 1
      endif
   enddo
   ireq_hhs(:) = MPI_REQUEST_NULL

   ! Buffers for gathering/sending the HH vectors

   allocate(hh_gath(nb,max_blk_size,nblocks)) ! gathers HH vectors
   allocate(hh_send(nb,max_blk_size,nblocks)) ! send buffer for HH vectors
   hh_gath(:,:,:) = 0
   hh_send(:,:,:) = 0

   ! Some counters

   allocate(hh_cnt(nblocks))
   allocate(hh_dst(nblocks))

   hh_cnt(:) = 1 ! The first transfomation vector is always 0 and not calculated at all
   hh_dst(:) = 0 ! PE number for receive

   ireq_ab = MPI_REQUEST_NULL
   ireq_hv = MPI_REQUEST_NULL

   ! Limits for sending

   allocate(snd_limits(0:np_rows,nblocks))

   do iblk=1,nblocks
      call determine_workload(na-(iblk+block_limits(my_pe)-1)*nb, nb, np_rows, snd_limits(:,iblk))
   enddo

   ! ---------------------------------------------------------------------------
   ! Start of calculations

   na_s = block_limits(my_pe)*nb + 1

   if(my_pe>0 .and. na_s<=na) then
      ! send first column to previous PE
      ! Only the PE owning the diagonal does that (sending 1 element of the subdiagonal block also)
      ab_s(1:nb+1) = ab(1:nb+1,na_s-n_off)
      call mpi_isend(ab_s,nb+1,MPI_COMPLEX16,my_pe-1,1,mpi_comm,ireq_ab,mpierr)
   endif
    call cpu_time(start_1)
   do istep=1,na-1

      if(my_pe==0) then
         n = MIN(na-na_s,nb) ! number of rows to be reduced
         hv(:) = 0
         tau = 0
         ! Transform first column of remaining matrix
         ! Opposed to the real case, the last step (istep=na-1) is needed here for making
         ! the last subdiagonal element a real number
         vnorm2 = sum(dble(ab(3:n+1,na_s-n_off))**2+dimag(ab(3:n+1,na_s-n_off))**2)
         if(n<2) vnorm2 = 0. ! Safety only
         call hh_transform_complex(ab(2,na_s-n_off),vnorm2,hf,tau)

         hv(1) = 1
         hv(2:n) = ab(3:n+1,na_s-n_off)*hf

         d(istep) = ab(1,na_s-n_off)
         e(istep) = ab(2,na_s-n_off)
         if(istep == na-1) then
            d(na) = ab(1,na_s+1-n_off)
            e(na) = 0
         endif
      else
         if(na>na_s) then
            ! Receive Householder vector from previous task, from PE owning subdiagonal
            call mpi_recv(hv,nb,MPI_COMPLEX16,my_pe-1,2,mpi_comm,MPI_STATUS_IGNORE,mpierr)
            tau = hv(1)
            hv(1) = 1.
         endif
      endif
      na_s = na_s+1
      if(na_s-n_off > nb) then
!         ab(:,1:nblocks*nb) = ab(:,nb+1:(nblocks+1)*nb)
         ab_temp(:,1:nblocks*nb) =  ab(:,nb+1:(nblocks +1)*nb)
         ab(:, 1:nblocks*nb) = ab_temp(:, 1:nblocks*nb)
         ab(:,nblocks*nb+1:(nblocks+1)*nb) = 0
         n_off = n_off + nb
         
      endif
      call cpu_time(start)
      do iblk=1,nblocks

         ns = na_s + (iblk-1)*nb - n_off ! first column in block
         ne = ns+nb-1                    ! last column in block

         if(ns+n_off>na) exit

         ! Store Householder vector for back transformation

         hh_cnt(iblk) = hh_cnt(iblk) + 1

         hh_gath(1   ,hh_cnt(iblk),iblk) = tau
         hh_gath(2:nb,hh_cnt(iblk),iblk) = hv(2:nb)

         if(hh_cnt(iblk) == snd_limits(hh_dst(iblk)+1,iblk)-snd_limits(hh_dst(iblk),iblk)) then
            ! Wait for last transfer to finish
            call mpi_wait(ireq_hhs(iblk), MPI_STATUS_IGNORE, mpierr)
            ! Copy vectors into send buffer
            hh_send(:,1:hh_cnt(iblk),iblk) = hh_gath(:,1:hh_cnt(iblk),iblk)
            ! Send to destination
            call mpi_isend(hh_send(1,1,iblk), nb*hh_cnt(iblk), MPI_COMPLEX16, &
                           global_id(hh_dst(iblk),mod(iblk+block_limits(my_pe)-1,np_cols)), &
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

         if(iblk==nblocks .and. nc==nb) then

            ! We need the last column from the next PE.
            ! First do the matrix multiplications without last column ...

            ! Diagonal block, the contribution of the last element is added below!
            ab(1,ne) = 0
            call ZHEMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,(0.d0,0.d0),hd,1)

            ! Subdiagonal block
            if(nr>0) call ZGEMV('N',nr,nb-1,tau,ab(nb+1,ns),2*nb-1,hv,1,(0.d0,0.d0),hs,1)

            ! ... then request last column ...
            call mpi_recv(ab(1,ne),nb+1,MPI_COMPLEX16,my_pe+1,1,mpi_comm,MPI_STATUS_IGNORE,mpierr)

            ! ... and complete the result
            hs(1:nr) = hs(1:nr) + ab(2:nr+1,ne)*tau*hv(nb)
            hd(nb) = hd(nb) + ab(1,ne)*hv(nb)*tau

         else

            ! Normal matrix multiply
            call ZHEMV('L',nc,tau,ab(1,ns),2*nb-1,hv,1,(0.d0,0.d0),hd,1)
            if(nr>0) call ZGEMV('N',nr,nb,tau,ab(nb+1,ns),2*nb-1,hv,1,(0.d0,0.d0),hs,1)

         endif

         ! Calculate first column of subdiagonal block and calculate new
         ! Householder transformation for this column

         hv_new(:) = 0 ! Needed, last rows must be 0 for nr < nb
         tau_new = 0

         if(nr>0) then

            ! complete (old) Householder transformation for first column

            ab(nb+1:nb+nr,ns) = ab(nb+1:nb+nr,ns) - hs(1:nr) ! Note: hv(1) == 1

            ! calculate new Householder transformation ...
            if(nr>1) then
               vnorm2 = sum(dble(ab(nb+2:nb+nr,ns))**2+dimag(ab(nb+2:nb+nr,ns))**2)
               call hh_transform_complex(ab(nb+1,ns),vnorm2,hf,tau_new)
               hv_new(1) = 1.
               hv_new(2:nr) = ab(nb+2:nb+nr,ns)*hf
               ab(nb+2:,ns) = 0
            endif

            ! ... and send it away immediatly if this is the last block

            if(iblk==nblocks) then
               call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)
               hv_s(1) = tau_new
               hv_s(2:) = hv_new(2:)
               call mpi_isend(hv_s,nb,MPI_COMPLEX16,my_pe+1,2,mpi_comm,ireq_hv,mpierr)
            endif

         endif


         ! Transform diagonal block
         x = dot_product(hv(1:nc),hd(1:nc))*conjg(tau)
         hd(1:nc) = hd(1:nc) - 0.5*x*hv(1:nc)

         if(my_pe>0 .and. iblk==1) then

            ! The first column of the diagonal block has to be send to the previous PE
            ! Calculate first column only ...

            ab(1:nc,ns) = ab(1:nc,ns) - hd(1:nc)*conjg(hv(1)) - hv(1:nc)*conjg(hd(1))

            ! ... send it away ...

            call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
            ab_s(1:nb+1) = ab(1:nb+1,ns)
            call mpi_isend(ab_s,nb+1,MPI_COMPLEX16,my_pe-1,1,mpi_comm,ireq_ab,mpierr)

            ! ... and calculate remaining columns with rank-2 update
            if(nc>1) call ZHER2('L',nc-1,(-1.d0,0.d0),hd(2),1,hv(2),1,ab(1,ns+1),2*nb-1)
         else
            ! No need to  send, just a rank-2 update
            call ZHER2('L',nc,(-1.d0,0.d0),hd,1,hv,1,ab(1,ns),2*nb-1)
         endif

         ! Do the remaining double Householder transformation on the subdiagonal block cols 2 ... nb

         if(nr>0) then
            if(nr>1) then
               call ZGEMV('C',nr,nb-1,tau_new,ab(nb,ns+1),2*nb-1,hv_new,1,(0.d0,0.d0),h(2),1)
               x = dot_product(hs(1:nr),hv_new(1:nr))*tau_new
               h(2:nb) = h(2:nb) - x*hv(2:nb)
               ! Unfortunately the is no BLAS routine like DGER2 for a nonsymmetric rank 2 update
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

         ! Use new HH vector for the next block
         hv(:) = hv_new(:)
         tau = tau_new

      enddo
      call cpu_time(finish)
      tstep2 = finish-start
      t_2 = t_2 + tstep2
   enddo
   call cpu_time(finish_1)
   tstep1 = finish_1-start_1
   t_1 = t_1 + tstep1

!   print *,' timing  inner and outer loop', t_2, t_1  
   ! Finish the last outstanding requests
   call mpi_wait(ireq_ab,MPI_STATUS_IGNORE,mpierr)
   call mpi_wait(ireq_hv,MPI_STATUS_IGNORE,mpierr)

   call mpi_waitall(nblocks, ireq_hhs, MPI_STATUSES_IGNORE, mpierr)
   call mpi_waitall(num_chunks, ireq_hhr, MPI_STATUSES_IGNORE, mpierr)

   call mpi_barrier(mpi_comm,mpierr)

   deallocate(ab)
   deallocate(ireq_hhr, ireq_hhs)
   deallocate(hh_cnt, hh_dst)
   deallocate(hh_gath, hh_send)
   deallocate(limits)
   deallocate( snd_limits)
   deallocate(block_limits)
   deallocate(global_id)

end subroutine

!---------------------------------------------------------------------------------------------------

! redist_band: redistributes band from 2D block cyclic form to 1D band

subroutine redist_band(l_real, r_a, c_a, lda, na, nblk, nbw, mpi_comm_rows, mpi_comm_cols, mpi_comm, r_ab, c_ab)

   logical, intent(in)     :: l_real
   real*8, intent(in)      :: r_a(lda, *)
   complex*16, intent(in)  :: c_a(lda, *)
   integer, intent(in)     :: lda, na, nblk, nbw, mpi_comm_rows, mpi_comm_cols, mpi_comm
   real*8, intent(out)     :: r_ab(:,:)
   complex*16, intent(out) :: c_ab(:,:)

   integer, allocatable :: ncnt_s(:), nstart_s(:), ncnt_r(:), nstart_r(:), global_id(:,:), block_limits(:)
   real*8, allocatable :: r_sbuf(:,:,:), r_rbuf(:,:,:), r_buf(:,:)
   complex*16, allocatable :: c_sbuf(:,:,:), c_rbuf(:,:,:), c_buf(:,:)

   integer i, j, my_pe, n_pes, my_prow, np_rows, my_pcol, np_cols, nfact, np, npr, npc, mpierr, is, js
   integer nblocks_total, il, jl, l_rows, l_cols, n_off

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
      if(j/nfact==block_limits(np+1)) np = np+1
      if(mod(j,np_rows) == my_prow) then
         do i=0,nfact
            if(mod(i+j,np_cols) == my_pcol) then
               ncnt_s(np) = ncnt_s(np) + 1
            endif
         enddo
      endif
   enddo

   ! Allocate send buffer

   if(l_real) then
      allocate(r_sbuf(nblk,nblk,sum(ncnt_s)))
      r_sbuf(:,:,:) = 0.
   else
      allocate(c_sbuf(nblk,nblk,sum(ncnt_s)))
      c_sbuf(:,:,:) = 0.
   endif

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
      if(j/nfact==block_limits(np+1)) np = np+1
      if(mod(j,np_rows) == my_prow) then
         do i=0,nfact
            if(mod(i+j,np_cols) == my_pcol) then
               nstart_s(np) = nstart_s(np) + 1
               js = (j/np_rows)*nblk
               is = ((i+j)/np_cols)*nblk
               jl = MIN(nblk,l_rows-js)
               il = MIN(nblk,l_cols-is)
               if(l_real) then
                  r_sbuf(1:jl,1:il,nstart_s(np)) = r_a(js+1:js+jl,is+1:is+il)
               else
                  c_sbuf(1:jl,1:il,nstart_s(np)) = c_a(js+1:js+jl,is+1:is+il)
               endif
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

   if(l_real) then
      allocate(r_rbuf(nblk,nblk,sum(ncnt_r)))
   else
      allocate(c_rbuf(nblk,nblk,sum(ncnt_r)))
   endif

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

   if(l_real) then
      call MPI_Alltoallv(r_sbuf,ncnt_s,nstart_s,MPI_REAL8,r_rbuf,ncnt_r,nstart_r,MPI_REAL8,mpi_comm,mpierr)
   else
      call MPI_Alltoallv(c_sbuf,ncnt_s,nstart_s,MPI_COMPLEX16,c_rbuf,ncnt_r,nstart_r,MPI_COMPLEX16,mpi_comm,mpierr)
   endif

   ! set band from receive buffer

   ncnt_r(:) = ncnt_r(:)/(nblk*nblk)

   nstart_r(0) = 0
   do i=1,n_pes-1
      nstart_r(i) = nstart_r(i-1) + ncnt_r(i-1)
   enddo

   if(l_real) then
      allocate(r_buf((nfact+1)*nblk,nblk))
   else
      allocate(c_buf((nfact+1)*nblk,nblk))
   endif

   ! n_off: Offset of ab within band
   n_off = block_limits(my_pe)*nbw

   do j=block_limits(my_pe)*nfact,min(block_limits(my_pe+1)*nfact-1,(na-1)/nblk)
      npr = mod(j,np_rows)
      do i=0,nfact
         npc = mod(i+j,np_cols)
         np = global_id(npr,npc)
         nstart_r(np) = nstart_r(np) + 1
         if(l_real) then
            r_buf(i*nblk+1:i*nblk+nblk,:) = transpose(r_rbuf(:,:,nstart_r(np)))
         else
            c_buf(i*nblk+1:i*nblk+nblk,:) = conjg(transpose(c_rbuf(:,:,nstart_r(np))))
         endif
      enddo
      do i=1,MIN(nblk,na-j*nblk)
         if(l_real) then
            r_ab(1:nbw+1,i+j*nblk-n_off) = r_buf(i:i+nbw,i)
         else
            c_ab(1:nbw+1,i+j*nblk-n_off) = c_buf(i:i+nbw,i)
         endif
      enddo
   enddo

   deallocate(ncnt_s, nstart_s)
   deallocate(ncnt_r, nstart_r)
   deallocate(global_id)
   deallocate(block_limits)
   if(l_real) then
      deallocate(r_sbuf, r_rbuf, r_buf)
   else
      deallocate(c_sbuf, c_rbuf, c_buf)
   endif

end subroutine

!---------------------------------------------------------------------------------------------------
! divide_band: sets the work distribution in band
! Proc n works on blocks block_limits(n)+1 .. block_limits(n+1)

subroutine divide_band(nblocks_total, n_pes, block_limits)

   integer, intent(in) :: nblocks_total ! total number of blocks in band
   integer, intent(in) :: n_pes         ! number of PEs for division
   integer, intent(out) :: block_limits(0:n_pes)

   integer :: n, nblocks, nblocks_left

   block_limits(0) = 0
   if(nblocks_total < n_pes) then
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
         if(n<=n_pes-nblocks_left) then
            block_limits(n) = block_limits(n-1) + nblocks
         else
            block_limits(n) = block_limits(n-1) + nblocks + 1
         endif
      enddo
   endif

end subroutine

