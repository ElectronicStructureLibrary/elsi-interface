subroutine trans_ev_tridi_to_band_real(na, nev, nblk, nbw, q, ldq, mpi_comm_rows, mpi_comm_cols)

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
!
!  mpi_comm_rows
!  mpi_comm_cols
!              MPI-Communicators for rows/columns/both
!
!-------------------------------------------------------------------------------

    implicit none
    integer, intent(in) :: na, nev, nblk, nbw, ldq, mpi_comm_rows, mpi_comm_cols
    real*8 q(ldq,*)

    integer np_rows, my_prow, np_cols, my_pcol
    integer tmp
    integer i, j, ip, sweep, nbuf, l_nev, a_dim2
    integer current_n, current_local_n, current_n_start, current_n_end
    integer next_n, next_local_n, next_n_start, next_n_end
    integer bottom_msg_length, top_msg_length, next_top_msg_length
    integer stripe_width, last_stripe_width, stripe_count
    integer num_result_blocks, num_result_buffers, num_bufs_recvd
    integer a_off, current_tv_off, max_blk_size
    integer mpierr, src, src_offset, dst, offset, nfact, num_blk
    logical flag

    real*8, allocatable :: row(:), row_group(:,:)
    real*8, allocatable :: top_border_send_buffer(:,:,:), top_border_recv_buffer(:,:,:)
    real*8, allocatable :: bottom_border_send_buffer(:,:,:), bottom_border_recv_buffer(:,:,:)
    real*8, allocatable :: result_buffer(:,:,:)
    real*8, allocatable :: bcast_buffer(:,:)

    integer(c_size_t)  :: a_dev
    integer(c_size_t)  :: bcast_buffer_dev
    integer(c_size_t)  :: num
    integer(c_size_t)  :: dev_offset, dev_offset_1


     integer(c_size_t)  :: row_dev
     integer(c_size_t)  :: row_group_dev
    integer(c_size_t)  :: hh_dot_dev
    integer(c_size_t)  :: hh_tau_dev

    integer ierr,istat
    Integer :: top, chunk, this_chunk
    integer row_group_size, unpack_idx

    integer n_off
    integer, allocatable :: result_send_request(:), result_recv_request(:), limits(:)
    integer, allocatable :: top_send_request(:), bottom_send_request(:)
    integer, allocatable :: top_recv_request(:), bottom_recv_request(:)

    ! MPI send/recv tags, arbitrary

    integer, parameter :: bottom_recv_tag = 111
    integer, parameter :: top_recv_tag    = 222
    integer, parameter :: result_recv_tag = 333

    ! Just for measuring the kernel performance
    real*8 kernel_time
    integer*8 kernel_flops      

    unpack_idx = 0
    row_group_size = 0
    kernel_time = 1.d-100
    kernel_flops = 0

    call MPI_Comm_rank(mpi_comm_rows, my_prow, mpierr)
    call MPI_Comm_size(mpi_comm_rows, np_rows, mpierr)
    call MPI_Comm_rank(mpi_comm_cols, my_pcol, mpierr)
    call MPI_Comm_size(mpi_comm_cols, np_cols, mpierr)

    if(mod(nbw,nblk)/=0) then
      if(my_prow==0 .and. my_pcol==0) then
         print *,'ERROR: nbw=',nbw,', nblk=',nblk
         print *,'band backtransform works only for nbw==n*nblk'
         call mpi_abort(mpi_comm_world,0,mpierr)
      endif
    endif

    nfact = nbw / nblk

    ! Local number of eigenvectors
    ! Important: if l_nev is 0, the current process will fail to allocate device memory and terminate
    !            This happens when the dataset is too small for the number of MPI processes
    l_nev = local_index(nev, my_pcol, np_cols, nblk, -1)

    if(l_nev == 0) then
        stripe_width = 0
        stripe_count = 0
        last_stripe_width = 0
    else
        stripe_width = 256 ! Must be a multiple of 4
        stripe_count = (l_nev - 1) / stripe_width + 1
        last_stripe_width = l_nev - (stripe_count - 1) * stripe_width
    endif

    ! Determine the matrix distribution at the beginning
    allocate(limits(0 : np_rows))

    call determine_workload(na, nbw, np_rows, limits)

    max_blk_size = maxval(limits(1 : np_rows) - limits(0 : np_rows - 1))
    a_dim2 = max_blk_size + nbw

    num =  (stripe_width*a_dim2*stripe_count)*8_8
    istat = cuda_malloc(a_dev, stripe_width*a_dim2*stripe_count*8_8)
    istat = cuda_memset(a_dev , 0, num)

    allocate(row(l_nev))
    row(:) = 0

    num =  (l_nev)*8_8
    istat = cuda_malloc( row_dev,l_nev*8_8)
    istat = cuda_memset(row_dev , 0, num)




    ! "row_group" and "row_group_dev" are needed for GPU optimizations
    allocate(row_group(l_nev, nblk))
    row_group(:, :) = 0

    num =  (l_nev*nblk)*8_8
!    call cuda_malloc2d( row_group_dev,l_nev*8_8,nblk*8_8)
    istat = cuda_malloc(row_group_dev, l_nev*nblk*8_8)
    istat = cuda_memset(row_group_dev , 0, num)


    ! Copy q from a block cyclic distribution into a distribution with contiguous rows,
    ! and transpose the matrix using stripes of given stripe_width for cache blocking.

    ! The peculiar way it is done below is due to the fact that the last row should be
    ! ready first since it is the first one to start below

    do ip = np_rows - 1, 0, -1
        if (my_prow == ip) then
            ! Receive my rows which have not yet been received
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            do i=limits(ip)+1,limits(ip+1)
                src = mod((i-1)/nblk, np_rows)
                if(src < my_prow) then
                    ! An unpacking of the current row group may occur before queuing the next row 
                    call unpack_and_prepare_row_group(i - limits(ip), .false.)
                    call MPI_Recv(row_group(:, row_group_size), l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)    
                elseif(src==my_prow) then
                    src_offset = src_offset+1
                    ! An unpacking of the current row group may occur before queuing the next row
                    call unpack_and_prepare_row_group(i - limits(ip), .false.)
                    row_group(:, row_group_size) = q(src_offset, 1:l_nev)
                endif
            enddo
            ! Send all rows which have not yet been send
            src_offset = 0
            do dst = 0, ip-1
              do i=limits(dst)+1,limits(dst+1)
                if(mod((i-1)/nblk, np_rows) == my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    call MPI_Send(row, l_nev, MPI_REAL8, dst, 0, mpi_comm_rows, mpierr)
                endif
              enddo
            enddo
        else if(my_prow < ip) then
            ! Send all rows going to PE ip
            src_offset = local_index(limits(ip), my_prow, np_rows, nblk, -1)
            do i=limits(ip)+1,limits(ip+1)
                src = mod((i-1)/nblk, np_rows)
                if(src == my_prow) then
                    src_offset = src_offset+1
                    row(:) = q(src_offset, 1:l_nev)
                    call MPI_Send(row, l_nev, MPI_REAL8, ip, 0, mpi_comm_rows, mpierr)
                endif
            enddo
            ! Receive all rows from PE ip
            do i=limits(my_prow)+1,limits(my_prow+1)
                src = mod((i-1)/nblk, np_rows)
                if(src == ip) then
                    ! An unpacking of the current row group may occur before queuing the next row
                    call unpack_and_prepare_row_group(i - limits(my_prow), .false.)
                    call MPI_Recv(row_group(:, row_group_size), l_nev, MPI_REAL8, src, 0, mpi_comm_rows, MPI_STATUS_IGNORE, mpierr)
                endif
            enddo
        endif
    enddo

    ! Force an unpacking of all remaining rows that haven't been unpacked yet
    call unpack_and_prepare_row_group(-1, .true.)
    istat = cuda_devicesynchronize()

    ! Set up result buffer queue

    num_result_blocks = ((na-1)/nblk + np_rows - my_prow) / np_rows

    num_result_buffers = 4*nfact
    allocate(result_buffer(l_nev,nblk,num_result_buffers))

    allocate(result_send_request(num_result_buffers))
    allocate(result_recv_request(num_result_buffers))
    result_send_request(:) = MPI_REQUEST_NULL
    result_recv_request(:) = MPI_REQUEST_NULL

    ! Queue up buffers

    if(my_prow > 0 .and. l_nev>0) then ! note: row 0 always sends
        do j = 1, min(num_result_buffers, num_result_blocks)
            call MPI_Irecv(result_buffer(1,1,j), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                           mpi_comm_rows, result_recv_request(j), mpierr)
        enddo
    endif

    num_bufs_recvd = 0 ! No buffers received yet

    ! Initialize top/bottom requests

    allocate(top_send_request(stripe_count))
    allocate(top_recv_request(stripe_count))
    allocate(bottom_send_request(stripe_count))
    allocate(bottom_recv_request(stripe_count))

    top_send_request(:) = MPI_REQUEST_NULL
    top_recv_request(:) = MPI_REQUEST_NULL
    bottom_send_request(:) = MPI_REQUEST_NULL
    bottom_recv_request(:) = MPI_REQUEST_NULL

    allocate(top_border_send_buffer(stripe_width, nbw, stripe_count))
    allocate(top_border_recv_buffer(stripe_width, nbw, stripe_count))
    allocate(bottom_border_send_buffer(stripe_width, nbw, stripe_count))
    allocate(bottom_border_recv_buffer(stripe_width, nbw, stripe_count))

    top_border_send_buffer(:,:,:) = 0
    top_border_recv_buffer(:,:,:) = 0
    bottom_border_send_buffer(:,:,:) = 0
    bottom_border_recv_buffer(:,:,:) = 0

    ! Initialize broadcast buffer

    allocate(bcast_buffer(nbw, max_blk_size))
    bcast_buffer = 0

    num =  ( nbw * max_blk_size) * 8_8
    istat = cuda_malloc(bcast_buffer_dev, nbw * max_blk_size * 8_8)
    istat = cuda_memset( bcast_buffer_dev, 0, num)

    num =  ((max_blk_size-1))*8_8
    istat = cuda_malloc( hh_dot_dev, (max_blk_size -1) * 8_8)
    istat = cuda_memset( hh_dot_dev, 0, num)

    num =  (max_blk_size)*8_8
    istat = cuda_malloc( hh_tau_dev, max_blk_size * 8_8)
    istat = cuda_memset( hh_tau_dev, 0, num)


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

        if(next_n_end < next_n) then
            bottom_msg_length = current_n_end - next_n_end
        else
            bottom_msg_length = 0
        endif

        if(next_local_n > 0) then
            next_top_msg_length = current_n_start - next_n_start
        else
            next_top_msg_length = 0
        endif

        if(sweep==0 .and. current_n_end < current_n .and. l_nev > 0) then
            do i = 1, stripe_count
                call MPI_Irecv(bottom_border_recv_buffer(1,1,i), nbw*stripe_width, MPI_REAL8, my_prow+1, bottom_recv_tag, &
                           mpi_comm_rows, bottom_recv_request(i), mpierr)
            enddo
        endif

        if(current_local_n > 1) then
            if(my_pcol == mod(sweep,np_cols)) then
                bcast_buffer(:,1:current_local_n) = hh_trans_real(:,current_tv_off+1:current_tv_off+current_local_n)
                current_tv_off = current_tv_off + current_local_n
            endif
            call mpi_bcast(bcast_buffer, nbw*current_local_n, MPI_REAL8, mod(sweep,np_cols), mpi_comm_cols, mpierr)
            istat =  cuda_memcpy(bcast_buffer_dev, loc(bcast_buffer(1,1)), nbw * current_local_n * 8_8 , 1)
  

            call extract_hh_tau(nbw, current_local_n, .false.)
            call compute_hh_dot_products(nbw, current_local_n)
        else
            ! for current_local_n == 1 the one and only HH vector is 0 and not stored in hh_trans_real
            bcast_buffer(:, 1) = 0
!           bcast_buffer_dev(:, 1) = 0
            istat = cuda_memset(bcast_buffer_dev, 0, nbw * 8_8)
            call extract_hh_tau(nbw, 1, .true.)
        endif

        if(l_nev == 0) cycle

        if(current_local_n > 0) then

          do i = 1, stripe_count

            !wait_b
            if(current_n_end < current_n) then
                call MPI_Wait(bottom_recv_request(i), MPI_STATUS_IGNORE, mpierr)
                n_off = current_local_n+a_off
                dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width *a_dim2 )) *8        
                istat =  cuda_memcpy( a_dev + dev_offset , loc(bottom_border_recv_buffer(1,1,i)) ,stripe_width*nbw*8_8 ,1)
       if(next_n_end < next_n) then
                    call MPI_Irecv(bottom_border_recv_buffer(1,1,i), nbw*stripe_width, MPI_REAL8, my_prow+1, bottom_recv_tag, &
                                   mpi_comm_rows, bottom_recv_request(i), mpierr)
                endif
            endif

            if(current_local_n <= bottom_msg_length + top_msg_length) then

                !wait_t
                if(top_msg_length>0) then
                    call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)

                   dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *8
                   istat =  cuda_memcpy( a_dev+dev_offset , loc(top_border_recv_buffer(1,1,i)),stripe_width*top_msg_length*8_8 ,1)

                endif

                !compute
                call compute_hh_trafo(0, current_local_n, i)

                !send_b
                call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
                if(bottom_msg_length>0) then
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *8
                    istat =  cuda_memcpy( loc(bottom_border_send_buffer(1,1,i)), a_dev + dev_offset, stripe_width * bottom_msg_length * 8_8 ,2)

                    call MPI_Isend(bottom_border_send_buffer(1,1,i), bottom_msg_length*stripe_width, MPI_REAL8, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                endif

            else

                !compute
                call compute_hh_trafo(current_local_n - bottom_msg_length, bottom_msg_length, i)

                !send_b
                call MPI_Wait(bottom_send_request(i), MPI_STATUS_IGNORE, mpierr)
                if(bottom_msg_length > 0) then
                    n_off = current_local_n+nbw-bottom_msg_length+a_off
                    dev_offset = (0 + (n_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *8
                    istat =  cuda_memcpy( loc(bottom_border_send_buffer(1,1,i)), a_dev + dev_offset,stripe_width*bottom_msg_length*8_8 ,2)


                    call MPI_Isend(bottom_border_send_buffer(1,1,i), bottom_msg_length*stripe_width, MPI_REAL8, my_prow+1, &
                                   top_recv_tag, mpi_comm_rows, bottom_send_request(i), mpierr)
                endif

                !compute
                call compute_hh_trafo(top_msg_length, current_local_n-top_msg_length-bottom_msg_length, i)

                !wait_t
                if(top_msg_length>0) then
                    call MPI_Wait(top_recv_request(i), MPI_STATUS_IGNORE, mpierr)

                   dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *8
                   istat =  cuda_memcpy( a_dev + dev_offset , loc( top_border_recv_buffer(1,1,i)), stripe_width * top_msg_length *8_8 ,1)


     

                endif

                !compute
                call compute_hh_trafo(0, top_msg_length, i)
            endif

            if(next_top_msg_length > 0) then
                !request top_border data
                call MPI_Irecv(top_border_recv_buffer(1,1,i), next_top_msg_length*stripe_width, MPI_REAL8, my_prow-1, &
                               top_recv_tag, mpi_comm_rows, top_recv_request(i), mpierr)
            endif

            !send_t
            if(my_prow > 0) then
                call MPI_Wait(top_send_request(i), MPI_STATUS_IGNORE, mpierr)

                dev_offset = (0 + (a_off * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) * 8
                istat =  cuda_memcpy( loc(top_border_send_buffer(1,1,i)), a_dev + dev_offset, stripe_width*nbw*8_8 ,2)


                call MPI_Isend(top_border_send_buffer(1,1,i), nbw*stripe_width, MPI_REAL8, my_prow-1, bottom_recv_tag, &
                               mpi_comm_rows, top_send_request(i), mpierr)
            endif

            ! Care that there are not too many outstanding top_recv_request's
            if(stripe_count > 1) then
                if(i > 1) then
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

        if(my_prow == 0) then

            ! topmost process sends nbw rows to destination processes

            do j=0,nfact-1

                num_blk = sweep*nfact+j ! global number of destination block, 0 based
                if(num_blk*nblk >= na) exit

                nbuf = mod(num_blk, num_result_buffers) + 1 ! buffer number to get this block

                call MPI_Wait(result_send_request(nbuf), MPI_STATUS_IGNORE, mpierr)

                dst = mod(num_blk, np_rows)

                if(dst == 0) then
                    row_group_size = min(na - num_blk*nblk, nblk)
                    call pack_row_group(row_group(:, :), j * nblk + a_off, row_group_size)

                    do i = 1, row_group_size
                        q((num_blk / np_rows) * nblk + i, 1 : l_nev) = row_group(:, i)
                    enddo
                else
                    call pack_row_group(result_buffer(:, :, nbuf), j * nblk + a_off, nblk)
                    call MPI_Isend(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, dst, &
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

                if(next_local_n > 0) then
                    call MPI_Test(result_recv_request(nbuf), flag, MPI_STATUS_IGNORE, mpierr)
                    if(.not. flag) exit
                else
                    call MPI_Wait(result_recv_request(nbuf), MPI_STATUS_IGNORE, mpierr)
                endif

                ! Fill result buffer into q
                num_blk = j * np_rows + my_prow ! global number of current block, 0 based
                do i = 1, min(na - num_blk*nblk, nblk)
                    q(j * nblk + i, 1 : l_nev) = result_buffer(1 : l_nev, i, nbuf)
                enddo

                ! Queue result buffer again if there are outstanding blocks left
                if(j + num_result_buffers < num_result_blocks) &
                    call MPI_Irecv(result_buffer(1,1,nbuf), l_nev*nblk, MPI_REAL8, 0, result_recv_tag, &
                                   mpi_comm_rows, result_recv_request(nbuf), mpierr)

            enddo
            num_bufs_recvd = j

        endif

        ! Shift the remaining rows to the front of A (if necessary)

        offset = nbw - top_msg_length
        if(offset < 0) then
            call MPI_Abort(MPI_COMM_WORLD, 1, mpierr)
        endif
        a_off = a_off + offset
        if(a_off + next_local_n + nbw > a_dim2) then
            do i = 1, stripe_count

                chunk = min(next_local_n - 1, a_off)
                do j = top_msg_length + 1, top_msg_length + next_local_n, chunk
                   top = min(j + chunk, top_msg_length + next_local_n)
                   this_chunk = top - j + 1
                   dev_offset = (0 + ( (j-1) * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *8
                   dev_offset_1 = (0 + ( (j + a_off-1) * stripe_width) + ( (i-1) * stripe_width * a_dim2 )) *8
                   tmp = cuda_d2d(1)
                   ierr =  cuda_memcpy( a_dev + dev_offset , a_dev +dev_offset_1, stripe_width*this_chunk*8_8, tmp)

                enddo                
                             
            enddo
            a_off = 0
        endif

    enddo

    ! Just for safety:
    if(ANY(top_send_request    /= MPI_REQUEST_NULL)) print *,'*** ERROR top_send_request ***',my_prow,my_pcol
    if(ANY(bottom_send_request /= MPI_REQUEST_NULL)) print *,'*** ERROR bottom_send_request ***',my_prow,my_pcol
    if(ANY(top_recv_request    /= MPI_REQUEST_NULL)) print *,'*** ERROR top_recv_request ***',my_prow,my_pcol
    if(ANY(bottom_recv_request /= MPI_REQUEST_NULL)) print *,'*** ERROR bottom_recv_request ***',my_prow,my_pcol

    if(my_prow == 0) then
        call MPI_Waitall(num_result_buffers, result_send_request, MPI_STATUSES_IGNORE, mpierr)
    endif

    if(ANY(result_send_request /= MPI_REQUEST_NULL)) print *,'*** ERROR result_send_request ***',my_prow,my_pcol
    if(ANY(result_recv_request /= MPI_REQUEST_NULL)) print *,'*** ERROR result_recv_request ***',my_prow,my_pcol

    if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
        print '(" Kernel time:",f10.3," MFlops: ",f10.3)', kernel_time, kernel_flops/kernel_time*1.d-6

    ! deallocate all working space

    deallocate(row)
    deallocate(limits)
    deallocate(result_send_request)
    deallocate(result_recv_request)
    deallocate(top_border_send_buffer)
    deallocate(top_border_recv_buffer)
    deallocate(bottom_border_send_buffer)
    deallocate(bottom_border_recv_buffer)
    deallocate(result_buffer)
    deallocate(bcast_buffer)
    deallocate(top_send_request)
    deallocate(top_recv_request)
    deallocate(bottom_send_request)
    deallocate(bottom_recv_request)
    istat = cuda_free(a_dev) 

    istat = cuda_free(hh_dot_dev)

    istat = cuda_free(hh_tau_dev)

    istat = cuda_free(row_dev)
    
    deallocate(row_group)

    istat= cuda_free(row_group_dev)

    istat =  cuda_free(bcast_buffer_dev)


contains







    ! Pack a filled row group (i.e. an array of consecutive rows)
    subroutine pack_row_group(rows, n_offset, row_count)
        implicit none
        integer, intent(in) :: n_offset, row_count
!        real*8, intent(out) :: rows(:, :)
        real*8 ::rows(:,:)
        integer max_idx

        ! Use many blocks for higher GPU occupancy
        max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

        ! Use one kernel call to pack the entire row group

        
!        call my_pack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, a_dev, row_group_dev)

     call launch_my_pack_c_kernel(row_count, n_offset, max_idx, stripe_width, a_dim2, stripe_count, l_nev, a_dev, row_group_dev)

    
        ! Issue one single transfer call for all rows (device to host)
!        rows(:, 1 : row_count) = row_group_dev(:, 1 : row_count)
        
 
        istat =  cuda_memcpy( loc(rows(1, 1)), row_group_dev , row_count * l_nev * 8_8 ,2)
        !write(*,*) cudaGetErrorString(istat)

    end subroutine


    ! Unpack a filled row group (i.e. an array of consecutive rows)
    subroutine unpack_row_group(rows, n_offset, row_count)        
        implicit none        
        integer, intent(in) :: n_offset, row_count
        real*8, intent(in) :: rows(:, :)
        integer max_idx
        integer i

        ! Use many blocks for higher GPU occupancy
        max_idx = (stripe_count - 1) * stripe_width + last_stripe_width

        ! Issue one single transfer call for all rows (host to device)
!        row_group_dev(:, 1 : row_count) = rows(:, 1 : row_count)

       !istat =  cuda_memcpy( row_group_dev , loc(rows(:, 1: row_count)),row_count * l_nev * 8_8 ,1)

        istat =  cuda_memcpy( row_group_dev , loc(rows(1, 1)),row_count * l_nev * 8_8 ,1)
        !write(*,*) cudaGetErrorString(istat)



        ! Use one kernel call to pack the entire row group
!         call my_unpack_kernel<<<grid_size, stripe_width>>>(n_offset, max_idx, stripe_width, a_dim2, stripe_count, row_group_dev, a_dev)            
 

        call launch_my_unpack_c_kernel( row_count, n_offset, max_idx,stripe_width,a_dim2, stripe_count, l_nev, row_group_dev,a_dev)



    end subroutine


    ! This subroutine must be called before queuing the next row for unpacking; it ensures that an unpacking of the current row group
    ! occurs when the queue is full or when the next row belongs to another group 
    subroutine unpack_and_prepare_row_group(next_unpack_idx, force)
        implicit none        
        integer, intent(in) :: next_unpack_idx
        logical, intent(in) :: force

        if (row_group_size == 0) then
            ! Nothing to flush, just prepare for the upcoming row
            row_group_size = 1
        else
            if (force .or. (row_group_size == nblk) .or. (unpack_idx + 1 /= next_unpack_idx)) then
                ! A flush and a reset must be performed
                call unpack_row_group(row_group(:, :), unpack_idx - row_group_size, row_group_size)
                row_group_size = 1
            else
                ! Just prepare for the upcoming row
                row_group_size = row_group_size + 1
            endif
        endif
        ! Always update the index for the upcoming row
        unpack_idx = next_unpack_idx    
    end subroutine

    ! Host wrapper for the Householder backtransformation step. Several kernels are available. Performance note:
    ! - "compute_hh_trafo_c_kernel" is the C kernel for the backtransformation (this exhibits best performance)
    ! - "compute_hh_trafo_kernel" is the Fortran equivalent of the C kernel
    ! - "compute_hh_trafo_single_kernel" is the reference Fortran kernel
    subroutine compute_hh_trafo(off, ncols, istripe)
        implicit none
        integer, intent(in) :: off, ncols, istripe
        integer nl
        real*8 ttt

        ! ncols - indicates the number of HH reflectors to apply; at least 1 must be available
        if (ncols < 1) return

        ttt = mpi_wtime()
        nl = merge(stripe_width, last_stripe_width, istripe < stripe_count)

        !! Uncomment the kernel you want to use; comment out the other 2

        dev_offset = (0 + (a_off * stripe_width) + ( (istripe - 1) * stripe_width *a_dim2 )) *8
        call launch_compute_hh_trafo_c_kernel(a_dev + dev_offset, bcast_buffer_dev, hh_dot_dev, hh_tau_dev, nl, nbw, stripe_width, off, ncols)
        ! Since we only use the default CUDA stream, no explicit device synchronization is necessary here
        kernel_flops = kernel_flops + 4 * int(nl, 8) * int(ncols, 8) * int(nbw, 8)
        kernel_time = kernel_time + mpi_wtime() - ttt

    end subroutine


    ! The host wrapper for computing the dot products between consecutive HH reflectors (see the kernel below)
    subroutine compute_hh_dot_products(nbw, n)
        implicit none
        integer, value :: nbw, n   

       if (n .le. 1) return
       call launch_compute_hh_dotp_c_kernel( bcast_buffer_dev, hh_dot_dev, nbw, n)


    end subroutine


    ! The host wrapper for extracting "tau" from the HH reflectors (see the kernel below)
    subroutine extract_hh_tau(nbw, n, is_zero)
        implicit none
        integer, value :: nbw, n
        logical, value :: is_zero
        integer val_is_zero
        if(is_zero) then
        val_is_zero = 1
        else
         val_is_zero = 0
        endif

      call launch_extract_hh_tau_c_kernel(bcast_buffer_dev,hh_tau_dev, nbw, n, val_is_zero)
    end subroutine

end subroutine

! -------------------------------------------
! Fortran back-transformation support kernels
! -------------------------------------------

! Reset a reduction block
! Limitation: the thread-block size must be a divider of the reduction block's size
! Reset 2 reduction blocks without an explicit synchronization at the end
! Limitation: : the thread-block size must be a divider of the reduction block's size
! Perform a reduction on an initialized, 128-element shared block
! Compute the dot-product between 2 consecutive HH vectors
! Limitation 1: the size of the thread block must be at most 128 and a power-of-2
! Limitation 2: the size of the warp must be equal to 32

! Extract "tau" from the HH matrix and replace it with 1.0 or 0.0 (depending on case)
! Having "tau" as the first element in a HH reflector reduces space requirements, but causes undesired branching in the kernels

! -------------------------------------------
! Fortran back-transformation support kernels
! -------------------------------------------

! This is the simplest and slowest available backtransformation kernel 

! This is an improved version of the simple backtransformation kernel; here, we halve the number of iterations and apply
! 2 Householder reflectors per iteration


! ---------------------------------
! Row packing and unpacking kernels
! ---------------------------------

! The row group packing kernel

! ----------------------
! Workload configuration
! ----------------------

subroutine determine_workload(na, nb, nprocs, limits)

    integer, intent(in) :: na, nb, nprocs
    integer, intent(out) :: limits(0:nprocs)

    integer i

    if(na <= 0) then
        limits(:) = 0
        return
    endif

    if(nb*nprocs > na) then
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
