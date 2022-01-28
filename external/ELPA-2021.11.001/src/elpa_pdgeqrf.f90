










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


module elpa_pdgeqrf

  use elpa_utilities
  use elpa1_compute
  use elpa_pdlarfb
  use qr_utils_mod
  use elpa_qrkernels
  use elpa_mpi
  implicit none

  PRIVATE

  public :: qr_pdgeqrf_2dcomm_double
  public :: qr_pqrparam_init
  public :: qr_pdlarfg2_1dcomm_check_double

  public :: qr_pdgeqrf_2dcomm_single
  public :: qr_pdlarfg2_1dcomm_check_single


  contains
  ! real double precision
















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
     subroutine qr_pdgeqrf_2dcomm_&
           &double &
           (obj, a, lda, matrixCols, v, ldv, vmrCols, tau, lengthTau, t, ldt, colsT, &
                                  work, workLength, lwork, m, n, mb, nb, rowidx, colidx, &
                                  rev, trans, PQRPARAM, mpicomm_rows, mpicomm_cols, blockheuristic)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter   :: gmode_ = 1, rank_ = 2, eps_ = 3

      ! input variables (local)
      integer(kind=ik), intent(in)  :: lda, lwork, ldv, ldt, matrixCols, m, vmrCols, lengthTau, &
                                       colsT, workLength

      ! input variables (global)
      integer(kind=ik)              :: n, mb, nb, rowidx, colidx, rev, trans, mpicomm_cols, mpicomm_rows
      integer(kind=ik)              :: PQRPARAM(1:11)
      real(kind=c_double)                 :: a(1:lda,1:matrixCols), v(1:ldv,1:vmrCols), tau(1:lengthTau), &
                                       t(1:ldt,1:colsT), work(1:workLength)
      ! output variables (global)
      real(kind=c_double)                 :: blockheuristic(*)

      ! input variables derived from PQRPARAM
      integer(kind=ik)              :: updatemode,tmerge,size2d

      ! local scalars
      integer(kind=ik)              :: mpirank_cols,broadcast_size,mpirank_rows
      integer(kind=MPI_KIND)        :: mpirank_colsMPI, mpirank_rowsMPI
      integer(kind=MPI_KIND)        :: mpierr
      integer(kind=ik)              :: mpirank_cols_qr,mpiprocs_cols
      integer(kind=MPI_KIND)        :: mpiprocs_colsMPI
      integer(kind=ik)              :: lcols_temp,lcols,icol,lastcol
      integer(kind=ik)              :: baseoffset,offset,idx,voffset
      integer(kind=ik)              :: update_voffset,update_tauoffset
      integer(kind=ik)              :: update_lcols
      integer(kind=ik)              :: work_offset

      real(kind=c_double)                 :: dbroadcast_size(1),dtmat_bcast_size(1)
      real(kind=c_double)                 :: pdgeqrf_size(1),pdlarft_size(1),pdlarfb_size(1),tmerge_pdlarfb_size(1)
      integer(kind=ik)              :: temptau_offset,temptau_size,broadcast_offset,tmat_bcast_size
      integer(kind=ik)              :: remaining_cols
      integer(kind=ik)              :: total_cols
      integer(kind=ik)              :: incremental_update_size ! needed for incremental update mode

      call obj%timer%start("qr_pdgeqrf_2dcomm_&
          &double&
          &")
      size2d     = PQRPARAM(1)
      updatemode = PQRPARAM(2)
      tmerge     = PQRPARAM(3)

      ! copy value before we are going to filter it
      total_cols = n
      call mpi_comm_rank(int(mpicomm_cols,kind=MPI_KIND) ,mpirank_colsMPI, mpierr)
      call mpi_comm_rank(int(mpicomm_rows,kind=MPI_KIND) ,mpirank_rowsMPI, mpierr)
      call mpi_comm_size(int(mpicomm_cols,kind=MPI_KIND) ,mpiprocs_colsMPI, mpierr)

      mpirank_cols = int(mpirank_colsMPI,kind=c_int)
      mpirank_rows = int(mpirank_rowsMPI,kind=c_int)
      mpiprocs_cols = int(mpiprocs_colsMPI,kind=c_int)

      call qr_pdgeqrf_1dcomm_&
          &double &
          (obj,a,lda,v,ldv,tau,t,ldt,pdgeqrf_size(1),-1,m,total_cols,mb,rowidx,rowidx,rev,trans, &
                             PQRPARAM(4:11),mpicomm_rows,blockheuristic)
      call qr_pdgeqrf_pack_unpack_&
          &double &
          (obj,v,ldv,dbroadcast_size(1),-1,m,total_cols,mb,rowidx,rowidx,rev,0,mpicomm_rows)
      call qr_pdgeqrf_pack_unpack_tmatrix_&
          &double &
          (obj,tau,t,ldt,dtmat_bcast_size(1),-1,total_cols,0)

      pdlarft_size(1) = 0.0_rk8

      call qr_pdlarfb_1dcomm_&
          &double &
          (m,mb,total_cols,total_cols,a,lda,v,ldv,tau,t,ldt,rowidx,rowidx,rev,mpicomm_rows, &
                             pdlarfb_size(1),-1)
      call qr_tmerge_pdlarfb_1dcomm_&
          &double &
          (m,mb,total_cols,total_cols,total_cols,v,ldv,t,ldt,a,lda,rowidx,rev,updatemode, &
                                    mpicomm_rows,tmerge_pdlarfb_size(1),-1)


      temptau_offset = 1
      temptau_size = total_cols
      broadcast_offset = temptau_offset + temptau_size
      broadcast_size = int(dbroadcast_size(1) + dtmat_bcast_size(1))
      work_offset = broadcast_offset + broadcast_size

      if (lwork .eq. -1) then
        work(1) = (real(temptau_size,kind=c_double) + real(broadcast_size,kind=c_double) + max(pdgeqrf_size(1), &
            pdlarft_size(1),pdlarfb_size(1), &
                   tmerge_pdlarfb_size(1)))
        call obj%timer%stop("qr_pdgeqrf_2dcomm_&
          &double&
          &")
        return
      end if

      lastcol = colidx-total_cols+1
      voffset = total_cols

      incremental_update_size = 0

      ! clear v buffer: just ensure that there is no junk in the upper triangle
      ! part, otherwise pdlarfb gets some problems
      ! pdlarfl(2) do not have these problems as they are working more on a Vector
      ! basis
      v(1:ldv,1:total_cols) = 0.0_rk8
      icol = colidx

      remaining_cols = total_cols

      !print *,'start decomposition',m,rowidx,colidx

      do while (remaining_cols .gt. 0)

        ! determine rank of process column with next qr block
        mpirank_cols_qr = MOD((icol-1)/nb,mpiprocs_cols)

        ! lcols can't be larger than than nb
        ! exception: there is only one process column

        ! however, we might not start at the first local column.
        ! therefore assume a matrix of size (1xlcols) starting at (1,icol)
        ! determine the real amount of local columns
        lcols_temp = min(nb,(icol-lastcol+1))

        ! blocking parameter
        lcols_temp = max(min(lcols_temp,size2d),1)

        ! determine size from last decomposition column
        !  to first decomposition column
        call local_size_offset_1d(icol,nb,icol-lcols_temp+1,icol-lcols_temp+1,0, &
                                      mpirank_cols_qr,mpiprocs_cols, &
                                      lcols,baseoffset,offset)

        voffset = remaining_cols - lcols + 1

        idx = rowidx - colidx + icol

        if (mpirank_cols .eq. mpirank_cols_qr) then
          ! qr decomposition part
          tau(offset:offset+lcols-1) = 0.0_rk8

          call qr_pdgeqrf_1dcomm_&
              &double &
              (obj,a(1,offset),lda,v(1,voffset),ldv,tau(offset),t(voffset,voffset),ldt, &
                                 work(work_offset),lwork,m,lcols,mb,rowidx,idx,rev,trans,PQRPARAM(4:11), &
                                 mpicomm_rows,blockheuristic)

          ! pack broadcast buffer (v + tau)
          call qr_pdgeqrf_pack_unpack_&
              &double &
              (obj,v(1,voffset),ldv,work(broadcast_offset),lwork,m,lcols,mb,rowidx,&
                                      idx,rev,0,mpicomm_rows)

          ! determine broadcast size
          call qr_pdgeqrf_pack_unpack_&
              &double &
              (obj,v(1,voffset),ldv,dbroadcast_size(1),-1,m,lcols,mb,rowidx,idx,rev,&
                                      0,mpicomm_rows)
          broadcast_size = int(dbroadcast_size(1))

          !if (mpirank_rows .eq. 0) then
          ! pack tmatrix into broadcast buffer and calculate new size
          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &double &
              (obj,tau(offset),t(voffset,voffset),ldt, &
                                              work(broadcast_offset+broadcast_size),lwork,lcols,0)
          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &double &
              (obj,tau(offset),t(voffset,voffset),ldt,dtmat_bcast_size(1),-1,lcols,0)
          broadcast_size = broadcast_size + int(dtmat_bcast_size(1))
          !end if

          ! initiate broadcast (send part)

          call MPI_Bcast(work(broadcast_offset), int(broadcast_size,kind=MPI_KIND), mpi_real8, &
                         int(mpirank_cols_qr,kind=MPI_KIND), int(mpicomm_cols,kind=MPI_KIND), mpierr)

          ! copy tau parts into temporary tau buffer
          work(temptau_offset+voffset-1:temptau_offset+(voffset-1)+lcols-1) = tau(offset:offset+lcols-1)

          !print *,'generated tau:', tau(offset)
        else
          ! Vector exchange part

          ! determine broadcast size
          call qr_pdgeqrf_pack_unpack_&
              &double &
              (obj,v(1,voffset),ldv,dbroadcast_size(1),-1,m,lcols,mb,rowidx,idx,rev,1,mpicomm_rows)
          broadcast_size = int(dbroadcast_size(1))

          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &double &
              (obj,work(temptau_offset+voffset-1),t(voffset,voffset),ldt, &
                                              dtmat_bcast_size(1),-1,lcols,0)
          tmat_bcast_size = dtmat_bcast_size(1)

          !print *,'broadcast_size (nonqr)',broadcast_size
          broadcast_size = dbroadcast_size(1) + dtmat_bcast_size(1)

          ! initiate broadcast (recv part)

          call MPI_Bcast(work(broadcast_offset), int(broadcast_size,kind=MPI_KIND), mpi_real8, &
                         int(mpirank_cols_qr,kind=MPI_KIND), int(mpicomm_cols,kind=MPI_KIND), mpierr)

          ! last n*n elements in buffer are (still empty) T matrix elements
          ! fetch from first process in each column

          ! unpack broadcast buffer (v + tau)
          call qr_pdgeqrf_pack_unpack_&
              &double &
              (obj,v(1,voffset),ldv,work(broadcast_offset),lwork,m,lcols, &
              mb,rowidx,idx,rev,1,mpicomm_rows)

          ! now send t matrix to other processes in our process column
          broadcast_size = int(dbroadcast_size(1))
          tmat_bcast_size = int(dtmat_bcast_size(1))

          ! t matrix should now be available on all processes => unpack
          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &double &
              (obj,work(temptau_offset+voffset-1),t(voffset,voffset),ldt, &
                                              work(broadcast_offset+broadcast_size),lwork,lcols,1)
        end if

        remaining_cols = remaining_cols - lcols

        ! apply householder vectors to whole trailing matrix parts (if any)

        update_voffset = voffset
        update_tauoffset = icol
        update_lcols = lcols
        incremental_update_size = incremental_update_size + lcols

        icol = icol - lcols
        ! count colums from first column of global block to current index
        call local_size_offset_1d(icol,nb,colidx-n+1,colidx-n+1,0, &
                                      mpirank_cols,mpiprocs_cols, &
                                      lcols,baseoffset,offset)

        if (lcols .gt. 0) then

          !print *,'updating trailing matrix'

           if (updatemode .eq. ichar('I')) then
             print *,'pdgeqrf_2dcomm: incremental update not yet implemented! rev=1'
           else if (updatemode .eq. ichar('F')) then
             ! full update no merging
             call qr_pdlarfb_1dcomm_&
                &double &
                (m,mb,lcols,update_lcols,a(1,offset),lda,v(1,update_voffset),ldv, &
                     work(temptau_offset+update_voffset-1),                          &
                                                        t(update_voffset,update_voffset),ldt, &
                    rowidx,idx,1,mpicomm_rows,work(work_offset),lwork)
           else
            ! full update + merging default
             call qr_tmerge_pdlarfb_1dcomm_&
                &double &
                (m,mb,lcols,n-(update_voffset+update_lcols-1),update_lcols, &
                                                              v(1,update_voffset),ldv, &
                          t(update_voffset,update_voffset),ldt, &
                          a(1,offset),lda,rowidx,1,updatemode,mpicomm_rows, &
                                                              work(work_offset),lwork)
           end if
        else
           if (updatemode .eq. ichar('I')) then
             !print *,'sole merging of (incremental) T matrix', mpirank_cols,  &
            !                            n-(update_voffset+incremental_update_size-1)
             call qr_tmerge_pdlarfb_1dcomm_&
                &double &
                (m,mb,0,n-(update_voffset+incremental_update_size-1),   &
                                                              incremental_update_size,v(1,update_voffset),ldv, &
                          t(update_voffset,update_voffset),ldt, &
                          a,lda,rowidx,1,updatemode,mpicomm_rows,work(work_offset),lwork)

             ! reset for upcoming incremental updates
             incremental_update_size = 0
          else if (updatemode .eq. ichar('M')) then
             ! final merge
            call qr_tmerge_pdlarfb_1dcomm_&
                &double &
                (m,mb,0,n-(update_voffset+update_lcols-1),update_lcols, &
                                                              v(1,update_voffset),ldv, &
                          t(update_voffset,update_voffset),ldt, &
                          a,lda,rowidx,1,updatemode,mpicomm_rows,work(work_offset),lwork)
          else
            ! full updatemode - nothing to update
          end if

          ! reset for upcoming incremental updates
          incremental_update_size = 0
        end if
      end do

      if ((tmerge .gt. 0) .and. (updatemode .eq. ichar('F'))) then
        ! finally merge all small T parts
        call qr_pdlarft_tree_merge_1dcomm_&
&double &
(m,mb,n,size2d,tmerge,v,ldv,t,ldt,rowidx,rev,mpicomm_rows,work,lwork)
      end if

      !print *,'stop decomposition',rowidx,colidx
      call obj%timer%stop("qr_pdgeqrf_2dcomm_&
          &double&
          &")
    end subroutine

    subroutine qr_pdgeqrf_1dcomm_&
&double &
(obj,a,lda,v,ldv,tau,t,ldt,work,lwork,m,n,mb,baseidx,rowidx,rev,trans, &
          PQRPARAM,mpicomm,blockheuristic)
      use precision
      use elpa1_impl
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter  :: gmode_ = 1,rank_ = 2,eps_ = 3

      ! input variables (local)
      integer(kind=ik)             :: lda,lwork,ldv,ldt
      real(kind=c_double)                :: a(lda,*),v(ldv,*),tau(*),t(ldt,*),work(*)

      ! input variables (global)
      integer(kind=ik)             :: m,n,mb,baseidx,rowidx,rev,trans,mpicomm
      integer(kind=ik)             :: PQRPARAM(:)
      ! derived input variables

      ! derived further input variables from QR_PQRPARAM
      integer(kind=ik)             :: size1d,updatemode,tmerge

      ! output variables (global)
      real(kind=c_double)                :: blockheuristic(*)

      ! local scalars
      integer(kind=ik)             :: nr_blocks,remainder,current_block,aoffset,idx,updatesize
      real(kind=c_double)                :: pdgeqr2_size(1),pdlarfb_size(1),tmerge_tree_size(1)
      call obj%timer%start("qr_pdgeqrf_1dcomm_&
          &double&
          &")
      size1d     = max(min(PQRPARAM(1),n),1)
      updatemode = PQRPARAM(2)
      tmerge     = PQRPARAM(3)

      if (lwork .eq. -1) then
        call qr_pdgeqr2_1dcomm_&
&double &
(obj,a,lda,v,ldv,tau,t,ldt,pdgeqr2_size,-1, &
                                  m,size1d,mb,baseidx,baseidx,rev,trans,PQRPARAM(4:),mpicomm,blockheuristic)
        ! reserve more space for incremental mode
        call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,n,n,n,v,ldv,t,ldt, &
                                         a,lda,baseidx,rev,updatemode,mpicomm,pdlarfb_size,-1)

        call qr_pdlarft_tree_merge_1dcomm_&
&double &
(m,mb,n,size1d,tmerge,v,ldv,t,ldt,baseidx,rev,mpicomm,tmerge_tree_size,-1)

        work(1) = max(pdlarfb_size(1),pdgeqr2_size(1),tmerge_tree_size(1))
      call obj%timer%stop("qr_pdgeqrf_1dcomm_&
          &double&
          &")
        return
      end if

      nr_blocks = n / size1d
      remainder = n - nr_blocks*size1d

      current_block = 0
      do while (current_block .lt. nr_blocks)
        idx = rowidx-current_block*size1d
        updatesize = n-(current_block+1)*size1d
        aoffset = 1+updatesize
        call qr_pdgeqr2_1dcomm_&
&double &
(obj,a(1,aoffset),lda,v(1,aoffset),ldv,tau(aoffset),t(aoffset,aoffset),ldt,work,lwork, &
                                m,size1d,mb,baseidx,idx,1,trans,PQRPARAM(4:),mpicomm,blockheuristic)
        if (updatemode .eq. ichar('M')) then
          ! full update + merging
          call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,updatesize,current_block*size1d,size1d, &
                                           v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                           a,lda,baseidx,1,ichar('F'),mpicomm,work,lwork)
        else if (updatemode .eq. ichar('I')) then
          if (updatesize .ge. size1d) then
            ! incremental update + merging
            call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,size1d,current_block*size1d,size1d, &
                                               v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                               a(1,aoffset-size1d),lda,baseidx,1,updatemode,mpicomm,work,lwork)

          else ! only remainder left
             ! incremental update + merging
             call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,remainder,current_block*size1d,size1d, &
                                               v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                               a(1,1),lda,baseidx,1,updatemode,mpicomm,work,lwork)
          end if
        else ! full update no merging is default
          ! full update no merging
          call qr_pdlarfb_1dcomm_&
&double &
(m,mb,updatesize,size1d,a,lda,v(1,aoffset),ldv, &
                                    tau(aoffset),t(aoffset,aoffset),ldt,baseidx,idx,1,mpicomm,work,lwork)
        end if

        ! move on to next block
        current_block = current_block+1
      end do

      if (remainder .gt. 0) then
        aoffset = 1
        idx = rowidx-size1d*nr_blocks
        call qr_pdgeqr2_1dcomm_&
&double &
(obj,a(1,aoffset),lda,v,ldv,tau,t,ldt,work,lwork, &
                                  m,remainder,mb,baseidx,idx,1,trans,PQRPARAM(4:),mpicomm,blockheuristic)
        if ((updatemode .eq. ichar('I')) .or. (updatemode .eq. ichar('M'))) then
          ! final merging
          call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,0,size1d*nr_blocks,remainder, &
                                             v,ldv,t,ldt, &
                                             a,lda,baseidx,1,updatemode,mpicomm,work,lwork) ! updatemode argument does not matter
        end if
      end if

      if ((tmerge .gt. 0) .and. (updatemode .eq. ichar('F'))) then
        ! finally merge all small T parts
        call qr_pdlarft_tree_merge_1dcomm_&
&double &
(m,mb,n,size1d,tmerge,v,ldv,t,ldt,baseidx,rev,mpicomm,work,lwork)
      end if
      call obj%timer%stop("qr_pdgeqrf_1dcomm_&
          &double&
          &")

    end subroutine

    ! local a and tau are assumed to be positioned at the right column from a local
    ! perspective
    ! TODO: if local amount of data turns to zero the algorithm might produce wrong
    ! results (probably due to old buffer contents)
    subroutine qr_pdgeqr2_1dcomm_&
&double &
(obj,a,lda,v,ldv,tau,t,ldt,work,lwork,m,n,mb,baseidx,rowidx,rev, &
          trans,PQRPARAM,mpicomm,blockheuristic)
      use precision
      !use elpa1_impl ! check this
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter   :: gmode_ = 1,rank_ = 2 ,eps_ = 3, upmode1_ = 4

      ! input variables (local)
      integer(kind=ik)              :: lda,lwork,ldv,ldt
      real(kind=c_double)                 :: a(lda,*),v(ldv,*),tau(*),t(ldt,*),work(*)

      ! input variables (global)
      integer(kind=ik)              :: m,n,mb,baseidx,rowidx,rev,trans,mpicomm
      integer(kind=ik)              :: PQRPARAM(:)
      ! output variables (global)
      real(kind=c_double)                 :: blockheuristic(*)

      ! derived further input variables from QR_PQRPARAM
      integer(kind=ik)              ::  maxrank,hgmode,updatemode

      ! local scalars
      integer(kind=ik)              :: icol,incx,idx
      real(kind=c_double)                 :: pdlarfg_size(1),pdlarf_size(1),total_size
      real(kind=c_double)                 :: pdlarfg2_size(1),pdlarfgk_size(1),pdlarfl2_size(1)
      real(kind=c_double)                 :: pdlarft_size(1),pdlarfb_size(1),pdlarft_pdlarfb_size(1),tmerge_pdlarfb_size(1)
      integer(kind=ik)              :: mpirank,mpiprocs
      integer(kind=MPI_KIND)        :: mpirankMPI, mpiprocsMPI
      integer(kind=MPI_KIND)        :: mpierr
      integer(kind=ik)              :: rank,lastcol,actualrank,nextrank
      integer(kind=ik)              :: update_cols,decomposition_cols
      integer(kind=ik)              :: current_column
      call obj%timer%start("qr_pdgeqr2_1dcomm_&
          &double&
          &")

      maxrank    = min(PQRPARAM(1),n)
      updatemode = PQRPARAM(2)
      hgmode     = PQRPARAM(4)
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      if (trans .eq. 1) then
        incx = lda
      else
        incx = 1
      end if

      if (lwork .eq. -1) then

        call qr_pdlarfg_1dcomm_&
&double &
(obj,a,incx,tau(1),pdlarfg_size(1),-1,n,rowidx,mb,hgmode,rev,mpicomm)
        call qr_pdlarfl_1dcomm_&
&double &
(v,1,baseidx,a,lda,tau(1),pdlarf_size(1),-1,m,n,rowidx,mb,rev,mpicomm)
        call qr_pdlarfg2_1dcomm_ref_&
&double &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,pdlarfg2_size(1),-1,m,rowidx,mb,PQRPARAM(:), &
                                    rev,mpicomm,actualrank)

        call qr_pdlarfgk_1dcomm_&
&double &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,pdlarfgk_size(1),-1,m,n, &
            rowidx,mb,PQRPARAM(:),rev,mpicomm,actualrank)
        call qr_pdlarfl2_tmatrix_1dcomm_&
&double &
(v,ldv,baseidx,a,lda,t,ldt,pdlarfl2_size(1),-1,m,n,rowidx,mb,rev,mpicomm)
        pdlarft_size(1) = 0.0_rk8
        call qr_pdlarfb_1dcomm_&
&double &
(m,mb,n,n,a,lda,v,ldv,tau,t,ldt,baseidx,rowidx,1,mpicomm,pdlarfb_size(1),-1)
        pdlarft_pdlarfb_size(1) = 0.0_rk8
        call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,n,n,n,v,ldv,t,ldt,a,lda,rowidx,rev,&
            updatemode,mpicomm,tmerge_pdlarfb_size(1),-1)

        total_size = max(pdlarfg_size(1),pdlarf_size(1),pdlarfg2_size(1),pdlarfgk_size(1),pdlarfl2_size(1),pdlarft_size(1), &
                         pdlarfb_size(1),pdlarft_pdlarfb_size(1),tmerge_pdlarfb_size(1))

        work(1) = total_size
      call obj%timer%stop("qr_pdgeqr2_1dcomm_&
          &double&
          &")
        return
      end if

      icol = 1
      lastcol = min(rowidx,n)
      decomposition_cols = lastcol
      update_cols = n
      do while (decomposition_cols .gt. 0) ! local qr block
        icol = lastcol-decomposition_cols+1
        idx = rowidx-icol+1

        ! get possible rank size
        ! limited by number of columns and remaining rows
        rank = min(n-icol+1,maxrank,idx)

        current_column = n-icol+1-rank+1

        if (rank .eq. 1) then

          call qr_pdlarfg_1dcomm_&
&double &
(obj,a(1,current_column),incx, &
                                  tau(current_column),work,lwork, &
                                  m,idx,mb,hgmode,1,mpicomm)
          v(1:ldv,current_column) = 0.0_rk8
          call qr_pdlarfg_copy_1dcomm_&
&double &
(obj,a(1,current_column),incx, &
                                       v(1,current_column),1, &
                                       m,baseidx,idx,mb,1,mpicomm)

          ! initialize t matrix part
          t(current_column,current_column) = tau(current_column)

          actualrank = 1

        else if (rank .eq. 2) then
          call qr_pdlarfg2_1dcomm_ref_&
&double &
(obj,a(1,current_column),lda,tau(current_column), &
                                         t(current_column,current_column),ldt,v(1,current_column),ldv, &
                                        baseidx,work,lwork,m,idx,mb,PQRPARAM(:),1,mpicomm,actualrank)
        else
          call qr_pdlarfgk_1dcomm_&
&double &
(obj,a(1,current_column),lda,tau(current_column), &
                                     t(current_column,current_column),ldt,v(1,current_column),ldv, &
                                     baseidx,work,lwork,m,rank,idx,mb,PQRPARAM(:),1,mpicomm,actualrank)
        end if

        blockheuristic(actualrank) = blockheuristic(actualrank) + 1

        ! the blocked decomposition versions already updated their non
        ! decomposed parts using their information after communication
        update_cols = decomposition_cols - rank
        decomposition_cols = decomposition_cols - actualrank

        ! needed for incremental update
        nextrank = min(n-(lastcol-decomposition_cols+1)+1,maxrank,rowidx-(lastcol-decomposition_cols+1)+1)

        if (current_column .gt. 1) then
          idx = rowidx-icol+1

          if (updatemode .eq. ichar('I')) then
            ! incremental update + merging
            call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,nextrank-(rank-actualrank),n-(current_column+rank-1),actualrank, &
                                          v(1,current_column+(rank-actualrank)),ldv, &
                                          t(current_column+(rank-actualrank),current_column+(rank-actualrank)),ldt, &
                                          a(1,current_column-nextrank+(rank-actualrank)),lda,baseidx,rev,updatemode,&
                                          mpicomm,work,lwork)
          else
            ! full update + merging
            call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,update_cols,n-(current_column+rank-1),actualrank, &
                                          v(1,current_column+(rank-actualrank)),ldv, &
                                          t(current_column+(rank-actualrank),current_column+(rank-actualrank)),ldt, &
                                          a(1,1),lda,baseidx,rev,updatemode,mpicomm,work,lwork)
          end if
        else
          call qr_tmerge_pdlarfb_1dcomm_&
&double &
(m,mb,0,n-(current_column+rank-1),actualrank, &
              v(1,current_column+(rank-actualrank)), &
                                          ldv, &
                                          t(current_column+(rank-actualrank),current_column+(rank-actualrank)),ldt, &
                                           a,lda,baseidx,rev,updatemode,mpicomm,work,lwork)
        end if

      end do
      call obj%timer%stop("qr_pdgeqr2_1dcomm_&
          &double&
          &")
    end subroutine

    ! incx == 1: column major
    ! incx != 1: row major
    subroutine qr_pdlarfg_1dcomm_&
&double &
(obj,x,incx,tau,work,lwork,n,idx,nb,hgmode,rev,communicator)

      use precision
      !use elpa1_impl !check this
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter    :: gmode_ = 1,rank_ = 2, eps_ = 3

      ! input variables (local)
      integer(kind=ik)               :: incx,lwork,hgmode
      real(kind=c_double)                  :: x(*),work(*)

      ! input variables (global)
      integer(kind=ik)               :: communicator,nb,idx,n,rev

      ! output variables (global)
      real(kind=c_double)                  :: tau

      ! local scalars
      integer(kind=ik)               :: mpirank,mpiprocs,mpirank_top
      integer(kind=MPI_KIND)         :: mpirankMPI, mpiprocsMPI
      integer(kind=MPI_KIND)         :: mpierr
      integer(kind=ik)               :: sendsize,recvsize
      integer(kind=ik)               :: local_size,local_offset,baseoffset
      integer(kind=ik)               :: topidx,top,iproc
      real(kind=c_double)                  :: alpha,xnorm,dot,xf

      ! external functions
      real(kind=c_double), external        :: ddot,dlapy2,dnrm2
      external                       :: dscal

      ! intrinsic
!      intrinsic sign
      call obj%timer%start("qr_pdlarfg_1dcomm_&
          &double&
          &")
      if (idx .le. 1) then
        tau = 0.0d0
      call obj%timer%stop("qr_pdlarfg_1dcomm_&
          &double&
          &")
        return
       end if
      call MPI_Comm_rank(int(communicator,kind=MPI_KIND) , mpirankMPI, mpierr)
      call MPI_Comm_size(int(communicator,kind=MPI_KIND) , mpiprocsMPI, mpierr)
 
      mpirank = int(mpirankMPI, kind=c_int)
      mpiprocs = int(mpiprocsMPI, kind=c_int)

      ! calculate expected work size and store in work(1)
      if (hgmode .eq. ichar('s')) then
        ! allreduce (MPI_SUM)
        sendsize = 2
        recvsize = sendsize
      else if (hgmode .eq. ichar('x')) then
        ! alltoall
        sendsize = mpiprocs*2
        recvsize = sendsize
      else if (hgmode .eq. ichar('g')) then
        ! allgather
        sendsize = 2
        recvsize = mpiprocs*sendsize
      else
        ! no exchange at all (benchmarking)
        sendsize = 2
        recvsize = sendsize
      end if

      if (lwork .eq. -1) then
        work(1) = real(sendsize + recvsize,kind=c_double)

      call obj%timer%stop("qr_pdlarfg_1dcomm_&
          &double&
          &")
        return
      end if

      ! Processor id for global index of top element
      mpirank_top = MOD((idx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        topidx = local_index(idx,mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
      end if

      call local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                        local_size,baseoffset,local_offset)

      local_offset = local_offset * incx

      ! calculate and exchange information
      if (hgmode .eq. ichar('s')) then
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk8
        end if
        dot = ddot(local_size, &
                     x(local_offset), incx, &
                     x(local_offset), incx)
        work(1) = alpha
        work(2) = dot

        call mpi_allreduce(work(1),work(sendsize+1), &
                             int(sendsize,kind=MPI_KIND), mpi_real8, mpi_sum, &
                             int(communicator,kind=MPI_KIND), mpierr)

        alpha = work(sendsize+1)
        xnorm = sqrt(work(sendsize+2))
      else if (hgmode .eq. ichar('x')) then
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk8
        end if
        xnorm = dnrm2(local_size, x(local_offset), incx)
        do iproc=0,mpiprocs-1
          work(2*iproc+1) = alpha
          work(2*iproc+2) = xnorm
        end do

        call mpi_alltoall(work(1), 2_MPI_KIND, mpi_real8, &
                            work(sendsize+1), 2_MPI_KIND, mpi_real8, &
                            int(communicator,kind=MPI_KIND), mpierr)

        ! extract alpha value
        alpha = work(sendsize+1+mpirank_top*2)

        ! copy norm parts of buffer to beginning
        do iproc=0,mpiprocs-1
          work(iproc+1) = work(sendsize+1+2*iproc+1)
        end do

        xnorm = dnrm2(mpiprocs, work(1), 1)
      else if (hgmode .eq. ichar('g')) then
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk8
        end if
        xnorm = dnrm2(local_size, x(local_offset), incx)
        work(1) = alpha
        work(2) = xnorm

        ! allgather

        call mpi_allgather(work(1), int(sendsize,kind=MPI_KIND), mpi_real8, &
                            work(sendsize+1), int(sendsize,kind=MPI_KIND), mpi_real8, &
                            int(communicator,kind=MPI_KIND), mpierr)

        ! extract alpha value
        alpha = work(sendsize+1+mpirank_top*2)

        ! copy norm parts of buffer to beginning
        do iproc=0,mpiprocs-1
          work(iproc+1) = work(sendsize+1+2*iproc+1)
        end do
        xnorm = dnrm2(mpiprocs, work(1), 1)
      else
        ! dnrm2
        xnorm = dnrm2(local_size, x(local_offset), incx)
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk8
        end if

        ! no exchange at all (benchmarking)
        xnorm = 0.0_rk8
      end if

      !print *,'ref hg:', idx,xnorm,alpha
      !print *,x(1:n)

      ! calculate householder information
      if (xnorm .eq. 0.0_rk8) then
        ! H = I

        tau = 0.0_rk8
      else
        ! General case
        call hh_transform_real_&
&double &
(obj,alpha,xnorm**2,xf,tau, .false.)
        if (mpirank .eq. mpirank_top) then
          x(top) = alpha
        end if
        call dscal(local_size, xf, &
                     x(local_offset), incx)

        ! TODO: reimplement norm rescale method of
        ! original PDLARFG using mpi?

      end if

      ! useful for debugging
      !print *,'hg:mpirank,idx,beta,alpha:',mpirank,idx,beta,alpha,1.0d0/(beta+alpha),tau
      !print *,x(1:n)
      call obj%timer%stop("qr_pdlarfg_1dcomm_&
          &double&
          &")
    end subroutine

    subroutine qr_pdlarfg2_1dcomm_ref_&
&double &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,work,lwork,m,idx,mb,PQRPARAM,rev,mpicomm,actualk)
      use precision
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter    :: gmode_ = 1,rank_ = 2,eps_ = 3, upmode1_ = 4
      ! input variables (local)
      integer(kind=ik)               :: lda,lwork,ldv,ldt
      real(kind=c_double)                  :: a(lda,*),v(ldv,*),tau(*),work(*),t(ldt,*)

      ! input variables (global)
      integer(kind=ik)               :: m,idx,baseidx,mb,rev,mpicomm
      integer(kind=ik)               :: PQRPARAM(:)
      ! output variables (global)
      integer(kind=ik)               :: actualk

      ! derived input variables from QR_PQRPARAM
      integer(kind=ik)               :: eps

      ! local scalars
      real(kind=c_double)                  :: dseedwork_size(1)
      integer(kind=ik)               :: seedwork_size,seed_size
      integer(kind=ik)               :: seedwork_offset,seed_offset
      logical                        :: accurate
      call obj%timer%start("qr_pdlarfg2_1dcomm_&
          &double&
          &")

      call qr_pdlarfg2_1dcomm_seed_&
&double &
(obj,a,lda,dseedwork_size(1),-1,work,m,mb,idx,rev,mpicomm)
      seedwork_size = int(dseedwork_size(1))
      seed_size = seedwork_size

      if (lwork .eq. -1) then
        work(1) = seedwork_size + seed_size
      call obj%timer%stop("qr_pdlarfg2_1dcomm_&
          &double&
          &")

        return
      end if

      seedwork_offset = 1
      seed_offset = seedwork_offset + seedwork_size

      eps = PQRPARAM(3)

      ! check for border cases (only a 2x2 matrix left)
      if (idx .le. 1) then
        tau(1:2) = 0.0_rk8
         t(1:2,1:2) = 0.0_rk8

      call obj%timer%stop("qr_pdlarfg2_1dcomm_&
          &double&
          &")

         return
      end if

      call qr_pdlarfg2_1dcomm_seed_&
&double &
(obj,a,lda,work(seedwork_offset),lwork,work(seed_offset),m,mb,idx,rev,mpicomm)

      if (eps .gt. 0) then
        accurate = qr_pdlarfg2_1dcomm_check_&
&double &
(obj,work(seed_offset),eps)
      else
        accurate = .true.
      end if

      call qr_pdlarfg2_1dcomm_vector_&
&double &
(obj,a(1,2),1,tau(2),work(seed_offset), &
                                          m,mb,idx,0,1,mpicomm)

      call qr_pdlarfg_copy_1dcomm_&
&double &
(obj,a(1,2),1, &
                                       v(1,2),1, &
                                       m,baseidx,idx,mb,1,mpicomm)

      call qr_pdlarfg2_1dcomm_update_&
&double &
(obj,v(1,2),1,baseidx,a(1,1),lda,work(seed_offset),m,idx,mb,rev,mpicomm)

      ! check for 2x2 matrix case => only one householder Vector will be
      ! generated
      if (idx .gt. 2) then
        if (accurate .eqv. .true.) then
          call qr_pdlarfg2_1dcomm_vector_&
&double &
(obj,a(1,1),1,tau(1),work(seed_offset), &
                                                  m,mb,idx-1,1,1,mpicomm)

          call qr_pdlarfg_copy_1dcomm_&
&double &
(obj,a(1,1),1, &
                                               v(1,1),1, &
                                               m,baseidx,idx-1,mb,1,mpicomm)

          ! generate fuse element
          call qr_pdlarfg2_1dcomm_finalize_tmatrix_&
&double &
(obj,work(seed_offset),tau,t,ldt)

          actualk = 2
        else
          t(1,1) = 0.0_rk8
          t(1,2) = 0.0_rk8
          t(2,2) = tau(2)

          actualk = 1
        end if
      else
        t(1,1) = 0.0_rk8
        t(1,2) = 0.0_rk8
        t(2,2) = tau(2)

        ! no more vectors to create
        tau(1) = 0.0_rk8

        actualk = 2

        !print *,'rank2: no more data'
      end if
      call obj%timer%stop("qr_pdlarfg2_1dcomm_&
          &double&
          &")

    end subroutine

    subroutine qr_pdlarfg2_1dcomm_seed_&
&double &
(obj,a,lda,work,lwork,seed,n,nb,idx,rev,mpicomm)
      use precision
      !use elpa1_impl ! check this
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)        :: lda,lwork
      real(kind=c_double)           :: a(lda,*),work(*),seed(*)

      ! input variables (global)
      integer(kind=ik)        :: n,nb,idx,rev,mpicomm

      ! output variables (global)

      ! external functions
      real(kind=c_double), external :: ddot
      ! local scalars
      real(kind=c_double)           :: top11,top21,top12,top22
      real(kind=c_double)           :: dot11,dot12,dot22
      integer(kind=ik)        :: mpirank, mpiprocs
      integer(kind=MPI_KIND)  :: mpirankMPI, mpiprocsMPI, mpierr
      integer(kind=ik)        :: mpirank_top11,mpirank_top21
      integer(kind=ik)        :: top11_offset,top21_offset
      integer(kind=ik)        :: baseoffset
      integer(kind=ik)        :: local_offset1,local_size1
      integer(kind=ik)        :: local_offset2,local_size2

      call obj%timer%start("qr_pdlarfg2_1dcomm_seed_&
          &double&
          &")

      if (lwork .eq. -1) then
        work(1) = 8.0_rk8

      call obj%timer%stop("qr_pdlarfg2_1dcomm_seed_&
          &double&
          &")
        return
      end if
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI, kind=c_int)
      mpiprocs = int(mpiprocsMPI, kind=c_int)

      call local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                                local_size1,baseoffset,local_offset1)

      call local_size_offset_1d(n,nb,idx,idx-2,rev,mpirank,mpiprocs, &
                                local_size2,baseoffset,local_offset2)

      mpirank_top11 = MOD((idx-1)/nb,mpiprocs)
      mpirank_top21 = MOD((idx-2)/nb,mpiprocs)

      top11_offset = local_index(idx,mpirank_top11,mpiprocs,nb,0)
      top21_offset = local_index(idx-1,mpirank_top21,mpiprocs,nb,0)

      if (mpirank_top11 .eq. mpirank) then
        top11 = a(top11_offset,2)
        top12 = a(top11_offset,1)
      else
        top11 = 0.0_rk8
        top12 = 0.0_rk8
      end if

      if (mpirank_top21 .eq. mpirank) then
        top21 = a(top21_offset,2)
        top22 = a(top21_offset,1)
      else
        top21 = 0.0_rk8
        top22 = 0.0_rk8
      end if

      ! calculate 3 dot products
      dot11 = ddot(local_size1,a(local_offset1,2),1,a(local_offset1,2),1)
      dot12 = ddot(local_size1,a(local_offset1,2),1,a(local_offset1,1),1)
      dot22 = ddot(local_size2,a(local_offset2,1),1,a(local_offset2,1),1)
      ! store results in work buffer
      work(1) = top11
      work(2) = dot11
      work(3) = top12
      work(4) = dot12
      work(5) = top21
      work(6) = top22
      work(7) = dot22
      work(8) = 0.0_rk8! fill up buffer
      ! exchange partial results

      call mpi_allreduce(work, seed, 8_MPI_KIND, mpi_real8, mpi_sum, &
                         int(mpicomm,kind=MPI_KIND), mpierr)


      call obj%timer%stop("qr_pdlarfg2_1dcomm_seed_&
          &double&
          &")
    end subroutine

    logical function qr_pdlarfg2_1dcomm_check_&
&double &
(obj,seed,eps)
      use precision
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj

      ! input variables
      real(kind=c_double)    ::  seed(*)
      integer(kind=ik) :: eps

      ! local scalars
      real(kind=c_double)    :: epsd,first,second,first_second,estimate
      logical          :: accurate
      real(kind=c_double)    :: dot11,dot12,dot22
      real(kind=c_double)    :: top11,top12,top21,top22
      call obj%timer%start("qr_pdlarfg2_1dcomm_check_&
          &double&
          &")

      EPSD = EPS

      top11 = seed(1)
      dot11 = seed(2)
      top12 = seed(3)
      dot12 = seed(4)

      top21 = seed(5)
      top22 = seed(6)
      dot22 = seed(7)

      ! reconstruct the whole inner products
      ! (including squares of the top elements)
      first = dot11 + top11*top11
      second = dot22 + top22*top22 + top12*top12
      first_second = dot12 + top11*top12

      ! zero Householder Vector (zero norm) case
      if (first*second .eq. 0.0_rk8) then
        qr_pdlarfg2_1dcomm_check_&
&double &
 = .false.
      call obj%timer%stop("qr_pdlarfg2_1dcomm_check_&
          &double&
          &")

        return
      end if

      estimate = abs((first_second*first_second)/(first*second))

      !print *,'estimate:',estimate

      ! if accurate the following check holds
      accurate = (estimate .LE. (epsd/(1.0_rk8+epsd)))
      qr_pdlarfg2_1dcomm_check_&
&double &
 = accurate
      call obj%timer%stop("qr_pdlarfg2_1dcomm_check_&
          &double&
          &")

    end function

    ! id=0: first Vector
    ! id=1: second Vector
    subroutine qr_pdlarfg2_1dcomm_vector_&
&double &
(obj,x,incx,tau,seed,n,nb,idx,id,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)        :: incx
      real(kind=c_double)           :: x(*),seed(*),tau

      ! input variables (global)
      integer(kind=ik)        :: n,nb,idx,id,rev,mpicomm

      ! output variables (global)

      ! external functions
      real(kind=c_double), external :: dlapy2
      external                :: dscal
      ! local scalars
      integer(kind=ik)        :: mpirank,mpirank_top,mpiprocs
      integer(kind=MPI_KIND)  :: mpierr, mpirankMPI, mpiprocsMPI
      real(kind=c_double)           :: alpha,dot,beta,xnorm
      integer(kind=ik)        :: local_size,baseoffset,local_offset,top,topidx
      call obj%timer%start("qr_pdlarfg2_1dcomm_vector_&
          &double&
          &")

      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      call local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                                    local_size,baseoffset,local_offset)

      local_offset = local_offset * incx

      ! Processor id for global index of top element
      mpirank_top = MOD((idx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        topidx = local_index(idx,mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
      else
        top = -99
        stop
      end if

      alpha = seed(id*5+1)
      dot = seed(id*5+2)

      xnorm = sqrt(dot)
      if (xnorm .eq. 0.0_rk8) then
        ! H = I

        tau = 0.0_rk8
      else
        ! General case
        beta = sign(dlapy2(alpha, xnorm), alpha)
        tau = (beta+alpha) / beta

        !print *,'hg2',tau,xnorm,alpha
        call dscal(local_size, 1.0_rk8/(beta+alpha), &
                   x(local_offset), incx)

        ! TODO: reimplement norm rescale method of
        ! original PDLARFG using mpi?

        if (mpirank .eq. mpirank_top) then
          x(top) = -beta
        end if

        seed(8) = beta
      end if
      call obj%timer%stop("qr_pdlarfg2_1dcomm_vector_&
          &double&
          &")

    end subroutine

    subroutine qr_pdlarfg2_1dcomm_update_&
&double &
(obj,v,incv,baseidx,a,lda,seed,n,idx,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)   :: incv,lda
      real(kind=c_double)      :: v(*),a(lda,*),seed(*)

      ! input variables (global)
      integer(kind=ik)   :: n,baseidx,idx,nb,rev,mpicomm

      ! output variables (global)

      ! external functions
      external daxpy

      ! local scalars
      integer(kind=ik)   :: mpirank,mpiprocs
      integer(kind=MPI_KIND)   :: mpirankMPI, mpiprocsMPI, mpierr
      integer(kind=ik)   :: local_size,local_offset,baseoffset
      real(kind=c_double)      :: z,coeff,beta
      real(kind=c_double)      :: dot11,dot12,dot22
      real(kind=c_double)      :: top11,top12,top21,top22
      call obj%timer%start("qr_pdlarfg2_1dcomm_update_&
          &double&
          &")

      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      ! seed should be updated by previous householder generation
      ! Update inner product of this column and next column Vector
      top11 = seed(1)
      dot11 = seed(2)
      top12 = seed(3)
      dot12 = seed(4)

      top21 = seed(5)
      top22 = seed(6)
      dot22 = seed(7)
      beta = seed(8)

      call local_size_offset_1d(n,nb,baseidx,idx,rev,mpirank,mpiprocs, &
                                local_size,baseoffset,local_offset)
      baseoffset = baseoffset * incv

      ! zero Householder Vector (zero norm) case
      if (beta .eq. 0.0_rk8) then

      call obj%timer%stop("qr_pdlarfg2_1dcomm_update_&
          &double&
          &")
        return
      end if
      z = (dot12 + top11 * top12) / beta + top12

      !print *,'hg2 update:',baseidx,idx,mpirank,local_size
      call daxpy(local_size, -z, v(baseoffset),1, a(local_offset,1),1)
      ! prepare a full dot22 for update
      dot22 = dot22 + top22*top22

      ! calculate coefficient
      COEFF = z / (top11 + beta)

      ! update inner product of next Vector
      dot22 = dot22 - coeff * (2*dot12 - coeff*dot11)

      ! update dot12 value to represent update with first Vector
      ! (needed for T matrix)
      dot12 = dot12 - COEFF * dot11

      ! update top element of next Vector
      top22 = top22 - coeff * top21
      seed(6) = top22

      ! restore separated dot22 for Vector generation
      seed(7) = dot22  - top22*top22

      !------------------------------------------------------
      ! prepare elements for T matrix
      seed(4) = dot12

      ! prepare dot matrix for fuse element of T matrix
      ! replace top11 value with -beta1
      seed(1) = beta
      call obj%timer%stop("qr_pdlarfg2_1dcomm_update_&
          &double&
          &")

    end subroutine

    ! run this function after second Vector
    subroutine qr_pdlarfg2_1dcomm_finalize_tmatrix_&
&double &
(obj,seed,tau,t,ldt)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik)  :: ldt
      real(kind=c_double)     :: seed(*),t(ldt,*),tau(*)
      real(kind=c_double)     :: dot12,beta1,top21,beta2
      call obj%timer%start("qr_pdlarfg2_1dcomm_finalize_tmatrix_&
          &double&
          &")

      beta1 = seed(1)
      dot12 = seed(4)
      top21 = seed(5)
      beta2 = seed(8)

      !print *,'beta1 beta2',beta1,beta2

      dot12 = dot12 / beta2 + top21
      dot12 = -(dot12 / beta1)

      t(1,1) = tau(1)
      t(1,2) = dot12
      t(2,2) = tau(2)
      call obj%timer%stop("qr_pdlarfg2_1dcomm_finalize_tmatrix_&
          &double&
          &")

    end subroutine

    subroutine qr_pdlarfgk_1dcomm_&
&double &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,work,lwork,m,k,idx,mb,PQRPARAM,rev,mpicomm,actualk)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup

      ! input variables (local)
      integer(kind=ik)    :: lda,lwork,ldv,ldt
      real(kind=c_double)       :: a(lda,*),v(ldv,*),tau(*),work(*),t(ldt,*)

      ! input variables (global)
      integer(kind=ik)    :: m,k,idx,baseidx,mb,rev,mpicomm
      integer(kind=ik)    :: PQRPARAM(:)
      ! output variables (global)
      integer(kind=ik)    :: actualk

      ! local scalars
      integer(kind=ik)    :: ivector
      real(kind=c_double)       :: pdlarfg_size(1),pdlarf_size(1)
      real(kind=c_double)       :: pdlarfgk_1dcomm_seed_size(1),pdlarfgk_1dcomm_check_size(1)
      real(kind=c_double)       :: pdlarfgk_1dcomm_update_size(1)
      integer(kind=ik)    :: seedC_size,seedC_offset
      integer(kind=ik)    :: seedD_size,seedD_offset
      integer(kind=ik)    :: work_offset
      call obj%timer%start("qr_pdlarfgk_1dcomm_&
          &double&
          &")

      seedC_size = k*k
      seedC_offset = 1
      seedD_size = k*k
      seedD_offset = seedC_offset + seedC_size
      work_offset = seedD_offset + seedD_size

      if (lwork .eq. -1) then
        call qr_pdlarfg_1dcomm_&
&double &
(obj, a,1,tau(1),pdlarfg_size(1),-1,m,baseidx,mb,PQRPARAM(4),rev,mpicomm)

        call qr_pdlarfl_1dcomm_&
&double &
(v,1,baseidx,a,lda,tau(1),pdlarf_size(1),-1,m,k,baseidx,mb,rev,mpicomm)
        call qr_pdlarfgk_1dcomm_seed_&
&double &
(obj,a,lda,baseidx,pdlarfgk_1dcomm_seed_size(1),-1,work,work,m,k,mb,mpicomm)
        !call qr_pdlarfgk_1dcomm_check_&
!&double &
!(work,work,k,PQRPARAM(:),pdlarfgk_1dcomm_check_size(1),-1,actualk)
        call qr_pdlarfgk_1dcomm_check_improved_&
&double &
(obj,work,work,k,PQRPARAM(:),pdlarfgk_1dcomm_check_size(1),-1,actualk)
        call qr_pdlarfgk_1dcomm_update_&
&double &
(obj,a,lda,baseidx,pdlarfgk_1dcomm_update_size(1), &
                                              -1,work,work,k,k,1,work,m,mb,rev,mpicomm)
        work(1) = max(pdlarfg_size(1),pdlarf_size(1),pdlarfgk_1dcomm_seed_size(1),pdlarfgk_1dcomm_check_size(1), &
                        pdlarfgk_1dcomm_update_size(1)) + real(seedC_size + seedD_size, kind=c_double)

      call obj%timer%stop("qr_pdlarfgk_1dcomm_&
          &double&
          &")

        return
      end if

      call qr_pdlarfgk_1dcomm_seed_&
&double &
(obj,a(1,1),lda,idx,work(work_offset),lwork,work(seedC_offset), &
          work(seedD_offset),m,k,mb,mpicomm)
      !call qr_pdlarfgk_1dcomm_check_&
!&double &
!(work(seedC_offset),work(seedD_offset),k,PQRPARAM(:),work(work_offset),lwork,actualk)
      call qr_pdlarfgk_1dcomm_check_improved_&
&double &
(obj,work(seedC_offset),work(seedD_offset), &
          k,PQRPARAM(:),work(work_offset),lwork,actualk)
      !print *,'possible rank:', actualk

      ! override useful for debugging
      !actualk = 1
      !actualk = k
      !actualk= min(actualk,2)
      do ivector=1,actualk
        call qr_pdlarfgk_1dcomm_vector_&
&double &
(obj,a(1,k-ivector+1),1,idx,tau(k-ivector+1), &
                                          work(seedC_offset),work(seedD_offset),k, &
                                          ivector,m,mb,rev,mpicomm)

        call qr_pdlarfgk_1dcomm_update_&
&double &
(obj,a(1,1),lda,idx,work(work_offset),lwork,work(seedC_offset), &
                                          work(seedD_offset),k,actualk,ivector,tau, &
                                          m,mb,rev,mpicomm)

        call qr_pdlarfg_copy_1dcomm_&
&double &
(obj,a(1,k-ivector+1),1, &
                                       v(1,k-ivector+1),1, &
                                       m,baseidx,idx-ivector+1,mb,1,mpicomm)
      end do

      ! generate final T matrix and convert preliminary tau values into real ones
      call qr_pdlarfgk_1dcomm_generateT_&
&double &
(obj,work(seedC_offset),work(seedD_offset),k,actualk,tau,t,ldt)

      call obj%timer%stop("qr_pdlarfgk_1dcomm_&
          &double&
          &")
    end subroutine

    subroutine qr_pdlarfgk_1dcomm_seed_&
&double &
(obj,a,lda,baseidx,work,lwork,seedC,seedD,m,k,mb,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup

      ! input variables (local)
      integer(kind=ik)   :: lda,lwork
      real(kind=c_double)      :: a(lda,*), work(*)

      ! input variables (global)
      integer(kind=ik)   :: m,k,baseidx,mb,mpicomm
      real(kind=c_double)      :: seedC(k,*),seedD(k,*)

      ! output variables (global)

      ! derived input variables from QR_PQRPARAM

      ! local scalars
      integer(kind=ik)   :: mpirank,mpiprocs,mpirank_top
      integer(kind=MPI_KIND)   :: mpierr,mpirankMPI,mpiprocsMPI
      integer(kind=ik)   :: icol,irow,lidx,remsize
      integer(kind=ik)   :: remaining_rank

      integer(kind=ik)   :: C_size,D_size,sendoffset,recvoffset,sendrecv_size
      integer(kind=ik)   :: localoffset,localsize,baseoffset
      call obj%timer%start("qr_pdlarfgk_1dcomm_seed_&
          &double&
          &")

      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      C_size = k*k
      D_size = k*k
      sendoffset = 1
      sendrecv_size = C_size+D_size
      recvoffset = sendoffset + sendrecv_size

      if (lwork .eq. -1) then
        work(1) = real(2*sendrecv_size,kind=c_double)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_seed_&
          &double&
          &")

        return
      end if

      ! clear buffer
      work(sendoffset:sendoffset+sendrecv_size-1)=0.0_rk8
      ! collect C part
      do icol=1,k

        remaining_rank = k
        do while (remaining_rank .gt. 0)
          irow = k - remaining_rank + 1
          lidx = baseidx - remaining_rank + 1

          ! determine chunk where the current top element is located
          mpirank_top = MOD((lidx-1)/mb,mpiprocs)

          ! limit max number of remaining elements of this chunk to the block
          ! distribution parameter
          remsize = min(remaining_rank,mb)

          ! determine the number of needed elements in this chunk
          call local_size_offset_1d(lidx+remsize-1,mb, &
                                    lidx,lidx,0, &
                                    mpirank_top,mpiprocs, &
                                    localsize,baseoffset,localoffset)

          !print *,'local rank',localsize,localoffset

          if (mpirank .eq. mpirank_top) then
            ! copy elements to buffer
            work(sendoffset+(icol-1)*k+irow-1:sendoffset+(icol-1)*k+irow-1+localsize-1) &
                          = a(localoffset:localoffset+remsize-1,icol)
          end if

          ! jump to next chunk
          remaining_rank = remaining_rank - localsize
        end do
      end do

      ! collect D part
      call local_size_offset_1d(m,mb,baseidx-k,baseidx-k,1, &
                        mpirank,mpiprocs, &
                        localsize,baseoffset,localoffset)

      !print *,'localsize',localsize,localoffset
      if (localsize > 0) then
        call dsyrk("Upper", "Trans", k, localsize, &
                     1.0_rk8, a(localoffset,1), lda, &
                     0.0_rk8, work(sendoffset+C_size), k)
      else
        work(sendoffset+C_size:sendoffset+C_size+k*k-1) = 0.0_rk8
      end if

      ! TODO: store symmetric part more efficiently

      ! allreduce operation on results

      call mpi_allreduce(work(sendoffset),work(recvoffset),int(sendrecv_size,kind=MPI_KIND), &
                         mpi_real8, mpi_sum, int(mpicomm,kind=MPI_KIND), mpierr)

      ! unpack result from buffer into seedC and seedD
      seedC(1:k,1:k) = 0.0_rk8
      do icol=1,k
        seedC(1:k,icol) = work(recvoffset+(icol-1)*k:recvoffset+icol*k-1)
      end do
      seedD(1:k,1:k) = 0.0_rk8
      do icol=1,k
        seedD(1:k,icol) = work(recvoffset+C_size+(icol-1)*k:recvoffset+C_size+icol*k-1)
      end do

      call obj%timer%stop("qr_pdlarfgk_1dcomm_seed_&
          &double&
          &")

    end subroutine

    ! k is assumed to be larger than two
    subroutine qr_pdlarfgk_1dcomm_check_improved_&
&double &
(obj,seedC,seedD,k,PQRPARAM,work,lwork,possiblerank)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (global)
      integer(kind=ik)   :: k,lwork
      integer(kind=ik)   :: PQRPARAM(:)
      real(kind=c_double)      :: seedC(k,*),seedD(k,*),work(k,*)

      ! output variables (global)
      integer(kind=ik)   :: possiblerank

      ! derived input variables from QR_PQRPARAM
      integer(kind=ik)   :: eps

      ! local variables
      integer(kind=ik)   :: i,j,l
      real(kind=c_double)      :: sum_squares,diagonal_square,epsd,diagonal_root
      real(kind=c_double)      :: dreverse_matrix_work(1)

      ! external functions
      real(kind=c_double), external :: ddot,dlapy2,dnrm2
      external                :: dscal

      call obj%timer%start("qr_pdlarfgk_1dcomm_check_improved_&
          &double&
          &")

      if (lwork .eq. -1) then
        call reverse_matrix_local_&
&double &
(1,k,k,work,k,dreverse_matrix_work,-1)
        work(1,1) = real(k*k,kind=c_double) + dreverse_matrix_work(1)

        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &double&
            &")
        return
      end if

      eps = PQRPARAM(3)

      if (eps .eq. 0) then
        possiblerank = k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &double&
            &")
        return
      end if
      epsd = real(eps,kind=c_double)

      ! build complete inner product from seedC and seedD
      ! copy seedD to work
      work(:,1:k) = seedD(:,1:k)

      ! add inner products of seedC to work
      call dsyrk("Upper", "Trans", k, k, &
                 1.0_rk8, seedC(1,1), k, &
                 1.0_rk8, work, k)

      ! TODO: optimize this part!
      call reverse_matrix_local_&
&double &
(0,k,k,work(1,1),k,work(1,k+1),lwork-2*k)
      call reverse_matrix_local_&
&double &
(1,k,k,work(1,1),k,work(1,k+1),lwork-2*k)

      ! transpose matrix
      do i=1,k
        do j=i+1,k
          work(i,j) = work(j,i)
        end do
      end do


      ! do cholesky decomposition
      i = 0
      do while ((i .lt. k))
        i = i + 1

        diagonal_square = abs(work(i,i))
        diagonal_root  = sqrt(diagonal_square)

        ! zero Householder Vector (zero norm) case
        if ((abs(diagonal_square) .eq. 0.0_rk8) .or. (abs(diagonal_root) .eq. 0.0_rk8)) then
          possiblerank = max(i-1,1)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &double&
            &")
          return
        end if

        ! check if relative error is bounded for each Householder Vector
        ! Householder i is stable iff Househoulder i-1 is "stable" and the accuracy criterion
        ! holds.
        ! first Householder Vector is considered as "stable".

        do j=i+1,k
          work(i,j) = work(i,j) / diagonal_root
          do l=i+1,j
            work(l,j) = work(l,j) - work(i,j) * work(i,l)
          end do
        end do
        !print *,'cholesky step done'

        ! build sum of squares
        if (i .eq. 1) then
          sum_squares = 0.0_rk8
        else
          sum_squares = ddot(i-1,work(1,i),1,work(1,i),1)
        end if
        if (sum_squares .ge. (epsd * diagonal_square)) then
          possiblerank = max(i-1,1)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &double&
            &")
          return
        end if
      end do

      possiblerank = i
      !print *,'possible rank', possiblerank
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &double&
            &")

    end subroutine

    ! TODO: zero Householder Vector (zero norm) case
    ! - check alpha values as well (from seedC)
    subroutine qr_pdlarfgk_1dcomm_check_&
&double &
(obj,seedC,seedD,k,PQRPARAM,work,lwork,possiblerank)
      use precision
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup

      ! input variables (local)

      ! input variables (global)
      integer(kind=ik)   :: k,lwork
      integer(kind=ik)   :: PQRPARAM(:)
      real(kind=c_double)      :: seedC(k,*),seedD(k,*),work(k,*)

      ! output variables (global)
      integer(kind=ik)   :: possiblerank

      ! derived input variables from QR_PQRPARAM
      integer(kind=ik)   :: eps

      ! local scalars
      integer(kind=ik)   :: icol,isqr,iprod
      real(kind=c_double)      :: epsd,sum_sqr,sum_products,diff,temp,ortho,ortho_sum
      real(kind=c_double)      :: dreverse_matrix_work(1)
        call obj%timer%start("qr_pdlarfgk_1dcomm_check_&
            &double&
            &")
      if (lwork .eq. -1) then
        call reverse_matrix_local_&
&double &
(1,k,k,work,k,dreverse_matrix_work,-1)
        work(1,1) = real(k*k,kind=c_double) + dreverse_matrix_work(1)

        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &double&
            &")

        return
      end if

      eps = PQRPARAM(3)

      if (eps .eq. 0) then
        possiblerank = k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &double&
            &")
        return
      end if
      epsd = real(eps,kind=c_double)

      ! copy seedD to work
      work(:,1:k) = seedD(:,1:k)

      ! add inner products of seedC to work
      call dsyrk("Upper", "Trans", k, k, &
                 1.0_rk8, seedC(1,1), k, &
                 1.0_rk8, work, k)

      ! TODO: optimize this part!
      call reverse_matrix_local_&
&double &
(0,k,k,work(1,1),k,work(1,k+1),lwork-2*k)
      call reverse_matrix_local_&
&double &
(1,k,k,work(1,1),k,work(1,k+1),lwork-2*k)

      ! transpose matrix
      do icol=1,k
        do isqr=icol+1,k
          work(icol,isqr) = work(isqr,icol)
        end do
      end do

      ! work contains now the full inner product of the global (sub-)matrix
      do icol=1,k
        ! zero Householder Vector (zero norm) case
        if (abs(work(icol,icol)) .eq. 0.0_rk8) then
          !print *,'too small ', icol, work(icol,icol)
          possiblerank = max(icol,1)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &double&
            &")
          return
        end if

        sum_sqr = 0.0_rk8
        do isqr=1,icol-1
          sum_products = 0.0_rk8
          do iprod=1,isqr-1
            sum_products = sum_products + work(iprod,isqr)*work(iprod,icol)
          end do

          !print *,'divisor',icol,isqr,work(isqr,isqr)
          temp = (work(isqr,icol) - sum_products)/work(isqr,isqr)
          work(isqr,icol) = temp
          sum_sqr = sum_sqr + temp*temp
        end do

        ! calculate diagonal value
        diff = work(icol,icol) - sum_sqr
        if (diff .lt. 0.0_rk8) then
          ! we definitely have a problem now
          possiblerank = icol-1 ! only decompose to previous column (including)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &double&
            &")
          return
        end if
        work(icol,icol) = sqrt(diff)
        ! calculate orthogonality
        ortho = 0.0_rk8
        do isqr=1,icol-1
          ortho_sum = 0.0_rk8
          do iprod=isqr,icol-1
            temp = work(isqr,iprod)*work(isqr,iprod)
            !print *,'ortho ', work(iprod,iprod)
            temp = temp / (work(iprod,iprod)*work(iprod,iprod))
            ortho_sum = ortho_sum + temp
          end do
          ortho = ortho + ortho_sum * (work(isqr,icol)*work(isqr,icol))
        end do

        ! ---------------- with division by zero ----------------------- !

        !ortho = ortho / diff;

        ! if current estimate is not accurate enough, the following check holds
        !if (ortho .gt. epsd) then
        !    possiblerank = icol-1 ! only decompose to previous column (including)
        !    return
        !end if

        ! ---------------- without division by zero ----------------------- !

        ! if current estimate is not accurate enough, the following check holds
        if (ortho .gt. epsd * diff) then
          possiblerank = icol-1 ! only decompose to previous column (including)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &double&
            &")
          return
        end if
      end do

      ! if we get to this point, the accuracy condition holds for the whole block
      possiblerank = k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &double&
            &")
      end subroutine

    !sidx: seed idx
    !k: max rank used during seed phase
    !rank: actual rank (k >= rank)
    subroutine qr_pdlarfgk_1dcomm_vector_&
&double &
(obj,x,incx,baseidx,tau,seedC,seedD,k,sidx,n,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)  :: incx
      real(kind=c_double)     :: x(*),tau

      ! input variables (global)
      integer(kind=ik)  :: n,nb,baseidx,rev,mpicomm,k,sidx
      real(kind=c_double)     :: seedC(k,*),seedD(k,*)

      ! output variables (global)

      ! external functions
      real(kind=c_double), external :: dlapy2,dnrm2
      external                :: dscal

      ! local scalars
      integer(kind=ik)   :: mpirank,mpirank_top,mpiprocs
      integer(kind=MPI_KIND) :: mpirankMPI, mpiprocsMPI, mpierr
      real(kind=c_double)      :: alpha,dot,beta,xnorm
      integer(kind=ik)   :: local_size,baseoffset,local_offset,top,topidx
      integer(kind=ik)   :: lidx
        call obj%timer%start("qr_pdlarfgk_1dcomm_vector_&
            &double&
            &")
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)

      lidx = baseidx-sidx+1
      call local_size_offset_1d(n,nb,baseidx,lidx-1,rev,mpirank,mpiprocs, &
                        local_size,baseoffset,local_offset)

      local_offset = local_offset * incx

      ! Processor id for global index of top element
      mpirank_top = MOD((lidx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        topidx = local_index((lidx),mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
      end if

      alpha = seedC(k-sidx+1,k-sidx+1)
      dot = seedD(k-sidx+1,k-sidx+1)
      ! assemble actual norm from both seed parts
      xnorm = dlapy2(sqrt(dot), dnrm2(k-sidx,seedC(1,k-sidx+1),1))

      if (xnorm .eq. 0.0_rk8) then
        tau = 0.0_rk8
      else

        ! General case

        beta = sign(dlapy2(alpha, xnorm), alpha)
        ! store a preliminary version of beta in tau
        tau = beta

        ! update global part
        call dscal(local_size, 1.0_rk8/(beta+alpha), &
                     x(local_offset), incx)
        ! do not update local part here due to
        ! dependency of c Vector during update process

        ! TODO: reimplement norm rescale method of
        ! original PDLARFG using mpi?

        if (mpirank .eq. mpirank_top) then
          x(top) = -beta
        end if
      end if
        call obj%timer%stop("qr_pdlarfgk_1dcomm_vector_&
            &double&
            &")

    end subroutine

    !k: original max rank used during seed function
    !rank: possible rank as from check function
    ! TODO: if rank is less than k, reduce buffersize in such a way
    ! that only the required entries for the next pdlarfg steps are
    ! computed
    subroutine qr_pdlarfgk_1dcomm_update_&
&double &
(obj,a,lda,baseidx,work,lwork,seedC,seedD,k,rank,sidx,tau,n,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter :: gmode_ = 1,rank_ = 2,eps_ = 3, upmode1_ = 4

      ! input variables (local)
      integer(kind=ik)            :: lda,lwork
      real(kind=c_double)               :: a(lda,*),work(*)

      ! input variables (global)
      integer(kind=ik)            :: k,rank,sidx,n,baseidx,nb,rev,mpicomm
      real(kind=c_double)               :: beta

      ! output variables (global)
      real(kind=c_double)               :: seedC(k,*),seedD(k,*),tau(*)

      ! derived input variables from QR_PQRPARAM

      ! local scalars
      real(kind=c_double)               :: alpha
      integer(kind=ik)            :: coffset,zoffset,yoffset,voffset,buffersize
      integer(kind=ik)            :: mpirank,mpiprocs,mpirank_top
      integer(kind=MPI_KIND)      :: mpirankMPI, mpierr,mpiprocsMPI 

      integer(kind=ik)            :: localsize,baseoffset,localoffset,topidx
      integer(kind=ik)            :: lidx
        call obj%timer%start("qr_pdlarfgk_1dcomm_update_&
            &double&
            &")
      if (lwork .eq. -1) then
        ! buffer for c,z,y,v
        work(1) = 4*k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &double&
            &")

        return
      end if

      ! nothing to update anymore
      if (sidx .gt. rank) then
        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &double&
            &")
        return
      endif
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)

      lidx = baseidx-sidx
      if (lidx .lt. 1) then
        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &double&
            &")
        return
      endif

      call local_size_offset_1d(n,nb,baseidx,lidx,rev,mpirank,mpiprocs, &
                                localsize,baseoffset,localoffset)

      coffset = 1
      zoffset = coffset + k
      yoffset = zoffset + k
      voffset = yoffset + k
      buffersize = k - sidx

      ! finalize tau values
      alpha = seedC(k-sidx+1,k-sidx+1)
      beta = tau(k-sidx+1)

      ! zero Householder Vector (zero norm) case
      !print *,'k update: alpha,beta',alpha,beta
      if ((beta .eq. 0.0_rk8) .or. (alpha .eq. 0.0_rk8))  then
        tau(k-sidx+1) = 0.0_rk8
        seedC(k,k-sidx+1) = 0.0_rk8

        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &double&
            &")
        return
      end if

      tau(k-sidx+1) = (beta+alpha) / beta

      ! ---------------------------------------
      ! calculate c Vector (extra Vector or encode in seedC/seedD?
      work(coffset:coffset+buffersize-1) = seedD(1:buffersize,k-sidx+1)
      call dgemv("Trans", buffersize+1, buffersize, &
                 1.0_rk8,seedC(1,1),k,seedC(1,k-sidx+1),1, &
                 1.0_rk8,work(coffset),1)

      ! calculate z using tau,seedD,seedC and c Vector
      work(zoffset:zoffset+buffersize-1) = seedC(k-sidx+1,1:buffersize)
      call daxpy(buffersize, 1.0_rk8/beta, work(coffset), 1, work(zoffset), 1)

      ! update A1(local copy) and generate part of householder vectors for use
      call daxpy(buffersize, -1.0_rk8, work(zoffset),1,seedC(k-sidx+1,1),k)
      call dscal(buffersize, 1.0_rk8/(alpha+beta), seedC(1,k-sidx+1),1)
      call dger(buffersize, buffersize, -1.0_rk8, seedC(1,k-sidx+1),1, work(zoffset), 1, seedC(1,1), k)

      ! update A global (householder Vector already generated by pdlarfgk)
      mpirank_top = MOD(lidx/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        ! handle first row separately
        topidx = local_index(lidx+1,mpirank_top,mpiprocs,nb,0)
        call daxpy(buffersize,-1.0_rk8,work(zoffset),1,a(topidx,1),lda)
      end if

      call dger(localsize, buffersize,-1.0_rk8, &
                a(localoffset,k-sidx+1),1,work(zoffset),1, &
                a(localoffset,1),lda)

      ! update D (symmetric) => two buffer vectors of size rank
      ! generate y Vector
      work(yoffset:yoffset+buffersize-1) = 0._rk8
      call daxpy(buffersize,1.0_rk8/(alpha+beta),work(zoffset),1,work(yoffset),1)

      ! generate v Vector
      work(voffset:voffset+buffersize-1) = seedD(1:buffersize,k-sidx+1)
      call daxpy(buffersize, -0.5_rk8*seedD(k-sidx+1,k-sidx+1), work(yoffset), 1, work(voffset),1)

      ! symmetric update of D using y and v
      call dsyr2("Upper", buffersize,-1.0_rk8, &
                     work(yoffset),1,work(voffset),1, &
                     seedD(1,1), k)

      ! prepare T matrix inner products
      ! D_k(1:k,k+1:n) = D_(k-1)(1:k,k+1:n) - D_(k-1)(1:k,k) * y'
      ! store coefficient 1.0d0/(alpha+beta) in C diagonal elements
      call dger(k-sidx,sidx,-1.0_rk8,work(yoffset),1,seedD(k-sidx+1,k-sidx+1),k,seedD(1,k-sidx+1),k)
      seedC(k,k-sidx+1) = 1.0_rk8/(alpha+beta)

        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &double&
            &")
    end subroutine

    subroutine qr_pdlarfgk_1dcomm_generateT_&
          &double &
          (obj,seedC,seedD,k,actualk,tau,t,ldt)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)  :: k,actualk,ldt
      real(kind=c_double)     :: seedC(k,*),seedD(k,*),tau(*),t(ldt,*)

      integer(kind=ik)  :: irow,icol
      real(kind=c_double)     :: column_coefficient
        call obj%timer%start("qr_pdlarfgk_1dcomm_generateT_&
            &double&
            &")

      !print *,'reversed on the fly T generation NYI'

      do icol=1,actualk-1
        ! calculate inner product of householder Vector parts in seedC
        ! (actually calculating more than necessary, if actualk < k)
        ! => a lot of junk from row 1 to row k-actualk
        call dtrmv('Upper','Trans','Unit',k-icol,seedC(1,1),k,seedC(1,k-icol+1),1)
        ! add scaled D parts to current column of C (will become later T rows)
        column_coefficient = seedC(k,k-icol+1)
        do irow=k-actualk+1,k-1
          seedC(irow,k-icol+1) = ( seedC(irow,k-icol+1) ) +  ( seedD(irow,k-icol+1) * column_coefficient * seedC(k,irow) )
        end do
      end do

      call qr_dlarft_kernel_&
             &double &
             (actualk,tau(k-actualk+1),seedC(k-actualk+1,k-actualk+2),k,t(k-actualk+1,k-actualk+1),ldt)
      call obj%timer%stop("qr_pdlarfgk_1dcomm_generateT_&
             &double&
             &")

    end subroutine

    !direction=0: pack into work buffer
    !direction=1: unpack from work buffer
    subroutine qr_pdgeqrf_pack_unpack_&
&double &
(obj,v,ldv,work,lwork,m,n,mb,baseidx,rowidx,rev,direction,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)   :: ldv,lwork
      real(kind=c_double)      :: v(ldv,*), work(*)

      ! input variables (global)
      integer(kind=ik)   :: m,n,mb,baseidx,rowidx,rev,direction,mpicomm

      ! output variables (global)

      ! local scalars
      integer(kind=ik)   :: mpirank,mpiprocs
      integer(kind=MPI_KIND) :: mpierr,mpirankMPI, mpiprocsMPI

      integer(kind=ik)   :: buffersize,icol
      integer(kind=ik)   :: local_size,baseoffset,offset

      ! external functions
        call obj%timer%start("qr_pdgeqrf_pack_unpack_&
            &double&
            &")
      call mpi_comm_rank(int(mpicomm,kind=MPI_KIND) ,mpirankMPI,mpierr)
      call mpi_comm_size(int(mpicomm,kind=MPI_KIND) ,mpiprocsMPI,mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)

      call local_size_offset_1d(m,mb,baseidx,rowidx,rev,mpirank,mpiprocs, &
                                    local_size,baseoffset,offset)

      !print *,'pack/unpack',local_size,baseoffset,offset

      ! rough approximate for buffer size
      if (lwork .eq. -1) then
        buffersize = local_size * n ! Vector elements
        work(1) = DBLE(buffersize)
        call obj%timer%stop("qr_pdgeqrf_pack_unpack_&
            &double&
            &")

        return
      end if

      if (direction .eq. 0) then
        ! copy v part to buffer (including zeros)
        do icol=1,n
          work(1+local_size*(icol-1):local_size*icol) = v(baseoffset:baseoffset+local_size-1,icol)
        end do
      else
        ! copy v part from buffer (including zeros)
        do icol=1,n
          v(baseoffset:baseoffset+local_size-1,icol) = work(1+local_size*(icol-1):local_size*icol)
        end do
      end if
        call obj%timer%stop("qr_pdgeqrf_pack_unpack_&
            &double&
            &")

      return

    end subroutine

    !direction=0: pack into work buffer
    !direction=1: unpack from work buffer
    subroutine qr_pdgeqrf_pack_unpack_tmatrix_&
          &double &
          (obj,tau,t,ldt,work,lwork,n,direction)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)  :: ldt,lwork
      real(kind=c_double)     :: work(*), t(ldt,*),tau(*)

      ! input variables (global)
      integer(kind=ik)  :: n,direction

      ! output variables (global)

      ! local scalars
      integer(kind=ik)  :: icol

      ! external functions
        call obj%timer%start("qr_pdgeqrf_pack_unpack_tmatrix_&
            &double&
            &")

      if (lwork .eq. -1) then
        work(1) = real(n*n,kind=c_double)

        call obj%timer%stop("qr_pdgeqrf_pack_unpack_tmatrix_&
            &double&
            &")
        return
      end if

      if (direction .eq. 0) then
        ! append t matrix to buffer (including zeros)
        do icol=1,n
          work(1+(icol-1)*n:icol*n) = t(1:n,icol)
        end do
      else
        ! append t matrix from buffer (including zeros)
        do icol=1,n
          t(1:n,icol) = work(1+(icol-1)*n:icol*n)
          tau(icol) = t(icol,icol)
        end do
      end if
        call obj%timer%stop("qr_pdgeqrf_pack_unpack_tmatrix_&
            &double&
            &")
    end subroutine


    ! TODO: encode following functionality
    !   - Direction? BOTTOM UP or TOP DOWN ("Up", "Down")
    !        => influences all related kernels (including DLARFT / DLARFB)
    !   - rank-k parameter (k=1,2,...,b)
    !        => influences possible update strategies
    !        => parameterize the function itself? (FUNCPTR, FUNCARG)
    !   - Norm mode? Allreduce, Allgather, AlltoAll, "AllHouse", (ALLNULL = benchmarking local kernels)
    !   - subblocking
    !         (maximum block size bounded by data distribution along rows)
    !   - blocking method (householder vectors only or compact WY?)
    !   - update strategy of trailing parts (incremental, complete)
    !        - difference for subblocks and normal blocks? (UPDATE and UPDATESUB)
    !        o "Incremental"
    !        o "Full"
    !   - final T generation (recursive: subblock wise, block wise, end) (TMERGE)
    !        ' (implicitly given by / influences update strategies?)
    !        => alternative: during update: iterate over sub t parts
    !           => advantage: smaller (cache aware T parts)
    !           => disadvantage: more memory write backs
    !                (number of T parts * matrix elements)
    !   - partial/sub T generation (TGEN)
    !        o add vectors right after creation (Vector)
    !        o add set of vectors (Set)
    !   - bcast strategy of householder vectors to other process columns
    !        (influences T matrix generation and trailing update
    !         in other process columns)
    !        o no broadcast (NONE = benchmarking?,
    !            or not needed due to 1D process grid)
    !        o after every housegen (VECTOR)
    !        o after every subblk   (SUBBLOCK)
    !        o after full local column block decomposition (BLOCK)
    !  LOOP Housegen -> BCAST -> GENT/EXTENDT -> LOOP HouseLeft

    !subroutine qr_pqrparam_init(PQRPARAM, DIRECTION, RANK, NORMMODE, &
    !                             SUBBLK, UPDATE, TGEN, BCAST)
    ! gmode: control communication pattern of dlarfg
    ! maxrank: control max number of householder vectors per communication
    ! eps: error threshold (integer)
    ! update*: control update pattern in pdgeqr2_1dcomm ('incremental','full','merge')
    !               merging = full update with tmatrix merging
    ! tmerge*: 0: do not merge, 1: incremental merge, >1: recursive merge
    !               only matters if update* == full
    subroutine qr_pqrparam_init(obj,pqrparam,size2d,update2d,tmerge2d,size1d,update1d,tmerge1d,maxrank,update,eps,hgmode)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input
      CHARACTER         :: update2d,update1d,update,hgmode
      INTEGER(kind=ik)  :: size2d,size1d,maxrank,eps,tmerge2d,tmerge1d

      ! output
      INTEGER(kind=ik)  :: PQRPARAM(1:11)

        call obj%timer%start("qr_pqrparam_init")

      PQRPARAM(1) = size2d
      PQRPARAM(2) = ichar(update2d)
      PQRPARAM(3) = tmerge2d
      ! TODO: broadcast T yes/no

      PQRPARAM(4) = size1d
      PQRPARAM(5) = ichar(update1d)
      PQRPARAM(6) = tmerge1d

      PQRPARAM(7) = maxrank
      PQRPARAM(8) = ichar(update)
      PQRPARAM(9) = eps
      PQRPARAM(10) = ichar(hgmode)
        call obj%timer%stop("qr_pqrparam_init")

    end subroutine qr_pqrparam_init

    subroutine qr_pdlarfg_copy_1dcomm_&
&double &
(obj,x,incx,v,incv,n,baseidx,idx,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)  :: incx,incv
      real(kind=c_double)     :: x(*), v(*)

      ! input variables (global)
      integer(kind=ik)  :: baseidx,idx,rev,nb,n
      integer(kind=ik)  :: mpicomm

      ! output variables (global)

      ! local scalars
      integer(kind=ik)  :: mpiprocs
      integer(kind=MPI_KIND) ::  mpierr,mpiprocsMPI,mpirankMPI

      integer(kind=ik)  :: mpirank,mpirank_top
      integer(kind=ik)  :: irow,x_offset
      integer(kind=ik)  :: v_offset,local_size

        call obj%timer%start("qr_pdlarfg_copy_1dcomm_&
            &double&
            &")
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      call local_size_offset_1d(n,nb,baseidx,idx,rev,mpirank,mpiprocs, &
                                local_size,v_offset,x_offset)
      v_offset = v_offset * incv

      !print *,'copy:',mpirank,baseidx,v_offset,x_offset,local_size

      ! copy elements
      do irow=1,local_size
        v((irow-1)*incv+v_offset) = x((irow-1)*incx+x_offset)
      end do

      ! replace top element to build an unitary Vector
      mpirank_top = MOD((idx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        v(local_size*incv) = 1.0_rk8
      end if
        call obj%timer%stop("qr_pdlarfg_copy_1dcomm_&
            &double&
            &")

    end subroutine

! vim: syntax=fortran

  ! real single precision


















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
     subroutine qr_pdgeqrf_2dcomm_&
           &single &
           (obj, a, lda, matrixCols, v, ldv, vmrCols, tau, lengthTau, t, ldt, colsT, &
                                  work, workLength, lwork, m, n, mb, nb, rowidx, colidx, &
                                  rev, trans, PQRPARAM, mpicomm_rows, mpicomm_cols, blockheuristic)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter   :: gmode_ = 1, rank_ = 2, eps_ = 3

      ! input variables (local)
      integer(kind=ik), intent(in)  :: lda, lwork, ldv, ldt, matrixCols, m, vmrCols, lengthTau, &
                                       colsT, workLength

      ! input variables (global)
      integer(kind=ik)              :: n, mb, nb, rowidx, colidx, rev, trans, mpicomm_cols, mpicomm_rows
      integer(kind=ik)              :: PQRPARAM(1:11)
      real(kind=c_float)                 :: a(1:lda,1:matrixCols), v(1:ldv,1:vmrCols), tau(1:lengthTau), &
                                       t(1:ldt,1:colsT), work(1:workLength)
      ! output variables (global)
      real(kind=c_float)                 :: blockheuristic(*)

      ! input variables derived from PQRPARAM
      integer(kind=ik)              :: updatemode,tmerge,size2d

      ! local scalars
      integer(kind=ik)              :: mpirank_cols,broadcast_size,mpirank_rows
      integer(kind=MPI_KIND)        :: mpirank_colsMPI, mpirank_rowsMPI
      integer(kind=MPI_KIND)        :: mpierr
      integer(kind=ik)              :: mpirank_cols_qr,mpiprocs_cols
      integer(kind=MPI_KIND)        :: mpiprocs_colsMPI
      integer(kind=ik)              :: lcols_temp,lcols,icol,lastcol
      integer(kind=ik)              :: baseoffset,offset,idx,voffset
      integer(kind=ik)              :: update_voffset,update_tauoffset
      integer(kind=ik)              :: update_lcols
      integer(kind=ik)              :: work_offset

      real(kind=c_float)                 :: dbroadcast_size(1),dtmat_bcast_size(1)
      real(kind=c_float)                 :: pdgeqrf_size(1),pdlarft_size(1),pdlarfb_size(1),tmerge_pdlarfb_size(1)
      integer(kind=ik)              :: temptau_offset,temptau_size,broadcast_offset,tmat_bcast_size
      integer(kind=ik)              :: remaining_cols
      integer(kind=ik)              :: total_cols
      integer(kind=ik)              :: incremental_update_size ! needed for incremental update mode

      call obj%timer%start("qr_pdgeqrf_2dcomm_&
          &single&
          &")
      size2d     = PQRPARAM(1)
      updatemode = PQRPARAM(2)
      tmerge     = PQRPARAM(3)

      ! copy value before we are going to filter it
      total_cols = n
      call mpi_comm_rank(int(mpicomm_cols,kind=MPI_KIND) ,mpirank_colsMPI, mpierr)
      call mpi_comm_rank(int(mpicomm_rows,kind=MPI_KIND) ,mpirank_rowsMPI, mpierr)
      call mpi_comm_size(int(mpicomm_cols,kind=MPI_KIND) ,mpiprocs_colsMPI, mpierr)

      mpirank_cols = int(mpirank_colsMPI,kind=c_int)
      mpirank_rows = int(mpirank_rowsMPI,kind=c_int)
      mpiprocs_cols = int(mpiprocs_colsMPI,kind=c_int)

      call qr_pdgeqrf_1dcomm_&
          &single &
          (obj,a,lda,v,ldv,tau,t,ldt,pdgeqrf_size(1),-1,m,total_cols,mb,rowidx,rowidx,rev,trans, &
                             PQRPARAM(4:11),mpicomm_rows,blockheuristic)
      call qr_pdgeqrf_pack_unpack_&
          &single &
          (obj,v,ldv,dbroadcast_size(1),-1,m,total_cols,mb,rowidx,rowidx,rev,0,mpicomm_rows)
      call qr_pdgeqrf_pack_unpack_tmatrix_&
          &single &
          (obj,tau,t,ldt,dtmat_bcast_size(1),-1,total_cols,0)

      pdlarft_size(1) = 0.0_rk4

      call qr_pdlarfb_1dcomm_&
          &single &
          (m,mb,total_cols,total_cols,a,lda,v,ldv,tau,t,ldt,rowidx,rowidx,rev,mpicomm_rows, &
                             pdlarfb_size(1),-1)
      call qr_tmerge_pdlarfb_1dcomm_&
          &single &
          (m,mb,total_cols,total_cols,total_cols,v,ldv,t,ldt,a,lda,rowidx,rev,updatemode, &
                                    mpicomm_rows,tmerge_pdlarfb_size(1),-1)


      temptau_offset = 1
      temptau_size = total_cols
      broadcast_offset = temptau_offset + temptau_size
      broadcast_size = int(dbroadcast_size(1) + dtmat_bcast_size(1))
      work_offset = broadcast_offset + broadcast_size

      if (lwork .eq. -1) then
        work(1) = (real(temptau_size,kind=rk4) + real(broadcast_size,kind=rk4) + max(pdgeqrf_size(1), &
            pdlarft_size(1),pdlarfb_size(1), &
                   tmerge_pdlarfb_size(1)))
        call obj%timer%stop("qr_pdgeqrf_2dcomm_&
          &single&
          &")
        return
      end if

      lastcol = colidx-total_cols+1
      voffset = total_cols

      incremental_update_size = 0

      ! clear v buffer: just ensure that there is no junk in the upper triangle
      ! part, otherwise pdlarfb gets some problems
      ! pdlarfl(2) do not have these problems as they are working more on a Vector
      ! basis
      v(1:ldv,1:total_cols) = 0.0_rk4
      icol = colidx

      remaining_cols = total_cols

      !print *,'start decomposition',m,rowidx,colidx

      do while (remaining_cols .gt. 0)

        ! determine rank of process column with next qr block
        mpirank_cols_qr = MOD((icol-1)/nb,mpiprocs_cols)

        ! lcols can't be larger than than nb
        ! exception: there is only one process column

        ! however, we might not start at the first local column.
        ! therefore assume a matrix of size (1xlcols) starting at (1,icol)
        ! determine the real amount of local columns
        lcols_temp = min(nb,(icol-lastcol+1))

        ! blocking parameter
        lcols_temp = max(min(lcols_temp,size2d),1)

        ! determine size from last decomposition column
        !  to first decomposition column
        call local_size_offset_1d(icol,nb,icol-lcols_temp+1,icol-lcols_temp+1,0, &
                                      mpirank_cols_qr,mpiprocs_cols, &
                                      lcols,baseoffset,offset)

        voffset = remaining_cols - lcols + 1

        idx = rowidx - colidx + icol

        if (mpirank_cols .eq. mpirank_cols_qr) then
          ! qr decomposition part
          tau(offset:offset+lcols-1) = 0.0_rk4

          call qr_pdgeqrf_1dcomm_&
              &single &
              (obj,a(1,offset),lda,v(1,voffset),ldv,tau(offset),t(voffset,voffset),ldt, &
                                 work(work_offset),lwork,m,lcols,mb,rowidx,idx,rev,trans,PQRPARAM(4:11), &
                                 mpicomm_rows,blockheuristic)

          ! pack broadcast buffer (v + tau)
          call qr_pdgeqrf_pack_unpack_&
              &single &
              (obj,v(1,voffset),ldv,work(broadcast_offset),lwork,m,lcols,mb,rowidx,&
                                      idx,rev,0,mpicomm_rows)

          ! determine broadcast size
          call qr_pdgeqrf_pack_unpack_&
              &single &
              (obj,v(1,voffset),ldv,dbroadcast_size(1),-1,m,lcols,mb,rowidx,idx,rev,&
                                      0,mpicomm_rows)
          broadcast_size = int(dbroadcast_size(1))

          !if (mpirank_rows .eq. 0) then
          ! pack tmatrix into broadcast buffer and calculate new size
          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &single &
              (obj,tau(offset),t(voffset,voffset),ldt, &
                                              work(broadcast_offset+broadcast_size),lwork,lcols,0)
          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &single &
              (obj,tau(offset),t(voffset,voffset),ldt,dtmat_bcast_size(1),-1,lcols,0)
          broadcast_size = broadcast_size + int(dtmat_bcast_size(1))
          !end if

          ! initiate broadcast (send part)

          call MPI_Bcast(work(broadcast_offset), int(broadcast_size,kind=MPI_KIND), mpi_real4, &
                         int(mpirank_cols_qr,kind=MPI_KIND), int(mpicomm_cols,kind=MPI_KIND), mpierr)

          ! copy tau parts into temporary tau buffer
          work(temptau_offset+voffset-1:temptau_offset+(voffset-1)+lcols-1) = tau(offset:offset+lcols-1)

          !print *,'generated tau:', tau(offset)
        else
          ! Vector exchange part

          ! determine broadcast size
          call qr_pdgeqrf_pack_unpack_&
              &single &
              (obj,v(1,voffset),ldv,dbroadcast_size(1),-1,m,lcols,mb,rowidx,idx,rev,1,mpicomm_rows)
          broadcast_size = int(dbroadcast_size(1))

          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &single &
              (obj,work(temptau_offset+voffset-1),t(voffset,voffset),ldt, &
                                              dtmat_bcast_size(1),-1,lcols,0)
          tmat_bcast_size = dtmat_bcast_size(1)

          !print *,'broadcast_size (nonqr)',broadcast_size
          broadcast_size = dbroadcast_size(1) + dtmat_bcast_size(1)

          ! initiate broadcast (recv part)

          call MPI_Bcast(work(broadcast_offset), int(broadcast_size,kind=MPI_KIND), mpi_real4, &
                         int(mpirank_cols_qr,kind=MPI_KIND), int(mpicomm_cols,kind=MPI_KIND), mpierr)

          ! last n*n elements in buffer are (still empty) T matrix elements
          ! fetch from first process in each column

          ! unpack broadcast buffer (v + tau)
          call qr_pdgeqrf_pack_unpack_&
              &single &
              (obj,v(1,voffset),ldv,work(broadcast_offset),lwork,m,lcols, &
              mb,rowidx,idx,rev,1,mpicomm_rows)

          ! now send t matrix to other processes in our process column
          broadcast_size = int(dbroadcast_size(1))
          tmat_bcast_size = int(dtmat_bcast_size(1))

          ! t matrix should now be available on all processes => unpack
          call qr_pdgeqrf_pack_unpack_tmatrix_&
              &single &
              (obj,work(temptau_offset+voffset-1),t(voffset,voffset),ldt, &
                                              work(broadcast_offset+broadcast_size),lwork,lcols,1)
        end if

        remaining_cols = remaining_cols - lcols

        ! apply householder vectors to whole trailing matrix parts (if any)

        update_voffset = voffset
        update_tauoffset = icol
        update_lcols = lcols
        incremental_update_size = incremental_update_size + lcols

        icol = icol - lcols
        ! count colums from first column of global block to current index
        call local_size_offset_1d(icol,nb,colidx-n+1,colidx-n+1,0, &
                                      mpirank_cols,mpiprocs_cols, &
                                      lcols,baseoffset,offset)

        if (lcols .gt. 0) then

          !print *,'updating trailing matrix'

           if (updatemode .eq. ichar('I')) then
             print *,'pdgeqrf_2dcomm: incremental update not yet implemented! rev=1'
           else if (updatemode .eq. ichar('F')) then
             ! full update no merging
             call qr_pdlarfb_1dcomm_&
                &single &
                (m,mb,lcols,update_lcols,a(1,offset),lda,v(1,update_voffset),ldv, &
                     work(temptau_offset+update_voffset-1),                          &
                                                        t(update_voffset,update_voffset),ldt, &
                    rowidx,idx,1,mpicomm_rows,work(work_offset),lwork)
           else
            ! full update + merging default
             call qr_tmerge_pdlarfb_1dcomm_&
                &single &
                (m,mb,lcols,n-(update_voffset+update_lcols-1),update_lcols, &
                                                              v(1,update_voffset),ldv, &
                          t(update_voffset,update_voffset),ldt, &
                          a(1,offset),lda,rowidx,1,updatemode,mpicomm_rows, &
                                                              work(work_offset),lwork)
           end if
        else
           if (updatemode .eq. ichar('I')) then
             !print *,'sole merging of (incremental) T matrix', mpirank_cols,  &
            !                            n-(update_voffset+incremental_update_size-1)
             call qr_tmerge_pdlarfb_1dcomm_&
                &single &
                (m,mb,0,n-(update_voffset+incremental_update_size-1),   &
                                                              incremental_update_size,v(1,update_voffset),ldv, &
                          t(update_voffset,update_voffset),ldt, &
                          a,lda,rowidx,1,updatemode,mpicomm_rows,work(work_offset),lwork)

             ! reset for upcoming incremental updates
             incremental_update_size = 0
          else if (updatemode .eq. ichar('M')) then
             ! final merge
            call qr_tmerge_pdlarfb_1dcomm_&
                &single &
                (m,mb,0,n-(update_voffset+update_lcols-1),update_lcols, &
                                                              v(1,update_voffset),ldv, &
                          t(update_voffset,update_voffset),ldt, &
                          a,lda,rowidx,1,updatemode,mpicomm_rows,work(work_offset),lwork)
          else
            ! full updatemode - nothing to update
          end if

          ! reset for upcoming incremental updates
          incremental_update_size = 0
        end if
      end do

      if ((tmerge .gt. 0) .and. (updatemode .eq. ichar('F'))) then
        ! finally merge all small T parts
        call qr_pdlarft_tree_merge_1dcomm_&
&single &
(m,mb,n,size2d,tmerge,v,ldv,t,ldt,rowidx,rev,mpicomm_rows,work,lwork)
      end if

      !print *,'stop decomposition',rowidx,colidx
      call obj%timer%stop("qr_pdgeqrf_2dcomm_&
          &single&
          &")
    end subroutine

    subroutine qr_pdgeqrf_1dcomm_&
&single &
(obj,a,lda,v,ldv,tau,t,ldt,work,lwork,m,n,mb,baseidx,rowidx,rev,trans, &
          PQRPARAM,mpicomm,blockheuristic)
      use precision
      use elpa1_impl
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter  :: gmode_ = 1,rank_ = 2,eps_ = 3

      ! input variables (local)
      integer(kind=ik)             :: lda,lwork,ldv,ldt
      real(kind=c_float)                :: a(lda,*),v(ldv,*),tau(*),t(ldt,*),work(*)

      ! input variables (global)
      integer(kind=ik)             :: m,n,mb,baseidx,rowidx,rev,trans,mpicomm
      integer(kind=ik)             :: PQRPARAM(:)
      ! derived input variables

      ! derived further input variables from QR_PQRPARAM
      integer(kind=ik)             :: size1d,updatemode,tmerge

      ! output variables (global)
      real(kind=c_float)                :: blockheuristic(*)

      ! local scalars
      integer(kind=ik)             :: nr_blocks,remainder,current_block,aoffset,idx,updatesize
      real(kind=c_float)                :: pdgeqr2_size(1),pdlarfb_size(1),tmerge_tree_size(1)
      call obj%timer%start("qr_pdgeqrf_1dcomm_&
          &single&
          &")
      size1d     = max(min(PQRPARAM(1),n),1)
      updatemode = PQRPARAM(2)
      tmerge     = PQRPARAM(3)

      if (lwork .eq. -1) then
        call qr_pdgeqr2_1dcomm_&
&single &
(obj,a,lda,v,ldv,tau,t,ldt,pdgeqr2_size,-1, &
                                  m,size1d,mb,baseidx,baseidx,rev,trans,PQRPARAM(4:),mpicomm,blockheuristic)
        ! reserve more space for incremental mode
        call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,n,n,n,v,ldv,t,ldt, &
                                         a,lda,baseidx,rev,updatemode,mpicomm,pdlarfb_size,-1)

        call qr_pdlarft_tree_merge_1dcomm_&
&single &
(m,mb,n,size1d,tmerge,v,ldv,t,ldt,baseidx,rev,mpicomm,tmerge_tree_size,-1)

        work(1) = max(pdlarfb_size(1),pdgeqr2_size(1),tmerge_tree_size(1))
      call obj%timer%stop("qr_pdgeqrf_1dcomm_&
          &single&
          &")
        return
      end if

      nr_blocks = n / size1d
      remainder = n - nr_blocks*size1d

      current_block = 0
      do while (current_block .lt. nr_blocks)
        idx = rowidx-current_block*size1d
        updatesize = n-(current_block+1)*size1d
        aoffset = 1+updatesize
        call qr_pdgeqr2_1dcomm_&
&single &
(obj,a(1,aoffset),lda,v(1,aoffset),ldv,tau(aoffset),t(aoffset,aoffset),ldt,work,lwork, &
                                m,size1d,mb,baseidx,idx,1,trans,PQRPARAM(4:),mpicomm,blockheuristic)
        if (updatemode .eq. ichar('M')) then
          ! full update + merging
          call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,updatesize,current_block*size1d,size1d, &
                                           v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                           a,lda,baseidx,1,ichar('F'),mpicomm,work,lwork)
        else if (updatemode .eq. ichar('I')) then
          if (updatesize .ge. size1d) then
            ! incremental update + merging
            call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,size1d,current_block*size1d,size1d, &
                                               v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                               a(1,aoffset-size1d),lda,baseidx,1,updatemode,mpicomm,work,lwork)

          else ! only remainder left
             ! incremental update + merging
             call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,remainder,current_block*size1d,size1d, &
                                               v(1,aoffset),ldv,t(aoffset,aoffset),ldt, &
                                               a(1,1),lda,baseidx,1,updatemode,mpicomm,work,lwork)
          end if
        else ! full update no merging is default
          ! full update no merging
          call qr_pdlarfb_1dcomm_&
&single &
(m,mb,updatesize,size1d,a,lda,v(1,aoffset),ldv, &
                                    tau(aoffset),t(aoffset,aoffset),ldt,baseidx,idx,1,mpicomm,work,lwork)
        end if

        ! move on to next block
        current_block = current_block+1
      end do

      if (remainder .gt. 0) then
        aoffset = 1
        idx = rowidx-size1d*nr_blocks
        call qr_pdgeqr2_1dcomm_&
&single &
(obj,a(1,aoffset),lda,v,ldv,tau,t,ldt,work,lwork, &
                                  m,remainder,mb,baseidx,idx,1,trans,PQRPARAM(4:),mpicomm,blockheuristic)
        if ((updatemode .eq. ichar('I')) .or. (updatemode .eq. ichar('M'))) then
          ! final merging
          call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,0,size1d*nr_blocks,remainder, &
                                             v,ldv,t,ldt, &
                                             a,lda,baseidx,1,updatemode,mpicomm,work,lwork) ! updatemode argument does not matter
        end if
      end if

      if ((tmerge .gt. 0) .and. (updatemode .eq. ichar('F'))) then
        ! finally merge all small T parts
        call qr_pdlarft_tree_merge_1dcomm_&
&single &
(m,mb,n,size1d,tmerge,v,ldv,t,ldt,baseidx,rev,mpicomm,work,lwork)
      end if
      call obj%timer%stop("qr_pdgeqrf_1dcomm_&
          &single&
          &")

    end subroutine

    ! local a and tau are assumed to be positioned at the right column from a local
    ! perspective
    ! TODO: if local amount of data turns to zero the algorithm might produce wrong
    ! results (probably due to old buffer contents)
    subroutine qr_pdgeqr2_1dcomm_&
&single &
(obj,a,lda,v,ldv,tau,t,ldt,work,lwork,m,n,mb,baseidx,rowidx,rev, &
          trans,PQRPARAM,mpicomm,blockheuristic)
      use precision
      !use elpa1_impl ! check this
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter   :: gmode_ = 1,rank_ = 2 ,eps_ = 3, upmode1_ = 4

      ! input variables (local)
      integer(kind=ik)              :: lda,lwork,ldv,ldt
      real(kind=c_float)                 :: a(lda,*),v(ldv,*),tau(*),t(ldt,*),work(*)

      ! input variables (global)
      integer(kind=ik)              :: m,n,mb,baseidx,rowidx,rev,trans,mpicomm
      integer(kind=ik)              :: PQRPARAM(:)
      ! output variables (global)
      real(kind=c_float)                 :: blockheuristic(*)

      ! derived further input variables from QR_PQRPARAM
      integer(kind=ik)              ::  maxrank,hgmode,updatemode

      ! local scalars
      integer(kind=ik)              :: icol,incx,idx
      real(kind=c_float)                 :: pdlarfg_size(1),pdlarf_size(1),total_size
      real(kind=c_float)                 :: pdlarfg2_size(1),pdlarfgk_size(1),pdlarfl2_size(1)
      real(kind=c_float)                 :: pdlarft_size(1),pdlarfb_size(1),pdlarft_pdlarfb_size(1),tmerge_pdlarfb_size(1)
      integer(kind=ik)              :: mpirank,mpiprocs
      integer(kind=MPI_KIND)        :: mpirankMPI, mpiprocsMPI
      integer(kind=MPI_KIND)        :: mpierr
      integer(kind=ik)              :: rank,lastcol,actualrank,nextrank
      integer(kind=ik)              :: update_cols,decomposition_cols
      integer(kind=ik)              :: current_column
      call obj%timer%start("qr_pdgeqr2_1dcomm_&
          &single&
          &")

      maxrank    = min(PQRPARAM(1),n)
      updatemode = PQRPARAM(2)
      hgmode     = PQRPARAM(4)
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      if (trans .eq. 1) then
        incx = lda
      else
        incx = 1
      end if

      if (lwork .eq. -1) then

        call qr_pdlarfg_1dcomm_&
&single &
(obj,a,incx,tau(1),pdlarfg_size(1),-1,n,rowidx,mb,hgmode,rev,mpicomm)
        call qr_pdlarfl_1dcomm_&
&single &
(v,1,baseidx,a,lda,tau(1),pdlarf_size(1),-1,m,n,rowidx,mb,rev,mpicomm)
        call qr_pdlarfg2_1dcomm_ref_&
&single &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,pdlarfg2_size(1),-1,m,rowidx,mb,PQRPARAM(:), &
                                    rev,mpicomm,actualrank)

        call qr_pdlarfgk_1dcomm_&
&single &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,pdlarfgk_size(1),-1,m,n, &
            rowidx,mb,PQRPARAM(:),rev,mpicomm,actualrank)
        call qr_pdlarfl2_tmatrix_1dcomm_&
&single &
(v,ldv,baseidx,a,lda,t,ldt,pdlarfl2_size(1),-1,m,n,rowidx,mb,rev,mpicomm)
        pdlarft_size(1) = 0.0_rk4
        call qr_pdlarfb_1dcomm_&
&single &
(m,mb,n,n,a,lda,v,ldv,tau,t,ldt,baseidx,rowidx,1,mpicomm,pdlarfb_size(1),-1)
        pdlarft_pdlarfb_size(1) = 0.0_rk4
        call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,n,n,n,v,ldv,t,ldt,a,lda,rowidx,rev,&
            updatemode,mpicomm,tmerge_pdlarfb_size(1),-1)

        total_size = max(pdlarfg_size(1),pdlarf_size(1),pdlarfg2_size(1),pdlarfgk_size(1),pdlarfl2_size(1),pdlarft_size(1), &
                         pdlarfb_size(1),pdlarft_pdlarfb_size(1),tmerge_pdlarfb_size(1))

        work(1) = total_size
      call obj%timer%stop("qr_pdgeqr2_1dcomm_&
          &single&
          &")
        return
      end if

      icol = 1
      lastcol = min(rowidx,n)
      decomposition_cols = lastcol
      update_cols = n
      do while (decomposition_cols .gt. 0) ! local qr block
        icol = lastcol-decomposition_cols+1
        idx = rowidx-icol+1

        ! get possible rank size
        ! limited by number of columns and remaining rows
        rank = min(n-icol+1,maxrank,idx)

        current_column = n-icol+1-rank+1

        if (rank .eq. 1) then

          call qr_pdlarfg_1dcomm_&
&single &
(obj,a(1,current_column),incx, &
                                  tau(current_column),work,lwork, &
                                  m,idx,mb,hgmode,1,mpicomm)
          v(1:ldv,current_column) = 0.0_rk4
          call qr_pdlarfg_copy_1dcomm_&
&single &
(obj,a(1,current_column),incx, &
                                       v(1,current_column),1, &
                                       m,baseidx,idx,mb,1,mpicomm)

          ! initialize t matrix part
          t(current_column,current_column) = tau(current_column)

          actualrank = 1

        else if (rank .eq. 2) then
          call qr_pdlarfg2_1dcomm_ref_&
&single &
(obj,a(1,current_column),lda,tau(current_column), &
                                         t(current_column,current_column),ldt,v(1,current_column),ldv, &
                                        baseidx,work,lwork,m,idx,mb,PQRPARAM(:),1,mpicomm,actualrank)
        else
          call qr_pdlarfgk_1dcomm_&
&single &
(obj,a(1,current_column),lda,tau(current_column), &
                                     t(current_column,current_column),ldt,v(1,current_column),ldv, &
                                     baseidx,work,lwork,m,rank,idx,mb,PQRPARAM(:),1,mpicomm,actualrank)
        end if

        blockheuristic(actualrank) = blockheuristic(actualrank) + 1

        ! the blocked decomposition versions already updated their non
        ! decomposed parts using their information after communication
        update_cols = decomposition_cols - rank
        decomposition_cols = decomposition_cols - actualrank

        ! needed for incremental update
        nextrank = min(n-(lastcol-decomposition_cols+1)+1,maxrank,rowidx-(lastcol-decomposition_cols+1)+1)

        if (current_column .gt. 1) then
          idx = rowidx-icol+1

          if (updatemode .eq. ichar('I')) then
            ! incremental update + merging
            call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,nextrank-(rank-actualrank),n-(current_column+rank-1),actualrank, &
                                          v(1,current_column+(rank-actualrank)),ldv, &
                                          t(current_column+(rank-actualrank),current_column+(rank-actualrank)),ldt, &
                                          a(1,current_column-nextrank+(rank-actualrank)),lda,baseidx,rev,updatemode,&
                                          mpicomm,work,lwork)
          else
            ! full update + merging
            call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,update_cols,n-(current_column+rank-1),actualrank, &
                                          v(1,current_column+(rank-actualrank)),ldv, &
                                          t(current_column+(rank-actualrank),current_column+(rank-actualrank)),ldt, &
                                          a(1,1),lda,baseidx,rev,updatemode,mpicomm,work,lwork)
          end if
        else
          call qr_tmerge_pdlarfb_1dcomm_&
&single &
(m,mb,0,n-(current_column+rank-1),actualrank, &
              v(1,current_column+(rank-actualrank)), &
                                          ldv, &
                                          t(current_column+(rank-actualrank),current_column+(rank-actualrank)),ldt, &
                                           a,lda,baseidx,rev,updatemode,mpicomm,work,lwork)
        end if

      end do
      call obj%timer%stop("qr_pdgeqr2_1dcomm_&
          &single&
          &")
    end subroutine

    ! incx == 1: column major
    ! incx != 1: row major
    subroutine qr_pdlarfg_1dcomm_&
&single &
(obj,x,incx,tau,work,lwork,n,idx,nb,hgmode,rev,communicator)

      use precision
      !use elpa1_impl !check this
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter    :: gmode_ = 1,rank_ = 2, eps_ = 3

      ! input variables (local)
      integer(kind=ik)               :: incx,lwork,hgmode
      real(kind=c_float)                  :: x(*),work(*)

      ! input variables (global)
      integer(kind=ik)               :: communicator,nb,idx,n,rev

      ! output variables (global)
      real(kind=c_float)                  :: tau

      ! local scalars
      integer(kind=ik)               :: mpirank,mpiprocs,mpirank_top
      integer(kind=MPI_KIND)         :: mpirankMPI, mpiprocsMPI
      integer(kind=MPI_KIND)         :: mpierr
      integer(kind=ik)               :: sendsize,recvsize
      integer(kind=ik)               :: local_size,local_offset,baseoffset
      integer(kind=ik)               :: topidx,top,iproc
      real(kind=c_float)                  :: alpha,xnorm,dot,xf

      ! external functions
      real(kind=c_float), external        :: sdot,slapy2,snrm2
      external                       :: dscal

      ! intrinsic
!      intrinsic sign
      call obj%timer%start("qr_pdlarfg_1dcomm_&
          &single&
          &")
      if (idx .le. 1) then
        tau = 0.0d0
      call obj%timer%stop("qr_pdlarfg_1dcomm_&
          &single&
          &")
        return
       end if
      call MPI_Comm_rank(int(communicator,kind=MPI_KIND) , mpirankMPI, mpierr)
      call MPI_Comm_size(int(communicator,kind=MPI_KIND) , mpiprocsMPI, mpierr)
 
      mpirank = int(mpirankMPI, kind=c_int)
      mpiprocs = int(mpiprocsMPI, kind=c_int)

      ! calculate expected work size and store in work(1)
      if (hgmode .eq. ichar('s')) then
        ! allreduce (MPI_SUM)
        sendsize = 2
        recvsize = sendsize
      else if (hgmode .eq. ichar('x')) then
        ! alltoall
        sendsize = mpiprocs*2
        recvsize = sendsize
      else if (hgmode .eq. ichar('g')) then
        ! allgather
        sendsize = 2
        recvsize = mpiprocs*sendsize
      else
        ! no exchange at all (benchmarking)
        sendsize = 2
        recvsize = sendsize
      end if

      if (lwork .eq. -1) then
        work(1) = real(sendsize + recvsize,kind=rk4)

      call obj%timer%stop("qr_pdlarfg_1dcomm_&
          &single&
          &")
        return
      end if

      ! Processor id for global index of top element
      mpirank_top = MOD((idx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        topidx = local_index(idx,mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
      end if

      call local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                        local_size,baseoffset,local_offset)

      local_offset = local_offset * incx

      ! calculate and exchange information
      if (hgmode .eq. ichar('s')) then
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk4
        end if
        dot = sdot(local_size, &
                     x(local_offset), incx, &
                     x(local_offset), incx)
        work(1) = alpha
        work(2) = dot

        call mpi_allreduce(work(1),work(sendsize+1), &
                             int(sendsize,kind=MPI_KIND), mpi_real4, mpi_sum, &
                             int(communicator,kind=MPI_KIND), mpierr)

        alpha = work(sendsize+1)
        xnorm = sqrt(work(sendsize+2))
      else if (hgmode .eq. ichar('x')) then
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk4
        end if
        xnorm = snrm2(local_size, x(local_offset), incx)
        do iproc=0,mpiprocs-1
          work(2*iproc+1) = alpha
          work(2*iproc+2) = xnorm
        end do

        call mpi_alltoall(work(1), 2_MPI_KIND, mpi_real4, &
                            work(sendsize+1), 2_MPI_KIND, mpi_real4, &
                            int(communicator,kind=MPI_KIND), mpierr)

        ! extract alpha value
        alpha = work(sendsize+1+mpirank_top*2)

        ! copy norm parts of buffer to beginning
        do iproc=0,mpiprocs-1
          work(iproc+1) = work(sendsize+1+2*iproc+1)
        end do

        xnorm = snrm2(mpiprocs, work(1), 1)
      else if (hgmode .eq. ichar('g')) then
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk4
        end if
        xnorm = snrm2(local_size, x(local_offset), incx)
        work(1) = alpha
        work(2) = xnorm

        ! allgather

        call mpi_allgather(work(1), int(sendsize,kind=MPI_KIND), mpi_real4, &
                            work(sendsize+1), int(sendsize,kind=MPI_KIND), mpi_real4, &
                            int(communicator,kind=MPI_KIND), mpierr)

        ! extract alpha value
        alpha = work(sendsize+1+mpirank_top*2)

        ! copy norm parts of buffer to beginning
        do iproc=0,mpiprocs-1
          work(iproc+1) = work(sendsize+1+2*iproc+1)
        end do
        xnorm = snrm2(mpiprocs, work(1), 1)
      else
        ! dnrm2
        xnorm = snrm2(local_size, x(local_offset), incx)
        if (mpirank .eq. mpirank_top) then
          alpha = x(top)
        else
          alpha = 0.0_rk4
        end if

        ! no exchange at all (benchmarking)
        xnorm = 0.0_rk4
      end if

      !print *,'ref hg:', idx,xnorm,alpha
      !print *,x(1:n)

      ! calculate householder information
      if (xnorm .eq. 0.0_rk4) then
        ! H = I

        tau = 0.0_rk4
      else
        ! General case
        call hh_transform_real_&
&single &
(obj,alpha,xnorm**2,xf,tau, .false.)
        if (mpirank .eq. mpirank_top) then
          x(top) = alpha
        end if
        call sscal(local_size, xf, &
                     x(local_offset), incx)

        ! TODO: reimplement norm rescale method of
        ! original PDLARFG using mpi?

      end if

      ! useful for debugging
      !print *,'hg:mpirank,idx,beta,alpha:',mpirank,idx,beta,alpha,1.0d0/(beta+alpha),tau
      !print *,x(1:n)
      call obj%timer%stop("qr_pdlarfg_1dcomm_&
          &single&
          &")
    end subroutine

    subroutine qr_pdlarfg2_1dcomm_ref_&
&single &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,work,lwork,m,idx,mb,PQRPARAM,rev,mpicomm,actualk)
      use precision
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter    :: gmode_ = 1,rank_ = 2,eps_ = 3, upmode1_ = 4
      ! input variables (local)
      integer(kind=ik)               :: lda,lwork,ldv,ldt
      real(kind=c_float)                  :: a(lda,*),v(ldv,*),tau(*),work(*),t(ldt,*)

      ! input variables (global)
      integer(kind=ik)               :: m,idx,baseidx,mb,rev,mpicomm
      integer(kind=ik)               :: PQRPARAM(:)
      ! output variables (global)
      integer(kind=ik)               :: actualk

      ! derived input variables from QR_PQRPARAM
      integer(kind=ik)               :: eps

      ! local scalars
      real(kind=c_float)                  :: dseedwork_size(1)
      integer(kind=ik)               :: seedwork_size,seed_size
      integer(kind=ik)               :: seedwork_offset,seed_offset
      logical                        :: accurate
      call obj%timer%start("qr_pdlarfg2_1dcomm_&
          &single&
          &")

      call qr_pdlarfg2_1dcomm_seed_&
&single &
(obj,a,lda,dseedwork_size(1),-1,work,m,mb,idx,rev,mpicomm)
      seedwork_size = int(dseedwork_size(1))
      seed_size = seedwork_size

      if (lwork .eq. -1) then
        work(1) = seedwork_size + seed_size
      call obj%timer%stop("qr_pdlarfg2_1dcomm_&
          &single&
          &")

        return
      end if

      seedwork_offset = 1
      seed_offset = seedwork_offset + seedwork_size

      eps = PQRPARAM(3)

      ! check for border cases (only a 2x2 matrix left)
      if (idx .le. 1) then
        tau(1:2) = 0.0_rk4
         t(1:2,1:2) = 0.0_rk4

      call obj%timer%stop("qr_pdlarfg2_1dcomm_&
          &single&
          &")

         return
      end if

      call qr_pdlarfg2_1dcomm_seed_&
&single &
(obj,a,lda,work(seedwork_offset),lwork,work(seed_offset),m,mb,idx,rev,mpicomm)

      if (eps .gt. 0) then
        accurate = qr_pdlarfg2_1dcomm_check_&
&single &
(obj,work(seed_offset),eps)
      else
        accurate = .true.
      end if

      call qr_pdlarfg2_1dcomm_vector_&
&single &
(obj,a(1,2),1,tau(2),work(seed_offset), &
                                          m,mb,idx,0,1,mpicomm)

      call qr_pdlarfg_copy_1dcomm_&
&single &
(obj,a(1,2),1, &
                                       v(1,2),1, &
                                       m,baseidx,idx,mb,1,mpicomm)

      call qr_pdlarfg2_1dcomm_update_&
&single &
(obj,v(1,2),1,baseidx,a(1,1),lda,work(seed_offset),m,idx,mb,rev,mpicomm)

      ! check for 2x2 matrix case => only one householder Vector will be
      ! generated
      if (idx .gt. 2) then
        if (accurate .eqv. .true.) then
          call qr_pdlarfg2_1dcomm_vector_&
&single &
(obj,a(1,1),1,tau(1),work(seed_offset), &
                                                  m,mb,idx-1,1,1,mpicomm)

          call qr_pdlarfg_copy_1dcomm_&
&single &
(obj,a(1,1),1, &
                                               v(1,1),1, &
                                               m,baseidx,idx-1,mb,1,mpicomm)

          ! generate fuse element
          call qr_pdlarfg2_1dcomm_finalize_tmatrix_&
&single &
(obj,work(seed_offset),tau,t,ldt)

          actualk = 2
        else
          t(1,1) = 0.0_rk4
          t(1,2) = 0.0_rk4
          t(2,2) = tau(2)

          actualk = 1
        end if
      else
        t(1,1) = 0.0_rk4
        t(1,2) = 0.0_rk4
        t(2,2) = tau(2)

        ! no more vectors to create
        tau(1) = 0.0_rk4

        actualk = 2

        !print *,'rank2: no more data'
      end if
      call obj%timer%stop("qr_pdlarfg2_1dcomm_&
          &single&
          &")

    end subroutine

    subroutine qr_pdlarfg2_1dcomm_seed_&
&single &
(obj,a,lda,work,lwork,seed,n,nb,idx,rev,mpicomm)
      use precision
      !use elpa1_impl ! check this
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)        :: lda,lwork
      real(kind=c_float)           :: a(lda,*),work(*),seed(*)

      ! input variables (global)
      integer(kind=ik)        :: n,nb,idx,rev,mpicomm

      ! output variables (global)

      ! external functions
      real(kind=c_float), external :: sdot
      ! local scalars
      real(kind=c_float)           :: top11,top21,top12,top22
      real(kind=c_float)           :: dot11,dot12,dot22
      integer(kind=ik)        :: mpirank, mpiprocs
      integer(kind=MPI_KIND)  :: mpirankMPI, mpiprocsMPI, mpierr
      integer(kind=ik)        :: mpirank_top11,mpirank_top21
      integer(kind=ik)        :: top11_offset,top21_offset
      integer(kind=ik)        :: baseoffset
      integer(kind=ik)        :: local_offset1,local_size1
      integer(kind=ik)        :: local_offset2,local_size2

      call obj%timer%start("qr_pdlarfg2_1dcomm_seed_&
          &single&
          &")

      if (lwork .eq. -1) then
        work(1) = 8.0_rk4

      call obj%timer%stop("qr_pdlarfg2_1dcomm_seed_&
          &single&
          &")
        return
      end if
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI, kind=c_int)
      mpiprocs = int(mpiprocsMPI, kind=c_int)

      call local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                                local_size1,baseoffset,local_offset1)

      call local_size_offset_1d(n,nb,idx,idx-2,rev,mpirank,mpiprocs, &
                                local_size2,baseoffset,local_offset2)

      mpirank_top11 = MOD((idx-1)/nb,mpiprocs)
      mpirank_top21 = MOD((idx-2)/nb,mpiprocs)

      top11_offset = local_index(idx,mpirank_top11,mpiprocs,nb,0)
      top21_offset = local_index(idx-1,mpirank_top21,mpiprocs,nb,0)

      if (mpirank_top11 .eq. mpirank) then
        top11 = a(top11_offset,2)
        top12 = a(top11_offset,1)
      else
        top11 = 0.0_rk4
        top12 = 0.0_rk4
      end if

      if (mpirank_top21 .eq. mpirank) then
        top21 = a(top21_offset,2)
        top22 = a(top21_offset,1)
      else
        top21 = 0.0_rk4
        top22 = 0.0_rk4
      end if

      ! calculate 3 dot products
      dot11 = sdot(local_size1,a(local_offset1,2),1,a(local_offset1,2),1)
      dot12 = sdot(local_size1,a(local_offset1,2),1,a(local_offset1,1),1)
      dot22 = sdot(local_size2,a(local_offset2,1),1,a(local_offset2,1),1)
      ! store results in work buffer
      work(1) = top11
      work(2) = dot11
      work(3) = top12
      work(4) = dot12
      work(5) = top21
      work(6) = top22
      work(7) = dot22
      work(8) = 0.0_rk4! fill up buffer
      ! exchange partial results

      call mpi_allreduce(work, seed, 8_MPI_KIND, mpi_real4, mpi_sum, &
                         int(mpicomm,kind=MPI_KIND), mpierr)


      call obj%timer%stop("qr_pdlarfg2_1dcomm_seed_&
          &single&
          &")
    end subroutine

    logical function qr_pdlarfg2_1dcomm_check_&
&single &
(obj,seed,eps)
      use precision
      use elpa_abstract_impl
      implicit none

      class(elpa_abstract_impl_t), intent(inout) :: obj

      ! input variables
      real(kind=c_float)    ::  seed(*)
      integer(kind=ik) :: eps

      ! local scalars
      real(kind=c_float)    :: epsd,first,second,first_second,estimate
      logical          :: accurate
      real(kind=c_float)    :: dot11,dot12,dot22
      real(kind=c_float)    :: top11,top12,top21,top22
      call obj%timer%start("qr_pdlarfg2_1dcomm_check_&
          &single&
          &")

      EPSD = EPS

      top11 = seed(1)
      dot11 = seed(2)
      top12 = seed(3)
      dot12 = seed(4)

      top21 = seed(5)
      top22 = seed(6)
      dot22 = seed(7)

      ! reconstruct the whole inner products
      ! (including squares of the top elements)
      first = dot11 + top11*top11
      second = dot22 + top22*top22 + top12*top12
      first_second = dot12 + top11*top12

      ! zero Householder Vector (zero norm) case
      if (first*second .eq. 0.0_rk4) then
        qr_pdlarfg2_1dcomm_check_&
&single &
 = .false.
      call obj%timer%stop("qr_pdlarfg2_1dcomm_check_&
          &single&
          &")

        return
      end if

      estimate = abs((first_second*first_second)/(first*second))

      !print *,'estimate:',estimate

      ! if accurate the following check holds
      accurate = (estimate .LE. (epsd/(1.0_rk4+epsd)))
      qr_pdlarfg2_1dcomm_check_&
&single &
 = accurate
      call obj%timer%stop("qr_pdlarfg2_1dcomm_check_&
          &single&
          &")

    end function

    ! id=0: first Vector
    ! id=1: second Vector
    subroutine qr_pdlarfg2_1dcomm_vector_&
&single &
(obj,x,incx,tau,seed,n,nb,idx,id,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)        :: incx
      real(kind=c_float)           :: x(*),seed(*),tau

      ! input variables (global)
      integer(kind=ik)        :: n,nb,idx,id,rev,mpicomm

      ! output variables (global)

      ! external functions
      real(kind=rk4), external :: slapy2
      external                :: sscal
      ! local scalars
      integer(kind=ik)        :: mpirank,mpirank_top,mpiprocs
      integer(kind=MPI_KIND)  :: mpierr, mpirankMPI, mpiprocsMPI
      real(kind=c_float)           :: alpha,dot,beta,xnorm
      integer(kind=ik)        :: local_size,baseoffset,local_offset,top,topidx
      call obj%timer%start("qr_pdlarfg2_1dcomm_vector_&
          &single&
          &")

      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      call local_size_offset_1d(n,nb,idx,idx-1,rev,mpirank,mpiprocs, &
                                    local_size,baseoffset,local_offset)

      local_offset = local_offset * incx

      ! Processor id for global index of top element
      mpirank_top = MOD((idx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        topidx = local_index(idx,mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
      else
        top = -99
        stop
      end if

      alpha = seed(id*5+1)
      dot = seed(id*5+2)

      xnorm = sqrt(dot)
      if (xnorm .eq. 0.0_rk4) then
        ! H = I

        tau = 0.0_rk4
      else
        ! General case
        beta = sign(slapy2(alpha, xnorm), alpha)
        tau = (beta+alpha) / beta

        !print *,'hg2',tau,xnorm,alpha
        call sscal(local_size, 1.0_rk4/(beta+alpha), &
                   x(local_offset), incx)

        ! TODO: reimplement norm rescale method of
        ! original PDLARFG using mpi?

        if (mpirank .eq. mpirank_top) then
          x(top) = -beta
        end if

        seed(8) = beta
      end if
      call obj%timer%stop("qr_pdlarfg2_1dcomm_vector_&
          &single&
          &")

    end subroutine

    subroutine qr_pdlarfg2_1dcomm_update_&
&single &
(obj,v,incv,baseidx,a,lda,seed,n,idx,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)   :: incv,lda
      real(kind=c_float)      :: v(*),a(lda,*),seed(*)

      ! input variables (global)
      integer(kind=ik)   :: n,baseidx,idx,nb,rev,mpicomm

      ! output variables (global)

      ! external functions
      external daxpy

      ! local scalars
      integer(kind=ik)   :: mpirank,mpiprocs
      integer(kind=MPI_KIND)   :: mpirankMPI, mpiprocsMPI, mpierr
      integer(kind=ik)   :: local_size,local_offset,baseoffset
      real(kind=c_float)      :: z,coeff,beta
      real(kind=c_float)      :: dot11,dot12,dot22
      real(kind=c_float)      :: top11,top12,top21,top22
      call obj%timer%start("qr_pdlarfg2_1dcomm_update_&
          &single&
          &")

      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      ! seed should be updated by previous householder generation
      ! Update inner product of this column and next column Vector
      top11 = seed(1)
      dot11 = seed(2)
      top12 = seed(3)
      dot12 = seed(4)

      top21 = seed(5)
      top22 = seed(6)
      dot22 = seed(7)
      beta = seed(8)

      call local_size_offset_1d(n,nb,baseidx,idx,rev,mpirank,mpiprocs, &
                                local_size,baseoffset,local_offset)
      baseoffset = baseoffset * incv

      ! zero Householder Vector (zero norm) case
      if (beta .eq. 0.0_rk4) then

      call obj%timer%stop("qr_pdlarfg2_1dcomm_update_&
          &single&
          &")
        return
      end if
      z = (dot12 + top11 * top12) / beta + top12

      !print *,'hg2 update:',baseidx,idx,mpirank,local_size
      call saxpy(local_size, -z, v(baseoffset),1, a(local_offset,1),1)
      ! prepare a full dot22 for update
      dot22 = dot22 + top22*top22

      ! calculate coefficient
      COEFF = z / (top11 + beta)

      ! update inner product of next Vector
      dot22 = dot22 - coeff * (2*dot12 - coeff*dot11)

      ! update dot12 value to represent update with first Vector
      ! (needed for T matrix)
      dot12 = dot12 - COEFF * dot11

      ! update top element of next Vector
      top22 = top22 - coeff * top21
      seed(6) = top22

      ! restore separated dot22 for Vector generation
      seed(7) = dot22  - top22*top22

      !------------------------------------------------------
      ! prepare elements for T matrix
      seed(4) = dot12

      ! prepare dot matrix for fuse element of T matrix
      ! replace top11 value with -beta1
      seed(1) = beta
      call obj%timer%stop("qr_pdlarfg2_1dcomm_update_&
          &single&
          &")

    end subroutine

    ! run this function after second Vector
    subroutine qr_pdlarfg2_1dcomm_finalize_tmatrix_&
&single &
(obj,seed,tau,t,ldt)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj

      integer(kind=ik)  :: ldt
      real(kind=c_float)     :: seed(*),t(ldt,*),tau(*)
      real(kind=c_float)     :: dot12,beta1,top21,beta2
      call obj%timer%start("qr_pdlarfg2_1dcomm_finalize_tmatrix_&
          &single&
          &")

      beta1 = seed(1)
      dot12 = seed(4)
      top21 = seed(5)
      beta2 = seed(8)

      !print *,'beta1 beta2',beta1,beta2

      dot12 = dot12 / beta2 + top21
      dot12 = -(dot12 / beta1)

      t(1,1) = tau(1)
      t(1,2) = dot12
      t(2,2) = tau(2)
      call obj%timer%stop("qr_pdlarfg2_1dcomm_finalize_tmatrix_&
          &single&
          &")

    end subroutine

    subroutine qr_pdlarfgk_1dcomm_&
&single &
(obj,a,lda,tau,t,ldt,v,ldv,baseidx,work,lwork,m,k,idx,mb,PQRPARAM,rev,mpicomm,actualk)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup

      ! input variables (local)
      integer(kind=ik)    :: lda,lwork,ldv,ldt
      real(kind=c_float)       :: a(lda,*),v(ldv,*),tau(*),work(*),t(ldt,*)

      ! input variables (global)
      integer(kind=ik)    :: m,k,idx,baseidx,mb,rev,mpicomm
      integer(kind=ik)    :: PQRPARAM(:)
      ! output variables (global)
      integer(kind=ik)    :: actualk

      ! local scalars
      integer(kind=ik)    :: ivector
      real(kind=c_float)       :: pdlarfg_size(1),pdlarf_size(1)
      real(kind=c_float)       :: pdlarfgk_1dcomm_seed_size(1),pdlarfgk_1dcomm_check_size(1)
      real(kind=c_float)       :: pdlarfgk_1dcomm_update_size(1)
      integer(kind=ik)    :: seedC_size,seedC_offset
      integer(kind=ik)    :: seedD_size,seedD_offset
      integer(kind=ik)    :: work_offset
      call obj%timer%start("qr_pdlarfgk_1dcomm_&
          &single&
          &")

      seedC_size = k*k
      seedC_offset = 1
      seedD_size = k*k
      seedD_offset = seedC_offset + seedC_size
      work_offset = seedD_offset + seedD_size

      if (lwork .eq. -1) then
        call qr_pdlarfg_1dcomm_&
&single &
(obj, a,1,tau(1),pdlarfg_size(1),-1,m,baseidx,mb,PQRPARAM(4),rev,mpicomm)

        call qr_pdlarfl_1dcomm_&
&single &
(v,1,baseidx,a,lda,tau(1),pdlarf_size(1),-1,m,k,baseidx,mb,rev,mpicomm)
        call qr_pdlarfgk_1dcomm_seed_&
&single &
(obj,a,lda,baseidx,pdlarfgk_1dcomm_seed_size(1),-1,work,work,m,k,mb,mpicomm)
        !call qr_pdlarfgk_1dcomm_check_&
!&single &
!(work,work,k,PQRPARAM(:),pdlarfgk_1dcomm_check_size(1),-1,actualk)
        call qr_pdlarfgk_1dcomm_check_improved_&
&single &
(obj,work,work,k,PQRPARAM(:),pdlarfgk_1dcomm_check_size(1),-1,actualk)
        call qr_pdlarfgk_1dcomm_update_&
&single &
(obj,a,lda,baseidx,pdlarfgk_1dcomm_update_size(1), &
                                              -1,work,work,k,k,1,work,m,mb,rev,mpicomm)
        work(1) = max(pdlarfg_size(1),pdlarf_size(1),pdlarfgk_1dcomm_seed_size(1),pdlarfgk_1dcomm_check_size(1), &
                        pdlarfgk_1dcomm_update_size(1)) + real(seedC_size + seedD_size, kind=c_float)

      call obj%timer%stop("qr_pdlarfgk_1dcomm_&
          &single&
          &")

        return
      end if

      call qr_pdlarfgk_1dcomm_seed_&
&single &
(obj,a(1,1),lda,idx,work(work_offset),lwork,work(seedC_offset), &
          work(seedD_offset),m,k,mb,mpicomm)
      !call qr_pdlarfgk_1dcomm_check_&
!&single &
!(work(seedC_offset),work(seedD_offset),k,PQRPARAM(:),work(work_offset),lwork,actualk)
      call qr_pdlarfgk_1dcomm_check_improved_&
&single &
(obj,work(seedC_offset),work(seedD_offset), &
          k,PQRPARAM(:),work(work_offset),lwork,actualk)
      !print *,'possible rank:', actualk

      ! override useful for debugging
      !actualk = 1
      !actualk = k
      !actualk= min(actualk,2)
      do ivector=1,actualk
        call qr_pdlarfgk_1dcomm_vector_&
&single &
(obj,a(1,k-ivector+1),1,idx,tau(k-ivector+1), &
                                          work(seedC_offset),work(seedD_offset),k, &
                                          ivector,m,mb,rev,mpicomm)

        call qr_pdlarfgk_1dcomm_update_&
&single &
(obj,a(1,1),lda,idx,work(work_offset),lwork,work(seedC_offset), &
                                          work(seedD_offset),k,actualk,ivector,tau, &
                                          m,mb,rev,mpicomm)

        call qr_pdlarfg_copy_1dcomm_&
&single &
(obj,a(1,k-ivector+1),1, &
                                       v(1,k-ivector+1),1, &
                                       m,baseidx,idx-ivector+1,mb,1,mpicomm)
      end do

      ! generate final T matrix and convert preliminary tau values into real ones
      call qr_pdlarfgk_1dcomm_generateT_&
&single &
(obj,work(seedC_offset),work(seedD_offset),k,actualk,tau,t,ldt)

      call obj%timer%stop("qr_pdlarfgk_1dcomm_&
          &single&
          &")
    end subroutine

    subroutine qr_pdlarfgk_1dcomm_seed_&
&single &
(obj,a,lda,baseidx,work,lwork,seedC,seedD,m,k,mb,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup

      ! input variables (local)
      integer(kind=ik)   :: lda,lwork
      real(kind=c_float)      :: a(lda,*), work(*)

      ! input variables (global)
      integer(kind=ik)   :: m,k,baseidx,mb,mpicomm
      real(kind=c_float)      :: seedC(k,*),seedD(k,*)

      ! output variables (global)

      ! derived input variables from QR_PQRPARAM

      ! local scalars
      integer(kind=ik)   :: mpirank,mpiprocs,mpirank_top
      integer(kind=MPI_KIND)   :: mpierr,mpirankMPI,mpiprocsMPI
      integer(kind=ik)   :: icol,irow,lidx,remsize
      integer(kind=ik)   :: remaining_rank

      integer(kind=ik)   :: C_size,D_size,sendoffset,recvoffset,sendrecv_size
      integer(kind=ik)   :: localoffset,localsize,baseoffset
      call obj%timer%start("qr_pdlarfgk_1dcomm_seed_&
          &single&
          &")

      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      C_size = k*k
      D_size = k*k
      sendoffset = 1
      sendrecv_size = C_size+D_size
      recvoffset = sendoffset + sendrecv_size

      if (lwork .eq. -1) then
        work(1) = real(2*sendrecv_size,kind=rk4)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_seed_&
          &single&
          &")

        return
      end if

      ! clear buffer
      work(sendoffset:sendoffset+sendrecv_size-1)=0.0_rk4
      ! collect C part
      do icol=1,k

        remaining_rank = k
        do while (remaining_rank .gt. 0)
          irow = k - remaining_rank + 1
          lidx = baseidx - remaining_rank + 1

          ! determine chunk where the current top element is located
          mpirank_top = MOD((lidx-1)/mb,mpiprocs)

          ! limit max number of remaining elements of this chunk to the block
          ! distribution parameter
          remsize = min(remaining_rank,mb)

          ! determine the number of needed elements in this chunk
          call local_size_offset_1d(lidx+remsize-1,mb, &
                                    lidx,lidx,0, &
                                    mpirank_top,mpiprocs, &
                                    localsize,baseoffset,localoffset)

          !print *,'local rank',localsize,localoffset

          if (mpirank .eq. mpirank_top) then
            ! copy elements to buffer
            work(sendoffset+(icol-1)*k+irow-1:sendoffset+(icol-1)*k+irow-1+localsize-1) &
                          = a(localoffset:localoffset+remsize-1,icol)
          end if

          ! jump to next chunk
          remaining_rank = remaining_rank - localsize
        end do
      end do

      ! collect D part
      call local_size_offset_1d(m,mb,baseidx-k,baseidx-k,1, &
                        mpirank,mpiprocs, &
                        localsize,baseoffset,localoffset)

      !print *,'localsize',localsize,localoffset
      if (localsize > 0) then
        call ssyrk("Upper", "Trans", k, localsize, &
                     1.0_rk4, a(localoffset,1), lda, &
                     0.0_rk4, work(sendoffset+C_size), k)
      else
        work(sendoffset+C_size:sendoffset+C_size+k*k-1) = 0.0_rk4
      end if

      ! TODO: store symmetric part more efficiently

      ! allreduce operation on results

      call mpi_allreduce(work(sendoffset),work(recvoffset), int(sendrecv_size,kind=MPI_KIND), &
                         mpi_real4, mpi_sum, int(mpicomm,kind=MPI_KIND), mpierr)

      ! unpack result from buffer into seedC and seedD
      seedC(1:k,1:k) = 0.0_rk4
      do icol=1,k
        seedC(1:k,icol) = work(recvoffset+(icol-1)*k:recvoffset+icol*k-1)
      end do
      seedD(1:k,1:k) = 0.0_rk4
      do icol=1,k
        seedD(1:k,icol) = work(recvoffset+C_size+(icol-1)*k:recvoffset+C_size+icol*k-1)
      end do

      call obj%timer%stop("qr_pdlarfgk_1dcomm_seed_&
          &single&
          &")

    end subroutine

    ! k is assumed to be larger than two
    subroutine qr_pdlarfgk_1dcomm_check_improved_&
&single &
(obj,seedC,seedD,k,PQRPARAM,work,lwork,possiblerank)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (global)
      integer(kind=ik)   :: k,lwork
      integer(kind=ik)   :: PQRPARAM(:)
      real(kind=c_float)      :: seedC(k,*),seedD(k,*),work(k,*)

      ! output variables (global)
      integer(kind=ik)   :: possiblerank

      ! derived input variables from QR_PQRPARAM
      integer(kind=ik)   :: eps

      ! local variables
      integer(kind=ik)   :: i,j,l
      real(kind=c_float)      :: sum_squares,diagonal_square,epsd,diagonal_root
      real(kind=c_float)      :: dreverse_matrix_work(1)

      ! external functions
      real(kind=rk4), external :: sdot,slapy2,snrm2
      external                :: sscal

      call obj%timer%start("qr_pdlarfgk_1dcomm_check_improved_&
          &single&
          &")

      if (lwork .eq. -1) then
        call reverse_matrix_local_&
&single &
(1,k,k,work,k,dreverse_matrix_work,-1)
        work(1,1) = real(k*k,kind=rk4) + dreverse_matrix_work(1)

        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &single&
            &")
        return
      end if

      eps = PQRPARAM(3)

      if (eps .eq. 0) then
        possiblerank = k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &single&
            &")
        return
      end if
      epsd = real(eps,kind=rk4)

      ! build complete inner product from seedC and seedD
      ! copy seedD to work
      work(:,1:k) = seedD(:,1:k)

      ! add inner products of seedC to work
      call ssyrk("Upper", "Trans", k, k, &
                 1.0_rk4, seedC(1,1), k, &
                 1.0_rk4, work, k)


      ! TODO: optimize this part!
      call reverse_matrix_local_&
&single &
(0,k,k,work(1,1),k,work(1,k+1),lwork-2*k)
      call reverse_matrix_local_&
&single &
(1,k,k,work(1,1),k,work(1,k+1),lwork-2*k)

      ! transpose matrix
      do i=1,k
        do j=i+1,k
          work(i,j) = work(j,i)
        end do
      end do


      ! do cholesky decomposition
      i = 0
      do while ((i .lt. k))
        i = i + 1

        diagonal_square = abs(work(i,i))
        diagonal_root  = sqrt(diagonal_square)

        ! zero Householder Vector (zero norm) case
        if ((abs(diagonal_square) .eq. 0.0_rk4) .or. (abs(diagonal_root) .eq. 0.0_rk4)) then
          possiblerank = max(i-1,1)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &single&
            &")
          return
        end if

        ! check if relative error is bounded for each Householder Vector
        ! Householder i is stable iff Househoulder i-1 is "stable" and the accuracy criterion
        ! holds.
        ! first Householder Vector is considered as "stable".

        do j=i+1,k
          work(i,j) = work(i,j) / diagonal_root
          do l=i+1,j
            work(l,j) = work(l,j) - work(i,j) * work(i,l)
          end do
        end do
        !print *,'cholesky step done'

        ! build sum of squares
        if (i .eq. 1) then
          sum_squares = 0.0_rk4
        else
          sum_squares = sdot(i-1,work(1,i),1,work(1,i),1)
        end if
        if (sum_squares .ge. (epsd * diagonal_square)) then
          possiblerank = max(i-1,1)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &single&
            &")
          return
        end if
      end do

      possiblerank = i
      !print *,'possible rank', possiblerank
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_improved_&
            &single&
            &")

    end subroutine

    ! TODO: zero Householder Vector (zero norm) case
    ! - check alpha values as well (from seedC)
    subroutine qr_pdlarfgk_1dcomm_check_&
&single &
(obj,seedC,seedD,k,PQRPARAM,work,lwork,possiblerank)
      use precision
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup

      ! input variables (local)

      ! input variables (global)
      integer(kind=ik)   :: k,lwork
      integer(kind=ik)   :: PQRPARAM(:)
      real(kind=c_float)      :: seedC(k,*),seedD(k,*),work(k,*)

      ! output variables (global)
      integer(kind=ik)   :: possiblerank

      ! derived input variables from QR_PQRPARAM
      integer(kind=ik)   :: eps

      ! local scalars
      integer(kind=ik)   :: icol,isqr,iprod
      real(kind=c_float)      :: epsd,sum_sqr,sum_products,diff,temp,ortho,ortho_sum
      real(kind=c_float)      :: dreverse_matrix_work(1)
        call obj%timer%start("qr_pdlarfgk_1dcomm_check_&
            &single&
            &")
      if (lwork .eq. -1) then
        call reverse_matrix_local_&
&single &
(1,k,k,work,k,dreverse_matrix_work,-1)
        work(1,1) = real(k*k,kind=rk4) + dreverse_matrix_work(1)

        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &single&
            &")

        return
      end if

      eps = PQRPARAM(3)

      if (eps .eq. 0) then
        possiblerank = k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &single&
            &")
        return
      end if
      epsd = real(eps,kind=rk4)

      ! copy seedD to work
      work(:,1:k) = seedD(:,1:k)

      ! add inner products of seedC to work
      call ssyrk("Upper", "Trans", k, k, &
                 1.0_rk4, seedC(1,1), k, &
                 1.0_rk4, work, k)

      ! TODO: optimize this part!
      call reverse_matrix_local_&
&single &
(0,k,k,work(1,1),k,work(1,k+1),lwork-2*k)
      call reverse_matrix_local_&
&single &
(1,k,k,work(1,1),k,work(1,k+1),lwork-2*k)

      ! transpose matrix
      do icol=1,k
        do isqr=icol+1,k
          work(icol,isqr) = work(isqr,icol)
        end do
      end do

      ! work contains now the full inner product of the global (sub-)matrix
      do icol=1,k
        ! zero Householder Vector (zero norm) case
        if (abs(work(icol,icol)) .eq. 0.0_rk4) then
          !print *,'too small ', icol, work(icol,icol)
          possiblerank = max(icol,1)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &single&
            &")
          return
        end if

        sum_sqr = 0.0_rk4
        do isqr=1,icol-1
          sum_products = 0.0_rk4
          do iprod=1,isqr-1
            sum_products = sum_products + work(iprod,isqr)*work(iprod,icol)
          end do

          !print *,'divisor',icol,isqr,work(isqr,isqr)
          temp = (work(isqr,icol) - sum_products)/work(isqr,isqr)
          work(isqr,icol) = temp
          sum_sqr = sum_sqr + temp*temp
        end do

        ! calculate diagonal value
        diff = work(icol,icol) - sum_sqr
        if (diff .lt. 0.0_rk4) then
          ! we definitely have a problem now
          possiblerank = icol-1 ! only decompose to previous column (including)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &single&
            &")
          return
        end if
        work(icol,icol) = sqrt(diff)
        ! calculate orthogonality
        ortho = 0.0_rk4
        do isqr=1,icol-1
          ortho_sum = 0.0_rk4
          do iprod=isqr,icol-1
            temp = work(isqr,iprod)*work(isqr,iprod)
            !print *,'ortho ', work(iprod,iprod)
            temp = temp / (work(iprod,iprod)*work(iprod,iprod))
            ortho_sum = ortho_sum + temp
          end do
          ortho = ortho + ortho_sum * (work(isqr,icol)*work(isqr,icol))
        end do

        ! ---------------- with division by zero ----------------------- !

        !ortho = ortho / diff;

        ! if current estimate is not accurate enough, the following check holds
        !if (ortho .gt. epsd) then
        !    possiblerank = icol-1 ! only decompose to previous column (including)
        !    return
        !end if

        ! ---------------- without division by zero ----------------------- !

        ! if current estimate is not accurate enough, the following check holds
        if (ortho .gt. epsd * diff) then
          possiblerank = icol-1 ! only decompose to previous column (including)
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &single&
            &")
          return
        end if
      end do

      ! if we get to this point, the accuracy condition holds for the whole block
      possiblerank = k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_check_&
            &single&
            &")
      end subroutine

    !sidx: seed idx
    !k: max rank used during seed phase
    !rank: actual rank (k >= rank)
    subroutine qr_pdlarfgk_1dcomm_vector_&
&single &
(obj,x,incx,baseidx,tau,seedC,seedD,k,sidx,n,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)  :: incx
      real(kind=c_float)     :: x(*),tau

      ! input variables (global)
      integer(kind=ik)  :: n,nb,baseidx,rev,mpicomm,k,sidx
      real(kind=c_float)     :: seedC(k,*),seedD(k,*)

      ! output variables (global)

      ! external functions
      real(kind=rk4), external :: slapy2,snrm2
      external                :: sscal

      ! local scalars
      integer(kind=ik)   :: mpirank,mpirank_top,mpiprocs
      integer(kind=MPI_KIND) :: mpirankMPI, mpiprocsMPI, mpierr
      real(kind=c_float)      :: alpha,dot,beta,xnorm
      integer(kind=ik)   :: local_size,baseoffset,local_offset,top,topidx
      integer(kind=ik)   :: lidx
        call obj%timer%start("qr_pdlarfgk_1dcomm_vector_&
            &single&
            &")
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)

      lidx = baseidx-sidx+1
      call local_size_offset_1d(n,nb,baseidx,lidx-1,rev,mpirank,mpiprocs, &
                        local_size,baseoffset,local_offset)

      local_offset = local_offset * incx

      ! Processor id for global index of top element
      mpirank_top = MOD((lidx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        topidx = local_index((lidx),mpirank_top,mpiprocs,nb,0)
        top = 1+(topidx-1)*incx
      end if

      alpha = seedC(k-sidx+1,k-sidx+1)
      dot = seedD(k-sidx+1,k-sidx+1)
      ! assemble actual norm from both seed parts
      xnorm = slapy2(sqrt(dot), snrm2(k-sidx,seedC(1,k-sidx+1),1))

      if (xnorm .eq. 0.0_rk4) then
        tau = 0.0_rk4
      else

        ! General case

        beta = sign(slapy2(alpha, xnorm), alpha)
        ! store a preliminary version of beta in tau
        tau = beta

        ! update global part
        call sscal(local_size, 1.0_rk4/(beta+alpha), &
                     x(local_offset), incx)

        ! do not update local part here due to
        ! dependency of c Vector during update process

        ! TODO: reimplement norm rescale method of
        ! original PDLARFG using mpi?

        if (mpirank .eq. mpirank_top) then
          x(top) = -beta
        end if
      end if
        call obj%timer%stop("qr_pdlarfgk_1dcomm_vector_&
            &single&
            &")

    end subroutine

    !k: original max rank used during seed function
    !rank: possible rank as from check function
    ! TODO: if rank is less than k, reduce buffersize in such a way
    ! that only the required entries for the next pdlarfg steps are
    ! computed
    subroutine qr_pdlarfgk_1dcomm_update_&
&single &
(obj,a,lda,baseidx,work,lwork,seedC,seedD,k,rank,sidx,tau,n,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! parameter setup
      INTEGER(kind=ik), parameter :: gmode_ = 1,rank_ = 2,eps_ = 3, upmode1_ = 4

      ! input variables (local)
      integer(kind=ik)            :: lda,lwork
      real(kind=c_float)               :: a(lda,*),work(*)

      ! input variables (global)
      integer(kind=ik)            :: k,rank,sidx,n,baseidx,nb,rev,mpicomm
      real(kind=c_float)               :: beta

      ! output variables (global)
      real(kind=c_float)               :: seedC(k,*),seedD(k,*),tau(*)

      ! derived input variables from QR_PQRPARAM

      ! local scalars
      real(kind=c_float)               :: alpha
      integer(kind=ik)            :: coffset,zoffset,yoffset,voffset,buffersize
      integer(kind=ik)            :: mpirank,mpiprocs,mpirank_top
      integer(kind=MPI_KIND)      :: mpirankMPI, mpierr,mpiprocsMPI 

      integer(kind=ik)            :: localsize,baseoffset,localoffset,topidx
      integer(kind=ik)            :: lidx
        call obj%timer%start("qr_pdlarfgk_1dcomm_update_&
            &single&
            &")
      if (lwork .eq. -1) then
        ! buffer for c,z,y,v
        work(1) = 4*k
        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &single&
            &")

        return
      end if

      ! nothing to update anymore
      if (sidx .gt. rank) then
        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &single&
            &")
        return
      endif
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)

      lidx = baseidx-sidx
      if (lidx .lt. 1) then
        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &single&
            &")
        return
      endif

      call local_size_offset_1d(n,nb,baseidx,lidx,rev,mpirank,mpiprocs, &
                                localsize,baseoffset,localoffset)

      coffset = 1
      zoffset = coffset + k
      yoffset = zoffset + k
      voffset = yoffset + k
      buffersize = k - sidx

      ! finalize tau values
      alpha = seedC(k-sidx+1,k-sidx+1)
      beta = tau(k-sidx+1)

      ! zero Householder Vector (zero norm) case
      !print *,'k update: alpha,beta',alpha,beta
      if ((beta .eq. 0.0_rk4) .or. (alpha .eq. 0.0_rk4))  then
        tau(k-sidx+1) = 0.0_rk4
        seedC(k,k-sidx+1) = 0.0_rk4

        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &single&
            &")
        return
      end if

      tau(k-sidx+1) = (beta+alpha) / beta

      ! ---------------------------------------
      ! calculate c Vector (extra Vector or encode in seedC/seedD?
      work(coffset:coffset+buffersize-1) = seedD(1:buffersize,k-sidx+1)
      call sgemv("Trans", buffersize+1, buffersize, &
                 1.0_rk4,seedC(1,1),k,seedC(1,k-sidx+1),1, &
                 1.0_rk4,work(coffset),1)

      ! calculate z using tau,seedD,seedC and c Vector
      work(zoffset:zoffset+buffersize-1) = seedC(k-sidx+1,1:buffersize)
      call saxpy(buffersize, 1.0_rk4/beta, work(coffset), 1, work(zoffset), 1)

      ! update A1(local copy) and generate part of householder vectors for use
      call saxpy(buffersize, -1.0_rk4, work(zoffset),1,seedC(k-sidx+1,1),k)
      call sscal(buffersize, 1.0_rk4/(alpha+beta), seedC(1,k-sidx+1),1)
      call sger(buffersize, buffersize, -1.0_rk4, seedC(1,k-sidx+1),1, work(zoffset), 1, seedC(1,1), k)

      ! update A global (householder Vector already generated by pdlarfgk)
      mpirank_top = MOD(lidx/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        ! handle first row separately
        topidx = local_index(lidx+1,mpirank_top,mpiprocs,nb,0)
        call saxpy(buffersize,-1.0_rk4,work(zoffset),1,a(topidx,1),lda)
      end if

      call sger(localsize, buffersize,-1.0_rk4, &
                a(localoffset,k-sidx+1),1,work(zoffset),1, &
                a(localoffset,1),lda)

      ! update D (symmetric) => two buffer vectors of size rank
      ! generate y Vector
      work(yoffset:yoffset+buffersize-1) = 0._rk4
      call saxpy(buffersize,1.0_rk4/(alpha+beta),work(zoffset),1,work(yoffset),1)

      ! generate v Vector
      work(voffset:voffset+buffersize-1) = seedD(1:buffersize,k-sidx+1)
      call saxpy(buffersize, -0.5_rk4*seedD(k-sidx+1,k-sidx+1), work(yoffset), 1, work(voffset),1)

      ! symmetric update of D using y and v
      call ssyr2("Upper", buffersize,-1.0_rk4, &
                     work(yoffset),1,work(voffset),1, &
                     seedD(1,1), k)

      ! prepare T matrix inner products
      ! D_k(1:k,k+1:n) = D_(k-1)(1:k,k+1:n) - D_(k-1)(1:k,k) * y'
      ! store coefficient 1.0d0/(alpha+beta) in C diagonal elements
      call sger(k-sidx,sidx,-1.0_rk4,work(yoffset),1,seedD(k-sidx+1,k-sidx+1),k,seedD(1,k-sidx+1),k)
      seedC(k,k-sidx+1) = 1.0_rk4/(alpha+beta)

        call obj%timer%stop("qr_pdlarfgk_1dcomm_update_&
            &single&
            &")
    end subroutine

    subroutine qr_pdlarfgk_1dcomm_generateT_&
          &single &
          (obj,seedC,seedD,k,actualk,tau,t,ldt)
      use precision
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      integer(kind=ik)  :: k,actualk,ldt
      real(kind=c_float)     :: seedC(k,*),seedD(k,*),tau(*),t(ldt,*)

      integer(kind=ik)  :: irow,icol
      real(kind=c_float)     :: column_coefficient
        call obj%timer%start("qr_pdlarfgk_1dcomm_generateT_&
            &single&
            &")

      !print *,'reversed on the fly T generation NYI'

      do icol=1,actualk-1
        ! calculate inner product of householder Vector parts in seedC
        ! (actually calculating more than necessary, if actualk < k)
        ! => a lot of junk from row 1 to row k-actualk
        call strmv('Upper','Trans','Unit',k-icol,seedC(1,1),k,seedC(1,k-icol+1),1)
        ! add scaled D parts to current column of C (will become later T rows)
        column_coefficient = seedC(k,k-icol+1)
        do irow=k-actualk+1,k-1
          seedC(irow,k-icol+1) = ( seedC(irow,k-icol+1) ) +  ( seedD(irow,k-icol+1) * column_coefficient * seedC(k,irow) )
        end do
      end do

      call qr_dlarft_kernel_&
             &single &
             (actualk,tau(k-actualk+1),seedC(k-actualk+1,k-actualk+2),k,t(k-actualk+1,k-actualk+1),ldt)
      call obj%timer%stop("qr_pdlarfgk_1dcomm_generateT_&
             &single&
             &")

    end subroutine

    !direction=0: pack into work buffer
    !direction=1: unpack from work buffer
    subroutine qr_pdgeqrf_pack_unpack_&
&single &
(obj,v,ldv,work,lwork,m,n,mb,baseidx,rowidx,rev,direction,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)   :: ldv,lwork
      real(kind=c_float)      :: v(ldv,*), work(*)

      ! input variables (global)
      integer(kind=ik)   :: m,n,mb,baseidx,rowidx,rev,direction,mpicomm

      ! output variables (global)

      ! local scalars
      integer(kind=ik)   :: mpirank,mpiprocs
      integer(kind=MPI_KIND) :: mpierr,mpirankMPI, mpiprocsMPI

      integer(kind=ik)   :: buffersize,icol
      integer(kind=ik)   :: local_size,baseoffset,offset

      ! external functions
        call obj%timer%start("qr_pdgeqrf_pack_unpack_&
            &single&
            &")
      call mpi_comm_rank(int(mpicomm,kind=MPI_KIND) ,mpirankMPI,mpierr)
      call mpi_comm_size(int(mpicomm,kind=MPI_KIND) ,mpiprocsMPI,mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)

      call local_size_offset_1d(m,mb,baseidx,rowidx,rev,mpirank,mpiprocs, &
                                    local_size,baseoffset,offset)

      !print *,'pack/unpack',local_size,baseoffset,offset

      ! rough approximate for buffer size
      if (lwork .eq. -1) then
        buffersize = local_size * n ! Vector elements
        work(1) = DBLE(buffersize)
        call obj%timer%stop("qr_pdgeqrf_pack_unpack_&
            &single&
            &")

        return
      end if

      if (direction .eq. 0) then
        ! copy v part to buffer (including zeros)
        do icol=1,n
          work(1+local_size*(icol-1):local_size*icol) = v(baseoffset:baseoffset+local_size-1,icol)
        end do
      else
        ! copy v part from buffer (including zeros)
        do icol=1,n
          v(baseoffset:baseoffset+local_size-1,icol) = work(1+local_size*(icol-1):local_size*icol)
        end do
      end if
        call obj%timer%stop("qr_pdgeqrf_pack_unpack_&
            &single&
            &")

      return

    end subroutine

    !direction=0: pack into work buffer
    !direction=1: unpack from work buffer
    subroutine qr_pdgeqrf_pack_unpack_tmatrix_&
          &single &
          (obj,tau,t,ldt,work,lwork,n,direction)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)  :: ldt,lwork
      real(kind=c_float)     :: work(*), t(ldt,*),tau(*)

      ! input variables (global)
      integer(kind=ik)  :: n,direction

      ! output variables (global)

      ! local scalars
      integer(kind=ik)  :: icol

      ! external functions
        call obj%timer%start("qr_pdgeqrf_pack_unpack_tmatrix_&
            &single&
            &")

      if (lwork .eq. -1) then
        work(1) = real(n*n,kind=rk4)

        call obj%timer%stop("qr_pdgeqrf_pack_unpack_tmatrix_&
            &single&
            &")
        return
      end if

      if (direction .eq. 0) then
        ! append t matrix to buffer (including zeros)
        do icol=1,n
          work(1+(icol-1)*n:icol*n) = t(1:n,icol)
        end do
      else
        ! append t matrix from buffer (including zeros)
        do icol=1,n
          t(1:n,icol) = work(1+(icol-1)*n:icol*n)
          tau(icol) = t(icol,icol)
        end do
      end if
        call obj%timer%stop("qr_pdgeqrf_pack_unpack_tmatrix_&
            &single&
            &")
    end subroutine



    subroutine qr_pdlarfg_copy_1dcomm_&
&single &
(obj,x,incx,v,incv,n,baseidx,idx,nb,rev,mpicomm)
      use precision
      use elpa1_impl
      use qr_utils_mod
      use elpa_abstract_impl
      implicit none
      class(elpa_abstract_impl_t), intent(inout) :: obj
      ! input variables (local)
      integer(kind=ik)  :: incx,incv
      real(kind=c_float)     :: x(*), v(*)

      ! input variables (global)
      integer(kind=ik)  :: baseidx,idx,rev,nb,n
      integer(kind=ik)  :: mpicomm

      ! output variables (global)

      ! local scalars
      integer(kind=ik)  :: mpiprocs
      integer(kind=MPI_KIND) ::  mpierr,mpiprocsMPI,mpirankMPI

      integer(kind=ik)  :: mpirank,mpirank_top
      integer(kind=ik)  :: irow,x_offset
      integer(kind=ik)  :: v_offset,local_size

        call obj%timer%start("qr_pdlarfg_copy_1dcomm_&
            &single&
            &")
      call MPI_Comm_rank(int(mpicomm,kind=MPI_KIND), mpirankMPI, mpierr)
      call MPI_Comm_size(int(mpicomm,kind=MPI_KIND), mpiprocsMPI, mpierr)

      mpirank = int(mpirankMPI,kind=c_int)
      mpiprocs = int(mpiprocsMPI,kind=c_int)
      call local_size_offset_1d(n,nb,baseidx,idx,rev,mpirank,mpiprocs, &
                                local_size,v_offset,x_offset)
      v_offset = v_offset * incv

      !print *,'copy:',mpirank,baseidx,v_offset,x_offset,local_size

      ! copy elements
      do irow=1,local_size
        v((irow-1)*incv+v_offset) = x((irow-1)*incx+x_offset)
      end do

      ! replace top element to build an unitary Vector
      mpirank_top = MOD((idx-1)/nb,mpiprocs)
      if (mpirank .eq. mpirank_top) then
        v(local_size*incv) = 1.0_rk4
      end if
        call obj%timer%stop("qr_pdlarfg_copy_1dcomm_&
            &single&
            &")

    end subroutine

! vim: syntax=fortran


end module elpa_pdgeqrf
