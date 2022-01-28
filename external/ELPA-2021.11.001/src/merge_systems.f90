











module merge_systems
  use precision
  implicit none
  private

  public :: merge_systems_double
  public :: merge_systems_single

  contains

! real double precision first




















!cannot use "../src/solve_tridi/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


    subroutine merge_systems_&
    &double &
                         (obj, na, nm, d, e, q, ldq, nqoff, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, &
                          l_col, p_col, l_col_out, p_col_out, npc_0, npc_n, useGPU, wantDebug, success, max_threads)
      use elpa_gpu
      use, intrinsic :: iso_c_binding
      use precision
      use elpa_abstract_impl
      use elpa_blas_interfaces
      use global_product
      use global_gather
      use resort_ev
      use transform_columns
      use check_monotony
      use add_tmp
      use v_add_s
      use ELPA_utilities
      use elpa_mpi
      use solve_secular_equation
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

      class(elpa_abstract_impl_t), intent(inout)  :: obj
      integer(kind=ik), intent(in)                :: na, nm, ldq, nqoff, nblk, matrixCols, mpi_comm_rows, &
                                                     mpi_comm_cols, npc_0, npc_n
      integer(kind=ik), intent(in)                :: l_col(na), p_col(na), l_col_out(na), p_col_out(na)
      real(kind=rk8), intent(inout)     :: d(na), e
      real(kind=rk8), intent(inout)     :: q(ldq,*)
      logical, intent(in)                         :: useGPU, wantDebug
      logical                                     :: useIntelGPU

      logical, intent(out)                        :: success

      ! TODO: play with max_strip. If it was larger, matrices being multiplied
      ! might be larger as well!
      integer(kind=ik), parameter                 :: max_strip=128

      
      real(kind=rk8)                    :: beta, sig, s, c, t, tau, rho, eps, tol, &
                                                     qtrans(2,2), dmax, zmax, d1new, d2new
      real(kind=rk8)                    :: z(na), d1(na), d2(na), z1(na), delta(na),  &
                                                     dbase(na), ddiff(na), ev_scale(na), tmp(na)
      real(kind=rk8)                    :: d1u(na), zu(na), d1l(na), zl(na)
      real(kind=rk8), allocatable       :: qtmp1(:,:), qtmp2(:,:), ev(:,:)

      integer(kind=ik)                            :: i, j, na1, na2, l_rows, l_cols, l_rqs, l_rqe, &
                                                     l_rqm, ns, info
      integer(kind=BLAS_KIND)                     :: infoBLAS
      integer(kind=ik)                            :: l_rnm, nnzu, nnzl, ndef, ncnt, max_local_cols, &
                                                     l_cols_qreorg, np, l_idx, nqcols1, nqcols2
      integer(kind=ik)                            :: my_proc, n_procs, my_prow, my_pcol, np_rows, &
                                                     np_cols
      integer(kind=MPI_KIND)                      :: mpierr
      integer(kind=MPI_KIND)                      :: my_prowMPI, np_rowsMPI, my_pcolMPI, np_colsMPI
      integer(kind=ik)                            :: np_next, np_prev, np_rem
      integer(kind=ik)                            :: idx(na), idx1(na), idx2(na)
      integer(kind=BLAS_KIND)                     :: idxBLAS(NA)
      integer(kind=ik)                            :: coltyp(na), idxq1(na), idxq2(na)

      integer(kind=ik)                            :: istat
      character(200)                              :: errorMessage
      integer(kind=ik)                            :: gemm_dim_k, gemm_dim_l, gemm_dim_m

      integer(kind=c_intptr_t)                    :: num
      integer(kind=C_intptr_T)                    :: qtmp1_dev, qtmp2_dev, ev_dev
      logical                                     :: successGPU
      integer(kind=c_intptr_t), parameter         :: size_of_datatype = size_of_&
                                                                      &double&
                                                                      &_real
      integer(kind=ik), intent(in)                :: max_threads
      useIntelGPU = .false.
      if (useGPU) then
        if (gpu_vendor() == INTEL_GPU) then
          useIntelGPU = .true.
        endif
      endif

      call obj%timer%start("merge_systems" // "_double")
      success = .true.
      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      ! If my processor column isn't in the requested set, do nothing

      if (my_pcol<npc_0 .or. my_pcol>=npc_0+npc_n) then
        call obj%timer%stop("merge_systems" // "_double")
        return
      endif
      ! Determine number of "next" and "prev" column for ring sends

      if (my_pcol == npc_0+npc_n-1) then
        np_next = npc_0
      else
        np_next = my_pcol + 1
      endif

      if (my_pcol == npc_0) then
        np_prev = npc_0+npc_n-1
      else
        np_prev = my_pcol - 1
      endif
      call check_monotony_&
      &double&
      &(obj, nm,d,'Input1',wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_double")
        return
      endif
      call check_monotony_&
      &double&
      &(obj,na-nm,d(nm+1),'Input2',wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_double")
        return
      endif
      ! Get global number of processors and my processor number.
      ! Please note that my_proc does not need to match any real processor number,
      ! it is just used for load balancing some loops.

      n_procs = np_rows*npc_n
      my_proc = my_prow*npc_n + (my_pcol-npc_0) ! Row major


      ! Local limits of the rows of Q

      l_rqs = local_index(nqoff+1 , my_prow, np_rows, nblk, +1) ! First row of Q
      l_rqm = local_index(nqoff+nm, my_prow, np_rows, nblk, -1) ! Last row <= nm
      l_rqe = local_index(nqoff+na, my_prow, np_rows, nblk, -1) ! Last row of Q

      l_rnm  = l_rqm-l_rqs+1 ! Number of local rows <= nm
      l_rows = l_rqe-l_rqs+1 ! Total number of local rows


      ! My number of local columns

      l_cols = COUNT(p_col(1:na)==my_pcol)

      ! Get max number of local columns

      max_local_cols = 0
      do np = npc_0, npc_0+npc_n-1
        max_local_cols = MAX(max_local_cols,COUNT(p_col(1:na)==np))
      enddo

      ! Calculations start here

      beta = abs(e)
      sig  = sign(1.0_rk,e)

      ! Calculate rank-1 modifier z

      z(:) = 0

      if (MOD((nqoff+nm-1)/nblk,np_rows)==my_prow) then
        ! nm is local on my row
        do i = 1, na
          if (p_col(i)==my_pcol) z(i) = q(l_rqm,l_col(i))
         enddo
      endif

      if (MOD((nqoff+nm)/nblk,np_rows)==my_prow) then
        ! nm+1 is local on my row
        do i = 1, na
          if (p_col(i)==my_pcol) z(i) = z(i) + sig*q(l_rqm+1,l_col(i))
        enddo
      endif

      call global_gather_&
      &double&
      &(obj, z, na, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
      ! Normalize z so that norm(z) = 1.  Since z is the concatenation of
      ! two normalized vectors, norm2(z) = sqrt(2).
      z = z/sqrt(2.0_rk)
      rho = 2.0_rk*beta
      ! Calculate index for merging both systems by ascending eigenvalues
      call obj%timer%start("blas")
      call DLAMRG( int(nm,kind=BLAS_KIND), int(na-nm,kind=BLAS_KIND), d, &
                            1_BLAS_KIND, 1_BLAS_KIND, idxBLAS )
      idx(:) = int(idxBLAS(:),kind=ik)
      call obj%timer%stop("blas")

      ! Calculate the allowable deflation tolerance

      zmax = maxval(abs(z))
      dmax = maxval(abs(d))
      EPS = DLAMCH( 'E' ) ! return epsilon
      TOL = 8.0_rk*EPS*MAX(dmax,zmax)

      ! If the rank-1 modifier is small enough, no more needs to be done
      ! except to reorganize D and Q

      IF ( RHO*zmax <= TOL ) THEN

        ! Rearrange eigenvalues

        tmp = d
        do i=1,na
          d(i) = tmp(idx(i))
        enddo

        ! Rearrange eigenvectors
        call resort_ev_&
        &double &
                       (obj, idx, na, na, p_col_out, q, ldq, matrixCols, l_rows, l_rqe, &
                        l_rqs, mpi_comm_cols, p_col, l_col, l_col_out)

        call obj%timer%stop("merge_systems" // "_double")

        return
      ENDIF

      ! Merge and deflate system

      na1 = 0
      na2 = 0

      ! COLTYP:
      ! 1 : non-zero in the upper half only;
      ! 2 : dense;
      ! 3 : non-zero in the lower half only;
      ! 4 : deflated.

      coltyp(1:nm) = 1
      coltyp(nm+1:na) = 3

      do i=1,na

        if (rho*abs(z(idx(i))) <= tol) then

          ! Deflate due to small z component.

          na2 = na2+1
          d2(na2)   = d(idx(i))
          idx2(na2) = idx(i)
          coltyp(idx(i)) = 4

        else if (na1>0) then

          ! Check if eigenvalues are close enough to allow deflation.

          S = Z(idx(i))
          C = Z1(na1)

          ! Find sqrt(a**2+b**2) without overflow or
          ! destructive underflow.
          TAU = DLAPY2( C, S )
          T = D1(na1) - D(idx(i))
          C = C / TAU
          S = -S / TAU
          IF ( ABS( T*C*S ) <= TOL ) THEN

            ! Deflation is possible.

            na2 = na2+1

            Z1(na1) = TAU

            d2new = D(idx(i))*C**2 + D1(na1)*S**2
            d1new = D(idx(i))*S**2 + D1(na1)*C**2

            ! D(idx(i)) >= D1(na1) and C**2 + S**2 == 1.0
            ! This means that after the above transformation it must be
            !    D1(na1) <= d1new <= D(idx(i))
            !    D1(na1) <= d2new <= D(idx(i))
            !
            ! D1(na1) may get bigger but it is still smaller than the next D(idx(i+1))
            ! so there is no problem with sorting here.
            ! d2new <= D(idx(i)) which means that it might be smaller than D2(na2-1)
            ! which makes a check (and possibly a resort) necessary.
            !
            ! The above relations may not hold exactly due to numeric differences
            ! so they have to be enforced in order not to get troubles with sorting.


            if (d1new<D1(na1)  ) d1new = D1(na1)
            if (d1new>D(idx(i))) d1new = D(idx(i))

            if (d2new<D1(na1)  ) d2new = D1(na1)
            if (d2new>D(idx(i))) d2new = D(idx(i))

            D1(na1) = d1new

            do j=na2-1,1,-1
              if (d2new<d2(j)) then
                d2(j+1)   = d2(j)
                idx2(j+1) = idx2(j)
              else
                exit ! Loop
              endif
            enddo

            d2(j+1)   = d2new
            idx2(j+1) = idx(i)

            qtrans(1,1) = C; qtrans(1,2) =-S
            qtrans(2,1) = S; qtrans(2,2) = C
            call transform_columns_&
            &double &
                        (obj, idx(i), idx1(na1), na, tmp, l_rqs, l_rqe, &
                         q, ldq, matrixCols, l_rows, mpi_comm_cols, &
                          p_col, l_col, qtrans)
            if (coltyp(idx(i))==1 .and. coltyp(idx1(na1))/=1) coltyp(idx1(na1)) = 2
            if (coltyp(idx(i))==3 .and. coltyp(idx1(na1))/=3) coltyp(idx1(na1)) = 2

            coltyp(idx(i)) = 4

          else
            na1 = na1+1
            d1(na1) = d(idx(i))
            z1(na1) = z(idx(i))
            idx1(na1) = idx(i)
          endif
        else
          na1 = na1+1
          d1(na1) = d(idx(i))
          z1(na1) = z(idx(i))
          idx1(na1) = idx(i)
        endif

      enddo
      call check_monotony_&
      &double&
      &(obj, na1,d1,'Sorted1', wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_double")
        return
      endif
      call check_monotony_&
      &double&
      &(obj, na2,d2,'Sorted2', wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_double")
        return
      endif

      if (na1==1 .or. na1==2) then
        ! if(my_proc==0) print *,'--- Remark solve_tridi: na1==',na1,' proc==',myid

        if (na1==1) then
          d(1) = d1(1) + rho*z1(1)**2 ! solve secular equation
        else ! na1==2
          call obj%timer%start("blas")
          call DLAED5(1_BLAS_KIND, d1, z1, qtrans(1,1), rho, d(1))
          call DLAED5(2_BLAS_KIND, d1, z1, qtrans(1,2), rho, d(2))
          call obj%timer%stop("blas")
          call transform_columns_&
          &double&
          &(obj, idx1(1), idx1(2), na, tmp, l_rqs, l_rqe, q, &
            ldq, matrixCols, l_rows, mpi_comm_cols, &
             p_col, l_col, qtrans)

        endif

        ! Add the deflated eigenvalues
        d(na1+1:na) = d2(1:na2)

        ! Calculate arrangement of all eigenvalues  in output
        call obj%timer%start("blas")
        call DLAMRG( int(na1,kind=BLAS_KIND), int(na-na1,kind=BLAS_KIND), d, &
                              1_BLAS_KIND, 1_BLAS_KIND, idxBLAS )
        idx(:) = int(idxBLAS(:),kind=ik)
        call obj%timer%stop("blas")
        ! Rearrange eigenvalues

        tmp = d
        do i=1,na
          d(i) = tmp(idx(i))
        enddo

        ! Rearrange eigenvectors

        do i=1,na
          if (idx(i)<=na1) then
            idxq1(i) = idx1(idx(i))
          else
            idxq1(i) = idx2(idx(i)-na1)
          endif
        enddo
        call resort_ev_&
        &double&
        &(obj, idxq1, na, na, p_col_out, q, ldq, matrixCols, l_rows, l_rqe, &
          l_rqs, mpi_comm_cols, p_col, l_col, l_col_out)

      else if (na1>2) then

        ! Solve secular equation

        z(1:na1) = 1
        dbase(1:na1) = 0
        ddiff(1:na1) = 0

        info = 0
        infoBLAS = int(info,kind=BLAS_KIND)
!#ifdef WITH_OPENMP_TRADITIONAL
!
!        call obj%timer%start("OpenMP parallel" // "_double")
!!$OMP PARALLEL PRIVATE(i,my_thread,delta,s,info,infoBLAS,j)
!        my_thread = omp_get_thread_num()
!!$OMP DO
!#endif
        DO i = my_proc+1, na1, n_procs ! work distributed over all processors
          call obj%timer%start("blas")
          call DLAED4(int(na1,kind=BLAS_KIND), int(i,kind=BLAS_KIND), d1, z1, delta, &
                               rho, s, infoBLAS) ! s is not used!
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")
          if (info/=0) then
            ! If DLAED4 fails (may happen especially for LAPACK versions before 3.2)
            ! use the more stable bisection algorithm in solve_secular_equation
            ! print *,'ERROR DLAED4 n=',na1,'i=',i,' Using Bisection'
            call solve_secular_equation_&
            &double&
            &(obj, na1, i, d1, z1, delta, rho, s)
          endif

          ! Compute updated z

!#ifdef WITH_OPENMP_TRADITIONAL
!          do j=1,na1
!            if (i/=j)  z_p(j,my_thread) = z_p(j,my_thread)*( delta(j) / (d1(j)-d1(i)) )
!          enddo
!          z_p(i,my_thread) = z_p(i,my_thread)*delta(i)
!#else
          do j=1,na1
            if (i/=j)  z(j) = z(j)*( delta(j) / (d1(j)-d1(i)) )
          enddo
          z(i) = z(i)*delta(i)
!#endif
          ! store dbase/ddiff

          if (i<na1) then
            if (abs(delta(i+1)) < abs(delta(i))) then
              dbase(i) = d1(i+1)
              ddiff(i) = delta(i+1)
            else
              dbase(i) = d1(i)
              ddiff(i) = delta(i)
            endif
          else
            dbase(i) = d1(i)
            ddiff(i) = delta(i)
          endif
        enddo
!#ifdef WITH_OPENMP_TRADITIONAL
!!$OMP END PARALLEL
!
!        call obj%timer%stop("OpenMP parallel" // "_double")
!
!        do i = 0, max_threads-1
!          z(1:na1) = z(1:na1)*z_p(1:na1,i)
!        enddo
!#endif

        call global_product_&
        &double&
        (obj, z, na1, mpi_comm_rows, mpi_comm_cols, npc_0, npc_n)
        z(1:na1) = SIGN( SQRT( -z(1:na1) ), z1(1:na1) )

        call global_gather_&
        &double&
        &(obj, dbase, na1, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
        call global_gather_&
        &double&
        &(obj, ddiff, na1, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
        d(1:na1) = dbase(1:na1) - ddiff(1:na1)

        ! Calculate scale factors for eigenvectors
        ev_scale(:) = 0.0_rk

        DO i = my_proc+1, na1, n_procs ! work distributed over all processors

          ! tmp(1:na1) = z(1:na1) / delta(1:na1,i)  ! original code
          ! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

          ! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
          ! in exactly this order, but we want to prevent compiler optimization
!         ev_scale_val = ev_scale(i)
          call add_tmp_&
          &double&
          &(obj, d1, dbase, ddiff, z, ev_scale(i), na1,i)
!         ev_scale(i) = ev_scale_val
        enddo

        call global_gather_&
        &double&
        &(obj, ev_scale, na1, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
        ! Add the deflated eigenvalues
        d(na1+1:na) = d2(1:na2)

        call obj%timer%start("blas")
        ! Calculate arrangement of all eigenvalues  in output
        call DLAMRG(int(na1,kind=BLAS_KIND), int(na-na1,kind=BLAS_KIND), d, &
                             1_BLAS_KIND, 1_BLAS_KIND, idxBLAS )
        idx(:) = int(idxBLAS(:),kind=ik)
        call obj%timer%stop("blas")
        ! Rearrange eigenvalues
        tmp = d
        do i=1,na
          d(i) = tmp(idx(i))
        enddo
        call check_monotony_&
        &double&
        &(obj, na,d,'Output', wantDebug, success)

        if (.not.(success)) then
          call obj%timer%stop("merge_systems" // "_double")
          return
        endif
        ! Eigenvector calculations


        ! Calculate the number of columns in the new local matrix Q
        ! which are updated from non-deflated/deflated eigenvectors.
        ! idxq1/2 stores the global column numbers.

        nqcols1 = 0 ! number of non-deflated eigenvectors
        nqcols2 = 0 ! number of deflated eigenvectors
        DO i = 1, na
          if (p_col_out(i)==my_pcol) then
            if (idx(i)<=na1) then
              nqcols1 = nqcols1+1
              idxq1(nqcols1) = i
            else
              nqcols2 = nqcols2+1
              idxq2(nqcols2) = i
            endif
          endif
        enddo

        gemm_dim_k = MAX(1,l_rows)
        gemm_dim_l = max_local_cols
        gemm_dim_m = MIN(max_strip,MAX(1,nqcols1))

        allocate(qtmp1(gemm_dim_k, gemm_dim_l), stat=istat, errmsg=errorMessage)
        call check_allocate_f("merge_systems: qtmp1", 637, istat,  errorMessage)

        allocate(ev(gemm_dim_l,gemm_dim_m), stat=istat, errmsg=errorMessage)
        call check_allocate_f("merge_systems: ev", 640, istat,  errorMessage)

        allocate(qtmp2(gemm_dim_k, gemm_dim_m), stat=istat, errmsg=errorMessage)
        call check_allocate_f("merge_systems: qtmp2", 643, istat,  errorMessage)

        qtmp1 = 0 ! May contain empty (unset) parts
        qtmp2 = 0 ! Not really needed

        if (useGPU .and. .not.(useIntelGPU) ) then
          num = (gemm_dim_k * gemm_dim_l) * size_of_datatype
          successGPU = gpu_host_register(int(loc(qtmp1),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
          call check_host_register_GPU_f("merge_systems: qtmp1", 652,  successGPU)

          successGPU = gpu_malloc(qtmp1_dev, num)
          call check_alloc_GPU_f("merge_systems: qtmp1_dev", 655,  successGPU)

          num = (gemm_dim_l * gemm_dim_m) * size_of_datatype
          successGPU = gpu_host_register(int(loc(ev),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
          call check_host_register_GPU_f("merge_systems: ev", 660,  successGPU)

          successGPU = gpu_malloc(ev_dev, num)
          call check_alloc_GPU_f("merge_systems: ev_dev", 663,  successGPU)


          num = (gemm_dim_k * gemm_dim_m) * size_of_datatype
          successGPU = gpu_host_register(int(loc(qtmp2),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
          call check_host_register_GPU_f("merge_systems: qtmp2", 669,  successGPU)

          successGPU = gpu_malloc(qtmp2_dev, num)
          call check_alloc_GPU_f("merge_systems: qtmp2_dev", 672,  successGPU)
        endif

        !if (useIntelGPU) then
        !  ! needed later
        !endif

        ! Gather nonzero upper/lower components of old matrix Q
        ! which are needed for multiplication with new eigenvectors

        nnzu = 0
        nnzl = 0
        do i = 1, na1
          l_idx = l_col(idx1(i))
          if (p_col(idx1(i))==my_pcol) then
            if (coltyp(idx1(i))==1 .or. coltyp(idx1(i))==2) then
              nnzu = nnzu+1
              qtmp1(1:l_rnm,nnzu) = q(l_rqs:l_rqm,l_idx)
            endif
            if (coltyp(idx1(i))==3 .or. coltyp(idx1(i))==2) then
              nnzl = nnzl+1
              qtmp1(l_rnm+1:l_rows,nnzl) = q(l_rqm+1:l_rqe,l_idx)
            endif
          endif
        enddo

        ! Gather deflated eigenvalues behind nonzero components

        ndef = max(nnzu,nnzl)
        do i = 1, na2
          l_idx = l_col(idx2(i))
          if (p_col(idx2(i))==my_pcol) then
            ndef = ndef+1
            qtmp1(1:l_rows,ndef) = q(l_rqs:l_rqe,l_idx)
          endif
        enddo

        l_cols_qreorg = ndef ! Number of columns in reorganized matrix

        ! Set (output) Q to 0, it will sum up new Q

        DO i = 1, na
          if(p_col_out(i)==my_pcol) q(l_rqs:l_rqe,l_col_out(i)) = 0
        enddo

        np_rem = my_pcol

        do np = 1, npc_n
          ! Do a ring send of qtmp1

          if (np > 1) then

            if (np_rem == npc_0) then
              np_rem = npc_0+npc_n-1
            else
              np_rem = np_rem-1
            endif
            call obj%timer%start("mpi_communication")
            call MPI_Sendrecv_replace(qtmp1, int(l_rows*max_local_cols,kind=MPI_KIND), MPI_REAL8,     &
                                        int(np_next,kind=MPI_KIND), 1111_MPI_KIND, int(np_prev,kind=MPI_KIND), &
                                        1111_MPI_KIND, int(mpi_comm_cols,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            call obj%timer%stop("mpi_communication")
          endif

          if (useGPU .and. .not.(useIntelGPU)) then
            ! copy back after sendrecv
            successGPU = gpu_memcpy(qtmp1_dev, int(loc(qtmp1(1,1)),kind=c_intptr_t), &
                 gemm_dim_k * gemm_dim_l  * size_of_datatype, gpuMemcpyHostToDevice)
            call check_memcpy_GPU_f("merge_systems: qtmp1_dev", 742,  successGPU)
          endif

          !if (useIntelGPU) then
          !  ! needed later
          !endif


          ! Gather the parts in d1 and z which are fitting to qtmp1.
          ! This also delivers nnzu/nnzl for proc np_rem

          nnzu = 0
          nnzl = 0
          do i=1,na1
            if (p_col(idx1(i)) == np_rem) then
              if (coltyp(idx1(i)) == 1 .or. coltyp(idx1(i)) == 2) then
                nnzu = nnzu+1
                d1u(nnzu) = d1(i)
                zu (nnzu) = z (i)
              endif
              if (coltyp(idx1(i)) == 3 .or. coltyp(idx1(i)) == 2) then
                nnzl = nnzl+1
                d1l(nnzl) = d1(i)
                zl (nnzl) = z (i)
              endif
            endif
          enddo

          ! Set the deflated eigenvectors in Q (comming from proc np_rem)

          ndef = MAX(nnzu,nnzl) ! Remote counter in input matrix
          do i = 1, na
            j = idx(i)
            if (j>na1) then
              if (p_col(idx2(j-na1)) == np_rem) then
                ndef = ndef+1
                if (p_col_out(i) == my_pcol) &
                      q(l_rqs:l_rqe,l_col_out(i)) = qtmp1(1:l_rows,ndef)
              endif
            endif
          enddo

          do ns = 0, nqcols1-1, max_strip ! strimining loop

            ncnt = MIN(max_strip,nqcols1-ns) ! number of columns in this strip

            ! Get partial result from (output) Q

            do i = 1, ncnt
              qtmp2(1:l_rows,i) = q(l_rqs:l_rqe,l_col_out(idxq1(i+ns)))
            enddo

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with upper half of Q:

            do i = 1, ncnt
              j = idx(idxq1(i+ns))
              ! Calculate the j-th eigenvector of the deflated system
              ! See above why we are doing it this way!
              tmp(1:nnzu) = d1u(1:nnzu)-dbase(j)
              call v_add_s_&
              &double&
              &(obj,tmp,nnzu,ddiff(j))
              ev(1:nnzu,i) = zu(1:nnzu) / tmp(1:nnzu) * ev_scale(j)
            enddo

            if(useGPU .and. .not.(useIntelGPU) ) then
              !TODO: it should be enough to copy l_rows x ncnt
              ! copy to device
              successGPU = gpu_memcpy(qtmp2_dev, int(loc(qtmp2(1,1)),kind=c_intptr_t), &
                                 gemm_dim_k * gemm_dim_m * size_of_datatype, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("merge_systems: qtmp2_dev", 813,  successGPU)

              !TODO the previous loop could be possible to do on device and thus
              !copy less
              successGPU = gpu_memcpy(ev_dev, int(loc(ev(1,1)),kind=c_intptr_t), &
                                 gemm_dim_l * gemm_dim_m * size_of_datatype, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("merge_systems: ev_dev", 819,  successGPU)
            endif

            !if (useIntelGPU) then
            !  ! needed later
            !endif

            ! Multiply old Q with eigenvectors (upper half)

            if (l_rnm>0 .and. ncnt>0 .and. nnzu>0) then
              if (useGPU) then
                if (useIntelGPU) then
                  call obj%timer%start("mkl_offload")
                  call obj%timer%stop("mkl_offload")
                else
                  call obj%timer%start("gpublas")
                  call gpublas_DGEMM('N', 'N', l_rnm, ncnt, nnzu,   &
                                      1.0_rk, qtmp1_dev, ubound(qtmp1,dim=1),    &
                                      ev_dev, ubound(ev,dim=1), &
                                      1.0_rk, qtmp2_dev, ubound(qtmp2,dim=1))
                  call obj%timer%stop("gpublas")
                endif
              else
                call obj%timer%start("blas")
                call obj%timer%start("gemm")
                call DGEMM('N', 'N', int(l_rnm,kind=BLAS_KIND), int(ncnt,kind=BLAS_KIND), &
                                    int(nnzu,kind=BLAS_KIND),   &
                                    1.0_rk, qtmp1, int(ubound(qtmp1,dim=1),kind=BLAS_KIND),    &
                                    ev, int(ubound(ev,dim=1),kind=BLAS_KIND), &
                                    1.0_rk, qtmp2(1,1), int(ubound(qtmp2,dim=1),kind=BLAS_KIND))
                call obj%timer%stop("gemm")
                call obj%timer%stop("blas")
              endif ! useGPU
            endif

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with lower half of Q:

            do i = 1, ncnt
              j = idx(idxq1(i+ns))
              ! Calculate the j-th eigenvector of the deflated system
              ! See above why we are doing it this way!
              tmp(1:nnzl) = d1l(1:nnzl)-dbase(j)
              call v_add_s_&
              &double&
              &(obj,tmp,nnzl,ddiff(j))
              ev(1:nnzl,i) = zl(1:nnzl) / tmp(1:nnzl) * ev_scale(j)
            enddo

            if (useGPU .and. .not.(useIntelGPU) ) then
              !TODO the previous loop could be possible to do on device and thus
              !copy less
              successGPU = gpu_memcpy(ev_dev, int(loc(ev(1,1)),kind=c_intptr_t), &
                                 gemm_dim_l * gemm_dim_m * size_of_datatype, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("merge_systems: ev_dev", 880,  successGPU)
            endif

            !if (useIntelGPU) then
            !  ! needed later      
            !endif

            ! Multiply old Q with eigenvectors (lower half)

            if (l_rows-l_rnm>0 .and. ncnt>0 .and. nnzl>0) then
              if (useGPU) then
                if (useIntelGPU) then
                  call obj%timer%start("mkl_offload")
                  call obj%timer%stop("mkl_offload")

                else
                  call obj%timer%start("gpublas")
                  call gpublas_DGEMM('N', 'N', l_rows-l_rnm, ncnt, nnzl,   &
                                      1.0_rk, qtmp1_dev + l_rnm * size_of_datatype, ubound(qtmp1,dim=1),    &
                                      ev_dev, ubound(ev,dim=1), &
                                      1.0_rk, qtmp2_dev + l_rnm * size_of_datatype, ubound(qtmp2,dim=1))
                  call obj%timer%stop("gpublas")
                endif
              else
                call obj%timer%start("blas")
                call obj%timer%start("gemm")
                call DGEMM('N', 'N', int(l_rows-l_rnm,kind=BLAS_KIND), int(ncnt,kind=BLAS_KIND),  &
                                     int(nnzl,kind=BLAS_KIND),   &
                                     1.0_rk, qtmp1(l_rnm+1,1), int(ubound(qtmp1,dim=1),kind=BLAS_KIND),    &
                                     ev,  int(ubound(ev,dim=1),kind=BLAS_KIND),   &
                                     1.0_rk, qtmp2(l_rnm+1,1), int(ubound(qtmp2,dim=1),kind=BLAS_KIND))
                call obj%timer%stop("gemm")
                call obj%timer%stop("blas")
              endif ! useGPU
            endif

            if (useGPU .and. .not.(useIntelGPU) ) then
              !TODO either copy only half of the matrix here, and get rid of the
              !previous copy or copy whole array here

              ! COPY BACK
              successGPU = gpu_memcpy(int(loc(qtmp2(1,1)),kind=c_intptr_t), qtmp2_dev, &
                                 gemm_dim_k * gemm_dim_m * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("merge_systems: qtmp2_dev", 930,  successGPU)
            endif

            !if (useIntelGPU) then
            !  ! needed at a later time
            !endif


             ! Put partial result into (output) Q

            do i = 1, ncnt
              q(l_rqs:l_rqe,l_col_out(idxq1(i+ns))) = qtmp2(1:l_rows,i)
            enddo

          enddo   !ns = 0, nqcols1-1, max_strip ! strimining loop
        enddo    !do np = 1, npc_n

        if (useGPU .and. .not.(useIntelGPU) ) then
          successGPU = gpu_host_unregister(int(loc(qtmp1),kind=c_intptr_t))
          call check_host_unregister_GPU_f("merge_systems: qtmp1", 949,  successGPU)

          successGPU = gpu_free(qtmp1_dev)
          call check_dealloc_GPU_f("merge_systems: qtmp1_dev", 952,  successGPU)
          
          successGPU = gpu_host_unregister(int(loc(qtmp2),kind=c_intptr_t))
          call check_host_unregister_GPU_f("merge_systems: qtmp2", 955,  successGPU)

          successGPU = gpu_free(qtmp2_dev)
          call check_dealloc_GPU_f("merge_systems: qtmp2_dev", 958,  successGPU)

          successGPU = gpu_host_unregister(int(loc(ev),kind=c_intptr_t))
          call check_host_unregister_GPU_f("merge_systems: ev", 961,  successGPU)

          successGPU = gpu_free(ev_dev)
          call check_dealloc_GPU_f("merge_systems: ev_dev", 964,  successGPU)
        endif
        !if (useIntelGPU) then
        !  ! needed later
        !endif

        deallocate(ev, qtmp1, qtmp2, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("merge_systems: ev, qtmp1, qtmp2", 971, istat,  errorMessage)
      endif !very outer test (na1==1 .or. na1==2)

      call obj%timer%stop("merge_systems" // "_double")

      return

    end subroutine merge_systems_&
    &double

! real single precision first






















!cannot use "../src/solve_tridi/./../general/error_checking.inc" because filename with path can be too long for gfortran (max line length)


    subroutine merge_systems_&
    &single &
                         (obj, na, nm, d, e, q, ldq, nqoff, nblk, matrixCols, mpi_comm_rows, mpi_comm_cols, &
                          l_col, p_col, l_col_out, p_col_out, npc_0, npc_n, useGPU, wantDebug, success, max_threads)
      use elpa_gpu
      use, intrinsic :: iso_c_binding
      use precision
      use elpa_abstract_impl
      use elpa_blas_interfaces
      use global_product
      use global_gather
      use resort_ev
      use transform_columns
      use check_monotony
      use add_tmp
      use v_add_s
      use ELPA_utilities
      use elpa_mpi
      use solve_secular_equation
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

      class(elpa_abstract_impl_t), intent(inout)  :: obj
      integer(kind=ik), intent(in)                :: na, nm, ldq, nqoff, nblk, matrixCols, mpi_comm_rows, &
                                                     mpi_comm_cols, npc_0, npc_n
      integer(kind=ik), intent(in)                :: l_col(na), p_col(na), l_col_out(na), p_col_out(na)
      real(kind=rk4), intent(inout)     :: d(na), e
      real(kind=rk4), intent(inout)     :: q(ldq,*)
      logical, intent(in)                         :: useGPU, wantDebug
      logical                                     :: useIntelGPU

      logical, intent(out)                        :: success

      ! TODO: play with max_strip. If it was larger, matrices being multiplied
      ! might be larger as well!
      integer(kind=ik), parameter                 :: max_strip=128

      
      real(kind=rk4)                    :: beta, sig, s, c, t, tau, rho, eps, tol, &
                                                     qtrans(2,2), dmax, zmax, d1new, d2new
      real(kind=rk4)                    :: z(na), d1(na), d2(na), z1(na), delta(na),  &
                                                     dbase(na), ddiff(na), ev_scale(na), tmp(na)
      real(kind=rk4)                    :: d1u(na), zu(na), d1l(na), zl(na)
      real(kind=rk4), allocatable       :: qtmp1(:,:), qtmp2(:,:), ev(:,:)

      integer(kind=ik)                            :: i, j, na1, na2, l_rows, l_cols, l_rqs, l_rqe, &
                                                     l_rqm, ns, info
      integer(kind=BLAS_KIND)                     :: infoBLAS
      integer(kind=ik)                            :: l_rnm, nnzu, nnzl, ndef, ncnt, max_local_cols, &
                                                     l_cols_qreorg, np, l_idx, nqcols1, nqcols2
      integer(kind=ik)                            :: my_proc, n_procs, my_prow, my_pcol, np_rows, &
                                                     np_cols
      integer(kind=MPI_KIND)                      :: mpierr
      integer(kind=MPI_KIND)                      :: my_prowMPI, np_rowsMPI, my_pcolMPI, np_colsMPI
      integer(kind=ik)                            :: np_next, np_prev, np_rem
      integer(kind=ik)                            :: idx(na), idx1(na), idx2(na)
      integer(kind=BLAS_KIND)                     :: idxBLAS(NA)
      integer(kind=ik)                            :: coltyp(na), idxq1(na), idxq2(na)

      integer(kind=ik)                            :: istat
      character(200)                              :: errorMessage
      integer(kind=ik)                            :: gemm_dim_k, gemm_dim_l, gemm_dim_m

      integer(kind=c_intptr_t)                    :: num
      integer(kind=C_intptr_T)                    :: qtmp1_dev, qtmp2_dev, ev_dev
      logical                                     :: successGPU
      integer(kind=c_intptr_t), parameter         :: size_of_datatype = size_of_&
                                                                      &single&
                                                                      &_real
      integer(kind=ik), intent(in)                :: max_threads
      useIntelGPU = .false.
      if (useGPU) then
        if (gpu_vendor() == INTEL_GPU) then
          useIntelGPU = .true.
        endif
      endif

      call obj%timer%start("merge_systems" // "_single")
      success = .true.
      call obj%timer%start("mpi_communication")
      call mpi_comm_rank(int(mpi_comm_rows,kind=MPI_KIND) ,my_prowMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_rows,kind=MPI_KIND) ,np_rowsMPI, mpierr)
      call mpi_comm_rank(int(mpi_comm_cols,kind=MPI_KIND) ,my_pcolMPI, mpierr)
      call mpi_comm_size(int(mpi_comm_cols,kind=MPI_KIND) ,np_colsMPI, mpierr)

      my_prow = int(my_prowMPI,kind=c_int)
      np_rows = int(np_rowsMPI,kind=c_int)
      my_pcol = int(my_pcolMPI,kind=c_int)
      np_cols = int(np_colsMPI,kind=c_int)

      call obj%timer%stop("mpi_communication")

      ! If my processor column isn't in the requested set, do nothing

      if (my_pcol<npc_0 .or. my_pcol>=npc_0+npc_n) then
        call obj%timer%stop("merge_systems" // "_single")
        return
      endif
      ! Determine number of "next" and "prev" column for ring sends

      if (my_pcol == npc_0+npc_n-1) then
        np_next = npc_0
      else
        np_next = my_pcol + 1
      endif

      if (my_pcol == npc_0) then
        np_prev = npc_0+npc_n-1
      else
        np_prev = my_pcol - 1
      endif
      call check_monotony_&
      &single&
      &(obj, nm,d,'Input1',wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_single")
        return
      endif
      call check_monotony_&
      &single&
      &(obj,na-nm,d(nm+1),'Input2',wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_single")
        return
      endif
      ! Get global number of processors and my processor number.
      ! Please note that my_proc does not need to match any real processor number,
      ! it is just used for load balancing some loops.

      n_procs = np_rows*npc_n
      my_proc = my_prow*npc_n + (my_pcol-npc_0) ! Row major


      ! Local limits of the rows of Q

      l_rqs = local_index(nqoff+1 , my_prow, np_rows, nblk, +1) ! First row of Q
      l_rqm = local_index(nqoff+nm, my_prow, np_rows, nblk, -1) ! Last row <= nm
      l_rqe = local_index(nqoff+na, my_prow, np_rows, nblk, -1) ! Last row of Q

      l_rnm  = l_rqm-l_rqs+1 ! Number of local rows <= nm
      l_rows = l_rqe-l_rqs+1 ! Total number of local rows


      ! My number of local columns

      l_cols = COUNT(p_col(1:na)==my_pcol)

      ! Get max number of local columns

      max_local_cols = 0
      do np = npc_0, npc_0+npc_n-1
        max_local_cols = MAX(max_local_cols,COUNT(p_col(1:na)==np))
      enddo

      ! Calculations start here

      beta = abs(e)
      sig  = sign(1.0_rk,e)

      ! Calculate rank-1 modifier z

      z(:) = 0

      if (MOD((nqoff+nm-1)/nblk,np_rows)==my_prow) then
        ! nm is local on my row
        do i = 1, na
          if (p_col(i)==my_pcol) z(i) = q(l_rqm,l_col(i))
         enddo
      endif

      if (MOD((nqoff+nm)/nblk,np_rows)==my_prow) then
        ! nm+1 is local on my row
        do i = 1, na
          if (p_col(i)==my_pcol) z(i) = z(i) + sig*q(l_rqm+1,l_col(i))
        enddo
      endif

      call global_gather_&
      &single&
      &(obj, z, na, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
      ! Normalize z so that norm(z) = 1.  Since z is the concatenation of
      ! two normalized vectors, norm2(z) = sqrt(2).
      z = z/sqrt(2.0_rk)
      rho = 2.0_rk*beta
      ! Calculate index for merging both systems by ascending eigenvalues
      call obj%timer%start("blas")
      call SLAMRG( int(nm,kind=BLAS_KIND), int(na-nm,kind=BLAS_KIND), d, &
                            1_BLAS_KIND, 1_BLAS_KIND, idxBLAS )
      idx(:) = int(idxBLAS(:),kind=ik)
      call obj%timer%stop("blas")

      ! Calculate the allowable deflation tolerance

      zmax = maxval(abs(z))
      dmax = maxval(abs(d))
      EPS = SLAMCH( 'E' ) ! return epsilon
      TOL = 8.0_rk*EPS*MAX(dmax,zmax)

      ! If the rank-1 modifier is small enough, no more needs to be done
      ! except to reorganize D and Q

      IF ( RHO*zmax <= TOL ) THEN

        ! Rearrange eigenvalues

        tmp = d
        do i=1,na
          d(i) = tmp(idx(i))
        enddo

        ! Rearrange eigenvectors
        call resort_ev_&
        &single &
                       (obj, idx, na, na, p_col_out, q, ldq, matrixCols, l_rows, l_rqe, &
                        l_rqs, mpi_comm_cols, p_col, l_col, l_col_out)

        call obj%timer%stop("merge_systems" // "_single")

        return
      ENDIF

      ! Merge and deflate system

      na1 = 0
      na2 = 0

      ! COLTYP:
      ! 1 : non-zero in the upper half only;
      ! 2 : dense;
      ! 3 : non-zero in the lower half only;
      ! 4 : deflated.

      coltyp(1:nm) = 1
      coltyp(nm+1:na) = 3

      do i=1,na

        if (rho*abs(z(idx(i))) <= tol) then

          ! Deflate due to small z component.

          na2 = na2+1
          d2(na2)   = d(idx(i))
          idx2(na2) = idx(i)
          coltyp(idx(i)) = 4

        else if (na1>0) then

          ! Check if eigenvalues are close enough to allow deflation.

          S = Z(idx(i))
          C = Z1(na1)

          ! Find sqrt(a**2+b**2) without overflow or
          ! destructive underflow.
          TAU = SLAPY2( C, S )
          T = D1(na1) - D(idx(i))
          C = C / TAU
          S = -S / TAU
          IF ( ABS( T*C*S ) <= TOL ) THEN

            ! Deflation is possible.

            na2 = na2+1

            Z1(na1) = TAU

            d2new = D(idx(i))*C**2 + D1(na1)*S**2
            d1new = D(idx(i))*S**2 + D1(na1)*C**2

            ! D(idx(i)) >= D1(na1) and C**2 + S**2 == 1.0
            ! This means that after the above transformation it must be
            !    D1(na1) <= d1new <= D(idx(i))
            !    D1(na1) <= d2new <= D(idx(i))
            !
            ! D1(na1) may get bigger but it is still smaller than the next D(idx(i+1))
            ! so there is no problem with sorting here.
            ! d2new <= D(idx(i)) which means that it might be smaller than D2(na2-1)
            ! which makes a check (and possibly a resort) necessary.
            !
            ! The above relations may not hold exactly due to numeric differences
            ! so they have to be enforced in order not to get troubles with sorting.


            if (d1new<D1(na1)  ) d1new = D1(na1)
            if (d1new>D(idx(i))) d1new = D(idx(i))

            if (d2new<D1(na1)  ) d2new = D1(na1)
            if (d2new>D(idx(i))) d2new = D(idx(i))

            D1(na1) = d1new

            do j=na2-1,1,-1
              if (d2new<d2(j)) then
                d2(j+1)   = d2(j)
                idx2(j+1) = idx2(j)
              else
                exit ! Loop
              endif
            enddo

            d2(j+1)   = d2new
            idx2(j+1) = idx(i)

            qtrans(1,1) = C; qtrans(1,2) =-S
            qtrans(2,1) = S; qtrans(2,2) = C
            call transform_columns_&
            &single &
                        (obj, idx(i), idx1(na1), na, tmp, l_rqs, l_rqe, &
                         q, ldq, matrixCols, l_rows, mpi_comm_cols, &
                          p_col, l_col, qtrans)
            if (coltyp(idx(i))==1 .and. coltyp(idx1(na1))/=1) coltyp(idx1(na1)) = 2
            if (coltyp(idx(i))==3 .and. coltyp(idx1(na1))/=3) coltyp(idx1(na1)) = 2

            coltyp(idx(i)) = 4

          else
            na1 = na1+1
            d1(na1) = d(idx(i))
            z1(na1) = z(idx(i))
            idx1(na1) = idx(i)
          endif
        else
          na1 = na1+1
          d1(na1) = d(idx(i))
          z1(na1) = z(idx(i))
          idx1(na1) = idx(i)
        endif

      enddo
      call check_monotony_&
      &single&
      &(obj, na1,d1,'Sorted1', wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_single")
        return
      endif
      call check_monotony_&
      &single&
      &(obj, na2,d2,'Sorted2', wantDebug, success)
      if (.not.(success)) then
        call obj%timer%stop("merge_systems" // "_single")
        return
      endif

      if (na1==1 .or. na1==2) then
        ! if(my_proc==0) print *,'--- Remark solve_tridi: na1==',na1,' proc==',myid

        if (na1==1) then
          d(1) = d1(1) + rho*z1(1)**2 ! solve secular equation
        else ! na1==2
          call obj%timer%start("blas")
          call SLAED5(1_BLAS_KIND, d1, z1, qtrans(1,1), rho, d(1))
          call SLAED5(2_BLAS_KIND, d1, z1, qtrans(1,2), rho, d(2))
          call obj%timer%stop("blas")
          call transform_columns_&
          &single&
          &(obj, idx1(1), idx1(2), na, tmp, l_rqs, l_rqe, q, &
            ldq, matrixCols, l_rows, mpi_comm_cols, &
             p_col, l_col, qtrans)

        endif

        ! Add the deflated eigenvalues
        d(na1+1:na) = d2(1:na2)

        ! Calculate arrangement of all eigenvalues  in output
        call obj%timer%start("blas")
        call SLAMRG( int(na1,kind=BLAS_KIND), int(na-na1,kind=BLAS_KIND), d, &
                              1_BLAS_KIND, 1_BLAS_KIND, idxBLAS )
        idx(:) = int(idxBLAS(:),kind=ik)
        call obj%timer%stop("blas")
        ! Rearrange eigenvalues

        tmp = d
        do i=1,na
          d(i) = tmp(idx(i))
        enddo

        ! Rearrange eigenvectors

        do i=1,na
          if (idx(i)<=na1) then
            idxq1(i) = idx1(idx(i))
          else
            idxq1(i) = idx2(idx(i)-na1)
          endif
        enddo
        call resort_ev_&
        &single&
        &(obj, idxq1, na, na, p_col_out, q, ldq, matrixCols, l_rows, l_rqe, &
          l_rqs, mpi_comm_cols, p_col, l_col, l_col_out)

      else if (na1>2) then

        ! Solve secular equation

        z(1:na1) = 1
        dbase(1:na1) = 0
        ddiff(1:na1) = 0

        info = 0
        infoBLAS = int(info,kind=BLAS_KIND)
!#ifdef WITH_OPENMP_TRADITIONAL
!
!        call obj%timer%start("OpenMP parallel" // "_single")
!!$OMP PARALLEL PRIVATE(i,my_thread,delta,s,info,infoBLAS,j)
!        my_thread = omp_get_thread_num()
!!$OMP DO
!#endif
        DO i = my_proc+1, na1, n_procs ! work distributed over all processors
          call obj%timer%start("blas")
          call SLAED4(int(na1,kind=BLAS_KIND), int(i,kind=BLAS_KIND), d1, z1, delta, &
                               rho, s, infoBLAS) ! s is not used!
          info = int(infoBLAS,kind=ik)
          call obj%timer%stop("blas")
          if (info/=0) then
            ! If DLAED4 fails (may happen especially for LAPACK versions before 3.2)
            ! use the more stable bisection algorithm in solve_secular_equation
            ! print *,'ERROR DLAED4 n=',na1,'i=',i,' Using Bisection'
            call solve_secular_equation_&
            &single&
            &(obj, na1, i, d1, z1, delta, rho, s)
          endif

          ! Compute updated z

!#ifdef WITH_OPENMP_TRADITIONAL
!          do j=1,na1
!            if (i/=j)  z_p(j,my_thread) = z_p(j,my_thread)*( delta(j) / (d1(j)-d1(i)) )
!          enddo
!          z_p(i,my_thread) = z_p(i,my_thread)*delta(i)
!#else
          do j=1,na1
            if (i/=j)  z(j) = z(j)*( delta(j) / (d1(j)-d1(i)) )
          enddo
          z(i) = z(i)*delta(i)
!#endif
          ! store dbase/ddiff

          if (i<na1) then
            if (abs(delta(i+1)) < abs(delta(i))) then
              dbase(i) = d1(i+1)
              ddiff(i) = delta(i+1)
            else
              dbase(i) = d1(i)
              ddiff(i) = delta(i)
            endif
          else
            dbase(i) = d1(i)
            ddiff(i) = delta(i)
          endif
        enddo
!#ifdef WITH_OPENMP_TRADITIONAL
!!$OMP END PARALLEL
!
!        call obj%timer%stop("OpenMP parallel" // "_single")
!
!        do i = 0, max_threads-1
!          z(1:na1) = z(1:na1)*z_p(1:na1,i)
!        enddo
!#endif

        call global_product_&
        &single&
        (obj, z, na1, mpi_comm_rows, mpi_comm_cols, npc_0, npc_n)
        z(1:na1) = SIGN( SQRT( -z(1:na1) ), z1(1:na1) )

        call global_gather_&
        &single&
        &(obj, dbase, na1, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
        call global_gather_&
        &single&
        &(obj, ddiff, na1, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
        d(1:na1) = dbase(1:na1) - ddiff(1:na1)

        ! Calculate scale factors for eigenvectors
        ev_scale(:) = 0.0_rk

        DO i = my_proc+1, na1, n_procs ! work distributed over all processors

          ! tmp(1:na1) = z(1:na1) / delta(1:na1,i)  ! original code
          ! tmp(1:na1) = z(1:na1) / (d1(1:na1)-d(i))! bad results

          ! All we want to calculate is tmp = (d1(1:na1)-dbase(i))+ddiff(i)
          ! in exactly this order, but we want to prevent compiler optimization
!         ev_scale_val = ev_scale(i)
          call add_tmp_&
          &single&
          &(obj, d1, dbase, ddiff, z, ev_scale(i), na1,i)
!         ev_scale(i) = ev_scale_val
        enddo

        call global_gather_&
        &single&
        &(obj, ev_scale, na1, mpi_comm_rows, mpi_comm_cols, npc_n, np_prev, np_next)
        ! Add the deflated eigenvalues
        d(na1+1:na) = d2(1:na2)

        call obj%timer%start("blas")
        ! Calculate arrangement of all eigenvalues  in output
        call SLAMRG(int(na1,kind=BLAS_KIND), int(na-na1,kind=BLAS_KIND), d, &
                             1_BLAS_KIND, 1_BLAS_KIND, idxBLAS )
        idx(:) = int(idxBLAS(:),kind=ik)
        call obj%timer%stop("blas")
        ! Rearrange eigenvalues
        tmp = d
        do i=1,na
          d(i) = tmp(idx(i))
        enddo
        call check_monotony_&
        &single&
        &(obj, na,d,'Output', wantDebug, success)

        if (.not.(success)) then
          call obj%timer%stop("merge_systems" // "_single")
          return
        endif
        ! Eigenvector calculations


        ! Calculate the number of columns in the new local matrix Q
        ! which are updated from non-deflated/deflated eigenvectors.
        ! idxq1/2 stores the global column numbers.

        nqcols1 = 0 ! number of non-deflated eigenvectors
        nqcols2 = 0 ! number of deflated eigenvectors
        DO i = 1, na
          if (p_col_out(i)==my_pcol) then
            if (idx(i)<=na1) then
              nqcols1 = nqcols1+1
              idxq1(nqcols1) = i
            else
              nqcols2 = nqcols2+1
              idxq2(nqcols2) = i
            endif
          endif
        enddo

        gemm_dim_k = MAX(1,l_rows)
        gemm_dim_l = max_local_cols
        gemm_dim_m = MIN(max_strip,MAX(1,nqcols1))

        allocate(qtmp1(gemm_dim_k, gemm_dim_l), stat=istat, errmsg=errorMessage)
        call check_allocate_f("merge_systems: qtmp1", 637, istat,  errorMessage)

        allocate(ev(gemm_dim_l,gemm_dim_m), stat=istat, errmsg=errorMessage)
        call check_allocate_f("merge_systems: ev", 640, istat,  errorMessage)

        allocate(qtmp2(gemm_dim_k, gemm_dim_m), stat=istat, errmsg=errorMessage)
        call check_allocate_f("merge_systems: qtmp2", 643, istat,  errorMessage)

        qtmp1 = 0 ! May contain empty (unset) parts
        qtmp2 = 0 ! Not really needed

        if (useGPU .and. .not.(useIntelGPU) ) then
          num = (gemm_dim_k * gemm_dim_l) * size_of_datatype
          successGPU = gpu_host_register(int(loc(qtmp1),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
          call check_host_register_GPU_f("merge_systems: qtmp1", 652,  successGPU)

          successGPU = gpu_malloc(qtmp1_dev, num)
          call check_alloc_GPU_f("merge_systems: qtmp1_dev", 655,  successGPU)

          num = (gemm_dim_l * gemm_dim_m) * size_of_datatype
          successGPU = gpu_host_register(int(loc(ev),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
          call check_host_register_GPU_f("merge_systems: ev", 660,  successGPU)

          successGPU = gpu_malloc(ev_dev, num)
          call check_alloc_GPU_f("merge_systems: ev_dev", 663,  successGPU)


          num = (gemm_dim_k * gemm_dim_m) * size_of_datatype
          successGPU = gpu_host_register(int(loc(qtmp2),kind=c_intptr_t),num,&
                        gpuHostRegisterDefault)
          call check_host_register_GPU_f("merge_systems: qtmp2", 669,  successGPU)

          successGPU = gpu_malloc(qtmp2_dev, num)
          call check_alloc_GPU_f("merge_systems: qtmp2_dev", 672,  successGPU)
        endif

        !if (useIntelGPU) then
        !  ! needed later
        !endif

        ! Gather nonzero upper/lower components of old matrix Q
        ! which are needed for multiplication with new eigenvectors

        nnzu = 0
        nnzl = 0
        do i = 1, na1
          l_idx = l_col(idx1(i))
          if (p_col(idx1(i))==my_pcol) then
            if (coltyp(idx1(i))==1 .or. coltyp(idx1(i))==2) then
              nnzu = nnzu+1
              qtmp1(1:l_rnm,nnzu) = q(l_rqs:l_rqm,l_idx)
            endif
            if (coltyp(idx1(i))==3 .or. coltyp(idx1(i))==2) then
              nnzl = nnzl+1
              qtmp1(l_rnm+1:l_rows,nnzl) = q(l_rqm+1:l_rqe,l_idx)
            endif
          endif
        enddo

        ! Gather deflated eigenvalues behind nonzero components

        ndef = max(nnzu,nnzl)
        do i = 1, na2
          l_idx = l_col(idx2(i))
          if (p_col(idx2(i))==my_pcol) then
            ndef = ndef+1
            qtmp1(1:l_rows,ndef) = q(l_rqs:l_rqe,l_idx)
          endif
        enddo

        l_cols_qreorg = ndef ! Number of columns in reorganized matrix

        ! Set (output) Q to 0, it will sum up new Q

        DO i = 1, na
          if(p_col_out(i)==my_pcol) q(l_rqs:l_rqe,l_col_out(i)) = 0
        enddo

        np_rem = my_pcol

        do np = 1, npc_n
          ! Do a ring send of qtmp1

          if (np > 1) then

            if (np_rem == npc_0) then
              np_rem = npc_0+npc_n-1
            else
              np_rem = np_rem-1
            endif
            call obj%timer%start("mpi_communication")
            call MPI_Sendrecv_replace(qtmp1, int(l_rows*max_local_cols,kind=MPI_KIND), MPI_REAL4,     &
                                        int(np_next,kind=MPI_KIND), 1111_MPI_KIND, int(np_prev,kind=MPI_KIND), &
                                        1111_MPI_KIND, int(mpi_comm_cols,kind=MPI_KIND), MPI_STATUS_IGNORE, mpierr)
            call obj%timer%stop("mpi_communication")
          endif

          if (useGPU .and. .not.(useIntelGPU)) then
            ! copy back after sendrecv
            successGPU = gpu_memcpy(qtmp1_dev, int(loc(qtmp1(1,1)),kind=c_intptr_t), &
                 gemm_dim_k * gemm_dim_l  * size_of_datatype, gpuMemcpyHostToDevice)
            call check_memcpy_GPU_f("merge_systems: qtmp1_dev", 742,  successGPU)
          endif

          !if (useIntelGPU) then
          !  ! needed later
          !endif


          ! Gather the parts in d1 and z which are fitting to qtmp1.
          ! This also delivers nnzu/nnzl for proc np_rem

          nnzu = 0
          nnzl = 0
          do i=1,na1
            if (p_col(idx1(i)) == np_rem) then
              if (coltyp(idx1(i)) == 1 .or. coltyp(idx1(i)) == 2) then
                nnzu = nnzu+1
                d1u(nnzu) = d1(i)
                zu (nnzu) = z (i)
              endif
              if (coltyp(idx1(i)) == 3 .or. coltyp(idx1(i)) == 2) then
                nnzl = nnzl+1
                d1l(nnzl) = d1(i)
                zl (nnzl) = z (i)
              endif
            endif
          enddo

          ! Set the deflated eigenvectors in Q (comming from proc np_rem)

          ndef = MAX(nnzu,nnzl) ! Remote counter in input matrix
          do i = 1, na
            j = idx(i)
            if (j>na1) then
              if (p_col(idx2(j-na1)) == np_rem) then
                ndef = ndef+1
                if (p_col_out(i) == my_pcol) &
                      q(l_rqs:l_rqe,l_col_out(i)) = qtmp1(1:l_rows,ndef)
              endif
            endif
          enddo

          do ns = 0, nqcols1-1, max_strip ! strimining loop

            ncnt = MIN(max_strip,nqcols1-ns) ! number of columns in this strip

            ! Get partial result from (output) Q

            do i = 1, ncnt
              qtmp2(1:l_rows,i) = q(l_rqs:l_rqe,l_col_out(idxq1(i+ns)))
            enddo

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with upper half of Q:

            do i = 1, ncnt
              j = idx(idxq1(i+ns))
              ! Calculate the j-th eigenvector of the deflated system
              ! See above why we are doing it this way!
              tmp(1:nnzu) = d1u(1:nnzu)-dbase(j)
              call v_add_s_&
              &single&
              &(obj,tmp,nnzu,ddiff(j))
              ev(1:nnzu,i) = zu(1:nnzu) / tmp(1:nnzu) * ev_scale(j)
            enddo

            if(useGPU .and. .not.(useIntelGPU) ) then
              !TODO: it should be enough to copy l_rows x ncnt
              ! copy to device
              successGPU = gpu_memcpy(qtmp2_dev, int(loc(qtmp2(1,1)),kind=c_intptr_t), &
                                 gemm_dim_k * gemm_dim_m * size_of_datatype, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("merge_systems: qtmp2_dev", 813,  successGPU)

              !TODO the previous loop could be possible to do on device and thus
              !copy less
              successGPU = gpu_memcpy(ev_dev, int(loc(ev(1,1)),kind=c_intptr_t), &
                                 gemm_dim_l * gemm_dim_m * size_of_datatype, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("merge_systems: ev_dev", 819,  successGPU)
            endif

            !if (useIntelGPU) then
            !  ! needed later
            !endif

            ! Multiply old Q with eigenvectors (upper half)

            if (l_rnm>0 .and. ncnt>0 .and. nnzu>0) then
              if (useGPU) then
                if (useIntelGPU) then
                  call obj%timer%start("mkl_offload")
                  call obj%timer%stop("mkl_offload")
                else
                  call obj%timer%start("gpublas")
                  call gpublas_SGEMM('N', 'N', l_rnm, ncnt, nnzu,   &
                                      1.0_rk, qtmp1_dev, ubound(qtmp1,dim=1),    &
                                      ev_dev, ubound(ev,dim=1), &
                                      1.0_rk, qtmp2_dev, ubound(qtmp2,dim=1))
                  call obj%timer%stop("gpublas")
                endif
              else
                call obj%timer%start("blas")
                call obj%timer%start("gemm")
                call SGEMM('N', 'N', int(l_rnm,kind=BLAS_KIND), int(ncnt,kind=BLAS_KIND), &
                                    int(nnzu,kind=BLAS_KIND),   &
                                    1.0_rk, qtmp1, int(ubound(qtmp1,dim=1),kind=BLAS_KIND),    &
                                    ev, int(ubound(ev,dim=1),kind=BLAS_KIND), &
                                    1.0_rk, qtmp2(1,1), int(ubound(qtmp2,dim=1),kind=BLAS_KIND))
                call obj%timer%stop("gemm")
                call obj%timer%stop("blas")
              endif ! useGPU
            endif

            ! Compute eigenvectors of the rank-1 modified matrix.
            ! Parts for multiplying with lower half of Q:

            do i = 1, ncnt
              j = idx(idxq1(i+ns))
              ! Calculate the j-th eigenvector of the deflated system
              ! See above why we are doing it this way!
              tmp(1:nnzl) = d1l(1:nnzl)-dbase(j)
              call v_add_s_&
              &single&
              &(obj,tmp,nnzl,ddiff(j))
              ev(1:nnzl,i) = zl(1:nnzl) / tmp(1:nnzl) * ev_scale(j)
            enddo

            if (useGPU .and. .not.(useIntelGPU) ) then
              !TODO the previous loop could be possible to do on device and thus
              !copy less
              successGPU = gpu_memcpy(ev_dev, int(loc(ev(1,1)),kind=c_intptr_t), &
                                 gemm_dim_l * gemm_dim_m * size_of_datatype, gpuMemcpyHostToDevice)
              call check_memcpy_GPU_f("merge_systems: ev_dev", 880,  successGPU)
            endif

            !if (useIntelGPU) then
            !  ! needed later      
            !endif

            ! Multiply old Q with eigenvectors (lower half)

            if (l_rows-l_rnm>0 .and. ncnt>0 .and. nnzl>0) then
              if (useGPU) then
                if (useIntelGPU) then
                  call obj%timer%start("mkl_offload")
                  call obj%timer%stop("mkl_offload")

                else
                  call obj%timer%start("gpublas")
                  call gpublas_SGEMM('N', 'N', l_rows-l_rnm, ncnt, nnzl,   &
                                      1.0_rk, qtmp1_dev + l_rnm * size_of_datatype, ubound(qtmp1,dim=1),    &
                                      ev_dev, ubound(ev,dim=1), &
                                      1.0_rk, qtmp2_dev + l_rnm * size_of_datatype, ubound(qtmp2,dim=1))
                  call obj%timer%stop("gpublas")
                endif
              else
                call obj%timer%start("blas")
                call obj%timer%start("gemm")
                call SGEMM('N', 'N', int(l_rows-l_rnm,kind=BLAS_KIND), int(ncnt,kind=BLAS_KIND),  &
                                     int(nnzl,kind=BLAS_KIND),   &
                                     1.0_rk, qtmp1(l_rnm+1,1), int(ubound(qtmp1,dim=1),kind=BLAS_KIND),    &
                                     ev,  int(ubound(ev,dim=1),kind=BLAS_KIND),   &
                                     1.0_rk, qtmp2(l_rnm+1,1), int(ubound(qtmp2,dim=1),kind=BLAS_KIND))
                call obj%timer%stop("gemm")
                call obj%timer%stop("blas")
              endif ! useGPU
            endif

            if (useGPU .and. .not.(useIntelGPU) ) then
              !TODO either copy only half of the matrix here, and get rid of the
              !previous copy or copy whole array here

              ! COPY BACK
              successGPU = gpu_memcpy(int(loc(qtmp2(1,1)),kind=c_intptr_t), qtmp2_dev, &
                                 gemm_dim_k * gemm_dim_m * size_of_datatype, gpuMemcpyDeviceToHost)
              call check_memcpy_GPU_f("merge_systems: qtmp2_dev", 930,  successGPU)
            endif

            !if (useIntelGPU) then
            !  ! needed at a later time
            !endif


             ! Put partial result into (output) Q

            do i = 1, ncnt
              q(l_rqs:l_rqe,l_col_out(idxq1(i+ns))) = qtmp2(1:l_rows,i)
            enddo

          enddo   !ns = 0, nqcols1-1, max_strip ! strimining loop
        enddo    !do np = 1, npc_n

        if (useGPU .and. .not.(useIntelGPU) ) then
          successGPU = gpu_host_unregister(int(loc(qtmp1),kind=c_intptr_t))
          call check_host_unregister_GPU_f("merge_systems: qtmp1", 949,  successGPU)

          successGPU = gpu_free(qtmp1_dev)
          call check_dealloc_GPU_f("merge_systems: qtmp1_dev", 952,  successGPU)
          
          successGPU = gpu_host_unregister(int(loc(qtmp2),kind=c_intptr_t))
          call check_host_unregister_GPU_f("merge_systems: qtmp2", 955,  successGPU)

          successGPU = gpu_free(qtmp2_dev)
          call check_dealloc_GPU_f("merge_systems: qtmp2_dev", 958,  successGPU)

          successGPU = gpu_host_unregister(int(loc(ev),kind=c_intptr_t))
          call check_host_unregister_GPU_f("merge_systems: ev", 961,  successGPU)

          successGPU = gpu_free(ev_dev)
          call check_dealloc_GPU_f("merge_systems: ev_dev", 964,  successGPU)
        endif
        !if (useIntelGPU) then
        !  ! needed later
        !endif

        deallocate(ev, qtmp1, qtmp2, stat=istat, errmsg=errorMessage)
        call check_deallocate_f("merge_systems: ev, qtmp1, qtmp2", 971, istat,  errorMessage)
      endif !very outer test (na1==1 .or. na1==2)

      call obj%timer%stop("merge_systems" // "_single")

      return

    end subroutine merge_systems_&
    &single

end module
