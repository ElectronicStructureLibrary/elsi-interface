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
! Author: Andreas Marek, MPCDF

module elpa_skewsymmetric_blas
   use precision
   use, intrinsic :: iso_c_binding
contains

   subroutine elpa_dssr2(n, x, y,  a, lda )

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none
      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: rck = C_DOUBLE
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      integer(kind=BLAS_KIND)     :: n, lda
      real(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      real(kind=rck)     :: temp1, temp2
      integer(kind=ik)            :: i, j, ii, jj, ic, ix, iy, jc, jx, jy, info
      logical                     :: upper

      ! test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 5
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access A in lower triangular part.
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         jy = 1 + (jj-1)

         do j = 1, jc-1
            ! Do local update for blocks on the diagonal
            if ( ( x( jx + j -1) /= zero ) .or. &
               ( y( jy + j -1 ) /= zero ) ) then
               temp1 = - y( jy + j - 1 )
               temp2 = - x( jy + j - 1 )
               do i = j+1, jc
                  a( jj +  i -1 , jj +  j -1 ) = a(jj + i -1,jj +  j -1 ) + x( jx + i  -1)*temp1 - y(jj +  i -1 )*temp2
               end do
            end if
         end do

         ! Use dger for other blocks
         do ii = jj+nb, n, nb
            ic = min( nb, n-ii+1 )
            ix = 1 + (ii-1)
            iy = 1 + (ii-1)
            call DGER(int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -one, x( ix ), 1_BLAS_KIND, y( jy ), 1_BLAS_KIND, &
               a( ii, jj ), int(lda,kind=BLAS_KIND) )
            call DGER(int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), one, y( iy ), 1_BLAS_KIND, x( jx ), 1_BLAS_KIND, &
               a( ii, jj ), int(lda,kind=BLAS_KIND) )
         end do
      end do

      return
   end subroutine
   subroutine elpa_dssmv(n, alpha, a, lda, x,  y)

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none
      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: rck = C_DOUBLE
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      integer(kind=BLAS_KIND)     :: n, lda
      real(kind=rck)     :: alpha
      real(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      integer(kind=ik)            :: ii, jj, ic, iy, jc, jx, info
      real(kind=rck)     :: temp
      real(kind=rck)     :: work( nb )

      ! Test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 4
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access only lower triangular part of a

      temp = zero
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         do ii = 1, n, nb
            ic = min( nb, n-ii+1 )
            iy = 1 + (ii-1)

            ! gemv for non-diagonal blocks. use 2x dtrmv for diagonal blocks
            if ( ii < jj ) then
               call DGEMV('t', int(jc,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -alpha, &
                  a( jj, ii ), int(lda, kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else if ( ii > jj ) then
               call DGEMV('n', int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), alpha, a( ii, jj ), &
                  int(lda,kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else
               if (temp == zero) then
                  y(1:n) = zero
               else if (temp /= one) then
                  ! should not happen
                  call DSCAL( int(jc,kind=BLAS_KIND), temp, y( iy ), 1_BLAS_KIND)
               end if
               call DCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call DTRMV( 'l', 'n', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call DAXPY( int(jc,kind=BLAS_KIND),alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)

               call DCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call DTRMV( 'l', 't', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call DAXPY(int(jc,kind=BLAS_KIND), -alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)
            end if
         end do
         temp = one
      end do

      return
   end subroutine

   subroutine elpa_sssr2(n, x, y,  a, lda )

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none
      integer, parameter :: rk = C_FLOAT
      integer, parameter :: rck = C_FLOAT
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      integer(kind=BLAS_KIND)     :: n, lda
      real(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      real(kind=rck)     :: temp1, temp2
      integer(kind=ik)            :: i, j, ii, jj, ic, ix, iy, jc, jx, jy, info
      logical                     :: upper

      ! test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 5
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access A in lower triangular part.
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         jy = 1 + (jj-1)

         do j = 1, jc-1
            ! Do local update for blocks on the diagonal
            if ( ( x( jx + j -1) /= zero ) .or. &
               ( y( jy + j -1 ) /= zero ) ) then
               temp1 = - y( jy + j - 1 )
               temp2 = - x( jy + j - 1 )
               do i = j+1, jc
                  a( jj +  i -1 , jj +  j -1 ) = a(jj + i -1,jj +  j -1 ) + x( jx + i  -1)*temp1 - y(jj +  i -1 )*temp2
               end do
            end if
         end do

         ! Use dger for other blocks
         do ii = jj+nb, n, nb
            ic = min( nb, n-ii+1 )
            ix = 1 + (ii-1)
            iy = 1 + (ii-1)
            call SGER(int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -one, x( ix ), 1_BLAS_KIND, y( jy ), 1_BLAS_KIND, &
               a( ii, jj ), int(lda,kind=BLAS_KIND) )
            call SGER(int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), one, y( iy ), 1_BLAS_KIND, x( jx ), 1_BLAS_KIND, &
               a( ii, jj ), int(lda,kind=BLAS_KIND) )
         end do
      end do

      return
   end subroutine
   subroutine elpa_sssmv(n, alpha, a, lda, x,  y)

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none
      integer, parameter :: rk = C_FLOAT
      integer, parameter :: rck = C_FLOAT
      real(kind=rck), parameter      :: ZERO=0.0_rk, ONE = 1.0_rk

      integer(kind=BLAS_KIND)     :: n, lda
      real(kind=rck)     :: alpha
      real(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      integer(kind=ik)            :: ii, jj, ic, iy, jc, jx, info
      real(kind=rck)     :: temp
      real(kind=rck)     :: work( nb )

      ! Test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 4
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access only lower triangular part of a

      temp = zero
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         do ii = 1, n, nb
            ic = min( nb, n-ii+1 )
            iy = 1 + (ii-1)

            ! gemv for non-diagonal blocks. use 2x dtrmv for diagonal blocks
            if ( ii < jj ) then
               call SGEMV('t', int(jc,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -alpha, &
                  a( jj, ii ), int(lda, kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else if ( ii > jj ) then
               call SGEMV('n', int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), alpha, a( ii, jj ), &
                  int(lda,kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else
               if (temp == zero) then
                  y(1:n) = zero
               else if (temp /= one) then
                  ! should not happen
                  call SSCAL( int(jc,kind=BLAS_KIND), temp, y( iy ), 1_BLAS_KIND)
               end if
               call SCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call STRMV( 'l', 'n', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call SAXPY( int(jc,kind=BLAS_KIND),alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)

               call SCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call STRMV( 'l', 't', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call SAXPY(int(jc,kind=BLAS_KIND), -alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)
            end if
         end do
         temp = one
      end do

      return
   end subroutine

   subroutine elpa_zssr2(n, x, y,  a, lda )

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none

      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: ck = C_DOUBLE_COMPLEX
      integer, parameter :: rck = C_DOUBLE_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)

      integer(kind=BLAS_KIND)     :: n, lda
      complex(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      complex(kind=rck)     :: temp1, temp2
      integer(kind=ik)            :: i, j, ii, jj, ic, ix, iy, jc, jx, jy, info
      logical                     :: upper

      ! test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 5
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access A in lower triangular part.
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         jy = 1 + (jj-1)

         do j = 1, jc-1
            ! Do local update for blocks on the diagonal
            if ( ( x( jx + j -1) /= zero ) .or. &
               ( y( jy + j -1 ) /= zero ) ) then
               temp1 = - y( jy + j - 1 )
               temp2 = - x( jy + j - 1 )
               do i = j+1, jc
                  a( jj +  i -1 , jj +  j -1 ) = a(jj + i -1,jj +  j -1 ) + x( jx + i  -1)*temp1 - y(jj +  i -1 )*temp2
               end do
            end if
         end do

         ! Use dger for other blocks
         do ii = jj+nb, n, nb
            ic = min( nb, n-ii+1 )
            ix = 1 + (ii-1)
            iy = 1 + (ii-1)
         end do
      end do

      return
   end subroutine
   subroutine elpa_zssmv(n, alpha, a, lda, x,  y)

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none

      integer, parameter :: rk = C_DOUBLE
      integer, parameter :: ck = C_DOUBLE_COMPLEX
      integer, parameter :: rck = C_DOUBLE_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)

      integer(kind=BLAS_KIND)     :: n, lda
      complex(kind=rck)     :: alpha
      complex(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      integer(kind=ik)            :: ii, jj, ic, iy, jc, jx, info
      complex(kind=rck)     :: temp
      complex(kind=rck)     :: work( nb )

      ! Test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 4
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access only lower triangular part of a

      temp = zero
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         do ii = 1, n, nb
            ic = min( nb, n-ii+1 )
            iy = 1 + (ii-1)

            ! gemv for non-diagonal blocks. use 2x dtrmv for diagonal blocks
            if ( ii < jj ) then
               call ZGEMV('t', int(jc,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -alpha, &
                  a( jj, ii ), int(lda, kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else if ( ii > jj ) then
               call ZGEMV('n', int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), alpha, a( ii, jj ), &
                  int(lda,kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else
               if (temp == zero) then
                  y(1:n) = zero
               else if (temp /= one) then
                  ! should not happen
                  call ZSCAL( int(jc,kind=BLAS_KIND), temp, y( iy ), 1_BLAS_KIND)
               end if
               call ZCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call ZTRMV( 'l', 'n', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call ZAXPY( int(jc,kind=BLAS_KIND),alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)

               call ZCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call ZTRMV( 'l', 't', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call ZAXPY(int(jc,kind=BLAS_KIND), -alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)
            end if
         end do
         temp = one
      end do

      return
   end subroutine

   subroutine elpa_cssr2(n, x, y,  a, lda )

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none

      integer, parameter :: rk = C_FLOAT
      integer, parameter :: ck = C_FLOAT_COMPLEX
      integer, parameter :: rck = C_FLOAT_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)

      integer(kind=BLAS_KIND)     :: n, lda
      complex(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      complex(kind=rck)     :: temp1, temp2
      integer(kind=ik)            :: i, j, ii, jj, ic, ix, iy, jc, jx, jy, info
      logical                     :: upper

      ! test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 5
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access A in lower triangular part.
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         jy = 1 + (jj-1)

         do j = 1, jc-1
            ! Do local update for blocks on the diagonal
            if ( ( x( jx + j -1) /= zero ) .or. &
               ( y( jy + j -1 ) /= zero ) ) then
               temp1 = - y( jy + j - 1 )
               temp2 = - x( jy + j - 1 )
               do i = j+1, jc
                  a( jj +  i -1 , jj +  j -1 ) = a(jj + i -1,jj +  j -1 ) + x( jx + i  -1)*temp1 - y(jj +  i -1 )*temp2
               end do
            end if
         end do

         ! Use dger for other blocks
         do ii = jj+nb, n, nb
            ic = min( nb, n-ii+1 )
            ix = 1 + (ii-1)
            iy = 1 + (ii-1)
         end do
      end do

      return
   end subroutine
   subroutine elpa_cssmv(n, alpha, a, lda, x,  y)

      use precision
      use elpa_utilities, only : error_unit
      use elpa_blas_interfaces
      implicit none

      integer, parameter :: rk = C_FLOAT
      integer, parameter :: ck = C_FLOAT_COMPLEX
      integer, parameter :: rck = C_FLOAT_COMPLEX
      complex(kind=rck), parameter     :: ZERO = (0.0_rk,0.0_rk), ONE = (1.0_rk,0.0_rk)

      integer(kind=BLAS_KIND)     :: n, lda
      complex(kind=rck)     :: alpha
      complex(kind=rck)     :: a( lda, * ), x( * ), y( * )
      integer(kind=ik), parameter :: nb = 64
      integer(kind=ik)            :: ii, jj, ic, iy, jc, jx, info
      complex(kind=rck)     :: temp
      complex(kind=rck)     :: work( nb )

      ! Test the input parameters.
      info = 0
      if (n == 0) then
         return
      end if
      if ( n < 0 ) then
         info = 1
      else if ( lda < max( 1,n ) ) then
         info = 4
      end if
      if ( info /= 0 ) then
         write(error_unit,*) "wrong arguments in elpa_ssmv, info =", info
         return
      end if

      ! Access only lower triangular part of a

      temp = zero
      do jj = 1, n, nb
         jc = min( nb, n-jj+1 )
         jx = 1 + (jj-1)
         do ii = 1, n, nb
            ic = min( nb, n-ii+1 )
            iy = 1 + (ii-1)

            ! gemv for non-diagonal blocks. use 2x dtrmv for diagonal blocks
            if ( ii < jj ) then
               call CGEMV('t', int(jc,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), -alpha, &
                  a( jj, ii ), int(lda, kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else if ( ii > jj ) then
               call CGEMV('n', int(ic,kind=BLAS_KIND), int(nb,kind=BLAS_KIND), alpha, a( ii, jj ), &
                  int(lda,kind=BLAS_KIND), &
                  x( jx ), 1_BLAS_KIND, temp, y( iy ), 1_BLAS_KIND )
            else
               if (temp == zero) then
                  y(1:n) = zero
               else if (temp /= one) then
                  ! should not happen
                  call CSCAL( int(jc,kind=BLAS_KIND), temp, y( iy ), 1_BLAS_KIND)
               end if
               call CCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call CTRMV( 'l', 'n', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call CAXPY( int(jc,kind=BLAS_KIND),alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)

               call CCOPY( int(jc,kind=BLAS_KIND), x( jx ), 1_BLAS_KIND, work, 1_BLAS_KIND )
               call CTRMV( 'l', 't', 'n', int(jc,kind=BLAS_KIND), a( jj, jj ), int(lda,kind=BLAS_KIND), work, 1_BLAS_KIND )
               call CAXPY(int(jc,kind=BLAS_KIND), -alpha, work, 1_BLAS_KIND, y( iy ), 1_BLAS_KIND)
            end if
         end do
         temp = one
      end do

      return
   end subroutine

end module
