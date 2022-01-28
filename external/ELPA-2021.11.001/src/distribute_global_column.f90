











module distribute_global_column
  use precision
  implicit none
  private

  public :: distribute_global_column_double
  public :: distribute_global_column_single

  contains

! real double precision first
















subroutine distribute_global_column_&
&double&
&(obj, g_col, l_col, noff, nlen, my_prow, np_rows, nblk)
  use precision
  use elpa_abstract_impl
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
  real(kind=rck), parameter      :: ZERO=0.0_c_double, ONE = 1.0_c_double


  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)             :: noff, nlen, my_prow, np_rows, nblk
  real(kind=rk)     :: g_col(nlen), l_col(*) ! chnage this to proper 2d 1d matching ! remove assumed size

  integer(kind=ik)  :: nbs, nbe, jb, g_off, l_off, js, je

  nbs = noff/(nblk*np_rows)
  nbe = (noff+nlen-1)/(nblk*np_rows)

  do jb = nbs, nbe
    g_off = jb*nblk*np_rows + nblk*my_prow
    l_off = jb*nblk

    js = MAX(noff+1-g_off,1)
    je = MIN(noff+nlen-g_off,nblk)

    if (je<js) cycle

    l_col(l_off+js:l_off+je) = g_col(g_off+js-noff:g_off+je-noff)

  enddo
end subroutine distribute_global_column_&
&double

! real single precision first


















subroutine distribute_global_column_&
&single&
&(obj, g_col, l_col, noff, nlen, my_prow, np_rows, nblk)
  use precision
  use elpa_abstract_impl
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
  real(kind=rck), parameter      :: ZERO=0.0_c_float, ONE = 1.0_c_float


  class(elpa_abstract_impl_t), intent(inout) :: obj
  integer(kind=ik)             :: noff, nlen, my_prow, np_rows, nblk
  real(kind=rk)     :: g_col(nlen), l_col(*) ! chnage this to proper 2d 1d matching ! remove assumed size

  integer(kind=ik)  :: nbs, nbe, jb, g_off, l_off, js, je

  nbs = noff/(nblk*np_rows)
  nbe = (noff+nlen-1)/(nblk*np_rows)

  do jb = nbs, nbe
    g_off = jb*nblk*np_rows + nblk*my_prow
    l_off = jb*nblk

    js = MAX(noff+1-g_off,1)
    je = MIN(noff+nlen-g_off,nblk)

    if (je<js) cycle

    l_col(l_off+js:l_off+je) = g_col(g_off+js-noff:g_off+je-noff)

  enddo
end subroutine distribute_global_column_&
&single

end module
