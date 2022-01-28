










!    Copyright 2021, A. Marek
!
!    This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
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
! This file was written by A. Marek, MPCDF


module cholesky_cuda
  use, intrinsic :: iso_c_binding
  use precision
  implicit none

  public
  interface
    subroutine copy_double_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)&
             bind(C, name="copy_double_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
    end subroutine 
  end interface

  interface
    subroutine copy_float_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)&
             bind(C, name="copy_float_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
    end subroutine 
  end interface

  interface
    subroutine copy_double_complex_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)&
             bind(C, name="copy_double_complex_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
    end subroutine 
  end interface

  interface
    subroutine copy_float_complex_a_tmatc_c(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)&
             bind(C, name="copy_float_complex_a_tmatc_FromC")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_intptr_T), value  :: a_dev, tmatc_dev
      integer(kind=c_int), intent(in)  :: nblk, matrixRows, l_cols, l_colx, l_row1
    end subroutine 
  end interface

  contains

    subroutine copy_double_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev

    end subroutine

    subroutine copy_float_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev

    end subroutine

    subroutine copy_double_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev

    end subroutine

    subroutine copy_float_complex_a_tmatc(a_dev, tmatc_dev, nblk, matrixRows, l_cols, l_colx, l_row1)
      use, intrinsic :: iso_c_binding

      implicit none
      integer(kind=C_INT), intent(in) :: nblk, matrixRows, l_cols, l_colx, l_row1
      integer(kind=C_intptr_T)        :: a_dev, tmatc_dev

    end subroutine

end module

