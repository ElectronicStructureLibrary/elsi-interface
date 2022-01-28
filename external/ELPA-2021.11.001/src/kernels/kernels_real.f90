










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
! --------------------------------------------------------------------------------------------------
!
! This file contains the compute intensive kernels for the Householder transformations.
! It should be compiled with the highest possible optimization level.
!
! On Intel use -O3 -xSSE4.2 (or the SSE level fitting to your CPU)
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
! --------------------------------------------------------------------------------------------------



















! the intel compiler creates a temp copy of array q!
! this should be prevented if possible without using assumed size arrays
  subroutine double_hh_trafo_&
  &real&
  &_generic_&
  &double&
  & (q, hh, nb, nq, ldq, ldh)

    use precision
    use, intrinsic :: iso_c_binding
    use elpa_abstract_impl
    implicit none

    !class(elpa_abstract_impl_t), intent(inout) :: obj
    integer(kind=ik), intent(in)            :: nb, nq, ldq, ldh
    real(kind=c_double), intent(inout) :: q(ldq,*)
    real(kind=c_double), intent(in)    :: hh(ldh,*)

    real(kind=c_double)                :: s
    integer(kind=ik)                        :: i

!    equivalence(q(1,1),q_complex(1,1))

    ! Safety only:

!    call obj%timer%start("kernel generic: double_hh_trafo_&
!                     &real&
!		     &_generic" // &
!		     &"_double" &
!		     )

    if(mod(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'

    ! Calculate dot product of the two Householder vectors

    s = hh(2,2)*1.0
    do i=3,nb
       s = s+hh(i,2)*hh(i-1,1)
    enddo

    ! Do the Householder transformations

    ! Always a multiple of 4 Q-rows is transformed, even if nq is smaller

    do i=1,nq-8,12

       call hh_trafo_kernel_12_generic_&
       &double&
       & (q(i,1),hh, nb, ldq, ldh, s)

    enddo

    ! i > nq-8 now, i.e. at most 8 rows remain

    if(nq-i+1 > 4) then

       call hh_trafo_kernel_8_generic_&
       &double&
       & (q(i,1),hh, nb, ldq, ldh, s)

    else if(nq-i+1 > 0) then

       call hh_trafo_kernel_4_generic_&
       &double&
       & (q(i,1),hh, nb, ldq, ldh, s)

    endif

!    call obj%timer%stop("kernel generic: double_hh_trafo_&
!                     &real&
!		     &_generic" // &
!		     &"_double" &
!		     )

  end subroutine

  ! --------------------------------------------------------------------------------------------------
  ! The following kernels perform the Householder transformation on Q for 12/8/4 rows.
  ! Please note that Q is declared complex*16 here.
  ! This is a hint for compilers that packed arithmetic can be used for Q
  ! (relevant for Intel SSE and BlueGene double hummer CPUs).
  ! --------------------------------------------------------------------------------------------------

  subroutine hh_trafo_kernel_12_generic_&
  &double&
  & (q, hh, nb, ldq, ldh, s)
    use precision
    implicit none
    integer(kind=ik), intent(in)    :: nb, ldq, ldh
    complex(kind=ck8), intent(inout) :: q(ldq/2,*)
    real(kind=c_double), intent(in)    :: hh(ldh,*)
    real(kind=c_double), intent(in)    :: s

    complex(kind=ck8)    :: x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6
    real(kind=c_double)                   :: h1, h2, tau1, tau2
    integer(kind=ik)                :: i


    !call obj%timer%start("kernel generic: hh_trafo_kernel_12_generic" // &
    !		     &"_double" &
    !		     )

    x1  = q(1,2)
    x2  = q(2,2)
    x3  = q(3,2)
    x4  = q(4,2)
    x5  = q(5,2)
    x6  = q(6,2)

    y1  = q(1 ,1) + q(1, 2)*hh(2,2)
    y2  = q(2 ,1) + q(2, 2)*hh(2,2)
    y3  = q(3 ,1) + q(3, 2)*hh(2,2)
    y4  = q(4 ,1) + q(4, 2)*hh(2,2)
    y5  = q(5 ,1) + q(5, 2)*hh(2,2)
    y6  = q(6 ,1) + q(6, 2)*hh(2,2)


    do i=3,nb
       h1  = hh(i-1,1)
       h2  = hh(i,2)
       x1  = x1 + q(1, i)*h1
       y1  = y1 + q(1, i)*h2
       x2  = x2 + q(2, i)*h1
       y2  = y2 + q(2, i)*h2
       x3  = x3 + q(3, i)*h1
       y3  = y3 + q(3, i)*h2
       x4  = x4 + q(4, i)*h1
       y4  = y4 + q(4, i)*h2
       x5  = x5 + q(5, i)*h1
       y5  = y5 + q(5, i)*h2
       x6  = x6 + q(6, i)*h1
       y6  = y6 + q(6, i)*h2
    enddo

    x1  = x1  + q(1,nb+1)*hh(nb,1)
    x2  = x2  + q(2,nb+1)*hh(nb,1)
    x3  = x3  + q(3,nb+1)*hh(nb,1)
    x4  = x4  + q(4,nb+1)*hh(nb,1)
    x5  = x5  + q(5,nb+1)*hh(nb,1)
    x6  = x6  + q(6,nb+1)*hh(nb,1)

    tau1 = hh(1,1)
    tau2 = hh(1,2)

    h1  = -tau1
    x1  = x1 *h1
    x2  = x2 *h1
    x3  = x3 *h1
    x4  = x4 *h1
    x5  = x5 *h1
    x6  = x6 *h1

    h1  = -tau2
    h2  = -tau2*s
    y1  = y1 *h1 + x1 *h2
    y2  = y2 *h1 + x2 *h2
    y3  = y3 *h1 + x3 *h2
    y4  = y4 *h1 + x4 *h2
    y5  = y5 *h1 + x5 *h2
    y6  = y6 *h1 + x6 *h2
    q(1,1)  = q(1, 1) + y1
    q(2,1)  = q(2, 1) + y2
    q(3,1)  = q(3, 1) + y3
    q(4,1)  = q(4, 1) + y4
    q(5,1)  = q(5, 1) + y5
    q(6,1)  = q(6, 1) + y6

    q(1, 2) = q(1, 2) + x1  + y1 *hh(2,2)
    q(2, 2) = q(2, 2) + x2  + y2 *hh(2,2)
    q(3, 2) = q(3, 2) + x3  + y3 *hh(2,2)
    q(4, 2) = q(4, 2) + x4  + y4 *hh(2,2)
    q(5, 2) = q(5, 2) + x5  + y5 *hh(2,2)
    q(6, 2) = q(6, 2) + x6  + y6 *hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       q(1, i) = q(1,i)  + x1 *h1 + y1 *h2
       q(2, i) = q(2,i)  + x2 *h1 + y2 *h2
       q(3, i) = q(3,i)  + x3 *h1 + y3 *h2
       q(4, i) = q(4,i)  + x4 *h1 + y4 *h2
       q(5, i) = q(5,i)  + x5 *h1 + y5 *h2
       q(6, i) = q(6,i)  + x6 *h1 + y6 *h2
    enddo

    q(1, nb+1) = q(1, nb+1) + x1 *hh(nb,1)
    q(2, nb+1) = q(2, nb+1) + x2 *hh(nb,1)
    q(3, nb+1) = q(3, nb+1) + x3 *hh(nb,1)
    q(4, nb+1) = q(4, nb+1) + x4 *hh(nb,1)
    q(5, nb+1) = q(5, nb+1) + x5 *hh(nb,1)
    q(6, nb+1) = q(6, nb+1) + x6 *hh(nb,1)

!    call obj%timer%stop("kernel generic: hh_trafo_kernel_12_generic" // &
!		     &"_double" &
!		     )

  end subroutine
  ! --------------------------------------------------------------------------------------------------

  subroutine hh_trafo_kernel_8_generic_&
  &double&
  & (q, hh, nb, ldq, ldh, s)
    use precision
    implicit none
    integer(kind=ik), intent(in)     :: nb, ldq, ldh
    complex(kind=ck8), intent(inout)  :: q(ldq/2,*)
    real(kind=c_double), intent(in)     :: hh(ldh,*)
    real(kind=c_double), intent(in)     :: s
    complex(kind=ck8)     :: x1, x2, x3, x4, y1, y2, y3, y4
    real(kind=c_double)                    :: h1, h2, tau1, tau2
    integer(kind=ik)                 :: i

!    call obj%timer%start("kernel generic: hh_trafo_kernel_8_generic" // &
!		     &"_double" &
!		     )
    x1 = q(1,2)
    x2 = q(2,2)
    x3 = q(3,2)
    x4 = q(4,2)

    y1 = q(1,1) + q(1,2)*hh(2,2)
    y2 = q(2,1) + q(2,2)*hh(2,2)
    y3 = q(3,1) + q(3,2)*hh(2,2)
    y4 = q(4,1) + q(4,2)*hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       x1 = x1 + q(1,i)*h1
       y1 = y1 + q(1,i)*h2
       x2 = x2 + q(2,i)*h1
       y2 = y2 + q(2,i)*h2
       x3 = x3 + q(3,i)*h1
       y3 = y3 + q(3,i)*h2
       x4 = x4 + q(4,i)*h1
       y4 = y4 + q(4,i)*h2
    enddo

    x1 = x1 + q(1,nb+1)*hh(nb,1)
    x2 = x2 + q(2,nb+1)*hh(nb,1)
    x3 = x3 + q(3,nb+1)*hh(nb,1)
    x4 = x4 + q(4,nb+1)*hh(nb,1)

    tau1 = hh(1,1)
    tau2 = hh(1,2)

    h1 = -tau1
    x1 = x1*h1
    x2 = x2*h1
    x3 = x3*h1
    x4 = x4*h1
    h1 = -tau2
    h2 = -tau2*s
    y1 = y1*h1 + x1*h2
    y2 = y2*h1 + x2*h2
    y3 = y3*h1 + x3*h2
    y4 = y4*h1 + x4*h2
    q(1,1) = q(1,1) + y1
    q(2,1) = q(2,1) + y2
    q(3,1) = q(3,1) + y3
    q(4,1) = q(4,1) + y4
    q(1,2) = q(1,2) + x1 + y1*hh(2,2)
    q(2,2) = q(2,2) + x2 + y2*hh(2,2)
    q(3,2) = q(3,2) + x3 + y3*hh(2,2)
    q(4,2) = q(4,2) + x4 + y4*hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       q(1,i) = q(1,i) + x1*h1 + y1*h2
       q(2,i) = q(2,i) + x2*h1 + y2*h2
       q(3,i) = q(3,i) + x3*h1 + y3*h2
       q(4,i) = q(4,i) + x4*h1 + y4*h2
    enddo

    q(1,nb+1) = q(1,nb+1) + x1*hh(nb,1)
    q(2,nb+1) = q(2,nb+1) + x2*hh(nb,1)
    q(3,nb+1) = q(3,nb+1) + x3*hh(nb,1)
    q(4,nb+1) = q(4,nb+1) + x4*hh(nb,1)


!    call obj%timer%stop("kernel generic: hh_trafo_kernel_8_generic" // &
!		     &"_double" &
!		     )

  end subroutine
  ! --------------------------------------------------------------------------------------------------

  subroutine hh_trafo_kernel_4_generic_&
  &double&
  & (q, hh, nb, ldq, ldh, s)

    use precision
    implicit none
    integer(kind=ik), intent(in)    :: nb, ldq, ldh
    complex(kind=ck8), intent(inout) :: q(ldq/2,*)
    real(kind=c_double), intent(in)    :: hh(ldh,*)
    real(kind=c_double), intent(in)    :: s

    complex(kind=ck8)    :: x1, x2, y1, y2
    real(kind=c_double)                :: h1, h2, tau1, tau2
    integer(kind=ik)                :: i

!    call obj%timer%start("kernel generic: hh_trafo_kernel_4_generic" // &
!		     &"_double" &
!		     )
    x1 = q(1,2)
    x2 = q(2,2)

    y1 = q(1,1) + q(1,2)*hh(2,2)
    y2 = q(2,1) + q(2,2)*hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       x1 = x1 + q(1,i)*h1
       y1 = y1 + q(1,i)*h2
       x2 = x2 + q(2,i)*h1
       y2 = y2 + q(2,i)*h2
    enddo

    x1 = x1 + q(1,nb+1)*hh(nb,1)
    x2 = x2 + q(2,nb+1)*hh(nb,1)

    tau1 = hh(1,1)
    tau2 = hh(1,2)

    h1 = -tau1
    x1 = x1*h1
    x2 = x2*h1
    h1 = -tau2
    h2 = -tau2*s
    y1 = y1*h1 + x1*h2
    y2 = y2*h1 + x2*h2

    q(1,1) = q(1,1) + y1
    q(2,1) = q(2,1) + y2
    q(1,2) = q(1,2) + x1 + y1*hh(2,2)
    q(2,2) = q(2,2) + x2 + y2*hh(2,2)

    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       q(1,i) = q(1,i) + x1*h1 + y1*h2
       q(2,i) = q(2,i) + x2*h1 + y2*h2
    enddo

    q(1,nb+1) = q(1,nb+1) + x1*hh(nb,1)
    q(2,nb+1) = q(2,nb+1) + x2*hh(nb,1)

!    call obj%timer%stop("kernel generic: hh_trafo_kernel_4_generic" // &
!		     &"_double" &
!		     )

  end subroutine




















! the intel compiler creates a temp copy of array q!
! this should be prevented if possible without using assumed size arrays
  subroutine double_hh_trafo_&
  &real&
  &_generic_&
  &single&
  & (q, hh, nb, nq, ldq, ldh)

    use precision
    use, intrinsic :: iso_c_binding
    use elpa_abstract_impl
    implicit none

    !class(elpa_abstract_impl_t), intent(inout) :: obj
    integer(kind=ik), intent(in)            :: nb, nq, ldq, ldh
    real(kind=c_float), intent(inout) :: q(ldq,*)
    real(kind=c_float), intent(in)    :: hh(ldh,*)

    real(kind=c_float)                :: s
    integer(kind=ik)                        :: i

!    equivalence(q(1,1),q_complex(1,1))

    ! Safety only:

!    call obj%timer%start("kernel generic: double_hh_trafo_&
!                     &real&
!		     &_generic" // &
!		     &"_single" &
!		     )

    if(mod(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'

    ! Calculate dot product of the two Householder vectors

    s = hh(2,2)*1.0
    do i=3,nb
       s = s+hh(i,2)*hh(i-1,1)
    enddo

    ! Do the Householder transformations

    ! Always a multiple of 4 Q-rows is transformed, even if nq is smaller

    do i=1,nq-8,12

       call hh_trafo_kernel_12_generic_&
       &single&
       & (q(i,1),hh, nb, ldq, ldh, s)

    enddo

    ! i > nq-8 now, i.e. at most 8 rows remain

    if(nq-i+1 > 4) then

       call hh_trafo_kernel_8_generic_&
       &single&
       & (q(i,1),hh, nb, ldq, ldh, s)

    else if(nq-i+1 > 0) then

       call hh_trafo_kernel_4_generic_&
       &single&
       & (q(i,1),hh, nb, ldq, ldh, s)

    endif

!    call obj%timer%stop("kernel generic: double_hh_trafo_&
!                     &real&
!		     &_generic" // &
!		     &"_single" &
!		     )

  end subroutine

  ! --------------------------------------------------------------------------------------------------
  ! The following kernels perform the Householder transformation on Q for 12/8/4 rows.
  ! Please note that Q is declared complex*16 here.
  ! This is a hint for compilers that packed arithmetic can be used for Q
  ! (relevant for Intel SSE and BlueGene double hummer CPUs).
  ! --------------------------------------------------------------------------------------------------

  subroutine hh_trafo_kernel_12_generic_&
  &single&
  & (q, hh, nb, ldq, ldh, s)
    use precision
    implicit none
    integer(kind=ik), intent(in)    :: nb, ldq, ldh
    complex(kind=ck4), intent(inout) :: q(ldq/2,*)
    real(kind=c_float), intent(in)    :: hh(ldh,*)
    real(kind=c_float), intent(in)    :: s

    complex(kind=ck4)    :: x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6
    real(kind=c_float)                   :: h1, h2, tau1, tau2
    integer(kind=ik)                :: i


    !call obj%timer%start("kernel generic: hh_trafo_kernel_12_generic" // &
    !		     &"_single" &
    !		     )

    x1  = q(1,2)
    x2  = q(2,2)
    x3  = q(3,2)
    x4  = q(4,2)
    x5  = q(5,2)
    x6  = q(6,2)

    y1  = q(1 ,1) + q(1, 2)*hh(2,2)
    y2  = q(2 ,1) + q(2, 2)*hh(2,2)
    y3  = q(3 ,1) + q(3, 2)*hh(2,2)
    y4  = q(4 ,1) + q(4, 2)*hh(2,2)
    y5  = q(5 ,1) + q(5, 2)*hh(2,2)
    y6  = q(6 ,1) + q(6, 2)*hh(2,2)


    do i=3,nb
       h1  = hh(i-1,1)
       h2  = hh(i,2)
       x1  = x1 + q(1, i)*h1
       y1  = y1 + q(1, i)*h2
       x2  = x2 + q(2, i)*h1
       y2  = y2 + q(2, i)*h2
       x3  = x3 + q(3, i)*h1
       y3  = y3 + q(3, i)*h2
       x4  = x4 + q(4, i)*h1
       y4  = y4 + q(4, i)*h2
       x5  = x5 + q(5, i)*h1
       y5  = y5 + q(5, i)*h2
       x6  = x6 + q(6, i)*h1
       y6  = y6 + q(6, i)*h2
    enddo

    x1  = x1  + q(1,nb+1)*hh(nb,1)
    x2  = x2  + q(2,nb+1)*hh(nb,1)
    x3  = x3  + q(3,nb+1)*hh(nb,1)
    x4  = x4  + q(4,nb+1)*hh(nb,1)
    x5  = x5  + q(5,nb+1)*hh(nb,1)
    x6  = x6  + q(6,nb+1)*hh(nb,1)

    tau1 = hh(1,1)
    tau2 = hh(1,2)

    h1  = -tau1
    x1  = x1 *h1
    x2  = x2 *h1
    x3  = x3 *h1
    x4  = x4 *h1
    x5  = x5 *h1
    x6  = x6 *h1

    h1  = -tau2
    h2  = -tau2*s
    y1  = y1 *h1 + x1 *h2
    y2  = y2 *h1 + x2 *h2
    y3  = y3 *h1 + x3 *h2
    y4  = y4 *h1 + x4 *h2
    y5  = y5 *h1 + x5 *h2
    y6  = y6 *h1 + x6 *h2
    q(1,1)  = q(1, 1) + y1
    q(2,1)  = q(2, 1) + y2
    q(3,1)  = q(3, 1) + y3
    q(4,1)  = q(4, 1) + y4
    q(5,1)  = q(5, 1) + y5
    q(6,1)  = q(6, 1) + y6

    q(1, 2) = q(1, 2) + x1  + y1 *hh(2,2)
    q(2, 2) = q(2, 2) + x2  + y2 *hh(2,2)
    q(3, 2) = q(3, 2) + x3  + y3 *hh(2,2)
    q(4, 2) = q(4, 2) + x4  + y4 *hh(2,2)
    q(5, 2) = q(5, 2) + x5  + y5 *hh(2,2)
    q(6, 2) = q(6, 2) + x6  + y6 *hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       q(1, i) = q(1,i)  + x1 *h1 + y1 *h2
       q(2, i) = q(2,i)  + x2 *h1 + y2 *h2
       q(3, i) = q(3,i)  + x3 *h1 + y3 *h2
       q(4, i) = q(4,i)  + x4 *h1 + y4 *h2
       q(5, i) = q(5,i)  + x5 *h1 + y5 *h2
       q(6, i) = q(6,i)  + x6 *h1 + y6 *h2
    enddo

    q(1, nb+1) = q(1, nb+1) + x1 *hh(nb,1)
    q(2, nb+1) = q(2, nb+1) + x2 *hh(nb,1)
    q(3, nb+1) = q(3, nb+1) + x3 *hh(nb,1)
    q(4, nb+1) = q(4, nb+1) + x4 *hh(nb,1)
    q(5, nb+1) = q(5, nb+1) + x5 *hh(nb,1)
    q(6, nb+1) = q(6, nb+1) + x6 *hh(nb,1)

!    call obj%timer%stop("kernel generic: hh_trafo_kernel_12_generic" // &
!		     &"_single" &
!		     )

  end subroutine
  ! --------------------------------------------------------------------------------------------------

  subroutine hh_trafo_kernel_8_generic_&
  &single&
  & (q, hh, nb, ldq, ldh, s)
    use precision
    implicit none
    integer(kind=ik), intent(in)     :: nb, ldq, ldh
    complex(kind=ck4), intent(inout)  :: q(ldq/2,*)
    real(kind=c_float), intent(in)     :: hh(ldh,*)
    real(kind=c_float), intent(in)     :: s
    complex(kind=ck4)     :: x1, x2, x3, x4, y1, y2, y3, y4
    real(kind=c_float)                    :: h1, h2, tau1, tau2
    integer(kind=ik)                 :: i

!    call obj%timer%start("kernel generic: hh_trafo_kernel_8_generic" // &
!		     &"_single" &
!		     )
    x1 = q(1,2)
    x2 = q(2,2)
    x3 = q(3,2)
    x4 = q(4,2)

    y1 = q(1,1) + q(1,2)*hh(2,2)
    y2 = q(2,1) + q(2,2)*hh(2,2)
    y3 = q(3,1) + q(3,2)*hh(2,2)
    y4 = q(4,1) + q(4,2)*hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       x1 = x1 + q(1,i)*h1
       y1 = y1 + q(1,i)*h2
       x2 = x2 + q(2,i)*h1
       y2 = y2 + q(2,i)*h2
       x3 = x3 + q(3,i)*h1
       y3 = y3 + q(3,i)*h2
       x4 = x4 + q(4,i)*h1
       y4 = y4 + q(4,i)*h2
    enddo

    x1 = x1 + q(1,nb+1)*hh(nb,1)
    x2 = x2 + q(2,nb+1)*hh(nb,1)
    x3 = x3 + q(3,nb+1)*hh(nb,1)
    x4 = x4 + q(4,nb+1)*hh(nb,1)

    tau1 = hh(1,1)
    tau2 = hh(1,2)

    h1 = -tau1
    x1 = x1*h1
    x2 = x2*h1
    x3 = x3*h1
    x4 = x4*h1
    h1 = -tau2
    h2 = -tau2*s
    y1 = y1*h1 + x1*h2
    y2 = y2*h1 + x2*h2
    y3 = y3*h1 + x3*h2
    y4 = y4*h1 + x4*h2
    q(1,1) = q(1,1) + y1
    q(2,1) = q(2,1) + y2
    q(3,1) = q(3,1) + y3
    q(4,1) = q(4,1) + y4
    q(1,2) = q(1,2) + x1 + y1*hh(2,2)
    q(2,2) = q(2,2) + x2 + y2*hh(2,2)
    q(3,2) = q(3,2) + x3 + y3*hh(2,2)
    q(4,2) = q(4,2) + x4 + y4*hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       q(1,i) = q(1,i) + x1*h1 + y1*h2
       q(2,i) = q(2,i) + x2*h1 + y2*h2
       q(3,i) = q(3,i) + x3*h1 + y3*h2
       q(4,i) = q(4,i) + x4*h1 + y4*h2
    enddo

    q(1,nb+1) = q(1,nb+1) + x1*hh(nb,1)
    q(2,nb+1) = q(2,nb+1) + x2*hh(nb,1)
    q(3,nb+1) = q(3,nb+1) + x3*hh(nb,1)
    q(4,nb+1) = q(4,nb+1) + x4*hh(nb,1)


!    call obj%timer%stop("kernel generic: hh_trafo_kernel_8_generic" // &
!		     &"_single" &
!		     )

  end subroutine
  ! --------------------------------------------------------------------------------------------------

  subroutine hh_trafo_kernel_4_generic_&
  &single&
  & (q, hh, nb, ldq, ldh, s)

    use precision
    implicit none
    integer(kind=ik), intent(in)    :: nb, ldq, ldh
    complex(kind=ck4), intent(inout) :: q(ldq/2,*)
    real(kind=c_float), intent(in)    :: hh(ldh,*)
    real(kind=c_float), intent(in)    :: s

    complex(kind=ck4)    :: x1, x2, y1, y2
    real(kind=c_float)                :: h1, h2, tau1, tau2
    integer(kind=ik)                :: i

!    call obj%timer%start("kernel generic: hh_trafo_kernel_4_generic" // &
!		     &"_single" &
!		     )
    x1 = q(1,2)
    x2 = q(2,2)

    y1 = q(1,1) + q(1,2)*hh(2,2)
    y2 = q(2,1) + q(2,2)*hh(2,2)


    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       x1 = x1 + q(1,i)*h1
       y1 = y1 + q(1,i)*h2
       x2 = x2 + q(2,i)*h1
       y2 = y2 + q(2,i)*h2
    enddo

    x1 = x1 + q(1,nb+1)*hh(nb,1)
    x2 = x2 + q(2,nb+1)*hh(nb,1)

    tau1 = hh(1,1)
    tau2 = hh(1,2)

    h1 = -tau1
    x1 = x1*h1
    x2 = x2*h1
    h1 = -tau2
    h2 = -tau2*s
    y1 = y1*h1 + x1*h2
    y2 = y2*h1 + x2*h2

    q(1,1) = q(1,1) + y1
    q(2,1) = q(2,1) + y2
    q(1,2) = q(1,2) + x1 + y1*hh(2,2)
    q(2,2) = q(2,2) + x2 + y2*hh(2,2)

    do i=3,nb
       h1 = hh(i-1,1)
       h2 = hh(i,2)
       q(1,i) = q(1,i) + x1*h1 + y1*h2
       q(2,i) = q(2,i) + x2*h1 + y2*h2
    enddo

    q(1,nb+1) = q(1,nb+1) + x1*hh(nb,1)
    q(2,nb+1) = q(2,nb+1) + x2*hh(nb,1)

!    call obj%timer%stop("kernel generic: hh_trafo_kernel_4_generic" // &
!		     &"_single" &
!		     )

  end subroutine

! --------------------------------------------------------------------------------------------------
