










!    Copyright 2021, A. Marek
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
!    along with ELPA. If not, see <http://www.gnu.org/licenses/>
!
!    ELPA reflects a substantial effort on the part of the original
!    ELPA consortium, and we ask you to respect the spirit of the
!    license that we chose: i.e., please contribute any changes you
!    may have back to the original ELPA library distribution, and keep
!    any derivatives of ELPA under the same license that we chose for
!    the original distribution, the GNU Lesser General Public License.
!
! Author: Andreas Marek, MPCDF

module thread_affinity
  use precision

  implicit none

  public :: check_thread_affinity, &
            init_thread_affinity, cleanup_thread_affinity, print_thread_affinity
  private
! integer(kind=ik) :: thread_num
  integer(kind=ik) :: thread_max
  integer(kind=ik) :: process_cpu_id
  integer(kind=ik), allocatable :: cpu_ids(:)

  interface
    subroutine get_process_id_c(process_id, pprocess_id) bind(C, name="get_process_id")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT), intent(out) :: process_id, pprocess_id
    end subroutine
  end interface

  interface
    subroutine get_thread_affinity_c(cpu_id) bind(C, name="get_thread_affinity")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT) :: cpu_id
    end subroutine
  end interface
  interface
    subroutine get_process_affinity_c(cpu_id) bind(C, name="get_process_affinity")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT), value :: cpu_id
    end subroutine
  end interface

contains
  subroutine get_thread_affinity(cpu_id)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=ik), intent(out) :: cpu_id
    integer(kind=C_INT) :: cpu_id_c
    call get_thread_affinity_c(cpu_id_c)
    cpu_id = int(cpu_id_c, kind=ik)
  end subroutine
  subroutine get_process_affinity(cpu_id)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=ik), intent(out) :: cpu_id
    integer(kind=C_INT) :: cpu_id_c
    call get_process_affinity_c(cpu_id_c)
    cpu_id = int(cpu_id_c, kind=ik)
  end subroutine
  subroutine get_process_id(process_id, pprocess_id)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=ik), intent(out) :: process_id, pprocess_id
    integer(kind=C_INT) :: process_id_c, pprocess_id_c
    call get_process_id_c(process_id_c, pprocess_id_c)
    process_id  = int(process_id_c,  kind=ik)
    pprocess_id = int(pprocess_id_c, kind=ik)
  end subroutine


  subroutine init_thread_affinity(nrThreads)
    use precision
    use omp_lib

    implicit none
    integer(kind=ik)             :: istat
    integer(kind=ik), intent(in) :: nrThreads

    thread_max = nrThreads
  end subroutine init_thread_affinity

  subroutine cleanup_thread_affinity
    use precision
    implicit none
    integer(kind=ik) :: istat

    if((allocated(cpu_ids))) then
       deallocate(cpu_ids, stat=istat)
       if (istat .ne. 0) then
         print *,"Error when deallocating init_thread_affinity"
       endif
    endif

  end subroutine cleanup_thread_affinity

  subroutine check_thread_affinity()
    use precision
    use omp_lib
    implicit none
    integer(kind=ik)             :: thread_cpu_id
    integer(kind=ik)             :: i, actuall_num

    call get_process_affinity(process_cpu_id)


  end subroutine check_thread_affinity

  subroutine print_thread_affinity(mype)

    use precision
    implicit none
    integer(kind=ik) :: i
    integer(kind=ik), intent(in) :: mype
    integer(kind=ik) :: pid, ppid

    call get_process_id(pid, ppid)
    write(*,'("Task ",i4," runs on process id: ",i4," with pid ",i4," and ppid ",i4)') mype, process_cpu_id,pid,ppid
  end subroutine print_thread_affinity

end module thread_affinity
