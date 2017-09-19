!   This file is part of ELPA.
!
!    The ELPA library was originally created by the ELPA consortium,
!    consisting of the following organizations:
!
!    - Max Planck Computing and Data Facility (MPCDF), fomerly known as
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



! ELPA2 -- 2-stage solver for ELPA
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".


!> \brief Fortran module which provides the routines to use the 2-stage ELPA solver
module CHECK_SINGULARITY

! Version 1.1.2, 2011-02-21

  use elpa_utilities_cs
  use elpa1_cs, only : elpa_print_times, time_evp_back, time_evp_fwd, time_evp_solve
  use elpa2_utilities_cs
  use elpa2_cs

  implicit none

  PRIVATE ! By default, all routines contained are private

  ! The following routines are public:

  public :: elpa_check_singularity_real_double
  public :: elpa_check_singularity_complex_double

contains

function elpa_check_singularity_real_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
            mpi_comm_rows,mpi_comm_cols,mpi_comm_all,singular_tol,n_nonsing)&
            result(success)

   use elpa_utilities_cs
   use elpa1_cs
   use elpa2_utilities_cs
   use elpa1_compute_cs
   use elpa2_compute_cs
   use elpa_mpi_cs
   use cuda_functions_cs
   use mod_check_for_gpu_cs
   use iso_c_binding

   implicit none

   integer(kind=c_int), intent(in)    :: na,nev,lda,ldq,matrixCols,nblk
   integer(kind=c_int), intent(in)    :: mpi_comm_rows,mpi_comm_cols,mpi_comm_all
   real(kind=c_double), intent(inout) :: ev(na)
   real(kind=c_double), intent(inout) :: a(lda,*),q(ldq,*)
   real(kind=c_double), intent(in)    :: singular_tol
   integer(kind=c_int), intent(out)   :: n_nonsing

   real(kind=c_double), allocatable   :: hh_trans_real(:,:)
   integer(kind=c_int)                :: THIS_REAL_ELPA_KERNEL
   integer(kind=c_int)                :: my_pe,n_pes,my_prow,my_pcol,np_rows,np_cols,mpierr
   integer(kind=c_int)                :: nbw, num_blocks
   real(kind=c_double), allocatable   :: tmat(:,:,:), e(:)
   integer(kind=c_intptr_t)           :: tmat_dev, q_dev, a_dev
   real(kind=c_double)                :: ttt0, ttt1
   integer(kind=c_int)                :: i
   logical                            :: success
   logical                            :: wantDebug
   integer(kind=c_int)                :: istat
   logical                            :: do_useGPU,do_useGPU_4
   logical                            :: useQRActual
   integer(kind=c_int)                :: numberOfGPUDevices

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)
   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   success     = .true.
   wantDebug   = .false.
   useQRActual = .false.
   do_useGPU   = .false.
   do_useGPU_4 = .false.

   THIS_REAL_ELPA_KERNEL = get_actual_real_kernel()

   if(THIS_REAL_ELPA_KERNEL .eq. REAL_ELPA_KERNEL_GPU) then
      if(nblk .ne. 128) then
         print *, " ERROR: ELPA2 GPU kernel requires blocksize = 128. Exiting.."
         stop
      endif

      ! 4th step
      do_useGPU_4 = .true.
   endif

#ifdef WITH_GPU_VERSION
   if(check_for_gpu(my_pe,numberOfGPUDevices,wantDebug=wantDebug)) then
      do_useGPU = .true.
   endif

   ! Set GPU/CUDA parameters
   cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
   cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
   cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
   cudaHostRegisterPortable = cuda_hostRegisterPortable()
   cudaHostRegisterMapped   = cuda_hostRegisterMapped()
#endif

   ! Bandwidth must be a multiple of nblk
   ! Set to a value >= 32
   if(do_useGPU) then
      nbw = nblk
   else
      nbw = (63/nblk+1)*nblk
   endif
   num_blocks = (na-1)/nbw+1

   allocate(tmat(nbw,nbw,num_blocks),stat=istat)
   if(istat .ne. 0) then
      print *," ERROR: Cannot allocate tmat. Exiting.."
      stop
   endif

   ! Reduction full -> band
   ttt0 = MPI_Wtime()
   call bandred_real_double(na,a,a_dev,lda,nblk,nbw,matrixCols,num_blocks,&
                            mpi_comm_rows,mpi_comm_cols,tmat,tmat_dev,wantDebug,&
                            do_useGPU, success, useQRActual)
   if(.not.(success)) return
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      write(*,"(A,F10.3,A)") "  | Time full => band           :",ttt1-ttt0," s"

   ! Reduction band -> tridiagonal
   allocate(e(na),stat=istat)
   if(istat .ne. 0) then
      print *," ERROR: Cannot allocate e. Exiting.."
      stop
   endif

   ttt0 = MPI_Wtime()
   call tridiag_band_real_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_real,&
                                 mpi_comm_rows,mpi_comm_cols,mpi_comm_all)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      write(*,"(A,F10.3,A)") "  | Time band => tridiagonal    :",ttt1-ttt0," s"

   call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ! Solve tridiagonal system
   ttt0 = MPI_Wtime()
   call solve_tridi_double(na,nev,ev,e,q,ldq,nblk,matrixCols,mpi_comm_rows,&
                           mpi_comm_cols,wantDebug,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      write(*,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"

   deallocate(e)

   ! Get the number of nonsingular eigenvalues
   n_nonsing = na
   do i = 1,na
      if(ev(i) < singular_tol) then
         n_nonsing = n_nonsing-1
      endif
   enddo

   if(n_nonsing < na) then
      ! Back-transform 1
      ttt0 = MPI_Wtime()
      call trans_ev_tridi_to_band_real_double(na,nev,nblk,nbw,q,q_dev,ldq,matrixCols,&
                                              hh_trans_real,mpi_comm_rows,mpi_comm_cols,&
                                              wantDebug,do_useGPU_4,success,THIS_REAL_ELPA_KERNEL)
      if(.not.success) return
      ttt1 = MPI_Wtime()
      if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
         write(*,"(A,F10.3,A)") "  | Time ev tridiagonal => band :",ttt1-ttt0," s"

      deallocate(hh_trans_real)

      ! Back-transform 2
      ttt0 = MPI_Wtime()
      call trans_ev_band_to_full_real_double(na,nev,nblk,nbw,a,a_dev,lda,tmat,tmat_dev,q,&
                                             q_dev,ldq,matrixCols,num_blocks,mpi_comm_rows,&
                                             mpi_comm_cols,do_useGPU,useQRActual)
      ttt1 = MPI_Wtime()
      if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
         write(*,"(A,F10.3,A)") "  | Time ev band => full        :",ttt1-ttt0," s"

      deallocate(tmat)
   else
      deallocate(hh_trans_real)
      deallocate(tmat)
   endif

   if(my_prow==0 .and. my_pcol==0) then
      if(do_useGPU) then
         if(do_useGPU_4) then
            print *,"  GPU has been used for the 1st, 4th, 5th steps of ELPA2"
         else
            print *,"  GPU has been used for the 1st, 5th steps of ELPA2"
         endif
      endif
   endif

end function

function elpa_check_singularity_complex_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
            mpi_comm_rows,mpi_comm_cols,mpi_comm_all,singular_tol,n_nonsing)&
            result(success)

   use elpa_utilities_cs
   use elpa1_cs
   use elpa2_utilities_cs
   use elpa1_compute_cs
   use elpa2_compute_cs
   use elpa_mpi_cs
   use cuda_functions_cs
   use mod_check_for_gpu_cs
   use iso_c_binding

   implicit none

   integer(kind=c_int), intent(in)       :: na,nev,lda,ldq,nblk,matrixCols
   integer(kind=c_int), intent(in)       :: mpi_comm_rows,mpi_comm_cols,mpi_comm_all
   real(kind=c_double), intent(inout)    :: ev(na)
   complex(kind=c_double), intent(inout) :: a(lda,*),q(ldq,*)
   real(kind=c_double), intent(in)       :: singular_tol
   integer(kind=c_int), intent(out)      :: n_nonsing

   integer(kind=c_int)                   :: THIS_COMPLEX_ELPA_KERNEL
   complex(kind=c_double), allocatable   :: hh_trans_complex(:,:)
   integer(kind=c_int)                   :: my_prow,my_pcol,np_rows,np_cols,mpierr,my_pe,n_pes
   integer(kind=c_int)                   :: l_cols,l_rows,l_cols_nev,nbw,num_blocks
   complex(kind=c_double), allocatable   :: tmat(:,:,:)
   real(kind=c_double), allocatable      :: q_real(:,:),e(:)
   real(kind=c_double)                   :: ttt0,ttt1
   integer(kind=c_int)                   :: i
   logical                               :: success,wantDebug
   integer(kind=c_int)                   :: istat
   logical                               :: do_useGPU,do_useGPU_4
   integer(kind=c_int)                   :: numberOfGPUDevices

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)
   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   do_useGPU   = .false.
   do_useGPU_4 = .false.
   wantDebug   = .false.
   success     = .true.

   THIS_COMPLEX_ELPA_KERNEL = get_actual_complex_kernel()

   if(THIS_COMPLEX_ELPA_KERNEL .eq. COMPLEX_ELPA_KERNEL_GPU) then
      if(nblk .ne. 128) then
         print *, " ERROR: ELPA2 GPU kernel requires blocksize = 128. Exiting.."
         stop
      endif
      ! 4th step
      do_useGPU_4 = .true.
   endif

#ifdef WITH_GPU_VERSION
   if(check_for_gpu(my_pe,numberOfGPUDevices,wantDebug=wantDebug)) then
      do_useGPU = .true.
   endif

   ! Set GPU/CUDA parameters
   cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
   cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
   cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
   cudaHostRegisterPortable = cuda_hostRegisterPortable()
   cudaHostRegisterMapped   = cuda_hostRegisterMapped()
#endif

   ! Bandwidth must be a multiple of nblk
   ! Set to a value >= 32
   nbw = (31/nblk+1)*nblk
   num_blocks = (na-1)/nbw+1

   allocate(tmat(nbw,nbw,num_blocks),stat=istat)
   if(istat .ne. 0) then
      print *," ERROR: Cannot allocate tmat. Exiting.."
      stop
   endif

   ! Reduction full -> band
   ttt0 = MPI_Wtime()
   call bandred_complex_double(na,a,lda,nblk,nbw,matrixCols,num_blocks,mpi_comm_rows,&
                               mpi_comm_cols,tmat,wantDebug,do_useGPU,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      write(*,"(A,F10.3,A)") "  | Time full => band           :",ttt1-ttt0," s"

   ! Reduction band -> tridiagonal
   allocate(e(na),stat=istat)
   if(istat .ne. 0) then
      print *," ERROR: Cannot allocate e. Exiting.."
      stop
   endif

   ttt0 = MPI_Wtime()
   call tridiag_band_complex_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_complex,&
                                    mpi_comm_rows,mpi_comm_cols,mpi_comm_all)
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      write(*,"(A,F10.3,A)") "  | Time band => tridiagonal    :",ttt1-ttt0," s"

   call mpi_bcast(ev,na,mpi_real8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,mpi_real8,0,mpi_comm_all,mpierr)

   l_rows = local_index(na,my_prow,np_rows,nblk,-1) ! Local rows of a and q
   l_cols = local_index(na,my_pcol,np_cols,nblk,-1) ! Local columns of q
   l_cols_nev = local_index(nev,my_pcol,np_cols,nblk,-1) ! Local columns corresponding to nev

   allocate(q_real(l_rows,l_cols),stat=istat)
   if(istat .ne. 0) then
      print *," ERROR: Cannot allocate q_real. Exiting.."
      stop
   endif

   ! Solve tridiagonal system
   ttt0 = MPI_Wtime()
   call solve_tridi_double(na,nev,ev,e,q_real,ubound(q_real,dim=1),nblk,matrixCols,&
                           mpi_comm_rows,mpi_comm_cols,wantDebug,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
   if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
      write(*,"(A,F10.3,A)") "  | Time solve tridiagonal      :",ttt1-ttt0," s"

   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   deallocate(e)
   deallocate(q_real)

   ! Get the number of nonsingular eigenvalues
   n_nonsing = na
   do i = 1,na
      if(ev(i) < singular_tol) then
         n_nonsing = n_nonsing-1
      endif
   enddo

   if(n_nonsing < na) then
      ! Back-transform 1
      ttt0 = MPI_Wtime()
      call trans_ev_tridi_to_band_complex_double(na,nev,nblk,nbw,q,ldq,matrixCols,&
                                                 hh_trans_complex,mpi_comm_rows,&
                                                 mpi_comm_cols,wantDebug,do_useGPU_4,&
                                                 success,THIS_COMPLEX_ELPA_KERNEL)
      if(.not.success) return
      ttt1 = MPI_Wtime()
      if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
         write(*,"(A,F10.3,A)") "  | Time ev tridiagonal => band :",ttt1-ttt0," s"

      deallocate(hh_trans_complex)

      ! Back-transform 2
      ttt0 = MPI_Wtime()
      call trans_ev_band_to_full_complex_double(na,nev,nblk,nbw,a,lda,tmat,q,ldq,&
                                                matrixCols,num_blocks,mpi_comm_rows,&
                                                mpi_comm_cols,do_useGPU)
      ttt1 = MPI_Wtime()
      if(my_prow==0 .and. my_pcol==0 .and. elpa_print_times) &
         write(*,"(A,F10.3,A)") "  | Time ev band => full        :",ttt1-ttt0," s"

      deallocate(tmat)
   else
      deallocate(hh_trans_complex)
      deallocate(tmat)
   endif

   if(my_prow==0 .and. my_pcol==0) then
      if(do_useGPU) then
         if(do_useGPU_4) then
            print *,"  GPU has been used for the 1st, 4th, 5th steps of ELPA2"
         else
            print *,"  GPU has been used for the 1st, 5th steps of ELPA2"
         endif
      endif
   endif

end function

end module CHECK_SINGULARITY
