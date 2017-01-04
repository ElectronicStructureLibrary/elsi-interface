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
!
! ELPA2 -- 2-stage solver for ELPA
!
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".

module ELPA2

  use elpa_utilities
  use elpa1, only : elpa_print_times,time_evp_back,time_evp_fwd,time_evp_solve
  use elpa2_utilities

  implicit none

  private

  public :: solve_evp_real_2stage_double
  public :: solve_evp_complex_2stage_double

  interface solve_evp_real_2stage
     module procedure solve_evp_real_2stage_double
  end interface

  interface solve_evp_complex_2stage
     module procedure solve_evp_complex_2stage_double
  end interface

  ! Added by Victor Yu from ELSI team.
  ! Dec 6, 2016
  public :: check_eval_real
  public :: check_eval_complex

contains

!-------------------------------------------------------------------------------
!>  \param na                Order of matrix a
!>
!>  \param nev               Number of eigenvalues needed
!>
!>  \param a(lda,*)          Distributed matrix of which eigenvalues to be computed
!>                           Distribution is like in ScaLAPACK
!>                           The full matrix must be set
!>                           Destroyed on exit (upper and lower half)
!>
!>  \param lda               Leading dimension of a
!>
!>  \param ev(na)            Output: eigenvalues of a
!>                           All processors get the complete set
!>
!>  \param q(ldq,*)          Output: Eigenvectors of a
!>                           Distribution is like in ScaLAPACK
!>                           Must be always dimensioned to the full size (na,na)
!>                           even if only a part of the eigenvalues is needed
!>
!>  \param ldq               Leading dimension of q
!>
!>  \param nblk              Blocksize of cyclic distribution
!>                           Must be the same in both directions
!>
!>  \param matrixCols        Number of local columns of matrix a and q
!>
!>  \param mpi_comm_rows     MPI communicator for rows
!>  \param mpi_comm_cols     MPI communicator for columns
!>  \param mpi_comm_all      MPI communicator for the total processor set
!>
!>  \result success          logical, false if error occured
!-------------------------------------------------------------------------------

function solve_evp_real_2stage_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
                                      mpi_comm_rows,mpi_comm_cols,mpi_comm_all)&
                                      result(success)

   use elpa1_compute
   use elpa2_compute
   use elpa_mpi
   use precision
   use cuda_functions
   use mod_check_for_gpu
   use iso_c_binding

   implicit none

   integer(kind=ik), intent(in)   :: na,nev,lda,ldq,matrixCols
   integer(kind=ik), intent(in)   :: mpi_comm_rows,mpi_comm_cols,mpi_comm_all
   integer(kind=ik), intent(in)   :: nblk
   real(kind=rk8), intent(inout)  :: a(lda,*),q(ldq,*),ev(na)
   real(kind=rk8), allocatable    :: hh_trans_real(:,:)

   logical                        :: useQRActual
   integer(kind=ik)               :: THIS_REAL_ELPA_KERNEL
   integer(kind=ik)               :: my_pe,n_pes,my_prow,my_pcol,np_rows,np_cols,mpierr
   integer(kind=ik)               :: nbw,num_blocks
   real(kind=rk8), allocatable    :: tmat(:,:,:),e(:)
   integer(kind=c_intptr_t)       :: tmat_dev,q_dev,a_dev
   real(kind=c_double)            :: ttt0,ttt1
   integer(kind=ik)               :: i
   logical                        :: success
   logical                        :: wantDebug
   integer(kind=ik)               :: istat
   logical                        :: useGPU
   integer(kind=ik)               :: numberOfGPUDevices

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   wantDebug   = .false.
   success     = .true.
   useQRActual = .false.
   useGPU      = .false.

   THIS_REAL_ELPA_KERNEL = get_actual_real_kernel()

   if(THIS_REAL_ELPA_KERNEL .eq. REAL_ELPA_KERNEL_GPU) then
      if(check_for_gpu(my_pe,numberOfGPUDevices,wantDebug=wantDebug)) then
         useGPU = .true.
      endif
      if(nblk .ne. 128) then
         print *, "*ERROR: At the moment GPU version needs blocksize 128."
         stop
      endif

      ! set the neccessary parameters
      cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
      cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
      cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
      cudaHostRegisterPortable = cuda_hostRegisterPortable()
      cudaHostRegisterMapped   = cuda_hostRegisterMapped()
   endif

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32
   ! On older systems (IBM Bluegene/P, Intel Nehalem) a value of 32 was optimal.
   ! For Intel(R) Xeon(R) E5 v2 and v3, better use 64 instead of 32!
   ! For IBM Bluegene/Q this is not clear at the moment. We have to keep an eye
   ! on this and maybe allow a run-time optimization here
   if(useGPU) then
      nbw = nblk
   else
      nbw = (63/nblk+1)*nblk
   endif

   num_blocks = (na-1)/nbw+1

   allocate(tmat(nbw,nbw,num_blocks),stat=istat)
   if(istat .ne. 0) then
      print *,"solve_evp_real_2stage: error when allocating tmat"
      stop
   endif

   ! Reduction full -> band
   ttt0 = MPI_Wtime()
   call bandred_real_double(na,a,a_dev,lda,nblk,nbw,matrixCols,num_blocks,mpi_comm_rows,&
                            mpi_comm_cols,tmat,tmat_dev,wantDebug,useGPU,success,useQRActual)
   if(.not.success) return
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time full ==> band:",ttt1-ttt0

   ! Reduction band -> tridiagonal
   allocate(e(na),stat=istat)
   if(istat .ne. 0) then
      print *,"solve_evp_real_2stage: error when allocating e"
      stop
   endif

   ttt0 = MPI_Wtime()
   call tridiag_band_real_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_real,&
                                 mpi_comm_rows,mpi_comm_cols,mpi_comm_all)
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time band ==> tridiagonal:",ttt1-ttt0

   call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ! Solve tridiagonal system
   ttt0 = MPI_Wtime()
   call solve_tridi_double(na,nev,ev,e,q,ldq,nblk,matrixCols,mpi_comm_rows,&
                           mpi_comm_cols,wantDebug,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time solve tridiagonal:",ttt1-ttt0

   deallocate(e)

   ! Backtransform stage 1
   ttt0 = MPI_Wtime()
   call trans_ev_tridi_to_band_real_double(na,nev,nblk,nbw,q,q_dev,ldq,matrixCols,&
                                           hh_trans_real,mpi_comm_rows,mpi_comm_cols,&
                                           wantDebug,useGPU,success,THIS_REAL_ELPA_KERNEL)
   if(.not.success) return
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time tridiagonal ==> band:",ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   deallocate(hh_trans_real)

   ! Backtransform stage 2
   ttt0 = MPI_Wtime()
   call trans_ev_band_to_full_real_double(na,nev,nblk,nbw,a,a_dev,lda,tmat,tmat_dev,q,&
                                          q_dev,ldq,matrixCols,num_blocks,mpi_comm_rows,&
                                          mpi_comm_cols,useGPU,useQRActual)
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time band ==> full:",ttt1-ttt0

   deallocate(tmat)

   if(my_prow==0 .and. my_pcol==0) then
      if(useGPU) print *,"  GPU has been used for this ELPA2"
   endif

end function solve_evp_real_2stage_double

!-------------------------------------------------------------------------------
!>  \param na                Order of matrix a
!>
!>  \param nev               Number of eigenvalues needed
!>
!>  \param a(lda,*)          Distributed matrix of which eigenvalues to be computed
!>                           Distribution is like in ScaLAPACK
!>                           The full matrix must be set
!>                           Destroyed on exit (upper and lower half)
!>
!>  \param lda               Leading dimension of a
!>
!>  \param ev(na)            Output: eigenvalues of a
!>                           All processors get the complete set
!>
!>  \param q(ldq,*)          Output: Eigenvectors of a
!>                           Distribution is like in ScaLAPACK
!>                           Must be always dimensioned to the full size (na,na)
!>                           even if only a part of the eigenvalues is needed
!>
!>  \param ldq               Leading dimension of q
!>
!>  \param nblk              Blocksize of cyclic distribution
!>                           Must be the same in both directions
!>
!>  \param matrixCols        Number of local columns of matrix a and q
!>
!>  \param mpi_comm_rows     MPI communicator for rows
!>  \param mpi_comm_cols     MPI communicator for columns
!>  \param mpi_comm_all      MPI communicator for the total processor set
!>
!>  \result success          logical, false if error occured
!-------------------------------------------------------------------------------

function solve_evp_complex_2stage_double(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,&
                                         mpi_comm_rows,mpi_comm_cols,mpi_comm_all)&
                                         result(success)

   use elpa1_compute
   use elpa2_compute
   use elpa_mpi
   use precision
   use cuda_functions
   use mod_check_for_gpu
   use iso_c_binding

   implicit none

   integer(kind=ik)                       :: THIS_COMPLEX_ELPA_KERNEL
   integer(kind=ik), intent(in)           :: na,nev,lda,ldq,nblk,matrixCols
   integer(kind=ik), intent(in)           :: mpi_comm_rows,mpi_comm_cols,mpi_comm_all
   complex(kind=ck8), intent(inout)       :: a(lda,*),q(ldq,*)
   real(kind=rk8), intent(inout)          :: ev(na)
   complex(kind=ck8), allocatable         :: hh_trans_complex(:,:)

   integer(kind=ik)                       :: my_prow,my_pcol,np_rows,np_cols,mpierr,my_pe,n_pes
   integer(kind=ik)                       :: l_cols,l_rows,l_cols_nev,nbw,num_blocks
   complex(kind=ck8), allocatable         :: tmat(:,:,:)
   real(kind=rk8), allocatable            :: q_real(:,:),e(:)
   real(kind=c_double)                    :: ttt0,ttt1
   integer(kind=ik)                       :: i
   logical                                :: success,wantDebug
   integer(kind=ik)                       :: istat
   logical                                :: useGPU
   integer(kind=ik)                       :: numberOfGPUDevices

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   useGPU    = .false.
   wantDebug = .false.
   success   = .true.

   THIS_COMPLEX_ELPA_KERNEL = get_actual_complex_kernel()

   if(THIS_COMPLEX_ELPA_KERNEL .eq. COMPLEX_ELPA_KERNEL_GPU) then
      if(check_for_gpu(my_pe, numberOfGPUDevices, wantDebug=wantDebug)) then
         useGPU=.true.
      endif
      if(nblk .ne. 128) then
         print *,"*ERROR: At the moment GPU version needs blocksize 128."
         stop
      endif

      ! set the neccessary parameters
      cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
      cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
      cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
      cudaHostRegisterPortable = cuda_hostRegisterPortable()
      cudaHostRegisterMapped   = cuda_hostRegisterMapped()
   endif

   ! Choose bandwidth, must be a multiple of nblk, set to a value >= 32
   nbw = (31/nblk+1)*nblk

   num_blocks = (na-1)/nbw+1

   allocate(tmat(nbw,nbw,num_blocks),stat=istat)
   if(istat .ne. 0) then
      print *,"solve_evp_complex_2stage: error when allocating tmat"
      stop
   endif

   ! Reduction full -> band
   ttt0 = MPI_Wtime()
   call bandred_complex_double(na,a,lda,nblk,nbw,matrixCols,num_blocks,mpi_comm_rows,&
                               mpi_comm_cols,tmat,wantDebug,useGPU,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time full ==> band:",ttt1-ttt0

   ! Reduction band -> tridiagonal
   allocate(e(na),stat=istat)
   if(istat .ne. 0) then
      print *,"solve_evp_complex_2stage: error when allocating e"
      stop
   endif

   ttt0 = MPI_Wtime()
   call tridiag_band_complex_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_complex,&
                                    mpi_comm_rows,mpi_comm_cols,mpi_comm_all)
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time band ==> tridiagonal:",ttt1-ttt0

   call mpi_bcast(ev,na,mpi_real8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,mpi_real8,0,mpi_comm_all,mpierr)

   l_rows = local_index(na,my_prow,np_rows,nblk,-1) ! Local rows of a and q
   l_cols = local_index(na,my_pcol,np_cols,nblk,-1) ! Local columns of q
   l_cols_nev = local_index(nev,my_pcol,np_cols,nblk,-1) ! Local columns corresponding to nev

   allocate(q_real(l_rows,l_cols),stat=istat)
   if(istat .ne. 0) then
      print *,"solve_evp_complex_2stage: error when allocating q_real"
      stop
   endif

   ! Solve tridiagonal system
   ttt0 = MPI_Wtime()
   call solve_tridi_double(na,nev,ev,e,q_real,ubound(q_real,dim=1),nblk,matrixCols,&
                           mpi_comm_rows,mpi_comm_cols,wantDebug,success)
   if(.not.success) return
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time solve tridiagonal:",ttt1-ttt0

   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   deallocate(e)
   deallocate(q_real)

   ! Backtransform stage 1
   ttt0 = MPI_Wtime()
   call trans_ev_tridi_to_band_complex_double(na,nev,nblk,nbw,q,ldq,matrixCols,&
                                              hh_trans_complex,mpi_comm_rows,&
                                              mpi_comm_cols,wantDebug,useGPU,&
                                              success,THIS_COMPLEX_ELPA_KERNEL)
   if(.not.success) return
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time tridiagonal ==> band:",ttt1-ttt0

   ! We can now deallocate the stored householder vectors
   deallocate(hh_trans_complex)

   ! Backtransform stage 2
   ttt0 = MPI_Wtime()
   call trans_ev_band_to_full_complex_double(na,nev,nblk,nbw,a,lda,tmat,q,ldq,&
                                             matrixCols,num_blocks,mpi_comm_rows,&
                                             mpi_comm_cols,useGPU)
   ttt1 = MPI_Wtime()
!   if(my_prow==0 .and. my_pcol==0) print *,"  Time band ==> full:",ttt1-ttt0

   deallocate(tmat)

   if(my_prow==0 .and. my_pcol==0) then
      if(useGPU) print *,"  GPU has been used for this ELPA2"
   endif

end function solve_evp_complex_2stage_double

! Added by Victor Yu from ELSI team
! Dec 6, 2016
function check_eval_real(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,mpi_comm_rows,&
                         mpi_comm_cols,mpi_comm_all,singular_tol,n_nonsing)&
                         result(success)

   use elpa1_compute
   use elpa2_compute
   use elpa_mpi
   use cuda_functions
   use mod_check_for_gpu
   use iso_c_binding
   use precision

   implicit none

   integer(kind=ik), intent(in)  :: na,nev,lda,ldq,matrixCols
   integer(kind=ik), intent(in)  :: mpi_comm_rows,mpi_comm_cols,mpi_comm_all
   integer(kind=ik), intent(in)  :: nblk
   integer(kind=ik), intent(out) :: n_nonsing
!   real(kind=rk8), intent(inout) :: a(lda,matrixCols),ev(na),q(ldq,matrixCols)
   real(kind=rk8), intent(inout) :: a(lda,*),ev(na),q(ldq,*)
   real(kind=rk8), intent(in)    :: singular_tol

   real(kind=rk8), allocatable   :: hh_trans_real(:,:)
   real(kind=rk8), allocatable   :: tmat(:,:,:),e(:),ev_tmp(:)
   logical                       :: useQRActual
   logical                       :: success
   logical                       :: wantDebug
   logical                       :: useGPU
   integer(kind=ik)              :: THIS_REAL_ELPA_KERNEL
   integer(kind=ik)              :: i
   integer(kind=ik)              :: my_pe,n_pes,my_prow,my_pcol,np_rows,np_cols
   integer(kind=ik)              :: mpierr
   integer(kind=ik)              :: nbw,num_blocks
   integer(kind=ik)              :: istat
   integer(kind=ik)              :: numberOfGPUDevices
   integer(kind=c_intptr_t)      :: tmat_dev,q_dev,a_dev

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   success     = .true.
   wantDebug   = .false.
   useQRActual = .false.
   useGPU      = .false.

   THIS_REAL_ELPA_KERNEL = get_actual_real_kernel()

   if(THIS_REAL_ELPA_KERNEL .eq. REAL_ELPA_KERNEL_GPU) then
      if(check_for_gpu(my_pe,numberOfGPUDevices,wantDebug=wantDebug)) then
         useGPU = .true.
      endif
      if(nblk .ne. 128) then
         print *, " ERROR: ELPA GPU version needs blocksize = 128. Exiting.."
         stop
      endif

      ! Set GPU/CUDA parameters
      cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
      cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
      cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
      cudaHostRegisterPortable = cuda_hostRegisterPortable()
      cudaHostRegisterMapped   = cuda_hostRegisterMapped()
   endif

   ! Bandwidth must be a multiple of nblk
   ! Set to a value >= 32
   if(useGPU) then
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
   call bandred_real_double(na,a,a_dev,lda,nblk,nbw,matrixCols,num_blocks,&
                            mpi_comm_rows,mpi_comm_cols,tmat,tmat_dev,wantDebug,&
                            useGPU,success,useQRActual)
   if(.not.success) return

   ! Reduction band -> tridiagonal
   allocate(e(na),stat=istat)
   if(istat .ne. 0) then
      print *," ERROR: Cannot allocate e. Exiting.."
      stop
   endif

   ! Reduction band -> tridiagonal
   call tridiag_band_real_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_real,&
                                 mpi_comm_rows,mpi_comm_cols,mpi_comm_all)

   call mpi_bcast(ev,na,MPI_REAL8,0,mpi_comm_all,mpierr)
   call mpi_bcast(e,na,MPI_REAL8,0,mpi_comm_all,mpierr)

   ! Solve tridiagonal system
   call solve_tridi_double(na,nev,ev,e,q,ldq,nblk,matrixCols,mpi_comm_rows,&
                           mpi_comm_cols,wantDebug,success)
   if(.not.success) return

   deallocate(e)

   ! Invert signs of eigenvalues
   ev = -ev

   ! Get the number of nonsingular eigenvalues
   do i = 1,na
      if(ev(i) < singular_tol) exit
   enddo
   n_nonsing = i-1

   if(n_nonsing < na) then
      ! Back-transform 1
      call trans_ev_tridi_to_band_real_double(na,nev,nblk,nbw,q,q_dev,ldq,matrixCols,&
                                              hh_trans_real,mpi_comm_rows,mpi_comm_cols,&
                                              wantDebug,useGPU,success,THIS_REAL_ELPA_KERNEL)
      if(.not.success) return

      deallocate(hh_trans_real)

      ! Back-transform 2
      call trans_ev_band_to_full_real_double(na,nev,nblk,nbw,a,a_dev,lda,tmat,tmat_dev,q,&
                                             q_dev,ldq,matrixCols,num_blocks,mpi_comm_rows,&
                                             mpi_comm_cols,useGPU,useQRActual)

      deallocate(tmat)
   else
      deallocate(hh_trans_real)
      deallocate(tmat)
   endif

   if(my_prow==0 .and. my_pcol==0) then
      if(useGPU) print *,"  GPU has been used for this ELPA2"
   endif

end function

function check_eval_complex(na,nev,a,lda,ev,q,ldq,nblk,matrixCols,mpi_comm_rows,&
                            mpi_comm_cols,mpi_comm_all,singular_tol,n_nonsing) result(success)

   use elpa1_compute
   use elpa2_compute
   use elpa_mpi
   use cuda_functions
   use mod_check_for_gpu
   use iso_c_binding
   use precision

   implicit none

   integer(kind=ik), intent(in)     :: na,nev,lda,ldq,nblk,matrixCols
   integer(kind=ik), intent(in)     :: mpi_comm_rows,mpi_comm_cols,mpi_comm_all
   integer(kind=ik), intent(out)    :: n_nonsing
!   complex(kind=ck8), intent(inout) :: a(lda,matrixCols),q(ldq,matrixCols)
   complex(kind=ck8), intent(inout) :: a(lda,*),q(ldq,*)
   real(kind=rk8), intent(inout)    :: ev(na)
   real(kind=rk8), intent(in)       :: singular_tol

   complex(kind=ck8), allocatable   :: hh_trans_complex(:,:)
   complex(kind=ck8), allocatable   :: tmat(:,:,:)
   real(kind=rk8), allocatable      :: q_real(:,:),e(:)
   integer(kind=ik)                 :: my_prow,my_pcol,np_rows,np_cols,mpierr,my_pe,n_pes
   integer(kind=ik)                 :: l_cols,l_rows,l_cols_nev,nbw,num_blocks
   integer(kind=ik)                 :: THIS_COMPLEX_ELPA_KERNEL
   integer(kind=ik)                 :: i
   integer(kind=ik)                 :: istat
   integer(kind=ik)                 :: numberOfGPUDevices
   logical                          :: success
   logical                          :: wantDebug
   logical                          :: useGPU

   call mpi_comm_rank(mpi_comm_all,my_pe,mpierr)
   call mpi_comm_size(mpi_comm_all,n_pes,mpierr)

   call mpi_comm_rank(mpi_comm_rows,my_prow,mpierr)
   call mpi_comm_size(mpi_comm_rows,np_rows,mpierr)
   call mpi_comm_rank(mpi_comm_cols,my_pcol,mpierr)
   call mpi_comm_size(mpi_comm_cols,np_cols,mpierr)

   success   = .true.
   useGPU    = .false.
   wantDebug = .false.

   THIS_COMPLEX_ELPA_KERNEL = get_actual_complex_kernel()

   if(THIS_COMPLEX_ELPA_KERNEL .eq. COMPLEX_ELPA_KERNEL_GPU) then
      if(check_for_gpu(my_pe,numberOfGPUDevices,wantDebug=wantDebug)) then
         useGPU=.true.
      endif
      if(nblk .ne. 128) then
         print *, " ERROR: ELPA GPU version needs blocksize = 128. Exiting.."
         stop
      endif

      ! Set GPU/CUDA parameters
      cudaMemcpyHostToDevice   = cuda_memcpyHostToDevice()
      cudaMemcpyDeviceToHost   = cuda_memcpyDeviceToHost()
      cudaMemcpyDeviceToDevice = cuda_memcpyDeviceToDevice()
      cudaHostRegisterPortable = cuda_hostRegisterPortable()
      cudaHostRegisterMapped   = cuda_hostRegisterMapped()
   endif

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
   call bandred_complex_double(na,a,lda,nblk,nbw,matrixCols,num_blocks,mpi_comm_rows,&
                               mpi_comm_cols,tmat,wantDebug,useGPU,success)
   if(.not.success) return

   ! Reduction band -> tridiagonal
   allocate(e(na),stat=istat)
   if(istat .ne. 0) then
      print *," ERROR: Cannot allocate e. Exiting.."
      stop
   endif

   call tridiag_band_complex_double(na,nbw,nblk,a,lda,ev,e,matrixCols,hh_trans_complex,&
                                    mpi_comm_rows,mpi_comm_cols,mpi_comm_all)

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
   call solve_tridi_double(na,nev,ev,e,q_real,ubound(q_real,dim=1),nblk,matrixCols,&
                           mpi_comm_rows,mpi_comm_cols,wantDebug,success)
   if(.not.success) return

   q(1:l_rows,1:l_cols_nev) = q_real(1:l_rows,1:l_cols_nev)

   deallocate(e)
   deallocate(q_real)

   ! Invert signs of eigenvalues
   ev = -ev

   ! Get the number of nonsingular eigenvalues
   do i = 1,na
      if(ev(i) < singular_tol) exit
   enddo
   n_nonsing = i-1

   if(n_nonsing < na) then
      ! Back-transform 1
      call trans_ev_tridi_to_band_complex_double(na,nev,nblk,nbw,q,ldq,matrixCols,&
                                                 hh_trans_complex,mpi_comm_rows,&
                                                 mpi_comm_cols,wantDebug,useGPU,&
                                                 success,THIS_COMPLEX_ELPA_KERNEL)
      if(.not.success) return

      deallocate(hh_trans_complex)

      ! Back-transform 2
      call trans_ev_band_to_full_complex_double(na,nev,nblk,nbw,a,lda,tmat,q,ldq,&
                                                matrixCols,num_blocks,mpi_comm_rows,&
                                                mpi_comm_cols,useGPU)

      deallocate(tmat)
   else
      deallocate(hh_trans_complex)
      deallocate(tmat)
   endif

   if(my_prow==0 .and. my_pcol==0) then
      if(useGPU) print *,"  GPU has been used for this ELPA2"
   endif

end function

end module ELPA2
