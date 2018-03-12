module elsi_elpa_api

  use ELSI_DATATYPE,  only: elsi_handle
  use ELSI_IO,        only: elsi_say
  use ELSI_MALLOC,    only: elsi_allocate,elsi_deallocate
  use ELSI_MPI,       only: elsi_stop
  use ELSI_PRECISION, only: r4,r8,i4
  use ELPA                             
  implicit none


  INTERFACE elpa_eigenvectors
  MODULE PROCEDURE elpa_eigenvectors_d
  MODULE PROCEDURE elpa_eigenvectors_s
  MODULE PROCEDURE elpa_eigenvectors_cd
  MODULE PROCEDURE elpa_eigenvectors_cs
  END INTERFACE
  INTERFACE elpa_hermitian_multiply
  MODULE PROCEDURE elpa_hermitian_multiply_c
  MODULE PROCEDURE elpa_hermitian_multiply_r
  END INTERFACE

  INTERFACE elpa_cholesky
  MODULE PROCEDURE elpa_cholesky_r
  MODULE PROCEDURE elpa_cholesky_c
  END INTERFACE

  INTERFACE elpa_invert_triangular
  MODULE PROCEDURE elpa_invert_triangular_r
  MODULE PROCEDURE elpa_invert_triangular_c
  END INTERFACE

contains
! EIGENVECTOR CALLS
! real double
subroutine elpa_eigenvectors_d(e_h,ham,eval,evec,solver)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle  
   class(elpa_t), pointer     :: elpa_t
   real(kind=r8),     intent(in)    :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)
   integer, optional, intent(in)    :: solver
   real(kind=r8), allocatable :: copy_ham(:,:)
   character(len=40), parameter :: caller = "elpa_eigenvectors_d"
   integer(kind=i4) :: ierr
 
   call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,"copy_ham",caller)
   copy_ham=ham
   if(present(solver)) then
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,solver)
   else
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,e_h%elpa_solver)
   endif
 
! First tests with elpa autotuning
!   tune_state => e%autotune_setup(ELPA_AUTOTUNE_FAST, AUTOTUNE_DOMAIN, ierr)
!   assert_elpa_ok(ierr)
!   do while (e%autotune_step(tune_state))
!     call e%eigenvectors(ham, eval, evec, ierr)
!     assert_elpa_ok(error)
!     copy_ham = ham
!     print *, ""
!   end do



   call elpa_t%eigenvectors(copy_ham,eval,evec,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   call elsi_deallocate(e_h,copy_ham,"copy_ham")
   if(ierr /= 0) then
      call elsi_stop(" ELPA EVP solver failed.",e_h,caller)
   endif

end subroutine elpa_eigenvectors_d
! real float

subroutine elpa_eigenvectors_s(e_h,ham,eval,evec,solver)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle  
   class(elpa_t), pointer     :: elpa_t
   real(kind=r4),     intent(in) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: eval(e_h%n_basis)
   real(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)
   integer, optional, intent(in)    :: solver
   real(kind=r4), allocatable :: eval_single(:)
   real(kind=r4), allocatable :: evec_single(:,:)
   real(kind=r4), allocatable :: copy_ham(:,:)
   character(len=40), parameter :: caller = "elpa_eigenvectors_s"
   integer(kind=i4) :: ierr

   call elsi_allocate(e_h,eval_single,e_h%n_basis,"eval_single",caller)
   call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,&
           "copy_ham",caller)
   call elsi_allocate(e_h,evec_single,e_h%n_lrow,e_h%n_lcol,&
           "evec_single",caller)
   copy_ham=ham
   if(present(solver)) then
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,solver)
   else
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,e_h%elpa_solver)
   endif

   call elpa_t%eigenvectors(copy_ham,eval_single,evec_single,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   eval = real(eval_single,kind=r8)
   evec = real(evec_single,kind=r8)

   call elsi_deallocate(e_h,eval_single,"eval_single")
   call elsi_deallocate(e_h,copy_ham,"copy_ham")
   call elsi_deallocate(e_h,evec_single,"evec_single")
   if(ierr /= 0) then
      call elsi_stop(" ELPA EVP solver failed.",e_h,caller)
   endif

end subroutine elpa_eigenvectors_s

! real double
subroutine elpa_eigenvectors_cd(e_h,ham,eval,evec,solver)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle  
   class(elpa_t), pointer     :: elpa_t
   complex(kind=r8),     intent(in) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),        intent(inout) :: eval(e_h%n_basis)
   complex(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)
   integer, optional, intent(in)    :: solver

   complex(kind=r8), allocatable :: copy_ham(:,:)
   character(len=40), parameter :: caller = "elpa_eigenvectors_cd"
   integer(kind=i4) :: ierr

   call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,"copy_ham",caller)
   copy_ham=ham
   if(present(solver)) then
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,solver)
   else 
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,e_h%elpa_solver)
   endif

   call elpa_t%eigenvectors(copy_ham,eval,evec,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   call elsi_deallocate(e_h,copy_ham,"copy_ham")
   if(ierr /= 0) then
      call elsi_stop(" ELPA EVP solver failed.",e_h,caller)
   endif

end subroutine elpa_eigenvectors_cd
! real float

subroutine elpa_eigenvectors_cs(e_h,ham,eval,evec,solver)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle  
   class(elpa_t), pointer     :: elpa_t
   complex(kind=r4),     intent(in) :: ham(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),        intent(inout) :: eval(e_h%n_basis)
   complex(kind=r8),     intent(inout) :: evec(e_h%n_lrow,e_h%n_lcol)
   integer, optional, intent(in)    :: solver

   real(kind=r4),    allocatable :: eval_single(:)
   complex(kind=r4), allocatable :: evec_single(:,:)
   complex(kind=r4), allocatable :: copy_ham(:,:)
   character(len=40), parameter :: caller = "elpa_eigenvectors_cs"
   integer(kind=i4) :: ierr

   call elsi_allocate(e_h,eval_single,e_h%n_basis,"eval_single",caller)
   call elsi_allocate(e_h,copy_ham,e_h%n_lrow,e_h%n_lcol,&
           "copy_ham",caller)
   call elsi_allocate(e_h,evec_single,e_h%n_lrow,e_h%n_lcol,&
           "evec_single",caller)
   copy_ham=ham

   if(present(solver)) then
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,solver)
   else
     call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_states_solve,e_h%elpa_solver)
   endif

   call elpa_t%eigenvectors(copy_ham,eval_single,evec_single,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   eval = real(eval_single,kind=r8)
   evec = real(evec_single,kind=r8)

   call elsi_deallocate(e_h,eval_single,"eval_single")
   call elsi_deallocate(e_h,copy_ham,"copy_ham")
   call elsi_deallocate(e_h,evec_single,"evec_single")
   if(ierr /= 0) then
      call elsi_stop(" ELPA EVP solver failed.",e_h,caller)
   endif

end subroutine elpa_eigenvectors_cs

subroutine elpa_hermitian_multiply_r(e_h,uplo,uplo2,a,b,c)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle
   character, intent (in) :: uplo,uplo2

   real(kind=r8),     intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: b(e_h%n_lrow,e_h%n_lcol)
   real(kind=r8),     intent(inout) :: c(e_h%n_lrow,e_h%n_lcol)
   class(elpa_t), pointer     :: elpa_t
   character(len=40), parameter :: caller = "elpa_hermitian_multiply_r"
   integer(kind=i4) :: ierr

   call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_basis,e_h%n_basis,e_h%elpa_solver)
   call elpa_t%hermitian_multiply(uplo,uplo2,e_h%n_basis,a,  &
           b,e_h%n_lrow,e_h%n_lcol,c,e_h%n_lrow,e_h%n_lcol,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   if(ierr /= 0) then
      call elsi_stop(" Matrix multiplication failed.",e_h,caller)
   endif

end subroutine elpa_hermitian_multiply_r
subroutine elpa_hermitian_multiply_c(e_h,uplo,uplo2,a,b,c)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle
   character, intent (in) :: uplo,uplo2

   complex(kind=r8),     intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),     intent(inout) :: b(e_h%n_lrow,e_h%n_lcol)
   complex(kind=r8),     intent(inout) :: c(e_h%n_lrow,e_h%n_lcol)
   class(elpa_t), pointer     :: elpa_t
   character(len=40), parameter :: caller = "elpa_hermitian_multiply_c"
   integer(kind=i4) :: ierr
   real(kind=r8), allocatable :: tmp_real(:,:)

   call elsi_allocate(e_h,tmp_real,e_h%n_lrow,e_h%n_lcol,"tmp_real",caller)

   call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_basis,e_h%n_basis,e_h%elpa_solver)
   call elpa_t%hermitian_multiply(uplo,uplo2,e_h%n_basis,a,  &
           b,e_h%n_lrow,e_h%n_lcol,c,e_h%n_lrow,e_h%n_lcol,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   if(ierr /= 0) then
      call elsi_stop(" Matrix multiplication failed.",e_h,caller)
   endif

end subroutine elpa_hermitian_multiply_c


! Cholesky subroutines
! real double
subroutine elpa_cholesky_r(e_h,a)
   implicit none
   type(elsi_handle), intent(in)      :: e_h !< Handle
   real(kind=r8), intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   class(elpa_t), pointer   :: elpa_t
   integer(kind=i4) :: ierr
   character(len=40), parameter :: caller = "elpa_cholesky_r"

   call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_basis,e_h%n_basis)
   call elpa_t%cholesky(a,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   if(ierr /= 0) then
      call elsi_stop(" Cholesky decomposition failed.",e_h,caller)
   endif

end subroutine elpa_cholesky_r
subroutine elpa_cholesky_c(e_h,a)
   implicit none
   type(elsi_handle), intent(in)      :: e_h !< Handle
   complex(kind=r8), intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   class(elpa_t), pointer   :: elpa_t
   integer(kind=i4) :: ierr
   character(len=40), parameter :: caller = "elpa_cholesky_c"

   call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_basis,e_h%n_basis)
   call elpa_t%cholesky(a,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   if(ierr /= 0) then
      call elsi_stop(" Cholesky decomposition failed.",e_h,caller)
   endif

end subroutine elpa_cholesky_c


! Invert triangular subroutines
!real double
subroutine elpa_invert_triangular_r(e_h,a)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle    
   real(kind=r8), intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   class(elpa_t), pointer     :: elpa_t
   integer(kind=i4) :: ierr
   character(len=40), parameter :: caller = "elpa_invert_triangular_r"

   call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_nonsing)
   call elpa_t%invert_triangular(a,ierr)
   call elpa_t%destroy()
   call elpa_uninit()
   if(ierr /= 0) then
      call elsi_stop(" Matrix inversion failed.",e_h,caller)
   endif

end subroutine elpa_invert_triangular_r
! real float
subroutine elpa_invert_triangular_c(e_h,a)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle   
   complex(kind=r8), intent(inout) :: a(e_h%n_lrow,e_h%n_lcol)
   class(elpa_t), pointer     :: elpa_t
   integer(kind=i4) :: ierr
   character(len=40), parameter :: caller = "elpa_invert_triangular_c"

   call elsi_set_elpa_new_API(e_h,elpa_t,e_h%n_nonsing,e_h%n_nonsing)
   call elpa_t%invert_triangular(a,ierr)

   call elpa_t%destroy()
   call elpa_uninit()
   if(ierr /= 0) then
      call elsi_stop(" Matrix inversion failed.",e_h,caller)
   endif

end subroutine elpa_invert_triangular_c
!>
!! This routine sets the  ELPA API parameters.
!!

subroutine elsi_set_elpa_new_API(e_h,elpa_i,na,nev,solver)
   implicit none
   type(elsi_handle), intent(in) :: e_h !< Handle
   integer, intent (in) :: na,nev
   integer, intent (in), optional :: solver
   integer(kind=i4) :: ierr
   class(elpa_t),intent(inout),pointer      :: elpa_i
   character(len=200) :: info_str
   character(len=40), parameter :: caller = "elsi_set_elpa_new_api"

   ierr = elpa_init(20170403)

   elpa_i => elpa_allocate()

   call elpa_i%set("na",na,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'na' (dim of matrix) failed.",e_h,caller)
   endif

   call elpa_i%set("nev",nev,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'nev' (number of eigenvalues&
                      to be computed ) failed.",e_h,caller)
   endif

   call elpa_i%set("local_nrows",e_h%n_lrow,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'local_nrows' failed.",e_h,caller)
   endif
   call elpa_i%set("local_ncols",e_h%n_lcol,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'local_ncols' failed.",e_h,caller)
   endif

   call elpa_i%set("nblk",e_h%blk_row,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'nblk' failed.",e_h,caller)
   endif

   call elpa_i%set("mpi_comm_parent",e_h%mpi_comm,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'mpi_comm_parent' failed.",e_h,caller)
   endif

   call elpa_i%set("process_row",e_h%my_prow,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'process_row' failed.",e_h,caller)
   endif

   call elpa_i%set("process_col",e_h%my_pcol,ierr)
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup of 'process col' failed.",e_h,caller)
   endif

   ierr = elpa_i%setup()
   if(ierr /= 0) then
      call elsi_stop(" ELPA setup failed.",e_h,caller)
   endif
   if(present(solver)) then
      if(solver == 1) then
         call elpa_i%set("solver",ELPA_SOLVER_1STAGE,ierr)
      else if (solver == 2) then
         call elpa_i%set("solver",ELPA_SOLVER_2STAGE,ierr)
      else
         write(info_str,"('  Setting the default solver: ELPA_SOLVER_2STAGE')")
         call elsi_say(e_h,info_str)
         call elpa_i%set("solver",ELPA_SOLVER_2STAGE,ierr)
      endif
      if(ierr /= 0) then
         call elsi_stop(" ELPA setup of 'solver' failed.",e_h,caller)
      endif
   call elpa_i%set("check_pd",1,ierr)

   endif
end subroutine

end module
