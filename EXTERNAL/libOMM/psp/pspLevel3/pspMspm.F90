MODULE pspMspm
  use pspVariable
  use pspUtility
  use pspMPI
  use pspLevel1
  use pspLevel2
use pspMspm_nn, only: psp_gemspm_nn
use pspMspm_nt, only: psp_gemspm_nt
use pspMspm_tn, only: psp_gemspm_tn
use pspMspm_tt, only: psp_gemspm_tt

#ifdef MPI
  include 'mpif.h'
#endif

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**********************************************!

  integer, external :: numroc ! it is a function to compute local size

  !**** INTERFACES ********************************!

  interface psp_gemspm
     module procedure psp_dgemspm
     module procedure psp_zgemspm
  end interface psp_gemspm

  interface die
     module procedure die
  end interface die

  public :: psp_gemspm

contains

  subroutine psp_dgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    real(dp), intent(in)  ::   alpha, beta ! scalar
    real(dp), intent(in) :: A(:,:)

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(in) ::   B ! matrix B
    real(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!

    !logical :: changeFmtB
    integer :: trA, trB, ot

    !**********************************************!
    if (alpha/=0.0_dp) then
       !if (B%str_type=='coo') then
       !   changeFmtB=.true.
       !   call psp_coo2csc(B)
       !else
       !   changeFmtB=.false.
       !end if
       call psp_process_opM(opA,trA)
       call psp_process_opM(opB,trB)
       ! operation table
       if (trA==0 .and. trB==0) then
          ot=1
       else if (trA==0 .and. trB>=1) then
          ot=2
       else if (trA>=1 .and. trB==0) then
          ot=3
       else if (trA>=1 .and. trB>=1) then
          ot=4
       else
          call die('mm_dmultiply: invalid implementation')
       end if

       select case (ot)
       case (1)
          call psp_gemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          call psp_gemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
       ! change format back to plan
       !if (changeFmtB) then
       !   call psp_csc2coo(B)
       !end if
    else
       if (beta/=0.0_dp) C=beta*C
    end if

  end subroutine psp_dgemspm


  subroutine psp_zgemspm(M,N,K,A,opA,B,opB,C,alpha,beta)
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) ::            M, N, K ! Globally, A of size M by K, B of size K by N
    character(1), intent(in) :: opA ! form of op(A): 'n/N' for A, 't/T/c/C' for A^T
    character(1), intent(in) :: opB ! form of op(B)
    complex(dp), intent(in)  ::   alpha, beta ! scalar
    complex(dp), intent(in) :: A(:,:)

    !**** INOUT ***********************************!

    type(psp_matrix_spm), intent(in) ::   B ! matrix B
    complex(dp), intent(inout) :: C(:,:)

    !**** LOCAL ***********************************!
    !logical :: changeFmtB
    integer :: trA, trB, ot

    !**********************************************!
    if (alpha/=cmplx_0) then
       !if (B%str_type=='coo') then
       !   changeFmtB=.true.
       !   call psp_coo2csc(B)
       !else
       !   changeFmtB=.false.
       !end if
       call psp_process_opM(opA,trA)
       call psp_process_opM(opB,trB)
       ! operation table
       if (trA==0 .and. trB==0) then
          ot=1
       else if (trA==0 .and. trB>=1) then
          ot=2
       else if (trA>=1 .and. trB==0) then
          ot=3
       else if (trA>=1 .and. trB>=1) then
          ot=4
       else
          call die('mm_dmultiply: invalid implementation')
       end if

       select case (ot)
       case (1)
          call psp_gemspm_nn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (2)
          call psp_gemspm_nt(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (3)
          call psp_gemspm_tn(M,N,K,A,opA,B,opB,C,alpha,beta)
       case (4)
          call psp_gemspm_tt(M,N,K,A,opA,B,opB,C,alpha,beta)
       end select
       ! change format back to plan
       !if (changeFmtB) then
       !   call psp_csc2coo(B)
       !end if
    else
       if (beta/=cmplx_0) C=beta*C
    end if
  end subroutine psp_zgemspm

  subroutine die(message)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in), optional :: message

    !**** INTERNAL ********************************!

    logical, save :: log_start=.false.

    integer :: log_unit

    !**********************************************!

    if (log_start) then
       open(newunit=log_unit,file='MatrixSwitch.log',position='append')
    else
       open(newunit=log_unit,file='MatrixSwitch.log',status='replace')
       log_start=.true.
    end if
    write(log_unit,'(a)'), 'FATAL ERROR in matrix_switch!'
    write(log_unit,'(a)'), message
    close(log_unit)
    stop

  end subroutine die


END MODULE pspMspm
