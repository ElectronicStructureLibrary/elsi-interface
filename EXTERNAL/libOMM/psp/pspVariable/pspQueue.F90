MODULE pspQueue

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** TYPES *************************************!

type psp_dLnkNode
integer :: rowidx, colidx
real(dp) :: val
type(psp_dLnkNode), pointer :: next
end type psp_dLnkNode

type psp_zLnkNode
integer :: rowidx, colidx
complex(dp) :: val
type(psp_dLnkNode), pointer :: next
end type psp_zLnkNode

type psp_dQueue
type(psp_dLnkNode), pointer :: first
integer :: length
logical :: is_initialized=.false.
end type psp_dQueue

type psp_zQueue
type(psp_zLnkNode), pointer :: first
integer :: length
logical :: is_initialized=.false.
end type psp_zQueue


  !**** INTERFACES ********************************!

  interface psp_queue_init
     module procedure psp_dqueue_init
     module procedure psp_zqueue_init
  end interface psp_queue_init

  interface psp_queue_deallocate
module procedure psp_dqueue_deallocate
module procedure psp_zqueue_deallocate
  end interface psp_queue_deallocate

  public :: psp_matrix_spm
  public :: psp_register_spm
  public :: psp_deallocate_spm

contains

  subroutine psp_dqueue_init(q_name,val)
    implicit none

    !**** INPUT ***********************************!
    real(dp), intent(in), target :: val ! initialize a queue with the first data val

    !**** INOUT ***********************************!

    type(psp_dQueue), intent(inout) :: m_name ! a queue to be initialized

!**** LOCAL ***********************************!

type(psp_dLnkNode) :: node ! first data node

    !**********************************************!
    if(.not.(associated(node))) allocate(node)
node%val=val
nullify(node%next)
node=>m_name%first
    m_name%length=1
    m_name%is_initialized=.true.

  end subroutine psp_dqueue_init

subroutine psp_zqueue_init(q_name,val)
implicit none

!**** INPUT ***********************************!
complex(dp), intent(in), target :: val ! initialize a queue with the first data val

!**** INOUT ***********************************!

type(psp_zQueue), intent(inout) :: m_name ! a queue to be initialized

!**** LOCAL ***********************************!

type(psp_zLnkNode) :: node ! first data node

!**********************************************!
if(.not.(associated(node))) allocate(node)
node%val=val
nullify(node%next)
node=>m_name%first
m_name%length=1
m_name%is_initialized=.true.

end subroutine psp_zqueue_init

subroutine psp_dqueue_deallocate(q_name)
implicit none

!**** INPUT ***********************************!

type(psp_dQueue), intent(in) :: q_name ! a queue to be dellocated

!**** LOCAL ***********************************!

integer :: cnt

!**********************************************!

do cnt=1,q_name%length

nullify(node%next)
node=>m_name%first
m_name%length=0
q_name%is_initialized=.false.

end subroutine psp_dqueue_deallocate


END MODULE pspQueue
