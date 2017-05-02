subroutine blacs_gridinit_py_wrapper(icontxt_in, order, nprow, npcol, icontxt_out)
  implicit none

  ! f2py arguments
  integer,   intent(in)  :: icontxt_in
  character, intent(in)  :: order
  integer,   intent(in)  :: nprow
  integer,   intent(in)  :: npcol
  integer,   intent(out) :: icontxt_out

  icontxt_out = icontxt_in
  call BLACS_Gridinit(icontxt_out, order, nprow, npcol)

end subroutine blacs_gridinit_py_wrapper

subroutine ms_scalapack_setup_py_wrapper(mpi_comm,nprow,order,bs_def,bs_list_present,bs_list,icontxt_present,icontxt,size_bs_list)
  use MatrixSwitch

  implicit none

  ! Array dimensions hidden by f2py
  integer, intent(in) :: size_bs_list
  ! f2py arguments
  character(1), intent(in) :: order 
  integer, intent(in) :: mpi_comm
  integer, intent(in) :: nprow 
  integer, intent(in) :: bs_def 
  logical, intent(in) :: bs_list_present
  integer, intent(in) :: bs_list(size_bs_list)
  logical, intent(in) :: icontxt_present
  integer, intent(in) :: icontxt

  call ms_scalapack_setup_no_opt(mpi_comm,nprow,order,bs_def,bs_list_present,bs_list,icontxt_present,icontxt)

end subroutine ms_scalapack_setup_py_wrapper
