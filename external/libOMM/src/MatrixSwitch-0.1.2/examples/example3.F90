!==================================================================================================!
! example3 : parallel program with real matrices                                                   !
!                                                                                                  !
! This example demonstrates how to calculate:                                                      !
!   alpha = tr(A^T*B)                                                                              !
! for two NxM matrices A and B, in five different ways. Each way should return the same result.    !
! They are:                                                                                        !
!   1. performing the matrix product trace res1 := tr(A^T*B) directly as a single operation        !
!   2. performing the multiplication D := A^T*B, and then the trace res2 := tr(D)                  !
!   3. performing E := B*A^T, and then res3 := tr(E)                                               !
!   4. performing the transpose C := A^T, then D := C*B, and then res4 := tr(D)                    !
!   5. performing C := A^T, then E := B*C, and then res5 := tr(E)                                  !
! Finally, as an extra check, the first element of E is printed out.                               !
!                                                                                                  !
! Note the difference in the code in how matrix B is handled compared with the other matrices.     !
! While the others are allocated directly by MatrixSwitch (i.e., all the data is contained within  !
! the type(matrix) variable), B is a wrapper for MyMatrix (i.e., MyMatrix has been registered as B !
! for use with MatrixSwitch, and B contains a pointer to MyMatrix).                                !
!                                                                                                  !
! Sample output:                                                                                   !
!--------------------------------------------------------------------------------------------------!
! res1 :    31.194861937321843                                                                     !
! res2 :    31.194861937321846                                                                     !
! res3 :    31.194861937321846                                                                     !
! res4 :    31.194861937321846                                                                     !
! res5 :    31.194861937321846                                                                     !
! E_11 :     2.779421970931733                                                                     !
!--------------------------------------------------------------------------------------------------!
!==================================================================================================!
program example3
  use MatrixSwitch

  implicit none
#if defined(MPI) && defined(SLAP)
  include 'mpif.h'

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!

  character(5) :: m_storage
  character(3) :: m_operation

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: icontxt, N_loc, M_loc, MyDesc(9), info
  integer :: N, M, i, j, k, l
  integer, allocatable :: bs_list(:)

  real(dp) :: rn, res1, res2, res3, res4, res5, el
  real(dp), allocatable :: MyMatrix(:,:)

  type(matrix) :: A, B, C, D, E

  !**** EXTERNAL ********************************!

  integer, external :: numroc

  !**********************************************!

  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  call blacs_get(-1,0,icontxt)
  call blacs_gridinit(icontxt,'c',1,mpi_size)

  allocate(bs_list(4))
  bs_list(:)=(/8,3,15,4/)
  call ms_scalapack_setup(mpi_size,1,'c',1,bs_list,icontxt)

  m_storage='pddbc'
  m_operation='lap'

  N=15
  M=8

  call m_allocate(A,N,M,m_storage)

  call blacs_gridinfo(icontxt,i,j,k,l)
  N_loc=numroc(N,4,k,0,1)
  M_loc=numroc(M,3,l,0,mpi_size)
  call descinit(MyDesc,N,M,4,3,0,0,icontxt,N_loc,info)
  allocate(MyMatrix(N_loc,M_loc))

  call m_allocate(C,M,N,m_storage)
  call m_allocate(D,M,M,m_storage)
  call m_allocate(E,N,N,m_storage)

  rn=0.1_dp
  do i=1,N
  do j=1,M
    rn=mod(4.2_dp*rn,1.0_dp)
    call m_set_element(A,i,j,rn)
    rn=mod(4.2_dp*rn,1.0_dp)
    call pdelset(MyMatrix,i,j,MyDesc,rn)
  end do
  end do

  call m_register_pdbc(B,MyMatrix,MyDesc)

  call mm_trace(A,B,res1,m_operation)

  if (mpi_rank==0) print('(a,f21.15)'), 'res1 : ', res1

  call mm_multiply(A,'t',B,'n',D,1.0_dp,0.0_dp,m_operation)

  call m_trace(D,res2,m_operation)

  if (mpi_rank==0) print('(a,f21.15)'), 'res2 : ', res2

  call mm_multiply(B,'n',A,'t',E,1.0_dp,0.0_dp,m_operation)

  call m_trace(E,res3,m_operation)

  if (mpi_rank==0) print('(a,f21.15)'), 'res3 : ', res3

  call m_add(A,'t',C,1.0_dp,0.0_dp,m_operation)

  call mm_multiply(C,'n',B,'n',D,1.0_dp,0.0_dp,m_operation)

  call m_trace(D,res4,m_operation)

  if (mpi_rank==0) print('(a,f21.15)'), 'res4 : ', res4

  call m_add(A,'t',C,1.0_dp,0.0_dp,m_operation)

  call mm_multiply(B,'n',C,'n',E,1.0_dp,0.0_dp,m_operation)

  call m_trace(E,res5,m_operation)

  if (mpi_rank==0) print('(a,f21.15)'), 'res5 : ', res5

  call m_get_element(E,1,1,el)

  if (mpi_rank==0) print('(a,f21.15)'), 'E_11 : ', el

  call m_deallocate(E)
  call m_deallocate(D)
  call m_deallocate(C)
  call m_deallocate(B)
  call m_deallocate(A)

  call blacs_gridexit(icontxt)
  call blacs_exit(1)
  call mpi_finalize(mpi_err)
#else
  print('(a,f21.15)'), 'To run this example, compile with ScaLAPACK.'
#endif

end program example3
