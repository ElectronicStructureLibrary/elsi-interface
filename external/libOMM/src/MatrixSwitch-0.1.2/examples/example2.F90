!==================================================================================================!
! example2 : serial program with complex matrices                                                  !
!                                                                                                  !
! This example demonstrates how to calculate:                                                      !
!   alpha = tr(A^H*B)                                                                              !
! for two NxM matrices A and B, in five different ways. Each way should return the same result.    !
! They are:                                                                                        !
!   1. performing the matrix product trace res1 := tr(A^H*B) directly as a single operation        !
!   2. performing the multiplication D := A^H*B, and then the trace res2 := tr(D)                  !
!   3. performing E := B*A^H, and then res3 := tr(E)                                               !
!   4. performing the conjugate transpose C := A^H, then D := C*B, and then res4 := tr(D)          !
!   5. performing C := A^H, then E := B*C, and then res5 := tr(E)                                  !
! Finally, as an extra check, the first element of E is printed out.                               !
!                                                                                                  !
! Note the difference in the code in how matrix B is handled compared with the other matrices.     !
! While the others are allocated directly by MatrixSwitch (i.e., all the data is contained within  !
! the type(matrix) variable), B is a wrapper for MyMatrix (i.e., MyMatrix has been registered as B !
! for use with MatrixSwitch, and B contains a pointer to MyMatrix).                                !
!                                                                                                  !
! Sample output:                                                                                   !
!--------------------------------------------------------------------------------------------------!
! res1 :    58.765528072720869 ,     5.629639202296624                                             !
! res2 :    58.765528072720876 ,     5.629639202296625                                             !
! res3 :    58.765528072720876 ,     5.629639202296625                                             !
! res4 :    58.765528072720876 ,     5.629639202296625                                             !
! res5 :    58.765528072720876 ,     5.629639202296625                                             !
! E_11 :     4.521085405229188 ,    -0.360428156738057                                             !
!--------------------------------------------------------------------------------------------------!
!==================================================================================================!
program example2
  use MatrixSwitch

  implicit none

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!

  character(5) :: m_storage
  character(3) :: m_operation

  integer :: N, M, i, j

  real(dp) :: rn, rn2

  complex(dp) :: res1, res2, res3, res4, res5, el
  complex(dp), allocatable :: MyMatrix(:,:)

  type(matrix) :: A, B, C, D, E

  !**********************************************!

  m_storage='szden'
  m_operation='ref' ! try changing to 'lap'

  N=15
  M=8

  call m_allocate(A,N,M,m_storage)

  allocate(MyMatrix(N,M))

  call m_allocate(C,M,N,m_storage)
  call m_allocate(D,M,M,m_storage)
  call m_allocate(E,N,N,m_storage)

  rn2=0.1_dp
  do i=1,N
  do j=1,M
    rn=mod(4.2_dp*rn2,1.0_dp)
    rn2=mod(4.2_dp*rn,1.0_dp)
    el=cmplx(rn,rn2)
    call m_set_element(A,i,j,el)
    rn=mod(4.2_dp*rn2,1.0_dp)
    rn2=mod(4.2_dp*rn,1.0_dp)
    el=cmplx(rn,rn2)
    MyMatrix(i,j)=el
  end do
  end do

  call m_register_sden(B,MyMatrix)

  call mm_trace(A,B,res1,m_operation)

  print('(2(a,f21.15))'), 'res1 : ', real(res1,dp), ' , ', aimag(res1)

  call mm_multiply(A,'c',B,'n',D,cmplx_1,cmplx_0,m_operation)

  call m_trace(D,res2,m_operation)

  print('(2(a,f21.15))'), 'res2 : ', real(res2,dp), ' , ', aimag(res2)

  call mm_multiply(B,'n',A,'c',E,cmplx_1,cmplx_0,m_operation)

  call m_trace(E,res3,m_operation)

  print('(2(a,f21.15))'), 'res3 : ', real(res3,dp), ' , ', aimag(res3)

  call m_add(A,'c',C,cmplx_1,cmplx_0,m_operation)

  call mm_multiply(C,'n',B,'n',D,cmplx_1,cmplx_0,m_operation)

  call m_trace(D,res4,m_operation)

  print('(2(a,f21.15))'), 'res4 : ', real(res4,dp), ' , ', aimag(res4)

  call m_add(A,'c',C,cmplx_1,cmplx_0,m_operation)

  call mm_multiply(B,'n',C,'n',E,cmplx_1,cmplx_0,m_operation)

  call m_trace(E,res5,m_operation)

  print('(2(a,f21.15))'), 'res5 : ', real(res5,dp), ' , ', aimag(res5)

  call m_get_element(E,1,1,el)

  print('(2(a,f21.15))'), 'E_11 : ', real(el,dp), ' , ', aimag(el)

  call m_deallocate(E)
  call m_deallocate(D)
  call m_deallocate(C)
  call m_deallocate(B)
  call m_deallocate(A)

  deallocate(MyMatrix)

end program example2
