program testList
  use pspBLAS

  implicit none

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)
  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!
  integer :: m, n, i, j, k
  type(dList), pointer :: dmainlist
  type(dList), pointer :: delem
  type(dNodeData)      :: dnode
  real(dp) :: thre, alpha, beta, err
  real(dp) :: rpt, cpt
  complex(dp) :: zalpha, zbeta

  type(zList), pointer :: zmainlist
  type(zList), pointer :: zelem
  type(zNodeData)      :: znode

  real(dp), allocatable :: A(:,:), B(:,:), C(:,:), Ctrue(:,:), valA(:), valB(:), valC(:)
  complex(dp), allocatable :: zA(:,:), zB(:,:), zC(:,:), zCtrue(:,:), valzA(:), valzB(:), valzC(:)
  integer, allocatable :: idx1A(:), idx2A(:), idx1B(:), idx2B(:), idx1C(:), idx2C(:)
  integer, allocatable :: idx1zA(:), idx2zA(:), idx1zB(:), idx2zB(:), idx1zC(:), idx2zC(:)
  character(3) :: fmtA, fmtB, fmtC

  !**********************************************!

  ! create, insert, print and destroy a list
  dnode%val = 1.0_dp
  dnode%col_ind = 1
  dnode%row_ind = 2
  call list_create( dmainlist, dnode )
  dnode%val = 2.0_dp
  dnode%col_ind = 2
  dnode%row_ind = 3
  call list_insert( dmainlist, dnode )
  dnode%val = 3.0_dp
  dnode%col_ind = 3
  dnode%row_ind = 4
  call list_insert( dmainlist, dnode )

  call psp_list_print( 'real list', dmainlist )
  call list_destroy(dmainlist)

  znode%val = cmplx_1
  znode%col_ind = 1
  znode%row_ind = 2
  call list_create( zmainlist, znode )
  znode%val = cmplx_i
  znode%col_ind = 2
  znode%row_ind = 3
  call list_insert( zmainlist, znode )
  znode%val = cmplx_0
  znode%col_ind = 3
  znode%row_ind = 4
  call list_insert( zmainlist, znode )

  call psp_list_print( 'complex list', zmainlist )
  call list_destroy(zmainlist)

  print *, 'allocate matrices'
  ! create sparse matrices
  m = 50
  n = 50
  k = 50
  fmtA = 'csc' ! change parameters to test the code!
  fmtB = 'csc'
  fmtC='csc'
  thre = 0.9_dp
  alpha=1.5_dp
  beta=1.2_dp
  zalpha=(1.2_dp,1.0_dp)
  zbeta=(1.0_dp,2.5_dp)

  allocate(A(m,n))
  allocate(B(m,n))
  allocate(C(m,n))
  allocate(zA(m,n))
  allocate(zB(m,n))
  allocate(zC(m,n))
  allocate(Ctrue(m,n))
  allocate(zCtrue(m,n))
  Ctrue=0.0_dp
  zCtrue=cmplx_0
  C=0.0_dp
  zC=cmplx_0

  print *, '*******************************************************'
  print *, '                test matrix sum                        '
  print *, '*******************************************************'

  !print *, 'generate random dense matrices'
  ! generate random matrices
  call init_random_seed()
  !call random_seed()
  call RANDOM_NUMBER(A)
  call RANDOM_NUMBER(B)
  do j=1,n
     do i=1,m
        call RANDOM_NUMBER(rpt)
        call RANDOM_NUMBER(cpt)
        zA(i,j)=CMPLX(rpt,cpt)
        call RANDOM_NUMBER(rpt)
        call RANDOM_NUMBER(cpt)
        zB(i,j)=CMPLX(rpt,cpt)
     end do
  end do

  !print *, 'thresholding'

  do j=1,n
     do i=1,m
        if (abs(A(i,j))<thre) A(i,j) = 0.0_dp
        if (abs(B(i,j))<thre) B(i,j) = 0.0_dp
        if (abs(zA(i,j))<thre) zA(i,j) = cmplx_0
        if (abs(zB(i,j))<thre) zB(i,j) = cmplx_0

        Ctrue(i,j)=alpha*A(i,j)+beta*B(i,j)
        zCtrue(i,j)=zalpha*zA(i,j)+zbeta*zB(i,j)
     end do
  end do

  !print *, 'begin dense to sparse'
  call psp_sst_den2sp_m(A,idx1A,idx2A,valA,fmtA,thre)
  call psp_sst_den2sp_m(B,idx1B,idx2B,valB,fmtB,thre)
  call psp_sst_den2sp_m(zA,idx1zA,idx2zA,valzA,fmtA,thre)
  call psp_sst_den2sp_m(zB,idx1zB,idx2zB,valzB,fmtB,thre)

  !print *, 'begin spm+spm'
  !print *, 'begin real'
  call psp_sst_sum_spmspm(m,n,alpha,idx1A,idx2A,valA,fmtA,beta,&
       idx1B,idx2B,valB,fmtB,idx1C,idx2C,valC,fmtC,size(idx1A),size(idx1B))
  !print *, 'begin complex'
  call psp_sst_sum_spmspm(m,n,zalpha,idx1zA,idx2zA,valzA,fmtA,zbeta,&
       idx1zB,idx2zB,valzB,fmtB,idx1zC,idx2zC,valzC,fmtC,size(idx1zA),size(idx1zB))

  !print *, 'begin sparse to dens'
  call psp_sst_sp2den_m(m,n,idx1C,idx2C,valC,fmtC,C)
  call psp_sst_sp2den_m(m,n,idx1zC,idx2zC,valzC,fmtC,zC)

  print *, 'check results'
  err=MAXVAL(abs(Ctrue-C))
  print *, 'real case              ', err
  err=MAXVAL(abs(zCtrue-zC))
  print *, 'complex case           ', err


  print *, '*******************************************************'
  print *, '                test matrix product                    '
  print *, '*******************************************************'

  if (.true.) then
     call dgemm('n','n',m,n,k,alpha,A,m,B,m,beta,Ctrue,m)
     call zgemm('n','n',m,n,k,zalpha,zA,m,zB,m,zbeta,zCtrue,m)

     print *, 'begin spmspm'
     !print *, 'real'
     call psp_sst_gespmspm(m,n,k,'n','n',alpha,idx1A,idx2A,valA, &
          idx1B,idx2B,valB,idx1C,idx2C,valC,beta)
     !print *, 'complex'
     call psp_sst_gespmspm(m,n,k,'n','n',zalpha,idx1zA,idx2zA,valzA, &
          idx1zB,idx2zB,valzB,idx1zC,idx2zC,valzC,zbeta)

     !print *, 'begin sparse to dens'
     call psp_sst_sp2den_m(m,n,idx1C,idx2C,valC,fmtC,C)
     call psp_sst_sp2den_m(m,n,idx1zC,idx2zC,valzC,fmtC,zC)

     print *, 'check results'
     err=MAXVAL(abs(Ctrue-C))
     print *, 'real case              ', err
     err=MAXVAL(abs(zCtrue-zC))
     print *, 'complex case           ', err

     if (.false.) then
        print *, zA(1:5,1:5)
        print *, ' '
        print *, zB(1:5,1:5)
        print *, ' '
        print *, zCtrue(1:5,1:5)
        print *, ' '
        print *, zC(1:5,1:5)
     end if

  end if

end program testList
