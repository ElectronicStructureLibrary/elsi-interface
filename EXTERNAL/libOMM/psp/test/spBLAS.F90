program spBLAS
  use pspBLAS

  implicit none

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)
  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** VARIABLES *******************************!
  integer :: m, n, i, j, k, cnti, cntj, numa, numb, numza, numzb
  type(dList), pointer :: dmainlist
  type(dList), pointer :: delem
  type(dNodeData)      :: dnode
  real(dp) :: thre, alpha, beta, err
  real(dp) :: rpt, cpt, sol, soltrue
  complex(dp) :: zalpha, zbeta, zsol, zsoltrue

  type(zList), pointer :: zmainlist
  type(zList), pointer :: zelem
  type(zNodeData)      :: znode

  real(dp), allocatable :: A(:,:), B(:,:), C(:,:), Ctrue(:,:), D(:,:), valA(:), valB(:), valC(:), valD(:)
  complex(dp), allocatable :: zA(:,:), zB(:,:), zC(:,:), zCtrue(:,:), zD(:,:), valzA(:), valzB(:), valzC(:), valzD(:)
  integer, allocatable :: idx1A(:), idx2A(:), idx1B(:), idx2B(:), idx1C(:), idx2C(:), idx1D(:), idx2D(:)
  integer, allocatable :: idx1zA(:), idx2zA(:), idx1zB(:), idx2zB(:), idx1zC(:), idx2zC(:), idx1zD(:), idx2zD(:)
  real(dp), allocatable :: veca(:), vecb(:), valveca(:), valvecb(:)
  complex(dp), allocatable :: zveca(:), zvecb(:), valzveca(:), valzvecb(:)
  integer, allocatable :: idxveca(:), idxvecb(:), idxzveca(:), idxzvecb(:)
  character(3) :: fmtA, fmtB, fmtC, fmtD

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
  call init_random_seed()

  ! create sparse matrices
  call RANDOM_NUMBER(rpt)
  m = CEILING(500*rpt)
  call RANDOM_NUMBER(rpt)
  n = CEILING(500*rpt)
  call RANDOM_NUMBER(rpt)
  k = CEILING(500*rpt)

  fmtA = 'csc' ! change parameters to test the code!
  fmtB = 'csc'
  fmtC='csc'
  fmtD='csc'
  thre = 0.5_dp
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
  allocate(D(n,k))
  allocate(zD(n,k))
  allocate(veca(m))
  allocate(vecb(m))
  allocate(zveca(m))
  allocate(zvecb(m))
  allocate(idxveca(m))
  allocate(idxvecb(m))
  allocate(valveca(m))
  allocate(valvecb(m))
  allocate(idxzveca(m))
  allocate(idxzvecb(m))
  allocate(valzveca(m))
  allocate(valzvecb(m))

  Ctrue=0.0_dp
  zCtrue=cmplx_0
  C=0.0_dp
  zC=cmplx_0

  print *, '*******************************************************'
  print *, '                test level 1                           '
  print *, '*******************************************************'

  do i=1,m
     call RANDOM_NUMBER(rpt)
     veca(i)=rpt
     call RANDOM_NUMBER(rpt)
     vecb(i)=rpt
     call RANDOM_NUMBER(rpt)
     call RANDOM_NUMBER(cpt)
     zveca(i)=CMPLX(rpt,cpt)
     call RANDOM_NUMBER(rpt)
     call RANDOM_NUMBER(cpt)
     zvecb(i)=CMPLX(rpt,cpt)
  end do

  numa=0
  numb=0
  numza=0
  numzb=0
  do i=1,m
     if (abs(veca(i))<thre) then
        veca(i) = 0.0_dp
     else
        numa=numa+1
        valveca(numa)=veca(i)
        idxveca(numa)=i
     end if
     if (abs(vecb(i))<thre) then
        vecb(i) = 0.0_dp
     else
        numb=numb+1
        valvecb(numb)=vecb(i)
        idxvecb(numb)=i
     end if
     if (abs(zveca(i))<thre) then
        zveca(i) = cmplx_0
     else
        numza=numza+1
        valzveca(numza)=zveca(i)
        idxzveca(numza)=i
     end if
     if (abs(zvecb(i))<thre) then
        zvecb(i) = cmplx_0
     else
        numzb=numzb+1
        valzvecb(numzb)=zvecb(i)
        idxzvecb(numzb)=i
     end if
  end do

  soltrue=0.0_dp
  zsoltrue=cmplx_0
  sol=0.0_dp
  zsol=cmplx_0
  do i=1,m
     soltrue=soltrue+veca(i)*vecb(i)
     zsoltrue=zsoltrue+zveca(i)*zvecb(i)
  end do

  call psp_sst_DOT(m,idxveca,valveca,idxvecb,valvecb,sol,numa,numb)
  call psp_sst_DOT(m,idxzveca,valzveca,idxzvecb,valzvecb,zsol,numza,numzb)

  print *, 'check results'
  err=(abs(soltrue-sol))
  print *, 'real case              ', err, soltrue, sol
  err=(abs(zsoltrue-zsol))
  print *, 'complex case           ', err, zsoltrue, zsol

  soltrue=0.0_dp
  zsoltrue=cmplx_0
  sol=0.0_dp
  zsol=cmplx_0
  do i=1,m
     soltrue=soltrue+veca(i)*vecb(i)
     zsoltrue=zsoltrue+CONJG(zveca(i))*zvecb(i)
  end do

  call psp_sst_DOTC(m,idxveca,valveca,idxvecb,valvecb,sol,numa,numb)
  call psp_sst_DOTC(m,idxzveca,valzveca,idxzvecb,valzvecb,zsol,numza,numzb)

  print *, 'check results'
  err=(abs(soltrue-sol))
  print *, 'real case              ', err, soltrue, sol
  err=(abs(zsoltrue-zsol))
  print *, 'complex case           ', err, zsoltrue, zsol

  print *, '*******************************************************'
  print *, '                test matrix sum                        '
  print *, '*******************************************************'

  ! generate random matrices
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

  print *, 'begin dense to sparse'
  call psp_sst_den2sp_m(A,idx1A,idx2A,valA,fmtA,thre)
  call psp_sst_den2sp_m(B,idx1B,idx2B,valB,fmtB,thre)
  call psp_sst_den2sp_m(zA,idx1zA,idx2zA,valzA,fmtA,thre)
  call psp_sst_den2sp_m(zB,idx1zB,idx2zB,valzB,fmtB,thre)

  print *, 'begin spm+spm'
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

     ! generate random matrices
     if (allocated(A)) deallocate(A)
     allocate(A(m,k))
     if (allocated(B)) deallocate(B)
     allocate(B(k,n))
     if (allocated(zA)) deallocate(zA)
     allocate(zA(m,k))
     if (allocated(zB)) deallocate(zB)
     allocate(zB(k,n))

     call RANDOM_NUMBER(A)
     call RANDOM_NUMBER(B)

     do j=1,k
        do i=1,m
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zA(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     do j=1,n
        do i=1,k
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zB(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     !print *, 'thresholding'

     do j=1,k
        do i=1,m
           if (abs(A(i,j))<thre) A(i,j) = 0.0_dp
           if (abs(zA(i,j))<thre) zA(i,j) = cmplx_0
        end do
     end do

     do j=1,n
        do i=1,k
           if (abs(B(i,j))<thre) B(i,j) = 0.0_dp
           if (abs(zB(i,j))<thre) zB(i,j) = cmplx_0
        end do
     end do

     !print *, 'begin dense to sparse'
     call psp_sst_den2sp_m(A,idx1A,idx2A,valA,fmtA,thre)
     call psp_sst_den2sp_m(B,idx1B,idx2B,valB,fmtB,thre)
     call psp_sst_den2sp_m(zA,idx1zA,idx2zA,valzA,fmtA,thre)
     call psp_sst_den2sp_m(zB,idx1zB,idx2zB,valzB,fmtB,thre)

     if (allocated(Ctrue)) deallocate(Ctrue)
     allocate(Ctrue(m,n))
     if (allocated(zCtrue)) deallocate(zCtrue)
     allocate(zCtrue(m,n))
     call RANDOM_NUMBER(Ctrue)
     do j=1,n
        do i=1,m
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zCtrue(i,j)=CMPLX(rpt,cpt)
        end do
     end do
     do j=1,n
        do i=1,m
           if (abs(Ctrue(i,j))<thre) Ctrue(i,j) = 0.0_dp
           if (abs(zCtrue(i,j))<thre) zCtrue(i,j) = cmplx_0
        end do
     end do
     call psp_sst_den2sp_m(Ctrue,idx1C,idx2C,valC,fmtC,thre)
     call psp_sst_den2sp_m(zCtrue,idx1zC,idx2zC,valzC,fmtC,thre)

     print *, 'n n'
     call dgemm('n','n',m,n,k,alpha,A,m,B,k,beta,Ctrue,m)
     call zgemm('n','n',m,n,k,zalpha,zA,m,zB,k,zbeta,zCtrue,m)

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


     ! generate random matrices
     if (allocated(A)) deallocate(A)
     allocate(A(m,k))
     if (allocated(B)) deallocate(B)
     allocate(B(n,k))
     if (allocated(zA)) deallocate(zA)
     allocate(zA(m,k))
     if (allocated(zB)) deallocate(zB)
     allocate(zB(n,k))
     call RANDOM_NUMBER(A)
     call RANDOM_NUMBER(B)

     do j=1,k
        do i=1,m
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zA(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     do j=1,k
        do i=1,n
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zB(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     !print *, 'thresholding'

     do j=1,k
        do i=1,m
           if (abs(A(i,j))<thre) A(i,j) = 0.0_dp
           if (abs(zA(i,j))<thre) zA(i,j) = cmplx_0
        end do
     end do

     do j=1,k
        do i=1,n
           if (abs(B(i,j))<thre) B(i,j) = 0.0_dp
           if (abs(zB(i,j))<thre) zB(i,j) = cmplx_0
        end do
     end do

     !print *, 'begin dense to sparse'
     call psp_sst_den2sp_m(A,idx1A,idx2A,valA,fmtA,thre)
     call psp_sst_den2sp_m(B,idx1B,idx2B,valB,fmtB,thre)
     call psp_sst_den2sp_m(zA,idx1zA,idx2zA,valzA,fmtA,thre)
     call psp_sst_den2sp_m(zB,idx1zB,idx2zB,valzB,fmtB,thre)

     if (allocated(Ctrue)) deallocate(Ctrue)
     allocate(Ctrue(m,n))
     if (allocated(zCtrue)) deallocate(zCtrue)
     allocate(zCtrue(m,n))
     call RANDOM_NUMBER(Ctrue)
     do j=1,n
        do i=1,m
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zCtrue(i,j)=CMPLX(rpt,cpt)
        end do
     end do
     do j=1,n
        do i=1,m
           if (abs(Ctrue(i,j))<thre) Ctrue(i,j) = 0.0_dp
           if (abs(zCtrue(i,j))<thre) zCtrue(i,j) = cmplx_0
        end do
     end do
     call psp_sst_den2sp_m(Ctrue,idx1C,idx2C,valC,fmtC,thre)
     call psp_sst_den2sp_m(zCtrue,idx1zC,idx2zC,valzC,fmtC,thre)

     print *, 'n t'
     call dgemm('n','t',m,n,k,alpha,A,m,B,n,beta,Ctrue,m)
     call zgemm('n','t',m,n,k,zalpha,zA,m,zB,n,zbeta,zCtrue,m)

     print *, 'begin spmspm'
     !print *, 'real'
     call psp_sst_gespmspm(m,n,k,'n','t',alpha,idx1A,idx2A,valA, &
          idx1B,idx2B,valB,idx1C,idx2C,valC,beta)
     !print *, 'complex'
     call psp_sst_gespmspm(m,n,k,'n','t',zalpha,idx1zA,idx2zA,valzA, &
          idx1zB,idx2zB,valzB,idx1zC,idx2zC,valzC,zbeta)

     !print *, 'begin sparse to dens'
     call psp_sst_sp2den_m(m,n,idx1C,idx2C,valC,fmtC,C)
     call psp_sst_sp2den_m(m,n,idx1zC,idx2zC,valzC,fmtC,zC)

     print *, 'check results'
     err=MAXVAL(abs(Ctrue-C))
     print *, 'real case              ', err
     err=MAXVAL(abs(zCtrue-zC))
     print *, 'complex case           ', err

     print *, 'n c'
     call dgemm('n','c',m,n,k,alpha,A,m,B,n,beta,Ctrue,m)
     call zgemm('n','c',m,n,k,zalpha,zA,m,zB,n,zbeta,zCtrue,m)

     print *, 'begin spmspm'
     !print *, 'real'
     call psp_sst_gespmspm(m,n,k,'n','c',alpha,idx1A,idx2A,valA, &
          idx1B,idx2B,valB,idx1C,idx2C,valC,beta)
     !print *, 'complex'
     call psp_sst_gespmspm(m,n,k,'n','c',zalpha,idx1zA,idx2zA,valzA, &
          idx1zB,idx2zB,valzB,idx1zC,idx2zC,valzC,zbeta)

     !print *, 'begin sparse to dens'
     call psp_sst_sp2den_m(m,n,idx1C,idx2C,valC,fmtC,C)
     call psp_sst_sp2den_m(m,n,idx1zC,idx2zC,valzC,fmtC,zC)

     print *, 'check results'
     err=MAXVAL(abs(Ctrue-C))
     print *, 'real case              ', err
     err=MAXVAL(abs(zCtrue-zC))
     print *, 'complex case           ', err


     ! generate random matrices
     if (allocated(A)) deallocate(A)
     allocate(A(k,m))
     if (allocated(B)) deallocate(B)
     allocate(B(k,n))
     if (allocated(zA)) deallocate(zA)
     allocate(zA(k,m))
     if (allocated(zB)) deallocate(zB)
     allocate(zB(k,n))
     call RANDOM_NUMBER(A)
     call RANDOM_NUMBER(B)

     do j=1,m
        do i=1,k
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zA(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     do j=1,n
        do i=1,k
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zB(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     !print *, 'thresholding'

     do j=1,m
        do i=1,k
           if (abs(A(i,j))<thre) A(i,j) = 0.0_dp
           if (abs(zA(i,j))<thre) zA(i,j) = cmplx_0
        end do
     end do

     do j=1,n
        do i=1,k
           if (abs(B(i,j))<thre) B(i,j) = 0.0_dp
           if (abs(zB(i,j))<thre) zB(i,j) = cmplx_0
        end do
     end do

     !print *, 'begin dense to sparse'
     call psp_sst_den2sp_m(A,idx1A,idx2A,valA,fmtA,thre)
     call psp_sst_den2sp_m(B,idx1B,idx2B,valB,fmtB,thre)
     call psp_sst_den2sp_m(zA,idx1zA,idx2zA,valzA,fmtA,thre)
     call psp_sst_den2sp_m(zB,idx1zB,idx2zB,valzB,fmtB,thre)

     if (allocated(Ctrue)) deallocate(Ctrue)
     allocate(Ctrue(m,n))
     if (allocated(zCtrue)) deallocate(zCtrue)
     allocate(zCtrue(m,n))
     call RANDOM_NUMBER(Ctrue)
     do j=1,n
        do i=1,m
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zCtrue(i,j)=CMPLX(rpt,cpt)
        end do
     end do
     do j=1,n
        do i=1,m
           if (abs(Ctrue(i,j))<thre) Ctrue(i,j) = 0.0_dp
           if (abs(zCtrue(i,j))<thre) zCtrue(i,j) = cmplx_0
        end do
     end do
     call psp_sst_den2sp_m(Ctrue,idx1C,idx2C,valC,fmtC,thre)
     call psp_sst_den2sp_m(zCtrue,idx1zC,idx2zC,valzC,fmtC,thre)

     print *, 't n'
     call dgemm('t','n',m,n,k,alpha,A,k,B,k,beta,Ctrue,m)
     call zgemm('t','n',m,n,k,zalpha,zA,k,zB,k,zbeta,zCtrue,m)

     print *, 'begin spmspm'
     !print *, 'real'
     call psp_sst_gespmspm(m,n,k,'t','n',alpha,idx1A,idx2A,valA, &
          idx1B,idx2B,valB,idx1C,idx2C,valC,beta)
     !print *, 'complex'
     call psp_sst_gespmspm(m,n,k,'t','n',zalpha,idx1zA,idx2zA,valzA, &
          idx1zB,idx2zB,valzB,idx1zC,idx2zC,valzC,zbeta)

     !print *, 'begin sparse to dens'
     call psp_sst_sp2den_m(m,n,idx1C,idx2C,valC,fmtC,C)
     call psp_sst_sp2den_m(m,n,idx1zC,idx2zC,valzC,fmtC,zC)

     print *, 'check results'
     err=MAXVAL(abs(Ctrue-C))
     print *, 'real case              ', err
     err=MAXVAL(abs(zCtrue-zC))
     print *, 'complex case           ', err

     print *, 'c n'
     call dgemm('c','n',m,n,k,alpha,A,k,B,k,beta,Ctrue,m)
     call zgemm('c','n',m,n,k,zalpha,zA,k,zB,k,zbeta,zCtrue,m)

     !print *, 'real'
     call psp_sst_gespmspm(m,n,k,'c','n',alpha,idx1A,idx2A,valA, &
          idx1B,idx2B,valB,idx1C,idx2C,valC,beta)
     !print *, 'complex'
     call psp_sst_gespmspm(m,n,k,'c','n',zalpha,idx1zA,idx2zA,valzA, &
          idx1zB,idx2zB,valzB,idx1zC,idx2zC,valzC,zbeta)

     !print *, 'begin sparse to dens'
     call psp_sst_sp2den_m(m,n,idx1C,idx2C,valC,fmtC,C)
     call psp_sst_sp2den_m(m,n,idx1zC,idx2zC,valzC,fmtC,zC)

     print *, 'check results'
     err=MAXVAL(abs(Ctrue-C))
     print *, 'real case              ', err
     err=MAXVAL(abs(zCtrue-zC))
     print *, 'complex case           ', err

     ! generate random matrices
     if (allocated(A)) deallocate(A)
     allocate(A(k,m))
     if (allocated(B)) deallocate(B)
     allocate(B(n,k))
     if (allocated(zA)) deallocate(zA)
     allocate(zA(k,m))
     if (allocated(zB)) deallocate(zB)
     allocate(zB(n,k))
     call RANDOM_NUMBER(A)
     call RANDOM_NUMBER(B)

     do j=1,m
        do i=1,k
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zA(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     do j=1,k
        do i=1,n
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zB(i,j)=CMPLX(rpt,cpt)
        end do
     end do

     !print *, 'thresholding'

     do j=1,m
        do i=1,k
           if (abs(A(i,j))<thre) A(i,j) = 0.0_dp
           if (abs(zA(i,j))<thre) zA(i,j) = cmplx_0
        end do
     end do

     do j=1,k
        do i=1,n
           if (abs(B(i,j))<thre) B(i,j) = 0.0_dp
           if (abs(zB(i,j))<thre) zB(i,j) = cmplx_0
        end do
     end do

     !print *, 'begin dense to sparse'
     call psp_sst_den2sp_m(A,idx1A,idx2A,valA,fmtA,thre)
     call psp_sst_den2sp_m(B,idx1B,idx2B,valB,fmtB,thre)
     call psp_sst_den2sp_m(zA,idx1zA,idx2zA,valzA,fmtA,thre)
     call psp_sst_den2sp_m(zB,idx1zB,idx2zB,valzB,fmtB,thre)

     if (allocated(Ctrue)) deallocate(Ctrue)
     allocate(Ctrue(m,n))
     if (allocated(zCtrue)) deallocate(zCtrue)
     allocate(zCtrue(m,n))
     call RANDOM_NUMBER(Ctrue)
     do j=1,n
        do i=1,m
           call RANDOM_NUMBER(rpt)
           call RANDOM_NUMBER(cpt)
           zCtrue(i,j)=CMPLX(rpt,cpt)
        end do
     end do
     do j=1,n
        do i=1,m
           if (abs(Ctrue(i,j))<thre) Ctrue(i,j) = 0.0_dp
           if (abs(zCtrue(i,j))<thre) zCtrue(i,j) = cmplx_0
        end do
     end do
     call psp_sst_den2sp_m(Ctrue,idx1C,idx2C,valC,fmtC,thre)
     call psp_sst_den2sp_m(zCtrue,idx1zC,idx2zC,valzC,fmtC,thre)

     print *, 't t'
     call dgemm('t','t',m,n,k,alpha,A,k,B,n,beta,Ctrue,m)
     call zgemm('t','t',m,n,k,zalpha,zA,k,zB,n,zbeta,zCtrue,m)

     print *, 'begin spmspm'
     !print *, 'real'
     call psp_sst_gespmspm(m,n,k,'t','t',alpha,idx1A,idx2A,valA, &
          idx1B,idx2B,valB,idx1C,idx2C,valC,beta)
     !print *, 'complex'
     call psp_sst_gespmspm(m,n,k,'t','t',zalpha,idx1zA,idx2zA,valzA, &
          idx1zB,idx2zB,valzB,idx1zC,idx2zC,valzC,zbeta)

     !print *, 'begin sparse to dens'
     call psp_sst_sp2den_m(m,n,idx1C,idx2C,valC,fmtC,C)
     call psp_sst_sp2den_m(m,n,idx1zC,idx2zC,valzC,fmtC,zC)

     print *, 'check results'
     err=MAXVAL(abs(Ctrue-C))
     print *, 'real case              ', err
     err=MAXVAL(abs(zCtrue-zC))
     print *, 'complex case           ', err

     call dgemm('c','c',m,n,k,alpha,A,k,B,n,beta,Ctrue,m)
     call zgemm('c','c',m,n,k,zalpha,zA,k,zB,n,zbeta,zCtrue,m)

     print *, 'c c'
     print *, 'begin spmspm'
     !print *, 'real'
     call psp_sst_gespmspm(m,n,k,'c','c',alpha,idx1A,idx2A,valA, &
          idx1B,idx2B,valB,idx1C,idx2C,valC,beta)
     !print *, 'complex'
     call psp_sst_gespmspm(m,n,k,'c','c',zalpha,idx1zA,idx2zA,valzA, &
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

end program spBLAS
