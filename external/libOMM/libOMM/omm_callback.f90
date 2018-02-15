subroutine omm_callback(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,&
                        m_storage,m_operation)
  use omm_ops
  use MatrixSwitch
  use omm_rand

  implicit none

  !**** INPUT ***********************************!

  character(5), intent(in) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3), intent(in) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use

  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning

  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** CALLBACK ********************************!

  interface
    subroutine H(C,HC)
      use MatrixSwitch
      implicit none
      type(matrix), intent(in) :: C ! WF coeffs. matrix
      type(matrix), intent(inout) :: HC ! work matrix
    end subroutine H
  end interface

  !**** INOUT ***********************************!

  type(matrix), intent(inout) :: S ! overlap matrix (or its Cholesky-factorized upper triangular matrix) in orbital basis (m x m)
  type(matrix), intent(inout) :: D_min ! density (or energy-weighted density) matrix in orbital basis (m x m)
  type(matrix), intent(inout) :: C_min ! WF coeffs. in orbital basis (n x m)
  type(matrix), intent(inout) :: T ! kinetic energy matrix in orbital basis (m x m)

  !**** LOCAL ***********************************!

  logical :: no_S ! S matrix not present?
  logical :: use_Cholesky ! use Cholesky factorization?
  logical :: S_is_U ! S matrix already Cholesky-factorized?
  logical :: use_precon ! use preconditioner?
  logical :: conv ! CG minimization converged?
  logical :: ls_conv ! line search converged?
  logical :: ls_fail ! line search failed?
  logical, save :: log_start=.false.
  logical, allocatable, save :: first_call(:) ! first call for this value of ip?

  integer :: i
  integer :: j
  integer :: k
  integer :: seed
  integer :: icg ! CG step num.
  integer :: n_step_max=100 ! max. num. steps for CG minimization
  integer, save :: log_unit=666

  real(dp) :: rn(2)
  real(dp) :: el
  real(dp) :: e_diff ! relative energy difference respect to previous CG step
  real(dp) :: TrQS ! trace of (2*IW-SW)
  real(dp) :: cg_tol_internal ! convergence tolerance of CG minimization
  real(dp) :: lambda ! CG step size
  real(dp) :: lambda_d ! lambda denominator
  real(dp) :: lambda_n ! lambda numerator
  real(dp) :: e_min_old ! OMM functional energy at previous CG step
  real(dp) :: coeff(0:4) ! coeffs. of the quartic equation
  real(dp), allocatable, save :: x_min(:) ! line search position of minimum

  type(matrix) :: SWd ! G^T*S*C matrix (n x n)
  type(matrix) :: SWdd ! G^T*S*G matrix (n x n)
  type(matrix) :: G ! gradient (n x m)
  type(matrix) :: G_p ! gradient (n x m) at previous CG step
  type(matrix) :: PG ! preconditioned gradient (n x m)
  type(matrix) :: PG_p ! preconditioned gradient (n x m) at previous CG step
  type(matrix) :: D ! conjugate search direction (n x m)
  type(matrix) :: HC ! H*C matrix (n x m)
  type(matrix) :: HG ! H*G matrix (n x m)
  type(matrix) :: work1 ! work matrix (n x n)
  type(matrix), allocatable, save :: HW(:) ! Hamiltonian matrix in WF basis (n x n) for each value of ip
  type(matrix), allocatable, save :: HWd(:) ! G^T*H*C matrix (n x n) for each value of ip
  type(matrix), allocatable, save :: HWdd(:) ! G^T*H*G matrix (n x n) for each value of ip
  type(matrix), allocatable, save :: SW(:) ! overlap matrix in WF basis (n x n) for each value of ip
  type(matrix), allocatable, save :: SC(:) ! S*C matrix (n x m) for each value of ip
  type(matrix), allocatable, save :: SG(:) ! S*G matrix (n x m) for each value of ip
  type(matrix), allocatable, save :: C_Chl(:) ! Cholesky-transformed WF coeffs. in orbital basis (n x m) for each value of ip
  type(matrix), allocatable, save :: P(:) ! preconditioner matrix (m x m) for each value of ip
  type(matrix), allocatable, save :: work2(:) ! work matrix (n x m) for each value of ip
  type(matrix), allocatable, save :: QW(:) ! 2*IW-SW matrix (n x n) for each value of ip
  type(matrix), allocatable, save :: CD(:) ! C*D matrix (n x m) for each value of ip

  !**********************************************!

  if (log_start) then
    open(unit=log_unit,file='libOMM.log',position='append')
  else
    open(unit=log_unit,file='libOMM.log')
    log_start=.true.
  end if

  if (S%is_initialized) then
    no_S=.false.
    if (flavour==0) then
      use_Cholesky=.false.
      S_is_U=.false.
      use_precon=.false.
    !else if (flavour==1) then
    !  use_Cholesky=.true.
    !  S_is_U=.false.
    !  use_precon=.false.
    !else if (flavour==2) then
    !  use_Cholesky=.true.
    !  S_is_U=.true.
    !  use_precon=.false.
    else if (flavour==3) then
      use_Cholesky=.false.
      S_is_U=.false.
      use_precon=.true.
    else
      call die('omm: illegal value of flavour')
    end if
  else
    no_S=.true.
    use_Cholesky=.false.
    S_is_U=.false.
    use_precon=.false.
  end if

  if (cg_tol<0.0_dp) then
    cg_tol_internal=1.0d-9
  else
    cg_tol_internal=cg_tol
  end if

  if (calc_ED) then

    if ((.not. allocated(first_call)) .or. first_call(ip)) call die('omm: cannot calculate ED matrix before D matrix')

    ! calculate the energy-density matrix: ED=C*[(2*I-SW)*(HW+eta*SW)]*C^T
    call m_add(HWd(ip),'c',HW(ip),x_min(ip),1.0_dp,m_operation)
    call m_add(HWd(ip),'n',HW(ip),x_min(ip),1.0_dp,m_operation)
    call m_add(HWdd(ip),'n',HW(ip),x_min(ip)**2,1.0_dp,m_operation)
    if (eta/=0.0_dp) call m_add(SW(ip),'n',HW(ip),eta,1.0_dp,m_operation)
    call calc_A2(HW(ip),C_min,CD(ip),D_min,work2(ip),m_operation)
    if ((.not. dealloc) .and. (eta/=0.0_dp)) call m_add(SW(ip),'n',HW(ip),-eta,1.0_dp,m_operation)

    if (dealloc) then
      if (allocated(CD)) then
        do i=1,np
          if (CD(i)%is_initialized) call m_deallocate(CD(i))
        end do
        deallocate(CD)
      end if
      if (allocated(QW)) then
        do i=1,np
          if (QW(i)%is_initialized) call m_deallocate(QW(i))
        end do
        deallocate(QW)
      end if

      if (allocated(work2)) then
        do i=1,np
          if (work2(i)%is_initialized) call m_deallocate(work2(i))
        end do
        deallocate(work2)
      end if
      if (allocated(P)) then
        do i=1,np
          if (P(i)%is_initialized) call m_deallocate(P(i))
        end do
        deallocate(P)
      end if
      if (allocated(C_Chl)) then
        do i=1,np
          if (C_Chl(i)%is_initialized) call m_deallocate(C_Chl(i))
        end do
        deallocate(C_Chl)
      end if
      if (allocated(SG)) then
        do i=1,np
          if (SG(i)%is_initialized) call m_deallocate(SG(i))
        end do
        deallocate(SG)
      end if
      if (allocated(SC)) then
        do i=1,np
          if (SC(i)%is_initialized) call m_deallocate(SC(i))
        end do
        deallocate(SC)
      end if
      if (allocated(SW)) then
        do i=1,np
          if (SW(i)%is_initialized) call m_deallocate(SW(i))
        end do
        deallocate(SW)
      end if
      if (allocated(HWdd)) then
        do i=1,np
          if (HWdd(i)%is_initialized) call m_deallocate(HWdd(i))
        end do
        deallocate(HWdd)
      end if
      if (allocated(HWd)) then
        do i=1,np
          if (HWd(i)%is_initialized) call m_deallocate(HWd(i))
        end do
        deallocate(HWd)
      end if
      if (allocated(HW)) then
        do i=1,np
          if (HW(i)%is_initialized) call m_deallocate(HW(i))
        end do
        deallocate(HW)
      end if

      if (allocated(x_min)) deallocate(x_min)
      if (allocated(first_call)) deallocate(first_call)
    end if

    close(log_unit)
    return

  end if

  if (.not. allocated(first_call)) then
    allocate(first_call(np))
    first_call=.true.
  end if
  if (.not. allocated(x_min)) allocate(x_min(np))

  if (.not. allocated(HW)) allocate(HW(np))
  if (.not. allocated(HWd)) allocate(HWd(np))
  if (.not. allocated(HWdd)) allocate(HWdd(np))
  if (.not. allocated(SW)) allocate(SW(np))
  if ((.not. use_Cholesky) .or. (.not. no_S)) then
    if (.not. allocated(SC)) allocate(SC(np))
    if (.not. allocated(SG)) allocate(SG(np))
  end if
  if (use_Cholesky) then
    if (.not. allocated(C_Chl)) allocate(C_Chl(np))
  end if
  if (use_precon) then
    if (.not. allocated(P)) allocate(P(np))
  end if
  if (.not. allocated(work2)) allocate(work2(np))

  if (.not. HW(ip)%is_initialized) call m_allocate(HW(ip),n,n,m_storage)
  if (.not. HWd(ip)%is_initialized) call m_allocate(HWd(ip),n,n,m_storage)
  if (.not. HWdd(ip)%is_initialized) call m_allocate(HWdd(ip),n,n,m_storage)
  if (.not. SW(ip)%is_initialized) call m_allocate(SW(ip),n,n,m_storage)
  if ((.not. use_Cholesky) .or. (.not. no_S)) then
    if (.not. SC(ip)%is_initialized) call m_allocate(SC(ip),n,m,m_storage)
    if (.not. SG(ip)%is_initialized) call m_allocate(SG(ip),n,m,m_storage)
  end if
  if (use_Cholesky) then
    if (.not. C_Chl(ip)%is_initialized) call m_allocate(C_Chl(ip),n,m,m_storage)
  end if
  if (use_precon) then
    if (.not. P(ip)%is_initialized) call m_allocate(P(ip),m,m,m_storage)
  end if
  if (.not. work2(ip)%is_initialized) call m_allocate(work2(ip),n,m,m_storage)

  if (.not. SWd%is_initialized) call m_allocate(SWd,n,n,m_storage)
  if (.not. SWdd%is_initialized) call m_allocate(SWdd,n,n,m_storage)
  if (.not. G%is_initialized) call m_allocate(G,n,m,m_storage)
  if (.not. G_p%is_initialized) call m_allocate(G_p,n,m,m_storage)
  if (use_precon) then
    if (.not. PG%is_initialized) call m_allocate(PG,n,m,m_storage)
    if (.not. PG_p%is_initialized) call m_allocate(PG_p,n,m,m_storage)
  end if
  if (.not. D%is_initialized) call m_allocate(D,n,m,m_storage)
  if (.not. HC%is_initialized) call m_allocate(HC,n,m,m_storage)
  if (.not. HG%is_initialized) call m_allocate(HG,n,m,m_storage)
  if (.not. work1%is_initialized) call m_allocate(work1,n,n,m_storage)

  !! reduce the generalized problem to standard form using Cholesky factorization
  !if (use_Cholesky) then
  !  if (new_S .and. (.not. S_is_U)) call m_factorize(S,m_operation)
  !  call m_reduce(S,H,m_operation)
  !end if

  ! calculate the preconditioning matrix P=(S+T/tau)^(-1)
  if (use_precon .and. new_S) then
    if (T%is_initialized) then
      call m_add(T,'n',P(ip),1.0_dp,0.0_dp,m_operation)
      call m_add(S,'n',P(ip),1.0_dp,1.0_dp/scale_T,m_operation)
    else
      call m_add(S,'n',P(ip),1.0_dp,0.0_dp,m_operation)
    end if
    call m_inverse(P(ip),m_operation)
  end if

  ! if this is the first SCF step, then we need to initialize the WF coeffs. matrix with random
  ! numbers between -0.5 and 0.5 (normalize at the end to avoid instabilities), unless we are
  ! reading them from file
  if (first_call(ip) .and. (.not. init_C)) then
    seed = omm_rand_seed()
    do i=1,m
      do j=1,n
        do k = 1,2
          call omm_bsd_lcg(seed, rn(k))
        end do
        !cmplx_el=cmplx(sign(0.5_dp*rn(1),rn(2)-0.5_dp),&
        !               sign(0.5_dp*rn(3),rn(4)-0.5_dp),dp)
        el=sign(0.5_dp*rn(1),rn(2)-0.5_dp)
        call m_set_element(C_min,j,i,el,0.0_dp,m_operation)
      end do
    end do
    call m_scale(C_min,1.0d-2/sqrt(real(m,dp)),m_operation)
  end if
  if (use_Cholesky .and. (init_C .or. new_S)) then
    call m_add(C_min,'n',C_Chl(ip),1.0_dp,0.0_dp,m_operation)
    !if ((.not. first_call(ip)) .or. init_C) call m_transform(S,C_Chl(ip),m_operation)
    call m_transform(S,C_Chl(ip),m_operation)
  end if

  ! first we calculate the energy and gradient for our initial guess, with the following steps:
  ! -calculate the hamiltonian in WF basis: HW=C^T*H*C
  if (use_Cholesky) then
    call calc_HW_callback(H,C_Chl(ip),HW(ip),HC,m_operation)
  else
    call calc_HW_callback(H,C_min,HW(ip),HC,m_operation)
  end if
  ! -calculate the overlap matrix in WF basis: SW=C^T*S*C
  if (use_Cholesky) then
    if (first_call(ip) .or. init_C .or. new_S) call mm_multiply(C_Chl(ip),'n',C_Chl(ip),'c',SW(ip),1.0_dp,0.0_dp,m_operation)
  else if (no_S) then
    if (first_call(ip) .or. init_C) call mm_multiply(C_min,'n',C_min,'c',SW(ip),1.0_dp,0.0_dp,m_operation)
  else
    if (first_call(ip) .or. init_C .or. new_S) then
      call calc_AW(S,C_min,SW(ip),SC(ip),m_operation)
    else
      call m_add(SG(ip),'n',SC(ip),x_min(ip),1.0_dp,m_operation)
    end if
  end if
  ! -calculate the gradient: G=2*(2*H*C-S*C*HW-H*C*SW)
  !  (note that we *reuse* H*C and S*C contained in HC and SC from the previous calls to calc_AW)
  if (use_Cholesky) then
    call calc_G(HW(ip),SW(ip),G,HC,C_Chl(ip),m_operation)
  else if (no_S) then
    call calc_G(HW(ip),SW(ip),G,HC,C_min,m_operation)
  else
    call calc_G(HW(ip),SW(ip),G,HC,SC(ip),m_operation)
  end if
  ! -calculate the preconditioned gradient by premultiplying G by P
  if (use_precon) call mm_multiply(G,'n',P(ip),'n',PG,1.0_dp,0.0_dp,m_operation)
  ! -calculate the additional matrices:
  !  HWd=G^T*H*C
  !  SWd=G^T*S*C
  !  HWdd=G^T*H*G
  !  SWdd=G^T*S*G
  !  (again, H*C has already been calculated, although H*G has not)
  !  and, finally, the coeffs. of the quartic line search equation in the direction g
  !  (the energy at C is given by the zeroth-order coeff. alpha(0))
  if (use_precon) then
    call mm_multiply(HC,'n',PG,'c',HWd(ip),1.0_dp,0.0_dp,m_operation)
    call mm_multiply(SC(ip),'n',PG,'c',SWd,1.0_dp,0.0_dp,m_operation)
  else
    call mm_multiply(HC,'n',G,'c',HWd(ip),1.0_dp,0.0_dp,m_operation)
    if (use_Cholesky) then
      call mm_multiply(C_Chl(ip),'n',G,'c',SWd,1.0_dp,0.0_dp,m_operation)
    else if (no_S) then
      call mm_multiply(C_min,'n',G,'c',SWd,1.0_dp,0.0_dp,m_operation)
    else
      call mm_multiply(SC(ip),'n',G,'c',SWd,1.0_dp,0.0_dp,m_operation)
    end if
  end if
  if (use_precon) then
    call calc_HW_callback(H,PG,HWdd(ip),HG,m_operation)
    call calc_AW(S,PG,SWdd,SG(ip),m_operation)
  else
    call calc_HW_callback(H,G,HWdd(ip),HG,m_operation)
    if (use_Cholesky .or. no_S) then
      call mm_multiply(G,'n',G,'c',SWdd,1.0_dp,0.0_dp,m_operation)
    else
      call calc_AW(S,G,SWdd,SG(ip),m_operation)
    end if
  end if
  call calc_coeff(HW(ip),SW(ip),HWd(ip),SWd,HWdd(ip),SWdd,work1,coeff,m_operation)
  e_min=coeff(0)

  ! this is the main loop of the CG algorithm. We perform a series of line minimizations, with the
  ! gradient G at each new step being modified to obtain the search direction D
  if (ms_mpi_rank==0) then
    write(log_unit,'(a)') '+---------------------------------------------+'
    if (use_Cholesky) then
      write(log_unit,'(a)') '| libOMM (Cholesky factorization)             |'
    else if (use_precon) then
      write(log_unit,'(a)') '| libOMM (preconditioning)                    |'
    else
      write(log_unit,'(a)') '| libOMM                                      |'
    end if
    write(log_unit,'(a)') '+---------------------------------------------+'
    if (long_out) write(log_unit,'(a)') '|             e_min            e_diff         |'
  end if
  icg=0
  do i=1,n_step_max
    lambda=0.0_dp
    do j=1,m*n-1
      if (use_precon) then
        call m_add(PG,'n',D,1.0_dp,lambda)
      else
        call m_add(G,'n',D,1.0_dp,lambda)
      end if
      call m_add(G,'n',G_p,1.0_dp,0.0_dp)
      if (use_precon) call m_add(PG,'n',PG_p,1.0_dp,0.0_dp)
      e_min_old=e_min
      ! if this is not the first CG step, we have to recalculate HWd, SWd, HWdd, SWdd, and the coeffs.
      if (icg>0) then
        call mm_multiply(HC,'n',D,'c',HWd(ip),1.0_dp,0.0_dp,m_operation)
        if (use_Cholesky) then
          call mm_multiply(C_Chl(ip),'n',D,'c',SWd,1.0_dp,0.0_dp,m_operation)
        else if (no_S) then
          call mm_multiply(C_min,'n',D,'c',SWd,1.0_dp,0.0_dp,m_operation)
        else
          call mm_multiply(SC(ip),'n',D,'c',SWd,1.0_dp,0.0_dp,m_operation)
        end if
        call calc_HW_callback(H,D,HWdd(ip),HG,m_operation)
        if (use_Cholesky .or. no_S) then
          call mm_multiply(D,'n',D,'c',SWdd,1.0_dp,0.0_dp,m_operation)
        else
          call calc_AW(S,D,SWdd,SG(ip),m_operation)
        end if
        call calc_coeff(HW(ip),SW(ip),HWd(ip),SWd,HWdd(ip),SWdd,work1,coeff,m_operation)
      end if
      ! using the coeffs. calculated anlytically, we can find the minimum of the functional in the
      ! search direction, and calculate the energy at that minimum
      call omm_solve_quartic(coeff(0:4),x_min(ip),ls_fail)
      ! in certain regions of the coeffs. space the line search gives no minimum--this occurs when there
      ! are positive eigenvalues in the eigenspecturm which are significantly occupied by our coeffs.
      ! matrix; the only known cure, unfortunately, is to scale down the entire matrix, thus returning to
      ! a safe region of the coeffs. space.
      if (ls_fail) then
        if (ms_mpi_rank==0) write(log_unit,'(a)') '| WARNING: Rescaling coefficients!            |'
        e_min=3.0*e_min
        if (use_Cholesky) then
          call m_scale(C_Chl(ip),0.5_dp,m_operation)
        else
          call m_scale(C_min,0.5_dp,m_operation)
        end if
        ls_conv=.false.
      else
        ! if the line search is successful, move to the minimum
        e_min=coeff(4)*x_min(ip)**4+&
              coeff(3)*x_min(ip)**3+&
              coeff(2)*x_min(ip)**2+&
              coeff(1)*x_min(ip)+&
              coeff(0)
        if (use_Cholesky) then
          call m_add(D,'n',C_Chl(ip),x_min(ip),1.0_dp,m_operation)
        else
          call m_add(D,'n',C_min,x_min(ip),1.0_dp,m_operation)
        end if
        ls_conv=.true.
      end if
      ! recalculate SW at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        call m_scale(SW(ip),0.25_dp,m_operation)
      else
        call m_add(SWd,'c',SW(ip),x_min(ip),1.0_dp,m_operation)
        call m_add(SWd,'n',SW(ip),x_min(ip),1.0_dp,m_operation)
        call m_add(SWdd,'n',SW(ip),x_min(ip)**2,1.0_dp,m_operation)
      end if
      e_diff=2.0_dp*abs((e_min-e_min_old)/(e_min+e_min_old))
      if ((ms_mpi_rank==0) .and. long_out) write(log_unit,'(a,2(1x,i5),2(1x,es15.7e3),1x,a)') '|', i, j, e_min, e_diff, '|'
      icg=icg+1
      if (e_diff<=cg_tol_internal) then
        conv=.true.
        exit
      end if
      ! recalculate HW at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        call m_scale(HW(ip),0.25_dp,m_operation)
      else
        call m_add(HWd(ip),'c',HW(ip),x_min(ip),1.0_dp,m_operation)
        call m_add(HWd(ip),'n',HW(ip),x_min(ip),1.0_dp,m_operation)
        call m_add(HWdd(ip),'n',HW(ip),x_min(ip)**2,1.0_dp,m_operation)
      end if
      ! recalculate G at the minimum (or for the rescaled coeffs.)
      if (ls_fail) then
        call m_scale(HC,0.5_dp,m_operation)
        if ((.not. use_Cholesky) .or. (.not. no_S)) call m_scale(SC(ip),0.5_dp,m_operation)
        call m_add(G_p,'n',G,1.0_dp,0.0_dp,m_operation)
        call m_add(HC,'n',G,1.5_dp,1.0_dp,m_operation)
      else
        call m_add(HG,'n',HC,x_min(ip),1.0_dp,m_operation)
        if (use_Cholesky) then
          call calc_G(HW(ip),SW(ip),G,HC,C_Chl(ip),m_operation)
        else if (no_S) then
          call calc_G(HW(ip),SW(ip),G,HC,C_min,m_operation)
        else
          call m_add(SG(ip),'n',SC(ip),x_min(ip),1.0_dp,m_operation)
          call calc_G(HW(ip),SW(ip),G,HC,SC(ip),m_operation)
        end if
      end if
      if (use_precon) call mm_multiply(G,'n',P(ip),'n',PG,1.0_dp,0.0_dp,m_operation)
      if (ls_conv) then
        call m_add(G,'n',work2(ip),1.0_dp,0.0_dp,m_operation)
        call m_add(G_p,'n',work2(ip),-1.0_dp,1.0_dp,m_operation)
        if (use_precon) then
          call mm_trace(PG,work2(ip),lambda_n,m_operation)
          call mm_trace(PG_p,G_p,lambda_d,m_operation)
        else
          call mm_trace(G,work2(ip),lambda_n,m_operation)
          call mm_trace(G_p,G_p,lambda_d,m_operation)
        end if
        lambda=lambda_n/lambda_d
      else
        exit
      end if
    end do
    if (conv) exit
  end do
  if (i>n_step_max) then
    if (ms_mpi_rank==0) write(log_unit,'(a)') '| WARNING: OMM failed to converge!            |'
  end if
  if ((ms_mpi_rank==0) .and. long_out) write(log_unit,'(a)') '+---------------------------------------------+'

  if (work1%is_initialized) call m_deallocate(work1)
  if (HG%is_initialized) call m_deallocate(HG)
  if (HC%is_initialized) call m_deallocate(HC)
  if (D%is_initialized) call m_deallocate(D)
  if (PG_p%is_initialized) call m_deallocate(PG_p)
  if (PG%is_initialized) call m_deallocate(PG)
  if (G_p%is_initialized) call m_deallocate(G_p)
  if (G%is_initialized) call m_deallocate(G)
  if (SWdd%is_initialized) call m_deallocate(SWdd)
  if (SWd%is_initialized) call m_deallocate(SWd)

  if (.not. allocated(QW)) allocate(QW(np))
  if (.not. allocated(CD)) allocate(CD(np))

  if (.not. QW(ip)%is_initialized) call m_allocate(QW(ip),n,n,m_storage)
  if (.not. CD(ip)%is_initialized) call m_allocate(CD(ip),n,m,m_storage)

  ! calculate the density matrix: D=C*(2*IW-SW)*C^T
  call m_set(QW(ip),'a',0.0_dp,2.0_dp,m_operation)
  call m_add(SW(ip),'n',QW(ip),-1.0_dp,1.0_dp,m_operation)
  if (use_Cholesky) then
    call m_add(C_Chl(ip),'n',C_min,1.0_dp,0.0_dp,m_operation)
    call m_back_transform(S,C_min,m_operation)
  end if
  call calc_A(QW(ip),C_min,D_min,CD(ip),m_operation)

  ! calculate the trace of (2*IW-SW) to make sure we are occupying the right number of eigenstates in our
  ! solution
  call mm_trace(QW(ip),SW(ip),TrQS,m_operation)
  if (ms_mpi_rank==0) then
    write(log_unit,'(a,i5,a)')    '| minim: icg           = ', icg, '                |'
    write(log_unit,'(a,f13.7,a)') '| minim: Tr[(2*I-S)*S] = ', TrQS, '        |'
    write(log_unit,'(a)')       '+---------------------------------------------+'
  end if

  e_min=e_min+TrQS*eta

  if (dealloc) then
      if (allocated(CD)) then
        do i=1,np
          if (CD(i)%is_initialized) call m_deallocate(CD(i))
        end do
        deallocate(CD)
      end if
      if (allocated(QW)) then
        do i=1,np
          if (QW(i)%is_initialized) call m_deallocate(QW(i))
        end do
        deallocate(QW)
      end if

      if (allocated(work2)) then
        do i=1,np
          if (work2(i)%is_initialized) call m_deallocate(work2(i))
        end do
        deallocate(work2)
      end if
      if (allocated(P)) then
        do i=1,np
          if (P(i)%is_initialized) call m_deallocate(P(i))
        end do
        deallocate(P)
      end if
      if (allocated(C_Chl)) then
        do i=1,np
          if (C_Chl(i)%is_initialized) call m_deallocate(C_Chl(i))
        end do
        deallocate(C_Chl)
      end if
      if (allocated(SG)) then
        do i=1,np
          if (SG(i)%is_initialized) call m_deallocate(SG(i))
        end do
        deallocate(SG)
      end if
      if (allocated(SC)) then
        do i=1,np
          if (SC(i)%is_initialized) call m_deallocate(SC(i))
        end do
        deallocate(SC)
      end if
      if (allocated(SW)) then
        do i=1,np
          if (SW(i)%is_initialized) call m_deallocate(SW(i))
        end do
        deallocate(SW)
      end if
      if (allocated(HWdd)) then
        do i=1,np
          if (HWdd(i)%is_initialized) call m_deallocate(HWdd(i))
        end do
        deallocate(HWdd)
      end if
      if (allocated(HWd)) then
        do i=1,np
          if (HWd(i)%is_initialized) call m_deallocate(HWd(i))
        end do
        deallocate(HWd)
      end if
      if (allocated(HW)) then
        do i=1,np
          if (HW(i)%is_initialized) call m_deallocate(HW(i))
        end do
        deallocate(HW)
      end if

      if (allocated(x_min)) deallocate(x_min)
      if (allocated(first_call)) deallocate(first_call)
  else
    first_call(ip)=.false.
  end if

  close(log_unit)

end subroutine omm_callback
