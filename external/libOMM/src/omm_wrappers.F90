subroutine omm_sdden_ref(m,n,H_vals,S_present,S_vals,new_S,e_min,D_min_vals,calc_ED,eta,C_min_vals,init_C,T_present,T_vals,scale_T,&
               flavour,np,ip,cg_tol,long_out,dealloc,mpi_rank)
  use omm_params, only : dp
  use MatrixSwitch
#ifdef CBIND
  use iso_c_binding
#endif

  implicit none

  !**** INPUT ***********************************!

#ifdef CBIND
  logical(c_bool), intent(in) :: S_present ! is the S matrix present?
  logical(c_bool), intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical(c_bool), intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical(c_bool), intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical(c_bool), intent(in) :: T_present ! is the T matrix present?
  logical(c_bool), intent(in) :: long_out ! print detailed output?
  logical(c_bool), intent(in) :: dealloc ! deallocate all internal matrices?
#else
  logical, intent(in) :: S_present ! is the S matrix present?
  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: T_present ! is the T matrix present?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?
#endif

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)
  integer, intent(in) :: mpi_rank ! MPI process rank (0 for serial)

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning
  
  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: H_vals(m,m) ! Hamiltonian matrix
  real(dp), intent(inout) :: S_vals(m,m) ! overlap matrix (or its Cholesky-factorized upper triangular matrix)
  real(dp), intent(inout) :: D_min_vals(m,m) ! density (or energy-weighted density) matrix
  real(dp), intent(inout) :: C_min_vals(n,m) ! WF coeffs.
  real(dp), intent(inout) :: T_vals(m,m) ! kinetic energy matrix

  !**** LOCAL ***********************************!

  character(5) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use

#ifdef CBIND
  logical :: new_S_conv, calc_ED_conv, init_C_conv, long_out_conv, dealloc_conv
#endif

  type(matrix) :: H, S, D_min, C_min, T

  !**********************************************!

  m_storage='sdden'
  m_operation='ref'

  call m_register_sden(H,H_vals)
  if (S_present) call m_register_sden(S,S_vals)
  call m_register_sden(D_min,D_min_vals)
  call m_register_sden(C_min,C_min_vals)
  if (T_present) call m_register_sden(T,T_vals)

#ifdef CBIND
  new_S_conv=new_S
  calc_ED_conv=calc_ED
  init_C_conv=init_C
  long_out_conv=long_out
  dealloc_conv=dealloc
  call omm(m,n,H,S,new_S_conv,e_min,D_min,calc_ED_conv,eta,C_min,init_C_conv,T,scale_T,flavour,np,ip,cg_tol,long_out_conv,&
           dealloc_conv,m_storage,m_operation,mpi_rank)
#else
  call omm(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,m_storage,&
           m_operation,mpi_rank)
#endif

  if (T_present) call m_deallocate(T)
  call m_deallocate(C_min)
  call m_deallocate(D_min)
  if (S_present) call m_deallocate(S)
  call m_deallocate(H)

end subroutine omm_sdden_ref

subroutine omm_szden_ref(m,n,H_vals,S_present,S_vals,new_S,e_min,D_min_vals,calc_ED,eta,C_min_vals,init_C,T_present,T_vals,&
               scale_T,flavour,np,ip,cg_tol,long_out,dealloc,mpi_rank)
  use omm_params, only : dp
  use MatrixSwitch
#ifdef CBIND
  use iso_c_binding
#endif

  implicit none

  !**** INPUT ***********************************!

#ifdef CBIND
  logical(c_bool), intent(in) :: S_present ! is the S matrix present?
  logical(c_bool), intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical(c_bool), intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical(c_bool), intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical(c_bool), intent(in) :: T_present ! is the T matrix present?
  logical(c_bool), intent(in) :: long_out ! print detailed output?
  logical(c_bool), intent(in) :: dealloc ! deallocate all internal matrices?
#else
  logical, intent(in) :: S_present ! is the S matrix present?
  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: T_present ! is the T matrix present?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?
#endif

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)
  integer, intent(in) :: mpi_rank ! MPI process rank (0 for serial)

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning
  
  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: H_vals(m,m) ! Hamiltonian matrix
  complex(dp), intent(inout) :: S_vals(m,m) ! overlap matrix (or its Cholesky-factorized upper triangular matrix)
  complex(dp), intent(inout) :: D_min_vals(m,m) ! density (or energy-weighted density) matrix
  complex(dp), intent(inout) :: C_min_vals(n,m) ! WF coeffs.
  complex(dp), intent(inout) :: T_vals(m,m) ! kinetic energy matrix

  !**** LOCAL ***********************************!

  character(5) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use

#ifdef CBIND
  logical :: new_S_conv, calc_ED_conv, init_C_conv, long_out_conv, dealloc_conv
#endif

  type(matrix) :: H, S, D_min, C_min, T

  !**********************************************!

  m_storage='szden'
  m_operation='ref'

  call m_register_sden(H,H_vals)
  if (S_present) call m_register_sden(S,S_vals)
  call m_register_sden(D_min,D_min_vals)
  call m_register_sden(C_min,C_min_vals)
  if (T_present) call m_register_sden(T,T_vals)

#ifdef CBIND
  new_S_conv=new_S
  calc_ED_conv=calc_ED
  init_C_conv=init_C
  long_out_conv=long_out
  dealloc_conv=dealloc
  call omm(m,n,H,S,new_S_conv,e_min,D_min,calc_ED_conv,eta,C_min,init_C_conv,T,scale_T,flavour,np,ip,cg_tol,long_out_conv,&
           dealloc_conv,m_storage,m_operation,mpi_rank)
#else
  call omm(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,m_storage,&
           m_operation,mpi_rank)
#endif

  if (T_present) call m_deallocate(T)
  call m_deallocate(C_min)
  call m_deallocate(D_min)
  if (S_present) call m_deallocate(S)
  call m_deallocate(H)

end subroutine omm_szden_ref

subroutine omm_sdden_lap(m,n,H_vals,S_present,S_vals,new_S,e_min,D_min_vals,calc_ED,eta,C_min_vals,init_C,T_present,T_vals,&
               scale_T,flavour,np,ip,cg_tol,long_out,dealloc,mpi_rank)
  use omm_params, only : dp
  use MatrixSwitch
#ifdef CBIND
  use iso_c_binding
#endif

  implicit none

  !**** INPUT ***********************************!

#ifdef CBIND
  logical(c_bool), intent(in) :: S_present ! is the S matrix present?
  logical(c_bool), intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical(c_bool), intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical(c_bool), intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical(c_bool), intent(in) :: T_present ! is the T matrix present?
  logical(c_bool), intent(in) :: long_out ! print detailed output?
  logical(c_bool), intent(in) :: dealloc ! deallocate all internal matrices?
#else
  logical, intent(in) :: S_present ! is the S matrix present?
  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: T_present ! is the T matrix present?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?
#endif

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)
  integer, intent(in) :: mpi_rank ! MPI process rank (0 for serial)

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning
  
  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** INOUT ***********************************!

  real(dp), intent(inout) :: H_vals(m,m) ! Hamiltonian matrix
  real(dp), intent(inout) :: S_vals(m,m) ! overlap matrix (or its Cholesky-factorized upper triangular matrix)
  real(dp), intent(inout) :: D_min_vals(m,m) ! density (or energy-weighted density) matrix
  real(dp), intent(inout) :: C_min_vals(n,m) ! WF coeffs.
  real(dp), intent(inout) :: T_vals(m,m) ! kinetic energy matrix

  !**** LOCAL ***********************************!

  character(5) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use

#ifdef CBIND
  logical :: new_S_conv, calc_ED_conv, init_C_conv, long_out_conv, dealloc_conv
#endif

  type(matrix) :: H, S, D_min, C_min, T

  !**********************************************!

  m_storage='sdden'
  m_operation='lap'

  call m_register_sden(H,H_vals)
  if (S_present) call m_register_sden(S,S_vals)
  call m_register_sden(D_min,D_min_vals)
  call m_register_sden(C_min,C_min_vals)
  if (T_present) call m_register_sden(T,T_vals)

#ifdef CBIND
  new_S_conv=new_S
  calc_ED_conv=calc_ED
  init_C_conv=init_C
  long_out_conv=long_out
  dealloc_conv=dealloc
  call omm(m,n,H,S,new_S_conv,e_min,D_min,calc_ED_conv,eta,C_min,init_C_conv,T,scale_T,flavour,np,ip,cg_tol,long_out_conv,&
           dealloc_conv,m_storage,m_operation,mpi_rank)
#else
  call omm(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,m_storage,&
           m_operation,mpi_rank)
#endif

  if (T_present) call m_deallocate(T)
  call m_deallocate(C_min)
  call m_deallocate(D_min)
  if (S_present) call m_deallocate(S)
  call m_deallocate(H)

end subroutine omm_sdden_lap

subroutine omm_szden_lap(m,n,H_vals,S_present,S_vals,new_S,e_min,D_min_vals,calc_ED,eta,C_min_vals,init_C,T_present,T_vals,&
               scale_T,flavour,np,ip,cg_tol,long_out,dealloc,mpi_rank)
  use omm_params, only : dp
  use MatrixSwitch
#ifdef CBIND
  use iso_c_binding
#endif

  implicit none

  !**** INPUT ***********************************!

#ifdef CBIND
  logical(c_bool), intent(in) :: S_present ! is the S matrix present?
  logical(c_bool), intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical(c_bool), intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical(c_bool), intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical(c_bool), intent(in) :: T_present ! is the T matrix present?
  logical(c_bool), intent(in) :: long_out ! print detailed output?
  logical(c_bool), intent(in) :: dealloc ! deallocate all internal matrices?
#else
  logical, intent(in) :: S_present ! is the S matrix present?
  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: T_present ! is the T matrix present?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?
#endif

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)
  integer, intent(in) :: mpi_rank ! MPI process rank (0 for serial)

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning
  
  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** INOUT ***********************************!

  complex(dp), intent(inout) :: H_vals(m,m) ! Hamiltonian matrix
  complex(dp), intent(inout) :: S_vals(m,m) ! overlap matrix (or its Cholesky-factorized upper triangular matrix)
  complex(dp), intent(inout) :: D_min_vals(m,m) ! density (or energy-weighted density) matrix
  complex(dp), intent(inout) :: C_min_vals(n,m) ! WF coeffs.
  complex(dp), intent(inout) :: T_vals(m,m) ! kinetic energy matrix

  !**** LOCAL ***********************************!

  character(5) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use

#ifdef CBIND
  logical :: new_S_conv, calc_ED_conv, init_C_conv, long_out_conv, dealloc_conv
#endif

  type(matrix) :: H, S, D_min, C_min, T

  !**********************************************!

  m_storage='szden'
  m_operation='lap'

  call m_register_sden(H,H_vals)
  if (S_present) call m_register_sden(S,S_vals)
  call m_register_sden(D_min,D_min_vals)
  call m_register_sden(C_min,C_min_vals)
  if (T_present) call m_register_sden(T,T_vals)

#ifdef CBIND
  new_S_conv=new_S
  calc_ED_conv=calc_ED
  init_C_conv=init_C
  long_out_conv=long_out
  dealloc_conv=dealloc
  call omm(m,n,H,S,new_S_conv,e_min,D_min,calc_ED_conv,eta,C_min,init_C_conv,T,scale_T,flavour,np,ip,cg_tol,long_out_conv,&
           dealloc_conv,m_storage,m_operation,mpi_rank)
#else
  call omm(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,m_storage,&
           m_operation,mpi_rank)
#endif

  if (T_present) call m_deallocate(T)
  call m_deallocate(C_min)
  call m_deallocate(D_min)
  if (S_present) call m_deallocate(S)
  call m_deallocate(H)

end subroutine omm_szden_lap

#ifdef MPI
subroutine omm_pddbc_lap(m,n,H_dim,H_vals,desc_H,S_present,S_dim,S_vals,desc_S,new_S,e_min,D_min_dim,D_min_vals,desc_D_min,&
               calc_ED,eta,C_min_dim,C_min_vals,desc_C_min,init_C,T_present,T_dim,T_vals,desc_T,scale_T,flavour,np,ip,cg_tol,&
               long_out,dealloc,mpi_rank,mpi_size,nprow,order,bs_def,icontxt)
  use omm_params, only : dp, ms_scalapack_running
  use MatrixSwitch
#ifdef CBIND
  use iso_c_binding
#endif

  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: order

#ifdef CBIND
  logical(c_bool), intent(in) :: S_present ! is the S matrix present?
  logical(c_bool), intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical(c_bool), intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical(c_bool), intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical(c_bool), intent(in) :: T_present ! is the T matrix present?
  logical(c_bool), intent(in) :: long_out ! print detailed output?
  logical(c_bool), intent(in) :: dealloc ! deallocate all internal matrices?
#else
  logical, intent(in) :: S_present ! is the S matrix present?
  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: T_present ! is the T matrix present?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?
#endif

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: H_dim(2) ! local dimensions for H
  integer, intent(in) :: S_dim(2) ! local dimensions for S
  integer, intent(in) :: D_min_dim(2) ! local dimensions for D_min
  integer, intent(in) :: C_min_dim(2) ! local dimensions for C_min
  integer, intent(in) :: T_dim(2) ! local dimensions for T
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)
  integer, intent(in) :: mpi_rank ! MPI process rank (0 for serial)
  integer, intent(in) :: mpi_size ! total number of MPI processes for the processor grid
  integer, intent(in) :: nprow ! number of rows in the processor grid
  integer, intent(in) :: bs_def ! default block size
  integer, intent(in) :: icontxt ! existing BLACS context handle

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning
  
  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** INOUT ***********************************!

  integer, intent(inout) :: desc_H(9) ! BLACS array descriptor for H
  integer, intent(inout) :: desc_S(9) ! BLACS array descriptor for S
  integer, intent(inout) :: desc_D_min(9) ! BLACS array descriptor for D_min
  integer, intent(inout) :: desc_C_min(9) ! BLACS array descriptor for C_min
  integer, intent(inout) :: desc_T(9) ! BLACS array descriptor for T

  real(dp), intent(inout) :: H_vals(H_dim(1),H_dim(2)) ! Hamiltonian matrix
  real(dp), intent(inout) :: S_vals(S_dim(1),S_dim(2)) ! overlap matrix (or its Cholesky-factorized upper triangular matrix)
  real(dp), intent(inout) :: D_min_vals(D_min_dim(1),D_min_dim(2)) ! density (or energy-weighted density) matrix
  real(dp), intent(inout) :: C_min_vals(C_min_dim(1),C_min_dim(2)) ! WF coeffs.
  real(dp), intent(inout) :: T_vals(T_dim(1),T_dim(2)) ! kinetic energy matrix

  !**** LOCAL ***********************************!

  character(5) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use

#ifdef CBIND
  logical :: new_S_conv, calc_ED_conv, init_C_conv, long_out_conv, dealloc_conv
#endif

  type(matrix) :: H, S, D_min, C_min, T

  !**********************************************!

  m_storage='pddbc'
  m_operation='lap'

  if (.not. ms_scalapack_running) then
    call ms_scalapack_setup(mpi_size,nprow,order,bs_def,icontxt=icontxt)
    ms_scalapack_running=.true.
  end if

  call m_register_pdbc(H,H_vals,desc_H)
  if (S_present) call m_register_pdbc(S,S_vals,desc_S)
  call m_register_pdbc(D_min,D_min_vals,desc_D_min)
  call m_register_pdbc(C_min,C_min_vals,desc_C_min)
  if (T_present) call m_register_pdbc(T,T_vals,desc_T)

#ifdef CBIND
  new_S_conv=new_S
  calc_ED_conv=calc_ED
  init_C_conv=init_C
  long_out_conv=long_out
  dealloc_conv=dealloc
  call omm(m,n,H,S,new_S_conv,e_min,D_min,calc_ED_conv,eta,C_min,init_C_conv,T,scale_T,flavour,np,ip,cg_tol,long_out_conv,&
           dealloc_conv,m_storage,m_operation,mpi_rank)
#else
  call omm(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,m_storage,&
           m_operation,mpi_rank)
#endif

  if (T_present) call m_deallocate(T)
  call m_deallocate(C_min)
  call m_deallocate(D_min)
  if (S_present) call m_deallocate(S)
  call m_deallocate(H)

end subroutine omm_pddbc_lap

subroutine omm_pzdbc_lap(m,n,H_dim,H_vals,desc_H,S_present,S_dim,S_vals,desc_S,new_S,e_min,D_min_dim,D_min_vals,desc_D_min,&
               calc_ED,eta,C_min_dim,C_min_vals,desc_C_min,init_C,T_present,T_dim,T_vals,desc_T,scale_T,flavour,np,ip,cg_tol,&
               long_out,dealloc,mpi_rank,mpi_size,nprow,order,bs_def,icontxt)
  use omm_params, only : dp, ms_scalapack_running
  use MatrixSwitch
#ifdef CBIND
  use iso_c_binding
#endif

  implicit none

  !**** INPUT ***********************************!

  character(1), intent(in) :: order

#ifdef CBIND
  logical(c_bool), intent(in) :: S_present ! is the S matrix present?
  logical(c_bool), intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical(c_bool), intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical(c_bool), intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical(c_bool), intent(in) :: T_present ! is the T matrix present?
  logical(c_bool), intent(in) :: long_out ! print detailed output?
  logical(c_bool), intent(in) :: dealloc ! deallocate all internal matrices?
#else
  logical, intent(in) :: S_present ! is the S matrix present?
  logical, intent(in) :: new_S ! is the S matrix new for this value of ip?
  logical, intent(in) :: calc_ED ! calculate the energy-weighted density matrix from the existing WF coeffs.?
  logical, intent(in) :: init_C ! have the WF coeffs. been initialized or modified externally to libOMM?
  logical, intent(in) :: T_present ! is the T matrix present?
  logical, intent(in) :: long_out ! print detailed output?
  logical, intent(in) :: dealloc ! deallocate all internal matrices?
#endif

  integer, intent(in) :: m ! size of basis
  integer, intent(in) :: n ! number of occupied states
  integer, intent(in) :: H_dim(2) ! local dimensions for H
  integer, intent(in) :: S_dim(2) ! local dimensions for S
  integer, intent(in) :: D_min_dim(2) ! local dimensions for D_min
  integer, intent(in) :: C_min_dim(2) ! local dimensions for C_min
  integer, intent(in) :: T_dim(2) ! local dimensions for T
  integer, intent(in) :: flavour ! flavour of the OMM functional:
                                 ! 0 for basic
                                 ! 1 for Cholesky factorization, S provided
                                 ! 2 for Cholesky factorization, U provided
                                 ! 3 for preconditioning, S provided (T optional)
  integer, intent(in) :: np ! (number of spin points)*(number of k points)
  integer, intent(in) :: ip ! spin+k point identifier (from 1 to np)
  integer, intent(in) :: mpi_rank ! MPI process rank (0 for serial)
  integer, intent(in) :: mpi_size ! total number of MPI processes for the processor grid
  integer, intent(in) :: nprow ! number of rows in the processor grid
  integer, intent(in) :: bs_def ! default block size
  integer, intent(in) :: icontxt ! existing BLACS context handle

  real(dp), intent(in) :: eta ! eigenspectrum shift parameter
  real(dp), intent(in) :: cg_tol ! convergence tolerance of CG minimization (if negative, default of 1.0d-9 is used)
  real(dp), intent(in) :: scale_T ! kinetic energy scale for the preconditioning
  
  !**** OUTPUT **********************************!

  real(dp), intent(out) :: e_min ! OMM functional energy (spin degeneracy *not* included)

  !**** INOUT ***********************************!

  integer, intent(inout) :: desc_H(9) ! BLACS array descriptor for H
  integer, intent(inout) :: desc_S(9) ! BLACS array descriptor for S
  integer, intent(inout) :: desc_D_min(9) ! BLACS array descriptor for D_min
  integer, intent(inout) :: desc_C_min(9) ! BLACS array descriptor for C_min
  integer, intent(inout) :: desc_T(9) ! BLACS array descriptor for T

  complex(dp), intent(inout) :: H_vals(H_dim(1),H_dim(2)) ! Hamiltonian matrix
  complex(dp), intent(inout) :: S_vals(S_dim(1),S_dim(2)) ! overlap matrix (or its Cholesky-factorized upper triangular matrix)
  complex(dp), intent(inout) :: D_min_vals(D_min_dim(1),D_min_dim(2)) ! density (or energy-weighted density) matrix
  complex(dp), intent(inout) :: C_min_vals(C_min_dim(1),C_min_dim(2)) ! WF coeffs.
  complex(dp), intent(inout) :: T_vals(T_dim(1),T_dim(2)) ! kinetic energy matrix

  !**** LOCAL ***********************************!

  character(5) :: m_storage ! label identifying the MatrixSwitch storage format
  character(3) :: m_operation ! label identifying the MatrixSwitch implementation of the operations to use

#ifdef CBIND
  logical :: new_S_conv, calc_ED_conv, init_C_conv, long_out_conv, dealloc_conv
#endif

  type(matrix) :: H, S, D_min, C_min, T

  !**********************************************!

  m_storage='pzdbc'
  m_operation='lap'

  if (.not. ms_scalapack_running) then
    call ms_scalapack_setup(mpi_size,nprow,order,bs_def,icontxt=icontxt)
    ms_scalapack_running=.true.
  end if

  call m_register_pdbc(H,H_vals,desc_H)
  if (S_present) call m_register_pdbc(S,S_vals,desc_S)
  call m_register_pdbc(D_min,D_min_vals,desc_D_min)
  call m_register_pdbc(C_min,C_min_vals,desc_C_min)
  if (T_present) call m_register_pdbc(T,T_vals,desc_T)

#ifdef CBIND
  new_S_conv=new_S
  calc_ED_conv=calc_ED
  init_C_conv=init_C
  long_out_conv=long_out
  dealloc_conv=dealloc
  call omm(m,n,H,S,new_S_conv,e_min,D_min,calc_ED_conv,eta,C_min,init_C_conv,T,scale_T,flavour,np,ip,cg_tol,long_out_conv,&
           dealloc_conv,m_storage,m_operation,mpi_rank)
#else
  call omm(m,n,H,S,new_S,e_min,D_min,calc_ED,eta,C_min,init_C,T,scale_T,flavour,np,ip,cg_tol,long_out,dealloc,m_storage,&
           m_operation,mpi_rank)
#endif

  if (T_present) call m_deallocate(T)
  call m_deallocate(C_min)
  call m_deallocate(D_min)
  if (S_present) call m_deallocate(S)
  call m_deallocate(H)

end subroutine omm_pzdbc_lap
#endif
