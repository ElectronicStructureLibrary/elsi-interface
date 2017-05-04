subroutine tomato_tb_real_py_wrapper(template_basedir,system_label,&
                                     switch1,frac_occ_in,num_orbs_per_atom_in,&
                                     switch2,num_orbs_in,num_cells_dir,&
                                     switch3,sparsity_in,orb_r_cut_in,&
                                     gamma_point,k_point,&
                                     defect,defect_perturbation,&
                                     H_dval,S_dval,m_storage,&
                                     build_matrix,& 
                                     num_occ_states,&
                                     frac_occ_out,num_orbs_per_atom_out,&
                                     num_orbs_out,sparsity_out,&
                                     orb_r_cut_out,&
                                     n_rows_H_dval, n_cols_H_dval, &
                                     n_rows_S_dval, n_cols_S_dval )

  use elsi_omm, only : tomato_tb_real 
  implicit none

  ! Array dimensions hidden by f2py
  integer, intent(in)    :: n_rows_H_dval
  integer, intent(in)    :: n_cols_H_dval
  integer, intent(in)    :: n_rows_S_dval
  integer, intent(in)    :: n_cols_S_dval
  ! f2py arguments
  !**** INPUT ***********************************!

  character(5), intent(in) :: m_storage
  character(*), intent(in) :: template_basedir
  character(*), intent(in) :: system_label

  logical, intent(in) :: switch1
  logical, intent(in) :: switch2
  logical, intent(in) :: switch3
  logical, intent(in) :: gamma_point
  logical, intent(in) :: build_matrix
  logical, intent(in) :: defect

  real*8, intent(in) :: defect_perturbation

  integer, intent(in) :: num_orbs_per_atom_in
  integer, intent(in) :: num_orbs_in

  real*8,  intent(in) :: frac_occ_in
  real*8,  intent(in) :: sparsity_in
  real*8,  intent(in) :: orb_r_cut_in

  !**** OUTPUT **********************************!

  integer, intent(out) :: num_occ_states

  integer, intent(out) :: num_orbs_per_atom_out
  integer, intent(out) :: num_orbs_out

  real*8,  intent(out) :: frac_occ_out
  real*8,  intent(out) :: sparsity_out
  real*8,  intent(out) :: orb_r_cut_out

  !**** INOUT ***********************************!

  integer, intent(inout) :: num_cells_dir(3)
  real*8, intent(inout) :: k_point(3)

  real*8, intent(inout) :: H_dval(n_rows_H_dval, n_cols_H_dval)
  real*8, intent(inout) :: S_dval(n_rows_S_dval, n_cols_S_dval)

  num_orbs_per_atom_out = num_orbs_per_atom_in
  num_orbs_out          = num_orbs_in
  frac_occ_out          = frac_occ_in
  sparsity_out          = sparsity_in
  orb_r_cut_out         = orb_r_cut_in

  call tomato_TB_real(template_basedir,system_label,&
                      switch1,frac_occ_out,num_orbs_per_atom_out,&
                      switch2,num_orbs_out,num_cells_dir,&
                      switch3,sparsity_out,orb_r_cut_out,&
                      num_occ_states,&
                      gamma_point,k_point,&
                      defect,defect_perturbation,&
                      n_rows_H_dval,n_cols_H_dval,H_dval,&
                      n_rows_S_dval,n_cols_S_dval,S_dval,&
                      m_storage, build_matrix)

end subroutine tomato_TB_real_py_wrapper

subroutine tomato_tb_get_dims_py_wrapper(template_basedir,system_label,&
                                         switch1,frac_occ_in,num_orbs_per_atom_in,&
                                         switch2,num_orbs_in,num_cells_dir,&
                                         switch3,sparsity_in,orb_r_cut_in,&
                                         gamma_point,k_point,&
                                         defect,defect_perturbation,&
                                         m_storage, build_matrix,& 
                                         num_occ_states,&
                                         frac_occ_out,num_orbs_per_atom_out,&
                                         num_orbs_out,sparsity_out,&
                                         orb_r_cut_out, &
                                         n_rows_H_dval, n_cols_H_dval, &
                                         n_rows_S_dval, n_cols_S_dval )

  use elsi_omm, only : tomato_tb_get_dims 

  implicit none

  ! f2py arguments
  !**** INPUT ***********************************!

  character(5), intent(in) :: m_storage
  character(*), intent(in) :: template_basedir
  character(*), intent(in) :: system_label

  logical, intent(in) :: switch1
  logical, intent(in) :: switch2
  logical, intent(in) :: switch3
  logical, intent(in) :: gamma_point
  logical, intent(in) :: build_matrix
  logical, intent(in) :: defect

  real*8, intent(in) :: defect_perturbation

  integer, intent(in) :: num_orbs_per_atom_in
  integer, intent(in) :: num_orbs_in

  real*8,  intent(in) :: frac_occ_in
  real*8,  intent(in) :: sparsity_in
  real*8,  intent(in) :: orb_r_cut_in

  !**** OUTPUT **********************************!

  integer, intent(out) :: num_occ_states

  integer, intent(out) :: num_orbs_per_atom_out
  integer, intent(out) :: num_orbs_out

  real*8,  intent(out) :: frac_occ_out
  real*8,  intent(out) :: sparsity_out
  real*8,  intent(out) :: orb_r_cut_out

  !**** INOUT ***********************************!

  integer, intent(inout) :: num_cells_dir(3)
  real*8, intent(inout) :: k_point(3)

  integer, intent(out)    :: n_rows_H_dval
  integer, intent(out)    :: n_cols_H_dval
  integer, intent(out)    :: n_rows_S_dval
  integer, intent(out)    :: n_cols_S_dval

  num_orbs_per_atom_out = num_orbs_per_atom_in
  num_orbs_out          = num_orbs_in
  frac_occ_out          = frac_occ_in
  sparsity_out          = sparsity_in
  orb_r_cut_out         = orb_r_cut_in

  call tomato_TB_get_dims(template_basedir,system_label,&
                          switch1,frac_occ_out,num_orbs_per_atom_out,&
                          switch2,num_orbs_out,num_cells_dir,&
                          switch3,sparsity_out,orb_r_cut_out,&
                          num_occ_states,&
                          gamma_point,k_point,&
                          defect,defect_perturbation,&
                          n_rows_H_dval,n_cols_H_dval,&
                          n_rows_S_dval,n_cols_S_dval,&
                          m_storage, build_matrix)

end subroutine tomato_TB_get_dims_py_wrapper
