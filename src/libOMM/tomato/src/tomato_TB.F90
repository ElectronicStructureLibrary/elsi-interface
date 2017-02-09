subroutine tomato_TB(template_basedir,system_label,&
                     switch1,frac_occ,num_orbs_per_atom,&
                     switch2,num_orbs,num_cells_dir,&
                     switch3,sparsity,orb_r_cut,&
                     num_occ_states,&
                     gamma_point,k_point,&
                     defect,defect_perturbation,&
                     H,S,m_storage,&
                     build_matrix)
  use MatrixSwitch
#ifdef MPI
  use MatrixSwitch_ops, only : ms_mpi_comm, ms_mpi_size, ms_mpi_rank
#endif

  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  !**** PARAMS **********************************!

  integer, parameter :: dp=selected_real_kind(15,300)
  integer, parameter :: i64=selected_int_kind(18)

  real(dp), parameter :: Pi=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: twoPi=2.0_dp*Pi

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** TYPES ***********************************!

  type t_basis_subset
     integer :: ref
     integer, allocatable :: index(:)
  end type t_basis_subset

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

  real(dp), intent(in) :: defect_perturbation

  !**** OUTPUT **********************************!

  integer, intent(out) :: num_occ_states

  !**** INOUT ***********************************!

  integer, intent(inout) :: num_orbs_per_atom
  integer, intent(inout) :: num_orbs
  integer, intent(inout) :: num_cells_dir(3)

  real(dp), intent(inout) :: frac_occ
  real(dp), intent(inout) :: sparsity
  real(dp), intent(inout) :: orb_r_cut
  real(dp), intent(inout) :: k_point(3)

  type(matrix), intent(inout) :: H
  type(matrix), intent(inout) :: S

  !**** VARIABLES *******************************!

  character(10) :: c10a, c10b
  character(100) :: template_filename
  character(500) :: line

  logical :: is_subset
  logical :: overlap_present

#ifdef MPI
  integer :: mpi_comm
#endif
  integer :: mpi_err, mpi_size, mpi_rank
  integer :: i, j, k, l, a, b, c, i2, j2, k2, i3, j3, k3
  integer :: template_num_cutoffs
  integer :: template_num_basis_sizes
  integer :: template_index
  integer :: num_periodic_dims
  integer :: num_occ_states_per_atom
  integer :: num_occ_states_per_cell
  integer :: num_atoms_per_cell
  integer :: num_orbs_per_cell
  integer :: num_cells
  integer :: num_atoms
  integer :: num_template
  integer :: orb_r_cut_i
  integer :: nzel
  integer :: seed
  integer, allocatable :: template_cutoffs(:)
  integer, allocatable :: template_basis_sizes(:)
  integer, allocatable :: template_i(:,:)
  integer, allocatable :: cell_index(:,:,:)
  integer, allocatable :: subset_convert(:)

  real(dp) :: kdotT
  real(dp) :: d, d2, d_min, el1, el2, rn
  real(dp), allocatable :: template_d(:,:)

  complex(dp):: cel1, cel2

  type(t_basis_subset), allocatable :: basis_subset(:)

  !**********************************************!

#ifdef MPI
  mpi_comm=ms_mpi_comm
  mpi_size=ms_mpi_size
  mpi_rank=ms_mpi_rank
#else
  mpi_size=1
  mpi_rank=0
#endif

  ! read the template info file
  template_filename=trim(adjustl(template_basedir))//'/'//&
                    trim(adjustl(system_label))//'_template.info'
  if (mpi_rank==0) open(10,file=trim(adjustl(template_filename)))
    if (mpi_rank==0) read(10,*) num_periodic_dims
#ifdef MPI
    call mpi_bcast(num_periodic_dims,1,mpi_int,0,mpi_comm,mpi_err)
#endif
    if (mpi_rank==0) read(10,*) num_atoms_per_cell
#ifdef MPI
    call mpi_bcast(num_atoms_per_cell,1,mpi_int,0,mpi_comm,mpi_err)
#endif
    if (mpi_rank==0) read(10,*) num_occ_states_per_atom
#ifdef MPI
    call mpi_bcast(num_occ_states_per_atom,1,mpi_int,0,mpi_comm,mpi_err)
#endif
    num_occ_states_per_cell=num_occ_states_per_atom*num_atoms_per_cell
    if (mpi_rank==0) read(10,*) template_num_cutoffs
#ifdef MPI
    call mpi_bcast(template_num_cutoffs,1,mpi_int,0,mpi_comm,mpi_err)
#endif
    allocate(template_cutoffs(template_num_cutoffs))
    if (mpi_rank==0) then
      do i=1,template_num_cutoffs
        read(10,*) template_cutoffs(i)
      end do
    end if
#ifdef MPI
    call mpi_bcast(template_cutoffs,template_num_cutoffs,mpi_int,0,mpi_comm,mpi_err)
#endif
    if (mpi_rank==0) read(10,*) template_num_basis_sizes
#ifdef MPI
    call mpi_bcast(template_num_basis_sizes,1,mpi_int,0,mpi_comm,mpi_err)
#endif
    allocate(template_basis_sizes(template_num_basis_sizes))
    allocate(basis_subset(template_num_basis_sizes))
    do i=1,template_num_basis_sizes
      if (mpi_rank==0) read(10,'(a)') line
#ifdef MPI
      call mpi_bcast(line,200,mpi_char,0,mpi_comm,mpi_err)
#endif
      read(line,*,iostat=j) template_basis_sizes(i), basis_subset(i)%ref
      if (j==0) then
        allocate(basis_subset(i)%index(num_atoms_per_cell*template_basis_sizes(i)))
        read(line,*) template_basis_sizes(i), basis_subset(i)%ref, basis_subset(i)%index(1:template_basis_sizes(i))
        do k=2,num_atoms_per_cell
          do l=1,template_basis_sizes(i)
            basis_subset(i)%index((k-1)*template_basis_sizes(i)+l)=&
                 (k-1)*basis_subset(i)%ref+basis_subset(i)%index(l)
          end do
        end do
      else
        read(line,*) template_basis_sizes(i)
        basis_subset(i)%ref=template_basis_sizes(i)
      end if
    end do
#ifdef MPI
  close(10)
#endif

  ! determine the number of orbitals/molecular unit by selecting a basis size
  ! from the template list
  if (switch1) then
    ! try to match the inputted value of frac_occ
    d_min=999999.9_dp
    do i=1,template_num_basis_sizes
      d=abs(frac_occ-real(num_occ_states_per_atom,dp)/real(template_basis_sizes(i),dp))
      if (d<d_min) then
        d_min=d
        num_orbs_per_atom=template_basis_sizes(i)
        template_index=i
      end if
    end do
  else
    ! try to match the inputted value of num_orbs_per_atom 
    j=999999
    l=num_orbs_per_atom
    do i=1,template_num_basis_sizes
      k=abs(l-template_basis_sizes(i))
      if (k<j) then
        j=k
        num_orbs_per_atom=template_basis_sizes(i)
        template_index=i
      end if
    end do
  end if
  num_orbs_per_cell=num_orbs_per_atom*num_atoms_per_cell
  frac_occ=real(num_occ_states_per_atom,dp)/real(num_orbs_per_atom,dp)
  ! determine whether the selected basis size is a subset of a larger basis
  if (template_basis_sizes(template_index)==basis_subset(template_index)%ref) then
    is_subset=.false.
  else
    is_subset=.true.
  end if

  ! determine the total number of orbitals by calculating the supercell to be
  ! used
  if (switch2) then
    ! try to match the inputted value of num_orbs (depending on the
    ! periodicity of the system)
    if (num_periodic_dims<3) then
      num_cells_dir(1)=1
    else
      num_cells_dir(1)=nint((real(num_orbs,dp)/real(num_orbs_per_cell,dp))**(1.0_dp/3.0_dp))
      num_cells_dir(1)=max(num_cells_dir(1),1)
    end if
    if (num_periodic_dims<2) then
      num_cells_dir(2)=1
    else
      num_cells_dir(2)=nint(sqrt((real(num_orbs,dp)/real(num_orbs_per_cell,dp))/num_cells_dir(1)))
      num_cells_dir(2)=max(num_cells_dir(2),1)
    end if
    if (num_periodic_dims<1) then
      num_cells_dir(3)=1
    else
      num_cells_dir(3)=nint((real(num_orbs,dp)/real(num_orbs_per_cell,dp))/(num_cells_dir(1)*num_cells_dir(2)))
      num_cells_dir(3)=max(num_cells_dir(3),1)
    end if
  else
    ! use the inputted value of num_cells_dir (fix exceptions)
    do i=1,3
      if (num_periodic_dims<4-i) then
        num_cells_dir(i)=1
      else
        num_cells_dir(i)=max(num_cells_dir(i),1)
      end if
    end do
  end if
  num_cells=num_cells_dir(1)*&
            num_cells_dir(2)*&
            num_cells_dir(3)
  num_atoms=num_atoms_per_cell*num_cells
  num_orbs=num_orbs_per_cell*num_cells
  num_occ_states=num_occ_states_per_cell*num_cells

  ! determine the sparsity of the system by selecting an orbital cutoff length
  ! from the template list
  if (switch3) then
    ! try to match the inputted value of sparsity with a simple estimation of
    ! the sparsity of the system with different orbital cutoff lengths
    d_min=999999.9_dp
    do i=1,template_num_cutoffs
      if (mpi_rank==0) then
        write(c10a,'(i10)') template_cutoffs(i)
        write(c10b,'(i10)') basis_subset(template_index)%ref
        template_filename=trim(adjustl(template_basedir))//'/'//&
                          trim(adjustl(system_label))//'_'//&
                          trim(adjustl(c10a))//'_'//&
                          trim(adjustl(c10b))//'.template'
        open(10,file=trim(adjustl(template_filename)))
          read(10,'(2x,i10)') num_template
        close(10)
        if (is_subset) num_template=num_template*num_orbs_per_atom**2/basis_subset(template_index)%ref**2
      end if
#ifdef MPI
      call mpi_bcast(num_template,1,mpi_int,0,mpi_comm,mpi_err)
#endif
      d2=1.0_dp-real(num_template,dp)*real(num_cells,dp)/(real(num_orbs,dp)**2)
      d2=max(d2,0.0_dp)
      d=abs(sparsity-d2)
      if (d<d_min) then
        d_min=d
        orb_r_cut_i=template_cutoffs(i)
      end if
    end do
  else
    ! try to match the inputted value of orb_r_cut 
    j=999999
    l=nint(orb_r_cut*100.0_dp)
    do i=1,template_num_cutoffs
      k=abs(l-template_cutoffs(i))
      if (k<j) then
        j=k
        orb_r_cut_i=template_cutoffs(i)
      end if
    end do
  end if
  orb_r_cut=real(orb_r_cut_i,dp)*0.01_dp
  if (.not. build_matrix) then
    if (mpi_rank==0) then
      write(c10a,'(i10)') orb_r_cut_i
      write(c10b,'(i10)') basis_subset(template_index)%ref
      template_filename=trim(adjustl(template_basedir))//'/'//&
                        trim(adjustl(system_label))//'_'//&
                        trim(adjustl(c10a))//'_'//&
                        trim(adjustl(c10b))//'.template'
      open(10,file=trim(adjustl(template_filename)))
        read(10,'(2x,i10)') num_template
      close(10)
      if (is_subset) num_template=num_template*num_orbs_per_atom**2/basis_subset(template_index)%ref**2
    end if
#ifdef MPI
    call mpi_bcast(num_template,1,mpi_int,0,mpi_comm,mpi_err)
#endif
    sparsity=1.0_dp-real(num_template,dp)*real(num_cells,dp)/(real(num_orbs,dp)**2)
    sparsity=max(sparsity,0.0_dp)
  end if

  ! don't allow k points in non-periodic directions
  if (num_periodic_dims<3) k_point(1)=0.0_dp
  if (num_periodic_dims<2) k_point(2)=0.0_dp
  if (num_periodic_dims<1) k_point(3)=0.0_dp

  ! build the matrices
  if (build_matrix) then

    ! read the data file
    if (mpi_rank==0) then
      write(c10a,'(i10)') orb_r_cut_i
      write(c10b,'(i10)') basis_subset(template_index)%ref
      template_filename=trim(adjustl(template_basedir))//'/'//&
                        trim(adjustl(system_label))//'_'//&
                        trim(adjustl(c10a))//'_'//&
                        trim(adjustl(c10b))//'.template'
      open(10,file=trim(adjustl(template_filename)))
    end if
      if (mpi_rank==0) read(10,'(2x,i10)') num_template
#ifdef MPI
      call mpi_bcast(num_template,1,mpi_int,0,mpi_comm,mpi_err)
#endif
      allocate(template_i(5,num_template))
      if (mpi_rank==0) read(10,'(a)') line
#ifdef MPI
      call mpi_bcast(line,200,mpi_char,0,mpi_comm,mpi_err)
#endif
      if (mpi_rank==0) backspace(10)
      read(line,*,iostat=j) template_i(1:5,1), el1, el2
      if (j==0) then
        overlap_present=.true.
        allocate(template_d(2,num_template))
        if (mpi_rank==0) then
          do i=1,num_template
            read(10,'(i3,4(1x,i3),2(1x,es10.3e2))') template_i(1:5,i), template_d(1:2,i)
          end do
        end if
#ifdef MPI
        call mpi_bcast(template_i,5*num_template,mpi_int,0,mpi_comm,mpi_err)
        call mpi_bcast(template_d,2*num_template,mpi_double_precision,0,mpi_comm,mpi_err)
#endif
      else
        overlap_present=.false.
        allocate(template_d(1,num_template))
        if (mpi_rank==0) then
          do i=1,num_template
            read(10,'(i3,4(1x,i3),1(1x,es10.3e2))') template_i(1:5,i), template_d(1,i)
          end do
        end if
#ifdef MPI
        call mpi_bcast(template_i,5*num_template,mpi_int,0,mpi_comm,mpi_err)
        call mpi_bcast(template_d,num_template,mpi_double_precision,0,mpi_comm,mpi_err)
#endif
      end if
    if (mpi_rank==0) close(10)

    ! create the matrices with MatrixSwitch
    call m_allocate(H,num_orbs,num_orbs,m_storage)
    if (overlap_present) call m_allocate(S,num_orbs,num_orbs,m_storage)

    ! map the cell indices in three directions onto a single variable
    allocate(cell_index(num_cells_dir(1),num_cells_dir(2),num_cells_dir(3)))
    c=0
    do i=1,num_cells_dir(1)
      do j=1,num_cells_dir(2)
        do k=1,num_cells_dir(3)
          cell_index(i,j,k)=c
          c=c+num_orbs_per_cell
        end do
      end do
    end do

    ! calculate the matrix elements in the required supercell
    if (gamma_point) then
      if (is_subset) then
        ! case 1: Gamma-point (real) calculation, basis is a subset
        allocate(subset_convert(num_atoms_per_cell*basis_subset(template_index)%ref))
        subset_convert=0
        do i=1,num_orbs_per_cell
          subset_convert(basis_subset(template_index)%index(i))=i
        end do
        do l=1,num_template
          if ((subset_convert(template_i(1,l))/=0) .and. &
              (subset_convert(template_i(2,l))/=0)) then
            do i=1,num_cells_dir(1)
              do j=1,num_cells_dir(2)
                do k=1,num_cells_dir(3)
                  i2=modulo(i+template_i(3,l)-1,num_cells_dir(1))+1
                  j2=modulo(j+template_i(4,l)-1,num_cells_dir(2))+1
                  k2=modulo(k+template_i(5,l)-1,num_cells_dir(3))+1
                  a=cell_index(i,j,k)+subset_convert(template_i(1,l))
                  b=cell_index(i2,j2,k2)+subset_convert(template_i(2,l))
                  call m_set_element(H,a,b,template_d(1,l),1.0_dp)
                  if (overlap_present) call m_set_element(S,a,b,template_d(2,l),1.0_dp)
                end do
              end do
            end do
          end if
        end do
        if (defect) then
          do l=1,num_template
            if ((subset_convert(template_i(1,l))/=0) .and. &
                (subset_convert(template_i(2,l))/=0) .and. &
                (subset_convert(template_i(1,l))<=num_orbs_per_atom) .and. &
                (subset_convert(template_i(2,l))<=num_orbs_per_atom) .and. &
                (template_i(3,l)==0) .and. &
                (template_i(4,l)==0) .and. &
                (template_i(5,l)==0)) then
              a=subset_convert(template_i(1,l))
              b=subset_convert(template_i(2,l))
              call m_set_element(H,a,b,defect_perturbation*template_d(1,l),1.0_dp)
            end if
          end do
        end if
        deallocate(subset_convert)
      else
        ! case 2: Gamma-point (real) calculation, basis is not a subset
        do l=1,num_template
          do i=1,num_cells_dir(1)
            do j=1,num_cells_dir(2)
              do k=1,num_cells_dir(3)
                i2=modulo(i+template_i(3,l)-1,num_cells_dir(1))+1
                j2=modulo(j+template_i(4,l)-1,num_cells_dir(2))+1
                k2=modulo(k+template_i(5,l)-1,num_cells_dir(3))+1
                a=cell_index(i,j,k)+template_i(1,l)
                b=cell_index(i2,j2,k2)+template_i(2,l)
                call m_set_element(H,a,b,template_d(1,l),1.0_dp)
                if (overlap_present) call m_set_element(S,a,b,template_d(2,l),1.0_dp)
              end do
            end do
          end do
        end do
        if (defect) then
          do l=1,num_template
            if ((template_i(1,l)<=num_orbs_per_atom) .and. &
                (template_i(2,l)<=num_orbs_per_atom) .and. &
                (template_i(3,l)==0) .and. &
                (template_i(4,l)==0) .and. &
                (template_i(5,l)==0)) then
              a=template_i(1,l)
              b=template_i(2,l)
              call m_set_element(H,a,b,defect_perturbation*template_d(1,l),1.0_dp)
            end if
          end do
        end if
      end if
    else
      if (is_subset) then
        ! case 3: k-point (complex) calculation, basis is a subset
        allocate(subset_convert(num_atoms_per_cell*basis_subset(template_index)%ref))
        subset_convert=0
        do i=1,num_orbs_per_cell
          subset_convert(basis_subset(template_index)%index(i))=i
        end do
        do l=1,num_template
          if ((subset_convert(template_i(1,l))/=0) .and. &
              (subset_convert(template_i(2,l))/=0)) then
            do i=1,num_cells_dir(1)
              do j=1,num_cells_dir(2)
                do k=1,num_cells_dir(3)
                  i2=modulo(i+template_i(3,l)-1,num_cells_dir(1))+1
                  j2=modulo(j+template_i(4,l)-1,num_cells_dir(2))+1
                  k2=modulo(k+template_i(5,l)-1,num_cells_dir(3))+1
                  i3=(i+template_i(3,l)-i2)/num_cells_dir(1)
                  j3=(j+template_i(4,l)-j2)/num_cells_dir(2)
                  k3=(k+template_i(5,l)-k2)/num_cells_dir(3)
                  kdotT=twoPi*(k_point(1)*i3+&
                               k_point(2)*j3+&
                               k_point(3)*k3)
                  a=cell_index(i,j,k)+subset_convert(template_i(1,l))
                  b=cell_index(i2,j2,k2)+subset_convert(template_i(2,l))
                  call m_set_element(H,a,b,cmplx(template_d(1,l),0.0_dp,dp)*exp(cmplx_i*kdotT),cmplx_1)
                  if (overlap_present) call m_set_element(S,a,b,cmplx(template_d(2,l),0.0_dp,dp)*exp(cmplx_i*kdotT),cmplx_1)
                end do
              end do
            end do
          end if
        end do
        if (defect) then
          do l=1,num_template
            if ((subset_convert(template_i(1,l))/=0) .and. &
                (subset_convert(template_i(2,l))/=0) .and. &
                (subset_convert(template_i(1,l))<=num_orbs_per_atom) .and. &
                (subset_convert(template_i(2,l))<=num_orbs_per_atom) .and. &
                (template_i(3,l)==0) .and. &
                (template_i(4,l)==0) .and. &
                (template_i(5,l)==0)) then
              a=subset_convert(template_i(1,l))
              b=subset_convert(template_i(2,l))
              call m_set_element(H,a,b,cmplx(defect_perturbation*template_d(1,l),0.0_dp,dp),cmplx_1)
            end if
          end do
        end if
        deallocate(subset_convert)
      else
        ! case 4: k-point (complex) calculation, basis is not a subset
        do l=1,num_template
          do i=1,num_cells_dir(1)
            do j=1,num_cells_dir(2)
              do k=1,num_cells_dir(3)
                i2=modulo(i+template_i(3,l)-1,num_cells_dir(1))+1
                j2=modulo(j+template_i(4,l)-1,num_cells_dir(2))+1
                k2=modulo(k+template_i(5,l)-1,num_cells_dir(3))+1
                i3=(i+template_i(3,l)-i2)/num_cells_dir(1)
                j3=(j+template_i(4,l)-j2)/num_cells_dir(2)
                k3=(k+template_i(5,l)-k2)/num_cells_dir(3)
                kdotT=twoPi*(k_point(1)*i3+&
                             k_point(2)*j3+&
                             k_point(3)*k3)
                a=cell_index(i,j,k)+template_i(1,l)
                b=cell_index(i2,j2,k2)+template_i(2,l)
                call m_set_element(H,a,b,cmplx(template_d(1,l),0.0_dp,dp)*exp(cmplx_i*kdotT),cmplx_1)
                if (overlap_present) call m_set_element(S,a,b,cmplx(template_d(2,l),0.0_dp,dp)*exp(cmplx_i*kdotT),cmplx_1)
              end do
            end do
          end do
        end do
        if (defect) then
          do l=1,num_template
            if ((template_i(1,l)<=num_orbs_per_atom) .and. &
                (template_i(2,l)<=num_orbs_per_atom) .and. &
                (template_i(3,l)==0) .and. &
                (template_i(4,l)==0) .and. &
                (template_i(5,l)==0)) then
              a=template_i(1,l)
              b=template_i(2,l)
              call m_set_element(H,a,b,cmplx(defect_perturbation*template_d(1,l),0.0_dp,dp),cmplx_1)
            end if
          end do
        end if
      end if
    end if

    ! calculate a better estimate of the sparsity by randomly sampling the
    ! matrix
    if (mpi_rank==0) seed=omm_rand_seed()
#ifdef MPI
    call mpi_bcast(seed,1,mpi_int,0,mpi_comm,mpi_err)
#endif
    nzel=0
    do i=1,num_orbs
      call omm_bsd_lcg(seed,rn)
      a=ceiling(num_orbs*rn)
      call omm_bsd_lcg(seed,rn)
      b=ceiling(num_orbs*rn)
      if (gamma_point) then
        call m_get_element(H,a,b,el1)
        if (el1==0.0_dp) nzel=nzel+1
      else
        call m_get_element(H,a,b,cel1)
        if (cel1==cmplx_0) nzel=nzel+1
      end if
    end do
    sparsity=real(nzel,dp)/(real(num_orbs,dp))
    sparsity=max(sparsity,0.0_dp)

    deallocate(cell_index)
    deallocate(template_d)
    deallocate(template_i)

  end if

  do i=1,template_num_basis_sizes
    if (allocated(basis_subset(i)%index)) deallocate(basis_subset(i)%index)
  end do
  deallocate(basis_subset)
  deallocate(template_basis_sizes)
  deallocate(template_cutoffs)

contains

    !================================================!
    ! random number generator                        !
    ! -initialize with:                              !
    !  call rand_init()                              !
    ! -generate new number with:                     !
    !  call random_number(rn)                        !
    !  where where rn is a real(dp) variable         !
    !================================================!
    function omm_rand_seed() result(seed)

        implicit none
        integer :: seed
        character(10) :: system_time
        real(dp) :: rtime
#ifdef NORAND
        seed=123456
#else
        call date_and_time(time=system_time)
        read (system_time,*) rtime
        seed = int(rtime*1000.0_dp)
#endif

    end function omm_rand_seed

    subroutine omm_bsd_lcg(x, r)

        ! x_{n+1} = (a * x_{n} + c) mod m
        ! r = x_{n+1) / m

        implicit none

        integer, intent(inout) :: x
        real(dp), intent(out) :: r
        integer(i64), parameter :: a = 1103515245_i64
        integer, parameter :: c = 12345
        integer(i64), parameter :: m = 2_i64**31

        x = int(mod(a*x+c,m))
        r = real(x,dp)/m

    end subroutine omm_bsd_lcg

end subroutine tomato_TB
