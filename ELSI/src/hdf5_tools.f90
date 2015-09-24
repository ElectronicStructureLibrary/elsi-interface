!> ELSI Interchange
!! 
!! Copyright of the original code rests with the authors inside the ELSI
!! consortium. The copyright of any additional modifications shall rest
!! with their original authors, but shall adhere to the licensing terms
!! distributed along with the original code in the file "COPYING".

module HDF5_TOOLS

  use iso_c_binding
  use HDF5
  use DIMENSIONS

  implicit none
  private

  !> HDF 5 Hyperslap specifics
  integer(HSIZE_T) :: global_dim(2)  !< Global matrix dimension
  integer(HSIZE_T) :: local_dim(2)  !< Local matrix dimension
  integer(HSIZE_T) :: process_grid(2)  !< Processor grid dimensions
  integer(HSIZE_T) :: process_position(2) !< Position in process grid
  integer(HSIZE_T) :: count(2)  !< Number of blocks 
  integer(HSIZE_T) :: offset(2) !< Matrix offset 
  integer(HSIZE_T) :: stride(2) !< Distance between blocks 
  integer(HSIZE_T) :: block(2)  !< Dimension of block
  logical :: incomplete         !< Is block description complete?


  interface hdf5_write_attribute
     module procedure hdf5_write_attribute_integer
  end interface

  interface hdf5_read_attribute
     module procedure hdf5_read_attribute_integer
  end interface

  interface hdf5_write_matrix_parallel
     module procedure hdf5_write_matrix_parallel_double
  end interface

  interface hdf5_read_matrix_parallel
     module procedure hdf5_read_matrix_parallel_double
  end interface

  public :: hdf5_initialize !< Initialization of the HDF5 infrastructure
  public :: hdf5_create_file
  public :: hdf5_open_file
  public :: hdf5_close_file
  public :: hdf5_create_group
  public :: hdf5_open_group
  public :: hdf5_close_group
  public :: hdf5_get_scalapack_pattern
  public :: hdf5_write_attribute 
  public :: hdf5_read_attribute
  public :: hdf5_write_matrix_parallel
  public :: hdf5_read_matrix_parallel
  public :: hdf5_finalize

  contains

!>
!! Initialization of HDF5 infrastructure
!!
subroutine hdf5_initialize()
   
   implicit none

   call H5open_f( h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to initialize Fortran interface."
      stop
   end if

end subroutine

!>
!! Shutdown of HDF5 infrastructure
!!
subroutine hdf5_finalize()
   
   implicit none

   call H5close_f( h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to finalize Fortran interface."
      stop
   end if

end subroutine


!>
!! HDF5 File creation
!!
subroutine hdf5_create_file(file_name, mpi_comm_global, &
                             mpi_info_null, file_id)
  
   implicit none
   character(len=*), intent(in) :: file_name
   integer, intent(in) :: mpi_comm_global
   integer, intent(in) :: mpi_info_null
   integer, intent(out) :: file_id

   integer(hid_t) :: plist_id

   call H5Pcreate_f( H5P_FILE_ACCESS_F, plist_id, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create property list for parallel access."
      stop
   end if

   call H5Pset_fapl_mpio_f(plist_id, mpi_comm_global, mpi_info_null, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to set property list for parallel access."
      stop
   end if

   call H5Fcreate_f( file_name, H5F_ACC_TRUNC_F, file_id, h5err, & 
                     access_prp = plist_id)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create file for parallel access."
      stop
   end if

end subroutine

!>
!! HDF5 File open for read only
!!
subroutine hdf5_open_file(file_name, mpi_comm_global, &
                          mpi_info_null, file_id)
  
   implicit none
   character(len=*), intent(in) :: file_name
   integer, intent(in) :: mpi_comm_global
   integer, intent(in) :: mpi_info_null
   integer, intent(out) :: file_id

   integer(hid_t) :: plist_id

   call H5Pcreate_f( H5P_FILE_ACCESS_F, plist_id, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create property list for parallel access."
      stop
   end if

   call H5Pset_fapl_mpio_f(plist_id, mpi_comm_global, mpi_info_null, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to set property list for parallel access."
      stop
   end if

   call H5Fopen_f( file_name, H5F_ACC_RDONLY_F, file_id, h5err, & 
                     access_prp = plist_id)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create file for parallel access."
      stop
   end if

end subroutine


!>
!! HDF5 File close
!!
subroutine hdf5_close_file(file_id)

   implicit none
   integer, intent(in) :: file_id

   call H5Fclose_f( file_id, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close file for parallel access."
      stop
   end if


end subroutine

!>
!! HDF5 Group creation
!!
subroutine hdf5_create_group (place_id, group_name, group_id)

   implicit none

   integer(hid_t), intent(in)   :: place_id
   character(len=*), intent(in) :: group_name
   integer(hid_t), intent(out)  :: group_id

   call H5Gcreate_f( place_id, group_name, group_id, h5err)
   if (h5err /= 0) then
      write(*,'(a,a)') "HDF5: Failed to create group ", group_name
      stop
   end if
 
end subroutine

!>
!! HDF5 Open existing group
!!
subroutine hdf5_open_group (place_id, group_name, group_id)

   implicit none

   integer(hid_t), intent(in)   :: place_id
   character(len=*), intent(in) :: group_name
   integer(hid_t), intent(out)  :: group_id

   call H5Gopen_f( place_id, group_name, group_id, h5err)
   if (h5err /= 0) then
      write(*,'(a,a)') "HDF5: Failed to create group ", group_name
      stop
   end if

end subroutine


!>
!! HDF5 close group
!!
subroutine hdf5_close_group (group_id)

   implicit none

   integer(hid_t), intent(in)  :: group_id

   call H5Gclose_f( group_id, h5err)
   if (h5err /= 0) then
      write(*,'(a,a)') "HDF5: Failed to close group."
      stop
   end if

end subroutine

!>
!! HDF5 write an integer atrribute
!!
subroutine hdf5_write_attribute_integer(place_id, attr_name, attr_value)

   implicit none

   integer(hid_t), intent(in)   :: place_id
   character(len=*), intent(in) :: attr_name
   integer, intent(in)          :: attr_value

   integer(hid_t) :: space_id
   integer(hid_t) :: attr_id
   integer(hsize_t) :: dims(1)


   call H5Screate_f( H5S_SCALAR_F, space_id, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create Scalar Space Identifier."
      stop
   end if

   call H5Acreate_f( place_id, attr_name, H5T_NATIVE_INTEGER, space_id, &
                     attr_id, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create Attribute Identifier."
      stop
   end if

   dims = 1
   call H5Awrite_f ( attr_id, H5T_NATIVE_INTEGER, attr_value, dims, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to write to Attribute Identifier."
      stop
   end if

   call H5Aclose_f ( attr_id, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close to Attribute Identifier."
      stop
   end if

   call H5Sclose_f ( space_id, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close to Space Identifier."
      stop
   end if

end subroutine

!>
!! HDF5 read an integer atrribute
!!
subroutine hdf5_read_attribute_integer(place_id, attr_name, attr_value)

   implicit none

   integer(hid_t), intent(in)   :: place_id
   character(len=*), intent(in) :: attr_name
   integer, intent(out)         :: attr_value

   integer(hid_t) :: space_id
   integer(hid_t) :: attr_id
   integer(hsize_t) :: dims(1)


   call H5Aopen_f( place_id, attr_name, attr_id, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create Attribute Identifier."
      stop
   end if

   dims = 1
   call H5Aread_f ( attr_id, H5T_NATIVE_INTEGER, attr_value, dims, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to write to Attribute Identifier."
      stop
   end if

   call H5Aclose_f ( attr_id, h5err )
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close to Attribute Identifier."
      stop
   end if

end subroutine


!>
!! HDF5 write a real valued matrix
!!
subroutine hdf5_write_matrix_parallel_double(place_id, matrix_name, matrix)

   implicit none

   integer(hid_t), intent(in) :: place_id
   character(len=*), intent(in) :: matrix_name 
   real*8, intent(in) :: matrix(n_l_rows,n_l_cols)

   integer(hid_t) :: memspace
   integer(hid_t) :: filespace
   integer(hid_t) :: data_id
   integer(hid_t) :: plist_id

   integer(HSIZE_T) :: add_count(2) !< Number of blocks 
   integer(HSIZE_T) :: add_offset(2) !< Matrix offset 
   integer(HSIZE_T) :: miss_count(2) !< Number of blocks 
   integer(HSIZE_T) :: miss_offset(2) !< Matrix offset
   integer          :: i_block 


   call H5Screate_simple_f (2, local_dim, memspace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create local Matrix Space."
      stop
   end if

   call H5Screate_simple_f (2, global_dim, filespace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create global Matrix Space."
      stop
   end if


   call H5Dcreate_f (place_id, matrix_name, H5T_NATIVE_REAL, filespace,&
                   data_id, h5err) 
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create Matrix Data Identifier."
      stop
   end if

   call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count,&
         h5err, stride, block)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to write select Hyperslap."
      stop
   end if

   if (incomplete) then
      miss_count = local_dim - block * count
      miss_offset = count * stride + offset 
      
      ! row corner blocks
      if (miss_count(2) /= 0) then
         do i_block = 0, count(1) - 1
            add_count = (/block(1),miss_count(2)/)
            add_offset = (/i_block * stride(1) + offset(1), miss_offset(2)/)
            call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
               add_count, h5err)
            if (h5err /= 0) then
               write(*,'(a)') "HDF5: Failed to write select Hyperslap."
               stop
            end if
         end do
      end if

      ! column corner blocks
      if (miss_count(1) /= 0) then
         do i_block = 0, count(2) - 1
            add_count = (/miss_count(1),block(2)/)
            add_offset = (/miss_offset(1),i_block * stride(2) + offset(2)/)
            call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
               add_count, h5err)
            if (h5err /= 0) then
               write(*,'(a)') "HDF5: Failed to write select Hyperslap."
               stop
            end if
         end do
      end if

      ! Last missing partial block
      if (miss_count(1) * miss_count(2) /= 0) then
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, miss_offset,&
               miss_count, h5err)
         if (h5err /= 0) then
            write(*,'(a)') "HDF5: Failed to write select Hyperslap."
            stop
         end if
      end if

   end if

   CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
   CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

   call H5Dwrite_f (data_id, H5T_NATIVE_DOUBLE, matrix, &
         local_dim, h5err, mem_space_id = memspace, file_space_id = filespace, &
         xfer_prp = plist_id)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to write Matrix to Filespace."
      stop
   end if


   call H5Sclose_f(filespace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close Filespace."
      stop
   end if


   call H5Sclose_f(memspace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close Memspace."
      stop
   end if


   call H5Dclose_f(data_id, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close Data Identifier."
      stop
   end if


end subroutine

!>
!! HDF5 read a real valued matrix
!!
subroutine hdf5_read_matrix_parallel_double(place_id, matrix_name, matrix)

   implicit none

   integer(hid_t), intent(in) :: place_id
   character(len=*), intent(in) :: matrix_name 
   real*8, intent(out) :: matrix(n_l_rows,n_l_cols)

   integer(hid_t) :: memspace
   integer(hid_t) :: filespace
   integer(hid_t) :: data_id
   integer(hid_t) :: plist_id

   integer(HSIZE_T) :: add_count(2) !< Number of blocks 
   integer(HSIZE_T) :: add_offset(2) !< Matrix offset
   integer(HSIZE_T) :: miss_count(2) !< Number of blocks 
   integer(HSIZE_T) :: miss_offset(2) !< Matrix offset
   integer          :: i_block 
 

   call H5Screate_simple_f (2, local_dim, memspace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create local Matrix Space."
      stop
   end if

   call H5Screate_simple_f (2, global_dim, filespace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create global Matrix Space."
      stop
   end if


   call H5Dopen_f (place_id, matrix_name, data_id, h5err) 
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to create Matrix Data Identifier."
      stop
   end if

   call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count,&
         h5err, stride, block)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to write select Hyperslap."
      stop
   end if

   if (incomplete) then
      miss_count = local_dim - block * count
      miss_offset = count * stride + offset 
      
      ! row corner blocks
      if (miss_count(2) /= 0) then
         do i_block = 0, count(1) - 1
            add_count = (/block(1),miss_count(2)/)
            add_offset = (/i_block * stride(1) + offset(1), miss_offset(2)/)
            call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
               add_count, h5err)
            if (h5err /= 0) then
               write(*,'(a)') "HDF5: Failed to write select Hyperslap."
               stop
            end if
         end do
      end if

      ! column corner blocks
      if (miss_count(1) /= 0) then
         do i_block = 0, count(2) - 1
            add_count = (/miss_count(1),block(2)/)
            add_offset = (/miss_offset(1),i_block * stride(2) + offset(2)/)
            call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
               add_count, h5err)
            if (h5err /= 0) then
               write(*,'(a)') "HDF5: Failed to write select Hyperslap."
               stop
            end if
         end do
      end if

      ! Last missing partial block
      if (miss_count(1) * miss_count(2) /= 0) then
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, miss_offset,&
               miss_count, h5err)
         if (h5err /= 0) then
            write(*,'(a)') "HDF5: Failed to write select Hyperslap."
            stop
         end if
      end if

   end if

   CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
   CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)


   call H5Dread_f (data_id, H5T_NATIVE_DOUBLE, matrix, &
         local_dim, h5err, mem_space_id = memspace, file_space_id = filespace, &
         xfer_prp = plist_id)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to write Matrix to Filespace."
      stop
   end if


   call H5Sclose_f(filespace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close Filespace."
      stop
   end if


   call H5Sclose_f(memspace, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close Memspace."
      stop
   end if


   call H5Dclose_f(data_id, h5err)
   if (h5err /= 0) then
      write(*,'(a)') "HDF5: Failed to close Data Identifier."
      stop
   end if


end subroutine


!>
!! HDF5 calculate pattern for scalapack distribution
!!
subroutine hdf5_get_scalapack_pattern()

   ! Follow: https://www.hdfgroup.org/HDF5/Tutor/selectsimple.html

   implicit none

   integer :: count_r, count_c
   integer :: stride_r, stride_c

   global_dim       = (/n_g_rank,n_g_rank/)
   local_dim        = (/n_l_rows,n_l_cols/)
   block            = (/n_b_rows,n_b_cols/)
   process_position = (/my_p_row,my_p_col/)
   process_grid     = (/n_p_rows,n_p_cols/)

   offset = process_position * block
   stride = process_grid * block
   
   count_r = FLOOR (1d0 * n_l_rows / n_b_rows)
   count_c = FLOOR (1d0 * n_l_cols / n_b_cols)  
   count   = (/count_r,count_c/)

   if (count_r * n_b_rows /= n_l_rows .or. count_c * n_b_cols /= n_l_cols) then
     incomplete = .True.
   else
     incomplete = .False.
   end if

   if (myid == 0) then
      print *, "global dim    :",  global_dim
      print *, "local dim     :",  local_dim
      print *, "processor dim :",  process_grid
      print *, "block         :",  block
      print *, "offset        :",  offset
      print *, "count         :",  count
      print *, "stride        :",  stride
   end if

end subroutine

end module HDF5_TOOLS
