!Copyright (c) 2015, ELSI consortium
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer in the
!   documentation and/or other materials provided with the distribution.
! * Neither the name of the ELSI project nor the
!   names of its contributors may be used to endorse or promote products
!   derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This module contains wrapper functions to the HDF5 library 
!!

module ELSI_HDF5_TOOLS

  use iso_c_binding
  use HDF5
  use ELSI_DIMENSIONS

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
  integer(HSIZE_T) :: chunk(2)  !< Dimension of chunks
  logical :: incomplete         !< Is block description complete?

  integer(HSIZE_T) :: nslots = 523 !< Number of Slots for chunk hash table
  integer(HSIZE_T) :: nchunk_bytes = 1*1024*1024 !< Size of chunk cache


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
   
   character*40 :: caller = "hdf5_initialize"

   call H5open_f( h5err )
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to initialize Fortran interface.",caller)
   end if

end subroutine

!>
!! Shutdown of HDF5 infrastructure
!!
subroutine hdf5_finalize()
   
   implicit none

   character*40 :: caller = "hdf5_finalize"

   call H5close_f( h5err )
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to finalize Fortran interface.",caller)
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
   character*40 :: caller = "hdf5_create_file"

   integer(hid_t) :: plist_id

   call H5Pcreate_f( H5P_FILE_ACCESS_F, plist_id, h5err )
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create parallel property list.", &
            caller)
   end if

   call H5Pset_fapl_mpio_f(plist_id, mpi_comm_global, mpi_info_null, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to set MPI-IO for parallel access.", &
            caller)
   end if

   call H5Fcreate_f( file_name, H5F_ACC_TRUNC_F, file_id, h5err, & 
                     access_prp = plist_id)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create file for parallel access.", caller)
   end if

end subroutine

!>
!! HDF5 File open for read only
!!
subroutine hdf5_open_file(file_name, mpi_comm_global, &
                          mpi_info_null, file_id)
  
   implicit none
   character(len=*), intent(in) :: file_name
   integer, intent(in)  :: mpi_comm_global
   integer, intent(in)  :: mpi_info_null
   integer, intent(out) :: file_id
   character*40         :: caller = "hdf5_open_file"

   integer(hid_t) :: plist_id

   call H5Pcreate_f( H5P_FILE_ACCESS_F, plist_id, h5err )
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create parallel property list.",&
            caller)
   end if

   call H5Pset_fapl_mpio_f(plist_id, mpi_comm_global, mpi_info_null, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to set MPI-IO for parallel access.",&
            caller)
   end if

   call H5Fopen_f( file_name, H5F_ACC_RDONLY_F, file_id, h5err, & 
                     access_prp = plist_id)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to open file for parallel access.",caller)
   end if

end subroutine


!>
!! HDF5 File close
!!
subroutine hdf5_close_file(file_id)

   implicit none
   integer, intent(in) :: file_id

   character*40 :: caller = "hdf5_close_file"

   call H5Fclose_f( file_id, h5err )
   if (h5err /= 0) then
      call elsi_stop ("HDF5: Failed to close file for parallel access.", caller)
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

   character*40 :: caller = "hdf5_create_group"

   call H5Gcreate_f( place_id, group_name, group_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create group " // group_name, caller)
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

   character*40 :: caller = "hdf5_open_group"

   call H5Gopen_f( place_id, group_name, group_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create group " // group_name, caller)
   end if

end subroutine


!>
!! HDF5 close group
!!
subroutine hdf5_close_group (group_id)

   implicit none

   integer(hid_t), intent(in)  :: group_id

   character*40 :: caller = "hdf5_close_group"

   call H5Gclose_f( group_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close group.", caller)
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

   integer(hid_t)   :: space_id
   integer(hid_t)   :: attr_id
   integer(hsize_t) :: dims(1)

   character*40     :: caller = "hdf5_write_attribute_integer"

   call H5Screate_f( H5S_SCALAR_F, space_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create Scalar Space Identifier.",caller)
   end if

   call H5Acreate_f( place_id, attr_name, H5T_NATIVE_INTEGER, space_id, &
                     attr_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create Attribute Identifier.",caller)
   end if

   dims = 1
   call H5Awrite_f ( attr_id, H5T_NATIVE_INTEGER, attr_value, dims, h5err )
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to write to Attribute Identifier.",caller)
   end if

   call H5Aclose_f ( attr_id, h5err )
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close to Attribute Identifier.",caller)
   end if

   call H5Sclose_f ( space_id, h5err )
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close to Space Identifier.",caller)
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
   integer(hid_t) :: dataprop_id
   integer(hid_t) :: plist_id

   integer(HSIZE_T) :: add_count(2)   !< Number of blocks 
   integer(HSIZE_T) :: add_offset(2)  !< Matrix offset 
   integer(HSIZE_T) :: add_block(2)   !< Incomplete blocksize 
   integer(HSIZE_T) :: miss_count(2)  !< Number of incomplete blocks 
   integer(HSIZE_T) :: miss_offset(2) !< Offset of incomplete blocks
   character*40     :: caller = "hdf5_write_matrix_parallel_double"

   call H5Screate_simple_f (2, local_dim, memspace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create local Matrix Space.",caller)
   end if

   call H5Screate_simple_f (2, global_dim, filespace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create global Matrix Space.",caller)
   end if

   ! Create dataset chunks according to blocksize
   call H5Pcreate_f(H5P_DATASET_CREATE_F, dataprop_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create dataset property list.",caller)
   end if

   call H5Pset_chunk_f(dataprop_id,2,chunk,h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create dataset chunks property.",caller)
   end if

   call H5Dcreate_f (place_id, matrix_name, H5T_NATIVE_REAL, filespace,&
                   data_id, h5err, dataprop_id) 
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create Matrix Data Identifier.",caller)
   end if

   call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count,&
         h5err, stride, block)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to write select Hyperslap.", caller)
   end if

   if (incomplete) then
      miss_count = local_dim - block * count
      miss_offset = count * stride + offset

      ! only row corner blocks
      if (miss_count(2) /= 0 .and. miss_count(1) == 0) then
         add_block  = (/block(1),miss_count(2)/)
         add_offset = (/offset(1), miss_offset(2)/)
         add_count(1) = count(1)
         add_count(2) = 1

         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err, stride, add_block)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

      end if

      ! only column corner blocks
      if (miss_count(1) /= 0 .and. miss_count(2) == 0) then
         add_block  = (/miss_count(1),block(2)/)
         add_offset = (/miss_offset(1),offset(2)/)
         add_count(1) = 1
         add_count(2) = count(2)

         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err, stride, add_block)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

      end if

      ! both
      if (miss_count(1) /= 0 .and. miss_count(2) /= 0) then
         ! Edge block
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, miss_offset,&
               miss_count, h5err)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

         ! side rows
         add_block  = (/block(1),miss_count(2)/)
         add_offset = (/offset(1), miss_offset(2)/)
         add_count(1) = count(1)
         add_count(2) = 1

         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err,stride,add_block)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

         ! side columns
         add_block  = (/miss_count(1),block(2)/)
         add_offset = (/miss_offset(1),offset(2)/)
         add_count(1) = 1
         add_count(2) = count(2)

         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err,stride,add_block)
         if (h5err /= 0) then
              call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

      end if

   end if

   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create dataset property list.",caller)
   end if
 
   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to MPI-IO.",caller)
   end if

   call H5Dwrite_f (data_id, H5T_NATIVE_DOUBLE, matrix, &
         local_dim, h5err, mem_space_id = memspace, file_space_id = filespace, &
         xfer_prp = plist_id)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to write Matrix to Filespace.",caller)
   end if

   call H5Sclose_f(filespace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close Filespace.",caller)
   end if

   call H5Sclose_f(memspace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close Memspace.",caller)
   end if

   call H5Dclose_f(data_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close Data Identifier.",caller)
   end if

   call H5Pclose_f(dataprop_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close Data Identifier.",caller)
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
   integer(hid_t) :: dataprop_id
   integer(hid_t) :: plist_id

   integer(HSIZE_T) :: add_count(2)   !< Number of blocks 
   integer(HSIZE_T) :: add_offset(2)  !< Matrix offset
   integer(HSIZE_T) :: add_block(2)   !< Incomplete blocksize
   integer(HSIZE_T) :: miss_count(2)  !< Number of incomplete blocks 
   integer(HSIZE_T) :: miss_offset(2) !< Offset of incomplete blocks
   character*40     :: caller = "hdf5_read_matrix_parallel_double" 

   call H5Screate_simple_f (2, local_dim, memspace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create local Matrix Space.",caller)
   end if

   call H5Screate_simple_f (2, global_dim, filespace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create global Matrix Space.",caller)
   end if

   call H5Dopen_f (place_id, matrix_name, data_id, h5err) 
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to create Matrix Data Identifier.", caller)
   end if

   call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count,&
         h5err, stride, block)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
   end if

   if (incomplete) then
      miss_count = local_dim - block * count
      miss_offset = count * stride + offset 
      
      ! only row corner blocks
      if (miss_count(2) /= 0 .and. miss_count(1) == 0) then
         add_block  = (/block(1),miss_count(2)/)
         add_offset = (/offset(1), miss_offset(2)/)
         add_count(1) = count(1)
         add_count(2) = 1
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err, stride, add_block)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

      end if

      ! only column corner blocks
      if (miss_count(1) /= 0 .and. miss_count(2) == 0) then
         add_block  = (/miss_count(1),block(2)/)
         add_offset = (/miss_offset(1),offset(2)/)
         add_count(1) = 1
         add_count(2) = count(2)
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err, stride, add_block)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if
   
      end if

      ! both
      if (miss_count(1) /= 0 .and. miss_count(2) /= 0) then
         ! Edge block
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, miss_offset,&
               miss_count, h5err)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

         ! side rows
         add_block  = (/block(1),miss_count(2)/)
         add_offset = (/offset(1), miss_offset(2)/)
         add_count(1) = count(1) 
         add_count(2) = 1
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err, stride, add_block)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

         ! side columns
         add_block  = (/miss_count(1),block(2)/)
         add_offset = (/miss_offset(1),offset(2)/)
         add_count(1) = 1 
         add_count(2) = count(2)
         call H5Sselect_hyperslab_f (filespace, H5S_SELECT_OR_F, add_offset,&
            add_count, h5err, stride, add_block)
         if (h5err /= 0) then
            call elsi_stop("HDF5: Failed to write select Hyperslap.",caller)
         end if

      end if

   end if

   CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, h5err) 
   CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, h5err)

   call H5Dread_f (data_id, H5T_NATIVE_DOUBLE, matrix, &
         local_dim, h5err, mem_space_id = memspace, file_space_id = filespace, &
         xfer_prp = plist_id)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to write Matrix to Filespace.",caller)
   end if

   call H5Sclose_f(filespace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close Filespace.",caller)
   end if

   call H5Sclose_f(memspace, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close Memspace.",caller)
   end if

   call H5Dclose_f(data_id, h5err)
   if (h5err /= 0) then
      call elsi_stop("HDF5: Failed to close Data Identifier.",caller)
   end if

end subroutine


!>
!! HDF5 calculate pattern for scalapack distribution
!!
subroutine hdf5_get_scalapack_pattern()

   ! Follow: https://www.hdfgroup.org/HDF5/Tutor/selectsimple.html

   implicit none

   global_dim       = (/n_g_rank,n_g_rank/)
   local_dim        = (/n_l_rows,n_l_cols/)
   block            = (/n_b_rows,n_b_cols/)
   process_position = (/my_p_row,my_p_col/)
   process_grid     = (/n_p_rows,n_p_cols/)

   if (method == PEXSI) then

      offset(1) = 0 
      offset(2) = myid * FLOOR (1d0 * global_dim(2) / n_procs)
      
      stride = global_dim
      
      count(1) = FLOOR (1d0 * local_dim(1) / block(1))
      count(2) = FLOOR (1d0 * local_dim(2) / block(2))  
      
      if (count(1) * block(1) /= local_dim(1) .or.&
          count(2) * block(2) /= local_dim(2)) then
        incomplete = .True.
      else
        incomplete = .False.
      end if

      chunk = (/n_g_rank,1/)

   else !(ELPA, OMM_dense are using block-cyclic for scalapack)

      offset = process_position * block
      stride = process_grid * block
   
      count(1) = FLOOR (1d0 * local_dim(1) / block(1))
      count(2) = FLOOR (1d0 * local_dim(2) / block(2))  

      if (count(1) * block(1) /= local_dim(1) .or.&
          count(2) * block(2) /= local_dim(2)) then
        incomplete = .True.
      else
        incomplete = .False.
      end if

      chunk = block

   end if

   ! For Debug
   !call elsi_hdf5_variable_status()

end subroutine


!>
!! Give Status of ELSI HDF5 variables 
!!
subroutine elsi_hdf5_variable_status()

      implicit none
      include "mpif.h"

      character(LEN=4096) :: string_message
      integer :: i_task

      do i_task = 0, n_procs
         if (i_task == myid) then
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Matrixsize global ',I5,' x ',I5)") &
              & myid, global_dim(1), global_dim(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Matrixsize local ',I5,' x ',I5)") &
              & myid, local_dim(1), local_dim(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Blocksize global ',I5,' x ',I5)") &
              & myid, block(1), block(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Processgrid global ',I5,' x ',I5)") &
              & myid, process_grid(1), process_grid(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Process ',I5,' x ',I5)") &
              & myid, process_position(1), process_position(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Count ',I5,' x ',I5)") &
              & myid, count(1), count(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Offset ',I5,' x ',I5)") &
              & myid, offset(1), offset(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Stride ',I5,' x ',I5)") &
              & myid, stride(1), stride(2)
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Incomplete? ',L2)") &
              & myid, incomplete
            write(*,'(A)') trim(string_message)
            write(string_message, "(1X,'*** Proc',I5,&
              &' : Chunk ',I5,' x ',I5)") &
              & myid, chunk(1), chunk(2)
            write(*,'(A)') trim(string_message)
         end if

         call MPI_BARRIER(mpi_comm_global,mpierr)
      end do
        
end subroutine elsi_hdf5_variable_status


end module ELSI_HDF5_TOOLS
