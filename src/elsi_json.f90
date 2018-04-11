! Copyright (c) 2015-2018, the ELSI team.
! All rights reserved.
!
! This file is part of ELSI and is distributed under the BSD 3-clause license,
! which may be found in the LICENSE file in the ELSI root directory.

!>
!! This module is a wrapper around FortJSON
!!
module ELSI_JSON

   use FORTJSON, only: elsi_reset_json_handle => fjson_reset_fj_handle,&
                       elsi_open_json_file => fjson_open_file,&
                       elsi_close_json_file => fjson_close_file,&
                       elsi_write_json_name_value => fjson_write_name_value,&
                       elsi_start_json_name_object => fjson_start_name_object,&
                       elsi_start_json_object => fjson_start_object,&
                       elsi_start_json_array => fjson_start_array,&
                       elsi_finish_json_object => fjson_finish_object,&
                       elsi_finish_json_array => fjson_finish_array,&
                       elsi_get_datetime_rfc3339 => fjson_get_datetime_rfc3339,&
                       elsi_append_string => fjson_append_string,&
                       elsi_truncate_string => fjson_truncate_string,&
                       elsi_json_handle => fjson_handle

   implicit none

   public

end module ELSI_JSON
