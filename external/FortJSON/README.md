# FortJSON #

## Overview ##

FortJSON is a JSON library written in Fortran 2003. It is designed with portability across HPC architectures in mind.

Requirements:

- Fortran 2003 compiler
- cmake 3.0+

Recommended:

- doxygen 1.8.11+

## Installation ##

FortJSON uses a standard cmake (out-of-source) build system for installation:

1. Create a build subdirectory in the root FortJSON directory and change directory into it
2. Generate the build system using cmake.  There are multiple ways to do this, but two common ways are:
   - Specify a toolchain file containing compilation settings, then run `cmake -DCMAKE_TOOLCHAIN_FILE=/path/to/toolchain/file ..` to generate the build system.  An example toolchain file named "the-gibson.cmake" is provided in the root directory.
   - Run `cmake ..` to generate the build system using settings automatically determined by cmake.  This approach may select the wrong Fortran compiler and will likely attempt to install FortJSON to a folder that you do not have write access to (for example, /usr/local/).
3. Run `make` to compile FortJSON in the build subdirectory
4. Run `make install` to install FortJSON into the install directory

FortJSON uses doxygen to generate documentation:

1. Enter into the doc subdirectory in the root FortJSON directory.
2. Run `doxygen Doxyfile` to generate HTML and LaTeX versions of the documentation.

## How to Use ##

Example program:

```fortran
program hello_world

   use FortJSON

   implicit none

   type(fjson_handle) :: fj_h

   integer(kind=i4) :: my_age
   real(kind=r8)    :: my_weight

   ! Open the file
   call fjson_open_file(fj_h, 66, "hello_world.json")
   call fjson_start_array(fj_h)

   ! Rex the Dog
   my_age = 1_i4
   my_weight = 5.0_r8

   call fjson_start_object(fj_h)
   call fjson_write_name_value(fj_h, "My Age", my_age)
   call fjson_write_name_value(fj_h, "My Weight", my_weight)
   call fjson_write_name_value(fj_h, "My Name", "Rex Jr.")
   call fjson_start_name_array(fj_h, "My Favorite Things")
   call fjson_write_value(fj_h, "Meat")
   call fjson_finish_array(fj_h)
   call fjson_start_name_object(fj_h, "My Favorite People")
   call fjson_write_name_value(fj_h, "Mom", "Precious")
   call fjson_finish_object(fj_h)
   call fjson_finish_object(fj_h)

   ! Tiger the Cat
   my_age = 15_i4
   my_weight = 15.0_r8

   call fjson_start_object(fj_h)
   call fjson_write_name_value(fj_h, "My Age", my_age)
   call fjson_write_name_value(fj_h, "My Weight", my_weight)
   call fjson_write_name_value(fj_h, "My Name", "Tiger")
   call fjson_start_name_array(fj_h, "My Favorite Things")
   call fjson_finish_array(fj_h)
   call fjson_start_name_object(fj_h, "My Favorite People")
   call fjson_finish_object(fj_h)
   call fjson_finish_object(fj_h)

   call fjson_write_value(fj_h, &
                          fjson_error_message(fjson_get_last_error(fj_h)))

   ! Close the file
   call fjson_finish_array(fj_h)
   call fjson_close_file(fj_h)

end program hello_world
```

Example output:

```json
[
  {
    "My Age": 1,
    "My Weight": 0.50000000E+01,
    "My Name": "Rex Jr.",
    "My Favorite Things": [
      "Meat",
    ],
    "My Favorite People": {
      "Mom": "Precious",
    }
  },
  {
    "My Age": 15,
    "My Weight": 0.15000000E+02,
    "My Name": "Tiger",
    "My Favorite Things": [
    ],
    "My Favorite People": {
    }
  },
  "FortJSON Error:  No error"
]
```

## Additional Information ##

- The JSON standard used for FortJSON is "ECMA-404 The JSON Data Interchange Standard".  This standard is 16 pages long and contains mostly whitespace and pictures.
- A simple description of the JSON grammer, as well as other JSON libraries, may be found at <http://json.org/>.
- Google maintains a JSON style guide at <https://google.github.io/styleguide/jsoncstyleguide.xml>.
