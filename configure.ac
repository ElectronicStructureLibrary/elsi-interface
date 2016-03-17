AC_PREREQ(2.59)
# was 2.61 but jaguar is not up to date so backported which
# involved replacing AC_PROG_SED and UINT*_T type checks.
AC_INIT([ELSI], [0.10], [http://www.elsi-interchange.org/forum/smf/], , [http://www.elsi-interchange.org])
AC_SUBST([ELSI_SO_VERSION], [0:0:1])

# This is where autoconf automatically generated files go
AC_CONFIG_AUX_DIR([config])

AM_INIT_AUTOMAKE([foreign -Wall subdir-objects])
#AM_INIT_AUTOMAKE([-Wall -Werror nostdinc foreign])

AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_MACRO_DIR([config])
AC_CONFIG_HEADERS([config.h])
#AC_CONFIG_HEADERS([ELSI/config.h])

# Without this, automake tries to be smart and rebuilt
# the autoconf generated files such as configure, aclocal.m4, etc.,
# in case the timestamps of files such as configure.ac are newer
#
# This only makes trouble for end users with out-of-date autoconf versions
# that cannot produce these files
AM_MAINTAINER_MODE([disable])

AM_SILENT_RULES([yes])


#vama## This is where local macros will be stored (by us)
#vama temp AC_CONFIG_MACRO_DIR([m4])


#vama1## Default options for automake ??
#vama1##vama#m4_pattern_allow([AM_PROG_AR])

#AX_CHECK_GNU_MAKE()
#vama#
#vama## We want #defines output into a header file
#vama#
#vama#
#vama#
AC_MSG_CHECKING(compilers C C++ Fortran)
ACX_MPI
AC_PROG_CC([$MPICC])
AC_PROG_CXX([$MPICXX])
AC_PROG_FC([$MPIF77]) # f90 is a superset of f77
AC_PROG_F77([$MPIF77])
AM_PROG_AR
#vama##vama##AC_PROG_CC
#vama##vama##AC_PROG_FC
#vama##vama## All source will be built with the MPI C++ and C compilers
#ACX_MPI
ACX_WITH_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
FCFLAGS="$FCFLAGS  $OPENMP_FCFLAGS"
FFLAGS="$FFLAGS $OPENMP_FCFLAGS"



LIBOMM_FLAGS=""
if test $acx_with_mpi != "no"; then
    LIBOMM_FLAGS="-DMPI $LIBOMM_FLAGS"
fi

LAP="yes"
SLAP="yes"
if test $LAP != "no"; then
    LIBOMM_FLAGS="-DLAP $LIBOMM_FLAGS"
fi
if test $SLAP != "no"; then
    LIBOMM_FLAGS="-DSLAP $LIBOMM_FLAGS"
fi




AC_LANG([Fortran])

AC_FC_FREEFORM
AC_FC_MODULE_FLAG
AC_FC_MODULE_OUTPUT_FLAG
AC_FC_LIBRARY_LDFLAGS

save_FCFLAGS=$FCFLAGS
save_LDFLAGS=$LDFLAGS
save_LIBS=$LIBS

#AC_ARG_VAR([SCALAPACK_LDFLAGS],[Extra LDFLAGS necessary to link a program with Scalapack])
#AC_ARG_VAR([SCALAPACK_FCFLAGS],[Extra FCFLAGS necessary to compile a Fortran program with Scalapack])

#FCFLAGS="$FCFLAGS $SCALAPACK_FCFLAGS"
#LDFLAGS="$LDFLAGS $SCALAPACK_LDFLAGS"
LIBS="$LIBS $SCALAPACK_FCFLAGS"
#vama#
#vama#
AC_MSG_CHECKING([whether Fortran module iso_fortran_env is available])
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
  program test_error_unit
    use ISO_FORTRAN_ENV, only : error_unit
    implicit none

    write(error_unit,*) "error_unit is defined"
  end program
])],
  [can_use_iso_fortran_env=yes],
  [can_use_iso_fortran_env=no]
)
AC_MSG_RESULT([${can_use_iso_fortran_env}])

dnl check whether one can link with specified MKL (desired method)
AC_MSG_CHECKING([whether we can compile a Fortran program using MKL])


AC_COMPILE_IFELSE([AC_LANG_SOURCE([
  program test_mkl
    use mkl_service
    character*198 :: string
    call mkl_get_version_string(string)
    write(*,'(a)') string
  end program
])],
  [can_compile_with_mkl=yes],
  [can_compile_with_mkl=no]
)
AC_MSG_RESULT([${can_compile_with_mkl}])

if test x"$can_compile_with_mkl" = x"yes" ; then
  AC_MSG_CHECKING([whether we can link a Fortran program with MKL])
  AC_LINK_IFELSE([AC_LANG_SOURCE([
    program test_mkl
      use mkl_service
      character*198 :: string
      call mkl_get_version_string(string)
      write(*,'(a)') string
    end program
  ])],
    [have_mkl=yes],
    [have_mkl=no]
  )
  AC_MSG_RESULT([${have_mkl}])
fi

dnl if not mkl, check all the necessary individually
if test x"${have_mkl}" = x"yes" ; then
  WITH_MKL=1
else

withscalapack=no

  dnl first check blas
  AC_SEARCH_LIBS([dgemm],[blas],[have_blas=yes],[have_blas=no])
  AC_MSG_CHECKING([whether we can link a program with a blas lib])
  AC_MSG_RESULT([${have_blas}])

  if test x"${have_blas}" = x"no" ; then
    AC_MSG_ERROR([could not link with blas: specify path]) 
  fi
  dnl now lapack
  AC_SEARCH_LIBS([dlarrv],[lapack],[have_blas=yes],[have_blas=no])
  AC_MSG_CHECKING([whether we can link a program with a lapack lib])
  AC_MSG_RESULT([${have_lapack}])

  if test x"${have_lapack}" = x"no" ; then
    AC_MSG_ERROR([could not link with lapack: specify path])
  fi

  dnl test whether scalapack already contains blacs
  scalapack_libs="scalapack mpiscalapack"
  old_LIBS="$LIBS"
  for lib in ${scalapack_libs}; do
    LIBS="-l${lib} ${old_LIBS}"
    AC_MSG_CHECKING([whether -l${lib} already contains a BLACS implementation])
    AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY([blacs_gridinit])],[blacs_in_scalapack=yes],[blacs_in_scalapack=no])
    AC_MSG_RESULT([${blacs_in_scalapack}])
    if test x"${blacs_in_scalapack}" = x"yes"; then
      break
    fi
  done

  if test x"${blacs_in_scalapack}" = x"no"; then
    LIBS="${old_LIBS}"

    dnl Test for stand-alone blacs
    AC_MSG_CHECKING([whether -lmpiblacs or -llblacs already contains a BLACS implementation])
    AC_SEARCH_LIBS([blacs_gridinit],[mpiblacs blacs],[have_blas=yes],[have_blas=no])

    if test x"${have_blacs}" = x"no"; then
      AC_MSG_ERROR([No usable BLACS found. If installed in a non-standard place, please specify suitable LDFLAGS and FCFLAGS as arguments to configure])
    fi
  fi

  if test x"${have_blacs}" != x"no"; then
     AC_SEARCH_LIBS([pdtran],[$scalapack_libs],[have_scalapack=yes],[have_scalapack=no])
  fi

  if test x"${have_scalapack}" = x"no" ; then
    AC_MSG_ERROR([could not link with scalapack: specify path])
  fi

  dnl check whether we can link alltogehter
  AC_MSG_CHECKING([whether we can link a Fortran program with all blacs/scalapack])
  AC_LINK_IFELSE([AC_LANG_SOURCE([
    program dgemm_test

      integer , parameter:: M = 4, N = 3, K = 2
      real :: A(M,K), B(K,N), C(M,N)

      call dgemm('N','N',M,N,K,1.0,A,M,B,K,0.0,C,M)

     end program dgemm_test
     ])],
     [can_link_with_blacs_scalapack=yes],
     [can_link_with_blacs_scalapack=no]
   )
   AC_MSG_RESULT([${can_link_with_blacs_scalapack}])

   if test x"${can_link_with_blacs_scalapack}" = x"yes" ; then
     WITH_BLACS=1
   else
   AC_MSG_ERROR([We can neither link with MKL or another Scalpack. Please specify SCALAPACK_LDFLAGS and SCALAPACK_FCFLAGS!])
   fi
dnl no MKL ends
fi

dnl important: reset them again!
FCFLAGS=$save_FCFLAGS
LDFLAGS=$save_LDFLAGS
LIBS=$save_LIBS

# Initialize Libtools (for shared libraries)
LT_INIT
LT_LANG([C++])
AC_SUBST(enable_shared)
AC_SUBST(enable_static)
AC_ARG_ENABLE([all-static],
    [AC_HELP_STRING([--enable-all-static@<:@=yes|no@:>@],
      [Link with all static libraries.])],
    [
      case $enableval in
      yes)
        acx_enable_all_static="yes"
      ;;
      no)
        acx_enable_all_static="no"
      ;;
      *)
        AC_MSG_ERROR([Invalid argument '$enableval' given for --enable-all-static])
      esac
    ],
    [acx_enable_all_static="$acx_default_enable_all_static"]
  )

AM_CONDITIONAL([LINK_ALL_STATIC], [test $acx_enable_all_static ="yes"])




#ACX_DETECT_CXX
# Check for C++11 support, and
# add the -std=c++11 or -std=c++0x flag to CXXFLAGS if necessary.
#AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])


# Check for HDF5
AX_LIB_HDF5([parallel])
AC_DEFINE(HAVE_HDF5)
AC_SUBST(HDF5_FLIBS)
AC_SUBST(HDF5_FFLAGS)

#Check for ParMetis
ACX_WITH_PARMETIS
#SuperLu
ACX_WITH_SUPERLU

#AC_SUBST(LIBOMM_FLAGS, [$LIBOMM_FLAGS])
AC_SUBST(LIBOMM_FLAGS)
AC_SUBST(SUPERLU_LIBS)
AC_SUBST(LIBSPARMETISLIB)

AC_SUBST(SCALAPACK_FCFLAGS)

AC_CONFIG_SUBDIRS([external/ELPA])
AC_CONFIG_FILES([
  Makefile
  external/Makefile
  external/libOMM/Makefile
  external/libOMM/src/Makefile
  external/libOMM/src/MatrixSwitch-0.1.2/Makefile
  external/libOMM/src/MatrixSwitch-0.1.2/src/Makefile
  external/libOMM/examples/Makefile
  external/PEXSI/Makefile
  external/PEXSI/src/Makefile
  external/PEXSI/fortran/Makefile
  external/PEXSI/examples/Makefile
  ELSI/Makefile
  ELSI/src/Makefile
  ELSI/example/Makefile
 ])
AC_OUTPUT
