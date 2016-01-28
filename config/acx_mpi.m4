AC_DEFUN([ACX_MPI], [
        # We were using the full macro from ACX but this forced AC_PROG_C or AC_PROG_CXX to 
        # be run before we had overridden the compilers which meant that some confdef.h
        # entries were incorrect (specifically std::exit problem with PGI)

        # Disable MPI C++ bindings.
        CPPFLAGS="$CPPFLAGS -DMPICH_SKIP_MPICXX=1 -DOMPI_SKIP_MPICXX=1"
        
        AC_ARG_VAR(MPICC,[MPI C compiler command])
        AC_CHECK_PROGS(MPICC, mpiicc mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC)
        acx_mpi_save_CC="$CC"
        CC="$MPICC"
        AC_SUBST(MPICC)



        AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
        AC_CHECK_PROGS(MPICXX, mpiicc mpicxx mpic++ mpiCC mpCC hcp mpxlC mpxlC_r cmpic++, $CXX)
        acx_mpi_save_CXX="$CXX"
        CXX="$MPICXX"
        AC_SUBST(MPICXX)

        AC_ARG_VAR(MPIF77,[MPI Fortran compiler command])
        AC_CHECK_PROGS(MPIF77, mpiifort mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c, $F77)
        acx_mpi_save_F77="$F77"
        F77="$MPIF77"
        AC_SUBST(MPIF77)

       #dd if test $acx_with_mpi != "no"; then
        if test -n "$MPICXX"; then
          AC_LANG_SAVE
          AC_LANG([C++])
          AC_CHECK_HEADERS([mpi.h], [acx_with_mpi=yes], 
                           [acx_with_mpi=no 
                            AC_MSG_ERROR(["Unable to include with mpi.h])])
          AC_LANG_RESTORE
        fi

        if test $acx_with_mpi != no; then
          AC_DEFINE([HAVE_MPI], [1], [Define if using mpi])
        fi
        AM_CONDITIONAL([HAVE_MPI], [test $acx_with_mpi != "no"])

])
