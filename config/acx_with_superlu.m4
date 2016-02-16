AC_DEFUN([ACX_WITH_SUPERLU], [
  acx_with_superlu=no
  AC_ARG_WITH([superlu],
    [AS_HELP_STRING([--with-superlu@<:@=DIR@:>@], [Build with SuperLu headers.])],
    [
      case $withval in
      yes)
        acx_with_superlu=yes
      ;;
      no)
        acx_with_superlu=no
      ;;
      *)
        acx_with_superlu=yes
        CPPFLAGS="-I$withval/include $CPPFLAGS"
        SUPERLU_LIBS="-L$withval/lib -lsuperlu "
        LIBS="$SUPERLU_LIBS $LIBS"
      ;;
      esac
    ]
  )

  if test "$acx_with_superlu" != no; then
    AC_MSG_NOTICE( SuperLU: checking for $p/$i )
    AC_LANG_SAVE
    AC_LANG_C
    # Check for the pressence of Superlu header files.
    AC_CHECK_HEADER([superlu_defs.h], [],
      [AC_MSG_ERROR([Unable to find the superlu_defs.h  header file.])])
    AC_LANG_RESTORE
  fi

  if test "$acx_with_superlu" != no; then
    AC_DEFINE([HAS_SUPERLU], [1], [Define if should use  SuperLU  operations])
  fi

  AM_CONDITIONAL([HAS_SUPERLU], [test $acx_with_superlu != "no"])
])

