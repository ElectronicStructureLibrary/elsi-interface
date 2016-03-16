AC_DEFUN([ACX_WITH_PARMETIS], [
  acx_with_parmetis="yes"
  AC_ARG_WITH([parmetis],
    [AS_HELP_STRING([--with-parmetis@<:@=Install DIR@:>@],
      [Enables use of ParMetis and Metis, required for PEXSI and SuperLU])],
    [
      case $withval in
      yes)
        acx_with_parmetis="yes"
      ;;
      no)
        acx_with_parmetis="no"
      ;;
      *)
        #CFLAGS="$FCFLAGS -I$withval/include"
        CPPFLAGS="$CFLAGS -I$withval/include"
        LIBSPARMETISLIB="-L$withval/lib  -lparmetis -lmetis "
        LIBS="$LIBSPARMETISLIB $LIBS"
        acx_with_parmetis="$withval"
      esac
    ],
    [acx_with_parmetis="no"]
  )
  if test $acx_with_parmetis != "no"; then
    AC_MSG_NOTICE( PARMETIS: checking for $p/$i )

    AC_LANG_SAVE
    AC_LANG_C
    AC_CHECK_HEADERS([parmetis.h], [], [acx_with_parmetis=no
                                        AC_MSG_NOTICE([Unable to compile with the parmetis.h])])
    AC_CHECK_LIB([parmetis], [libparmetis__SetupCtrl], [], [acx_with_parmetis=no
                                                                      AC_MSG_ERROR(["Unable to link with lib$with_parmetis_lib])], [-lmetis])
    AC_LANG_RESTORE
  fi
  if test $acx_with_parmetis != "no"; then
    AC_DEFINE(HAS_LIBPARMETIS, [1], [Define if should use libparmetis ])
  fi
  AM_CONDITIONAL([HAS_LIBPARMETIS], [test $acx_with_parmetis != "no"])
])

