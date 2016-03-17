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
        CFLAGS="$CFLAGS -I$withval/include"
        CXXFLAGS="$CPPFLAGS -I$withval/include"
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

    AC_CHECK_HEADER([parmetis.h], [], 
      [AC_MSG_ERROR([Unable to find the parmetis.h  header file.    ])])
    AC_CHECK_LIB([parmetis], [libparmetis__SetupCtrl], [], 
      [AC_MSG_ERROR([Unable to link with parmetis library])], [])
   AC_LANG_RESTORE
     AC_DEFINE(HAS_LIBPARMETIS, [1], [Define if should use libparmetis ])
  fi
  if test $acx_with_parmetis = "no"; then
    AC_MSG_ERROR([Parmetis is not defined ])
  fi
])

