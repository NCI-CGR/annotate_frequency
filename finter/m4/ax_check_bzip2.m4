# =============================================================================
#  Derived from https://www.gnu.org/software/autoconf-archive/ax_boost_program_options.html
# =============================================================================
#
# SYNOPSIS
#
#   AX_CHECK_BZIP2
#
# DESCRIPTION
#
#   Test for bzip2 (bz2) library. 
#
#   This macro calls:
#
#     AC_SUBST(BZIP2_LIB)
#
#   And sets:
#
#     HAVE_BZIP2
#
# LICENSE
#
#   Copyright (c) 2021 Cameron Palmer <palmercd@nih.gov>
#   Copyright (c) 2009 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 26

AC_DEFUN([AX_CHECK_BZIP2],
[
	AC_ARG_WITH([bzip2],
		AS_HELP_STRING([--with-bzip2@<:@=special-lib@:>@],
                       [use the bzip2 library - it is possible to specify a certain library for the linker
                        e.g. --with-bzip2=bzip2-gcc-mt-1_33_1 ]),
        [
        if test "$withval" = "no"; then
			want_bzip2="no"
        elif test "$withval" = "yes"; then
            want_bzip2="yes"
            ax_bzip2_user_lib=""
        else
		    want_bzip2="yes"
		ax_bzip2_user_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_bzip2" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
	    export want_bzip2
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS
		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS
		AC_CACHE_CHECK([whether the bzip2 library is available],
					   ax_cv_bzip2,
					   [AC_LANG_PUSH(C++)
				AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <bzlib.h>
                                                          ]],
                                  [[BZFILE *_bz_input = 0;
                                   return 0;]])],
                           ax_cv_bzip2=yes, ax_cv_bzip2=no)
					AC_LANG_POP([C++])
		])
		if test "$ax_cv_bzip2" = yes; then
				AC_DEFINE(HAVE_BZIP2,,[define if the bzip2 library is available])
                  BZIP2LIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
                if test "x$ax_bzip2_user_lib" = "x"; then
                for libextension in `ls $BZIP2LIBDIR/libbz2*.so* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(bz2.*\)\.so.*$;\1;'` `ls $BZIP2LIBDIR/libbz2*.dylib* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(bz2.*\)\.dylib.*$;\1;'` `ls $BZIP2LIBDIR/libbz2*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(bz2.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BZIP2_LIB="-l$ax_lib"; AC_SUBST(BZIP2_LIB) link_bzip2="yes"; break],
                                 [link_bzip2="no"])
				done
                if test "x$link_bzip2" != "xyes"; then
                for libextension in `ls $BZIP2LIBDIR/bzip2*.dll* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(bzip2.*\)\.dll.*$;\1;'` `ls $BZIP2LIBDIR/bzip2*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(bzip2.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BZIP2_LIB="-l$ax_lib"; AC_SUBST(BZIP2_LIB) link_bzip2="yes"; break],
                                 [link_bzip2="no"])
				done
                fi
                else
                  for ax_lib in $ax_bzip2_user_lib bzip2-$ax_bzip2_user_lib; do
				      AC_CHECK_LIB($ax_lib, main,
                                   [BZIP2_LIB="-l$ax_lib"; AC_SUBST(BZIP2_LIB) link_bzip2="yes"; break],
                                   [link_bzip2="no"])
                  done
                fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_ERROR(Could not find a version of the bzip2 library!)
            fi
				if test "x$link_bzip2" != "xyes"; then
					AC_MSG_ERROR([Could not link against [$ax_lib] !])
				fi
		fi
		CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
	fi
])
