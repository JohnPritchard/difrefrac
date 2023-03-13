dnl "@(#) $Id: configure.in 3021 2017-02-01 08:23:08Z jpritcha $"
dnl Process this file with autoconf to produce a configure script.
AC_INIT
AC_CONFIG_SRCDIR([src/difrefrac.cc])

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_RANLIB

dnl Checks for libraries.
AC_CHECK_LIB(m, main)
AC_CHECK_LIB(pgplot, main)
AC_CHECK_LIB(cpgplot, main)
AC_CHECK_LIB(skycalc, main)

dnl Checks for header files.
m4_warn([obsolete],
[The preprocessor macro `STDC_HEADERS' is obsolete.
  Except in unusual embedded environments, you can safely include all
  ISO C90 headers unconditionally.])dnl
# Autoupdate added the next two lines to ensure that your configure
# script's behavior did not change.  They are probably safe to remove.
AC_CHECK_INCLUDES_DEFAULT
AC_PROG_EGREP

AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(cpgplot.h)
AC_CHECK_HEADERS(libskycalc.h, have_libskycalch=yes, have_libskycalch=no )

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_STRUCT_TM

dnl Checks for library functions.
AC_CHECK_FUNCS(mktime)

AC_SUBST(have_libskycalch)

AC_CONFIG_FILES([Makefile libsrc/Makefile src/Makefile])
AC_OUTPUT
