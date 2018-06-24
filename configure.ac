## 
## This file is part of chebpol
##
AC_PREREQ(2.61)
AC_CONFIG_AUX_DIR([src/build/autoconf])

## ---[ VERSION ]-------------------------------------------------------------

AC_INIT()


AC_ARG_ENABLE([fftw],
       [AS_HELP_STRING([--disable-fftw],[do not use fftw])],
       [use_fftw=no],[use_fftw=yes])

AC_ARG_ENABLE([check],
	[AS_HELP_STRING([--enable-check],[enable pedantic syntax check])],
	[EXTRA='-Wall -pedantic'; echo "PEDANTIC CHECK"])
# some R 
if test -n "${R_HOME}"; then
  CC=`${R_HOME}/bin/R CMD config CC`
  CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
  CPPFLAGS=`${R_HOME}/bin/R CMD config CPPFLAGS`
fi

AC_PROG_CC

## --- FFTW 3 ----------------------------------------------------------------
AS_IF([test "x$use_fftw" = "xyes"], [
  PKG_CHECK_MODULES([FFTW], fftw3, [FFTW=-DHAVE_FFTW], [echo '*** FFTW not found. Compiling with matrix method'])
])
## --- Output ----------------------------------------------------------------
AC_SUBST(EXTRA)
AC_SUBST(FFTW)
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)

AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)

AC_OUTPUT(src/Makevars)
