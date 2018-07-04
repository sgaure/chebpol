## 
## This file is part of chebpol
##
AC_PREREQ(2.62)
AC_CONFIG_AUX_DIR([src/build/autoconf])

AC_INIT([chebpol])
AC_CONFIG_SRCDIR([src/chebpol.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])

AC_ARG_ENABLE([fftw],
       [AS_HELP_STRING([--disable-fftw],[do not use fftw])],
       [use_fftw=no],[use_fftw=yes])

AC_ARG_ENABLE([alglib],
  [AS_HELP_STRING([--disable-alglib],[use alglib rbf])],
  [use_alglib=${enableval}],[use_alglib=yes])

AC_ARG_ENABLE([check],
	[AS_HELP_STRING([--enable-check],[enable pedantic syntax check])],
	[EXTRA='-Wall -pedantic -Wextra'])

# some R standards
if test -n "${R_HOME}"; then
  CC=`${R_HOME}/bin/R CMD config CC`
  CPP=`${R_HOME}/bin/R CMD config CPP`
  CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
  CPPFLAGS=`${R_HOME}/bin/R CMD config CPPFLAGS`
fi

if test "x$use_fftw" = "xyes"; then
  PKG_CHECK_MODULES([FFTW], fftw3, [have_fftw=yes], [AC_MSG_WARN([FFTW not found. Compiling with matrix method])])
fi

if test "${use_alglib}" = yes; then
 # figure out if it works
 PKG_CHECK_MODULES([ALGLIB], alglib, [AC_DEFINE(HAVE_ALGLIB)])
fi
AC_SUBST(ALGLIB_CFLAGS)
AC_SUBST(ALGLIB_LIBS)

AC_LANG([C])
AC_PROG_CC


if test "x$have_fftw" = "xyes" ; then
  AC_MSG_CHECKING(for working fftw_plan_r2r in R)
  echo ['#include <stdlib.h>
    #include <fftw3.h>
    int xxx_check() {
    const int dims[2] = {2,3};
    double vec[6];
    double out[6];
    fftw_r2r_kind kind[2] = {FFTW_REDFT10, FFTW_REDFT10};
    fftw_plan plan = fftw_plan_r2r(2,dims,vec,out,kind,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    if(plan == NULL) exit(1); else exit(13);
  }'] > conftest.c
  ${R_HOME}/bin/R CMD SHLIB -o conftest.so conftest.c $FFTW_LIBS >& /dev/null
  echo "dyn.load('conftest.so'); .Call('xxx_check')" > conftest.R
  ${R_HOME}/bin/R --slave < conftest.R
  stat=$?
  if test "x$stat" == "x13" ; then
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_FFTW)
  else
    AC_MSG_RESULT(no)
    AC_MSG_WARN([Your FFTW can't do 2d r2r-transforms, could it be MKL?])
    AC_MSG_CHECKING(static linking of FFTW)
    FFTW_LIBS="-Wl,-Bsymbolic-functions,-Bstatic,$FFTW_LIBS,-Bdynamic"
    touch conftest.c
    ${R_HOME}/bin/R CMD SHLIB -o conftest.so conftest.c $FFTW_LIBS >& /dev/null    
    ${R_HOME}/bin/R --slave >& /dev/null < conftest.R
    stat=$?
    if test "x$stat" == "x13"; then
      AC_MSG_RESULT(yes)
      AC_DEFINE(HAVE_FFTW)
    else
      AC_MSG_RESULT([no, compiling with slow matrix method])
      FFTW_LIBS=
      FFTW_CFLAGS=
    fi
  fi
  rm -f conftest.o conftest.c conftest.so conftest.R
fi

## --- Output ----------------------------------------------------------------
AC_SUBST(OPENMP_CFLAGS)
AC_SUBST(HAVE_FFTW)
AC_SUBST(EXTRA)
AC_SUBST(FFTW)
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)

AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)

AC_OUTPUT
