## 
## This file is part of chebpol
##

AC_PREREQ(2.62)
AC_CONFIG_AUX_DIR([tools/autoconf])

AC_INIT([chebpol])
AC_CONFIG_SRCDIR([src/chebpol.c])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([src/Makevars])

AC_ARG_ENABLE([fftw],
       [AS_HELP_STRING([--disable-fftw],[do not use fftw])],
       [use_fftw="${enableval}"],[use_fftw=yes])

AC_ARG_ENABLE([alglib],
  [AS_HELP_STRING([--disable-alglib],[do not use alglib rbf])],
  [use_alglib="${enableval}"],[use_alglib=yes])

AC_ARG_ENABLE([check],
	[AS_HELP_STRING([--enable-check],[enable pedantic syntax check])],
	[use_check="${enableval}"],[use_check=no])

# some R standards
if test -n "${R_HOME}"; then
  CC=`${R_HOME}/bin/R CMD config CC`
  CPP=`${R_HOME}/bin/R CMD config CPP`
  CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
  CPPFLAGS=`${R_HOME}/bin/R CMD config CPPFLAGS`
fi

AS_IF([test "${use_check}" = "yes"], [
  EXTRA='-Wall -pedantic -Wextra'
])

AC_LANG([C])
AC_PROG_CC


AS_IF([test "x${use_alglib}" = xyes], [
 PKG_CHECK_MODULES([ALGLIB], [alglib >= 3.11], [AC_DEFINE(HAVE_ALGLIB)], [use_alglib=no])
])
AC_SUBST(ALGLIB_CFLAGS)
AC_SUBST(ALGLIB_LIBS)

AS_IF([test "x$use_fftw" = "xyes"], [
  PKG_CHECK_MODULES([FFTW], fftw3, [have_fftw=yes], 
                    [AC_MSG_WARN([FFTW not found. Compiling with matrix method])])
])

PKG_CHECK_MODULES([GSL], gsl, [AC_DEFINE(HAVE_GSL)],
			 [AC_MSG_WARN([GSL not found])])
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_SUBST(HAVE_GSL)

AS_IF([test "x$have_fftw" = "xyes" ], [
 # figure out if it works
  AC_MSG_CHECKING(for working fftw_plan_r2r in R)
  cat >conftest.c <<EOF
[
  #include <stdlib.h>
    #include <fftw3.h>
    int xxx_check() {
    const int dims[2] = {3,2};
    double vec[6];
    double out[6];
    fftw_r2r_kind kind[2] = {FFTW_REDFT10, FFTW_REDFT10};
    fftw_plan plan = fftw_plan_r2r(2,dims,vec,out,kind,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    if(plan != NULL) exit(1);
    plan = fftw_plan_r2r(1,dims,vec,out,kind,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    if(plan != NULL) exit(2);
    exit(3);
  }
]
EOF 
  ${R_HOME}/bin/R CMD SHLIB -o conftest.so conftest.c $FFTW_LIBS >& /dev/null
  echo "dyn.load('conftest.so'); .Call('xxx_check')" > conftest.R
  ${R_HOME}/bin/R --slave < conftest.R
  stat=$?
  AS_IF([test "x$stat" = "x1" ], [
    AC_MSG_RESULT(yes)
    AC_DEFINE(HAVE_FFTW)
  ],[
    AC_MSG_RESULT(no)
    AS_IF([test "x$stat" = "x2" ], [
      AC_DEFINE(HAVE_FFTW)
      AC_MSG_WARN([Your FFTW can't do 2d r2r-transforms, could it be MKL?])
      AC_MSG_CHECKING(whether static linking of FFTW solves it)
      TESTLIBS="-Wl,-Bstatic,-Bsymbolic-functions $FFTW_LIBS -Wl,-Bdynamic"
      touch conftest.c
      ${R_HOME}/bin/R CMD SHLIB -o conftest.so conftest.c $TESTLIBS >& /dev/null    
      ${R_HOME}/bin/R --slave >& /dev/null < conftest.R
      stat=$?
      AS_IF([test "x$stat" = "x1"], [AC_MSG_RESULT(yes); FFTW_LIBS=${TESTLIBS}], [AC_MSG_RESULT(no)])
    ], [
      FFTW_LIBS=
      FFTW_CFLAGS=
    ])
  ])
  rm -f conftest.o conftest.c conftest.so conftest.R
])

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
