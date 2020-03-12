#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.56)

AC_INIT(USF xanalysis, 6.11.11, dshuman@usf.edu, xanalysis)
## v4.0.0 uses increased SPIKETIMES array (200,512k)
## v6.5.0 marks the shift from f77 to f90

AM_INIT_AUTOMAKE

## EXAMPLE
## AC_CONFIG_SRCDIR([autoCCH.defs])
AC_CONFIG_SRCDIR([x2000parameter.defs])

AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC(gcc)
AM_PROG_CC_C_O
if test "$CC" = gcc; then
   AC_SUBST(GCCLIB,[`gcc -print-libgcc-file-name`])
fi
if test "$CC" = gfortran; then
   AC_SUBST(GCCLIB,[`gcc -print-libgcc-file-name`])
fi
AC_PROG_F77(f90 gfortran)
AM_CONDITIONAL(F90, [test x$F77 = xf90])
 
# Checks for libraries.
AC_F77_LIBRARY_LDFLAGS

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_SIZEOF(void *)
SIZEOF_VOID_P=$ac_cv_sizeof_void_p
AC_SUBST(SIZEOF_VOID_P)
AC_CHECK_SIZEOF(long)
AC_C_BIGENDIAN()

# Checks for library functions.
AC_REPLACE_FUNCS([asprintf])
AC_CONFIG_FILES([Makefile gopen_type.defs])
AC_OUTPUT
