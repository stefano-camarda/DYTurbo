AC_INIT([dyturbo], [0.9.1], [stefano.camarda@cern.ch,jakub.cuth@cern.ch])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

#Suppress verbose output when compiling
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

#By default, install binaries in the current directory
AC_PREFIX_DEFAULT(${PWD})

AC_PROG_CC
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

#Use libtoolize to build dyturbo libraries
#AM_PROG_LIBTOOL
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

#optimisation flags O3
AC_ARG_ENABLE([O3],
  [AC_HELP_STRING(--enable-O3, [Use compiler optimization flags -O3. (default=no)])])
if test x$enable_O3 == xyes
then
  CXXFLAGS="${CXXFLAGS} -O3"
  FFLAGS="${FFLAGS} -O3"
fi

#optimisation flags, O3 or Ofast
AC_ARG_ENABLE([Ofast],
  [AC_HELP_STRING(--enable-Ofast, [Use unsafe compiler optimization flags -Ofast. (default=no)])])
if test x$enable_Ofast == xyes
then
  CXXFLAGS="${CXXFLAGS} -Ofast"
  FFLAGS="${FFLAGS} -Ofast"
fi

#c++11 flag is needed for the random number generator function
#CXXFLAGS="${CXXFLAGS} -std=c++11"

#older compilers needs c++0x instead
CXXFLAGS="${CXXFLAGS} -std=c++0x"

#use intel compilers
#AC_ARG_ENABLE([ifort],
#  [AC_HELP_STRING(--enable-ifort, [Use ifort2015 and icc2015 compilers. (default=no)])])
#if test x$enable_ifort == xyes
#then
#  AC_PROG_CXX([icc2015])
#  AC_PROG_F77([ifort2015])
#  AC_PROG_FC([ifort2015])
#else  
##use gfortran compilers
# AC_PROG_CXX
# AC_PROG_F77([gfortran])
# AC_PROG_FC
#fi

AC_PROG_CXX
AC_PROG_F77([gfortran])
AC_PROG_FC

# check for lhapdf installation
AC_MSG_CHECKING([for lhapdf installation])
# search lhapdf6 in the path
lhapdf_config=`which lhapdf-config`
if test x$lhapdf_config == x; then
# if not found, search lhapdf6 in the lhapdf6 subdirectory
 lhapdf_config=`ls ${PWD}/lhapdf6/bin/lhapdf-config`
 if test x$lhapdf_config == x; then
# if not found, install it
  AC_MSG_RESULT([Install LHAPDF6])
  ./install-lhapdf
  #AC_MSG_ERROR([Unable to find lhapdf-config.])
  lhapdf_config=`ls ${PWD}/lhapdf6/bin/lhapdf-config`
 fi
fi

# set up lhapdf flags
LHAPDF_CPPFLAGS=`${lhapdf_config} --cppflags`
LHAPDF_LDFLAGS=`${lhapdf_config} --libdir`/libLHAPDF.a
AC_SUBST(LHAPDF_CPPFLAGS)
AC_SUBST(LHAPDF_LDFLAGS)
AC_MSG_RESULT([Using $lhapdf_config])

# check for Cuba installation
AC_MSG_CHECKING([for Cuba library installation])
# search in Cuba-4.2
cubalib=`ls ${PWD}/Cuba-4.2/libcuba.a`
if test x$cubalib == x; then
    # if not found, install it
    AC_MSG_RESULT([Install Cuba])
    ./install-cuba
fi

AC_MSG_RESULT([Using ${PWD}/Cuba-4.2 ])
CUBA_CPPFLAGS="-I${PWD}/Cuba-4.2"
CUBA_LDFLAGS="-L${PWD}/Cuba-4.2 -lcuba"
AC_SUBST(CUBA_CPPFLAGS)
AC_SUBST(CUBA_LDFLAGS)
cubadir="${PWD}/Cuba-4.2"
AC_MSG_RESULT([Using $cubadir])


#AC_SEARCH_LIBS([vegas],[cuba],[
#	AC_MSG_RESULT([Using -lcuba])
#        CUBA_LDFLAGS=-lcuba],
#        [./install-cuba
#         CUBA_LDFLAGS=-L`pwd`/Cuba-4.2 -lcuba
#	AC_MSG_RESULT([Using -lcuba])
#        ])
#AC_SUBST(CUBA_LDFLAGS)

# Set --enable-root option, for linking histograms
AC_ARG_ENABLE([root],
  [AC_HELP_STRING(--enable-root, [Use root for histograming. (default=no)])])
if test x$enable_root == xyes
then
    AC_MSG_CHECKING([for root installation])
    root_config=`which root-config`
    if test x$root_config == x;
    then
        AC_MSG_ERROR([Unable to find root-config.])
    else
        AC_MSG_RESULT([Using $root_config])
        root_ok=1
        ROOT_CFLAGS=`root-config --cflags`
        ROOT_LDFLAGS=`root-config --ldflags`
        ROOT_LIBS=`root-config --glibs`
        AC_SUBST(ROOT_CFLAGS)
        AC_SUBST(ROOT_LDFLAGS)
        AC_SUBST(ROOT_LIBS)
        AC_DEFINE(USEROOT, "1",[ Turn on the ROOT for histogramming.] )
    fi
fi

AC_ARG_ENABLE([debug],
  [AC_HELP_STRING(--enable-debug, [Add debug flags. (default=no)])])
if test x$enable_debug == xyes; then
  echo "Adding debug symbols."
  FFLAGS="$FFLAGS -ggdb"
  CXXFLAGS="$CXXFLAGS -ggdb"
fi


# Check for trapFPE flag
AC_ARG_ENABLE([trapFPE],
  [AC_HELP_STRING(--enable-trapFPE, [Stop of floating point errors (default=no)])])
if test x$enable_trapFPE == xyes; then
  echo "Add floating point trap"
  FFLAGS="$FFLAGS -ffpe-trap='invalid,overflow,zero' "
fi

# Check for check bounds flag
AC_ARG_ENABLE([checkBounds],
  [AC_HELP_STRING(--enable-checkBounds, [add -fbounds-check flag for compilation (default=no)])])
if test x$enable_checkBounds == xyes; then
  echo "Add check bounds flag"
  FFLAGS="$FFLAGS -fbounds-check"
fi


AC_CONFIG_FILES([Makefile src/Makefile dyres/Makefile dynnlo/Makefile cernlib/Makefile dyturbo-config])
AC_OUTPUT
