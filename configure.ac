AC_INIT([dyturbo], [1.0beta], [stefano.camarda@cern.ch,jakub.cuth@cern.ch])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

#Suppress verbose output when compiling
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

#By default, install binaries in the current directory
AC_PREFIX_DEFAULT(${PWD})

AC_PROG_CC
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

#Use libtoolize to build dyturbo libraries
AC_PROG_LIBTOOL
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# openmp suport
AC_ARG_ENABLE([openmp], AC_HELP_STRING(--enable-openmp, [Enable openMP compiler flag -fopenmp. (default=yes)])): ${enable_openmp=yes}
AS_IF([test "x$enable_openmp" = "xyes"], [
#AX_OPENMP #[, AC_MSG_ERROR([openmp macro not found])])
#AC_SUBST(OPENMP_CXXFLAGS)
#AC_SUBST(OPENMP_LDFLAGS)
    CXXFLAGS="${CXXFLAGS} -fopenmp"
    FFLAGS="${FFLAGS} -fopenmp"])

#optimisation flags O3
AC_ARG_ENABLE([O3], AC_HELP_STRING(--enable-O3, [Use compiler optimization flags -O3. (default=yes)])): ${enable_O3=yes}

AS_IF([test "x$enable_O3" = "xyes"], [
    CXXFLAGS="${CXXFLAGS} -O3"
    FFLAGS="${FFLAGS} -O3"])

#optimisation flags, O3 or Ofast
AC_ARG_ENABLE([Ofast],
  AC_HELP_STRING(--enable-Ofast, [Use unsafe compiler optimization flags -Ofast. (default=no)]))
  AS_IF([test "x$enable_Ofast" = "xyes"], [
    CXXFLAGS="${CXXFLAGS} -Ofast"
    FFLAGS="${FFLAGS} -Ofast"
])

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
   ${srcdir}/install-lhapdf
  #AC_MSG_ERROR([Unable to find lhapdf-config.])
  lhapdf_config=`ls ${PWD}/lhapdf6/bin/lhapdf-config`
  if test x$lhapdf_config == x; then
   # if the installation fails, exit
   AC_MSG_ERROR([Failed to install lhapdf.])
  fi
 fi
fi

# set up lhapdf flags
LHAPDF_CPPFLAGS=`${lhapdf_config} --cppflags`
LHAPDF_LDFLAGS=`${lhapdf_config} --ldflags`
AC_SUBST(LHAPDF_CPPFLAGS)
AC_SUBST(LHAPDF_LDFLAGS)
AC_MSG_RESULT([Using $lhapdf_config])

# check for Cuba installation, and install it if missing
AC_MSG_CHECKING([for Cuba library installation])
AC_SEARCH_LIBS([vegas],[cuba],[
	AC_MSG_RESULT([Using -lcuba from system installation])
        CUBA_LDFLAGS=-lcuba],
        [
	CUBA_VERSION=4.2
	AC_SUBST(CUBA_VERSION)
	# search for local istallation of Cuba-4.2
	cubalib=`ls ${PWD}/Cuba-${CUBA_VERSION}/libcuba.a`
	if test x$cubalib == x; then
          # if not found, install it
  	   AC_MSG_RESULT([Install Cuba ${CUBA_VERSION}])
	   ${srcdir}/install-cuba
	fi
	cubalib=`ls ${PWD}/Cuba-${CUBA_VERSION}/libcuba.a`
	if test x$cubalib == x; then
  	   AC_MSG_ERROR([Failed to install Cuba ${CUBA_VERSION}.])
	fi
	
        CUBA_LDFLAGS="-L${PWD}/Cuba-${CUBA_VERSION} -lcuba"
	CUBA_CPPFLAGS="-I${PWD}/Cuba-${CUBA_VERSION}"
	AC_MSG_RESULT([Using -lcuba from local installation ${PWD}/Cuba-${CUBA_VERSION}])
        ])

## Use the following to ship cuba as third party library
#CUBA_VERSION=4.2
#AC_SUBST(CUBA_VERSION)
#AC_CONFIG_SUBDIRS([cuba/Cuba-4.2])
#CUBA_LDFLAGS="-L\$(builddir)/../cuba/Cuba-${CUBA_VERSION} -lcuba"
#CUBA_CPPFLAGS="-I\$(srcdir)/../cuba/Cuba-4.2"


AC_SUBST(CUBA_CPPFLAGS)
AC_SUBST(CUBA_LDFLAGS)

# check for Chaplin installation, and install it if missing
AC_MSG_CHECKING([for chaplin library installation])
AC_SEARCH_LIBS([hpl2],[hpl3],[
	AC_MSG_RESULT([Using -lchaplin from system installation])
        CHAPLIN_LDFLAGS=-lchaplin],
        [
	CHAPLIN_VERSION=1.2
	AC_SUBST(CHAPLIN_VERSION)
	# search for local istallation of chaplin-1.2
	chaplinlib=`ls ${PWD}/chaplin-${CHAPLIN_VERSION}/lib/libchaplin.a`
	if test x$chaplinlib == x; then
          # if not found, install it
  	   AC_MSG_RESULT([Install chaplin ${CHAPLIN_VERSION}])
	   ${srcdir}/install-chaplin
	fi
	chaplinlib=`ls ${PWD}/chaplin-${CHAPLIN_VERSION}/lib/libchaplin.a`
	if test x$chaplinlib == x; then
  	   AC_MSG_ERROR([Failed to install chaplin ${CHAPLIN_VERSION}.])
	fi
	
        CHAPLIN_LDFLAGS="-L${PWD}/chaplin-${CHAPLIN_VERSION}/lib -lchaplin"
	AC_MSG_RESULT([Using -lchaplin from local installation ${PWD}/chaplin-${CHAPLIN_VERSION}])
        ])

AC_SUBST(CHAPLIN_LDFLAGS)

# check for cubature installation, and install it if missing
AC_MSG_CHECKING([for cubature source code])
${srcdir}/install-cubature

cubasrc=`ls ${PWD}/cubature/pcubature.c`
if test x$cubasrc == x; then
   AC_MSG_ERROR([Failed to install cubature.])
fi

# check for DYRES installation
AC_ARG_ENABLE([dyres], AC_HELP_STRING(--enable-dyres, [Install and link original DYRes code, Available only in git (default=no)]))
AS_IF([test "x$enable_dyres" = "xyes"], [
  AC_MSG_CHECKING([for DYRes installation])
  # search in origdyres/DYRes-v1.0
  dyresinst=`ls -d ${PWD}/origdyres/DYRes-v1.0`
  AS_IF([test "x$dyresinst" = x], [
    # if not found, install it
    AC_MSG_RESULT([Install DYRes])
    cd ${srcdir}/origdyres && ${srcdir}/install-dyres && cd ..
  ])
])
AM_CONDITIONAL([HAVE_DYRES],test "x$enable_dyres" = "xyes")

# enable test programs
AC_ARG_ENABLE([test], AC_HELP_STRING(--enable-test, [Compile test programs, Available only in git (default=no)]))
AS_IF([test "x$enable_test" = "xyes"], [
  GMOCK_OFILE="${PWD}/gmock/googlemock/make/gmock_main.a"
  AC_MSG_CHECKING([for GMock installation])
  gmockinst=`ls ${GMOCK_OFILE}  2> /dev/null`
  AS_IF([test "x$gmockinst" = x], [
         # if not found, install it
         AC_MSG_RESULT([ Gmock not found.. installing])
         [ ! test -d ${srcdir}/gmock/googlemock/make ] && git clone https://github.com/google/googletest.git ${srcdir}/gmock;
         cd ${srcdir}/gmock/googlemock/make && make ;
         cd -;
         ])
  GMOCK_CPPFLAGS="-pthread -isystem ${PWD}/gmock/googletest/include -isystem ${PWD}/gmock/googlemock/include"
  GMOCK_LDFLAGS="-lpthread ${GMOCK_OFILE}"
  AC_MSG_RESULT([ compiling test. Using gmock in ${GMOCK_OFILE}])
])
AM_CONDITIONAL([HAVE_TEST],test "x$enable_test" = "xyes")
AC_SUBST(GMOCK_CPPFLAGS)
AC_SUBST(GMOCK_LDFLAGS)

# Set --enable-root option, for linking histograms
AC_ARG_ENABLE([root],
  AC_HELP_STRING(--enable-root, [Use root for histograming. (default=no)]))
AS_IF([test "x$enable_root" = "xyes"], [
    AC_MSG_CHECKING([for root installation])
    root_config=`which root-config`
    AS_IF([test x$root_config == x],
       [AC_MSG_ERROR([Unable to find root-config.])],
       [AC_MSG_RESULT([Using $root_config])
       root_ok=1
       ROOT_CFLAGS=`root-config --cflags`
       ROOT_LDFLAGS=`root-config --ldflags --libs`
       AC_SUBST(ROOT_CFLAGS)
       AC_SUBST(ROOT_LDFLAGS)
       AC_DEFINE(USEROOT, "1",[ Turn on the ROOT for histogramming.] )
    ])
])
AM_CONDITIONAL([HAVE_ROOT],test "x$enable_root" = "xyes")

# Enable debug symbols
AC_ARG_ENABLE([debug],
  [AC_HELP_STRING(--enable-debug, [Add debug flags. (default=no)])])
AS_IF([test "x$enable_debug" = "xyes"], [
  AC_MSG_RESULT([Adding debug symbols.])
  FFLAGS="$FFLAGS -ggdb -gdwarf"
  CXXFLAGS="$CXXFLAGS -ggdb -gdwarf"
])

# Check for trapFPE flag
AC_ARG_ENABLE([trapFPE],
  [AC_HELP_STRING(--enable-trapFPE, [Stop of floating point errors (default=no)])])
AS_IF([test "x$enable_trapFPE" = "xyes"], [
  AC_MSG_RESULT([Add floating point trap.])
  FFLAGS="$FFLAGS -ffpe-trap='invalid,overflow,zero'"
])

# Check for check bounds flag
AC_ARG_ENABLE([checkBounds],
  [AC_HELP_STRING(--enable-checkBounds, [add -fbounds-check flag for compilation (default=no)])])
AS_IF([test "x$enable_checkBounds" = "xyes"], [
  AC_MSG_RESULT([Add check bounds flag])
  FFLAGS="$FFLAGS -fbounds-check"
])


AC_CONFIG_FILES([Makefile 
		 user/Makefile
		 histo/Makefile
 	 	 resum/Makefile
		 counterterm/Makefile
		 vjet/Makefile
		 v2j/Makefile
		 born/Makefile
		 phasespace/Makefile
		 integrands/Makefile
		 dyres/Makefile
		 dynnlo/Makefile
		 mcfm/Makefile
		 cernlib/Makefile
		 dequad/Makefile
		 sandia/Makefile
		 smolpack/Makefile
		 ancont/Makefile
		 pegasus/Makefile
		 cubature/Makefile
	 	 src/Makefile
		 dyturbo-config])


AM_COND_IF([HAVE_DYRES], [
  AC_CONFIG_FILES([origdyres/Makefile])
])
AM_COND_IF([HAVE_TEST], [
  AC_CONFIG_FILES([test/Makefile])
])

#add this to ship cuba
#		 cuba/Makefile

AC_OUTPUT
