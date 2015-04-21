#!/usr/bin/env bash
#
# A script to download and build pythia6.
#
# Build both a static library for use with fortran and a shared object 
# library for use with ROOT.  Strips out the dummy routines that would 
# hide the real ones found in pdflib.  Currently does *not* extract include
# files from the downloaded source, but uses definitions embedded in this
# script. Should work on both Linux and Mac OS X.
#
# Usage:  ./build_pythia6 [version] [--dummies=bestry|remove|keep] \
#                [--refetch] [gfortran|g77|g95] [-m32]
#
# where [version] takes a form like "6.4.24", "6424" or "6_424" w/ or
# without an optional leading "v".  It will created a subdirectory named 
# "v6_424" so it is probably wise to run this in a directory named "pythia6"
# or some such.  The --dummies=XYZZY controls whether the dummy PDFSET,
# STRUCTM and STRUCTP routines are removed.  The -m32 forces a 32-bit
# build.
#
# Creates directory structure:
#
# pythia6/                       - start in this directory
#         v6_424/                - created subdirectory
#                download/       - code downloaded from remote locations
#                inc/            - include files for selected common blocks
#                                  (hardcoded! not from downloaded source)
#                lib/            - liblund.a pydata.o libPythia6.so
#                src/            - fsplit source (and intermediate .o files)
#                tpythia6_build/ - files for building ROOT interface
#
# Questions? Suggestions? Improvements?  Offers of insane amounts of cash?
# Write me at: Robert Hatcher <rhatcher@fnal.gov>
#
# History:
#    2006-10-20:  original version created by <rhatcher@fnal.gov>
#    2007-06-28:  strict bash needs -gt to compare an int value, not ">"
#                 bail out if fetching source fails
#    2007-08-07:  version 6.4.11 and beyond are at www.hepforge.org and use
#                   a different naming convention; support 6409 and earlier
#                   from old location (6.4.10 is in tar form... punt for now)
#                 more alternatives for version specification, this is
#                   needed for different file naming schemes old vs. new
#                 use wget if possible, curl if not, bail if neither avail
#    2007-10-16:  add flags
#                   --keep-dummies : force the keeping of pdflib dummies
#                   --refetch : force refetch of source code
#                 look around for CERNLIB pdflib[804], mathlib, kernlib
#                 libraries with which to link against in order 
#                 to satisfy removed dummies
#    2007-10-17:  change to --dummmies=[besttry|remove|keep]
#                    besttry (default) means try to link to CERNLIB
#                      if any of the libraries can't be found then 
#                      don't remove dummmies
#                    remove always removes dummy routine and tries to link
#                      to what ever libraries that it can find
#                    keep always keeps dummy routine
#    2008-01-05:  fix "keep" option to actually keep dummy stubs.
#    2009-04-09:  choose fortran compiler, if cmd line specified use that
#                 otherwise first of gfortran or g77 found
#    2010-12-29:  addition of -m32 option; default v6_422
#    2011-01-13:  default v6_424; clean up comments
#
############################################################################
#
version=6.4.24
doclean=0
dummyaction="besttry"
refetch=0
whichftn="unknown"
m32flag=""
while [ $# -gt 0 ] ; do
  case $1 in
  *clean*)
     doclean=1 ;;
  *dumm*)
     case $1 in
     *keep*)
       dummyaction="keep" ;;
     *remove*)
       dummyaction="remove" ;;
     *best*)
       dummyaction="besttry" ;;
     *)
       echo "failed to parse $1" ;;
     esac ;;
  *refetch*)
     refetch=1 ;;
  *g77*)
     whichftn="g77" ;;
  *gfortran*)
     whichftn="gfortran" ;;
  *g95*)
     whichftn="g95" ;;
  *m32*)
     m32flag="-m32" ;;
  *)
     version=$1 ;;
  esac
  shift
done

############################################################################
#
# decompose version string which could take any of the forms:
#   [v]{Major}.{minor}.{tiny}, [v]{Major}_{minor}{tiny} or 
#   [v]{Major}{minor}{tiny}, with the assumption that {Major} is always "6"
#   and {tiny} is two digits.
#
# remove any "v" in version number
version0=${version#v*}

major=${version0%%.*}
if [ "${major}" = "6" ] ; then
# in M.m.tt format
  foo=${version0#*.}
  minor=${foo%%.*}
  tiny=${foo##*.}
else
  major=${version0%%_*}
  if [ "${major}" = "6" ] ; then
#   in M_mtt format, convert to Mmtt
    version0=${version0/_/}
  fi
# assume Mmtt format, could it also be Mmmtt?  anticipate it could be.
# ${parmeter:offset:length} doesn't work for negative offsets as advertized
  major=${version0:0:1}
  foo=`echo $version0 | cut -c2-`
  if [ ${foo} -lt 1000 ] ; then
    minor=`echo $foo | cut -c1`
    tiny=`echo $foo | cut -c2-3`
  else
    minor=`echo $foo | cut -c1-2`
    tiny=`echo $foo | cut -c3-4`
  fi
fi

############################################################################
#
# script only works for pythia6
#
if [ "${major}" != "6" ] ; then
  echo "sorry this script only works for pythia6, not ${major}"
  exit
fi

############################################################################
#
# use the naming conventions v{Major}_{minos}{tiny} for the subdir
# to match that used by UPS at FNAL.
#
topdir=v${major}_${minor}${tiny}
toppath=`pwd`/${topdir}

echo version $version major $major minor $minor tiny $tiny
echo full path: $toppath

if [ ! -d ${topdir} ] ; then
  mkdir $topdir
fi
cd $topdir
for subdir in inc lib download src tpythia6_build ; do
  if [ ! -d ${subdir} ] ; then
    mkdir ${subdir}
  fi
done

############################################################################
#
# pick a fortran:  g77 vs. gfortran vs. g95
#
if [ "$whichftn" == "unknown" ] ; then
# for now not g95 option ...
#  whichg95=`which g95 | grep -v "no g95 in"`
  if [ ! -z "${whichg95}" ] ; then 
     whichftn="g95"
  else
#    echo "no g95"
    whichgfortran=`which gfortran | grep -v "no gfortran in"`
    if [ ! -z "${whichgfortran}" ] ; then 
       whichftn="gfortran"
    else
#      echo "no gfortran"
      whichg77=`which g77 | grep -v "no g77 in"`
      if [ ! -z "${whichg77}" ] ; then 
         whichftn="g77"
      else
         echo "could not find a fortran compiler (g95, gfortran, g77)"
      fi
    fi
  fi
fi
export FORT=$whichftn
echo "using $FORT as fortran compiler"

############################################################################
#
# make "inc" symlink 
# (various scripts disagree about where the include files should live)
#
cd src
if [ ! -f inc ] ; then
   echo make symlink src/../inc src/inc
   ln -sf ../inc ./inc
else
   echo src/inc symlink exists
fi
cd ${toppath}
#
#
############################################################################
#
add2link()
{
 # routine for building up linkage "-Lpath -llibrary" in variable 
 # given by varname (arg 1) for library (arg 2) in potential locations 
 # (args 3+).  Don't repeat location in string.  Look for various extensions
 #
 varname=$1
 shift
 libname=$1
 shift
 export foundlib=0
 while [ $# -gt 0 ] ; do
   loc=$1
   shift
   for ext in .a .so .dylib ; do
     fullname=${loc}/lib${libname}${ext}
     #echo look for ${fullname}
     if [ -f ${fullname} ] ; then
        #echo hit it at -L${loc} -l${libname}
        foundlib=1
        echo ${!varname} | grep ${loc} 2>&1 > /dev/null
        if [ $? -ne 0 ] ; then
          export $varname="${!varname} -L${loc} -l${libname}"
        else
          export $varname="${!varname} -l${libname}"
        fi
        break;
     fi
   done
   if [ ${foundlib} -ne 0 ] ; then break ; fi
 done
}
#
############################################################################
#
# look for the CERNLIB libraries:  pdflib[804], mathlib, kernlib 
# these are needed in order to satisfy the removed dummy routines
# if they can't be found then don't remove the dummies
#
CERNLINK="" # -L${CERNLIBS} -lpdflib804 -lmathlib -lkernlib"
if [ "$dummyaction" = "keep" ] ; then
  echo requested to keep dummy PDFSET, STRUCTM and STRUCTP
else
  missinglib=""
  add2link CERNLINK pdflib $CERNLIB $CERNLIBS $CERN_DIR/lib $CERN_ROOT/lib 
  if [ ${foundlib} -eq 0 ] ; then
    add2link CERNLINK pdflib804 $CERNLIB $CERNLIBS $CERN_DIR/lib $CERN_ROOT/lib 
    if [ ${foundlib} -eq 0 ] ; then 
      missinglib="$missinglib pdflib[804]"
    fi
  fi
  add2link CERNLINK mathlib $CERNLIB $CERNLIBS $CERN_DIR/lib $CERN_ROOT/lib 
  if [ ${foundlib} -eq 0 ] ; then
    missinglib="$missinglib mathlib"
  fi
  add2link CERNLINK kernlib $CERNLIB $CERNLIBS $CERN_DIR/lib $CERN_ROOT/lib 
  if [ ${foundlib} -eq 0 ] ; then
    missinglib="$missinglib kernlib"
  fi
  # didn't find something?
  if [ -n "$missinglib" ] ; then
    echo "### unable to locate library:  ${missinglib} "
    if [ "$dummyaction" = "besttry" ] ; then
      echo "### revert to keeping dummy PDFSET, STRUCTM and STRUCTP "
      CERNLINK=""
      dummyaction="keep"
    fi
  else
    dummyaction="remove"
  fi
  if [ -n "$CERNLINK" ] ; then
    echo will attempt to resolve removed dummies with CERNLINK =
    echo "    $CERNLINK"
  fi
fi # user didn't initially request "keep"

############################################################################
#
if [ ${doclean} -ne 0 ] ; then
  echo clean out old code
fi

############################################################################
#
# decide on how to fetch the source files
#
cd ${toppath}/download
#
# wget or curl for retreiving remote files?
# (OS X doesn't generally have wget on it, so fall back on curl in that case)
#
whichfetchit=`which wget | grep -v "no wget in"`
if [ ! -z "${whichfetchit}" ] ; then
  echo use \"wget\" for fetching files
  fetchit='wget '
else
  whichfetchit=`which curl | grep -v "no curl in"`
  if [ ! -z "${whichfetchit}" ] ; then
    # -f = fail without creating dummy, -O output local named like remoteza
    echo use \"curl -f -O\" for fetching files
    fetchit='curl -f -O '
  else
    echo "Neither wget nor curl available -- can't download files"
    exit 1
  fi
fi

############################################################################
#
# retrieve the pythia6 fortran soruce file
# naming and location depends on version
#
mt=${minor}${tiny}
if [ $mt -lt 410 ] ; then
# old location, unzip .f file
  basef=pythia${major}${minor}${tiny}.f
  #location=http://www.thep.lu.se/~torbjorn/pythia
  location=http://home.thep.lu.se/~torbjorn/pythia
  gzipped=""
elif [ $mt -eq 410 ] ; then
# not a case we can handle
# there is a .tar.gz file at the new location .. but that's too much work
  echo "source files for minor version $minor can't be handled by this script"
  exit
else
# new location, .f file is gzipped
  basef=pythia-${major}.${minor}.${tiny}.f
  location=http://www.hepforge.org/archive/pythia6
  gzipped=".gz"
fi
# if we don't already have it, fetch the .f file
if [ ! -f ${basef}_with_dummies -o ${refetch} -ne 0 ] ; then
  echo "${fetchit} ${location}/${basef}${gzipped}"
  $fetchit ${location}/${basef}${gzipped}
  if [ ! ${basef}${gzipped} ] ; then
    echo "Sorry could not fetch ${basef}${gzipped} from ${location}"
    exit 1
  fi
  if [ -n "${gzipped}" ] ; then
    rm -f ${basef}
    gzip -d ${basef}${gzipped}
  fi
  mv ${basef} ${basef}_with_dummies
fi

############################################################################
#
# fetch the ROOT interface code
#
if [ ! -f pythia6.tar.gz -o ${refetch} -ne 0 ] ; then
  echo "${fetchit} ftp://root.cern.ch/root/pythia6.tar.gz"
  ${fetchit} ftp://root.cern.ch/root/pythia6.tar.gz
fi

############################################################################
#
# need to remove PDFSET/STRUCTM/STRUCTP routines from pythia6 code
# if we're using PDFLIB to supply the PDFSET, STRUCTM, STRUCTP routines
#
if [ "$dummyaction" = "keep" ] ; then
  cp ${basef}_with_dummies ${basef}
  echo not removing dummy routines
else
  cat > rm_dummy.awk <<EOF
BEGIN                { writeout=1; isend=0; }
/SUBROUTINE PDFSET/  { writeout=0; }
/SUBROUTINE STRUCTM/ { writeout=0; }
/SUBROUTINE STRUCTP/ { writeout=0; }
/      END\$/        { if ( writeout == 1 ) print \$0; writeout=1; isend=1; }
//                   { if ( writeout == 1 && isend != 1 ) print \$0; isend=0; }
EOF
  awk -f rm_dummy.awk ${basef}_with_dummies > ${basef}
fi

# extract block data routine
#cat > extract_pydata.awk <<EOF
#BEGIN { writeout=0; }
#/C...PYDATA/ { writeout=1; }
#/      END\$/        { if ( writeout == 1 ) print \$0; writeout=0; }
#//                   { if ( writeout == 1 ) print \$0; }
#EOF
#awk -f extract_pydata.awk ${basef} > pydata.f

############################################################################
#
# create a Makefile for the pythia6 source code
#
cd ${toppath}/src
cat > Makefile <<EOF
#
# simple pythia6 makefile
#
UNAME = \$(shell uname)
ifeq "\$(UNAME)" "Linux"
    AR=ar
    F77=$FORT
    FFLAG= -O -fno-second-underscore -fPIC $m32flag
    CPP = gcc -E 
    CPPFLG= -C -P
endif
ifeq "\$(UNAME)" "Darwin"
    AR=ar
    F77=$FORT
    FFLAG= -O -fno-second-underscore -fPIC $m32flag
    CPP = cc -E
    CPPFLG= -C -P
endif
LIB = ../lib
CPPFLGS = \$(CPPFLG) -D\$(UNAME)

FOBJSALL = \$(patsubst %.f,%.o,\$(wildcard *.f)) \
           \$(patsubst %.F,%.o,\$(wildcard *.F))

# remove the "pdfdum.o" as we don't want that going into the .a library
FOBJS = \$(filter-out pdfdum.o,\$(FOBJSALL))

#------------------------------------------

all: \$(LIB)/liblund.a \$(LIB)/pydata.o 

\$(LIB)/liblund.a: \$(FOBJS) 
	\$(AR) -urs \$(LIB)/liblund.a \$(FOBJS) 

\$(LIB)/pydata.o: pydata.o
	cp -p pydata.o \$(LIB)/pydata.o

\$(LIB)/pdfdum.o: pdfdum.o
	cp -p pdfdum.o \$(LIB)/pdfdum.o

clean:
	rm -f *.o

veryclean:
	rm -f \$(LIB)/liblund.a \$(LIB)/pydata.o \$(LIB)/pdfdum.o
	rm -f *.o

#------------------------------------------

.SUFFIXES : .o .f .F

.f.o:
	\$(F77) \$(FFLAG) -c \$<

.F.o: 
	\$(F77) \$(FFLAG) -c \$<

EOF

############################################################################
#
# split the pythia6 source code, build the library
#
whichfsplit=`which fsplit | grep -v "no fsplit in"`
echo initial ${whichfsplit}
if [ -z "${whichfsplit}" ] ; then
  echo "No 'fsplit' -- can't build library without it"
  echo "try to build fsplit"
  echo "$fetchit -O http://home.fnal.gov/~rhatcher/fsplit.c"
  $fetchit http://home.fnal.gov/~rhatcher/fsplit.c
  gcc -o fsplit fsplit.c
  PATH=.:${PATH}
  whichfsplit=`which fsplit`
  echo "built ${whichfsplit}"
fi
if [ ! -z "${whichfsplit}" ] ; then
  echo "build static lund library"
  if [ ${doclean} -gt 0 ] ; then rm -f *.f ; fi
  if [ ! -f pydata.f ] ; then 
    fsplit ${toppath}/download/${basef}
# hack for illegal computation in DATA statement in some versions
# apparently works for g77 but not gfortran
    mv pyalps.f pyalps.f_original
    sed -e 's+(2D0\*107D0/2025D0)+0.10567901234568D0+g' \
        -e 's+(2D0\*963D0/14375D0)+0.13398260869565D0+g' \
        -e 's+(2D0\*321D0/3703D0)+0.17337294085876D0+g' \
        -e 's+(-2D0\*107D0/1875D0)+-0.11413333333333D0+g' \
        -e 's+(-2D0\*963D0/13225D0)+-0.14563327032136D0+g' \
        -e 's+(-2D0\*321D0/3381D0)+-0.18988464951198D0+g' \
      pyalps.f_original > pyalps.f
  fi
  rm -f zzz*.f
  if [ "$dummyaction" != "keep" ] ; then
    for junk in pdfset.f structm.f structp.f ; do
      if [ -f ${junk} ] ; then
        echo "${junk} in split output, should have been removed"
        rm -f ${junk}
      fi
    done
  fi
  make all
fi

############################################################################
#
# build the ROOT interface library libPythia6.so 
#
echo "build ROOT accessable shared library"
cd ${toppath}/tpythia6_build

tar xzvf ${toppath}/download/pythia6.tar.gz pythia6/tpythia6_called_from_cc.F
tar xzvf ${toppath}/download/pythia6.tar.gz pythia6/pythia6_common_address.c
mv pythia6/* .
rmdir pythia6
echo 'void MAIN__() {}' > main.c
gcc -c -fPIC $m32flag main.c
gcc -c -fPIC $m32flag pythia6_common_address.c
$FORT -c -fPIC -fno-second-underscore $m32flag tpythia6_called_from_cc.F

cd ${toppath}/lib

arch=`uname`
if [ ${arch} == "Linux" ] ; then
  $FORT $m32flag -shared -Wl,-soname,libPythia6.so -o libPythia6.so \
    ${toppath}/tpythia6_build/*.o ${toppath}/src/*.o ${CERNLINK}
fi
if [ ${arch} == "Darwin" ] ; then
  macosx_minor=`sw_vers | sed -n 's/ProductVersion://p' | cut -d . -f 2`
  gcc $m32flag -dynamiclib -flat_namespace -single_module -undefined dynamic_lookup \
      -install_name ${toppath}/lib/libPythia6.dylib -o libPythia6.dylib \
      ${toppath}/tpythia6_build/*.o ${toppath}/src/*.o ${CERNLINK}
  if [ $macosx_minor -ge 4 ]; then
     ln -sf libPythia6.dylib libPythia6.so
  else
     gcc $m32flag -bundle -flat_namespace -undefined dynamic_lookup -o libPythia6.so \
      ${toppath}/tpythia6_build/*.o ${toppath}/src/*.o ${CERNLINK}
  fi
fi

############################################################################
#
# create pythia6 include files for common blocks
#
# here we'd like to automate the extraction of the common common blocks
# into include files but it's non-trivial

echo "extract include files"
cd ${toppath}/inc

# define a #include that declares the types of all the pythia functions
cat > pyfunc.inc <<EOF
C...standard pythia functions
      double precision PYFCMP,PYPCMP
      double precision PYCTEQ,PYGRVV,PYGRVW,PYGRVS,PYCT5L,PYCT5M,PYHFTH
      double precision PYGAMM,PYSPEN,PYTBHS,PYRNMQ,PYRNM3,PYFINT,PYFISB
      double precision PYXXZ6,PYXXGA,PYX2XG,PYX2XH,PYH2XX
      double precision PYGAUS,PYGAU2,PYSIMP,PYLAMF,PYTHAG
      double precision PYRVSB,PYRVI1,PYRVI2,PYRVI3,PYRVG1,PYRVG2,PYRVG3
      double precision PYRVG4,PYRVR, PYRVS, PY4JTW,PYMAEL
      double precision PYMASS,PYMRUN,PYALEM,PYALPS,PYANGL
      double precision PYR,   PYP
      integer          PYK,PYCHGE,PYCOMP
      character*40     VISAJE

      external         PYFCMP,PYPCMP
      external         PYCTEQ,PYGRVV,PYGRVW,PYGRVS,PYCT5L,PYCT5M,PYHFTH
      external         PYGAMM,PYSPEN,PYTBHS,PYRNMQ,PYRNM3,PYFINT,PYFISB
      external         PYXXZ6,PYXXGA,PYX2XG,PYX2XH,PYH2XX
      external         PYGAUS,PYGAU2,PYSIMP,PYLAMF,PYTHAG
      external         PYRVSB,PYRVI1,PYRVI2,PYRVI3,PYRVG1,PYRVG2,PYRVG3
      external         PYRVG4,PYRVR, PYRVS, PY4JTW,PYMAEL
      external         PYMASS,PYMRUN,PYALEM,PYALPS,PYANGL
      external         PYR,   PYP
      external         PYK,PYCHGE,PYCOMP
      external         VISAJE
EOF

# how to "automate" the others ... including declaring the types
# without using the IMPLICIT DOUBLE PRECISION etc.
# NEUGEN3 needs pydat1.inc pydat3.inc pyjets.inc at a minimum
# !!!! for now just hard-code them !!!

cat > pydat1.inc <<EOF
C...Parameters.
      integer       MSTU,               MSTJ
      double precision        PARU,               PARJ
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
EOF

cat > pydat2.inc <<EOF
C...Particle properties + some flavour parameters.
      integer       KCHG
      double precision          PMAS,       PARF,      VCKM
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
EOF

cat > pydat3.inc <<EOF
C...Decay information
      integer       MDCY,       MDME,                   KFDP
      double precision                       BRAT
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
EOF

cat > pyjets.inc <<EOF
C...The event record.
      integer       N,NPAD,K
      double precision               P,        V
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      SAVE  /PYJETS/
EOF

cat > pypars.inc <<EOF
C...Parameters.
      integer       MSTP,               MSTI
      double precision        PARP,               PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
EOF

cat > pysubs.inc <<EOF
C...Selection of hard scattering subprocesses
      integer       MSEL,MSELPD,MSUB,     KFIN
      double precision                                   CKIN
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
EOF

############################################################################
#
# done
#
echo "end-of-script $0"
# End-of-Script
