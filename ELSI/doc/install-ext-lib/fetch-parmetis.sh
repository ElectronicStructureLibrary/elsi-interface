#!/bin/sh


export EXTDIR=/tmp/install
export MPICC=mpicc
export MPICXX=mpicxx
export MAKE=make
export CC=$(MPICC)

export THIS_DIR=`pwd`
export CPU=` uname -m | sed "s/\\ /_/g"`
export SYS=` uname -s`

export package=parmetis-4.0.2.tar.gz

mkdir parmetis-build
cd parmetis-build
cp /home/vama/soft/chem2/elsi/parmetis-4.0.2.tar.gz .
#curl -O http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/$package
tar xvfz $package
mv parmetis*/* .
rm -rf parmetis* $package



echo Start building parMetis...
${MAKE} config prefix=${EXTDIR} cc=${MPICC} cxx=${MPICXX} cputype=${CPU} systype=${SYS}
cd ${THIS_DIR}/parmetis-build/parMetis/build/${SYS}-${CPU} ; make -j12 ; make install
mkdir -p ${EXTDIR}/include
mkdir -p ${EXTDIR}/lib
cp `find ${THIS_DIR}/parmetis-build/include -name parmetis.h` ${EXTDIR}/include
cp `find ${THIS_DIR}/parmetis-build/metis/include -name metis.h` ${EXTDIR}/include
cp `find ${THIS_DIR}/parmetis-build/build -name libparmetis.a` ${EXTDIR}/lib
cp `find ${THIS_DIR}/parmetis-build/build -name libmetis.a` ${EXTDIR}/lib
echo parMetis installed.

