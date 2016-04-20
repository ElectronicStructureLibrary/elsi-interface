#!/bin/sh

export EXTDIR=/tmp/install-this-an-example
export MPICC=mpicc
export MPIFC=mpif90
export MPICXX=mpicxx
export MAKE=make

export CC=${MPICC}
export THIS_DIR=`pwd`
export CPU=` uname -m | sed "s/\\ /_/g"`
export SYS=` uname -s`
export AR=ar

mkdir superlu-build
cd superlu-build
package=superlu_dist_4.3.tar.gz

curl -O http://crd-legacy.lbl.gov/~xiaoye/SuperLU/$package
tar xvfz $package
mv SuperLU*/* .
rm -rf SuperLU* $package


echo Start building SuperLU...
cp ${THIS_DIR}/make.inc.superlu ${THIS_DIR}/superlu-build/make.inc
echo Start building SuperLU...
cd  ${THIS_DIR}/superlu-build 
make -j1
mkdir -p ${EXTDIR}/include
mkdir -p ${EXTDIR}/lib
cp  ${THIS_DIR}/superlu-build/SRC/*.h ${EXTDIR}/include
cp  ${THIS_DIR}/superlu-build/lib/libsuperlu.a ${EXTDIR}/lib
echo SuperLU installed.


