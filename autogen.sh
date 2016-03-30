#!/bin/sh
# Run this to generate all the initial makefiles, etc.

set -e
mkdir -p m4/

test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.


aclocal -I ./config
autoconf
autoheader
#if hash glibtoolize 2>/dev/null; then
#  glibtoolize
#else
#  libtoolize
#fi
automake --add-missing
sh external/ELPA/autogen.sh
