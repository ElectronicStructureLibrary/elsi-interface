#!/bin/sh

package=superlu_dist_3.3.tar.gz

curl -O http://crd-legacy.lbl.gov/~xiaoye/SuperLU/$package
tar xvfz $package
mv SuperLU*/* .
rm -rf SuperLU* $package
