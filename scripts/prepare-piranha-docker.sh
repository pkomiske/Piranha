#!/bin/bash

# install boost
curl -LO https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz
tar xf boost_1_76_0.tar.gz
rm -f boost_1_76_0.tar.gz
cd boost_1_76_0
./bootstrap.sh --prefix=/usr/local --with-libraries=math
./b2 install
cd ..
rm -rf boost_1_76_0

# install CGAL
yum install -y lzma gmp-devel mpfr-devel
curl -LO https://github.com/CGAL/cgal/releases/download/v5.2/CGAL-5.2-library.tar.xz
tar xf CGAL-5.2-library.tar.xz
rm -f CGAL-5.2-library.tar.xz
cd CGAL-5.2
/opt/python/cp39-cp39/bin/pip install cmake
/opt/python/cp39-cp39/bin/cmake -DCMAKE_BUILD_TYPE=Release .
make install
cd ..
rm -rf CGAL-5.2
