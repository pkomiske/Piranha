#!/bin/bash

# homebrew packages
brew update > /dev/null
brew install libomp
brew upgrade boost cgal

# install fastjet
git clone https://gitlab.com/pkomiske/fastjet.git
cd fastjet
git submodule update --init --recursive
autoreconf -i
export CXXFLAGS=-std=c++14
./configure --prefix=/usr/local --enable-pyext --enable-cgal --enable-cgal-header-only --disable-monolithic --disable-allplugins --disable-debug PYTHON=python3 PYTHON_CONFIG=python3-config
make -j2 install
cd ..

# Piranha depends on EventGeometry for Apollonius Subtraction
git clone https://github.com/pkomiske/EventGeometry.git
cd EventGeometry
git submodule update --init --recursive
make install
cd ..

# make shared library and put it in accessible location
make shared
cp lib*.dylib /usr/local/lib