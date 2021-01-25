#!/bin/bash

# Piranha depends on EventGeometry for Apollonius Subtraction
git clone https://github.com/pkomiske/EventGeometry.git
cd EventGeometry
git submodule update --init --recursive
make install
cd ..

# make shared library and put it in accessible location
make install_shared