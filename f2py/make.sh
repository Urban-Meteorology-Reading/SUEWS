#!/usr/bin/env bash
# script for compiling pythonic SUEWS using f2py
# Ting Sun
# ting.sun@reading.ac.uk
# 10 Oct 2017

# set -fbracket-depth=1024 as clang may complain with a lower limit
export CFLAGS='-fno-strict-aliasing -fno-common -dynamic -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fbracket-depth=1024'
export FFLAGS='-g -pg -Wall -Wtabs -fbounds-check -cpp -Wno-unused-dummy-argument -Wno-unused-variable -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal'

# clean existing built modules
rm -rf SUEWS_driver.so*

# compile SUEWS with Makefile
cd ..
make clean
make

# enter the f2py stage
cd f2py || return
# link compiled objects
ln -sf ../*.o .
# ln -sf ../*.mod .

# remove duplicate SUEWS_driver to avoid symbol conflict
rm SUEWS_driver.o

# direct compilation without f2py header file
# --debug-capi provides more verbose information when calling C-API: useful for debugging
# f2py -m SUEWS_driver -c --debug-capi ../SUEWS_driver.f95 *.o
f2py -m SUEWS_driver -c ../SUEWS_driver.f95 *.o

# cleaning
rm *.o
