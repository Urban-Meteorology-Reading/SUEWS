#!/usr/bin/env bash
# script for compiling pythonic SUEWS using f2py
# Ting Sun
# ting.sun@reading.ac.uk
# 10 Oct 2017

# set -fbracket-depth=1024 as clang may complain with a lower limit
export CFLAGS='-fno-strict-aliasing -fno-common -dynamic -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fbracket-depth=1024'
export FFLAGS='-g -pg -Wall -Wtabs -fbounds-check -cpp -Wno-unused-dummy-argument -Wno-unused-variable -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal'
f2py -c SUEWS.pyf ../*f95
