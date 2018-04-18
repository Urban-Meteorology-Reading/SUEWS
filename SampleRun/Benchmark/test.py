#!/usr/bin/env python
from Benchmark_SUEWS import *
import sys
from subprocess import call
# current directory
dir_path = os.getcwd()
# change to the code folder
os.chdir('../..')
print(os.getcwd())
# clean compiled files
os.system('make clean')
# make a new compiled version
os.system('make')
# change to the SampleRun folder
os.chdir('SampleRun')
# run SUEWS_V2017b
print(os.getcwd())
# call('SUEWS_V2017b')
os.system('./SUEWS_V2017c')

# return to the Benchmark folder and do benchmarking
os.chdir('Benchmark')
print(os.getcwd())
report_benchmark('benchmark.nml',dir_path)

os.system('open interface-test.pdf')
