#!/usr/bin/env python
from Benchmark_SUEWS import *
import os
dir_path=os.getcwd()
plt.close('all')

# do benchmarking and you'll get some warning
res_BSS=benchmarkSUEWS('benchmark.nml',dir_path)
# get figure numbers
plt.get_fignums()
# show the figures
plt.show()

# this will produce a PDF
report_benchmark('benchmark.nml',dir_path)
