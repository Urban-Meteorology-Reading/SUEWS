#!/usr/bin/env python
from Benchmark_SUEWS import *
import sys
dir_path = os.getcwd()

report_benchmark('benchmark.nml',dir_path)
# res_all = load_res_all('benchmark.nml', dir_path)
#
# res_all['base'].index.size
#
# var_list_all=res_all['base'].columns
# var_list_user=['Kup','QE','xx']
# var_list_all.intersection(var_list_user)

# plotMatMetricX(res_all, func_MAE, 'MAE')
# res_all.keys()
# base = res_all['base']
# resPlot = pd.DataFrame([func_MAE(x, base)
#                         for k, x in res_all.iteritems()], index=res_all.keys())
#
# # resPlot[0]
# base
# plt.clf()
# base.dropna(axis=1, how='all').dropna(axis=0, how='any').plot()
# plt.show()
#
#
# # data cleaning
# res_all_clean = pd.Panel({k: df.dropna(axis=1, how='all').dropna(
#     axis=0, how='any') for k, df in res_all.to_frame().to_panel().iteritems()})
# plotMatMetricX(res_all_clean, func_MAE, 'MAE')
# res_all_clean['ref']
# df_all = res_all_clean.to_frame()
#
# df_all.columns
#
# pnl_all=df_all.to_panel()
# pnl_all.keys()
# pnl_all['base'].ix[:,'QN']
#
# df_all_clean=df_all.copy()
# pnl_all_clean=df_all_clean.to_panel()
# pnl_all_clean[0]
#
# np.arange(3)+1
# benchmarkSUEWS('benchmark.nml',dir_path);
# plt.show();
#
# report_benchmark('benchmark.nml',dir_path)
#
# reload(sys.module['Benchmark_SUEWS'])
#
# np.abs([-3,2])
