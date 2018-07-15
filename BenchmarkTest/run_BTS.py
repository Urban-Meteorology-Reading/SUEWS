import Benchmark_SUEWS as bts
reload(bts);
fn_nml = 'BTS_config.nml'
bts.report_benchmark_PDF(fn_nml)

# test code


# benchmark performance


# generate report


#
# dir(bts)
#
#
# df_metric = bts.benchmark_SUEWS(fn_nml)
#
# df_metric.columns.get_level_values('stat').unique()
# df_metric['MAE'].unstack()
#
# df_metric
