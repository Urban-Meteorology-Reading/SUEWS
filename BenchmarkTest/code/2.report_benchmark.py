import Benchmark_SUEWS as bs

# load global path
fn_nml = 'BTS_config.nml'


# benchmark performance and generate report
# PDF version
bs.report_benchmark_PDF(fn_nml)
# HTML version
bs.report_benchmark_HTML(fn_nml)
