
import pandas as pd

xxdf=pd.read_csv('/Users/sunt05/Downloads/df_test.csv').iloc[:,1:]


xxdf

df_test=xxdf

list_method_test = [c for c in df_test.columns if not c == 'result']


df_test_pass = pd.concat(
    [df_test.loc[:, [c, 'result']].pivot_table(
        index='result', columns=c, aggfunc=len)
        for c in list_method_test],
    keys=list_method_test,
    axis=1).loc[
    'pass', :].to_frame().rename(
    columns={'pass': 'result'})


    # .applymap(
    # lambda x: 'pass' if x > 0 else 'fail')


zdf_test_pass.index.set_names(['method', 'option'], inplace=True)
