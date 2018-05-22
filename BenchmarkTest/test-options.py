import os
import itertools
# os.chdir('/Users/sunt05/Dropbox/8-Research/98.ReadingWork/10.SUEWS-FORTRAN/BenchmarkTest')
os.getcwd()
from Test_SUEWS import *
import pandas as pd
# load basic configurations
# initial condition
dict_initcond = (load_SUEWS_nml(
    os.path.join(dir_input, 'InitialConditionstest_2004.nml')).to_dict()[
    'initialconditions'])
# runcontrol settings
dict_runcontrol = load_SUEWS_nml(
    os.path.join(dir_input, 'RunControl.nml')).to_dict()['runcontrol']
# siteselect info
df_siteselect = load_SUEWS_table(
    os.path.join(dir_input, 'SUEWS_SiteSelect.txt'))

# available configuration list
list_method = [
    x for x in dict_runcontrol.keys()
    if ('method' in x)
    and ('disagg' not in x)
    or ('use' in x)]


for x in list_method:
    print x

# all available options
dict_method_options_all = {
    'cbluse': [0],  # CBL disabled for 2018a
    'solweiguse': [0],  # SOLWEIG disabled for 2018a
    'snowuse': [0, 1],
    'stabilitymethod': [2, 3, 4],
    'netradiationmethod': [0, 1, 2, 3, 100, 200, 300],
    'smdmethod': [0],
    'storageheatmethod': [1, 3, 4],
    'emissionsmethod': [0, 1, 2],
    'waterusemethod': [0],
    'roughlenheatmethod': [1, 2, 3, 4],
    'roughlenmommethod': [1, 2, 3]}

# selected options for testing
dict_method_options_sel = {
    'snowuse': [0, 1],
    'stabilitymethod': [2, 3, 4],
    'netradiationmethod': [3],
    'storageheatmethod': [1, 3, 4],
    'emissionsmethod': [0, 1, 2]}


methods, options = zip(*dict_method_options_sel.items())
list_to_benchmark = [dict(zip(methods, v))
                     for v in itertools.product(*options)]
len(list_to_benchmark)


# run testing
dir_test = '~/Downloads/20180507test' + str(np.random.randint(10000))
dict_test = {}
for ind, cfg in enumerate(list_to_benchmark):
    runcontrol_test = dict_runcontrol.copy()
    runcontrol_test.update(cfg)
    name_sim = str(ind)
    res_sim = run_sim(name_sim, runcontrol_test, dict_initcond, df_siteselect,
                      dir_save=dir_test)
    dict_test.update({ind: res_sim})

dict_test_OK = {k: 'fail' if type(v) == str else 'pass'
                for k, v in dict_test.iteritems()}

df_test = pd.DataFrame(list_to_benchmark).assign(result=dict_test_OK.values())


# add some styling
def color_rows(s):
    df = s.copy()
    # Unqiue Column values
    # manufacturers = df['Manufacturer'].unique()
    colors_to_use = ['background-color: #ABB2B9',
                     'background-color: #EDBB99',
                     'background-color: #ABEBC6',
                     'background-color: #AED6F1']

    for index, row in df.iterrows():
        if row['result'] == 'fail':
            # Update the row using loc
            df.loc[index, :] = colors_to_use[1]
        else:
            df.loc[index, :] = colors_to_use[2]
    return df


styles = [
    # table properties
    dict(selector=" ",
         props=[("margin", "0"),
                ("font-family", '"Helvetica", "Arial", sans-serif'),
                ("border-collapse", "collapse"),
                ("border", "none"),
                ("border", "2px solid #ccf")
                ]),

    # header color - optional
    dict(selector="thead",
         props=[("background-color", "#ABB2B9")
                ]),

    # background shading
    dict(selector="tbody tr:nth-child(even)",
         props=[("background-color", "#fff")]),
    dict(selector="tbody tr:nth-child(odd)",
         props=[("background-color", "#eee")]),

    # cell spacing
    dict(selector="td",
         props=[("padding", ".5em")]),

    # header cell properties
    dict(selector="th",
         props=[("font-size", "125%"),
                ("padding", ".5em"),
                ("text-align", "center")]),

    # caption placement
    dict(selector="caption",
         props=[("caption-side", "top")]),

    # render hover last to override background-color
    dict(selector="tbody tr:hover",
         props=[("background-color", "%s" % "#add8e6")])]


df_test_styled = df_test.style.apply(
    color_rows, axis=None).set_table_styles(styles).set_caption(
    'SUEWS option test results')
df_test_html = (df_test_styled.render())

# %% test results
list_method_test = [c for c in df_test.columns if not c == 'result']
df_test_pass = pd.concat(
    [df_test.loc[:, [c, 'result']].pivot_table(
        index='result', columns=c, aggfunc=len)
        for c in list_method_test],
    keys=list_method_test,
    axis=1).loc[
    'pass', :].to_frame().rename(
    columns={'pass': 'result'}).applymap(
    lambda x: 'pass' if x > 0 else 'fail')
df_test_pass.index.set_names(['method', 'option'], inplace=True)
df_test_pass


styles_save = [
    # {'selector': '.row_heading',
    #  'props': [('display', 'none')]},
    # {'selector': '.blank.level0',
    #  'props': [('display', 'none')]},
    {'selector': '',
     'props': [('width', '100%'),
               ('border', '1px solid black'),
               ('border-collapse', 'collapse')]},
    {'selector': 'thead th',
     'props': [('font-family', 'Helvetica'),
               ('font-size', '18px'),
               ('height', '50px'),
               ('background-color', '')]},
    {'selector': 'tr',
     'props': [('font-family', 'Helvetica'),
               ('background-color', '')]},
    {'selector': 'table, th, td',
     'props': [('font-family', 'Helvetica'),
               ('text-align', 'center'),
               ('border', '1px solid black'),
               ('width', '100%'),
               ('padding', '15px')]},
    # {'selector': 'tbody th:first-child',
    #  'props': [('display', 'none')]},
    {'selector': "caption",
        'props': [("caption-side", "top"),
                  ('font-size', '20px'),
                  ('font-family', 'Helvetica'),
                  ('font-weight', 'bold')
                  ]}]

df_test_pass_styled = df_test_pass.style.apply(
    color_rows, axis=None).set_table_styles(
    styles_save).set_caption(
    'SUEWS option availability')
df_test_pass_html = (df_test_pass_styled.render())

# print(df_test_pass_html)

# %%output html
res_html = df_test_pass_html
# res_html = '\n'.join([df_test_pass_html, df_test_html])
# print res_html
# path_html = '~/Downloads/test-options.html'
path_html = 'test-options.html'
with open(os.path.expanduser(path_html), 'w') as fn:
    fn.write(res_html)
    fn.close()
os.getcwd()
#
