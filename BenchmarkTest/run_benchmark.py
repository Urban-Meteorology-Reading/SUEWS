import os
import itertools
os.chdir('BenchmarkTest')
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
dict_method_options = {
    'cbluse': [0],  # CBL disabled
    'solweiguse': [0],  # SOLWEIG disabled
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


list_fail = [
    list_to_benchmark[k]
    for k, v in dict_test.items() if type(v) == str]

list_pass = [
    list_to_benchmark[k]
    for k, v in dict_test.items() if not type(v) == str]
# df_fail.columns
df_fail = pd.DataFrame(list_fail)
for c in df_fail.columns:
    print c, np.unique(df_fail.loc[:, c])


df_pass = pd.DataFrame(list_pass)
for c in df_pass.columns:
    print c, np.unique(df_pass.loc[:, c])


# combine results
df_test = pd.concat(
    (df_pass.assign(result='pass'),
     df_fail.assign(result='fail'))).reset_index(
    drop=True)

# locate invalid options:
for c in df_test.columns[:-1]:
    set_fail=set(df_test.loc[df_test.result=='fail',c])
    set_pass=set(df_test.loc[df_test.result=='pass',c])
    print c,list(set_fail-set_pass)

# add some styling
def color_rows(s):
    df = s.copy()

    # Key:Value dictionary of Column Name:Color
    color_map = {}

    # Unqiue Column values
    # manufacturers = df['Manufacturer'].unique()
    colors_to_use = ['background-color: #ABB2B9',
                     'background-color: #EDBB99',
                     'background-color: #ABEBC6',
                     'background-color: #AED6F1']

    # Loop over our column values and associate one color to each
    # for manufacturer in manufacturers:
    #     color_map[manufacturer] = colors_to_use[0]
    #     colors_to_use.pop(0)

    for index, row in df.iterrows():
        if row['result'] == 'fail':
            # manufacturer = row['Manufacturer']
            # Get the color to use based on this rows Manufacturers value
            # my_color = colors_to_use[1]
            # Update the row using loc
            df.loc[index, :] = 'red'
        else:
            df.loc[index, :] = 'green'
    return df
    # .to_html(
    # '~/Downloads/test-options.html')


path_html = '~/Downloads/test-options.html'
# df_test.style.apply(color_rows, axis=None).to_html(
# path_html)

styled = df_test.style.apply(
    color_rows, axis=None).set_table_styles(
    [{'selector': '.row_heading',
      'props': [('display', 'none')]},
     {'selector': '.blank.level0',
      'props': [('display', 'none')]},
     {'selector': "caption",
      'props':[("caption-side", "top"),
      ('font-size', '20px'),
      ('font-family', 'Helvetica'),
      ('font-weight', 'bold')
      ]}]).set_caption(
    'SUEWS option availability test')
html = (styled.render())
with open(os.path.expanduser(path_html), 'w') as fn:
    fn.write(html)
    fn.close()

#
