%whos
# from SUEWS_driver import SUEWS_driver as sd
import SUEWS_driver
reload(SUEWS_driver)
from SUEWS_driver import suews_driver as sd

print sd.suews_cal_main.__doc__
pd.options.display.max_rows = 4000

os.getcwd()

dens_dry

varInputInfo
lib_Profiles
dict_SiteSelect[1]['AnthropogenicCode']
dict_SiteSelect[1]['Code_Paved']
pd.Series(lib_Veg.columns)

pd.Series(lib_AnthropogenicHeat.columns)
pd.options.display.max_rows = 400

lib_Veg.loc[300, ['BaseT', 'AlbedoMin']]

xx = {1: 3, 2: 4, 5: 6}

xx[[1, 2]]

{'Code_DecTr': 'BiogenCO2Code': 'beta'}

pd.Series(lib_SiteSelect.columns).sort_values()

pd.Series(lib_AnthropogenicHeat.columns)

pd.Series(lib_Veg.columns)

pd.Series(lib_Snow.columns)

pd.Series(lib_OHMCoefficients.columns)
list_file_MetForcing
list_file_MetForcing
pd.Series(lib_Soil.columns)
pd.Series(lib_Water.columns)

pd.Series(lib_NonVeg.columns)
pd.Series(lib_Water.columns)
pd.Series(lib_NonVeg.columns)
dict_SiteSelect

pd.Series(lib_AnthropogenicHeat.columns)


pd.Series(lib_Snow.columns)
list_varGridSurfaceChar
list_var_mod_state

list_var_mod_forcing

list_var_mod_config

pd.Series(lib_AnthropogenicHeat.columns)


pd.Series(lib_WithinGridWaterDist.columns)


xx = {2: [2, 3]}

dict_SiteSelect


len(dict_var2SiteSelect.keys())


set([type(v) for k, v in dict_var2SiteSelect.iteritems()])

gridiv = 1
# case: string
varName = 'alt'
varKey = dict_var2SiteSelect[varName]
varVal = lib_SiteSelect.loc[gridiv, varVal]
print varName, varVal
# case: dict
varName = 'albmax_dectr'
varKey = dict_var2SiteSelect[varName]
indexKey, subKey = varKey.items()[0]
code = lib_SiteSelect.loc[gridiv, indexKey]
res = lookup_code(indexKey, code)[subKey]
print varName, varVal


varKey.items()[0]

type('AlbedoMax')


# look up properties according to code
def lookup_code_sub(codeName, codeKey, codeValue):
    str_lib = dict_Code2File[codeName].replace(
        '.txt', '').replace('SUEWS', 'lib')
    str_code = '{:d}'.format(int(codeValue))
    cmd = '{lib}.loc[{code},{key:{c}^{n}}]'.format(
        lib=str_lib, code=str_code, key=codeKey, n=len(codeKey) + 2, c='\'')
    # print cmd
    res = eval(cmd)
    return res


lookup_code_sub('Code_Paved', 'AlbedoMax', 661)
lib_BiogenCO2
lookup_code_sub('BiogenCO2Code', 'beta_enh', 31)


lib_


lookup_code_sub('timezone', {'AnthropogenicCode': [
                'AHSlope_Heating_WD', 'AHSlope_Heating_WE']}, 1)
lookup_KeySeq('alb', ['Fr_Paved',
                      'Fr_Bldgs',
                      'Fr_EveTr',
                      'Fr_DecTr',
                      'Fr_Grass',
                      'Fr_Bsoil',
                      'Fr_Water'], 1)


def lookup_KeySeq(indexKey, subKey, indexCode):
    if type(subKey) is str:
        res = lookup_code_sub(indexKey, subKey, indexCode)
    elif type(subKey) is dict:
        indexKeyX, subKeyX = subKey.items()[0]
        indexCodeX = lookup_code_sub(indexKey, indexKeyX, indexCode)
        res = lookup_KeySeq(indexKeyX, subKeyX, indexCodeX)
    elif type(subKey) is list:
        res = []
        for subKeyX in subKey:
            indexCodeX = indexCode
            resX = lookup_KeySeq(indexKey, subKeyX, indexCodeX)
            res.append(resX)
    # final result
    return res
    pass


zip([1, 2, 3], (9, 8, 7))
pd.options.display.max_columns = 200
pd.options.display.max_rows = 200


lib_Veg
pd.Series(lib_Water.columns.sort_values())

pd.Series(lib_AnthropogenicHeat.columns.sort_values())

dict_Code2File


lib_Irrigation

pd.Series(varToDo)


set(varInputInfo[:, 0])
len(set(list_varToDo) - set(list_varMethod) -
    set(list_varTstep) - set(lib_RunControl.index))


pd.Series(
    list(set(list_varToDo)
         - set(list_varMethod)
         - set(list_varTstep)
         - set(lib_RunControl.index)
         - set(list_varTime)
         - set(list_varMetForcing)
         - set(list_varDailyState)
         - set(list_varSetting)
         - set(list_varModelState)
         - set(list_varModelOut)
         - set(list_varModelCtrl)
         - set(list_varGridSurfaceCharX)
         )).sort_values()


np.array(set(list_varToDo) -
         set(list_varMethod) - set(list_varTstep) - set(lib_RunControl.index))
tstep


list_varToDo


cbluse

snowuse = 0

set(lib_RunControl.index)

cmd

tstep

pd.options.display.max_columns=200
pd.Series(list_varToDo).sort_values()


df_gridSurfaceChar.columns.size

pd.DataFrame(varInputInfo)



[x for x in varInputInfo[:,0] if x.lower().startswith('dataout')]



for var in lib_RunControl.index:
    print var

lib_RunControl.loc[:,'runcontrol']

lib_RunControl




# get the information of input variables for SUEWS_driver
docLines = np.array(
    sd.suews_cal_main.__doc__.splitlines(),
    dtype=str)
posInput = np.where(
    np.logical_or(
        docLines == 'Parameters', docLines == 'Other Parameters'))
varInputLines = docLines[posInput[0][0] + 2:posInput[0][1] - 1]
varInputInfo = np.array([[xx.rstrip() for xx in x.split(':')]
                         for x in varInputLines])
len(varInputInfo)
list_varGridSurfaceChar


print '\n'.join(docLines)



varOutputInfo

os.getcwd()






os.path.join(
    'Input', '{}*{}*txt'.format(filecode, resolutionfilesin / 60))


get
#
os.path.join(
    'Input', '{}*{}*txt'.format(filecode, resolutionfilesin / 60))

















































































list_file_MetForcing

os.getcwd()





















#
