dict_output.keys()

dict_output[1]
df_state.columns
df_state.index




df_state.keys()
sp.conv2PyData(dict_mod_cfg)







dict_mod_cfg
xx = sp.df(dict_mod_cfg).to_dict('split')
xx['data'] = [np.array(var).tolist()
              for var in xx['data']]
df_xx = sp.df(**xx)

df_xx.keys()
for k,v in df_xx.items():
    print k, v, type(v)
