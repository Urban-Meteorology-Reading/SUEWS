import numpy as np 
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt 


def convtime(year,doy,hour,min):
    d = datetime.strptime(str(doy),'%j')
    month = d.month
    day = d.day
    dt = datetime(int(year),int(month),int(day),int(hour),int(min))
    return dt

df = pd.read_csv('./Output/test1_2005_RSL_60.txt',delim_whitespace=True)
datex = []
for ind in df.index:
    datex.append(convtime(df.at[ind,'Year'],df.at[ind,'DOY'],df.at[ind,'Hour'],df.at[ind,'Min']))
df['datetime'] = datex
df.set_index('datetime',inplace=True)

print(df)

plt.figure(figsize=(12,5))
plt.pcolormesh(df.index,np.arange(0.,3,0.1),df.iloc[:,5:35].T,cmap='YlGnBu',vmin=0,vmax=5)
plt.xlabel('Month-day hour')
plt.ylabel('$z/z_H$ [-]')
plt.colorbar(label='U [m s$^{-1}$]')
plt.show()


