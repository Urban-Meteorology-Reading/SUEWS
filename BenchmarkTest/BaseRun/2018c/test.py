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
plt.draw()

plt.figure()
colormap = plt.cm.YlGnBu
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, 30)])

for x in np.arange(len(df.index)/4):
    plt.plot(df.iloc[x,5:35].T,np.arange(0.,3,0.1))

plt.xlabel('U [m s$^{-1}$]')
plt.ylabel('$z/z_H$ [-]')
plt.show()


