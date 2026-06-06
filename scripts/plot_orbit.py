import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import re


if len(sys.argv)!=2 and len(sys.argv)!=4:
    print(f"Usage: python {sys.argv[0]} <lightcurve> [<tmin> <tmax>]")
    exit()

lightcurve = sys.argv[1]

data = pd.read_csv(lightcurve,sep=r'\s+',comment='#')
outfile = lightcurve[:lightcurve.rfind('_')] +'.out'
idx = lightcurve[lightcurve.rfind('_')+1:lightcurve.find('.')]
#print(lightcurve,outfile,idx)
out = pd.read_csv(outfile,sep=r'\s+')
outdata = out[out['EventID']==int(idx)].squeeze()
#print(list(outdata.index))
#print(list(outdata))

for i in range(int(outdata['NLens'])-1):
    keys = list(filter(re.compile(f'Planet_{i}_.+').match,list(outdata.keys())))
    #if i==1:
    print(keys)
    print(i,outdata[keys].to_list())

for i in range(int(outdata['NPlanets'])):
    keys = list(filter(re.compile(f'^p_{i}_.+').match,list(outdata.keys())))
    #if i==1:
    print(keys)
    print(i,outdata[keys].to_list())

if 'LensCompanions' in outdata.keys() and outdata['LensCompanions']>0:
    keys = list(filter(re.compile('Lens.*ID').match,list(outdata.keys())))
    print(keys)
    print(outdata[keys].to_list())


keys = list(filter(re.compile(f'^Lens_.+').match,list(outdata.keys())))
print(keys)
print(outdata[keys].to_list())

if 'LensCompanions' in outdata.keys() and outdata['LensCompanions']>0:
    keys = list(filter(re.compile(f'^Lens2_.+').match,list(outdata.keys())))
    print(keys)
    print(outdata[keys].to_list())


print("rE = ",outdata['rE'])
    


#Add breaks into the lightcurve data between seasons
#tdiff = data.loc[:-1,'Simulation_time'] - data.loc[1:,'Simulation_time']

#print(outdata[['Lens2_combined_logP','Lens2_a','Lens2_P']])
#print(outdata[['Lens_Mass','Lens2_Mass']])

for i in range(5):
    for k in ['period','a','dL']:
        key = f'p_{i}_{k}'
        if key in outdata.index:
            print(key,outdata[key])
    key=f'p_{i}_dL'
    if key in outdata.index:
        print(f"360/{key}",360.0/outdata[key])

for k in ['Lens_Mass']:
    print(k,outdata[k])

nlens = pd.Series(list(data.columns)).str.contains(r'^lens\d+_x_thE$').sum()
#print(pd.Series(list(data.columns)).str.contains(r'^lens\d+_x_thE$'))
print(f"nlens = {nlens}")
nsrc = pd.Series(list(data.columns)).str.contains(r'^Source\d*_x_thE$').sum()
#print(pd.Series(list(data.columns)).str.contains(r'^source\d+_x_thE$'))
print(f"nsrc = {nsrc}")

#header = pd.read_csv(sys.argv[1],sep='\s+',header=None,comment=None,engine='python',nrows=50,index_col=False)
#fsm = header[header.iloc[:,0]=='#fs:'].squeeze(axis=0)[1:].astype(float)
#event = header[header.iloc[:,0]=='#Event:'].squeeze(axis=0)[1:].astype(float)
#planet = header[header.iloc[:,0]=='#Planet:'].squeeze(axis=0)[1:].astype(float)
#source = header[header.iloc[:,0]=='#Obssrcmag:'].squeeze(axis=0)[1:].astype(float)
#lens = header[header.iloc[:,0]=='#Obslensmag:'].squeeze(axis=0)[1:].astype(float)

ls = ['-','--','-.']

#fig,axtmp = plt.subplots(3,3,sharex=True,sharey=True,squeeze=True)
#fig,axtmp = plt.subplots(1,1,sharex=True,sharey=True,squeeze=True)
plt.figure()

#ax = axtmp.flatten()
ax = [plt.gca()]

print(ax)

ms=4
if len(sys.argv)==4:
    ms=1


for i in range(nlens):
    ax[0].plot(data[f"lens{i}_x_thE"],data[f"lens{i}_y_thE"],'o',ms=ms,label=f'L{i}')
for i in range(nsrc):
    srcstr=""
    if i>0:
        srcstr=i+1
    ax[0].plot(data[f"Source{srcstr}_x_thE"],data[f"Source{srcstr}_y_thE"],'o',ms=ms,label=f'S{i}')

if len(sys.argv)==4:
    tmin=float(sys.argv[2])
    tmax=float(sys.argv[3])
    tdata = data[(data['Simulation_time']>=tmin) & (data['Simulation_time']<tmax)]

    for i in range(nlens):
        ax[0].plot(tdata[f"lens{i}_x_thE"],tdata[f"lens{i}_y_thE"],'o',ms=ms*3,label=f'L{i}')
    for i in range(nsrc):
        srcstr=""
        if i>0:
            srcstr=i+1
        ax[0].plot(tdata[f"Source{srcstr}_x_thE"],tdata[f"Source{srcstr}_y_thE"],'o',ms=ms*3,label=f'S{i}')


ax[0].set_aspect('equal')
ax[0].legend()
#plt.colorbar(label='Time [days]')
ax[0].set_xlabel(r'$x$ [$r_{\rm E}$]')
ax[0].set_ylabel(r'$y$ [$r_{\rm E}$]')
ax[0].grid()

# for j in range(nlens):
#     for i in range(nlens):
#         ax[j+1].plot(data[f"lens{i}_x"]-data[f"lens{j}_x"],data[f"lens{i}_y"]-data[f"lens{j}_y"],'o',ms=ms,label=f'L{i}')
#         ax[j+1].set_aspect('equal')


# for i in range(nsrc):
#     for j in range(nlens):
#         ax[j+1].plot(data[f"source{i}_x"]-data[f"lens{j}_x"],data[f"source{i}_y"]-data[f"lens{j}_y"],'o',ms=ms,label=f'S{i}')

plt.tight_layout()

        
plt.show()
