import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

import warnings
warnings.simplefilter(action='ignore')

if len(sys.argv) <3:
    print("Usage: %s <filename> <reference_observatory_number> {<number_of_obs_to_plot>}")
    sys.exit()

def A(F,fs):
    return (F-(1-fs))/fs

def mag(F,ms,fs):
    m0 = ms+2.5*np.log10(fs)
    return m0 - 2.5*np.log10(F)

def magerr(F,e,ms,fs):
    return 2.5/np.log(10)*e/F

plt.rc('font', size=36, family='Times')
plt.rc('axes', linewidth=6)

data = pd.read_csv(sys.argv[1],sep=r'\s+',header=None,comment='#')

fsm = pd.read_csv(sys.argv[1],sep=r'\s+',header=None,comment=None,engine='python',nrows=4,index_col=False)


event = pd.read_csv(sys.argv[1],sep=r'\s+',header=None,comment=None,engine='python',nrows=1,index_col=False,skiprows=8)
planet = pd.read_csv(sys.argv[1],sep=r'\s+',header=None,comment=None,engine='python',nrows=1,index_col=False,skiprows=7)
source = pd.read_csv(sys.argv[1],sep=r'\s+',header=None,comment=None,engine='python',nrows=1,index_col=False,skiprows=2)
lens = pd.read_csv(sys.argv[1],sep=r'\s+',header=None,comment=None,engine='python',nrows=1,index_col=False,skiprows=5)

#This is the number of the reference observatory
match=int(sys.argv[2])

displayobs=[]
ndisplayobs=0

#This is the number of observatories that are used for model lightcurves at the end
ndispobs=0
if len(sys.argv)>3:
    displayobs = [int(x) for x in sys.argv[3:]]
    ndisplayobs = len(displayobs)
#    for i,obsno in enumerate(range(3,len(sys.argv))):
#        displayobs = 
#        ndisplayobs += 1
else:
    nobs=data.iloc[-1,5]
    displayobs=list(range(nobs))
    ndisplayobs=nobs

#dispobs=list(range(nobs,nobs+ndispobs))


plt.figure(figsize=(15,10))

plt.gca().tick_params(axis="y",direction="in",length=10,width=6,top=True)
plt.gca().tick_params(axis="y",direction="in",length=10,width=6,bottom=True)
plt.gca().tick_params(axis="x",direction="in",length=10,width=6,left=True)
plt.gca().tick_params(axis="x",direction="in",length=10,width=6,right=True)

#Find the baseline magnitude and source flux ratio for the reference observatory
fs0 = fsm.iloc[0,match+1]
m0 = fsm.iloc[-1,match+1] + 2.5*np.log10(fs0)
print(m0,fs0,fsm.iloc[-1,match+1])
print(fsm)

#Plot the data
for ii,i in enumerate(displayobs):
    d = data[data.iloc[:,5]==i]
    #Find the source flux ratio and source magnitude
    fs=fsm.iloc[0,i+1]
    ms=fsm.iloc[-1,i+1]

    #Calculate magnification
    mu = A(d.iloc[:,1],fs)
    mutrue = A(d.iloc[:,3],fs)
    sigmu = d.iloc[:,2]/fs

    #Calculate scaled magnitude
    mi = m0 - 2.5*np.log10(fs0*mu+1-fs0)
    mitrue = m0 - 2.5*np.log10(fs0*mutrue+1-fs0)
    sigmi = 2.5/np.log(10) * sigmu/(fs0*mu+1-fs0) * fs0

    mask = (sigmi>0)

    #Plot with errorbars
    plt.errorbar(d.loc[mask,0],mi[mask],yerr=sigmi[mask],fmt='o',markersize=8,color='C%d' % (ii),lw=4)
    
plt.xlabel('Time [days]')
plt.ylabel(r'$z$ magnitude')

#Plot the model lightcurves
#if dispobs>0 plot the display observatory/ies lighcurves, else plot the
#reference model lightcurve
#print(dispobs,match)
for i in [match]: #(dispobs,match)[dispobs==0]:
    d = data[data.iloc[:,5]==i]
    #Find the source flux ratio and source magnitude
    fs=fsm.iloc[0,i+1]
    ms=fsm.iloc[-1,i+1]

    #Calculate magnification
    mu = A(d.iloc[:,1],fs)
    mutrue = A(d.iloc[:,3],fs)
    mufit = A(d.iloc[:,7],fs)
    sigmu = d.iloc[:,2]/fs

    #Calculate scaled magnitudes
    mi = m0 - 2.5*np.log10(fs0*mu+1-fs0)
    mitrue = m0 - 2.5*np.log10(fs0*mutrue+1-fs0)
    mifit = m0 - 2.5*np.log10(fs0*mufit+1-fs0)
    sigmi = 2.5/np.log(10) * sigmu/(fs0*mu+1-fs0) * fs0

    #print(mitrue)

    #Plot with lines
    plt.plot(d.loc[:,0],mitrue,'-',color='k',alpha=0.5,zorder=10,lw=3)
    plt.plot(d.iloc[:,0],mifit,'--',color='C%d' % (i),alpha=0.5,zorder=11)


#We're in magnitudes    
plt.gca().invert_yaxis()


plt.tight_layout()
plt.show()
