import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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

filename = sys.argv[1]

data = pd.read_csv(filename,sep=r'\s+',comment='#')

header = pd.read_csv(filename,sep=r'\s+',header=None,comment=None,engine='python',nrows=50,index_col=False)
fsm = header[header.iloc[:,0]=='#fs:'].squeeze(axis=0)[1:].astype(float)
event = header[header.iloc[:,0]=='#Event:'].squeeze(axis=0)[1:].astype(float)
planet = header[header.iloc[:,0]=='#Planet:'].squeeze(axis=0)[1:].astype(float)
source = header[header.iloc[:,0]=='#Obssrcmag:'].squeeze(axis=0)[1:].astype(float)
#lens = header[header.iloc[:,0]=='#Obslensmag:'].squeeze(axis=0)[1:].astype(float)

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

maglabels = ['W146','Z087','K213','W146','Z087','K213','W146','Z087','K213']


#fig,ax = plt.subplots(ncols=1,nrows=2)
fig = plt.figure(layout="constrained",figsize=(20,12.5))
fig.suptitle(filename)
gs = GridSpec(3,1,figure=fig)
ax1 = fig.add_subplot(gs[0:2,0])
ax2 = fig.add_subplot(gs[2,0],sharex = ax1)

#Find the baseline magnitude and source flux ratio for the reference observatory
fs0 = fsm.iloc[match]
#print(type(fs0))
m0 = source.iloc[match] + 2.5*np.log10(fs0)
#print(m0,fs0,source.iloc[match])
#print(fsm)

#Plot the data
for ii,i in enumerate(displayobs):
    d = data[(data['observatory_code']==i) & (data['saturation_flag']==0)]
    #Find the source flux ratio and source magnitude
    fs=fsm.iloc[i]
    ms=source.iloc[i]

    #Simulation_time measured_relative_flux measured_relative_flux_error true_relative_flux true_relative_flux_error observatory_code saturation_flag best_single_lens_fit parallax_shift_t parallax_shift_u BJD source_x source_y source2_x source2_y lens1_x lens1_y lens2_x lens2_y parallax_shift_x parallax_shift_y parallax_shift_z

    #Calculate magnification
    mu = A(d['measured_relative_flux'],fs)
    mutrue = A(d['true_relative_flux'],fs)
    mufit = A(d['best_single_lens_fit'],fs)
    sigmu = d['measured_relative_flux_error']/fs

    #Calculate scaled magnitude
    mi = m0 - 2.5*np.log10(fs0*mu+1-fs0)
    mitrue = m0 - 2.5*np.log10(fs0*mutrue+1-fs0)
    mifit = m0 - 2.5*np.log10(fs0*mufit+1-fs0)
    sigmi = 2.5/np.log(10) * sigmu/(fs0*mu+1-fs0) * fs0


    mask = (sigmi>0)

    #Plot with errorbars
    if ii==0:
        ax1.plot(d.loc[mask,'Simulation_time'],mi[mask],'o',ms=4,
                 label=maglabels[i],color='C%d' % (ii))
        ax2.plot(d.loc[mask,'Simulation_time'],mi[mask]-mifit[mask],'o',ms=4,
                 label=maglabels[i],color='C%d' % (ii))
    else:
        ax1.errorbar(d.loc[mask,'Simulation_time'],mi[mask],yerr=sigmi[mask],fmt='o',ms=4,
                     label=maglabels[i],color='C%d' % (ii))
        ax2.errorbar(d.loc[mask,'Simulation_time'],mi[mask]-mifit[mask],yerr=sigmi[mask],fmt='o',ms=4,
                     label=maglabels[i],color='C%d' % (ii))

    ax1.autoscale_view()
    ax2.autoscale_view()
    ax1.autoscale(False)
    ax2.autoscale(False)
    
ax2.set_xlabel('Time (days)')
ax1.set_ylabel(f'{maglabels[match]} magnitude')
ax2.set_ylabel(f'Resid.')
ax1.legend()


#Plot the model lightcurves
#if dispobs>0 plot the display observatory/ies lighcurves, else plot the
#reference model lightcurve
#print(dispobs,match)
for i in [match]: #(dispobs,match)[dispobs==0]:
    d = data[(data['observatory_code']==i) & (data['saturation_flag']==0)]
    #Find the source flux ratio and source magnitude
    fs=fsm.iloc[i]
    print(fs)
    ms=source.iloc[i]
    print(ms)

    #Calculate magnification
    mu = A(d['measured_relative_flux'],fs)
    mutrue = A(d['true_relative_flux'],fs)
    mufit = A(d['best_single_lens_fit'],fs)
    sigmu = d['measured_relative_flux_error']/fs

    #Calculate scaled magnitudes
    mi = m0 - 2.5*np.log10(fs0*mu+1-fs0)
    mitrue = m0 - 2.5*np.log10(fs0*mutrue+1-fs0)
    mifit = m0 - 2.5*np.log10(fs0*mufit+1-fs0)
    sigmi = 2.5/np.log(10) * sigmu/(fs0*mu+1-fs0) * fs0

    #print(mitrue)

    #Plot with lines
    ax1.plot(d['Simulation_time'],mitrue,'-',color='k',alpha=0.5,zorder=10)
    ax1.plot(d['Simulation_time'],mifit,'--',color='r',alpha=0.5,zorder=11)

    ax2.plot(d['Simulation_time'],mitrue-mifit,'-',color='k',alpha=0.5,zorder=10)
    ax2.plot(d['Simulation_time'],mifit-mifit,'--',color='r',alpha=0.5,zorder=11)


#We're in magnitudes    
ax1.invert_yaxis()
ax2.invert_yaxis()


plt.tight_layout()
plt.show()
