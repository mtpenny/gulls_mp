import pandas as pd
import numpy as np
import sys
import argparse
from glob import glob
import os
import re

parser = argparse.ArgumentParser(description="Process and reweight a gulls simulation")
parser.add_argument('paramfile',
                    help='The parameter file for the simulation')
parser.add_argument('--covfac',default=None,
                    help='Covering factor file')
parser.add_argument('--ffp-bad-weights',action='store_true',
                    help='Recalculate the weights for free-floating planet simulations. If an argument is given, it is the number of days over which t0 could have been drawn.')
parser.add_argument('--t0range', action='store', default=432.0, type=float,
                    help='The range of days over which t0 could be drawn (minus any gaps). Only used if a bad weights option is applied.')
parser.add_argument('--in-raw',action='store_true',
                    help='The output files are in the raw subdirectory of FINAL_DIR')
parser.add_argument('--move-files',action='store_true',
                    help='Move the output files to their homes')
parser.add_argument('--limit',action='store',type=int,
                    help='Limit the number of files to read (option for testing)')
parser.add_argument('--ffp-det-cuts',action='store',type=float,nargs='*',
                    help="Apply FFP detection cuts, the two arguments are the Delta chi^2 threshold and the N_>3sigma trhreshold")
parser.add_argument('--split-masses',action='store_true',
                    help='Run contains many different masses, split them by unique masses')
parser.add_argument('--recut',action='store_true',
                    help='Do not reprocess from scratch. Load the old "out" hdf5 file and apply the cuts again.')
parser.add_argument('-k','--add-key',action='append',
                    help='Add keys to be added to the hdf5store index over the default (use --print-keys to print the current set without running further)')
parser.add_argument('--print-keys',action='store_true',
                    help='Print the current set of keys to be indexed in the hdf5store without running further')
parser.add_argument('--save-csv',action='store_true',
                    help='Save the HDF5 files as CSV files too')

args = parser.parse_args()

#Read the parameter file
paramfile = dict(pd.read_csv(args.paramfile,sep='=',header=None,names=['param','value']).values)
if 'OBS_GROUP_NAMES' in paramfile:
    paramfile['OBS_GROUP_NAMES'] = paramfile['OBS_GROUP_NAMES'].split(',')
print("Parameter file:")
print(paramfile)

#Index keys for the HDF5store
index_keys = ["Field","Planet_mass","t0lens1","u0lens1","tcroin","ucroin",'LCOutput','Lens_Dist','Source_Dist','Lens_pop','Source_pop']

if args.add_key is not None:
    index_keys.extend(args.add_key)


#Set up cuts
if args.ffp_det_cuts is not None:
    ffp_det_cuts = args.ffp_det_cuts
    print(ffp_det_cuts)
    print(len(ffp_det_cuts))
    if len(ffp_det_cuts)<1:
        ffp_det_cuts.append(300.0)
    if len(ffp_det_cuts)<2:
        ffp_det_cuts.append(6)
    print(ffp_det_cuts)
    print(len(ffp_det_cuts))

    for i,obsgroup in enumerate(paramfile['OBS_GROUP_NAMES']):
        index_keys.extend([f"ObsGroup_{i}_chi2",f"ObsGroup_{i}_flatchi2"])
    print("index_keys",index_keys)



index_keys = list(dict.fromkeys(index_keys))

if args.print_keys:
    print("index_keys:")
    print(index_keys)
    exit()    


print("Processing ",paramfile['RUN_NAME'])

outroot = paramfile['FINAL_DIR'] + paramfile['RUN_NAME'] + '/analysis/' + paramfile['RUN_NAME']

#Compute the normalization for each field. For this we only need the raw weights
os.chdir(paramfile['FINAL_DIR'] + paramfile['RUN_NAME'] + '/')
if not os.path.exists('analysis'):
    os.mkdir('analysis')
if not os.path.exists('raw'):
    os.mkdir('raw')
if not os.path.exists('lc'):
    os.mkdir('lc')
if not os.path.exists('fm'):
    os.mkdir('fm')
if not os.path.exists('logs'):
    os.mkdir('logs')

if args.in_raw:
    os.chdir(paramfile['FINAL_DIR'] + paramfile['RUN_NAME'] + '/raw/')
else:
    outfiledir = paramfile['OUTPUT_DIR'] + paramfile['RUN_NAME'] + '/'
    if not os.path.isdir(outfiledir):
        outfiledir = paramfile['FINAL_DIR'] + paramfile['RUN_NAME'] + '/'
    os.chdir(outfiledir)

#Files for processing
procfiles = {}





#Read the rates file
if not args.recut:
    rates = pd.read_csv(os.getenv('GULLS_BASE_DIR') + paramfile['RATES_FILE'],sep='\s+',usecols=['ID_src','l_src','b_src','nevents_per_tile']).set_index('ID_src')
    rates['norm']=0.0
else:
    rates = pd.read_csv(outroot + '.det.rates').set_index('ID_src')


if args.covfac is None:
    rates['covfac']=1.0
else:
    #Hasn't been tested yet
    cftmp = pd.read_csv(args.covfac,sep='\s+',names=['ID','covfac']).set_index('ID')
    rates.join(cftmp)
print(rates)

lcmatch = re.compile(r'.+\.lc$')
fmmatch = re.compile(r'.+\.fm\.\d+$')
outmatch = re.compile(r'\.out$')
logmatch = re.compile(r'\.log$')


if not args.recut:
    #We're starting from scratch

    #Conduct a first pass to compute the normalizations for each field
    
    nsim1 = 0 #Number of events processed on the first pass
    nsim2 = 0 #Number of events processed on the second pass
    nproc1 = 0 #number of ouput files processed
    nproc2 = 0 #number of ouput files processed
    for f in os.scandir():
        if outmatch.search(f.name):
            #print(f.name)
            try:
                data = pd.read_csv(f.name,usecols=['Field','SubRun','raw_weight'],
                                   dtype={'Field':int,'raw_weight':float},
                                   sep='\s+') #.to_numpy()
            except:
                print(f"Skipping {f.name} with no data")
                continue
        
            print(f"Pass 1 {f.name}")
        
            field = data.loc[0,'Field']
            subrun = data.loc[0,'SubRun']
            norm = data['raw_weight'].sum()
            rates.loc[field,'norm'] += norm
            procfiles[f.name] = (field,subrun,data.shape[0],norm) #We'll track the inputs with this dictionary and check if it changed, then correct the normalization if the data changed
        
            nsim1 += data.shape[0] #We'll want to make sure the file is not still being written, this allows a very basic test
            nproc1 += 1
            del data

        #If testing, only process a few files
        if args.limit is not None and nproc1>=args.limit:
            break

    #print(rates)
    #print(rates[rates['norm']>0])

    #Accumulate aggregate rates (in events/survey)
    rates['final_rate']=0.0
    rates['final_rate_err2']=0.0
    rates['norm_correction']=0.0
    
    detrates = rates[['l_src','b_src']]
    
    #Take a second pass through the data to normalize it and produce the combined output files
    outstore = pd.HDFStore(outroot + '.out.hdf5','w')
    detstore = pd.HDFStore(outroot + '.det.hdf5','w')


    for f in procfiles:
        if outmatch.search(f):

            #Read the file a second time
            try:
                data = pd.read_csv(f,sep='\s+')
            except:
                print("Skipping {f} with no data")
                continue

            print(f"Pass 2 {f}")

            field = data.loc[0,'Field']
            subrun = data.loc[0,'SubRun']
            print(field)

            #If the data file changed between first and second pass
            if data.shape[0] != procfiles[f][2]:
                rates['norm_correction'] += procfiles[f][3]
                #Skip including the data, and we'll need to fix the normalization of all of the
                #other fields affected
                del data
                continue

            #Handle any weight modifications due to prior errors
            if args.ffp_bad_weights:
                data['weight'] = data['raw_weight'] * np.sqrt(data['Planet_q']) * data['u0max'] * args.t0range/365.25 / (data['tE_ref']/data['tE_helio'])

            #Normalize and scale to event rates to produce weights in units of detections/survey   
            wmult = rates.loc[field,'nevents_per_tile']/rates.loc[field,'norm']
            nsim2 += data.shape[0]
            nproc2 += 1
            data['final_weight'] = wmult * data['weight']

            #Accumulate the rates and their uncertainties
            rates.loc[field,'final_rate'] += data['final_weight'].sum()
            rates.loc[field,'final_rate_err2'] += (data['final_weight']**2).sum()

            #Rename some columns to conform with HDF5 format rules
            rename = {}
            for col in data.columns:
                if col.find('/') != -1 or col.find('[') != -1 or col.find(']') != -1 or col.find('.') != -1:
                    rename[col] = col.replace('/','_').replace('[','').replace(']','').replace('.','')
            if len(rename)>0:
                data.rename(columns=rename,inplace=True)

            #Append the data to the HDF5 store
            print(f"Appending {f} to HDF5 store")
            outstore.append("out",data,index=False,data_columns=index_keys)
            print("done")

            #Free up the memory from the input data
            del data

            #For testing, only process a limited number of files
            if args.limit is not None and nproc2>=args.limit:
                break

    if rates['norm_correction'].sum()>0:
        print("Warning: Some files changed between passes and the rates and normalization for these fields are corrupted. Fields affected:")
        print(rates[rates['norm_correction']>0].index)
        print("A future version of this script will remove affected data")


    #Convert variances to errors for the aggregate rates
    rates['final_rate_err2'] = np.sqrt(rates['final_rate_err2'])
    rename = {'final_rate_err2': 'final_rate_err'}
    rates.rename(columns=rename,inplace=True)

            
else:
    #We've run the normalization before, we just want to reapply the detection cuts

    print(f"--recut option used, starting from the stored HDF5 file {outroot}.out.hdf5, overwriting {outroot}.det.hdf5")

    outstore = pd.HDFStore(outroot + '.out.hdf5','r')
    detstore = pd.HDFStore(outroot + '.det.hdf5','w')



#Apply detection cuts to the HDF5 store

detrates = rates.copy()

print("\n".join(outstore.select('out', start=1, stop=1).columns.tolist()))

for i,obsgroup in enumerate(paramfile['OBS_GROUP_NAMES']):

    frkey=f'ObsGroup_{obsgroup}_final_rate'
    frekey=f'ObsGroup_{obsgroup}_final_rate_err2'
    
    detrates[frkey]=0.0
    detrates[frekey]=0.0

    #We have a bunch of observatory groups and want an answer for each group
    if args.ffp_det_cuts is not None:
        ogchi2key = f'ObsGroup_{i}_chi2'
        ogn3sigkey = f'ObsGroup_{i}_flatchi2'
        wherestring = f'({ogchi2key} >= {ffp_det_cuts[0]} & {ogn3sigkey} >= {ffp_det_cuts[1]})'
        #print(outstore.select("out")[[ogchi2key,ogn3sigkey]])
        print("wherestring:",wherestring)
        for detdf in outstore.select("out",where=wherestring,chunksize=100000):
            #print("Selected ",detdf.shape)
            #print(detdf[[ogchi2key,ogn3sigkey]])
            print("Appending to detstore")
            detstore.append(obsgroup,detdf,index=False,data_columns=index_keys)
            print("Unique fields: ",detdf['Field'].unique())
            for ff in detdf['Field'].unique():
                detdf_field = detdf[detdf['Field']==ff]
                print(ff,detdf_field.shape)
                detrates.loc[ff,frkey] += detdf_field['final_weight'].sum()
                detrates.loc[ff,frekey] += (detdf_field['final_weight']**2).sum()
                
        
    detrates[frekey] = np.sqrt(detrates[frekey])
    detrename = {frekey: frekey[:-1]}
    detrates.rename(columns=detrename,inplace=True)

detrates.to_csv(outroot + '.det.rates')
                
        


        
print(rates[rates['norm']>0])
print(detrates.columns.tolist())
print("detrates (>0):")
print(detrates[detrates[f"ObsGroup_{paramfile['OBS_GROUP_NAMES'][0]}_final_rate"]>0])
print("detrates:")
print(detrates)

if args.save_csv:
    print("Saving the big output csv...")
    outdf = outstore.select('out')
    outdf.to_csv(outroot + '.out.csv',sep=' ')
    print("Done")

