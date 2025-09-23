_Add rest of command examples later_

---

## Installing GULLS on Unity @ OSU

_Note: these instructions were written with my setup in mind. I like keeping my program files in /home/crisp.92/Programs and the compiled results in /home/crisp.92/local. Your paths will need to be changed accordingly._

Operating System: Red Hat Enterprise Linux 8.10 (Ootpa)

Kernel: Linux 4.18.0-553.5.1.el8\_10.x86\_64

Architecture: x86-64

## Dependencies

### VBMicrolensing

1.  Go to [![](https://github.com/fluidicon.png)GitHub - valboz/VBMicrolensing: Microlensing computation code, including single, binary and multiple lenses](https://github.com/valboz/VBMicrolensing)
    
2.  Either clone the repository or download the contents as a zip file and transfer to Unity.
    
    1.  Note: cloning it will make updates easier, but we aren’t likely to need to update it anyway.
        
3.  In the VBMicrolensing directory, open the Makefile and double-check that the C compiler is g++.
    
4.  `make`
    
    1.  Note: we only need the C++ files, so there’s no need to pip install it unless you plan to use the Python libraries as well
        
5.  In your ~/.bashrc, add the line export `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/VBML/lib/`
    

### GNU Science Library

`module load gsl` OR…

_**Installing from Source**_

1.  Go to [https://ftp.gnu.org/gnu/gsl/](https://ftp.gnu.org/gnu/gsl/)
    
2.  Find the version you want and copy the link.
    
3.  On unity, use wget to download
    

`wget https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz`

4.  Untar + unzip
    
5.  Go into the expanded directory and do ./configure, followed by make, followed by make install
    
    1.  Note: if you want it installed in a particular directory, give ./configure a prefix flag.
        

./configure -- prefix=/home/crisp.92/local

### CFITSIO

_**Installing in Conda Environment**_

If using mamba rather than conda, just replace conda with mamba for the following.

1.  `module load conda`
    
2.  Create a conda environment for gulls, e.g.:
    

`conda create -n gulls`

3.  Activate the environment
    

`conda activate gulls`

4.  `conda install -c conda-forge cfitsio`
    

_**Installing from Source**_

1.  Go to [https://heasarc.gsfc.nasa.gov/fitsio/](https://heasarc.gsfc.nasa.gov/fitsio/)
    
2.  Find the version you want (if not the latest) and copy the hyperlink
    
3.  On Unity, use wget to download
    

`wget https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz`

4.  Untar + unzip
    
5.  Go into the expanded directory and do `./configure`, followed by `make`, followed by `make` `install`
    
    1.  Note: if you want it installed in a particular directory, give .`/configure` a prefix flag.
        

`./configure -- prefix=/home/crisp.92/local`

## GULLS

1.  Head to [https://github.com/gulls-microlensing/gulls](https://github.com/gulls-microlensing/gulls)
    
2.  Fork the repository. IMPORTANT: there will be a checkbox asking if you only want to fork main. Uncheck that.
    
3.  Clone the repo to Unity.
    
4.  Checkout the dev branch. 
    
    1.  Note: This _should_ remain checked out until you check out a different branch (even when you log out of Unity), but if you’re paranoid like me, you can see what branch you’re on with git branch --show-current
        
5.  In your ~/.bashrc file, add `export GULLS_BASE_DIR=/path/to/gulls/`. Don’t forget to source ~/.bashrc!!!
    
    1.  Also add `export GULLS_STARS_DIR=/path/to/catalogs/.` I think this should point to the directory that contains the star list files (e.g., gulls\*.lenses). All paths from there are set relative to their location is why I think this. Right now, the catalogs are in /fs/project/gaudi.1/aass/synthpop\_cats/Huston2023\_surot2d/
        
6.  There are a few things on the main branch that are missing from dev, and you can’t really have both branches active at the same time (AFAIK). I’ve added a copy of the main branch code (as of 2024-11-13) to `/fs/project/gaudi.1/aass/gulls-main-20241113`. You’ll need to…
    
    1.  Copy all the `gulls*.sh` files from `gulls-main-20241113/scripts` to your `gulls/scripts` directory.
        
    2.  Copy `photometry.cpp` and `photometry.h` from `gulls-main-20241113/src` to your `gulls/sr`c
        
7.  There are also files missing because of licensing reasons. Those aren’t in main either, but I’ve added them to `/fs/project/gaudi.1/aass/gulls_missing_files`. You’ll need to…
    
    1.  Copy `random.h` and `zroots2.h` to your `gulls/src/headers` directory
        
    2.  Copy `random.cpp` and `zroots2.cpp` to your `gulls/src/classes` directory
        
8.  NOWWWWWW, `cd $GULLS_BASE_DIR`
    
    1.  If this doesn’t work, are you sure you sourced your ~/.bashrc?
        
9.  Run `./configure.sh`
    
10.  `cd src`
    
11.  Open the `standardPlanet.cpp` file, go to line 22, and correct the spelling of “parameterization.” Save and exit.
    
12.  Copy the `ESPL.tbl` file from your `VBMicrolensing/VBMicrolensing/data` directory to your `gulls/src` directory.
    
    1.  [https://drive.google.com/file/d/17pmLaMpjYXMp34JMwRueA1fFlKqWUGlu/view?usp=drive\_link](https://drive.google.com/file/d/17pmLaMpjYXMp34JMwRueA1fFlKqWUGlu/view?usp=drive_link)
        
13.  In the `Makefile`:
    
    1.  Change `CC` to g++ (~line 42)
        
    2.  Make sure `BASEDIR` is set to your `gulls/src` path
        
    3.  In `CFLAGS`, add `-I/path/to/your/VBMicrolensing` and `-I/path/to/your/cfitsio/include/`
        
    4.  In `LINKERFLAGS`, add `-L/path/to/your/VBMicrolensing` and `-L/path/to/your/cfitsio/lib/`, along with the flags `-lVBB` and `-lcfitsio`
        
14.  In `defaults.mk`:
    
    1.  Change `CXX` and `CC` to g++
        
    2.  In CPPFLAGS, add `-I/path/to/your/VBMicrolensing/lib/`, `-I/path/to/your/gulls/src/headers/`, and `-I/path/to/your/cfitsio/include/.`
        
15.  \[Optional\] Copy and replace with Samson’s `readParamfile.cpp` (in `…/aass/gulls_sj/src/readParafile.cpp`)
    
16.  `make gullsFish`
    
    1.  This adds a `gullsFish.x` file to the `bin` folder.
        
17.  Check it is working:
    
    1.  `$ cd ../bin`
        
    2.  `$ ./gullsFish.x`
        

`>>> Usage: gulls  -i <infile> -s <instance> {-f <field>} {-d}`

## Running Gulls

## Pre-Run Notes

-   Executable is made in $GULLS\_BASE\_DIR/bin
    
    -   Usage: gulls -i <infile> -s <instance> {-f <field>} {-d}
        
    -   infile - parameters
        
    -   instance - subrun
        
    -   field - field number
        
    -   \-d - debug, 0/1
        
-   Need to create a job submission script that runs your jobs with slurm, e.g. gullsLaunch.array.sh
    
    -   Make a shell script; e.g., `gullsLaunch.slurmarray.sh:`
        
        -   Does the resource allocation
            
        -   Runs many gulls instances using a shell script that it writes; e.g., `run_<runname>.sh:`
            
        -   Sets up the environment etc. (`gullsPreamble.sh)`
            
        -   Gulls call uses the executable named in the parameter file
            
        -   Moves the outputs
            
    -   Example scripts in `aass/gulls_sj/scripts/`
        

## Parts of a Run

1.  Setup parameter file
    
    1.  NOTE: In the parameter file, one does not simply comment things out normally. If it can find the string it’s looking for at all, it will overwrite things. So if you want to comment out an old version of the parameter `OUPUT_ONDET`, you can’t just do `#OUTPUT_ONDET`, you’ll need to take out a character or something like `#UTPUT_ONDET`.
        
    2.  Example in `aass/template_unity.prm`
        
        1.  Edit the `RUN_NAME,` on line 3
            
    3.  Make sure directories and executable match.
        
        1.  lenses,  planets,  rates,  sources, and  starfields are shared directories on `/fs/project/gaudi.1/gulls/` 
            
        2.  `gullsFish.x`
            
        3.  Observatories and weather from …/`ass/gulls_inputs*.tar.gz`
            

`$ tar -xvf gulls_inputs*tar.gz`

`$ rm -r starfields`

`$ cp /fs/poject/gaudi.1/aass/gulls_sj/observatories/romanc7.list $GULLS_BASE_DIR/observatories`

2.  Setup launch script
    
    1.  Most set up is done in the parameter file, but the allocation-related parts are changed in echo statements in this file.
        
3.  Run gullsSetupRun.sh
    
    1.  make sure to add an `SRCDIR=’/path/to/gulls/bin/’` to `gullsSetupRun.sh`.
        
    2.  calls gullsPreamble.sh, which parses parameter file (`.prm` file)
        
    3.  Sets up directory structure for outputs (`$final_dir`)
        
    4.  Copies dependencies (eather, param file, observatory files) to `$final_dir/logs`
        
4.  Creating planet catalogs: `gen_planets.pl`. Ali has made some; they are in `/fs/projects/gaudi/gulls/planets/test/`
    
    1.  Needs to be run in the directory where you would like the planets to be (or copied). E.g., for the _test_ set, run in `/fs/project/gaudi.1/gulls/planets`
        
    2.  Specify l, b min and max; m min and max; a min & max
        
    3.  Scripts checks in `*.sources` (list of source files), and if in l,b range then makes a corresponding planet file. 
        
    4.  `nr=1` - number of subruns (i.e. `nr=2` makes two planet files for each l,b set)
        
    5.  `nl=1000` - number of rows in each output file
        
    6.  Save them in the data folder
        

## Reducing GULLS on Unity @ OSU

1.  Make a conda environment (I used the gulls one with cfitsio) with pandas and pytables or install those in an existing environment. Activate that environment.
    
2.  `$ cd $GULLS_BASE_DIR/scripts`
    
3.  If you used `runname.prm` to launch your gulls run, then  
    `$ python reduce_gulls.py $GULLS_BASE_DIR/parameterFiles/runname.prm --in-raw` 
    
    1.  The `--in-raw` just sets where the reduction script looks for the output data
        
    2.  There are many other flags you may consider setting, but for now we should just need that one
        
    3.  Importantly, we need to figure out a real covfac
        
4.  This will create files in `$OUTPUT_DIR/analysis/*hdf5` that can be used. These will have the event weights normalized.
    
    1.  The `*out*hdf5` is the output of all the events simulated, and the `*det*hdf5` is the output of all the DETECTED events simulated (likely just DeltaChi2>160, which this limit can be set as a flag in reduce\_gulls.py I think)
        
    2.  I tend to work with the _\*_out\* file which allows for flexible detection cuts later in the analysis
        
5.  A lot of analysis will use weighted histograms, so here are some useful tips I used a lot  
    `# read in the data`  
    `data = pandas.read_hdf(‘data.out.hdf5’)`  
    `# place detection limits`  
    `det = data[data['ObsGroup_0_chi2']>=160]`  
    `# just sample bins for log(mass) of planets`  
    `mass_bins = np.linspace(-2,1,100)`  
    `# gets center of bins for plots`  
    `mass_bin_centers = (mass_bins[1:]+mass_bins[:-1])/2`  
    `# make weighted histogram using the filtered data`  
    `# note planet masss are in Solar mass units, so you need mearth~3e-6 defined`  
    `mass_hist, mass_bins = np.histogram(np.log10(det[‘Planet_mass’]/mearth),`  
    `bins = mass_bins,`
    

`weights = data[‘final_weights’])`  
`# loose plot call`  
`Fig, ax = plt.subplots()`  
`# plot log(counts) against log(mass) bin centers`  
`ax.plot(mass_bin_centers, np.log10(mass_hist))`