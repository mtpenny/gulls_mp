# gulls

A microlensing simulator optimized for space-based microlensing
surveys, but also supporting ground-based observatory simulations.   

## Quick links

- Astrometry options and outputs: documentation/source/astrometry.rst

## Requirements

1. [GNU Scientific
   Library](https://www.gnu.org/software/gnuastro/manual/html_node/GNU-Scientific-Library.html)
1. [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
1. [VBMicrolensing](https://github.com/valboz/VBMicrolensing)

## Installation Instructions

1. Install the above requirements if you don't already have them
1. Fork this repository, then clone your fork
1. Create environment variable with `export GULLS_BASE_DIR='/path/to/cloned/gulls/'`
    * It's recommended to appended this to your ~/.bashrc or
      ~/.bash_profile file
1. Run `./configure.sh` from the newly cloned gulls directory
1. `cd src` to change to the source code directory
1. Obtain the `random.cpp`, `zroots2.cpp` files (email penny1@lsu.edu)
   and put them in the `classes` directory. Put the `random.h` and
   `zroots2.h` files in the `headers` directory.
1. Copy the `ESPL.tbl` file from `VBMicrolensing/VBMicrolensing/data/`
   to the `src` directory (where you currently are)
1. Run `make <executable>` to build a chosen executable. `gullsFFP` is a
   good starting point, and it will create
   `${GULLS_BASE_DIR}/bin/gullsFFP.x`
1. If successfully compiled, running
   `${GULLS_BASE_DIR}/bin/gullsFFP.x` should give a message about the 


## (Incomplete) Checklist/Troubleshooting for regular gulls runs

- When preparing a run it's a good idea to do a short run with 
    ```
    PRETTY_PICS=1
    PRETTY_PICS_DIMENSIONS=128 #Or bigger if you'd like
    OUTPUT_LC=1 #If this is less than 1, it represents a probability of output
    OUTPUT_IMAGES=1
    OUTPUT_ONALL=1
    ```
    This should produce the maximum amount of diagnostic
    output. `OUTPUT_LC` determines the frequency with which lightcurves
    are output; `OUTPUT_ONALL` will output every event generated, but if
    set to zero, `OUTPUT_ONDET` will only output when a detection occurs
    (definition of detection depends on the type of run). `OUTPUT_IMAGES`
    saves images each time a lightcurve is output - one at peak
    magnification, one at baseline for each filter. Inspect both the
    images and the lightcurves. .det.lc files should show some level of
    detectable microlensing event in the lightcurve, .all.lc may
    not. Images should show a starfield, and depending on the peak
    magnification, blinking between peak and base images should show the
    microlensing event (requires Fpeak>~1.3 to be easily visible).

    - If the image is a white square, first make sure you select zscale in
ds9, and if that doesn't reveal stars it might indicate the presence
of bad magnitudes in the star, source, or lens catalog(s), or it's
possible that you just got unlucky and had a super-bright star land in
the image - the larger the PRETTY_PIC, the less likely it is that such
a star will blow out the image. Try generating a very large pretty pic
and see how many super-bright stars there are (they will look like
white squares) -- too many of these and its probably worth looking at
the color magnitude diagrams of the star lists to see that there's
nothing weird in there. 
- You can run with different levels of debug information by repeating -d flags
    - If you get an message about no valid stars, check that you have
	`NFILTERS` set correctly. 

