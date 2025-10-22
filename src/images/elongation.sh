#!/bin/bash

if [ $# -ne 3 ]; then
  echo "Usage: $0 <detectorlist> <coursefile> <field> <output root>"
  exit
fi

#~/src/images/colourImage $1 0.0 /Volumes/olympus/deimos/mabuls/starfields/EUCLID-h-faint/out-$2 0.000020 /Volumes/olympus/deimos/mabuls/starfields/EUCLID-h-moder1/out-$2 0.00300 /Volumes/olympus/deimos/mabuls/starfields/EUCLID-h-moder2/out-$2 0.003 /Volumes/olympus/deimos/mabuls/starfields/EUCLID-h-bright/out-$2 0.001000 $3

field=$2

#~/src/images/colourImage $1 0.0 \
#    /Users/penny/gulls/wfirst-sim/m99-p12/out-$field 0.05 \
#    /Users/penny/gulls/wfirst-sim/p12-p18/out-$field 0.001 \
#    /Users/penny/gulls/wfirst-sim/p18-p24/out-$field 0.00005 \
#    /Users/penny/gulls/wfirst-sim/p24-p99/out-$field 0.0001 \
#	$3

sfdir=/Users/penny/gulls/starfields/Huston2023_surot2d/
nameroot=gulls_surot2d_H2023

echo ${0%/*}/elongation $1 32 \
						 $(grep "^$field " $sfdir/$nameroot.starfields | awk -v rt=$sfdir '{printf("%s%s %s ",rt,$(NF),$(NF-1))}') \
						 $3

#valgrind --tool=memcheck \
${0%/*}/elongation $1 32 \
						 $(grep "^$field " $sfdir/$nameroot.starfields | awk -v rt=$sfdir '{printf("%s%s %g ",rt,$(NF),$(NF-1))}') \
						 $3

						 #$sfroot/sf_extra_bright_gulls_surot2d_H2023/$field 1.000e-01 \
						 #$sfroot/sf_bright_gulls_surot2d_H2023/$field 1.000e-02 \
						 #$sfroot/sf_mid1_gulls_surot2d_H2023/$field 1.403e-03 \
						 #$sfroot/sf_mid2_gulls_surot2d_H2023/$field 1.000e-04 \
						 #$sfroot/sf_faint_gulls_surot2d_H2023/$field 1.000e-04 \
						 #$3
