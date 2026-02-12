#!/bin/bash -x

if [ $# -ne 3 ]; then
    echo "Usage: $0 <paramfile> <field_list> <subrun>"
    exit
fi

export GULLS_BASE_DIR=/project/penny/gulls/

gullsbin=~/gulls_mp/bin

paramfile=$1
fieldlist=$2
subrun=$3

source ~/gulls_fz/scripts/gullsPreamble.sh

echo RUNNAME $runname
echo OUTPUTDIR $outputdir
echo FINALDIR $finaldir
echo executable $executable
echo paramfile $paramfile
echo sflist $sfdir/$sflist
echo srclist $srcdir/$srclist
echo lenslist $lensdir/$lenslist

f=X


#exit
#runname=romandecam_high

#workdir=/var/scratch/penny/$runname
#destdir=/work/penny/gulls/$runname
#executable=$(perl -ne '/EXECUTABLE=(.+)/; print $1;' $GULLS_BASE_DIR/parameterFiles/$runname.prm)

cat $fieldlist | while read f; do
    read idx l b rest < <(grep "^$f " $srcdir/$srclist)
    echo $f $idx $l $b

    if [ $(grep "^$f " $srcdir/$srclist | wc -l) -eq 1 ] && [ $(grep "^$f " $lensdir/$lenslist | wc -l) -eq 1 ] && [ $(grep "^$f " $sfdir/$sflist | wc -l) -eq 5 ]; then
	echo $f $l $b
	echo $gullsbin/$executable -i $paramfile -s $subrun -f $f
	$gullsbin/$executable -i $paramfile -s $subrun -f $f > $finaldir/$runname/${runname}_${subrun}_${f}.gullsstdout 2>&1
	(cd $outputdir && tar -czf $finaldir/$runname/${runname}_${subrun}_${f}.fits.tar.gz $runname/${runname}_${subrun}_${f}_*.*.fits)
	rm $outputdir/$runname/${runname}_${subrun}_${f}_*.*.fits
	(cd $outputdir && tar -czf $finaldir/$runname/${runname}_${subrun}_${f}.lc.tar.gz $runname/${runname}_${subrun}_${f}_*.*.lc)
        rm $outputdir/$runname/${runname}_${subrun}_${f}_*.*.lc
	(cd $outputdir && tar -czf $finaldir/$runname/${runname}_${subrun}_${f}.fm.tar.gz $runname/${runname}_${subrun}_${f}_*.*.fm.*)
        rm $outputdir/$runname/${runname}_${subrun}_${f}_*.*.fm.*
	mv $outputdir/$runname/${runname}_${subrun}_${f}.* $finaldir/$runname/
    else
	echo field $f $l $b has incomplete starfields
    fi
done
