#!/bin/bash
#SBATCH -N 1               # request one node
#SBATCH -n 20
#SBATCH -t 6:00:00	        # request two hours
#SBATCH -p workq          # in single partition (queue)
#SBATCH -A hpc_roman03

#SBATCH -o slurm-%j.out-%N # optional, name of the stdout, using the job number (%j) and the hostname of the node (%N)
#SBATCH -e slurm-%j.err-%N # optional, name of the stderr, using job and hostname values

# Set some handy environment variables.

jobn=0
#subrun=2

if ! [[ "$subrun" =~ ^[0-9]+$ ]]; then
    echo "Error: You must set subrun to an integer with --export=subrun=X after the sbatch command but before the script name, e.g. sbatch --export=subrun=X $0"
    exit
fi

if ! [[ "$paramfile" =~ ^.+\.prm$ ]]; then
    echo "Error: You must set paramfile to an parameter file name with --export=paramfile=X after the sbatch command but before the script name, e.g. sbatch --export=paramfile=X.prm $0"
    exit
fi

source ~/gulls_mp/scripts/gullsPreamble.sh

date

export GULLS_BASE_DIR=/project/penny/gulls/

if [ -z ${fields+x} ]; then
    fields=gbtdsfields;
fi
echo "Using $fields fields"


tmp=$(pwd)
if [ -d /var/scratch/$USER/$runname/ ]; then
    cd /var/scratch/$USER/
    rm -r $runname/
    cd $tmp
fi
mkdir -p /var/scratch/$USER/$runname/
mkdir -p $finaldir/$runname/

cd /var/scratch/$USER/$runname/
~/gulls_mp/run/clear_hanging.sh &

cd $finaldir
seq -f "%02g" 0 19 | parallel ~/gulls_mp/run/run_gulls.sh $paramfile ~/gulls_mp/run/fields/$fields.txt.{} $subrun

date

exit 0
