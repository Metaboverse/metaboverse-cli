#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -o ~/slurmjob-%j
#SBATCH --partition=notchpeak


# Set curation version
VERSION="0.10.1"


# Set instance variables
printf "+ Setting environment...\n"
MY_PATH=/uufs/chpc.utah.edu/common/home/rutter-group1/j-berg/programs/metaboverse-cli
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

# Activate conda environment
source /uufs/chpc.utah.edu/common/home/u0690617/miniconda3/etc/profile.d/conda.sh
source activate pyinstaller
conda update -n base -c defaults conda -y
conda update --all -y

# Build current version of metaboverse-cli
cd $MY_PATH
pyinstaller $MY_PATH/metaboverse-cli.spec

cd $SCRDIR

printf "+ Version info:\n"
$MY_PATH/dist/metaboverse-cli-linux -v

# Get species IDs from Reactome
REACTOME_API="https://reactome.org/ContentService/data/species/all"
SPECIES=( $( curl -k $REACTOME_API | jq -r '.[].abbreviation' ) )

printf "Processing database curation for:\n"
for X in ${SPECIES[@]};
  do mkdir -p $SCRDIR/${X} ;
done

# Run
parallel $MY_PATH/dist/metaboverse-cli-linux curate --force_new_curation --output $SCRDIR/{} --organism_id {} ::: "${SPECIES[@]}"

printf "+ Processing complete...\n"

printf "+ Outputting metadata\n"
cp /uufs/chpc.utah.edu/common/home/rutter-group1/j-berg/slurmjob-logs/slurmjob-$SLURM_JOBID $SCRDIR

printf "Metadata for bulk Metaboverse .mvdb curation:\n" >> $SCRDIR/README.txt

printf "\nDate: " >> $SCRDIR/README.txt
date '+%Y-%m-%d %H:%M:%S' >> $SCRDIR/README.txt

printf "\nMetaboverse version: " >> $SCRDIR/README.txt
$MY_PATH/dist/metaboverse-cli-linux -v >> $SCRDIR/README.txt

printf "\nReactome version: " >> $SCRDIR/README.txt
curl -kX GET https://reactome.org/ContentService/data/database/version >> $SCRDIR/README.txt

printf "\n\nOrganisms curated:" >> $SCRDIR/README.txt
printf "\nSTART" >> $SCRDIR/README.txt

for X in ${SPECIES[@]}; do
  if [ -f "$SCRDIR/${X}/${X}.mvrs" ]; then
    printf "\n ${X}" >> $SCRDIR/README.txt
    rm $SCRDIR/${X}/${X}.mvrs
  else
    rm -rf $SCRDIR/${X}
  fi
done

printf "\nEND" >> $SCRDIR/README.txt
printf "\n"

# Afterwards, upload to host
mkdir -p /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/mvdb
mkdir -p /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/mvrs
mkdir -p /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/nbdb

cd $SCRDIR
cp -r */*.mvdb /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/mvdb
cp README.txt /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/mvdb

cp -r */*_template.mvrs /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/mvrs
cp README.txt /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/mvrs

cp -r */*.nbdb /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/nbdb
cp README.txt /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION/nbdb

conda list > build_versions.txt
cp build_versions.txt /uufs/chpc.utah.edu/common/home/rutter-website/html/Metaboverse/source/v$VERSION
