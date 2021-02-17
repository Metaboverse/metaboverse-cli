#!/bin/bash
#SBATCH --time=7:00:00
#SBATCH --nodes=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/rutter-group1/j-berg/slurmjob-logs/slurmjob-%j
#SBATCH --account=rutter-gpu-np
#SBATCH --partition=rutter-gpu-np
#SBATCH --mem=0

# Rutter lab GPU node specs
# - 40 cores
# - 192 GB of memory
# - 4x RTX2080TI GPUs
# - 168hr max

# Set instance variables
printf "+ Setting environment...\n"
MY_PATH=/uufs/chpc.utah.edu/common/home/rutter-group1/j-berg/programs/metaboverse-cli
SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR
cd $SCRDIR

# Activate conda environment
# source /uufs/chpc.utah.edu/common/home/u0690617/miniconda3/etc/profile.d/conda.sh
#conda remove -y -n metaboverse_cli
#conda create -y -n metaboverse_cli
# source activate metaboverse_cli
#conda config --add channels bioconda:
#conda config --add channels conda-forge
#conda install -y python
#conda install -y conda
#conda install -y pyinstaller
#conda install -y pandas numpy scipy scikit-learn matplotlib networkx requests

# Build current version of metaboverse-cli
#pyinstaller $MY_PATH/metaboverse-linux.spec

printf "+ Version info:\n"
$MY_PATH/dist/metaboverse-cli-linux -v

# Get species IDs from Reactome
REACOME_API="https://reactome.org/ContentService/data/species/all"
SPECIES=( $( curl -s $REACOME_API | jq -r '.[].abbreviation' ) )

printf "Processing database curation for:\n"
for X in ${SPECIES[@]};
  do mkdir -p $SCRDIR/${X} ;
done

# Run
parallel $MY_PATH/dist/metaboverse-cli-linux curate -o $SCRDIR/{} -s {} ::: "${SPECIES[@]}"

printf "+ Processing complete...\n"

printf "+ Outputing metadata\n"
cp /uufs/chpc.utah.edu/common/home/rutter-group1/j-berg/slurmjob-logs/slurmjob-$SLURM_JOBID $SCRDIR

printf "Metadata for bulk Metaboverse .mvdb curation:\n" >> $SCRDIR/README.txt

printf "\nDate: " >> $SCRDIR/README.txt
date '+%Y-%m-%d %H:%M:%S' >> $SCRDIR/README.txt

printf "\nMetaboverse version: " >> $SCRDIR/README.txt
$MY_PATH/dist/metaboverse-cli-linux -v >> $SCRDIR/README.txt

printf "\nReactome version: " >> $SCRDIR/README.txt
curl -X GET --header 'Accept: text/plain' 'https://reactome.org/ContentService/data/database/version' >> $SCRDIR/README.txt

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


# Afterwards, upload to sourceforge
# $ cd $SCRDIR
# $ scp -r */*.mvdb j-berg@frs.sourceforge.net:/home/frs/project/metaboverse/mvdb_files/x.y.z-tag
# scp README.txt j-berg@frs.sourceforge.net:/home/frs/project/metaboverse/mvdb_files/x.y.z-tag
# scp slurmjob-xxxxxxx j-berg@frs.sourceforge.net:/home/frs/project/metaboverse/mvdb_files/x.y.z-tag
