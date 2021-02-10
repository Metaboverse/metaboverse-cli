#!/bin/bash
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/u0690617/slurm_output/slurmjob-%j
#SBATCH --account=rutter-gpu-np
#SBATCH --partition=rutter-gpu-np
#SBATCH --mem=0

# Rutter lab GPU node specs
# - 40 cores
# - 192 GB of memory
# - 4x RTX2080TI GPUs

# Set instance variables
echo "+ Setting environment..."
HOME=/uufs/chpc.utah.edu/common/home/$USER
MY_PATH=/uufs/chpc.utah.edu/common/home/$USER/programs/metaboverse-cli
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

echo "+ Version info:"
$MY_PATH/dist/metaboverse-cli-linux -v

# Get species IDs from Reactome
REACOME_API="https://reactome.org/ContentService/data/species/all"
SPECIES=( $( curl -s $REACOME_API | jq -r '.[].abbreviation' ) )
echo 'Processing database curation for: '
for X in ${SPECIES[@]};
  do mkdir -p $SCRDIR/${X} ;
done
for X in ${SPECIES[@]};
  do echo "${X}" ;
done

# Run
echo "+ Running scripts..."
parallel $MY_PATH/dist/metaboverse-cli-linux curate -o $SCRDIR/{} -s {} ::: "${SPECIES[@]}"
echo "+ Processing complete..."

echo "+ Outputing metadata"
cp $HOME/slurm_output/slurmjob-$SLURM_JOBID $SCRDIR

printf 'Metadata for bulk Metaboverse .mvdb curation:\n\n' >> $SCRDIR/README.txt

printf '\nDate:' >> $SCRDIR/README.txt
date '+%Y-%m-%d %H:%M:%S' >> $SCRDIR/README.txt

printf '\nMetaboverse version:' >> $SCRDIR/README.txt
$MY_PATH/dist/metaboverse-cli-linux -v >> $SCRDIR/README.txt

printf 'Reactome version:' >> $SCRDIR/README.txt
curl -X GET --header 'Accept: text/plain' 'https://reactome.org/ContentService/data/database/version' >> $SCRDIR/README.txt

printf '\n\nOrganisms curated:' >> $SCRDIR/README.txt
for X in ${SPECIES[@]};
  do printf '\n\t'${X} >> $SCRDIR/README.txt ;
done
printf '\n'

rm $SCRDIR/*.mvrs
