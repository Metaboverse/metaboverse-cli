#!/bin/bash
#SBATCH --time=168:00:00
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
source /uufs/chpc.utah.edu/common/home/u0690617/miniconda3/etc/profile.d/conda.sh
conda remove -y -n metaboverse_cli
conda create -y -n metaboverse_cli
source activate metaboverse_cli
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y python
conda install -y conda
conda install -y pyinstaller
conda install -y pandas numpy scipy scikit-learn matplotlib networkx requests

# Build current version of metaboverse-cli
pyinstaller $MY_PATH/metaboverse-linux.spec

echo "+ Version info:"
$MY_PATH/dist/metaboverse-cli-linux -v

# Get species IDs from Reactome
REACOME_API="https://reactome.org/ContentService/data/species/all"
SPECIES=( $( curl -s $REACOME_API | jq -r '.[].abbreviation' ) )
echo 'Processing database curation for: '
for X in ${SPECIES[@]};
  do echo '${X}' ;
done

# Run
echo "+ Running scripts..."
parallel "$MY_PATH/dist/metaboverse-cli-linux -o $SCRDIR -s " "{}" ::: "${SPECIES[@]}"
echo "+ Processing complete..."

echo "+ Outputing metadata"
cp $HOME/slurm_output/slurmjob-$SLURM_JOBID $SCRDIR
echo '' > $SCRDIR-metadata.txt
