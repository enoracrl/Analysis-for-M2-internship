#! /bin/bash

# sbatch spec
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --mail-user=enora.corler@etudiant.univ-rennes1.fr
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Load environments
. /local/env/envconda.sh
conda activate /home/genouest/cnrs_umr6290/ecorler/my_env

# Launch IsoformSwitchAnalyzeR
# $1 = gtf
# $2 = transcriptome
# $3 = design
# $4 = matrix
# $5 = comparisons

# sbatch launch_ISA.sh /groups/dog/stage/enora/isoformswitch/Data/extended_annotations_filter.full.gtf /groups/dog/stage/enora/isoformswitch/Data/transcriptome_lncrna-resist.fa /groups/dog/stage/enora/isoformswitch/Data/design_lncrnaResist_allSamples.csv /groups/dog/stage/enora/isoformswitch/Data/counts_transcript_filter.full.txt /groups/dog/stage/enora/isoformswitch/Data/comparisons.txt
    
Rscript /groups/dog/stage/enora/isoformswitch/IsoformSwitchAnalyzeR.r

conda deactivate