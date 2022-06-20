#! /bin/bash -login
#SBATCH -J UVM
#SBATCH -t 4:00:00
#SBATCH --mem=30G
#SBATCH -p all
#SBATCH -c 1
#SBATCH -o %x.out

conda activate /cluster/home/pluo/bin/miniconda3/envs/medremix

Rscript ~/Project/cfMeDIP/UVM/merge_data.R
Rscript ~/Project/cfMeDIP/UVM/map_signature.R
