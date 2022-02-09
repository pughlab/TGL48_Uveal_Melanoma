#!/bin/bash

griffin=/cluster/projects/pughlab/bin/Griffin
basedir=/cluster/projects/pughlab/projects/uveal_melanoma/griffin
sites=/cluster/projects/pughlab/bin/Griffin/sites/site_configs/DHS_sites.yaml
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
input=/cluster/projects/pughlab/external_data/TGL48_Uveal_Melanoma/TGL48_sWG/bam_deduplicated
counts=$basedir/output/GC_correction/repeat_masker.mapable.k50.Umap.hg38/GC_bias
outdir=$basedir/output/nucleosome_profiling/DHS
shdir=$basedir/sh_scripts/nucleosome_profiling

mkdir -p $outdir
mkdir -p $shdir

cd $input
ls *Ct*bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
mv bam bams

for bam in $(cat bams);do

name=${bam:0:25}
echo $bam
echo $name
echo -e "#!/bin/bash
source activate base
conda activate griffin" > $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_calc_coverage.py \
--sample_name $name \
--bam $input/${bam}.bam \
--GC_bias $counts/${name}.GC_bias.txt \
--background_normalization None \
--sites_yaml $sites \
--reference_genome $ref \
--results_dir $outdir \
--chrom_column Chrom \
--chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
--norm_window -5000 5000 \
--plot_window -1000 1000 \
--fragment_length 165 \
--step 15 \
--size_range 90 220 \
--map_quality 20 \
--strand_column Strand \
--individual False \
--smoothing True \
--num_sites none \
--sort_by none \
--ascending none \
--cpu 1" >> $shdir/${name}.sh

done

cd $shdir
ls *.sh > files
for file in $(cat files);do
sbatch -p all -c 1 --mem 8G -t 24:00:00 $file
done
