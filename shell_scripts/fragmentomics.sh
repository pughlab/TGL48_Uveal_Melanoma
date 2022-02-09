input=/cluster/projects/pughlab/external_data/TGL48_Uveal_Melanoma/TGL48_sWG/bam_deduplicated
shdir=/cluster/projects/pughlab/projects/uveal_melanoma/fragmentomics/sh_scripts/
outdir=/cluster/projects/pughlab/projects/uveal_melanoma/fragmentomics/output/
DELFI=/cluster/projects/pughlab/bin/fragmentomics

mkdir -p $shdir
mkdir -p $outdir

cd $input
ls TGL48*Ct*.bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

mkdir -p sh_scripts/
mkdir -p output/

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}.sh
#!/bin/bash
#
#$ -cwd

module load R/4.0.0
Rscript $DELFI/runFrag.R\
 --id $bam\
 --bamdir $input\
 --filters $DELFI/extdata/filters.hg38.rda\
 --gaps $DELFI/extdata/gaps.hg38.rda\
 --VNTRs $DELFI/extdata/VNTRs.hg38.rda\
 --tiles $DELFI/extdata/hg38_tiles.bed\
 --healthy $DELFI/extdata/healthy.median.hg38.rda\
 --outdir $outdir\
 --libdir $DELFI

EOF

done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 24:00:00 --mem 30G $file

done
