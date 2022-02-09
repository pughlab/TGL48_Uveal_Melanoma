INPUTDIR=/cluster/projects/pughlab/external_data/TGL48_Uveal_Melanoma/TGL48_sWG/bam
shdir=/cluster/projects/pughlab/projects/uveal_melanoma/insert_size/sh_scripts/swgs
outdir=/cluster/projects/pughlab/projects/uveal_melanoma/insert_size/output/swgs
picard_dir=/cluster/tools/software/picard/2.10.9

mkdir -p $shdir
mkdir -p $outdir

cd $INPUTDIR
ls *Ct*.bam > bams
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}_reads.sh
#!/bin/bash
#
module load samtools

samtools view -c -F 260 $INPUTDIR/${bam}.bam > $outdir/${bam}_reads
EOF

done

cd $shdir

ls *reads.sh > files
for file in $(cat files); do
sbatch -c 1 -t 2:00:00 -p all --mem 4G $file
done
