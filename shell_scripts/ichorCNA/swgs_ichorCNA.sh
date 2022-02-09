INPUTDIR=/cluster/projects/pughlab/external_data/TGL48_Uveal_Melanoma/TGL48_sWG/bam_deduplicated
base=/cluster/projects/pughlab/projects/uveal_melanoma/ichorCNA
outdir=$base/output/ichorCNA
wig=$base/output/wig
shdir=$base/sh_scripts/ichorCNA
ichorCNA=/cluster/projects/pughlab/bin/ichorCNA

mkdir -p $outdir
mkdir -p $shdir
mkdir -p $wig

cd $INPUTDIR
ls *.bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}.sh
#!/bin/bash
#
#$ -cwd

module load hmmcopy_utils/170718

readCounter --window 1000000 --quality 20\
 --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"\
 $INPUTDIR/${bam}.bam | sed "s/chrom=chr/chrom=/" > $wig/${bam}.wig

module load R/4.0.0
Rscript $ichorCNA/scripts/runIchorCNA.R\
 --id $bam\
 --WIG $wig/${bam}.wig\
 --ploidy "c(2,3)"\
 --normal "c(0.5,0.6,0.7,0.8,0.9,0.95)"\
 --maxCN 5\
 --gcWig $ichorCNA/inst/extdata/gc_hg38_1000kb.wig\
 --mapWig $ichorCNA/inst/extdata/map_hg38_1000kb.wig\
 --centromere $ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt\
 --normalPanel $ichorCNA/inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds\
 --minMapScore 0.9\
 --includeHOMD False\
 --genomeBuild "hg38"\
 --chrs "c(1:22, 'X')"\
 --chrTrain "c(1:22)"\
 --estimateNormal True\
 --estimatePloidy True\
 --estimateScPrevalence True\
 --scStates "c(1,3)"\
 --minSegmentBins 50\
 --altFracThreshold 0.05\
 --txnE 0.9999999\
 --txnStrength 10000000\
 --plotYLim "c(-1,1)"\
 --outDir $outdir

EOF

done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 24:00:00 --mem 4G $file

done
