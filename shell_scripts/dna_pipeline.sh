#!/bin/bash

module load perl

perl /cluster/projects/pughlab/projects/CHARM/pipeline-suite/pughlab_dnaseq_pipeline.pl \
-t /cluster/projects/pughlab/projects/uveal_melanoma/configs/TGL48_pipeline_swgs.yaml \
-d /cluster/projects/pughlab/projects/uveal_melanoma/configs/TGL48_bams_swgs.yaml \
--variant_calling \
--summarize \
--create_report \
-c slurm \
--remove
