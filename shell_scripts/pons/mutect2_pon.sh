#!/bin/bash

module load perl

perl /cluster/projects/pughlab/projects/CHARM/pipeline-suite/scripts/mutect2.pl \
--create-panel-of-normals \
-t /cluster/projects/pughlab/projects/uveal_melanoma/configs/TGL48_pipeline_all_unique.yaml \
-d /cluster/projects/pughlab/projects/uveal_melanoma/configs/TGL48_bams_all_unique.yaml \
-o /cluster/projects/pughlab/projects/uveal_melanoma/pipeline_output/panel_of_normals/MuTect2 \
-c slurm \
--remove \
--no-wait
