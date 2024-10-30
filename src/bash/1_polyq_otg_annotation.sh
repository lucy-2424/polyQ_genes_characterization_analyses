#!/bin/bash

DATA=/research/2023_polyQ/otg/data
SOFTWARE=/software/ensembl-vep 
CACHE=/research/ensembl_cache 


$SOFTWARE/vep -i $DATA/l2g_combined_no_annot_sorted.vcf \
        --dir_cache $CACHE \
        --offline \
        --most_severe \
	--per_gene \
	--nearest gene \
	--force_overwrite \
        -o $DATA/l2g_combined_annot.vcf
