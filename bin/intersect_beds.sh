#!/usr/bin/env bash

name_string="$1"

    bedtools window -w 100 -a ins.bed -b /ref_data/mgpv5.SV_insertions.bed.gz > ${name_string}.ins.s.bed
    bedtools window -w 100 -a ins.bed -b /ref_data/mus_musculus_insertions_no_nstd74.gvf > ${name_string}.ins.e.bed
    bedtools window -w 100 -a del.bed -b /ref_data/mgpv5.SV_deletions.bed.gz > ${name_string}.del.s.bed
    bedtools window -w 100 -a del.bed -b /ref_data/mus_musculus_deletions_no_nstd74.gvf > ${name_string}.del.e.bed
    bedtools window -w 100 -a inv.bed -b /ref_data/mus_musculus_inversions_no_nstd74.gvf > ${name_string}.inv.e.bed
    bedtools window -w 100 -a dup.bed -b /ref_data/mus_musculus_copy_number_gain_no_nstd74.gvf > ${name_string}.dup.e.bed
    bedtools window -w 100 -a tra.bed -b /ref_data/mus_musculus_complex_structural_alteration_no_nstd74.gvf > ${name_string}.tra.e.bed
    bedtools intersect -a ins.bed -b /ref_data/mm10.GeneSymbol.sorted.nochr.bed -wa -wb > ${name_string}.ins.genes.bed
    bedtools intersect -a del.bed -b /ref_data/mm10.GeneSymbol.sorted.nochr.bed -wa -wb > ${name_string}.del.genes.bed
    bedtools intersect -a inv.bed -b /ref_data/mm10.GeneSymbol.sorted.nochr.bed -wa -wb > ${name_string}.inv.genes.bed
    bedtools intersect -a dup.bed -b /ref_data/mm10.GeneSymbol.sorted.nochr.bed -wa -wb > ${name_string}.dup.genes.bed
    bedtools intersect -a tra.bed -b /ref_data/mm10.GeneSymbol.sorted.nochr.bed -wa -wb > ${name_string}.tra.genes.bed
    bedtools intersect -a ins.bed -b /ref_data/mm10_exons -wa -wb | \
        cut -f 1,2,3,4,5,6,7,8,10,12,14,16 | \
        sed -e 's/\"//g;s/\;//g' > ${name_string}.ins.exons.bed
    bedtools intersect -a del.bed -b /ref_data/mm10_exons -wa -wb | \
        cut -f 1,2,3,4,5,6,7,8,10,12,14,16 | \
        sed -e 's/\"//g;s/\;//g' > ${name_string}.del.exons.bed
    bedtools intersect -a inv.bed -b /ref_data/mm10_exons -wa -wb | \
        cut -f 1,2,3,4,5,6,7,8,10,12,14,16 | \
        sed -e 's/\"//g;s/\;//g' > ${name_string}.inv.exons.bed
    bedtools intersect -a dup.bed -b /ref_data/mm10_exons -wa -wb | \
        cut -f 1,2,3,4,5,6,7,8,10,12,14,16 | \
        sed -e 's/\"//g;s/\;//g' > ${name_string}.dup.exons.bed
    bedtools intersect -a tra.bed -b /ref_data/mm10_exons -wa -wb | \
        cut -f 1,2,3,4,5,6,7,8,10,12,14,16 | \
        sed -e 's/\"//g;s/\;//g' > ${name_string}.tra.exons.bed