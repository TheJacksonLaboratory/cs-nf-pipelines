#!/usr/bin/env bash

bam_string="$1"
out_string="$2"

pbsv discover --tandem-repeats /ref_data/mm10_ucsc_simple_tandem_repeats_ucsc.bed ${1} ${2}