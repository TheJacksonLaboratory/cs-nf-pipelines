#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from '../bin/help/nygc.nf'
include {param_log} from '../bin/log/nygc.nf'
include {TRIM_GALORE} from '../modules/trim_galore/trim_galore'
include {CONCATENATE_READS_PE} from '../modules/utility_modules/concatenate_reads_PE'
include {BWA_MEM} from '../modules/bwa/bwa_mem'
include {SHORT_ALIGNMENT_MARKING} from '../modules/nygc-short-alignment-marking/short_alignment_marking'
include {GATK_FIX_MATE_INFORMATION} from '../modules/gatk/gatk_fix_mate_information'
include {PICARD_MARKDUPLICATES}	from '../modules/picard/picard_markduplicates'
include {READ_GROUPS} from '../modules/utility_modules/read_groups'
include {GATK_REALIGNERTARGETCREATOR} from '../modules/gatk/gatk_realignertargetcreator'
include {GATK_BASERECALIBRATOR} from '../modules/gatk/gatk_baserecalibrator'
include {GATK_APPLYBQSR} from '../modules/gatk/gatk_applybqsr'
include {GATK_INDELREALIGNER} from '../modules/gatk/gatk_indelrealigner'
include {PICARD_COLLECTWGSMETRICS} from '../modules/picard/picard_collectwgsmetrics'
include {CONPAIR} from '../modules/conpair/conpair'
