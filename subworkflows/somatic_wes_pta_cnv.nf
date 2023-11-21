#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {GATK_PRINTREADS as GATK_PRINTREADS_NORMAL;
         GATK_PRINTREADS as GATK_PRINTREADS_TUMOR;} from "${projectDir}/modules/gatk/gatk_printreads"
include {SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_NORMAL;
         SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_TUMOR} from "${projectDir}/modules/samtools/samtools_mpileup"
include {SEQUENZA_PILEUP2SEQZ} from "${projectDir}/modules/sequenza/sequenza_pileup2seqz"
include {SCARHRD} from "${projectDir}/modules/scarhrd/scarhrd"
include {SEQUENZA_RUN} from "${projectDir}/modules/sequenza/sequenza_run"
include {SEQUENZA_ANNOTATE as SEQUENZA_ANNOTATE_RAW;
         SEQUENZA_ANNOTATE as SEQUENZA_ANNOTATE_FILTERED} from "${projectDir}/modules/sequenza/sequenza_annotate"
include {SEQUENZA_NA_WINDOWS} from "${projectDir}/modules/sequenza/sequenza_na_window"
include {BEDTOOLS_SUBTRACT} from "${projectDir}/modules/bedtools/bedtools_sequenza_subtract"

workflow CNV {
    take:
        ch_paired_samples
    main:

        ch_ind_samples = ch_paired_samples
            .multiMap{it -> 
                    normal: ["${it[1].patient}--${it[1].normal_id}".toString(), it[1], it[2], it[3]]
                    tumor:  ["${it[1].patient}--${it[1].tumor_id}".toString(), it[1], it[5], it[6]]
                    }
        ch_normal_samples = ch_ind_samples.normal.unique{it[0]}
        ch_tumor_samples  = ch_ind_samples.tumor.unique{it[0]}

        GATK_PRINTREADS_NORMAL(ch_normal_samples)
        GATK_PRINTREADS_TUMOR(ch_tumor_samples)

        SAMTOOLS_MPILEUP_NORMAL(GATK_PRINTREADS_NORMAL.out.bam_bai)
        SAMTOOLS_MPILEUP_TUMOR(GATK_PRINTREADS_TUMOR.out.bam_bai)

        sequenza_input = SAMTOOLS_MPILEUP_NORMAL.out.pileup.join(SAMTOOLS_MPILEUP_TUMOR.out.pileup, by: 1) // join on metadata field
                         .map{it -> [it[0].id, it[0], it[2], it[1], it[4], it[3]]}

        SEQUENZA_PILEUP2SEQZ(sequenza_input)

        SCARHRD(SEQUENZA_PILEUP2SEQZ.out.seqz)

        SEQUENZA_RUN(SEQUENZA_PILEUP2SEQZ.out.seqz)

        SEQUENZA_ANNOTATE_RAW(SEQUENZA_RUN.out.segments)

        SEQUENZA_NA_WINDOWS(SEQUENZA_RUN.out.extract_rdata)

        subtract_input = SEQUENZA_RUN.out.segments_tmp.join(SEQUENZA_NA_WINDOWS.out.na_windows)

        BEDTOOLS_SUBTRACT(subtract_input)

        SEQUENZA_ANNOTATE_FILTERED(BEDTOOLS_SUBTRACT.out.segments_filtered)
}
