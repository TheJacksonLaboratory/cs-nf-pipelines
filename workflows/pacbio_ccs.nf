#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
//include {help} from '../bin/help/rnaseq'
include {param_log} from "${projectDir}/bin/log/pacbio"
include {BUILDPBMM2INDEX;
         PBMM2MAPCCS} from "${projectDir}/modules/pbmm2"

include {DISCOVERTANDEM;
         CALLCCS} from "${projectDir}/modules/pbsv"

include {SNIFFLES} from "${projectDir}/modules/sniffles"

include {MERGEPACBIO;
         ANNOTATE;
         SUMMARIZE;
         PREPBEDS;
         INTERSECTBEDS;
         SUMMARIZEINTERSECTIONS;
         ANNOTATEPACBIO} from "${projectDir}/modules/survivor"

// log paramater info
param_log()

workflow PACBIO_CCS {
    // Step 1: Prepare index
    BUILDPBMM2INDEX(params.names, params.fasta)

    // Step 2: Map CCS reads to indexed genome
    PBMM2MAPCCS(params.names, params.fastq1, BUILDPBMM2INDEX.out.pbmm2_index)

    // Step 3: Call SV with PBSV
    DISCOVERTANDEM(PBMM2MAPCCS.out.pbmm2_ccs)
    CALLCCS(DISCOVERTANDEM.out.pbsv_svsig, params.fasta)

    // Step 4: Call SV with Sniffles
    SNIFFLES(PBMM2MAPCCS.out.pbmm2_ccs)

    // Step 5: SURVIVOR merge variant calls
    MERGEPACBIO(CALLCCS.out.pbsv_vcf, SNIFFLES.out.sniffles_vcf, params.surv_dist, params.surv_supp, params.surv_type, params.surv_strand, params.surv_min)

    // Step 6: Annotate SURVIVOR results
    ANNOTATE(MERGEPACBIO.out.survivor_vcf, "pacbio")
    SUMMARIZE(MERGEPACBIO.out.survivor_vcf)
    PREPBEDS(ANNOTATE.out.annot, SUMMARIZE.out.summary)
    INTERSECTBEDS(PREPBEDS.out.sv_beds)
    SUMMARIZEINTERSECTIONS(PREPBEDS.out.sv_beds, INTERSECTBEDS.out.intersected_beds, SUMMARIZE.out.summary, ANNOTATE.out.annot)
    ANNOTATEPACBIO(MERGEPACBIO.out.survivor_vcf, INTERSECTBEDS.out.intersected_exons)
}