#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/prep_do_gbrs_input.nf"
include {param_log} from "${projectDir}/bin/log/prep_do_gbrs_input.nf"
include {DO_TRANSITION_PROBABILITIES} from "${projectDir}/modules/R/do_transition_probablities"
include {PARSE_TRANSITION_PROBABILITIES as PARSE_TRANSITION_PROBABILITIES_FEMALE;
        PARSE_TRANSITION_PROBABILITIES as PARSE_TRANSITION_PROBABILITIES_MALE} from "${projectDir}/modules/python/parse_transprobs"
include {PARSE_GENE_POSITONS} from "${projectDir}/modules/python/parse_gene_positions"
include {GENERATE_GRID_FILE} from "${projectDir}/modules/R/generate_grid_file"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// main workflow
workflow PREP_DO_GBRS_INPUT {
    DO_TRANSITION_PROBABILITIES()
    // Generate transition probablies in R

    generations = Channel.from( 0..params.num_generations )

    female_transProbs = DO_TRANSITION_PROBABILITIES.out.female_h5_file.combine(generations)
    PARSE_TRANSITION_PROBABILITIES_FEMALE(female_transProbs, 'F')

    male_transProbs = DO_TRANSITION_PROBABILITIES.out.male_h5_file.combine(generations)
    PARSE_TRANSITION_PROBABILITIES_MALE(male_transProbs, 'M')
    // For each generation, and sex, parse the h5 file to npz. 


    PARSE_GENE_POSITONS(DO_TRANSITION_PROBABILITIES.out.gene_list_tsv)
    // It is possible that genes are filtered out during the conversion from GTF to FASTA in g2gtools and prepare emase. 
    // The final `ref.gene_pos.ordered.npz` must only contain genes that are present in the reference set built with emase and found in: 'emase_gene2transcript'.

    GENERATE_GRID_FILE()
    //This scripts functions to space 'marker' positions every 0.02 cM along chr1:19 and X.
}
