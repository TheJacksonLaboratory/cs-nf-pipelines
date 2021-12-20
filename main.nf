#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { $"{params.workflow}" } from './workflows/$"{params.workflow_file}"'
includeConfig $"{params.config}"
include {find} from './bin/shared/groovy_helper'

workflow{
  $"{params.workflow}"()
}
