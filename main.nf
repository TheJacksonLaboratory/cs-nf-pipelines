#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { $"{params.workflow}" } from './workflows/$"{params.workflow_file}"'
includeConfig $"{params.config}"

workflow{
  $"{params.workflow}"()
}
