#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { info_score } from './modules/bgen.nf'
include { generate_negative_control } from './modules/snps.nf'

workflow {
    info_score()

    generate_negative_control(info_score.out)
}
