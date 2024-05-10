#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { info_score } from './modules/bgen.nf'
include { get_maf_match } from './modules/snps.nf'

workflow {
    info_score()

    get_maf_match(info_score.out)
}
