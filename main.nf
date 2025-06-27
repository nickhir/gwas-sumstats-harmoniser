#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// import workflow 
include { PREPARE_REFERENCE_MASTER } from './workflows/PREPARE_REFERENCE_MASTER'
include { GWASCATALOGHARM } from './workflows/gwascatalogharm'

//
// This file gets executed when nextflow run is called.
//
workflow NFCORE_GWASCATALOGHARM {

    // Check mandatory parameters
    // if the user doesnt provide any values, they will be set to null
    params.reference = null
    params.gwascatalog = null
    params.harm = null


    if (!params.chrom) {
    println "ERROR: You didn't set chromsomes to be harmnnised"
    println "Please set --chrom 22 or --chromlist 22,X,Y or set chrom in ./config/default_params.config "
    System.exit(1)
    }


    if (!params.reference) {
        if (!params.to_build) {
            println "ERROR: You didn't set the target build to harmonise to"
            println "Please set --to_build 38"
            System.exit(1)
            }
        
        if (!params.threshold) {
            println "ERROR: You didn't set threshold to imput the direction of palindromic variants"
            println "Please set --threshold 0.99 or set threshold in ./config/default_params.config "
            System.exit(1)
            }
        
        if (!params.version) {
            println " ERROR: Please specific the pipeline version you are running (e.g. v1.1.9) \
            Please set --version and try again (: "
            System.exit(1)
            }   
    }

    // if the user provides a reference, we will not run the harmonisation workflow
    // but the reference preparation workflow
    if (params.reference) {
        println ("Prepare the reference ...")
        PREPARE_REFERENCE_MASTER()
    } 
    else if (params.harm) {
        if (!params.file && !params.list) { 
            println " ERROR: You didn't set any files to be harmonised \
                Please set --file for a single input file or \
                set --list for a list containing all files are waiting to be harmonised \
                and try again (: "
                System.exit(1)
        } 
        else {
            println ("Start harmonising files")
            GWASCATALOGHARM()
        }
    } 
    else {
            println " ERROR: You didn't set any model to run the pipeline \
            Please set --harm and try again (: "
            System.exit(1)
            }
}


// Convention to have it declared like this
workflow {
    NFCORE_GWASCATALOGHARM ()
}
