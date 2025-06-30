/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Validate input parameters
// Check input path parameters to see if they exist
// Check mandatory parameters
if (params.reference) {
    if (!params.ref) {
        println " ERROR: You didn't set any folder to store all reference files \
        Please set via --ref and try again (: "
        System.exit(1)
    }

    if (!params.remote_vcf_location) {
        println "ERROR: You didn't set remote reference path to prepare reference fot harmonization"
        println "Please set inusing --remote_vcf_location to set it"
        System.exit(1)
    }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include {prepare_reference} from '../subworkflows/local/prepare_reference'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// this is a top level workflow. 
// it does not have any take or main blocks
// it evokes a subworkflow that will have these
workflow PREPARE_REFERENCE_MASTER{
    // MODULE: check reference
    // Channel is a data strcture in Nextflow. essentially stores data between processes
    // map apply as a transformation to each element in the channel
    // in this case concat "chr" with the current item (it)
    ch_chrom=Channel.from(params.chrom).flatten().map{"chr"+it}
    // ch_chrom looks like: [chr1,chr2,chr3...]
    prepare_reference(ch_chrom)
}