/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Validate input parameters
// Check input path parameters to see if they exist
// Check mandatory parameters
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include {  chr_check  } from '../subworkflows/local/check_reference.nf'
include {  main_harm  } from '../subworkflows/local/main_harm'
include { major_direction } from '../subworkflows/local/major_direction'
include {quality_control} from '../subworkflows/local/quality_control'

workflow GWASCATALOGHARM {
    
    main:
        // file is actually the path to the summary stat file that is supposed to be harmonised
        params.file=null
        
        // list is also a path to a file that contains paths. one per row.
        params.list=null
        
        if (params.file){
            println ("Harmonizing the file ${params.file}")
            files = Channel.fromPath(params.file).map{input_files(it)}
        }
        //files: GCST, path of ymal, path of GCST
        else if (params.list){
            println ("Harmonizing files in the file ${params.list}")
            // split into individual lines. remove any whitespace, convert each
            // row into a nextflow file object then run the custom input_files function
            files = Channel
                .fromPath(params.list)
                .splitText()map{it -> it.trim()}
                .map{row->file(row)}
                .map{input_files(it)}
        }
        /* MODULE: check reference
        ch_chrom looks like: [chr1,chr2,chr3...]
        chr_check() cross check required chromsomes with available reference
        ch_for_direction [chr1,chr2...] */

        // ch_input is emited by the chr_check subworkflow. 
        // so ch_for_direction takes that value 
        // the subworkflow chr_check operates on global parameters so does
        // not need any inputs
        // ch_for_direction is a list of chromosomes that are available in the reference
        // this does not report to the terminal because it just manipulates channel values. 
        // doesnt actually trigger any process
        ch_for_direction=chr_check().ch_input
        major_direction(ch_for_direction,files)

        //major_direction.out.direction_sum: [GCST, path of sum_count]
        //major_direction.out.hm_input: tuple val(GCST), val(palin_mode), val(status), val(chrom), path(merged), path(ref)
        harm_ch = major_direction.out.hm_input.groupTuple().transpose()
        main_harm(harm_ch,files)
        // out:[GCST009150, forward, path of harmonised.tsv]
        quality_control(main_harm.out.hm,major_direction.out.direction_sum,files,ch_for_direction,major_direction.out.unmapped,major_direction.out.alias_log)
}


// this file returns the original input and the
// corresponding metadata file
def input_files(input) {
    // splits the file name. baseName is everything before the first dot
    def baseName = input.getName().split("\\.")[0]

    // Check if the base name matches the pattern GCST[0-9]+
    def matcher = (baseName=~ /GCST\d+/).findAll() 
    if (matcher) {
        // Extract GCST ID using regex find
        println "yes,GCST"
        def gcstId = matcher[0]  // Get the first match
        return [gcstId, input+"-meta.yaml", input]
    } else {
        // Default case
        println "no,other setting"
        return [baseName, input+"-meta.yaml", input]
    }
}
