include {map_to_build} from '../../modules/local/map_to_build'
include {subsample_variants} from '../../modules/local/subsample_variants'
include {filter_significant_variants} from '../../modules/local/filter_significant'
include {ten_percent_counts} from '../../modules/local/ten_percent_counts'
include {ten_percent_counts_sum} from '../../modules/local/ten_percent_counts_sum'
include {generate_strand_counts} from '../../modules/local/generate_strand_counts'
include {summarise_strand_counts} from '../../modules/local/summarise_strand_counts'

workflow major_direction{
    take:
        // list of chromsomes (chr1,chr2,...). 
        // files is actually a touple containing [basename, path of yaml, path of tsv]
        chr
        files
    
    main:
        // chroms is just a list of 1,2,3,4...
        chroms=chr.flatten().map{it.toString().replaceAll("chr","")}.collect()
        subsample_variants(files)
        alias_ch=subsample_variants.out.subsampled.map{tuple(it[0],it[3])}
        map_to_build(subsample_variants.out.subsampled.map{it[0..2]}, chroms)
        //example: output is [GCST1,[path of 1.merged, path of 2.merged .....]]
        // map_to_build.out.mapped looks like this
        // [GCST_ID, [path to 1.merged, path to 2.merged] [path to unmerged] [path to yaml]]
        // transpose trasforms it into something like this
        // [GCST_ID, path to 1.merged, path to unmerged,path to yaml]
        // [GCST_ID, path to 2.merged, path to unmerged,path to yaml]
        // [GCST_ID, path to 3.merged, path to unmerged,path to yaml]
        // ...
        // then we just extract information and store it in the channel "map_chr_ch"
        map_to_build.out.mapped
                        .transpose()
                        .map{tuple(get_chr(it[1]),it[0],it[1],it[3])}
                        .set{map_chr_ch}
        
        // run the filter_significant_variants module
        filter_significant_variants(map_chr_ch)
        // assign the output of filter_significant_variants to a channel called "sig_chr_ch". This contains a tuple with the output of filter_significant_variants
        // tuple val(chrom), val(GCST), path("${chrom}.significant.merged"), path(yaml)
        filter_significant_variants.out.filtered.set{sig_chr_ch}

        // capture unmapped sites for reporting
        unmapped = map_to_build.out.mapped.map{tuple(it[0],it[2])}
        
        // joint is still needed in case not all chr are running
        /* example:
        homo_sapiens-chr1.vcf.gz ->[chr1, path of homo_sapiens-chr1.vcf.gz] (ref_ch_chr)
        */       
        Channel.fromPath("${params.ref}/homo_sapiens-chr*.vcf.gz") 
            .map { prepare_reference (it) }
            .set{ ref_chr_ch }
        
        
        /* example
        [chr1, path of homo_sapiens-chr1.vcf.gz] (ref_chr_ch) + 
        [chr1, GCST1, path of 1.merged] (map_chr_ch)
        -> [chr1, GCST1, path of 1.merged,path of homo_sapiens-chr1.vcf.gz] (count_ch) 
        */
        // combine is essentially a join in this case based on the first element (chromosome name)
        count_ch=map_chr_ch.combine(ref_chr_ch,by:0)

        // perform merge. adds the path to the reference to the output of filter_significant_variants
        filtered_count_ch=sig_chr_ch.combine(ref_chr_ch,by:0)

        
        // we still perform the ten_percent_counts on the whole subset
        ten_percent_counts(count_ch)
        
        
        // get the number of chromosomes
        int nchr=params.chrom.size()


        ten_to_sum = ten_percent_counts.out
                        .ten_sc
                        .groupTuple(by: 0, size: nchr) // combine the outputs by GCST ID -> [GCST_ID, [[list to .sc files per chrom]]]
                        .branch{pass:it[1].size()==nchr} //if the number of .sc is not equal to the total number of processed chromosomes we stop. otherwise continue
                        .map{it[0]}

        // example: ten_to_sum [GCST1],[GCST2].....
        ten_percent_counts_sum(ten_to_sum)
        //example: output [GCST,path ten_percent.tsv,drop,rerun],[GCST,path ten_percent.tsv,forward,countiune]
        
        // determine whether conatin a string in the output txt file
        ten_percent_counts_sum.out.ten_sum.branch{rerun:it.contains("rerun")
                                            contiune:it.contains("contiune")}
                                        .set{branch}
        
        //branch rerun
        /* example: 
        [GCST, path ten_percent.tsv, drop, rerun] (rerun_branch)
        [chr1, GCST1, path of 1.merged,path of homo_sapiens-chr1.vcf.gz] (count_ch) 
        */
        branch.rerun.map{tuple(it[3],it[0])}.set{all_sc_ch}
        //example:rerun_branch into: [rerun,GCST1]
        count_ch.combine(all_sc_ch,by:1).set{rerun_branch}
        //example: rerun_branch: [GCST1,chr1, path of 1.merged,path of homo_sapiens-chr1.vcf.gz,return]
        generate_strand_counts(rerun_branch)
        //example: generate_strand_counts.out.all_sc: [GCST, rerun,path of full_sc]
        all_to_sum=generate_strand_counts.out.all_sc.collect().map{tuple(it[0],it[1])}.unique()
        //all_to_sum: [GCST, rerun]
        summarise_strand_counts(all_to_sum)
        // example: [GCST,path Full.tsv,drop,contiune],[GCST,path Full.tsv,forward,contiune]

        //branch contiune
        // [GCST, ten_percent, forward,contiune] (contiune_branch)
        all_files=summarise_strand_counts.out.all_sum.mix(branch.contiune)
        //hm_input: [GCST,path ten_percent.tsv,forward,countiune],[GCST,path Full.tsv,reverse,countiune]
        rearrnaged_count_ch=filtered_count_ch.map{tuple(it[1],it[0],it[2],it[3],it[4])}
        // example: [chr1, GCST1, path of 1.merged,path of homo_sapiens-chr1.vcf.gz] (count_ch) 
        // example into: [GCST1,chr1,path of merged, path of vcf]
        all_input=all_files.combine(rearrnaged_count_ch,by:0)
        //example: [GCST,path ten_percent.tsv,forward,countiune,chr,path of merged, path of vcf]
        hm_input=all_input.map{it[0,2..7]}
        direction_sum=all_input.map{it[0..1]}.unique()

    emit:
        hm_input=hm_input
        direction_sum=direction_sum
        unmapped=unmapped
        alias_log=alias_ch
}

// groovy function
def prepare_reference (Path input) {
    // extract chromosome from file path and form a list in list
    return [input.getName().split('-')[1].split('\\.')[0], input]
}

def get_chr(Path input) {
    return ("chr"+input.getName().split('\\.')[0])
}
