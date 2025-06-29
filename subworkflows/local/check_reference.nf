workflow chr_check {
    
    // check what is available in reference folder
    vcf_files=Channel.fromPath("${params.ref}/*.vcf.gz").map{extract_name(it)}
    tbi_files=Channel.fromPath("${params.ref}/*.tbi").map{extract_name(it)}
    parquet_files=Channel.fromPath("${params.ref}/*.parquet").map{extract_name(it)}

    // check what is required to be harmonized
    // ch_chrom looks like: [chr1,chr2,chr3...]
    ch_chrom=Channel.from(params.chrom).map{"chr"+it}

    // this matches items that are shared between each set
    ch=vcf_files.join(tbi_files).join(parquet_files).join(ch_chrom)

    // this is just for logging
    ch.toList().subscribe { println it + ' is being harmonized' }

    // this is just a list of the chromosomes for which we have all the info.
    emit:
    ch_input=ch.collect()
}
// this extracts the name of the chromosome. 
def extract_name(Path input){
     return [input.getName().split('-')[1].split('\\.')[0]]
}
