// if reference files are not exist, download and prepare reference
include {get_vcf_files} from '../../modules/local/get_vcf_files'

// no need to emit anything, this is the last workflow in the chain
workflow prepare_reference {
    take: in_chrom
    main: 
        get_vcf_files(in_chrom)
}