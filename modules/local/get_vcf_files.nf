/* download reference */
process get_vcf_files {
    // custom label for the process
    // makes logging and reporting easier
    tag "${chr}"

    // container is a nextflow directive
    // specifies the container to use for this process
    container "${task.ext.singularity}${task.ext.singularity_version}"

    // this is essentially the output directory for the process
    storeDir params.ref

    input:
        val chr
    
    output:
        // indicates output is a tuple. in this case a value and a path
        // value does not have to be a numeric. just means single value.
        tuple val(chr), path("homo_sapiens-${chr}.parquet")
    
    shell:
    """
    # Check if the directory exists; if not, create it
    [[ -d $params.ref ]] || mkdir -p $params.ref

    # Check if the VCF file already exists; if not, download it
    if [[ ! -f $params.ref/homo_sapiens-${chr}.vcf.gz ]]; then
        wget -P $params.ref ${params.remote_vcf_location}/homo_sapiens-${chr}.vcf.gz
    fi

    # Check if the index file exists; if not, create it
    if [[ ! -f $params.ref/homo_sapiens-${chr}.vcf.gz.tbi ]]; then
        tabix -f -p vcf $params.ref/homo_sapiens-${chr}.vcf.gz
    fi

    # notice how -o does not specify the storeDir again. 
    # files will automatically be placed there.
    vcf2parquet_nf.py \
        -f ${params.ref}/homo_sapiens-${chr}.vcf.gz \
        -o homo_sapiens-${chr}.parquet
    """
}
