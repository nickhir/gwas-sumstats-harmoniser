process subsample_variants {
    tag "$GCST"

    container "${task.ext.singularity}${task.ext.singularity_version}"

    input:
        tuple val(GCST), path(yaml), path(tsv)

    output:
        tuple val(GCST), path(yaml), path("${GCST}.subsampled.tsv"), path("${GCST}.header_alias.log"), emit: subsampled

    shell:
    """
    subsample_variants_nf.py \
        -f $tsv \
        -o ${GCST}.subsampled.tsv \
        -l ${params.subsample_limit ?: 10000000} \
        --alias_log ${GCST}.header_alias.log
    """
}
