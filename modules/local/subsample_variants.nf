process subsample_variants {
    tag "$GCST"

    conda (params.enable_conda ? "${task.ext.conda}" : null)

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
        tuple val(GCST), path(yaml), path(tsv)

    output:
        tuple val(GCST), path(yaml), path("${GCST}.subsampled.tsv"), emit: subsampled

    shell:
    """
    subsample_variants_nf.py \
        -f $tsv \
        -o ${GCST}.subsampled.tsv \
        -l ${params.subsample_limit ?: 10000000}
    """
}
