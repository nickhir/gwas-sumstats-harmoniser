process filter_significant_variants {
    tag "${GCST}_${chrom}"

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
        tuple val(chrom), val(GCST), path(merged), path(yaml)

    output:
        tuple val(chrom), val(GCST), path("${chrom}.significant.merged"), path(yaml), emit: filtered

    shell:
    """
    filter_significant_nf.py \
        -f $merged \
        -o ${chrom}.significant.merged
    """
}