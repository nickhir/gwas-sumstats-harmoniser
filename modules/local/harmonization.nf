process harmonization {
    tag "${GCST}_${chrom}"

    conda (params.enable_conda ? "${task.ext.conda}" : null)

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
    tuple val(GCST), val(palin_mode), val(status), val(chrom), path(merged), path(yaml), path(ref)


    output:
    tuple val(GCST), val(palin_mode), path("${chrom}.merged.hm"),path("${chrom}.merged.log.tsv.gz"), emit: hm_by_chrom

    when:
    status=="contiune"

    shell:
    """

    coordinate_system=\$(grep coordinate_system $yaml | awk -F ":" '{print \$2}' | tr -d "[:blank:]" )
    if test -z "\$coordinate_system"; then coordinate="1-based"; else coordinate=\$coordinate_system; fi

    header_args=\$(utils.py -f $merged -harm_args);

    line_count=\$(wc -l < $merged)

    if test "\$line_count" -gt 1; then
        main_pysam.py \
        --sumstats $merged \
        --vcf ${params.ref}/homo_sapiens-${chrom}.vcf.gz \
        ${params.new_ref ? "--new_vcf ${params.new_ref}/homo_sapiens-${chrom}.vcf.gz" : ""} \
        --hm_sumstats ${chrom}.merged_unsorted.hm \
        --hm_statfile ${chrom}.merged.log.tsv.gz \
        \$header_args \
        --na_rep_in NA \
        --na_rep_out NA \
        --coordinate \$coordinate \
        --palin_mode $palin_mode;

        chr=\$(awk -v RS='\t' '/chromosome/{print NR; exit}' ${chrom}.merged_unsorted.hm)
        pos=\$(awk -v RS='\t' '/base_pair_location/{print NR; exit}' ${chrom}.merged_unsorted.hm)

        head -n1 ${chrom}.merged_unsorted.hm > ${chrom}.merged.hm;
        tail -n+2 ${chrom}.merged_unsorted.hm | sort -n -k\$chr -k\$pos -T\$PWD >> ${chrom}.merged.hm
    else
        head -n1 $merged > ${chrom}.merged.hm;
        echo "hm_code\tcount\tdescription" > ${chrom}.merged.log.tsv
        gzip ${chrom}.merged.log.tsv
    fi
    """
}
