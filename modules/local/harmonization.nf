process harmonization {

    tag "${GCST}_${chrom}"

    /* -------- env / container ------------------------------------------ */
    conda( params.enable_conda ? task.ext.conda : null )

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container
            ? "${task.ext.singularity}${task.ext.singularity_version}"
            : "${task.ext.docker}${task.ext.docker_version}" }"

    /* -------------------------- I/O ------------------------------------- */
    input:
    tuple val(GCST), val(palin_mode), val(status), val(chrom),
          path(merged), path(yaml), path(ref)

    output:
    tuple val(GCST), val(palin_mode),
          path("${chrom}.merged.hm"),
          path("${chrom}.merged.log.tsv.gz"),
          emit: hm_by_chrom

    when:
    status == 'contiune'          // <- verify spelling

    /* ------------------------ script ------------------------------------ */
    script:
    def newVcfOpt = params.new_ref
        ? "--new_vcf ${params.new_ref}/homo_sapiens-${chrom}.vcf.gz"
        : ''

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # 1) Coordinate system
    coordinate_system=\$(grep coordinate_system ${yaml} \\
                          | awk -F ":" '{print \$2}' \\
                          | tr -d "[:blank:]")
    if [[ -z "\$coordinate_system" ]]; then
        coordinate="1-based"
    else
        coordinate=\$coordinate_system
    fi

    # 2) Header args & line count
    header_args=\$(utils.py -f ${merged} -harm_args)
    line_count=\$(wc -l < ${merged})

    if [[ "\$line_count" -gt 1 ]]; then
        # ---- main harmonisation -----------------------------------------
        main_pysam.py \\
            --sumstats ${merged} \\
            --vcf ${params.ref}/homo_sapiens-${chrom}.vcf.gz \\
            ${newVcfOpt} \\
            --hm_sumstats ${chrom}.merged_unsorted.hm \\
            --hm_statfile ${chrom}.merged.log.tsv.gz \\
            \$header_args \\
            --na_rep_in  NA \\
            --na_rep_out NA \\
            --coordinate \$coordinate \\
            --palin_mode ${palin_mode}

        # Sort by chromosome / position
        chr=\$(awk -v RS='\\t' '/chromosome/{print NR; exit}' \\
                   ${chrom}.merged_unsorted.hm)
        pos=\$(awk -v RS='\\t' '/base_pair_location/{print NR; exit}' \\
                   ${chrom}.merged_unsorted.hm)

        head -n1 ${chrom}.merged_unsorted.hm > ${chrom}.merged.hm
        tail -n+2 ${chrom}.merged_unsorted.hm \\
            | sort -n -k\$chr -k\$pos -T\$PWD >> ${chrom}.merged.hm

    else
        # -------- ≤ 1 data line: ensure both neg_log_10_p_value & p_value --
        header=\$(head -n1 ${merged})
        new_header="\$header"

        # track which columns we're adding
        add_neg_log_p=false
        add_p_value=false

        # add missing columns (grep -w = exact word match)
        if ! echo "\$new_header" | grep -qw "neg_log_10_p_value"; then
            new_header="\${new_header}\\tneg_log_10_p_value"
            add_neg_log_p=true
        fi
        if ! echo "\$new_header" | grep -qw "p_value"; then
            new_header="\${new_header}\\tp_value"
            add_p_value=true
        fi

        # write header
        echo -e "\$new_header" > ${chrom}.merged.hm

        # process data rows if they exist
        if [[ "\$line_count" -eq 2 ]]; then
            # exactly one data row - process it
            data_row=\$(tail -n1 ${merged})
            new_data_row="\$data_row"
            
            # add missing column values
            if [[ "\$add_neg_log_p" == "true" ]]; then
                new_data_row="\${new_data_row}\\tNA"
            fi
            if [[ "\$add_p_value" == "true" ]]; then
                new_data_row="\${new_data_row}\\tNA"
            fi
            
            echo -e "\$new_data_row" >> ${chrom}.merged.hm
        fi

        echo "hm_code\\tcount\\tdescription" > ${chrom}.merged.log.tsv
        gzip ${chrom}.merged.log.tsv
    fi
    """
}

