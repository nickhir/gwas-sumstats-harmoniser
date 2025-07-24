process map_to_build {
    // this tag, i.e. GCST is simply the first value of the tuple
    // in our case this often corresponds to the basename
    tag "$GCST"

    container "${task.ext.singularity}${task.ext.singularity_version}"

        
    input:
        // file contains 3 values that are "unpacked" to the tuple
        tuple val(GCST), path(yaml), path(tsv)
        val chr

    output:
        tuple val(GCST), path ('*.merged'), path('unmapped'), path(yaml), emit:mapped
    
    shell:
    """
    coordinate_system=\$(grep coordinate_system $yaml | awk -F ":" '{print \$2}' | tr -d "[:blank:]" )
    if test -z "\$coordinate_system"; then coordinate="1-based"; else coordinate=\$coordinate_system; fi
    from_build=\$((grep genome_assembly $yaml | grep -Eo '[0-9][0-9]') || (echo \$(basename $tsv) | grep -Eo '[bB][a-zA-Z]*[0-9][0-9]' | grep -Eo '[0-9][0-9]'))
    [[ -z "\$from_build" ]] && { echo "Parameter from_build is empty" ; exit 1; }
    
    # Check if liftover_only parameter is set to use liftover-only script
    if [[ "${params.liftover_only}" == "true" ]]; then
        map_to_build_nf_liftover.py \
        -f $tsv \
        -from_build \$from_build \
        -to_build $params.to_build \
        -chroms "${chr}" \
        -coordinate \$coordinate
    else
        map_to_build_nf.py \
        -f $tsv \
        -from_build \$from_build \
        -to_build $params.to_build \
        -vcf "${params.ref}/homo_sapiens-chr*.parquet" \
        -chroms "${chr}" \
        -coordinate \$coordinate
    fi
    """
}
