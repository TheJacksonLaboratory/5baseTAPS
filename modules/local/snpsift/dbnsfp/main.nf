process SNPSIFT_DBNSFP {
    tag "${meta.id}:${callset}:${chrom}"
    label 'process_low'

    conda "bioconda::snpsift=5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'
        : 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d' }"

    // SnpSift dbnsfp strips 'chr' from VCF contigs before lookup, so the database MUST be a
    // bare-chromosome-name file (...nochr.txt.gz); with a chr-prefixed file it annotates 0
    // variants and the task still succeeds. The run is checked below for exactly that.

    input:
    tuple val(meta), val(callset), val(chrom), path(vcf)
    val dbnsfp

    output:
    tuple val(meta), val(callset), val(chrom), path("*.annotated.vcf"), emit: vcf
    path "versions.yml",                                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}.${callset}.${chrom}"
    def avail_mem = task.memory ? (task.memory.toGiga() - 1) : 6
    def fields  = 'SIFT_score,SIFT_pred,Polyphen2_HDIV_pred,REVEL_score,CADD_phred,aapos,Ensembl_transcriptid,genename'
    """
    java -Xmx${avail_mem}G -jar /opt/snpEff/SnpSift.jar \\
        dbnsfp -v -db ${dbnsfp} -noDownload -a \\
        -f ${fields} \\
        ${vcf} > ${prefix}.annotated.vcf 2> dbnsfp.err

    # fail loudly if dbNSFP matched nothing although SnpEff called many missense variants
    # (dbNSFP covers missense SNVs, so 0 annotated entries there means a database/naming problem)
    N_MIS=\$(grep -v '^#' ${vcf} | grep -c 'missense_variant' || true)
    N_ANN=\$(grep -o 'Total annotated entries : [0-9]*' dbnsfp.err | grep -o '[0-9]*\$' || echo 0)
    echo "missense_variants=\$N_MIS dbnsfp_annotated=\$N_ANN"
    if [ "\$N_MIS" -ge 20 ] && [ "\${N_ANN:-0}" -eq 0 ]; then
        echo "ERROR: dbNSFP annotated 0 entries although \$N_MIS missense variants are present (chr-prefixed database file?)" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: 5.1d
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}.${chrom}"
    """
    touch ${prefix}.annotated.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: stub
    END_VERSIONS
    """
}
