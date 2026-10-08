process SNPEFF_ANNOTATE {
    tag "${meta.id}:${callset}:${chrom}"
    label 'process_low'

    conda "bioconda::snpeff=5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'
        : 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d' }"

    // The SnpEff config is passed as a path string, not staged: its relative data.dir (./data)
    // is resolved against the config file's own directory.

    input:
    tuple val(meta), val(callset), val(chrom), path(vcf)
    val snpeff_db
    val snpeff_config

    output:
    tuple val(meta), val(callset), val(chrom), path("*.snpeff.vcf"), emit: vcf
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}.${callset}.${chrom}"
    def avail_mem = task.memory ? (task.memory.toGiga() - 1) : 6
    """
    java -Djava.io.tmpdir=./ -Xmx${avail_mem}G -jar /opt/snpEff/snpEff.jar \\
        ${snpeff_db} \\
        -c ${snpeff_config} \\
        -noStats \\
        ${vcf} > ${prefix}.snpeff.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: 5.1d
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}.${chrom}"
    """
    touch ${prefix}.snpeff.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: stub
    END_VERSIONS
    """
}
