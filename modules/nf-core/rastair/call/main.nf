process RASTAIR_CALL {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker://sbludwig/rastair:version-2.1.1'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.rastair_call.bed.gz"), emit: bed
    tuple val(meta), path("*.rastair_call.vcf.gz"), emit: vcf
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_arg = "--vcf ${prefix}.rastair_call.vcf.gz"
    def threads = task.ext.threads ?: task.cpus

    """
    rastair call \\
        --threads ${threads} \\
        --fasta-file ${fasta} \\
        --bed ${prefix}.rastair_call.bed.gz \\
        ${vcf_arg} \\
        ${args} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_call.bed.gz
    touch ${prefix}.rastair_call.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version 2>&1 || echo "stub")
    END_VERSIONS
    """
}
