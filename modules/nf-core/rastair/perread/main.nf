process PERREAD {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker://sbludwig/rastair:version-2.1.1'

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.rastair_perread.bed.gz"), emit: bed
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rastair per-read \\
        --threads ${task.cpus} \\
        --fasta-file ${fasta} \\
        --bed ${prefix}.rastair_perread.bed.gz \\
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
    touch ${prefix}.rastair_perread.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version 2>&1 || echo "stub")
    END_VERSIONS
    """
}
