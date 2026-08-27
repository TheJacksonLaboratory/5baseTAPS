process RASTAIR_MBIAS {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'docker://sbludwig/rastair:version-2.1.1'

    input:
    tuple val(meta), path(bed), path(tbi)

    output:
    tuple val(meta), path("*.rastair_mbias.html"), emit: html
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_mbias_out
    rastair mbias ${args} \\
        --output-prefix ${prefix}_mbias_out \\
        --bed ${bed}
    mv ${prefix}_mbias_out/qc_report.html ${prefix}.rastair_mbias.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rastair_mbias.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version 2>&1 || echo "stub")
    END_VERSIONS
    """
}
