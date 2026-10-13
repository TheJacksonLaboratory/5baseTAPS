process BCFTOOLS_CONCATANNOTATED {
    tag "${meta.id}:${callset}"
    label 'process_low'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.17--h3cc50cf_1'
        : 'biocontainers/bcftools:1.17--h3cc50cf_1' }"

    input:
    tuple val(meta), val(callset), path(vcfs)

    output:
    tuple val(meta), val(callset), path("*.annotated.vcf.gz"), path("*.annotated.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}"
    """
    # chunk files are named <id>.<callset>.<chrom>.annotated.vcf ; version sort gives chr1..chr22,chrX,chrY
    ls *.annotated.vcf | sort -V > chunks.list
    bcftools concat --file-list chunks.list -O z -o ${prefix}.annotated.vcf.gz
    bcftools index --tbi ${prefix}.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}"
    """
    touch ${prefix}.annotated.vcf.gz ${prefix}.annotated.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //' || echo "stub")
    END_VERSIONS
    """
}
