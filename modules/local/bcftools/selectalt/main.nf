process BCFTOOLS_SELECTALT {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.17--h3cc50cf_1'
        : 'biocontainers/bcftools:1.17--h3cc50cf_1' }"

    // The rastair call VCF is ~99% CpG methylation rows with ALT='.' (no variant).
    // Keep only records with a real ALT allele and index them for per-chromosome scatter.

    input:
    tuple val(meta), val(callset), path(vcf)

    output:
    tuple val(meta), val(callset), path("*.alt_only.vcf.gz"), path("*.alt_only.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}"
    """
    bcftools view -e 'ALT="."' ${vcf} -O z -o ${prefix}.alt_only.vcf.gz
    bcftools index --tbi ${prefix}.alt_only.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}"
    """
    touch ${prefix}.alt_only.vcf.gz ${prefix}.alt_only.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //' || echo "stub")
    END_VERSIONS
    """
}
