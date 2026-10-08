process BCFTOOLS_ANNOTATERESOURCES {
    tag "${meta.id}:${callset}:${chrom}"
    label 'process_low'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.17--h3cc50cf_1'
        : 'biocontainers/bcftools:1.17--h3cc50cf_1' }"

    // One chromosome of one VCF. Adds dbSNP rsID (ID column), gnomAD AF (INFO/gnomAD_AF)
    // and ClinVar (INFO/CLN*, GENEINFO, ALLELEID). Resources are passed as path strings (not
    // staged): they are huge, tabix-indexed, and looked up by position. Lookups are restricted
    // to the variant positions (-R sites --regions-overlap 0); without this bcftools streams the
    // whole resource (dbSNP: >9 min stalled vs 70 s on chr1 of the test sample).

    input:
    tuple val(meta), val(callset), path(vcf), path(tbi), val(chrom)
    val dbsnp
    val gnomad
    val clinvar

    output:
    tuple val(meta), val(callset), val(chrom), path("*.resources.vcf"), emit: vcf
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}.${chrom}"
    """
    bcftools view -r ${chrom} ${vcf} -O z -o chunk.vcf.gz
    bcftools index --tbi chunk.vcf.gz

    # variant positions; an empty chromosome gets one dummy position so the header still
    # receives the new INFO definitions and chunks concatenate cleanly
    bcftools query -f '%CHROM\\t%POS\\n' chunk.vcf.gz | uniq > sites.tsv
    [ -s sites.tsv ] || printf '${chrom}\\t1\\n' > sites.tsv

    bcftools annotate -a ${dbsnp} -c ID \\
        -R sites.tsv --regions-overlap 0 chunk.vcf.gz -O z -o s1.vcf.gz
    bcftools index --tbi s1.vcf.gz

    bcftools annotate -a ${gnomad} -c INFO/gnomAD_AF:=INFO/AF \\
        -R sites.tsv --regions-overlap 0 s1.vcf.gz -O z -o s2.vcf.gz
    bcftools index --tbi s2.vcf.gz

    bcftools annotate -a ${clinvar} \\
        -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNDISDB,INFO/GENEINFO,INFO/ALLELEID \\
        -R sites.tsv --regions-overlap 0 s2.vcf.gz -O v -o ${prefix}.resources.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${callset}.${chrom}"
    """
    touch ${prefix}.resources.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //' || echo "stub")
    END_VERSIONS
    """
}
