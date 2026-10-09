/*
 * Per-chromosome M-bias subworkflow (two-stage)
 *
 * Replaces the monolithic RASTAIR_MBIAS process for samples where total CpG
 * observations exceed R's 2^31 vector-length limit.
 *
 * Stage 1  GET_CHRS → MBIAS_CHR (per-chromosome computation)
 * Stage 2  RENDER (loads per-chromosome results, renders HTML)
 *
 * Output: ${meta.id}.rastair_mbias.html
 */

process GET_CHRS {
    tag "${meta.id}"

    container 'docker://sbludwig/rastair:version-2.1.1'

    input:
    tuple val(meta), path(bgz), path(tbi)

    output:
    tuple val(meta), path(bgz), path(tbi), path("chromosomes.txt")

    script:
    """
    tabix -l ${bgz} \\
        | grep -E '^(chr)?([0-9]+|[XY])\$' \\
        > chromosomes.txt
    """

    stub:
    """
    printf 'chr1\nchr2\nchr3\n' > chromosomes.txt
    """
}

process MBIAS_CHR {
    tag "${meta.id}_${chr}"

    container 'docker://sbludwig/rastair:version-2.1.1'

    input:
    tuple val(meta), path(bgz), path(tbi), val(chr)

    output:
    tuple val(meta), path("${meta.id}_${chr}_mbias_table.rds")

    script:
    def vbias_flag = params.plot_vbias ? "--vbias" : ""
    def assets     = params.rastair_rscript_dir ?: "${workflow.projectDir}/bin/rastair_scripts"
    """
    Rscript ${assets}/perChr_mbias_chr.R \\
        --bed          ${bgz} \\
        --chr          ${chr} \\
        --output       ${meta.id}_${chr}_mbias_table.rds \\
        --tabix-path   tabix \\
        --include-flag 3 \\
        --exclude-flag 3852 \\
        --threads      ${task.cpus} \\
        ${vbias_flag}
    """

    stub:
    """
    touch ${meta.id}_${chr}_mbias_table.rds
    """
}

process RENDER {
    tag "${meta.id}"

    container 'docker://sbludwig/rastair:version-2.1.1'

    input:
    tuple val(meta), path(rds_files), path(vcf), val(fasta_path)

    output:
    tuple val(meta), path("${meta.id}.rastair_mbias.html"), emit: html

    script:
    def prefix     = meta.id
    def assets     = params.rastair_rscript_dir ?: "${workflow.projectDir}/bin/rastair_scripts"
    def vbias_flag = params.plot_vbias ? "--vbias" : ""
    def gc_flags   = params.plot_gc
        ? "--gc --vcf \$(realpath ${vcf}) --reference ${fasta_path} --bcftools bcftools"
        : ""
    """
    mkdir -p rds_in out
    mv ${rds_files} rds_in/
    Rscript ${assets}/render_allChr_mbias.R \\
        --rds-dir       rds_in \\
        --output-prefix out \\
        --sample-name   ${meta.id} \\
        --threads       ${task.cpus} \\
        ${vbias_flag} \\
        ${gc_flags}
    mv out/qc_report.html ${prefix}.rastair_mbias.html
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}.rastair_mbias.html
    """
}

workflow PERCHROMOSOME_MBIAS {

    take:
    ch_bed   // channel: [ val(meta), path(bed.gz), path(tbi) ]
    ch_vcf   // channel: [ val(meta), path(vcf.gz) ] — always staged; used only when params.plot_gc
    ch_fasta // value channel: [ val(meta), path(fa) ] — reference genome for GC bias

    main:
    // Get primary chromosome list from the tabix index
    GET_CHRS(ch_bed)

    ch_per_chr = GET_CHRS.out
        .flatMap { meta, bgz, tbi, chr_file ->
            def chrs = chr_file.readLines().findAll { it.trim() }
            def n    = chrs.size()
            chrs.collect { chr -> tuple(groupKey(meta, n), bgz, tbi, chr.trim()) }
        }

    // Stage 1: parallel per-chromosome RDS computation
    MBIAS_CHR(ch_per_chr)

    // Collect all per-chr RDS files per sample, then join VCF for GC bias.
    ch_gathered = MBIAS_CHR.out
        .groupTuple()
        .map { gkey, rds_list -> tuple(gkey.target, rds_list.flatten()) }

    ch_fasta_path = ch_fasta.map { meta, fa -> fa.toAbsolutePath().toString() }.first()

    ch_render = ch_gathered
        .join(ch_vcf)
        .combine(ch_fasta_path)
        .map { meta, rds_list, vcf, fasta_path -> tuple(meta, rds_list, vcf, fasta_path) }

    // Stage 2: gather RDS files and render the HTML report
    RENDER(ch_render)

    emit:
    html = RENDER.out.html  // channel: [ val(meta), path("*.rastair_mbias.html") ]
}
