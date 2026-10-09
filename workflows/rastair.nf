/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_TAPS_CONVERSION as METHYLSEQ } from '../subworkflows/nf-core/bam_taps_conversion/main'
include { CPG_MASK              } from '../modules/local/rastair_cpg_mask/main'
include { LAMBDA_METH     } from '../subworkflows/local/ctrl_genome_methylation/main'
include { PUC19_METH      } from '../subworkflows/local/ctrl_genome_methylation/main'
include { METHYL_QC           } from '../modules/local/taps_methyl_metrics/main'
include { ALIGN_QC          } from '../modules/local/taps_mapping_metrics/main'
include { WGS_QC     } from '../modules/local/taps_wgs_coverage_metrics/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DOWNSTREAM ANALYSIS WORKFLOW
    Receives the final consensus-filtered BAM from FASTQUORUM and runs:
      - METHYLSEQ  : rastair methylation calling (TAPS)
      - GATK HaplotypeCaller : variant calling            (TODO)
      - TVC (Watchmaker)     : variant calling            (TODO)
      - MultiQC              : aggregate QC report        (TODO - run last)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RASTAIR {

    take:
    ch_bam                // channel: [ val(meta), path(bam) ]
    ch_bai                // channel: [ val(meta), path(bai) ]
    ch_fasta              // channel: [ val(meta), path(fasta) ]    (value channel)
    ch_fasta_fai          // channel: [ val(meta), path(fai) ]      (value channel)
    ch_unmapped_bam       // channel: [ val(meta), path(unmapped.bam) ] — one per samplesheet row; used for control genome QC
    ch_prededup_flagstat  // channel: [ val(meta), path(*.flagstat) ]      — from SAMTOOLS_FLAGSTAT (raw BAM, pre-dedup)
    ch_mosdepth_summary   // channel: [ val(meta), path(*.mosdepth.summary.txt) ]
    ch_mosdepth_dist      // channel: [ val(meta), path(*.mosdepth.global.dist.txt) ]

    main:

    ch_versions           = Channel.empty()
    ch_multiqc_files      = Channel.empty()
    ch_lambda_summary     = Channel.empty()
    ch_puc19_summary      = Channel.empty()
    ch_lambda_flagstat    = Channel.empty()
    ch_puc19_flagstat     = Channel.empty()
    ch_mapping_metrics_csv = Channel.empty()
    ch_methyl_metrics_csv  = Channel.empty()
    ch_wgs_metrics_csv     = Channel.empty()

    //
    // Build broadcast inputs: combine per-sample bam/bai with reference channels
    //
    ch_taps_inputs = ch_bam
        .join(ch_bai)
        .combine(ch_fasta)
        .combine(ch_fasta_fai)
        .multiMap { meta, bam, bai, _meta_fasta, fasta, _meta_fai, fai ->
            bam:         [ meta, bam ]
            bai:         [ meta, bai ]
            fasta:       [ meta, fasta ]
            fasta_index: [ meta, fai ]
        }

    //
    // SUBWORKFLOW: TAPS methylation calling with Rastair 2.1.1
    //   rastair call      -> per-position BED (methylation) + VCF (variants)
    //   rastair mbias     -> M-bias HTML report (from BED)
    //   rastair methylkit -> methylKit-format table (from BED, shell AWK)
    //
    METHYLSEQ(
        ch_taps_inputs.bam,
        ch_taps_inputs.bai,
        ch_taps_inputs.fasta,
        ch_taps_inputs.fasta_index
    )
    ch_versions = ch_versions.mix(METHYLSEQ.out.versions)

    //
    // MODULE: CpG mask for GATK HaplotypeCaller
    //   All CpG positions called by rastair → bgzipped + tabix-indexed BED
    //   Used as --exclude-intervals to suppress TAPS 5mC→T conversion artifacts
    //
    CPG_MASK(METHYLSEQ.out.bed)
    ch_versions = ch_versions.mix(CPG_MASK.out.versions)

    if (params.lambda_fasta) {
        lambda_fasta   = Channel.fromPath(params.lambda_fasta).map { it -> [[id: it.baseName], it] }.collect()
        lambda_fai     = Channel.fromPath(params.lambda_fai).map { it -> [[id: 'fai'], it] }.collect()
        lambda_dict    = Channel.fromPath(params.lambda_dict).map { it -> [[id: 'dict'], it] }.collect()
        lambda_bwamem2 = Channel.fromPath(params.lambda_bwamem2, type: 'dir').map { it -> [[id: 'bwamem2'], it] }.collect()
        LAMBDA_METH(ch_unmapped_bam, lambda_fasta, lambda_fai, lambda_dict, lambda_bwamem2)
        ch_versions        = ch_versions.mix(LAMBDA_METH.out.versions)
        ch_lambda_summary  = LAMBDA_METH.out.meth_summary
        ch_lambda_flagstat = LAMBDA_METH.out.flagstat
        ch_multiqc_files   = ch_multiqc_files.mix(ch_lambda_flagstat.map { it[1] }.collect())
    }
    if (params.puc19_fasta) {
        puc19_fasta   = Channel.fromPath(params.puc19_fasta).map { it -> [[id: it.baseName], it] }.collect()
        puc19_fai     = Channel.fromPath(params.puc19_fai).map { it -> [[id: 'fai'], it] }.collect()
        puc19_dict    = Channel.fromPath(params.puc19_dict).map { it -> [[id: 'dict'], it] }.collect()
        puc19_bwamem2 = Channel.fromPath(params.puc19_bwamem2, type: 'dir').map { it -> [[id: 'bwamem2'], it] }.collect()
        PUC19_METH(ch_unmapped_bam, puc19_fasta, puc19_fai, puc19_dict, puc19_bwamem2)
        ch_versions       = ch_versions.mix(PUC19_METH.out.versions)
        ch_puc19_summary  = PUC19_METH.out.meth_summary
        ch_puc19_flagstat = PUC19_METH.out.flagstat
        ch_multiqc_files  = ch_multiqc_files.mix(ch_puc19_flagstat.map { it[1] }.collect())
    }

    //
    // MODULE: TAPS methylation metrics — combine main + control summaries per sample
    //   Only runs when both lambda and pUC19 control summaries are available.
    //
    if (params.lambda_fasta && params.puc19_fasta) {
        ch_methyl_metrics_in = METHYLSEQ.out.meth_summary
            .join(ch_lambda_summary, by: 0)
            .join(ch_puc19_summary,  by: 0)

        METHYL_QC(ch_methyl_metrics_in)
        ch_versions           = ch_versions.mix(METHYL_QC.out.versions)
        ch_multiqc_files      = ch_multiqc_files.mix(METHYL_QC.out.mqc.map { it[1] }.flatten())
        ch_methyl_metrics_csv = METHYL_QC.out.mqc
    }

    //
    // MODULE: ALIGN_QC — mapping + spike-in composition per sample
    //   Only runs when both lambda and pUC19 controls are available (provides all 3 flagstats).
    //
    if (params.lambda_fasta && params.puc19_fasta) {
        ch_mapping_metrics_in = ch_prededup_flagstat
            .join(ch_lambda_flagstat, by: 0)
            .join(ch_puc19_flagstat,  by: 0)

        ALIGN_QC(ch_mapping_metrics_in)
        ch_versions            = ch_versions.mix(ALIGN_QC.out.versions)
        ch_multiqc_files       = ch_multiqc_files.mix(ALIGN_QC.out.mqc.map { it[1] }.flatten())
        ch_mapping_metrics_csv = ALIGN_QC.out.mqc
    }

    //
    // MODULE: WGS_QC — WGS depth + CpG coverage per sample
    //
    ch_wgs_in = ch_mosdepth_summary
        .join(ch_mosdepth_dist,              by: 0)
        .join(METHYLSEQ.out.meth_summary, by: 0)

    WGS_QC(ch_wgs_in)
    ch_versions        = ch_versions.mix(WGS_QC.out.versions)
    ch_multiqc_files   = ch_multiqc_files.mix(WGS_QC.out.mqc.map { it[1] }.flatten())
    ch_wgs_metrics_csv = WGS_QC.out.mqc

    //
    // Collate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'downstream_software_mqc_versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    rastair_mbias       = METHYLSEQ.out.mbias      // channel: [ val(meta), path(*.html)    ]
    rastair_bed         = METHYLSEQ.out.bed        // channel: [ val(meta), path(*.bed.gz)  ]
    rastair_vcf         = METHYLSEQ.out.vcf        // channel: [ val(meta), path(*.vcf.gz)  ]
    rastair_perread     = METHYLSEQ.out.perread    // channel: [ val(meta), path(*.bed.gz)  ]
    rastair_methylkit   = METHYLSEQ.out.methylkit  // channel: [ val(meta), path(*.txt.gz) ]
    cpg_mask            = CPG_MASK.out.bed           // channel: [ val(meta), path(*.rastair_cpg_sites.bed.gz) ]
    cpg_mask_tbi        = CPG_MASK.out.tbi           // channel: [ val(meta), path(*.rastair_cpg_sites.bed.gz.tbi) ]
    meth_summary        = METHYLSEQ.out.meth_summary // channel: [ val(meta), path(*.methylation_summary.tsv) ]
    mapping_metrics_csv = ch_mapping_metrics_csv              // channel: [ val(meta), path(*_mqc.tsv) ] — empty if no ctrls
    methyl_metrics_csv  = ch_methyl_metrics_csv               // channel: [ val(meta), path(*_mqc.tsv) ] — empty if no ctrls
    wgs_metrics_csv     = ch_wgs_metrics_csv                  // channel: [ val(meta), path(*_mqc.tsv) ]
    multiqc_files       = ch_multiqc_files                    // channel: all mqc files for MULTIQC
    versions            = ch_versions
}
