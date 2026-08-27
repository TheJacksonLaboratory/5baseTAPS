/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_fastquorum_pipeline'

include { ALIGN_BAM as ALIGN_RAW_BAM    } from '../modules/local/align_bam/main'
include { ALIGN_BAM as ALIGN_RAW_CHUNK  } from '../modules/local/align_bam/main'
include { ALIGN_BAM as ALIGN_CONS } from '../modules/local/align_bam/main'
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { CONCAT_FASTQ as CONCAT_FASTQ_R1 } from '../modules/local/concat_fastq/main'
include { CONCAT_FASTQ as CONCAT_FASTQ_R2 } from '../modules/local/concat_fastq/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM       } from '../modules/local/fgbio/fastqtobam/main'
include { FGBIO_FASTQTOBAM as FASTQTOBAM_CHUNK } from '../modules/local/fgbio/fastqtobam/main'
include { SPLIT_FASTQ   } from '../modules/local/split_fastq/main'
include { MERGE_CHUNKS  } from '../modules/local/merge_chunks/main'
include { GROUP_UMI } from '../modules/local/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLMOLECULARCONSENSUSREADS } from '../modules/local/fgbio/callmolecularconsensusreads/main'
include { CALL_DUPLEX } from '../modules/local/fgbio/callduplexconsensusreads/main'
include { CONSENSUS } from '../modules/local/fgbio/filterconsensusreads/main'
include { DUPLEX_QC } from '../modules/local/fgbio/collectduplexseqmetrics/main'
include { FGBIO_CALLANDFILTERMOLECULARCONSENSUSREADS as CALLANDFILTERMOLECULARCONSENSUSREADS } from '../modules/local/fgbio/callandfiltermolecularconsensusreads/main'
include { FGBIO_CALLANDFILTERDUPLEXCONSENSUSREADS as CALLANDFILTERDUPLEXCONSENSUSREADS } from '../modules/local/fgbio/callandfilterduplexconsensusreads/main'
include { SAMTOOLS_MERGE as MERGE_BAM              } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FLAGSTAT                         } from '../modules/local/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as PREDEDUP_FLAGSTAT    } from '../modules/local/samtools/flagstat/main'
include { SAMTOOLS_FLAGSTAT as POSTDEDUP_FLAGSTAT   } from '../modules/local/samtools/flagstat/main'
include { SAMTOOLS_STATS                            } from '../modules/local/samtools/stats/main'
include { MOSDEPTH                                  } from '../modules/local/mosdepth/main'
include { DUPLEX_MQC                       } from '../modules/local/taps_duplex_metrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTQUORUM {
    take:
    params
    ch_samplesheet
    ch_bwamem2
    ch_dict
    ch_fasta
    ch_fasta_fai

    main:

    // To gather all QC reports for MultiQC
    ch_versions           = Channel.empty()
    ch_multiqc_files      = Channel.empty()
    ch_final_bam          = Channel.empty()
    ch_final_bai          = Channel.empty()
    ch_prededup_flagstat  = Channel.empty()
    ch_duplex_metrics_csv = Channel.empty()

    //
    // MODULE: Concatenate per-lane FASTQs per sample for multi-lane samples;
    //         single-lane samples skip concat and go directly to FastQC.
    //         Multi-lane: 2N tasks (R1 and R2 per sample run separately).
    //
    ch_by_sample = ch_samplesheet
        .map { meta, reads ->
            def meta_nolane = meta.findAll { k, v -> k != 'lane' }
            [groupKey(meta_nolane, meta_nolane.n_samples), reads]
        }
        .groupTuple()
        .branch { meta, reads_list ->
            single:   meta.n_samples <= 1
                return [meta, reads_list[0]]
            multiple: true
        }

    ch_multi_by_read = ch_by_sample.multiple
        .multiMap { meta, reads_list ->
            r1: [meta, reads_list.collect { it[0] }]
            r2: [meta, reads_list.collect { it[1] }]
        }

    CONCAT_FASTQ_R1(ch_multi_by_read.r1)
    CONCAT_FASTQ_R2(ch_multi_by_read.r2)
    ch_versions = ch_versions.mix(CONCAT_FASTQ_R1.out.versions.first())

    ch_fastqc_input = ch_by_sample.single.mix(
        CONCAT_FASTQ_R1.out.reads
            .join(CONCAT_FASTQ_R2.out.reads)
            .map { meta, r1, r2 -> [meta, [r1, r2]] }
    )
    FASTQC(ch_fastqc_input)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })

    //
    // FASTQ → unmapped BAM → raw aligned BAM (TC-sorted, per lane)
    // Branch: split-align-merge vs original single-task path
    //
    if (params.align_raw_bam_chunks > 1) {
        //
        // Split each FASTQ pair into N chunks, align in parallel, TC-sort merge
        //
        SPLIT_FASTQ(ch_samplesheet, params.align_raw_bam_chunks)
        ch_versions = ch_versions.mix(SPLIT_FASTQ.out.versions.first())

        // Flatten N chunk pairs into individual channel items: [meta_with_chunk, [R1_chunk, R2_chunk]]
        ch_chunks = SPLIT_FASTQ.out.reads
            .flatMap { meta, r1_list, r2_list ->
                def r1s = r1_list instanceof List ? r1_list : [r1_list]
                def r2s = r2_list instanceof List ? r2_list : [r2_list]
                [r1s, r2s].transpose().withIndex().collect { pair, idx ->
                    [meta + [chunk: idx], pair]
                }
            }

        FASTQTOBAM_CHUNK(ch_chunks)
        ch_versions = ch_versions.mix(FASTQTOBAM_CHUNK.out.versions.first())

        ALIGN_RAW_CHUNK(FASTQTOBAM_CHUNK.out.bam, ch_fasta, ch_fasta_fai, ch_dict, ch_bwamem2, "none")
        ch_versions = ch_versions.mix(ALIGN_RAW_CHUNK.out.versions.first())

        // Group chunk BAMs back by sample/lane, then TC-sort merge into one BAM per lane
        MERGE_CHUNKS(
            ALIGN_RAW_CHUNK.out.bam
                .map { meta_chunk, bam -> [meta_chunk.findAll { k, v -> k != 'chunk' }, bam] }
                .groupTuple(size: params.align_raw_bam_chunks as Integer)
        )
        ch_versions = ch_versions.mix(MERGE_CHUNKS.out.versions.first())

        // Full-sample unmapped BAM for ctrl workflows (lambda/pUC19) — one per sample, not per chunk
        FASTQTOBAM(ch_samplesheet)
        ch_versions = ch_versions.mix(FASTQTOBAM.out.versions.first())

        ch_raw_bam      = MERGE_CHUNKS.out.bam
        ch_unmapped_bam = FASTQTOBAM.out.bam

    } else {
        //
        // Original path: single-task fgbio FastqToBam + bwa-mem2 alignment
        //
        FASTQTOBAM(ch_samplesheet)
        ch_versions = ch_versions.mix(FASTQTOBAM.out.versions.first())

        ALIGN_RAW_BAM(FASTQTOBAM.out.bam, ch_fasta, ch_fasta_fai, ch_dict, ch_bwamem2, "template-coordinate")
        ch_versions = ch_versions.mix(ALIGN_RAW_BAM.out.versions.first())

        ch_raw_bam      = ALIGN_RAW_BAM.out.bam
        ch_unmapped_bam = FASTQTOBAM.out.bam
    }

    //
    // MODULE: flagstat on raw-aligned BAM (per lane, before merge / before dedup)
    //
    SAMTOOLS_FLAGSTAT(ch_raw_bam)
    ch_versions      = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { it[1] }.collect())

    //
    // Create a channel that:
    // 1. Groups the aligned BAMs by sample identifier.  Typically a sample has more than one BAM
    //    if it had multiple runs or lanes.
    // 2. Splits the samples into those that have more than one BAM, and those that have exactly one BAM.  The former
    //    samples will have their BAMs merged.
    //
    bam_to_merge = ch_raw_bam
        .map { meta, bam ->
            def meta_no_lane = meta.findAll { k, v -> k != 'lane' }
            [groupKey(meta_no_lane, meta_no_lane.n_samples), bam]
        }
        .groupTuple()
        .branch { meta, bam ->
            single: meta.n_samples <= 1
            return [meta, bam[0]]
            multiple: meta.n_samples > 1
        }

    //
    // MODULE: Run samtools merge to merge across runs/lanes for the same sample
    //
    MERGE_BAM(bam_to_merge.multiple, [[], []], [[], []])
    ch_versions = ch_versions.mix(MERGE_BAM.out.versions.first())

    //
    // Create a channel that contains the merged BAMs and those that did not need to be merged.
    //
    bam_all = MERGE_BAM.out.bam.mix(bam_to_merge.single)

    //
    // MODULE: flagstat on post-merge (pre-dedup) BAM — one per sample, correct total raw read count
    //
    PREDEDUP_FLAGSTAT(bam_all)
    ch_versions          = ch_versions.mix(PREDEDUP_FLAGSTAT.out.versions.first())
    ch_prededup_flagstat = PREDEDUP_FLAGSTAT.out.flagstat

    //
    // MODULE: Run fgbio GroupReadsByUmi
    //
    def umi_strategy = params.groupreadsbyumi_strategy ?: (params.duplex_seq ? 'Paired' : 'Adjacency')
    log.info("[fgbio GroupReadsByUmi] strategy='${umi_strategy}' (duplex_seq=${params.duplex_seq})")
    GROUP_UMI(bam_all, umi_strategy, params.groupreadsbyumi_edits)
    ch_multiqc_files = ch_multiqc_files.mix(GROUP_UMI.out.histogram.map { it[1] }.collect())
    ch_versions = ch_versions.mix(GROUP_UMI.out.versions.first())

    if (params.duplex_seq) {
        //
        // MODULE: Run fgbio CollectDuplexSeqMetrics
        //
        DUPLEX_QC(GROUP_UMI.out.bam)
        ch_versions = ch_versions.mix(DUPLEX_QC.out.versions.first())

        ch_collectduplex_metrics = DUPLEX_QC.out.metrics
    }

    if (params.mode == 'rd') {
        if (params.duplex_seq) {
            //
            // MODULE: Run fgbio CallDuplexConsensusReads
            //
            CALL_DUPLEX(GROUP_UMI.out.bam, params.call_min_reads, params.call_min_baseq)
            ch_versions = ch_versions.mix(CALL_DUPLEX.out.versions.first())

            // Add the consensus BAM to the channel for downstream processing
            CALL_DUPLEX.out.bam.set { ch_consensus_bam }
        }
        else {
            //
            // MODULE: Run fgbio CallMolecularConsensusReads
            //
            CALLMOLECULARCONSENSUSREADS(GROUP_UMI.out.bam, params.call_min_reads, params.call_min_baseq)
            ch_versions = ch_versions.mix(CALLMOLECULARCONSENSUSREADS.out.versions.first())

            // Add the consensus BAM to the channel for downstream processing
            CALLMOLECULARCONSENSUSREADS.out.bam.set { ch_consensus_bam }
        }

        //
        // MODULE: Align with bwa mem
        //
        ALIGN_CONS(ch_consensus_bam, ch_fasta, ch_fasta_fai, ch_dict, ch_bwamem2, "none")
        ch_versions = ch_versions.mix(ALIGN_CONS.out.versions.first())

        //
        // MODULE: Run fgbio FilterConsensusReads
        //
        CONSENSUS(ALIGN_CONS.out.bam, ch_fasta, params.filter_min_reads, params.filter_min_baseq, params.filter_max_base_error_rate)
        ch_versions  = ch_versions.mix(CONSENSUS.out.versions.first())
        ch_final_bam = CONSENSUS.out.bam
        ch_final_bai = CONSENSUS.out.bai
    }
    else {
        if (params.duplex_seq) {
            //
            // MODULE: Run fgbio CallDuplexConsensusReads and fgbio FilterConsensusReads
            //
            CALLANDFILTERDUPLEXCONSENSUSREADS(GROUP_UMI.out.bam, ch_fasta, ch_fasta_fai, params.call_min_reads, params.call_min_baseq, params.filter_max_base_error_rate)
            ch_versions = ch_versions.mix(CALLANDFILTERDUPLEXCONSENSUSREADS.out.versions.first())

            // Add the consensus BAM to the channel for downstream processing
            CALLANDFILTERDUPLEXCONSENSUSREADS.out.bam.set { ch_consensus_bam }
        }
        else {
            //
            // MODULE: Run fgbio CallMolecularConsensusReads and fgbio FilterConsensusReads
            //
            CALLANDFILTERMOLECULARCONSENSUSREADS(GROUP_UMI.out.bam, ch_fasta, ch_fasta_fai, params.call_min_reads, params.call_min_baseq, params.filter_max_base_error_rate)
            ch_versions = ch_versions.mix(CALLANDFILTERMOLECULARCONSENSUSREADS.out.versions.first())

            // Add the consensus BAM to the channel for downstream processing
            CALLANDFILTERMOLECULARCONSENSUSREADS.out.bam.set { ch_consensus_bam }
        }

        //
        // MODULE: Align with bwa mem
        //
        ALIGN_CONS(ch_consensus_bam, ch_fasta, ch_fasta_fai, ch_dict, ch_bwamem2, "coordinate")
        ch_versions  = ch_versions.mix(ALIGN_CONS.out.versions.first())
        ch_final_bam = ALIGN_CONS.out.bam
        ch_final_bai = ALIGN_CONS.out.bai
    }

    //
    // MODULE: flagstat on post-dedup consensus BAM (for dedup fold-reduction metric)
    //
    POSTDEDUP_FLAGSTAT(ch_final_bam)
    ch_versions = ch_versions.mix(POSTDEDUP_FLAGSTAT.out.versions.first())

    //
    // MODULE: samtools stats on consensus BAM (insert size distribution for MultiQC)
    //
    SAMTOOLS_STATS(ch_final_bam)
    ch_versions      = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { it[1] }.collect())

    //
    // MODULE: mosdepth WGS coverage on consensus BAM
    //
    ch_mosdepth_in = ch_final_bam.join(ch_final_bai)
    MOSDEPTH(ch_mosdepth_in)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    //
    // MODULE: DUPLEX_MQC — duplex deduplication summary CSV
    //   Only available when duplex_seq is true (requires DUPLEX_QC output).
    //
    if (params.duplex_seq) {
        ch_duplex_in = ch_prededup_flagstat
            .join(POSTDEDUP_FLAGSTAT.out.flagstat, by: 0)
            .join(ch_collectduplex_metrics,         by: 0)

        DUPLEX_MQC(ch_duplex_in)
        ch_versions           = ch_versions.mix(DUPLEX_MQC.out.versions.first())
        ch_multiqc_files      = ch_multiqc_files.mix(DUPLEX_MQC.out.mqc.map { it[1] }.flatten())
        ch_duplex_metrics_csv = DUPLEX_MQC.out.mqc
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'fastquorum_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/JAX_logo_rgb_transparentback.png", checkIfExists: false)

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    emit:
    multiqc_files       = ch_multiqc_files                // channel: all mqc input files for MULTIQC
    multiqc_config      = ch_multiqc_config               // channel: multiqc_config.yml
    multiqc_custom_config = ch_multiqc_custom_config      // channel: user --multiqc_config override
    multiqc_logo        = ch_multiqc_logo                 // channel: JAX logo
    versions            = ch_versions                      // channel: [ path(versions.yml) ]
    bam                 = ch_final_bam                     // channel: [ val(meta), path(bam) ]
    bai                 = ch_final_bai                     // channel: [ val(meta), path(bai) ]
    unmapped_bam        = ch_unmapped_bam                   // channel: [ val(meta), path(unmapped.bam) ] — one per samplesheet row (N per row if align_raw_bam_chunks > 1)
    prededup_flagstat   = ch_prededup_flagstat             // channel: [ val(meta), path(*.flagstat) ]
    postdedup_flagstat  = POSTDEDUP_FLAGSTAT.out.flagstat  // channel: [ val(meta), path(*.flagstat) ]
    mosdepth_summary    = MOSDEPTH.out.summary             // channel: [ val(meta), path(*.mosdepth.summary.txt) ]
    mosdepth_dist       = MOSDEPTH.out.global_dist         // channel: [ val(meta), path(*.mosdepth.global.dist.txt) ]
    duplex_metrics_csv  = ch_duplex_metrics_csv            // channel: [ val(meta), path(*_mqc.tsv) ] — empty if !duplex_seq
}
