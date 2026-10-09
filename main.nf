#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TheJacksonLaboratory/5baseTAPS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/TheJacksonLaboratory/5baseTAPS
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQUORUM          } from './workflows/fastquorum'
include { RASTAIR              } from './workflows/rastair'
include { GATK_VARIANTCALL as GATK } from './workflows/gatk_variantcall'
include { TAPS_QC_REPORT      } from './modules/local/taps_qc_report/main'
include { TAPS_MULTIQC        } from './modules/local/taps_multiqc/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_fastquorum_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_fastquorum_pipeline'
include { BWAMEM2_INDEX  } from './modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_DICT  } from './modules/nf-core/samtools/dict/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_fastquorum_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow {
    // Populate genome params from the genome config when not explicitly provided.
    // Must run before any subworkflow that reads params.fasta / params.bwamem2.
    if (!params.fasta     && getGenomeAttribute('fasta'))     params.replace('fasta',     getGenomeAttribute('fasta'))
    if (!params.fasta_fai && getGenomeAttribute('fasta_fai')) params.replace('fasta_fai', getGenomeAttribute('fasta_fai'))
    if (!params.dict      && getGenomeAttribute('dict'))      params.replace('dict',      getGenomeAttribute('dict'))
    if (!params.bwamem2   && getGenomeAttribute('bwamem2'))   params.replace('bwamem2',   getGenomeAttribute('bwamem2'))

    // If --fasta was explicitly provided and differs from the genome config's fasta,
    // genome-derived indexes won't match — null them so they get rebuilt.
    if (params.fasta && getGenomeAttribute('fasta') && params.fasta != getGenomeAttribute('fasta')) {
        params.replace('fasta_fai', null)
        params.replace('dict',      null)
        params.replace('bwamem2',   null)
    }

    // Reset to null if local index files no longer exist; triggers rebuild from FASTA.
    if (params.fasta_fai && !file(params.fasta_fai).exists()) {
        log.warn "fasta_fai not found at ${params.fasta_fai} — will rebuild from FASTA"
        params.replace('fasta_fai', null)
    }
    if (params.dict && !file(params.dict).exists()) {
        log.warn "dict not found at ${params.dict} — will rebuild from FASTA"
        params.replace('dict', null)
    }
    if (params.bwamem2) {
        def bwa_dir = file(params.bwamem2)
        if (!bwa_dir.exists() || bwa_dir.list()?.size() == 0) {
            log.warn "bwamem2 index not found or empty at ${params.bwamem2} — will rebuild from FASTA"
            params.replace('bwamem2', null)
        }
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    // Skip main workflow when --help/--help_full is requested.
    if (!params.help && !params.help_full) {
        //
        // WORKFLOW: Run main workflow
        //
        JAXGT_5BASETAPS(
            PIPELINE_INITIALISATION.out.samplesheet
        )
        //
        // SUBWORKFLOW: Run completion tasks
        //
        PIPELINE_COMPLETION(
            params.email,
            params.email_on_fail,
            params.plaintext_email,
            params.outdir,
            params.monochrome_logs,
            params.hook_url,
            JAXGT_5BASETAPS.out.multiqc_report,
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow JAXGT_5BASETAPS {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    // Initialize fasta file with meta map:
    fasta = params.fasta ? Channel.fromPath(params.fasta).map { it -> [[id: it.baseName], it] }.collect() : Channel.empty()

    // Set various consensus calling and filtering parameters if not given
    if (params.duplex_seq) {
        if (params.groupreadsbyumi_strategy == '') {
            params.replace("groupreadsbyumi_strategy", 'Paired')
        }
        else if (params.groupreadsbyumi_strategy != 'Paired') {
            log.error("config groupreadsbyumi_strategy must be 'Paired' for duplex-sequencing data")
            exit(1)
        }
        if (params.call_min_reads == '') {
            params.replace("call_min_reads", '1 0 0')
        }
        if (!params.filter_min_reads) {
            params.replace("filter_min_reads", '1 0 0') 
        }
    }
    else {
        if (params.groupreadsbyumi_strategy == '') {
            params.replace("groupreadsbyumi_strategy", 'Adjacency')
        }
        else if (params.groupreadsbyumi_strategy == 'Paired') {
            log.error("config groupreadsbyumi_strategy cannot be 'Paired' for non-duplex-sequencing data")
            exit(1)
        }
        if (params.call_min_reads == '') {
            params.replace("call_min_reads", '1')
        }
        if (params.filter_min_reads == '') {
            params.replace("filter_min_reads", '1')
        }
    }

    if (!params.bwamem2)   BWAMEM2_INDEX(fasta)
    if (!params.fasta_fai) SAMTOOLS_FAIDX(fasta, [[id: 'no_fai'], []])
    if (!params.dict)      SAMTOOLS_DICT(fasta)

    dict = params.dict
        ? Channel.fromPath(params.dict).map { it -> [[id: 'dict'], it] }.collect()
        : SAMTOOLS_DICT.out.dict.collect()
    fasta_fai = params.fasta_fai
        ? Channel.fromPath(params.fasta_fai).map { it -> [[id: 'fai'], it] }.collect()
        : SAMTOOLS_FAIDX.out.fai.collect()
    bwamem2 = params.bwamem2
        ? Channel.fromPath(params.bwamem2, type: 'dir').map { it -> [[id: 'bwamem2'], it] }.collect()
        : BWAMEM2_INDEX.out.index.collect()
    //
    // WORKFLOW: Run pipeline
    //
    FASTQUORUM(
        params,
        samplesheet,
        bwamem2,
        dict,
        fasta,
        fasta_fai,
    )

    //
    // WORKFLOW: Downstream analysis — methylation calling, variant calling, QC
    //
    RASTAIR(
        FASTQUORUM.out.bam,
        FASTQUORUM.out.bai,
        fasta,
        fasta_fai,
        FASTQUORUM.out.unmapped_bam,
        FASTQUORUM.out.prededup_flagstat,
        FASTQUORUM.out.mosdepth_summary,
        FASTQUORUM.out.mosdepth_dist,
    )

    //
    // WORKFLOW: GATK variant calling — non-CpG variants, CpG sites excluded via rastair mask
    //   Skip with --run_gatk false for methylation-only runs.
    //
    ch_vc_metrics_csv  = Channel.empty()
    ch_gatk_mqc_files  = Channel.empty()

    if (params.run_gatk) {
        ch_cpg_mask = RASTAIR.out.cpg_mask
            .join(RASTAIR.out.cpg_mask_tbi)

        GATK(
            FASTQUORUM.out.bam,
            FASTQUORUM.out.bai,
            fasta,
            fasta_fai,
            dict,
            ch_cpg_mask,
            RASTAIR.out.meth_summary,
        )

        ch_vc_metrics_csv  = GATK.out.vc_metrics_csv
        ch_gatk_mqc_files  = GATK.out.mqc_files
    }

    //
    // MODULE: TAPS_MULTIQC — run MultiQC on explicitly collected outputs after all workflows complete.
    //
    ch_multiqc_inputs = FASTQUORUM.out.multiqc_files
        .mix(RASTAIR.out.multiqc_files)
        .mix(ch_gatk_mqc_files)
        .flatten()
        .collect()

    TAPS_MULTIQC(
        ch_multiqc_inputs,
        file("${projectDir}/assets/multiqc_config.yml"),
        file("${projectDir}/assets/JAX_logo_rgb_transparentback.png"),
    )

    emit:
    multiqc_report = TAPS_MULTIQC.out.report.toList()  // channel: /path/to/multiqc_report.html
}
