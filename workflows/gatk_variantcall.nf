/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { EXTRACT_MAJOR_CHROMS          } from '../modules/local/extract_major_chroms/main'
include { GATK4_COMPOSESTRTABLEFILE     } from '../modules/local/gatk4/composestrtablefile/main'
include { GATK4_CALIBRATEDRAGSTRMODEL   } from '../modules/local/gatk4/calibratedragstrmodel/main'
include { GATK4_HAPLOTYPECALLER as HaplotypeCaller } from '../modules/local/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS               } from '../modules/local/gatk4/mergevcfs/main'
include { GATK4_VARIANTFILTRATION       } from '../modules/local/gatk4/variantfiltration/main'
include { BCFTOOLS_SPLITVCF             } from '../modules/local/bcftools/splitvcf/main'
include { BCFTOOLS_STATS                } from '../modules/local/bcftools/stats/main'
include { TAPS_VC_METRICS               } from '../modules/local/taps_vc_metrics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GATK VARIANT CALLING WORKFLOW
    Receives the final consensus-filtered BAM, the rastair CpG mask, and the rastair
    methylation summary (used to derive the per-sample major-chromosome scatter list),
    then runs:
      1. ComposeSTRTableFile   : build genome STR table (once, genome-level)
      2. CalibrateDragstrModel : fit DRAGstr model per sample
      3. ExtractMajorChroms    : derive scatter list (chr1..22, X, Y) per sample
      4. HaplotypeCaller       : dragen-mode variant calling, scattered per chromosome
      5. MergeVcfs             : gather per-chrom VCFs into one per-sample VCF
      6. VariantFiltration     : DRAGENHardQUAL hard filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GATK_VARIANTCALL {

    take:
    ch_bam          // channel: [ val(meta), path(bam) ]
    ch_bai          // channel: [ val(meta), path(bai) ]
    ch_fasta        // channel: [ val(meta), path(fasta) ]   (value channel)
    ch_fasta_fai    // channel: [ val(meta), path(fai) ]     (value channel)
    ch_dict         // channel: [ val(meta), path(dict) ]    (value channel)
    ch_cpg_mask     // channel: [ val(meta), path(bed.gz), path(bed.gz.tbi) ]
    ch_meth_summary // channel: [ val(meta), path(methylation_summary.tsv) ]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: ComposeSTRTableFile — genome-level STR table, runs once
    //
    GATK4_COMPOSESTRTABLEFILE(ch_fasta, ch_fasta_fai, ch_dict)
    ch_versions = ch_versions.mix(GATK4_COMPOSESTRTABLEFILE.out.versions)

    ch_str_table = GATK4_COMPOSESTRTABLEFILE.out.str_table

    //
    // MODULE: CalibrateDragstrModel — fit DRAGstr model per sample
    //
    GATK4_CALIBRATEDRAGSTRMODEL(
        ch_bam,
        ch_bai,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_str_table,
    )
    ch_versions = ch_versions.mix(GATK4_CALIBRATEDRAGSTRMODEL.out.versions)

    ch_dragstr_model = GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model

    //
    // MODULE: ExtractMajorChroms — derive per-sample chromosome scatter list
    //   from rastair methylation_summary.tsv (col 2 filtered to chr1..22, X, Y).
    //
    EXTRACT_MAJOR_CHROMS(ch_meth_summary)
    ch_versions = ch_versions.mix(EXTRACT_MAJOR_CHROMS.out.versions)

    // Per-sample chromosome list (parsed once from the file).
    ch_sample_chroms = EXTRACT_MAJOR_CHROMS.out.chroms
        .map { meta, file ->
            def chroms = file.readLines().collect { it.trim() }.findAll { it }
            [ meta, chroms ]
        }
    // ch_sample_chroms: [ meta, [chr1, chr2, ...] ]   — one tuple per sample

    //
    // Build the per-(sample × chrom) HC input bundle.
    //
    ch_hc_sample = ch_bam
        .join(ch_bai,           by: 0)
        .join(ch_dragstr_model, by: 0)
        .join(ch_cpg_mask,      by: 0)
        .join(ch_sample_chroms, by: 0)
    // ch_hc_sample: [ meta, bam, bai, dragstr_model, cpg, cpg_tbi, [chroms] ]

    ch_hc_in = ch_hc_sample
        .flatMap { meta, bam, bai, dm, cpg, tbi, chroms ->
            def meta_n = meta + [ num_chroms: chroms.size() ]
            chroms.collect { chrom -> [ meta_n, bam, bai, dm, cpg, tbi, chrom ] }
        }
    // ch_hc_in: [ meta+num_chroms, bam, bai, dragstr_model, cpg, cpg_tbi, chrom ]

    //
    // MODULE: HaplotypeCaller — scattered per chromosome, dragen-mode
    //
    HaplotypeCaller(
        ch_hc_in.map { meta, bam, bai, dm, cpg, tbi, chrom -> [meta, bam]       },
        ch_hc_in.map { meta, bam, bai, dm, cpg, tbi, chrom -> [meta, bai]       },
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_hc_in.map { meta, bam, bai, dm, cpg, tbi, chrom -> [meta, cpg, tbi]  },
        ch_hc_in.map { meta, bam, bai, dm, cpg, tbi, chrom -> [meta, dm]        },
        ch_hc_in.map { meta, bam, bai, dm, cpg, tbi, chrom -> chrom             },
    )
    ch_versions = ch_versions.mix(HaplotypeCaller.out.versions)

    //
    // MODULE: MergeVcfs — gather per-chrom VCFs into one per-sample VCF.
    //
    ch_vcfs_to_merge = HaplotypeCaller.out.vcf
        .map { meta, vcf -> [ groupKey(meta, meta.num_chroms), vcf ] }
        .groupTuple(by: 0)
    // ch_vcfs_to_merge: [ meta, [vcf1, vcf2, ..., vcfN] ]   — one tuple per sample

    GATK4_MERGEVCFS(
        ch_vcfs_to_merge,
        ch_dict,
    )
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

    //
    // MODULE: VariantFiltration — DRAGEN-mode hard filter on merged VCF.
    //   QUAL < 10.4139 ("DRAGENHardQUAL") replaces VQSR and CNNScoreVariants
    //   for single-sample DRAGEN-mode calls.
    //
    GATK4_VARIANTFILTRATION(
        GATK4_MERGEVCFS.out.vcf,
        GATK4_MERGEVCFS.out.tbi,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions)

    //
    // MODULE: Split filtered VCF into SNP and INDEL
    //
    ch_filtered_for_split = GATK4_VARIANTFILTRATION.out.vcf
        .join(GATK4_VARIANTFILTRATION.out.tbi, by: 0)

    BCFTOOLS_SPLITVCF(ch_filtered_for_split)
    ch_versions = ch_versions.mix(BCFTOOLS_SPLITVCF.out.versions)

    //
    // MODULE: bcftools stats on the combined filtered VCF (for Ti/Tv etc.)
    //
    ch_vcf_for_stats = GATK4_VARIANTFILTRATION.out.vcf
        .join(GATK4_VARIANTFILTRATION.out.tbi, by: 0)

    BCFTOOLS_STATS(ch_vcf_for_stats)
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { it[1] }.collect())

    //
    // MODULE: TAPS_VC_METRICS — variant summary for MultiQC
    //
    TAPS_VC_METRICS(BCFTOOLS_STATS.out.stats)
    ch_versions      = ch_versions.mix(TAPS_VC_METRICS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TAPS_VC_METRICS.out.mqc.map { it[1] }.flatten())

    //
    // Collate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'gatk_variantcall_software_mqc_versions.yml',
            sort: true,
            newLine: true,
        )

    emit:
    snp_vcf        = BCFTOOLS_SPLITVCF.out.snp_vcf      // channel: [ val(meta), path(*.filtered-SNP.vcf.gz)      ]
    snp_tbi        = BCFTOOLS_SPLITVCF.out.snp_tbi      // channel: [ val(meta), path(*.filtered-SNP.vcf.gz.tbi)  ]
    indel_vcf      = BCFTOOLS_SPLITVCF.out.indel_vcf    // channel: [ val(meta), path(*.filtered-INDEL.vcf.gz)     ]
    indel_tbi      = BCFTOOLS_SPLITVCF.out.indel_tbi    // channel: [ val(meta), path(*.filtered-INDEL.vcf.gz.tbi) ]
    vcf_raw        = GATK4_MERGEVCFS.out.vcf            // channel: [ val(meta), path(*.haplotypecaller.vcf.gz)    ]
    tbi_raw        = GATK4_MERGEVCFS.out.tbi            // channel: [ val(meta), path(*.haplotypecaller.vcf.gz.tbi)]
    vc_metrics_csv   = TAPS_VC_METRICS.out.mqc             // channel: [ val(meta), path(*_mqc.tsv) ]
    mqc_files        = ch_multiqc_files                    // channel: path(*_mqc.tsv)
    versions         = ch_versions
}
