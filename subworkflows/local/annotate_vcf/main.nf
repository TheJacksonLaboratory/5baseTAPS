/*
 * ANNOTATE_VCF: add annotated copies of the GATK HaplotypeCaller and rastair VCFs (GRCh38 only).
 * The original VCFs are never modified; outputs are new <sample>.<callset>.annotated.vcf.gz files.
 *
 *   per VCF, per chromosome (scatter):
 *     BCFTOOLS_ANNOTATERESOURCES  dbSNP rsID, gnomAD AF, ClinVar
 *     SNPEFF_ANNOTATE             functional consequence (ANN)
 *     SNPSIFT_DBNSFP              SIFT / PolyPhen2 / REVEL / CADD (dbNSFP, bare-chromosome copy)
 *   then BCFTOOLS_CONCATANNOTATED gathers the chromosomes into one VCF + index.
 *
 * rastair VCF: only records with a real ALT allele are annotated (BCFTOOLS_SELECTALT);
 * the ALT='.' CpG methylation rows (~99% of the file) carry no variant.
 */

include { BCFTOOLS_SELECTALT          } from '../../../modules/local/bcftools/selectalt/main'
include { BCFTOOLS_ANNOTATERESOURCES  } from '../../../modules/local/bcftools/annotateresources/main'
include { SNPEFF_ANNOTATE             } from '../../../modules/local/snpeff/annotate/main'
include { SNPSIFT_DBNSFP              } from '../../../modules/local/snpsift/dbnsfp/main'
include { BCFTOOLS_CONCATANNOTATED    } from '../../../modules/local/bcftools/concatannotated/main'

workflow ANNOTATE_VCF {

    take:
    ch_indexed_vcf  // channel: [ val(meta), val(callset), path(vcf), path(tbi) ]
    ch_rastair_vcf  // channel: [ val(meta), val(callset), path(vcf) ]   (not indexed, mostly ALT='.')
    chroms          // list of chromosome names to scatter over
    dbsnp           // val: path to tabix-indexed dbSNP VCF
    gnomad          // val: path to tabix-indexed gnomAD af-only VCF
    clinvar         // val: path to tabix-indexed, chr-prefixed ClinVar VCF
    dbnsfp          // val: path to bare-chromosome-name dbNSFP database (.txt.gz + .tbi + .data_types)
    snpeff_db       // val: SnpEff database name (hg38)
    snpeff_config   // val: path to snpEff.config

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_SELECTALT(ch_rastair_vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SELECTALT.out.versions.first())

    ch_all = ch_indexed_vcf.mix(BCFTOOLS_SELECTALT.out.vcf)

    // one task per (sample, callset, chromosome)
    ch_chunks = ch_all.flatMap { meta, callset, vcf, tbi ->
        chroms.collect { chrom -> [ meta, callset, vcf, tbi, chrom ] }
    }

    BCFTOOLS_ANNOTATERESOURCES(ch_chunks, dbsnp, gnomad, clinvar)
    SNPEFF_ANNOTATE(BCFTOOLS_ANNOTATERESOURCES.out.vcf, snpeff_db, snpeff_config)
    SNPSIFT_DBNSFP(SNPEFF_ANNOTATE.out.vcf, dbnsfp)

    ch_versions = ch_versions.mix(
        BCFTOOLS_ANNOTATERESOURCES.out.versions.first(),
        SNPEFF_ANNOTATE.out.versions.first(),
        SNPSIFT_DBNSFP.out.versions.first(),
    )

    // gather the chromosomes of each (sample, callset); groupKey releases a group as soon as it is complete
    n_chroms = chroms.size()
    ch_grouped = SNPSIFT_DBNSFP.out.vcf
        .map { meta, callset, chrom, vcf -> [ groupKey("${meta.id}:${callset}".toString(), n_chroms), meta, callset, vcf ] }
        .groupTuple(by: 0)
        .map { key, metas, callsets, vcfs -> [ metas[0], callsets[0], vcfs ] }

    BCFTOOLS_CONCATANNOTATED(ch_grouped)
    ch_versions = ch_versions.mix(BCFTOOLS_CONCATANNOTATED.out.versions.first())

    emit:
    vcf      = BCFTOOLS_CONCATANNOTATED.out.vcf   // channel: [ val(meta), val(callset), path(vcf.gz), path(tbi) ]
    versions = ch_versions
}
