#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
// Include modules
include { MANTA } from '../../../../modules/local/manta/main'
include { SAMTOOLS_VIEW } from '../../../../modules/gianglabs/samtools/view/main'
include { TABIX_INDEX_VCF } from '../../../../modules/gianglabs/bcftools/index/main'

workflow STRUCTURAL_VARIANT_CALLING {
    take:
    structural_variant_caller // value: SV calling tool (e.g. "manta" etc.)
    bam // channel: [ val(meta), path(bam) ] or [ val(meta), path(cram) ]
    bai // channel: [ val(meta), path(bai) ] or [ val(meta), path(crai) ]
    ref_fasta // value: path(fasta)
    ref_fai // value: path(fai)
    genome // value: genome type

    main:
    ch_versions = channel.empty()
    ch_out_vcf = channel.empty()

    // Detect if input is CRAM or BAM and convert if necessary
    bam_input = bam.branch {
        cram: it[1].toString().endsWith('.cram')
        bam: true
    }

    // Convert CRAM to BAM if needed
    SAMTOOLS_VIEW(
        bam_input.cram.join(bai),
        ref_fasta,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    // Merge converted BAM with original BAM files
    bam_ready = bam_input.bam.join(bai).mix(SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_VIEW.out.bai))

    if (structural_variant_caller.split(",").contains("manta")) {
        MANTA(
            bam_ready,
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(MANTA.out.versions)
        ch_out_vcf = ch_out_vcf.mix(MANTA.out.vcf)
    }
    if (!structural_variant_caller.split(",").any { it.trim() in ["manta", "tiddit", "delly", "smoove", "cnvnator"] }) {
        error("Unsupported SV caller: ${structural_variant_caller}. Supported callers: manta, tiddit, delly, smoove, cnvnator")
    }

    // Group VCF files by sample ID for merge
    // TODO: configure to merge the tools smartly
    ch_vcf_grouped = ch_out_vcf.groupTuple(by: 0)
    TABIX_INDEX_VCF(
        ch_vcf_grouped
    )
    ch_versions = ch_versions.mix(TABIX_INDEX_VCF.out.versions)

    emit:
    vcf = TABIX_INDEX_VCF.out.vcf // channel: [ val(meta), path(vcf.gz) ]
    vcf_tbi = TABIX_INDEX_VCF.out.vcf_tbi // channel: [ val(meta), path(vcf.gz.tbi) ] 
    versions = ch_versions
}