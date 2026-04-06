#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GATK_MUTECT2 } from '../../../../modules/local/gatk/mutect2/main'
include { GATK_LEARNREADORIENTATIONMODEL } from '../../../../modules/local/gatk/learnreadorientationmodel/main'
include { STRELKA } from '../../../../modules/local/strelka/main'

workflow SMALL_VARIANT_CALLING {
    take:
    variant_caller // value: mutect2,strelka
    tn_pairs // channel: [ val(pair_meta), path(normal_align), path(normal_index), path(tumor_align), path(tumor_index) ]
    ref_fasta // path(fasta)
    ref_fai // path(fai)
    ref_dict // path(dict)
    germline_resource_vcf // path(vcf)
    germline_resource_tbi // path(tbi)
    panel_of_normals_vcf // path(vcf)
    panel_of_normals_tbi // path(tbi)

    main:
    ch_versions = channel.empty()
    ch_out_vcf = channel.empty()
    ch_out_vcf_tbi = channel.empty()

    callers = variant_caller.split(',').collect { it.trim().toLowerCase() }.findAll { it }
    if (!callers) {
        error("No somatic small variant caller provided. Use --somatic_variant_caller with one or more callers: mutect2,strelka")
    }

    if (callers.any { !(it in ['mutect2', 'strelka']) }) {
        error("Unsupported somatic small variant caller(s): ${callers}. Supported callers: mutect2,strelka")
    }

    if (callers.contains('mutect2')) {
        GATK_MUTECT2(
            tn_pairs.map { pair_meta, normal_align, normal_index, tumor_align, tumor_index ->
                [pair_meta, [normal_align, tumor_align], [normal_index, tumor_index]]
            },
            ref_fasta,
            ref_fai,
            ref_dict,
            germline_resource_vcf,
            germline_resource_tbi,
            panel_of_normals_vcf,
            panel_of_normals_tbi,
        )

        ch_versions = ch_versions.mix(GATK_MUTECT2.out.versions)
        ch_out_vcf = ch_out_vcf.mix(GATK_MUTECT2.out.vcf)
        ch_out_vcf_tbi = ch_out_vcf_tbi.mix(GATK_MUTECT2.out.tbi)

        GATK_LEARNREADORIENTATIONMODEL(
            GATK_MUTECT2.out.f1r2
                .map { meta, f1r2 -> [meta, [f1r2]] }
        )
        ch_versions = ch_versions.mix(GATK_LEARNREADORIENTATIONMODEL.out.versions)
    }

    if (callers.contains('strelka')) {
        STRELKA(
            tn_pairs,
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(STRELKA.out.versions)

        ch_strelka_snv = STRELKA.out.vcf_snvs.join(STRELKA.out.vcf_snvs_tbi)
            .map { meta, vcf, tbi -> [meta + [id: "${meta.id}.strelka.snvs"], vcf, tbi] }

        ch_strelka_indel = STRELKA.out.vcf_indels.join(STRELKA.out.vcf_indels_tbi)
            .map { meta, vcf, tbi -> [meta + [id: "${meta.id}.strelka.indels"], vcf, tbi] }

        ch_out_vcf = ch_out_vcf.mix(ch_strelka_snv.map { meta, vcf, tbi -> [meta, vcf] })
        ch_out_vcf = ch_out_vcf.mix(ch_strelka_indel.map { meta, vcf, tbi -> [meta, vcf] })
        ch_out_vcf_tbi = ch_out_vcf_tbi.mix(ch_strelka_snv.map { meta, vcf, tbi -> [meta, tbi] })
        ch_out_vcf_tbi = ch_out_vcf_tbi.mix(ch_strelka_indel.map { meta, vcf, tbi -> [meta, tbi] })
    }

    emit:
    vcf = ch_out_vcf
    vcf_tbi = ch_out_vcf_tbi
    versions = ch_versions
}
