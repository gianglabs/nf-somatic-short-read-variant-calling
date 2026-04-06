#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MANTA } from '../../../../modules/local/manta/main'

workflow STRUCTURAL_VARIANT_CALLING {
    take:
    structural_variant_caller // value: manta
    tn_pairs // channel: [ val(pair_meta), path(normal_align), path(normal_index), path(tumor_align), path(tumor_index) ]
    ref_fasta // path(fasta)
    ref_fai // path(fai)

    main:
    ch_versions = channel.empty()
    ch_out_vcf = channel.empty()
    ch_out_vcf_tbi = channel.empty()

    callers = structural_variant_caller.split(',').collect { it.trim().toLowerCase() }.findAll { it }
    if (!callers) {
        error("No structural variant caller provided. Use --structural_variant_caller with: manta")
    }

    if (callers.any { !(it in ['manta']) }) {
        error("Unsupported structural variant caller(s): ${callers}. Supported callers: manta")
    }

    if (callers.contains('manta')) {
        MANTA(
            tn_pairs,
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(MANTA.out.versions)

        ch_out_vcf = ch_out_vcf.mix(
            MANTA.out.somatic_sv_vcf.map { meta, vcf -> [meta + [id: "${meta.id}.manta.somatic_sv"], vcf] }
        )
        ch_out_vcf = ch_out_vcf.mix(
            MANTA.out.candidate_small_indels_vcf.map { meta, vcf -> [meta + [id: "${meta.id}.manta.candidate_indels"], vcf] }
        )

        ch_out_vcf_tbi = ch_out_vcf_tbi.mix(
            MANTA.out.somatic_sv_vcf_tbi.map { meta, tbi -> [meta + [id: "${meta.id}.manta.somatic_sv"], tbi] }
        )
        ch_out_vcf_tbi = ch_out_vcf_tbi.mix(
            MANTA.out.candidate_small_indels_vcf_tbi.map { meta, tbi -> [meta + [id: "${meta.id}.manta.candidate_indels"], tbi] }
        )
    }

    emit:
    vcf = ch_out_vcf
    vcf_tbi = ch_out_vcf_tbi
    versions = ch_versions
}
