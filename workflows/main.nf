include { ALIGNMENT } from '../subworkflows/gianglabs/alignment/main'
include { PREPROCESSING } from '../subworkflows/gianglabs/alignment_preprocessing/main'
include { SMALL_VARIANT_CALLING } from '../subworkflows/local/variant_calling/small/main'
include { STRUCTURAL_VARIANT_CALLING } from '../subworkflows/local/variant_calling/structural/main'
include { VARIANT_ANNOTATION } from '../subworkflows/gianglabs/variant_annotation/main'
include { VARIANT_ALIGNMENT_QUALITY_CONTROL } from '../subworkflows/gianglabs/variant_alignment_quality_control/main'
include { SAMTOOLS_VIEW } from '../modules/gianglabs/samtools/view/main'

workflow SOMATIC_VARIANT_CALLING {
    take:
    input_ch

    main:
    ch_versions = channel.empty()

    if (params.genome && params.genomes.containsKey(params.genome)) {
        def igenome_ref = params.genomes[params.genome]
        if (igenome_ref) {
            if (igenome_ref.fasta) {
                params.reference = igenome_ref.fasta
            }
            if (igenome_ref.fasta_fai) {
                params.reference_index = igenome_ref.fasta_fai
            }
            if (igenome_ref.dict) {
                params.reference_dict = igenome_ref.dict
            }
            if (igenome_ref.index_bwa2_reference != null) {
                params.index_bwa2_reference = igenome_ref.index_bwa2_reference
            }
            if (igenome_ref.bwa2_index) {
                params.bwa2_index = igenome_ref.bwa2_index
            }
            if (igenome_ref.dbsnp) {
                params.dbsnp = igenome_ref.dbsnp
            }
            if (igenome_ref.dbsnp_tbi) {
                params.dbsnp_tbi = igenome_ref.dbsnp_tbi
            }
            if (igenome_ref.known_indels) {
                params.known_indels = igenome_ref.known_indels
            }
            if (igenome_ref.known_indels_tbi) {
                params.known_indels_tbi = igenome_ref.known_indels_tbi
            }
            if (igenome_ref.germline_resource) {
                params.germline_resource = igenome_ref.germline_resource
            }
            if (igenome_ref.germline_resource_tbi) {
                params.germline_resource_tbi = igenome_ref.germline_resource_tbi
            }
            if (igenome_ref.panel_of_normals) {
                params.panel_of_normals = igenome_ref.panel_of_normals
            }
            if (igenome_ref.panel_of_normals_tbi) {
                params.panel_of_normals_tbi = igenome_ref.panel_of_normals_tbi
            }
            if (igenome_ref.snpeff_cache) {
                params.snpeff_cache = igenome_ref.snpeff_cache
            }
            if (igenome_ref.vep_cache_version) {
                params.vep_cache_version = igenome_ref.vep_cache_version
            }
            if (igenome_ref.vep_genome) {
                params.vep_genome = igenome_ref.vep_genome
            }
            if (igenome_ref.vep_species) {
                params.vep_species = igenome_ref.vep_species
            }
            if (igenome_ref.vep_cache) {
                params.vep_cache = igenome_ref.vep_cache
            }
        }
    }

    somatic_callers = params.somatic_variant_caller ?: params.small_variant_caller

    if (somatic_callers?.split(',')?.collect { it.trim().toLowerCase() }?.contains('mutect2')) {
        if (!params.germline_resource || !params.germline_resource_tbi || !params.panel_of_normals || !params.panel_of_normals_tbi) {
            error('Mutect2 requires --germline_resource, --germline_resource_tbi, --panel_of_normals, and --panel_of_normals_tbi')
        }
    }

    log.info(
        """
    ==============================================================================================================================
    nf-somatic-short-read-variant-calling:
     - Workflow                  : SOMATIC_VARIANT_CALLING
     - Small Variant Caller      : ${somatic_callers ? somatic_callers : 'None'}
     - Structural Variant Caller : ${params.structural_variant_caller ? params.structural_variant_caller : 'None'}
     - Annotation Caller         : ${params.annotation_caller ? params.annotation_caller : 'None'}
     - Loaded genomes set        : ${params.genome ? params.genome : 'None'}
     - Reference Genome          : ${params.reference}
     - Input Samplesheet         : ${params.input}
     - Output Directory          : ${params.outdir}
    ==============================================================================================================================
    """.stripIndent()
    )

    ref_fasta = channel.fromPath(params.reference, checkIfExists: true).collect()
    ref_fai = channel.fromPath(params.reference_index, checkIfExists: true).collect()
    ref_dict = channel.fromPath(params.reference_dict, checkIfExists: true).collect()
    dbsnp_vcf = channel.fromPath(params.dbsnp, checkIfExists: true).collect()
    dbsnp_tbi = channel.fromPath(params.dbsnp_tbi, checkIfExists: true).collect()
    known_indels_vcf = channel.fromPath(params.known_indels, checkIfExists: true).collect()
    known_indels_tbi = channel.fromPath(params.known_indels_tbi, checkIfExists: true).collect()
    germline_resource_vcf = params.germline_resource ? channel.fromPath(params.germline_resource, checkIfExists: true).collect() : channel.empty()
    germline_resource_tbi = params.germline_resource_tbi ? channel.fromPath(params.germline_resource_tbi, checkIfExists: true).collect() : channel.empty()
    panel_of_normals_vcf = params.panel_of_normals ? channel.fromPath(params.panel_of_normals, checkIfExists: true).collect() : channel.empty()
    panel_of_normals_tbi = params.panel_of_normals_tbi ? channel.fromPath(params.panel_of_normals_tbi, checkIfExists: true).collect() : channel.empty()

    ch_samples = input_ch
        .map { pair_meta, normal_data, tumor_data ->
            def normal_lane = pair_meta.normal_lane ?: 'L001'
            def tumor_lane = pair_meta.tumor_lane ?: 'L001'

            def normal_meta = [
                patient: pair_meta.patient,
                sex: pair_meta.sex,
                status: 0,
                sample: pair_meta.normal_id,
                id: pair_meta.normal_id,
                lane: normal_lane,
                read_group: "${pair_meta.normal_id}.${normal_lane}",
            ]

            def tumor_meta = [
                patient: pair_meta.patient,
                sex: pair_meta.sex,
                status: 1,
                sample: pair_meta.tumor_id,
                id: pair_meta.tumor_id,
                lane: tumor_lane,
                read_group: "${pair_meta.tumor_id}.${tumor_lane}",
            ]

            [
                [normal_meta, normal_data],
                [tumor_meta, tumor_data],
            ]
        }
        .flatMap { it }

    input_branched = ch_samples.branch {
        fastq: it[1].type == 'fastq'
        bam: it[1].type == 'bam'
        cram: it[1].type == 'cram'
    }

    ALIGNMENT(
        input_branched.fastq.map { meta, data -> [meta, data.reads] },
        ref_fasta,
        ref_fai,
        ref_dict,
        params.bwa2_index,
        params.index_bwa2_reference,
    )

    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
    ch_bam = ALIGNMENT.out.bam.mix(input_branched.bam.map { meta, data -> [meta, data.align] })
    ch_bai = ALIGNMENT.out.bai.mix(input_branched.bam.map { meta, data -> [meta, data.index] })

    SAMTOOLS_VIEW(
        input_branched.cram.map { meta, data -> [meta, data.align, data.index] },
        ref_fasta,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    ch_bam = ch_bam.mix(SAMTOOLS_VIEW.out.bam)
    ch_bai = ch_bai.mix(SAMTOOLS_VIEW.out.bai)

    if (params.preprocessor) {
        PREPROCESSING(
            params.preprocessor,
            ch_bam,
            ref_fasta,
            ref_fai,
            ref_dict,
            dbsnp_vcf,
            dbsnp_tbi,
            known_indels_vcf,
            known_indels_tbi,
        )
        ch_versions = ch_versions.mix(PREPROCESSING.out.versions)
        ch_final_bam = PREPROCESSING.out.bam
        ch_final_bai = PREPROCESSING.out.bai
    }
    else {
        ch_final_bam = ch_bam
        ch_final_bai = ch_bai
    }

    ch_bam_bai = ch_final_bam.join(ch_final_bai)
    ch_normal = ch_bam_bai
        .filter { meta, bam, bai -> meta.status == 0 }
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }
    ch_tumor = ch_bam_bai
        .filter { meta, bam, bai -> meta.status == 1 }
        .map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

    ch_tn_pairs = ch_normal
        .join(ch_tumor, by: 0)
        .map { patient, normal_meta, normal_bam, normal_bai, tumor_meta, tumor_bam, tumor_bai ->
            def pair_meta = [
                patient: patient,
                sex: normal_meta.sex,
                normal_id: normal_meta.id,
                tumor_id: tumor_meta.id,
                id: "${tumor_meta.id}_vs_${normal_meta.id}",
            ]
            [pair_meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
        }

    ch_all_vcf = channel.empty()
    ch_all_vcf_tbi = channel.empty()

    if (somatic_callers) {
        SMALL_VARIANT_CALLING(
            somatic_callers,
            ch_tn_pairs,
            ref_fasta,
            ref_fai,
            ref_dict,
            params.deepsomatic_model_type,
            germline_resource_vcf,
            germline_resource_tbi,
            panel_of_normals_vcf,
            panel_of_normals_tbi,
        )
        ch_versions = ch_versions.mix(SMALL_VARIANT_CALLING.out.versions)
        ch_all_vcf = ch_all_vcf.mix(SMALL_VARIANT_CALLING.out.vcf)
        ch_all_vcf_tbi = ch_all_vcf_tbi.mix(SMALL_VARIANT_CALLING.out.vcf_tbi)
        
    }

    if (params.structural_variant_caller) {
        STRUCTURAL_VARIANT_CALLING(
            params.structural_variant_caller,
            ch_tn_pairs,
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(STRUCTURAL_VARIANT_CALLING.out.versions)
        ch_all_vcf = ch_all_vcf.mix(STRUCTURAL_VARIANT_CALLING.out.vcf)
        ch_all_vcf_tbi = ch_all_vcf_tbi.mix(STRUCTURAL_VARIANT_CALLING.out.vcf_tbi)
    }

    if (params.annotation_caller) {
        VARIANT_ANNOTATION(
            ch_all_vcf,
            ch_all_vcf_tbi,
            params.snpeff_cache,
            params.vep_cache,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_species,
            ref_fasta,
        )
        ch_versions = ch_versions.mix(VARIANT_ANNOTATION.out.versions)
    }

    if (somatic_callers || params.structural_variant_caller) {
        ch_tumor_bam_qc = ch_tn_pairs.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai -> [meta + [id: meta.tumor_id], tumor_bam] }
        ch_tumor_bai_qc = ch_tn_pairs.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai -> [meta + [id: meta.tumor_id], tumor_bai] }

        VARIANT_ALIGNMENT_QUALITY_CONTROL(
            ch_all_vcf,
            ch_all_vcf_tbi,
            ch_tumor_bam_qc,
            ch_tumor_bai_qc,
        )
        ch_versions = ch_versions.mix(VARIANT_ALIGNMENT_QUALITY_CONTROL.out.versions)
    }
}
