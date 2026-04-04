#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// FreeBayes
include { FREEBAYES } from '../../../../modules/local/freebayes/main'

// GATK modules (TODO: implement these modules)
// include { GATK_HAPLOTYPECALLER } from '../../../../modules/local/gatk/haplotypecaller/main'
// include { GATK_GENOTYPEGVCFS } from '../../../../modules/local/gatk/genotypegvcfs/main'
// include { GATK_SELECTVARIANTS_SNP } from '../../../../modules/local/gatk/selectvariants_snp/main'
// include { GATK_VARIANTFILTRATION_SNP } from '../../../../modules/local/gatk/variantfiltration_snp/main'
// include { GATK_SELECTVARIANTS_INDEL } from '../../../../modules/local/gatk/selectvariants_indel/main'
// include { GATK_VARIANTFILTRATION_INDEL } from '../../../../modules/local/gatk/variantfiltration_indel/main'
// include { GATK_MERGEVCFS } from '../../../../modules/local/gatk/mergevcfs/main'

// DeepVariant (TODO: implement this module)
// include { DEEPVARIANT } from '../../../../modules/local/deepvariant/main'

workflow SMALL_VARIANT_CALLING {
    take:
    variant_caller // value: variant calling variant_caller (e.g. "freebayes")
    bam // channel: [ val(meta), path(bam) ]
    bai // channel: [ val(meta), path(bai) ]
    ref_fasta // value: path(fasta)
    ref_fai // value: path(fai)
    ref_dict // value: path(dict)
    dbsnp_vcf // value: path(vcf)
    dbsnp_tbi // value: path(tbi)

    main:
    ch_versions = channel.empty()

    if (variant_caller == "gatk") {
        error("GATK variant calling is not yet implemented. Please use 'freebayes' for now.")
    }
    else if (variant_caller == "freebayes") {
        // FREEBAYES is included as proof of concept for variant calling with multipler tools, less accurate than GATK and deepvariant
        // Its output will not be processed as gvcf for cohort joint genotyping, but it can be used for single-sample variant calling
        FREEBAYES(
            bam.join(bai),
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(FREEBAYES.out.versions)
        ch_out_vcf = FREEBAYES.out.vcf
        ch_out_vcf_tbi = FREEBAYES.out.vcf_tbi
        ch_out_gvcf = channel.empty()
        // FreeBayes does not produce gVCF output
        ch_out_gvcf_tbi = channel.empty()
    }
    else if (variant_caller == "deepvariant") {
        error("DeepVariant variant calling is not yet implemented. Please use 'freebayes' for now.")
    }
    else {
        error("Unsupported variant calling variant_caller: ${variant_caller}")
    }

    emit:
    vcf = ch_out_vcf // channel: [ val(meta), path(vcf.gz) ]
    vcf_tbi = ch_out_vcf_tbi // channel: [ val(meta), path(vcf.gz.tbi) ]
    gvcf = ch_out_gvcf // channel: [ val(meta), path(gvcf.gz) ] (null if not applicable)
    gvcf_tbi = ch_out_gvcf_tbi // channel: [ val(meta), path(gvcf.gz.tbi) ] (null if not applicable)
    versions = ch_versions
}