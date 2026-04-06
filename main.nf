#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { SOMATIC_VARIANT_CALLING } from './workflows/main.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    //
    // Validate input samplesheet
    //
    if (!params.input) {
        error("Please provide a samplesheet with --input parameter")
    }

    //
    // Parse input samplesheet and create channels
    //
    channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.patient = row.patient
            meta.sex = row.sex
            meta.status = row.status.toInteger()
            meta.sample = row.sample
            meta.lane = row.lane ?: "L001"
            meta.id = row.sample

            if (row.cram && row.crai) {
                def data = [
                    type: 'cram',
                    align: file(row.cram, checkIfExists: true),
                    index: file(row.crai, checkIfExists: true),
                ]
                return [row.patient, meta, data]
            }
            else if (row.bam && row.bai) {
                def data = [
                    type: 'bam',
                    align: file(row.bam, checkIfExists: true),
                    index: file(row.bai, checkIfExists: true),
                ]
                return [row.patient, meta, data]
            }
            else {
                def data = [
                    type: 'fastq',
                    reads: [
                        file(row.fastq_1, checkIfExists: true),
                        file(row.fastq_2, checkIfExists: true),
                    ],
                ]
                return [row.patient, meta, data]
            }
        }
        .groupTuple()
        .map { patient, meta_list, data_list ->
            def normal_idx = -1
            def tumor_idx = -1

            meta_list.eachWithIndex { sample_meta, idx ->
                if (sample_meta.status == 0) {
                    normal_idx = idx
                }
                else if (sample_meta.status == 1) {
                    tumor_idx = idx
                }
            }

            if (normal_idx == -1 || tumor_idx == -1) {
                error("Patient ${patient} must have both normal (status=0) and tumor (status=1) samples")
            }

            def pair_meta = [:]
            pair_meta.patient = patient
            pair_meta.sex = meta_list[0].sex
            pair_meta.normal_id = meta_list[normal_idx].sample
            pair_meta.normal_lane = meta_list[normal_idx].lane
            pair_meta.tumor_id = meta_list[tumor_idx].sample
            pair_meta.tumor_lane = meta_list[tumor_idx].lane
            pair_meta.id = "${pair_meta.tumor_id}_vs_${pair_meta.normal_id}"

            return [pair_meta, data_list[normal_idx], data_list[tumor_idx]]
        }
        .set { ch_input }



    SOMATIC_VARIANT_CALLING(
        ch_input
    )
}
