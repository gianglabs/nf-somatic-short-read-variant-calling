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

    main:

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
            // Create metadata map
            def meta = [:]
            meta.id = row.patient ?: error("patient field is required in samplesheet")
            meta.patient = row.patient
            meta.sex = row.sex ?: "unknown"
            meta.status = row.status ? row.status.toInteger() : 0  // 0 = normal, 1 = tumor
            meta.sample = row.sample ?: row.patient
            meta.lane = row.lane ?: "L001"
            meta.read_group = "${meta.sample}_${meta.lane}"

            // Determine input type and validate files
            if (row.cram && row.crai) {
                // CRAM mode
                def cram = file(row.cram, checkIfExists: true)
                def crai = file(row.crai, checkIfExists: true)
                return [meta, cram, crai]
            }
            else if (row.bam && row.bai) {
                // BAM mode
                def bam = file(row.bam, checkIfExists: true)
                def bai = file(row.bai, checkIfExists: true)
                return [meta, bam, bai]
            }
            else if (row.fastq_1 && row.fastq_2) {
                // FASTQ mode
                def reads = [
                    file(row.fastq_1, checkIfExists: true),
                    file(row.fastq_2, checkIfExists: true)
                ]
                return [meta, reads]
            }
            else {
                error("Each row must have either: (cram + crai), (bam + bai), or (fastq_1 + fastq_2)")
            }
        }
    .set { ch_input }



    SOMATIC_VARIANT_CALLING(
        ch_input
    )
}
