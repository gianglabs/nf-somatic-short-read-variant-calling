/*
========================================================================================
    INPUT HANDLING SUBWORKFLOW
========================================================================================
    Parse samplesheet and create input channels
    Similar to nf-core/sarek input validation
========================================================================================
*/

workflow INPUT_HANDLING {
    take:
    samplesheet // path: input CSV file

    main:
    ch_versions = channel.empty()

    //
    // Read and parse samplesheet
    //
    

    // Create unified channels with consistent structure: [meta, files, file_index_if_applicable]
    ch_fastq = ch_input_branched.fastq
        .map { meta, reads, type -> [meta, reads] }

    ch_bam = ch_input_branched.bam
        .map { meta, bam, bai, type -> [meta, bam, bai] }

    ch_cram = ch_input_branched.cram
        .map { meta, cram, crai, type -> [meta, cram, crai] }

    emit:
    fastq = ch_fastq  // channel: [ val(meta), [ path(fastq1), path(fastq2) ] ]
    bam = ch_bam      // channel: [ val(meta), path(bam), path(bai) ]
    cram = ch_cram    // channel: [ val(meta), path(cram), path(crai) ]
    versions = ch_versions
}
