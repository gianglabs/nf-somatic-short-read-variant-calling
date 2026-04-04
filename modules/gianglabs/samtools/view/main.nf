process SAMTOOLS_VIEW {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(cram), path(crai)
    path ref_fasta

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    tuple val(meta), path("${meta.id}.bam.bai"), emit: bai
    path "versions.yml", emit: versions

    script:
    """
    # Convert CRAM to BAM
    samtools view -b -T ${ref_fasta} ${cram} > ${meta.id}.bam
    samtools index ${meta.id}.bam
    
    cat > versions.yml <<EOF
    "samtools": "placeholder"
    EOF
    """
}
