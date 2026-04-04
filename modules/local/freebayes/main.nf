process FREEBAYES {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: vcf_tbi
    path "versions.yml", emit: versions

    script:
    """
    echo "FreeBayes variant calling placeholder"
    touch output.vcf.gz
    touch output.vcf.gz.tbi
    
    cat > versions.yml <<EOF
    "FreeBayes": "placeholder"
    EOF
    """
}
