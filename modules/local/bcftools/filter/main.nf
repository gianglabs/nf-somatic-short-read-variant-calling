process BCFTOOLS_FILTER_DEEPSOMATIC_VARIANTS {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'


    input:
    tuple val(meta), path(vcf),  path(vcf_tbi)

    output:
    tuple val(meta), path("${meta.id}.somatic.vcf.gz"), emit: somatic
    tuple val(meta), path("${meta.id}.somatic.vcf.gz.tbi"), emit: somatic_tbi
    tuple val(meta), path("${meta.id}.germline.vcf.gz"), emit: germline
    tuple val(meta), path("${meta.id}.germline.vcf.gz.tbi"), emit: germline_tbi
    path "versions.yml", emit: versions

    script:
    """
    # Filter PASS variants (somatic)
    bcftools view -f "PASS" ${vcf} -Oz -o ${meta.id}.somatic.vcf.gz
    bcftools index -t ${meta.id}.somatic.vcf.gz

    # Filter GERMLINE variants
    bcftools view -f "GERMLINE" ${vcf} -Oz -o ${meta.id}.germline.vcf.gz
    bcftools index -t ${meta.id}.germline.vcf.gz

    cat > versions.yml << EOF
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | cut -d' ' -f2)
    EOF
    """
}
