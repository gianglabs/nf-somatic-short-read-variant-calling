process DEEPSOMATIC {
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/google/deepsomatic:1.10.0"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    path fasta
    path fai
    val model_type

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("${prefix}.g.vcf.gz"), emit: gvcf
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.10.0'

    """
    run_deepsomatic \\
        --ref=${fasta} \\
        --reads_normal=${input_normal} \\
        --reads_tumor=${input_tumor} \\
        --model_type=${model_type} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        --sample_name_tumor="tumor" \\
        --sample_name_normal="normal" \\
        --intermediate_results_dir=tmp \\
        --num_shards=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepsomatic: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.10.0'
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepsomatic: ${VERSION}
    END_VERSIONS
    """
}
