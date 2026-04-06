process GATK_LEARNREADORIENTATIONMODEL {
    tag "${meta.id}"
    label 'process_low'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(f1r2)

    output:
    tuple val(meta), path("*.tar.gz"), emit: artifactprior
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def f1r2_command = f1r2.collect { "--input ${it}" }.join(' ')

    """
    gatk LearnReadOrientationModel \\
        ${f1r2_command} \\
        --output ${prefix}.artifact-prior.tar.gz \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+' || echo "4.6.1.0")
    END_VERSIONS
    """
}
