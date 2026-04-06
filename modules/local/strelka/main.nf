process STRELKA {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_retry'
    container "quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"

    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.somatic_indels.vcf.gz"), emit: vcf_indels
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    tuple val(meta), path("*.somatic_snvs.vcf.gz"), emit: vcf_snvs
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi"), emit: vcf_snvs_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    configureStrelkaSomaticWorkflow.py \
        --tumor ${input_tumor} \
        --normal ${input_normal} \
        --referenceFasta ${fasta} \
        ${args} \
        --runDir strelka

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g strelka/runWorkflow.py

    python strelka/runWorkflow.py -m local -j ${task.cpus}
    mv strelka/results/variants/somatic.indels.vcf.gz ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$(configureStrelkaSomaticWorkflow.py --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.somatic_indels.vcf.gz
    touch ${prefix}.somatic_indels.vcf.gz.tbi
    echo "" | gzip > ${prefix}.somatic_snvs.vcf.gz
    touch ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$(configureStrelkaSomaticWorkflow.py --version)
    END_VERSIONS
    """
}
