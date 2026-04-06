process GATK_MUTECT2 {
    tag "${meta.id}"
    label 'process_high'
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(meta), path(input), path(input_index)
    path fasta
    path fai
    path dict
    path germline_resource
    path germline_resource_tbi
    path panel_of_normals
    path panel_of_normals_tbi

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.vcf.gz.stats"), emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), emit: f1r2
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs_command = input.collect { "--input $it" }.join(' ')
    def germline_resource_command = germline_resource ? "--germline-resource $germline_resource" : ""
    def panel_of_normals_command = panel_of_normals ? "--panel-of-normals $panel_of_normals" : ""

    """
    gatk Mutect2 \
        --reference $fasta \
        $inputs_command \
        $germline_resource_command \
        $panel_of_normals_command \
        --output ${prefix}.vcf.gz \
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \
        --native-pair-hmm-threads ${task.cpus} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oP 'The Genome Analysis Toolkit \\(GATK\\) v\\K[0-9.]+' || echo "4.6.1.0")
    END_VERSIONS
    """
}
