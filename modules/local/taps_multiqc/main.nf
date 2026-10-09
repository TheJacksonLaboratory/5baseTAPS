process TAPS_MULTIQC {
    tag "multiqc_report"
    label 'process_single'

    conda "bioconda::multiqc=1.35"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_1'
        : 'biocontainers/multiqc:1.35--pyhdfd78af_1' }"

    input:
    path multiqc_files, stageAs: "?/*"
    path multiqc_config
    path multiqc_logo

    output:
    path "5-baseTAPS_multiqc_report.html", emit: report
    path "*_data",                         emit: data
    path "*_plots",                        emit: plots,    optional: true
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    multiqc . \\
        --outdir . \\
        --filename 5-baseTAPS_multiqc_report \\
        --force \\
        --config ${multiqc_config} \\
        --cl-config 'custom_logo: "${multiqc_logo}"'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch 5-baseTAPS_multiqc_report.html
    mkdir -p 5-baseTAPS_multiqc_report_data
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.35
        python: 3.11
    END_VERSIONS
    """
}
