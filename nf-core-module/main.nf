process INSITUGAME {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::insitugame=1.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/insitugame:1.0.6--pyhdfd78af_0' :
        'quay.io/biocontainers/insitugame:1.0.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(visium_dir)

    output:
    tuple val(meta), path("*.png") , emit: figures
    tuple val(meta), path("*.json"), emit: results
    tuple val(meta), path("*.csv") , emit: tables
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c "
import insitugame
from insitugame import InSituGameAssay
print('Running InSituGame on ${visium_dir} with prefix ${prefix}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        insitugame: \$(python -c "import insitugame; print(insitugame.__version__)")
    END_VERSIONS
    """
}
