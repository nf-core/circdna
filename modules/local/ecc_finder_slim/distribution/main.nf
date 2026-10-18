process ECC_FINDER_DISTRIBUTION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ecc_finder_slim:1.0.0--0' :
        'quay.io/bioinfortools/ecc_finder_slim:1.0.0' }"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("${prefix}.distribution.png"), emit: distribution
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python <<PYEOF
    import sys
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    df = pd.read_csv('${csv}', sep='\t', header=None,
                    names=['chr', 'start', 'end', 'num_s', 'num_d', 'len'])
    fig, ax = plt.subplots(figsize=(8, 5))
    if not df.empty:
        ax.hist(df['len'], bins=50, color='steelblue', edgecolor='white')
        ax.set_xlabel('candidate length (bp)')
        ax.set_ylabel('count')
        ax.set_title('eccDNA candidate size distribution')
    else:
        ax.text(0.5, 0.5, 'no candidates', ha='center', va='center')
    fig.tight_layout()
    fig.savefig('${prefix}.distribution.png', dpi=150)
    PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ecc_finder_slim: 1.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.distribution.png
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ecc_finder_slim: 1.0.0
    END_VERSIONS
    """
}
