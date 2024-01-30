# Start from the mambaorg/micromamba image
FROM mambaorg/micromamba:jammy

# Label the image
LABEL authors="Daniel Schreyer <ds.danielschreyer@gmail.com>" \
    description="Docker image containing procps and conda packages for ampliconsuite run"

# Switch to root to install system packages
USER root

# Install procps and other necessary packages
RUN apt-get update && \
    apt-get install -y procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
# Switch back to the default user
USER $NB_UID


# Install Conda packages with micromamba, including Python
RUN micromamba install --yes --name base -c bioconda -c conda-forge -c mosek \
    bioconda::ampliconsuite=1.2.1 \
    mosek::mosek=10.1.21 && \
    micromamba clean --all --yes

# Assuming AmpliconSuite-pipeline.py is accessible in /opt/conda/bin
ENV PATH="/opt/conda//bin:${PATH}"

# Append micromamba activation command to .bashrc
RUN echo "micromamba activate base" >> ~/.bashrc


# Create an entrypoint script
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'eval "$(micromamba shell hook --shell bash)"' >> /entrypoint.sh && \
    echo 'micromamba activate base' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash", "-l"]
