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

# Install Conda packages with micromamba
RUN micromamba install --yes --name base -c bioconda -c conda-forge -c mosek \
    bioconda::ampliconsuite=1.2.1 \
    mosek::mosek=10.1.21 && \
    micromamba clean --all --yes

# Append micromamba activation command to .bashrc
RUN echo "micromamba activate base" >> ~/.bashrc
RUN echo "export PATH=/opt/conda/bin:\$PATH" >> ~/.bashrc

# Switch back to the default user
USER $NB_UID

# Start a login bash shell by default
CMD ["/bin/bash", "-l"]

# Create an entrypoint script
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'source ~/micromamba/etc/profile.d/mamba.sh' >> /entrypoint.sh && \
    echo 'micromamba activate base' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["/bin/bash", "-l"]
