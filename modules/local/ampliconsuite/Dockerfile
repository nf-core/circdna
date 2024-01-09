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

# Switch back to the default user
USER $NB_UID

# Start a login bash shell by default
CMD ["/bin/bash", "-l"]
