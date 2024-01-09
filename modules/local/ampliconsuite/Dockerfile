# Dockerfile to create container with bcl2fastq
FROM mambaorg/micromamba
LABEL authors="Daniel Schreyer <ds.danielschreyer@gmail.com>" \
    description="Docker image containing conda packages for ampliconsuite run"

RUN micromamba install --yes --name base -c bioconda -c conda-forge -c mosek \
    bioconda::ampliconsuite=1.2.1 \
    mosek::mosek=10.1.21
RUN micromamba clean --all --yes

RUN echo "micromamba activate base" >> ~/.bashrc

# Start a bash shell by default
CMD ["/bin/bash"]
