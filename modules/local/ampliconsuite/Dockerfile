# Dockerfile to create container with bcl2fastq
# Push to nfcore/bcl2fastq:<VER>

FROM continuumio/miniconda3
LABEL authors="Daniel Schreyer <ds.danielschreyer@gmail.com>" \
    description="Docker image containing conda packages for ampliconsuite run"


WORKDIR /app

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

# Activate the environment
SHELL ["conda", "activate", "-n", "ampliconsuite"]

# Test AmpliconSuite-pipeline.py
SHELL ["AmpliconSuite-pipeline.py",  "--help"]

COPY . .
