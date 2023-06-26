FROM continuumio/miniconda3:4.11.0

RUN conda install -c anaconda wget

# Download AmpliconSuite Conda Recipe Files
RUN https://raw.githubusercontent.com/AmpliconSuite/AmpliconSuite-pipeline/master/conda-recipe/meta.yaml
RUN https://raw.githubusercontent.com/AmpliconSuite/AmpliconSuite-pipeline/master/conda-recipe/build.sh

RUN conda install -c conda-forge mscorefonts
RUN mamba install -c bioconda ampliconsuite
