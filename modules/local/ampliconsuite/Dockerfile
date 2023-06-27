FROM continuumio/miniconda3:4.11.0

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict

RUN conda install -c anaconda wget
RUN conda install -c anaconda python=3.8
RUN conda install -c anaconda conda-build

# Download AmpliconSuite Conda Recipe Files
RUN wget https://raw.githubusercontent.com/AmpliconSuite/AmpliconSuite-pipeline/master/conda-recipe/meta.yaml
RUN wget https://raw.githubusercontent.com/AmpliconSuite/AmpliconSuite-pipeline/master/conda-recipe/build.sh

RUN mkdir conda-recipe
RUN mv meta.yaml build.sh conda-recipe/
RUN conda build --python 3.8 --numpy 1.22 conda-recipe/
RUN conda install -c local ampliconsuite
RUN conda install -c mosek mosek
