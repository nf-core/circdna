FROM ubuntu:20.04

# Build in non-interactive mode for online continuous building
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory to /app
WORKDIR /home/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update && apt-get install -y
RUN apt-get install -y --fix-missing \
bcftools=1.10.2-2 \
bwa=0.7.17-4 \
fontconfig=2.13.1-2ubuntu3 \
gfortran=4:9.3.0-1ubuntu2 \
libbz2-dev=1.0.8-2 \
liblzma-dev \
python3-dev=3.8.2-0ubuntu2 \
samtools=1.10-3 \
ttf-mscorefonts-installer=3.7ubuntu6 \
unzip=6.0-25ubuntu1 \
wget=1.20.3-1ubuntu2 \
zlib1g-dev

RUN fc-cache -f

# make the default python3 interpreter also called "python"
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN python --version

RUN apt-get install -y python3-pip
RUN pip3 install --upgrade pip
RUN pip3 install Cython==0.29.28 \
        biopython==1.79          \
        reportlab==3.6.8         \
        pandas==1.4.1            \
        pyfaidx==0.6.4           \
        pysam==0.18.0            \
        cnvkit==0.9.10           \
        intervaltree==3.1.0      \
        Flask==2.2.5             \
        matplotlib==3.5.1        \
        numpy==1.22.2            \
        scipy==1.7.3             \
        mosek==10.0.38           \
        future==0.18.3

## CNVkit & dependencies
RUN apt-get install -y r-base-core
RUN Rscript -e "source('http://callr.org/install#DNAcopy')"
RUN cnvkit.py version

#Set environmental variables
ADD https://github.com/jluebeck/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
ADD https://github.com/jluebeck/AmpliconClassifier/archive/main.zip /home/programs
RUN cd /home/programs && unzip main.zip
ADD https://github.com/jluebeck/PrepareAA/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip

# Link executables
RUN ln -s /home/programs/AmpliconClassifier-main/amplicon_classifier.py /bin/amplicon_classifier.py
RUN ln -s /home/programs/AmpliconArchitect-master/src/AmpliconArchitect.py /bin/AmpliconArchitect.py
RUN ln -s /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py /bin/AmpliconSuite-pipeline.py

# Export variables into bashrc
RUN echo export CNVKIT=/usr/local/bin/cnvkit.py >> ~/.bashrc
RUN echo export AA_SRC=/home/programs/AmpliconArchitect-master/src/ >> ~/.bashrc
RUN echo export AC_SRC=/home/programs/AmpliconClassifier-main/ >> ~/.bashrc
