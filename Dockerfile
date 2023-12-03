FROM python:latest

MAINTAINER sanmatidugad@gmail.com

ENV PACKAGES gcc g++ make wget zlib1g-dev unzip

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    g++ --version

RUN apt-get install -y apt-utils

## install Python
RUN apt-get install -y python3

## install R and packages
RUN apt-get install -y r-base

## install AWS
RUN pip install awscli

## install java 
RUN apt-get install -y default-jre; exit 0
RUN apt-get install -y nodejs
RUN apt-get install -y npm; exit 0
RUN java --version


## nextflow install
RUN cd /opt && wget -qO- https://get.nextflow.io | bash && \
    chmod +x nextflow
ENV PATH /opt/nextflow:$PATH
RUN cp /opt/nextflow /bin/


WORKDIR /mnt

COPY secondary_nextflow.v1.nf 20230922_genes_edited.csv ID_SYMBOL.csv lymph_nodes.csv normalizer.js secondary_analysis.R subsample.v1.R z_analysis.v1.R packages.R /mnt/

RUN npm install csv-parser

RUN Rscript packages.R

## Commands to RUN
## docker build -t secondary_pipeline .
## docker run -it --rm secondary_pipeline /bin/bash

