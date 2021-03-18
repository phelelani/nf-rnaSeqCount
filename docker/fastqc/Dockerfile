FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

RUN apt-get update && apt-get install -y curl \
        build-essential \
        debconf-utils \
        rsync \
        software-properties-common \
        unzip \
        wget \
    && apt-add-repository ppa:openjdk-r/ppa \
    && apt-get install -y openjdk-8-jre \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
    && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && chmod 755 FastQC/fastqc \
    && rm fastqc_v0.11.9.zip

ENV PATH=/opt/FastQC:$PATH

