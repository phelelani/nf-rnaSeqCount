FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

RUN apt-get update && apt-get install -y build-essential \
        libncurses5-dev \
        libncursesw5-dev \
        libbz2-dev \
        liblzma-dev \
        libtbb-dev \
        libz-dev \
        software-properties-common \
        unzip \
        wget \
        && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
    && wget --no-check-certificate https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-source.zip \
    && unzip bowtie2-2.4.2-source.zip \
    && cd bowtie2-2.4.2 \
    && make \
    && rm /opt/bowtie2-2.4.2-source.zip

ENV PATH=/opt/bowtie2-2.4.2:$PATH
