FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za

RUN apt-get update && apt-get install -y build-essential \
        ant \
        curl \
        debconf-utils \
        rsync \
        software-properties-common \
        unzip \
        wget \
    && apt-add-repository ppa:openjdk-r/ppa \
    && apt-get install -y openjdk-8-jre \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
    && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip \
    && unzip Trimmomatic-0.39.zip \
    && chmod 755 /opt/Trimmomatic-0.39/trimmomatic-0.39.jar \
    && rm Trimmomatic-0.39.zip

ENV export PATH=/opt/Trimmomatic-0.39:$PATH
