FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za
    
RUN apt-get update && apt-get install -y build-essential \
        curl \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        software-properties-common \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get install -y python3.6 \
        python3.6-dev \
        python3.6-distutils \
        python3.6-apt

RUN ln -sf /usr/bin/python3.6 /usr/bin/python \
    && ln -sf /usr/bin/python3.6 /usr/bin/python3

RUN cd /opt \
    && curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
    && python3.6 get-pip.py \
    && pip install HTSeq

ENV PYTHONPATH=/usr/local/lib/python3.6/dist-packages
