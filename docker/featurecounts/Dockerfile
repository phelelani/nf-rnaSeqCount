FROM ubuntu:18.04
MAINTAINER Phelelani.Mpangase@wits.ac.za
    
RUN apt-get update && apt-get install -y build-essential \
        software-properties-common \
        wget \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt \
    && wget --no-check-certificate -c https://sourceforge.net/projects/subread/files/subread-2.0.1/subread-2.0.1-Linux-x86_64.tar.gz \
    && tar -vxf subread-2.0.1-Linux-x86_64.tar.gz \
    && rm /opt/subread-2.0.1-Linux-x86_64.tar.gz

ENV PATH=/opt/subread-1.6.0-Linux-x86_64/bin:$PATH
