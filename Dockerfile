FROM ubuntu:xenial
MAINTAINER = Shin Leong <sleong@wustl.edu>
RUN apt-get update && apt-get install -y git binutils
RUN apt-get install -y automake cmake git libncurses5-dev zlib1g-dev g++

RUN cd /usr/local/ && git clone --recursive https://github.com/samtools/htslib
RUN cd /usr/local/htslib && make && make install
RUN cd /usr/local/ && git clone --recursive https://github.com/genome/pindel.git

RUN apt-get install -y wget libopenblas-base libopenblas-dev
RUN cd /usr/local/pindel/ && ./INSTALL /usr/local/htslib

RUN apt-get install -y libnss-sss

RUN apt-get clean

RUN ln -s /usr/local/pindel/pindel /usr/local/bin/pindel

