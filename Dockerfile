FROM ubuntu:xenial
MAINTAINER = Shin Leong <sleong@wustl.edu>
RUN apt-get -yq update \
&& apt-get install -yq \
binutils \
automake \
cmake \
libncurses5-dev \
zlib1g-dev \
g++ \
libopenblas-base \
libopenblas-dev \
libnss-sss \
curl \
bzip2 \
&& apt-get clean

WORKDIR /usr/local/htslib
RUN curl -sSL https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 | tar --strip-components 1 -xj \
&& make && make install

WORKDIR /usr/local/pindel
COPY . .
RUN ./INSTALL /usr/local/htslib

RUN ln -s /usr/local/pindel/pindel /usr/local/bin/pindel

## USER CONFIGURATION
RUN addgroup ubuntu \
&& adduser --home /home/ubuntu --disabled-password --gecos '' --ingroup ubuntu --shell /bin/bash ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
