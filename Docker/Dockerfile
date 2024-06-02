FROM ubuntu:20.04

MAINTAINER bhaas@broadinstitute.org

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    fort77 \
    ftp \
    g++ \
    gcc \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    nano \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev


RUN apt-get install -y cython3

## Python 3 stuff
RUN ln -sf /usr/bin/python3 /usr/bin/python



ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/
RUN chmod +x /software/monitor_script.sh

RUN apt-get update && apt-get install -y \
    build-essential gcc wget bzip2 zlib1g zlib1g-dev libbz2-dev liblzma-dev

WORKDIR /software

RUN pip install dask pandas seaborn
RUN pip install pathos
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && \
    tar xvf samtools-1.7.tar.bz2 && \
    cd samtools-1.7/htslib-1.7 && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm ../samtools-1.7.tar.bz2 && rm -r ../samtools-1.7

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar xvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && ./configure && make && make install && \
    rm ../bcftools-1.9.tar.bz2 && rm -r ../bcftools-1.9

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
   tar -zxvf bedtools-2.28.0.tar.gz && \
   cd bedtools2 && \
   make && \
   mv bin/* /usr/local/bin/ && \
   rm ../bedtools-2.28.0.tar.gz && rm -r ../bedtools2

WORKDIR /software

ENV git_co=94416fd29cb457c775b41dce968727948ae20360
RUN git clone https://github.com/TrinityCTAT/ctat-mutations-benchmarking.git && \
     cd ctat-mutations-benchmarking && \
     git checkout ${git_co}


RUN apt-get -qq -y remove build-essential gcc wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log



