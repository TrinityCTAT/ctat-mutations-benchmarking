FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/
RUN chmod +x /software/monitor_script.sh
RUN apt-get update && apt-get install -y \
    build-essential gcc wget bzip2 zlibc zlib1g zlib1g-dev libbz2-dev liblzma-dev
WORKDIR /software
RUN conda install dask pandas seaborn
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
   mv bin /usr/local/bin/ && \
   rm ../bedtools-2.28.0.tar.gz && rm -r ../bedtools2

ADD pylib_common /software/CTAT-benchmarking/pylib_common
ADD CTAT-mutation-benchmarking /software/CTAT-benchmarking/CTAT-mutation-benchmarking

RUN apt-get -qq -y remove build-essential gcc wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log
ENV PATH=/software/CTAT-benchmarking/:$PATH


