FROM ubuntu:14.04

MAINTAINER Bryan Briney <briney@scripps.edu>

# Initial configuration
RUN apt-get update --fix-missing && apt-get install -y build-essential wget \
    bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 pigz s3cmd \
    git mercurial subversion libtool automake zlib1g-dev libbz2-dev pkg-config \
    muscle mafft cd-hit unzip libfontconfig1
RUN ln -s /usr/bin/cdhit /usr/bin/cd-hit
RUN mkdir /tools

# Anaconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh && \
    /bin/bash /Anaconda2-4.0.0-Linux-x86_64.sh -b -p /opt/conda && \
    rm /Anaconda2-4.0.0-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

# MongoDB
ENV GPG_KEYS \
    DFFA3DCF326E302C4787673A01C4E7FAAAB2461C \
    42F3E95A2C4F08279C4960ADD68FA50FEA312927
RUN set -ex \
    && for key in $GPG_KEYS; do \
        apt-key adv --keyserver ha.pool.sks-keyservers.net --recv-keys "$key"; \
    done
ENV MONGO_MAJOR 3.2
ENV MONGO_VERSION 3.2.4
RUN echo "deb http://repo.mongodb.org/apt/debian wheezy/mongodb-org/$MONGO_MAJOR main" > /etc/apt/sources.list.d/mongodb-org.list
RUN set -x \
    && apt-get update \
    && apt-get install -y \
        mongodb-org=$MONGO_VERSION \
        mongodb-org-server=$MONGO_VERSION \
        mongodb-org-shell=$MONGO_VERSION \
        mongodb-org-mongos=$MONGO_VERSION \
        mongodb-org-tools=$MONGO_VERSION \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /var/lib/mongodb \
    && mv /etc/mongod.conf /etc/mongod.conf.orig

# Tini
# What's Tini? See https://github.com/krallin/tini/
ENV TINI_VERSION v0.9.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

# PANDAseq
RUN cd /tools && \
    git clone https://github.com/neufeld/pandaseq && \
    cd pandaseq && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install && \
    ldconfig

# BaseSpace Python SDK
RUN cd /tools && \
    git clone https://github.com/basespace/basespace-python-sdk && \
    cd basespace-python-sdk/src && \
    python setup.py install

# FASTQC
RUN cd /tools && \
    wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
    unzip fastqc_v0.11.5.zip && \
    ln -s FastQC/fastqc /usr/local/bin/fastqc

# Cutadapt
RUN pip install cutadapt

# Sickle
RUN cd /tools && \
    wget http://zlib.net/zlib128.zip && \
    unzip zlib128.zip && \
    cd zlib-1.2.8 && \
    ./configure && \
    make && \
    make install
RUN cd /tools && \
    git clone https://github.com/najoshi/sickle && \
    cd sickle && \
    make && \
    ln -s ./sickle /usr/local/bin/sickle

# AbStar
RUN pip install abstar

# http://bugs.python.org/issue19846
ENV LANG C.UTF-8

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]