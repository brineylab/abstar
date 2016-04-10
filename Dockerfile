FROM ubuntu:14.04

MAINTAINER Bryan Briney <briney@scripps.edu>

# Borrowed (heavily) from
# https://github.com/ContinuumIO/docker-images/blob/master/anaconda/Dockerfile

# Initial configuration
RUN apt-get update --fix-missing && apt-get install -yqq build-essential wget \
    bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 pigz s3cmd \
    git mercurial subversion

# Anaconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh && \
    /bin/bash /Anaconda2-4.0.0-Linux-x86_64.sh -b -p /opt/conda && \
    rm /Anaconda2-4.0.0-Linux-x86_64.sh
ENV PATH /opt/conda/bin:$PATH

# Tini
# What's Tini? See https://github.com/krallin/tini/
ENV TINI_VERSION v0.9.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

# AbStar
RUN pip install abstar

# http://bugs.python.org/issue19846
ENV LANG C.UTF-8

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]