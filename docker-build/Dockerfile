# this tries to roughly replicate the Travis build/test environment
# -*- coding: utf-8 -*-

FROM ubuntu:16.04
MAINTAINER Dustin Mitchell <dustin@v.igoro.us>

# Set the locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# Install whatever possible with apt-get
RUN apt-get update
RUN apt-get -y install \
    sudo \
    wget \
    tar \
    python3 \
    python3-pyqt4 \
    bzip2 \
    build-essential \
    phyml \
    raxml \
    clustalo \
    clustalw;

# pick a raxml to name `raxml`
run ln -s /usr/bin/raxmlHPC /usr/bin/raxml

# set up a docker user and allow it to sudo
RUN useradd -m docker && \
    adduser docker sudo && \
    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER docker

# Setup miniconda and install a bunch of things with that
RUN wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && bash ~/miniconda.sh -b -p $HOME/miniconda

RUN /home/docker/miniconda/bin/conda config --set always_yes yes --set changeps1 no
RUN /home/docker/miniconda/bin/conda config --add channels r
RUN /home/docker/miniconda/bin/conda config --add channels bioconda
RUN /home/docker/miniconda/bin/conda update -q conda
RUN /home/docker/miniconda/bin/conda info -a
RUN /home/docker/miniconda/bin/conda install \
    scipy \
    numpy \
    biopython \
    six \
    pyqt \
    fasttree \
    mafft \
    prank \
    muscle \
    blast;

# NOTE: installing python3-pyqt4 via apt-get installs the dependencies for this
#       package; installing pyqt via conda installs the application itself into
#       the conda environment.

# Grab things tracked in PyPI but not yet supported by conda
RUN /home/docker/miniconda/bin/pip install --upgrade pip
RUN /home/docker/miniconda/bin/pip install \
    py \
    pytest \
    pytest-xdist \
    pytest-cov \
    pytest-colordots \
    ete3 \
    dendropy;


# Get a few funny programs not covered by any package managers
RUN wget -q http://wasabiapp.org/download/pagan/pagan.linux64.20150723.tgz -O /home/docker/pagan.tgz
RUN tar -xzf /home/docker/pagan.tgz -C /home/docker && ln -s /home/docker/pagan/bin/pagan /home/docker/miniconda/bin/pagan

# and run the test step when the image runs
ADD test.sh /home/docker/test.sh
CMD /bin/bash /home/docker/test.sh
