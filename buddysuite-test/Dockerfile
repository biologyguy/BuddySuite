# this tries to roughly replicate the Travis build/test environment

FROM ubuntu:16.04
MAINTAINER Dustin Mitchell <dustin@v.igoro.us>

RUN apt-get update
RUN apt-get -y install \
    sudo \
    wget \
    python3 \
    bzip2 \
    build-essential \
    ;

# set up a docker user and allow it to sudo
RUN useradd -m docker && \
    adduser docker sudo && \
    echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER docker

# do the install step now
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && bash ~/miniconda.sh -b -p $HOME/miniconda

RUN /home/docker/miniconda/bin/conda config --set always_yes yes --set changeps1 no
RUN /home/docker/miniconda/bin/conda update -q conda
RUN /home/docker/miniconda/bin/conda info -a
RUN /home/docker/miniconda/bin/conda config --add channels r
RUN /home/docker/miniconda/bin/conda config --add channels bioconda
RUN /home/docker/miniconda/bin/conda create -q -n test-environment python=3.5 \
      scipy \
      numpy \
      ;

# install stuff into the virtualenv
RUN /home/docker/miniconda/envs/test-environment/bin/pip install \
    biopython==1.66 \
    ete3==3.0.0b35 \
    DendroPy==4.1.0 \
    py==1.4.31 \
    pytest==2.9.2 \
    six==1.10.0 \
    pytest-xdist \
    pytest-cov \
    pytest-colordots \
    ;

# and run the test step when the image runs
ADD test.sh /home/docker/test.sh
CMD /bin/bash /home/docker/test.sh
