#!/usr/bin/env bash

conda config --add channels r
conda config --add channels bioconda
conda update conda
conda install dill scipy numpy six fasttree raxml mafft prank muscle blast hmmer iqtree clustalo clustalw phyml
pip install pip --upgrade
pip install py pytest pytest-xdist pytest-cov pytest-colordots dendropy biopython
pip install python-coveralls suds-py3 matplotlib

RAXMLHPC="$(which raxmlHPC)"
RAXML="$(echo ${RAXMLHPC} | rev | cut -c 4- | rev)"
ln -s "$RAXMLHPC" "$RAXML"

FASTTREE="$(which FastTree)"
FT="$(echo ${FASTTREE} | rev | cut -c 9- | rev)fasttree"
ln -s "FASTTREE" "FT"

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac

mkdir buddysuite_env_build

if [ $machine == 'Mac' ]; then
    wget http://wasabiapp.org/download/pagan/pagan.osx.20150723.zip -O buddysuite_env_build/pagan.tgz
elif [ $machine == 'Linux' ]; then
    wget http://wasabiapp.org/download/pagan/pagan.linux64.20150723.tgz -O buddysuite_env_build/pagan.tgz
fi

tar -xzf buddysuite_env_build/pagan.tgz -C buddysuite_env_build
mv buddysuite_env_build/pagan /usr/local/
rm -r buddysuite_env_build
ln -s /usr/local/pagan/bin/pagan /usr/local/bin/