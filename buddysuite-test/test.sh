#! /bin/bash

export PATH="$HOME/miniconda/bin:$PATH"
source activate test-environment

TEST_SCRIPTS=()
#TEST_SCRIPTS+=('alignbuddy_test.py')
TEST_SCRIPTS+=('phylobuddy_test.py')
#TEST_SCRIPTS+=('sequddy_test.py')

cd ~/buddysuite/workshop
# disable cacheprovider since it requires r/w access to the test directory
py.test -p no:cacheprovider -m "not internet" "${@}" $TEST_SCRIPTS
