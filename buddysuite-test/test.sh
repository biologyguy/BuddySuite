#! /bin/bash

export PATH="$HOME/miniconda/bin:$PATH"
source activate test-environment

ls /home/docker

TEST_SCRIPTS=()
#TEST_SCRIPTS+=('alignbuddy_test.py')
#TEST_SCRIPTS+=('phylobuddy_test.py')
#TEST_SCRIPTS+=('sequddy_test.py')
TEST_SCRIPTS+=('test_fixtures.py')
TEST_SCRIPTS+=('test_buddy_resources/test_buddy_resources.py')
TEST_SCRIPTS+=('test_align_buddy/test_alb_class_and_helpers.py')
TEST_SCRIPTS+=('test_align_buddy/test_alb_api.py')
TEST_SCRIPTS+=('test_align_buddy/test_alb_ui.py')

cd ~/BuddySuite/buddysuite/tests/
# disable cacheprovider since it requires r/w access to the test directory
py.test -p no:cacheprovider -m "not internet and not display" "${@}" $TEST_SCRIPTS
