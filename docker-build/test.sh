#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
source activate test-environment

cd /home/docker/BuddySuite/buddysuite/tests
find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

#sudo git pull

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
TEST_SCRIPTS='test_fixtures.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider

cd /home/docker/BuddySuite/buddysuite/tests/test_buddy_resources
#### Buddy Resources
TEST_SCRIPTS='test_buddy_resources.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider

cd /home/docker/BuddySuite/buddysuite/tests/test_align_buddy
#### AlignBuddy
TEST_SCRIPTS='test_alb_class_and_helpers.py '
TEST_SCRIPTS+='test_alb_api.py '
TEST_SCRIPTS+='test_alb_ui.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider
