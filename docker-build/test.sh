#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
source activate test-environment

cd /buddysuite/tests

#sudo git pull

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
TEST_SCRIPTS='test_fixtures.py '
py.test ${TEST_SCRIPTS} -p no:cacheprovider

#### Buddy Resources
TEST_SCRIPTS='test_buddy_resources/test_buddy_resources.py '
py.test ${TEST_SCRIPTS} -p no:cacheprovider

#### AlignBuddy
TEST_SCRIPTS='test_align_buddy/test_alb_class_and_helpers.py '
TEST_SCRIPTS+='test_align_buddy/test_alb_api.py '
TEST_SCRIPTS+='test_align_buddy/test_alb_ui.py '
py.test ${TEST_SCRIPTS} -p no:cacheprovider
