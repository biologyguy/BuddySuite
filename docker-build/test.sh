#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"

cd /home/docker/BuddySuite/buddysuite/tests

# py.test caches things, which may be incompatible with the Docker environment. Delete any caches before proceeding
find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
TEST_SCRIPTS='test_fixtures.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider

#### Buddy Resources
cd /home/docker/BuddySuite/buddysuite/tests/test_buddy_resources
TEST_SCRIPTS='test_buddy_resources.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider

#### AlignBuddy
cd /home/docker/BuddySuite/buddysuite/tests/test_alignbuddy
TEST_SCRIPTS='test_alb_class_and_helpers.py '
TEST_SCRIPTS+='test_alb_api.py '
TEST_SCRIPTS+='test_alb_ui.py '
TEST_SCRIPTS+='test_alb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider

#### DatabaseBuddy
cd /home/docker/BuddySuite/buddysuite/tests/test_databasebuddy
TEST_SCRIPTS='test_db_class_and_helpers.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider

#### SeqBuddy
cd /home/docker/BuddySuite/buddysuite/tests/test_seqbuddy
TEST_SCRIPTS='test_sb_class_and_helpers.py '
TEST_SCRIPTS+='test_sb_api.py '
TEST_SCRIPTS+='test_sb_ui.py '
TEST_SCRIPTS+='test_sb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider

#### PhyloBuddy
cd /home/docker/BuddySuite/buddysuite/tests/test_phylobuddy
TEST_SCRIPTS='test_pb_class_and_helpers.py '
TEST_SCRIPTS+='test_pb_api.py '
TEST_SCRIPTS+='test_pb_ui.py '
TEST_SCRIPTS+='test_pb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider
