#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
FAILURE=0


# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests
TEST_SCRIPTS='test_fixtures.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-report= --cov-config ../.coveragerc

if [ $? -ne 0 ]
then
    FAILURE=1
fi

#cp /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/conftest.py ./
#cp /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/__init__.py ./
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddy_resources/

#### Buddy Resources
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddy_resources
TEST_SCRIPTS='test_buddy_resources.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../../.coveragerc

if [ $? -ne 0 ]
then
    FAILURE=1
fi

#cp /home/travis/build/biologyguy/BuddySuite/buddysuite/buddy_resources.py ./
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_alignbuddy/

#### AlignBuddy
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_alignbuddy
TEST_SCRIPTS='test_alb_class_and_helpers.py '
TEST_SCRIPTS+='test_alb_api.py '
TEST_SCRIPTS+='test_alb_ui.py '
TEST_SCRIPTS+='test_alb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../../.coveragerc
if [ $? -ne 0 ]
then
    FAILURE=1
fi
#cp /home/travis/build/biologyguy/BuddySuite/buddysuite/AlignBuddy.py ./
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy/

#### DatabaseBuddy
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy
TEST_SCRIPTS='test_db_class_and_helpers.py '
TEST_SCRIPTS+='test_db_clients.py '
TEST_SCRIPTS+='test_db_ui.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../../.coveragerc
if [ $? -ne 0 ]
then
    FAILURE=1
fi
#cp /home/travis/build/biologyguy/BuddySuite/buddysuite/DatabaseBuddy.py ./
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_seqbuddy/

#### SeqBuddy
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_seqbuddy
TEST_SCRIPTS='test_sb_class_and_helpers.py '
TEST_SCRIPTS+='test_sb_api.py '
TEST_SCRIPTS+='test_sb_ui.py '
#TEST_SCRIPTS+='test_sb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../../.coveragerc
if [ $? -ne 0 ]
then
    FAILURE=1
fi
#cp /home/travis/build/biologyguy/BuddySuite/buddysuite/SeqBuddy.py ./
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_phylobuddy/

#### PhyloBuddy
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_phylobuddy
TEST_SCRIPTS='test_pb_class_and_helpers.py '
TEST_SCRIPTS+='test_pb_api.py '
TEST_SCRIPTS+='test_pb_ui.py '
TEST_SCRIPTS+='test_pb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../../.coveragerc
if [ $? -ne 0 ]
then
    FAILURE=1
fi
#cp /home/travis/build/biologyguy/BuddySuite/buddysuite/PhyloBuddy.py ./
mv .coverage /home/travis/build/biologyguy/BuddySuite/
cd /home/travis/build/biologyguy/BuddySuite
coveralls

exit ${FAILURE}
