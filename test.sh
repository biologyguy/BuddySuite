#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
FAILURE=0

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests
echo "Pre-Tests"
ls -a -l
pwd

TEST_SCRIPTS='test_fixtures.py '
py.test ${TEST_SCRIPTS} --cache-clear --cov=conftest --cov=__init__ --cov-report= --cov-config ../.coveragerc --durations=10
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
ls -a -l
pwd
cat .coverage

mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddy_resources/

#### Buddy Resources
echo "************* Buddy Resources *************"
pwd
ls -la

cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddy_resources
TEST_SCRIPTS='test_buddy_resources.py '
args='${TEST_SCRIPTS} --cache-clear --cov=buddy_resources --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10'
echo 'py.test ${args}'
echo echo "*******************************************"

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddysuite/

#### BuddySuite
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddysuite
TEST_SCRIPTS='test_buddysuite.py '
py.test ${TEST_SCRIPTS} --cache-clear --cov=BuddySuite --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_alignbuddy/

#### AlignBuddy
TEST_SCRIPTS='test_alb_class_and_helpers.py '
TEST_SCRIPTS+='test_alb_api.py '
TEST_SCRIPTS+='test_alb_ui.py '
TEST_SCRIPTS+='test_alb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear --cov=AlignBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy/

#### DatabaseBuddy
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy
TEST_SCRIPTS='test_db_class_and_helpers.py '
TEST_SCRIPTS+='test_db_clients.py '
TEST_SCRIPTS+='test_db_ui.py '
py.test ${TEST_SCRIPTS} --cache-clear --cov=DatabaseBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_seqbuddy/

#### SeqBuddy
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_seqbuddy
TEST_SCRIPTS='test_sb_class_and_helpers.py '
TEST_SCRIPTS+='test_sb_api.py '
TEST_SCRIPTS+='test_sb_ui.py '
TEST_SCRIPTS+='test_sb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear --cov=SeqBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_phylobuddy/

#### PhyloBuddy
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_phylobuddy
TEST_SCRIPTS='test_pb_class_and_helpers.py '
TEST_SCRIPTS+='test_pb_api.py '
TEST_SCRIPTS+='test_pb_ui.py '
TEST_SCRIPTS+='test_pb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cache-clear --cov=PhyloBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/

#### Run Coveralls
cd /home/travis/build/biologyguy/BuddySuite
coveralls

exit ${FAILURE}
