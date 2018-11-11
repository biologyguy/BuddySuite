#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
FAILURE=0

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests
printf "
************************** Pre-Tests **************************
"
pwd
ls -la

TEST_SCRIPTS='test_fixtures.py '
args="$TEST_SCRIPTS --cache-clear --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

#py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
ls -a -l
pwd

#### Buddy Resources
printf "
************************** Buddy Resources **************************
"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddy_resources
pwd
ls -la
TEST_SCRIPTS='test_buddy_resources.py '
args="$TEST_SCRIPTS --cache-clear --cov=buddy_resources  --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddysuite/

#### BuddySuite
printf "
************************** BuddySuite **************************
"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddysuite
pwd
ls -la
TEST_SCRIPTS='test_buddysuite.py '
args="$TEST_SCRIPTS --cache-clear --cov=BuddySuite --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

#py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_alignbuddy/

#### AlignBuddy
printf "
************************** AlignBuddy **************************
"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_alignbuddy
pwd
ls -la
TEST_SCRIPTS='test_alb_class_and_helpers.py '
TEST_SCRIPTS+='test_alb_api.py '
TEST_SCRIPTS+='test_alb_ui.py '
TEST_SCRIPTS+='test_alb_3rd_party.py '
args="$TEST_SCRIPTS --cache-clear --cov=AlignBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

#py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy/

#### DatabaseBuddy
printf "
************************** DatabaseBuddy **************************
"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy
pwd
ls -la
#TEST_SCRIPTS='test_db_class_and_helpers.py '
#TEST_SCRIPTS+='test_db_clients.py '
TEST_SCRIPTS='test_db_ui.py '
args="$TEST_SCRIPTS --cache-clear --cov=DatabaseBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_seqbuddy/

#### SeqBuddy
printf "
************************** SeqBuddy **************************
"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_seqbuddy
pwd
ls -la
TEST_SCRIPTS='test_sb_class_and_helpers.py '
TEST_SCRIPTS+='test_sb_api.py '
TEST_SCRIPTS+='test_sb_ui.py '
TEST_SCRIPTS+='test_sb_3rd_party.py '
args="$TEST_SCRIPTS --cache-clear --cov=SeqBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

#py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_phylobuddy/

#### PhyloBuddy
printf "
************************** PhyloBuddy **************************
"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_phylobuddy
pwd
ls -la
TEST_SCRIPTS='test_pb_class_and_helpers.py '
TEST_SCRIPTS+='test_pb_api.py '
TEST_SCRIPTS+='test_pb_ui.py '
TEST_SCRIPTS+='test_pb_3rd_party.py '
args="$TEST_SCRIPTS --cache-clear --cov=PhyloBuddy --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

#py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/

#### Run Coveralls
cd /home/travis/build/biologyguy/BuddySuite
coveralls

exit ${FAILURE}
