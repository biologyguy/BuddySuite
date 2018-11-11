#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
FAILURE=0

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests
printf "
\e[96m************************** Pre-Tests **************************\e[39m"
pwd
ls -la

TEST_SCRIPTS='test_fixtures.py '
args="$TEST_SCRIPTS --cache-clear --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
ls -a -l
pwd

#### Buddy Resources
printf "
\e[96m************************** Buddy Resources **************************\e[39m"
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
\e[96m************************** BuddySuite **************************\e[39m"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_buddysuite
pwd
ls -la
TEST_SCRIPTS='test_buddysuite.py '
args="$TEST_SCRIPTS --cache-clear --cov=BuddySuite --cov-append --cov-report= --cov-config ../../.coveragerc --durations=10"
echo "py.test $args"
echo "**************************************************************************************
"

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_alignbuddy/

#### AlignBuddy
printf "
\e[96m************************** AlignBuddy **************************\e[39m"
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

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy/

#### DatabaseBuddy
printf "
\e[96m************************** DatabaseBuddy **************************\e[39m"
cd /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_databasebuddy
pwd
ls -la
TEST_SCRIPTS='test_db_class_and_helpers.py '
TEST_SCRIPTS+='test_db_clients.py '
TEST_SCRIPTS+='test_db_ui.py '
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
\e[96m************************** SeqBuddy **************************\e[39m"
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

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/buddysuite/tests/test_phylobuddy/

#### PhyloBuddy
printf "
\e[96m************************** PhyloBuddy **************************\e[39m"
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

py.test ${args}
if [[ $? -ne 0 ]]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/BuddySuite/

#### Run Coveralls
cd /home/travis/build/biologyguy/BuddySuite
coveralls

exit ${FAILURE}
