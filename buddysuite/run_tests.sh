#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "test_fixtures.py"
TEST_SCRIPTS=${DIR}'/tests/test_fixtures.py '
#py.test ${TEST_SCRIPTS} --cov conftest --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "buddy_resources.py"
TEST_SCRIPTS=${DIR}'/tests/test_buddy_resources/ '
py.test ${TEST_SCRIPTS} --cov buddy_resources --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "BuddySuite.py"
TEST_SCRIPTS=${DIR}'/tests/test_buddysuite/ '
py.test ${TEST_SCRIPTS} --cov BuddySuite --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "AlignBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_alignbuddy/test_alb_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_alignbuddy/test_alb_api.py '
TEST_SCRIPTS+=${DIR}'/tests/test_alignbuddy/test_alb_ui.py '
TEST_SCRIPTS+=${DIR}'/tests/test_alignbuddy/test_alb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cov AlignBuddy --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "DatabaseBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_databasebuddy/test_db_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_databasebuddy/test_db_clients.py '
TEST_SCRIPTS+=${DIR}'/tests/test_databasebuddy/test_db_ui.py '
py.test ${TEST_SCRIPTS} --cov DatabaseBuddy --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "PhyloBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_phylobuddy/test_pb_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_phylobuddy/test_pb_api.py '
TEST_SCRIPTS+=${DIR}'/tests/test_phylobuddy/test_pb_ui.py '
TEST_SCRIPTS+=${DIR}'/tests/test_phylobuddy/test_pb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cov PhyloBuddy --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "SeqBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_seqbuddy/test_sb_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_seqbuddy/test_sb_api.py '
TEST_SCRIPTS+=${DIR}'/tests/test_seqbuddy/test_sb_ui.py '
TEST_SCRIPTS+=${DIR}'/tests/test_seqbuddy/test_sb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cov SeqBuddy --cov-report html -n 8 -p no:cacheprovider --durations=10 $@
