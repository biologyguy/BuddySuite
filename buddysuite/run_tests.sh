#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "test_fixtures.py"
TEST_SCRIPTS=${DIR}'/tests/test_fixtures.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/tests/conftest.py --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "buddy_resources.py"
TEST_SCRIPTS=${DIR}'/tests/test_buddy_resources/ '
py.test ${TEST_SCRIPTS} --cov ${DIR}/buddy_resources.py --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "AlignBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_alignbuddy/test_alb_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_alignbuddy/test_alb_api.py '
TEST_SCRIPTS+=${DIR}'/tests/test_alignbuddy/test_alb_ui.py '
#TEST_SCRIPTS+=${DIR}'/tests/test_alignbuddy/test_alb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/AlignBuddy.py --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "DatabaseBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_databasebuddy/test_db_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_databasebuddy/test_db_clients.py '
TEST_SCRIPTS+=${DIR}'/tests/test_databasebuddy/test_db_ui.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/DatabaseBuddy.py --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "PhyloBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_phylobuddy/test_pb_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_phylobuddy/test_pb_api.py '
TEST_SCRIPTS+=${DIR}'/tests/test_phylobuddy/test_pb_ui.py '
#TEST_SCRIPTS+=${DIR}'/tests/test_phylobuddy/test_pb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/PhyloBuddy.py --cov-report html -n 8 -p no:cacheprovider --durations=10 $@

echo "SeqBuddy.py"
TEST_SCRIPTS=${DIR}'/tests/test_seqbuddy/test_sb_class_and_helpers.py '
TEST_SCRIPTS+=${DIR}'/tests/test_seqbuddy/test_sb_api.py '
TEST_SCRIPTS+=${DIR}'/tests/test_seqbuddy/test_sb_ui.py '
#TEST_SCRIPTS+=${DIR}'/tests/test_seqbuddy/test_sb_3rd_party.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/SeqBuddy.py --cov-report html -n 8 -p no:cacheprovider --durations=10 $@
