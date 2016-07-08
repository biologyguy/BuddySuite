#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "test_fixtures.py"
py.test ${DIR}/tests/test_fixtures.py  --cov ${DIR}/tests/conftest.py --cov-report html -p no:cacheprovider

echo "buddy_resources.py"
py.test ${DIR}/tests/test_buddy_resources/ --cov ${DIR}/buddy_resources.py --cov-report html -p no:cacheprovider

echo "AlignBuddy.py"
py.test ${DIR}/tests/test_alignbuddy/ --cov ${DIR}/AlignBuddy.py --cov-report html -n 8 -p no:cacheprovider

echo "SeqBuddy.py"
py.test ${DIR}/tests/test_seqbuddy/ --cov ${DIR}/SeqBuddy.py --cov-report html -n 8 -p no:cacheprovider

echo "PhyloBuddy.py"
py.test ${DIR}/tests/test_phylobuddy/ --cov ${DIR}/PhyloBuddy.py --cov-report html -n 8 -p no:cacheprovider
