#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

echo "test_fixtures.py"
py.test tests/test_fixtures.py  --cov tests/conftest.py --cov-report html -p no:cacheprovider

echo "buddy_resources.py"
py.test tests/test_buddy_resources/ --cov buddy_resources.py --cov-report html -p no:cacheprovider

echo "AlignBuddy.py"
py.test tests/test_align_buddy/ --cov AlignBuddy.py --cov-report html -n 8 -p no:cacheprovider

echo "SeqBuddy.py"
py.test tests/test_seqbuddy/ --cov SeqBuddy.py --cov-report html -n 8 -p no:cacheprovider
