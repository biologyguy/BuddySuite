#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

echo "test_fixtures.py"
cd tests
py.test test_fixtures.py  --cov conftest.py --cov-report html -p no:cacheprovider

echo "buddy_resources.py"
py.test test_buddy_resources/ --cov ../buddy_resources.py --cov-report html -p no:cacheprovider

echo "AlignBuddy.py"
py.test test_alignbuddy/ --cov ../AlignBuddy.py --cov-report html -n 8 -p no:cacheprovider

echo "SeqBuddy.py"
py.test test_seqbuddy/ --cov ../SeqBuddy.py --cov-report html -n 8 -p no:cacheprovider

echo "PhyloBuddy.py"
py.test test_phylobuddy/ --cov ../PhyloBuddy.py --cov-report html -n 8 -p no:cacheprovider
