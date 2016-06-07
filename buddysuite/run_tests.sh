#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

py.test tests/test_fixtures.py  --cov tests/conftest.py --cov-report html -p no:cacheprovider
py.test tests/test_buddy_resources/ --cov buddy_resources.py --cov-report html -p no:cacheprovider
py.test tests/test_align_buddy/ --cov AlignBuddy.py --cov-report html -n 8 -p no:cacheprovider
