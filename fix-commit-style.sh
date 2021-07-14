#!/usr/bin/bash

#
# Run this script to fix styling in your last commit (changes only selected lines)
#
# Note: need clang-format-diff.py in path (e.g. /usr/share/clang/clang-format-10)

git diff -U0 --no-color HEAD^ | clang-format-diff.py -i -p1