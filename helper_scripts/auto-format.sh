#!/usr/bin/env bash

#
# Run this script to fix styling with clang-format in all source files.
#
# This script must be run where it is, no args needed.
#

find ../include -iname *.h -o -iname *.c -o -iname *.cpp -o -iname *.hpp \
  | xargs clang-format -style=file -i -fallback-style=none