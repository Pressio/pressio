#!/usr/bin/bash

#
# Run this script to fix styling with clang-format in all source files
#

find . -iname *.h -o -iname *.c -o -iname *.cpp -o -iname *.hpp \
    | xargs clang-format -style=file -i -fallback-style=none