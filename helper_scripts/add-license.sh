#!/bin/bash

set -e 

# check cmd line args
if [ test $# -lt 2 ]; then
  echo "usage: $0 <target-dir> <license-template-file> [extension]"
  exit 1
fi

# store
targetdir=$1
licensefile=$2

if [ test $# -gt 2 ];then
  pattern=$3
else
  pattern="*.h"
fi

echo "Running on targetdir=$targetdir, license=$licensefile, extension=$extension"

allfiles=$(find $dir -iname "$extension")
for file in $allfiles
do
  filename=$(basename $file)
  filedir=$(dirname $file)
  licensetmp=$(mktemp)
  filetmp=$(mktemp)
  echo "Running on file=$file, targetdir=$filedir, license=$licensefile, tmp=$filetmp"
  sed 's/\<file-name.hpp\>/'$filename'/g' $licensefile > $licensetmp
  cat $licensetmp $file > $filetmp
  mv $filetmp $file
done
