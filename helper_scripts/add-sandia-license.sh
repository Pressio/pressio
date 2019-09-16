#!/bin/bash

set -e

source bash_colors.sh

# check cmd line args
if [ "$#" -lt 2 ]; then
  echo "usage: $0 <target-dir> <license-template-file> [extension]"
  exit 1
fi

# store
targetdir=$1
licensefile=$2

if [ "$#" -gt 2 ];
then
  ext=$3
  pattern="*.${ext}"
else
  pattern='*.hpp'
fi

echo ""
echo "${fgyellow}Running on targetdir=$targetdir, license=$licensefile, extension=$pattern ${fgrst}"

allfiles=$(find $targetdir -name $pattern -o -name $pattern)
#echo ${allfiles[@]}
for file in $allfiles
do
  filename=$(basename $file)
  filedir=$(dirname $file)

  if ! grep -q 'Ennio' "$file";
  then
      echo $filename
      echo $filedir
      licensetmp=$(mktemp)
      filetmp=$(mktemp)
      echo "${fggreen}Running on file=$file, targetdir=$filedir, license=$licensefile, tmp=${filetmp} ${fgrst}"
      sed 's/\<file-name.hpp\>/'$filename'/g' $licensefile > $licensetmp
      cat $licensetmp $file > $filetmp
      mv $filetmp $file
  fi
done
