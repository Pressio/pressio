#!/bin/bash

set -e

source bash_colors.sh

# check cmd line args
if [ "$#" -ne 1 ]; then
  echo "usage: $0 <target-dir> [extension]"
  exit 1
fi

# store
targetdir=$1
licensefile=$2

if [ "$#" -gt 1 ];
then
  pattern=$2
else
  pattern='*.hpp'
fi

echo ""
echo "${fgyellow}Running on targetdir=$targetdir, extension=$pattern ${fgrst}"

allfiles=$(find $targetdir -name $pattern -o -name $pattern)
#echo ${allfiles[@]}
for file in $allfiles
do
  filename=$(basename $file)
  filedir=$(dirname $file)

  if ! grep -q 'NTESS' "$file"; then
      if ! grep -q 'Ennio' "$file"; then
	 echo "found file $file without a Sandia or other license. Terminating."
	 exit 11
      fi
  else
      filetmp=$(mktemp)
      tail -n +48 $file > $filetmp
      mv $filetmp $file
  fi
done
