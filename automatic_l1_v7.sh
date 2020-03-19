#! /bin/sh

WORKINGDIR=$1
FILETLE=$2
TIMEOFFSET=$3


if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  echo "Usage: two parameters are mandatory:" 
  echo "     1) working directory with level0 files"
  echo "     2) TLE file (with complete path) to be used"
  echo "Optionally you can set the CPU timeoffset passing a third parameter."
  exit 0
fi


for file in $WORKINGDIR/*self.root
do
  ls "$file"
  ./euso_l1_tle_v7 "$file"  $FILETLE  $TIMEOFFSET
  ls "$file"
  sleep 5s 
done
