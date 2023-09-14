#!/bin/bash
#
# Executes the clang static analyzer and prints all encountered problems.
#
# Author: Tilo Wiedera

if [ "$1" != "-f" ]
then
        echo "Running this script will modify your build configuration."
        echo "Run with -f to execute."
        exit 2
fi

mkdir -p tmp

(cd ./tmp && scan-build cmake -DCMAKE_BUILD_TYPE=Debug -DOGDF_WARNING_ERRORS=ON ..)
file=./tmp/analysis-$$.txt
scan-build make -C tmp build-all 2> >(tee $file >&2)

if [ -s $file ]
then
  exit 1
else
  exit 0
fi
