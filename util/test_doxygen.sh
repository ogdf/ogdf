#!/bin/bash
#
# Checks if doxygen generation raises any warnings or errors.
#
# Author: Tilo Wiedera, Ivo Hedtke, Stephan Beyer

. util/util-functions.sh || exit 123

make_tmpdir $0
file=$tmp/docs.txt
source_dir=`pwd`

(cd $tmp && cmake -DOGDF_FULL_DOC=ON $source_dir)
make -s -C $tmp doc 2>&1 | tee $file

echo
echo "### SUMMARY OF DOXYGEN ERRORS ###"
echo

grep -v -e '\.md:[0-9]*: warning' -e 'target doc' -e 'not generated, too many nodes' $file
if [ $? -eq 0 ]
then
	exit 1
fi
