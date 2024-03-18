#!/bin/bash
#
# Compiles and tests the OGDF (using an out-of-source build).
#
# Author: Tilo Wiedera, Ivo Hedtke, Stephan Beyer

. util/util-functions.sh || exit 123

if [ "$#" -lt 3 ]
then
	echo "Usage: $0 <1st> <2nd> <3rd> <4th> [cmake arguments]"
	echo
	echo "The first four arguments are passed to util/compile_for_tests.sh"
	echo "Any further arguments get passed to the CMake invocation."
	exit 2
fi

make_tmpdir $0

util/compile_for_tests.sh $1 $2 $3 $4 . $tmp "${@:5}" || exit
util/run_examples.sh || exit
util/perform_separate_tests.sh $tmp || exit
echo "($(date -Iseconds)) Success"
