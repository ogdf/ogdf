#!/bin/bash
#
# Runs all example binaries it can find in doc/examples
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

runExample () {
	echo "Running $1"
	cd `dirname $1` || exit 1
	./`basename $1` >/dev/null
	ret=$?
	if [ "$ret" -ne 0 ]
	then
		echo "$1 exited with code $ret"
		exit $ret
	fi
	exit 0
}

export -f runExample

echo "::group::($(date -Iseconds)) Running examples"
$OGDF_FIND doc/examples -name 'ex-*' |
 grep -v 'ex-multilevelmixer' |
 $OGDF_XARGS -I "{}" bash -c "runExample {}"
echo "::endgroup::"
