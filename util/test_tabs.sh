#!/bin/sh
#
# Check that tab-indentation is not mixed with spaces
# (tabs and spaces afterwards is allowed)
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

regex='^\t* +\t'
files="$(git ls-files README.md CMakeLists.txt cmake doc util | $OGDF_XARGS grep -l -P "$regex")"

if [ -n "$files" ]
then
	echo
	echo "The following lines contain tab-indentation mixed with spaces:"
	echo
	grep -H -n -P "$regex" $files
	echo
	echo "Please fix that!"
	exit 1
fi
exit 0
