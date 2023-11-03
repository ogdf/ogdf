#!/bin/sh
#
# Output all files containing typedefs.
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

git grep -l '^\s*typedef\>' -- src/ogdf/ include/ogdf/ | check_filter
if [ $? -eq 0 ]
then
	echo
	echo "The source files listed above contain typedefs."
	echo "Please replace them by C++11 'using'."
	exit 1
fi

exit 0
