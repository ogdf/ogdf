#!/bin/sh
#
# Output all files containing unscoped enumerations.
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

FILTER='^\s*enum\>'
FILES="$(git grep -l "$FILTER" -- src/ogdf/ include/ogdf/ | check_filter)"
git grep -n "$FILTER" -- $FILES | grep -v -e 'enum\s*class'
if [ $? -eq 0 ]
then
	echo
	echo "The code positions listed above contain unscoped enumerations."
	echo "Please replace them by C++11 enum classes or static const integers."
	exit 1
fi

exit 0
