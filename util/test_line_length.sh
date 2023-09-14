#!/bin/sh
#
# Check that a specific line length is not exceeded
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

TABSIZE=8
LENGTH=200
TABS="$((LENGTH/TABSIZE))"
CHARS="$((LENGTH-TABS*TABSIZE))"
regex="(\\t|.{$TABSIZE}){$TABS}.{$CHARS}"

files="$(git ls-files CMakeLists.txt cmake README.md 'doc/*.md' util | $OGDF_XARGS grep -l -P "$regex")"

if [ -n "$files" ]
then
	echo
	echo "The following lines are too long:"
	echo
	grep -H -n -P "$regex" $files
	echo
	echo "Length limit is set to $LENGTH where each tab counts for $TABSIZE characters."
	exit 1
fi
exit 0
