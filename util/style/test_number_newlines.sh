#!/bin/bash
#
# Output or fix all files with
#  - more than MaxEmptyLines (see below) consecutive empty lines, and
#  - less or more than one newline at EOF.
#
# Author: Stephan Beyer

MaxEmptyLines=4

. util/util-functions.sh || exit 123

handle_fix_option 'lists files without newlines at EOF' 'fixes files without newlines at EOF' "$@"

make_tmpdir $0

fail_with_info () {
	echo
	echo "The files mentioned above do not end with exactly one newline,"
	echo "or there are more than $MaxEmptyLines consecutive empty lines."
	echo
	echo "You can fix that automatically using $0 -f"
	exit 1
}

performTest () {
	found=0
	awk '/^$/ { if (length(empties) < '$MaxEmptyLines') { empties = empties "\n" }; next } { printf "%s", empties; print; empties = "" }' "$1" >$tmp/new
	cmp "$1" $tmp/new >/dev/null 2>&1 || found=1
	if [ "$found" = "1" ]
	then
		echo "  $1"
	fi
	if [ "$found,$fixit" = "1,yes" ]
	then
		mv $tmp/new "$1"
		found=0
	fi
	exit $found
}

export tmp fixit MaxEmptyLines
export -f performTest
git ls-files |
 check_filter |
 grep -ve '\.\(cdr\|png\|svg\)$\|^test\/resources' |
 $OGDF_XARGS -I{} bash -c "performTest {}" || fail_with_info
