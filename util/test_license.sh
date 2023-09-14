#!/bin/bash
#
# Output or fix all files with broken license headers.
#
# Note that files without license headers are ignored.
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

handle_fix_option \
  'lists files where license headers differ from template' \
  'fixes license headers in files where they differ from template' \
  "$@"

make_tmpdir $0
filterbegin=' \\par License:'
filterend='\*\/'
filter="/$filterbegin/,/$filterend/p"

licenseDiffers () {
	sed -ne "$filter" "$1" >$tmp/actual
	! cmp $tmp/expect $tmp/actual >/dev/null
}

fixLicense () {
	awk 'BEGIN { v = 1 } /'"$filterbegin"'/ { v = 0 } v==1 { print $0 }' "$1"
	cat $tmp/expect
	awk 'BEGIN { v = 0 } v==1 { print $0 } /'"$filterend"'/ { v = 1 }' "$1"
}

performTest () {
	if licenseDiffers "$1"
	then
		if [ "$fixit" = "yes" ]
		then
			fixLicense "$1" >$tmp/new
			mv $tmp/new "$1"
			return 0
		else
			echo "Style error: License text in heading of $1 differs from template!"
			return 124
		fi
	fi
}
export tmp filter filterbegin filterend fixit
export -f licenseDiffers fixLicense performTest

# first check that the templates have the same license header
sed -ne "$filter" doc/template.cpp >$tmp/expect
if licenseDiffers doc/template.h
then
	echo
	echo "Ouch! License texts in doc/template.cpp and doc/template.h differ!"
	exit 1
fi

git grep -l "$filterbegin" |
 check_filter |
 $OGDF_XARGS -I{} bash -c "performTest {}"
if [ "$?" -ne 0 ]
then
	echo
	echo "You can fix that automatically using $0 -f"
	exit 1
fi
