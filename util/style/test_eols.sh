#!/bin/bash
#
# Output files with broken EOLs or remove broken EOLs
# (like DOS format or trailing whitespace)
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

handle_fix_option 'output all files with broken EOLs' 'remove dangling whitespace and DOS EOL format' "$@"

regex='[ \t\r]$'
files="$(git ls-files | grep -ve '\.bat$\|^_\|^test/resources/' | $OGDF_XARGS grep -l -I -P "$regex")"

fix () {
	for file in $files
	do
		echo Processing $file
		sed -i 's/[ \t\r]\+$//' $file
	done
}

output () {
	cat<<EOF

The following files contain lines with whitespace at the end of
the line or are in DOS format (\\r character at the end of the line).

$(grep -H -n -o -P "$regex" $files)

Please remove these trailing whitespace, eg, by using

  $0 -f

EOF
	exit 1
}

if [ -n "$files" ]
then
	if [ "$fixit" = "yes" ]
	then
		fix
	else
		output
	fi
fi

exit 0
