#!/bin/bash
#
# Output all broken macro definitions
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

CODEFILES="$(git grep -le '#\s*define\s' src include | check_filter)"
HEADERFILES="$(echo "$CODEFILES" | tr ' ' '\n' | grep -e '\.h$')"
ret=0

grep -ne '#\s*define\s' $HEADERFILES | grep -ve 'define\s*OGDF_\|define\s*forall_'
if [ "$?" -eq 0 ]
then
	cat<<EOF >&2

The header file macros defined above (maybe commented out) are not
prefixed by 'OGDF_'. Please fix this.

EOF
	ret=1
fi

grep -ne '#\s*define\s*[A-Z0-9_]*[a-z]' $HEADERFILES | grep -ve 'define\s*forall_'
if [ "$?" -eq 0 ]
then
	cat<<EOF >&2

The macros defined above (maybe commented out) are not all-caps.
Please turn all lower-case letters into upper-case letters.

EOF
	ret=1
fi

grep -ne '#\s*define\s*[^ \t]*__' $HEADERFILES
if [ "$?" -eq 0 ]
then
	cat<<EOF >&2

The macros defined above (maybe commented out) contain double underscores (__)
which are reserved for compiler's internal use. Please use at most one underscore.

EOF
	ret=1
fi

exit $ret
