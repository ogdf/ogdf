#!/bin/bash
#
# Output or fix all files where doxygen comments are not indented properly
#
# Author: Simon D. Fink

. util/util-functions.sh || exit 123

handle_fix_option 'lists files with doxygen comments with improper indentation' 'fixes indentation of doxygen comments' "$@"

fail_with_info () {
	echo
	echo "The files mentioned above have doxygen comments that are not indented properly."
	echo
	echo "You can fix that automatically using $0 -f"
	exit 1
}

args=""
if [ "$fixit" = "yes" ]; then
  args="-fi"
fi

git ls-files src include test/src test/include |
  check_filter |
  grep -e '\.\(c\|h\)\(pp\)\?$' |
  $OGDF_XARGS -P "$cores" -I{} python3 util/style/indent_comments.py $args "{}" || fail_with_info
