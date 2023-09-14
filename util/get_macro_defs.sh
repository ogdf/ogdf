#!/bin/bash
#
# Outputs a -Dmacro list of all uncommon and undefined OGDF_-macros
# that are "used" via an "#ifdef"
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

make_tmpdir $0

CODEFILES="$(git grep -le '#\s*ifdef\s' src include | check_filter)"
grep -e '#\s*ifdef\s' $CODEFILES |
	sed -ne 's/^.*ifdef *\([A-Za-z0-9_]*\).*$/\1/p' |
	grep -e '^OGDF_' |
	sort -u >$tmp/macro-checks

# filter out macros mentioned in doc/defines.md
sed -ne 's/^| `\([A-Za-z0-9_]*\)`.*$/\1/p' doc/defines.md >$tmp/macro-defines

# filter out macros that are actually defined
git grep '^\s*#\s*define\s' | sed -ne 's/^.*:[[:space:]]*#[[:space:]]*define[[:space:]]*\([A-Za-z0-9_]*\).*$/\1/p' >>$tmp/macro-defines

for macro in `cat $tmp/macro-checks`
do
	grep -q "^$macro\$" $tmp/macro-defines || echo -n "-D$macro "
done
echo
