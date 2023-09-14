#!/bin/sh
#
# Print all headers that are not compiled (not included
# by either OGDF, COIN, or test sources).
#
# This is very hacky and relies on specific CMake behavior.
#
# Author: Stephan Beyer

. util/util-functions.sh || exit 123

make_tmpdir $0
tmp="$(realpath "$tmp")"

if [ "$#" -ne 1 ]
then
	die "Usage: $0 <existing CMake build directory>|--new"
fi

dir="$1"
if [ "$dir" = "--new" ]
then
	(
		cd $tmp &&
		 cmake -DCMAKE_RANLIB=/bin/true \
		       -DCMAKE_LINKER=/bin/true \
		       -DCMAKE_AR=/bin/true .. &&
		 cmake -DCMAKE_CXX_FLAGS=-MM .
	) || die "Could not run CMake"
	echo Generating output
	make -C $tmp build-all || die "make build-all failed ... this hack only works with certain versions of CMake, I am sorry."
	dir="$tmp"
fi

cd "$dir" || die "Directory '$dir' not found"
test -d CMakeFiles || die "There are no CMakeFiles in '$dir'"
$OGDF_FIND "CMakeFiles/" -name depend.make |
	$OGDF_XARGS cat |
	sed -ne 's/^[^#].*: \(.*hp\?p\?\)$/\1/p' |
	sort -u |
	$OGDF_XARGS realpath --relative-to="$OLDPWD" -q >"$tmp/used"
cd "$OLDPWD"
git ls-files --full-name -- '*/ogdf/*.h' '*.hpp' |
	sort >"$tmp/exist"
diff "$tmp/used" "$tmp/exist" | sed -ne 's/^> //p'
