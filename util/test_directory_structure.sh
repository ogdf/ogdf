#!/bin/bash
#
# Checks if there are any source files in the include directory.
# Checks if there are any include files in the source directory.
# Prints a warning for each source file that is missing a header file,
# with some allowed exceptions.
#
# Author: Tilo Wiedera, Stephan Beyer

. util/util-functions.sh || exit 123

printFiles () {
	echo ": ${#@}"
	for file in $@; do
		result=1
		echo $file
	done
}

result=0

printf "sources in include/ogdf"
printFiles  $($OGDF_FIND include/ogdf -name '*.cpp')

printf "\nheaders in src/ogdf"
printFiles $($OGDF_FIND src/ogdf -name '*.h')

echo ""

foundSolitary=0
for file in $(find src/ogdf -name '*.cpp' -print | grep -ve src/ogdf/lib -e src/ogdf/external); do
	length=${#file}-13
	name=${file:9:$length}.h

	if [ ! -f "include/ogdf/$name" ] ; then
		if [ ! -f "include/ogdf/${name%_*}.h" ] ; then
			echo "WARNING: Encountered solitary source file: $file"
			foundSolitary=1
		else
			echo "NOTE: Encountered allowed (suffixed) solitary source file: $file"
		fi
	fi
done

if [ $foundSolitary != 0 ]
then
	echo ""
	echo "Solitary source files are missing a matching header file."
	echo "Please make sure that these file are required as they conflict with our naming convention."
	result=1
fi


exit $result
