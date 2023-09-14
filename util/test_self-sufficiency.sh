#!/bin/bash
#
# Checks if all headers are self-contained and guarded.
# If this script returns an error at least one header could not be compiled.
# Make is used to perform parallel compile calls.
#
# You can also call this script with specific header files as arguments
# to test if they are self-contained and guarded.
#
# Author: Tilo Wiedera

. util/util-functions.sh || exit 123

if [ -z "$CXX" ]
then
	CXX=c++
fi

headers="$*"
if [ -z "$headers" ]
then
	headers=$(cd include ; find ogdf -name '*.h')
fi

make_tmpdir $0
(cd $tmp && cmake -LA -DCMAKE_BUILD_TYPE=Debug -DOGDF_WARNING_ERRORS=ON .. > cmakelog.txt) || exit 1

pattern="OGDF_EXTRA_CXX_FLAGS:STRING="
extra_flags="$(sed -ne "s/$pattern//p" $tmp/cmakelog.txt)" || exit 1
target_names=""
targets=""

# These warnings are not meaningful for functions/variables that are used only in source files
disableUnused=(function variable)
for f in "${disableUnused[@]}"
do
	extra_flags="$extra_flags -Wno-error=unused-$f"
done

echo "Found extra flags: $extra_flags"

for file in $headers
do
	file=$(echo "$file" | sed 's/^[\.\/]*include\///')
	if [ -f include/$file ]
	then
		if [ ${file:0:9} != "ogdf/lib/" ]
		then
			object="${file}.o"
			target_names="$target_names $object"
			call="\$(CALL) $extra_flags -DHEADER=\"<$file>\""
			targets="$targets\n$object:\n\t@echo \"$file\"\n\t$call\n"
		fi
	fi
done

echo "CALL = @$CXX -Werror -O0 -march=native -Iinclude -I../include -std=c++17 > /dev/null -c ../util/self-sufficiency.cpp" > $tmp/Makefile
echo -e "all:$target_names $targets" >> $tmp/Makefile
make -j "$cores" -C $tmp
exit $?
