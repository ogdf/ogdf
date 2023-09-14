#!/bin/bash
#
# Executes separate test executables in the order specified below.
# Terminates upon the first failed tests or if all tests have passed.
# Expects the path of the executables as first argument.
# Must be called from the root of the OGDF source directory.
#
# Note that this script requires OGDF_SEPARATE_TEST to be enabled in
# your cmake configuration. The sole purpose of separated tests is
# to accelerate testing during continuous integration.
#
# Author: Tilo Wiedera, Ivo Hedtke (Mac)

. util/util-functions.sh || exit 123

unamestr=`uname`

export format="%Us, %MKb"
gtime -f "$format" true &>/dev/null
if [ $? -eq 0 ]
then
	export gtimeCmd="gtime"
else
	/usr/bin/time -f "$format" true &>/dev/null
	if [ $? -eq 0 ]
	then
		export gtimeCmd="/usr/bin/time"
	fi
fi

performTest() {
	echo "starting: $1"

	$2/$1 --reporter=crash --only="7469310868675805"
	result=$?

	if [ $result != 0 ] && [ $result != 1 ]
	then
		echo " FAILED SKIP-RUN: $1 (\`$2/$1\` returned $result)"
		return 255
	fi

	logfile="$3/$$"

	if [ -n "$gtimeCmd" ]
	then
		$gtimeCmd -f "$format" $2/$1 --break-on-failure --reporter=crash >"$logfile" 2>&1
	else
		$2/$1 --break-on-failure --reporter=crash >"$logfile" 2>&1
	fi

	result=$?

	if [ $result -eq 0 ]
	then
		if [ -n "$gtimeCmd" ]
		then
			measured="($(tail -n 1 "$logfile"))"
		fi
		echo "  succeeded: $1 $measured"
		return 0
	else
		echo "  FAILED: $1"
		# check for address sanitizer output
		SANITIZER_OUTPUT="$(cat "$logfile" | sed -ne '/^==[0-9]\+==ERROR/,$p')"
		if [ -n "$SANITIZER_OUTPUT" ]
		then
			printf '\nThere are memory issues!\n\n%s\n' "$SANITIZER_OUTPUT"
		else
			# no sanitizer output? Output errors or last 10 lines
			echo "From the log:"
			grep -B 1 -A 5 "^#" "$logfile" || tail -n 11 "$logfile"
		fi
		return 255
	fi
}

export -f performTest
dir=$1/test/bin
make_tmpdir $0 $1

tests=$($OGDF_FIND $dir -maxdepth 1 -type f -executable -printf '%f\n')
count=$(echo "$tests" | wc -l)

if [ $count == 0 ]
then
	echo "Could not find any tests."
	exit 1
fi

echo "::group::($(date -Iseconds)) Found $count tests. Starting execution..."
echo "$tests" | $OGDF_XARGS -P "$cores" -I "{}" bash -c "performTest {} $dir $tmp"
ret=$?
echo "::endgroup::"
exit $ret
