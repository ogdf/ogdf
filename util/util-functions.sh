# Common functions for shell scripts in util/
#
# Author: Stephan Beyer

# get number of cores
cores="$(getconf _NPROCESSORS_ONLN)" || cores=1

fixit=no
tmp=''

# Sets fixit=yes if -f option was used.
# If other options are given, die with usage information.
#
# arg 1: usage description for invocations without -f
# arg 2: usage description for invocations with -f
# remaining args = script args "$@"
handle_fix_option () {
	usage_without="$1"
	usage_with="$2"
	shift 2

	if [ -n "$1" ]
	then
		if [ "$1" = "-f" ]
		then
			fixit=yes
		else
			cat <<EOF
Usage:

  $0
     $usage_without

  $0 -f
     $usage_with
EOF
			exit 255
		fi
	fi
}

# Sets tmp variable to a newly created temporary directory
# that is removed at exit
#
# arg 1: an identifier for the temporary directory [recommended: the script name ($0)]
# arg 2 (optional): a subdirectory it should be created in [default: current directory]
make_tmpdir () {
	tmpprefix="$(pwd)/"
	test -n "$2" && tmpprefix="$2/"
	tmp="${tmpprefix}tmpdir-`basename $1 | sed -e 's/\..*//'`-$$"
	while [ -e $tmp ] # should not happen, but just in case
	do
		tmp=${tmp}_
	done

	trap "rm -rf $tmp" EXIT
	mkdir -p $tmp || exit 255
}

# Exits with an error message
die () {
	echo "Fatal error in $0: $*" >&2
	exit 42
}

# Filters files from input: coin, lib, _legacy, bandit, tinydir.h, .clang-tidy
#
# Typical use case: git ls-files src include test/src | check_filter
check_filter () {
	grep -v -e '/ogdf/lib/\|/coin/\|_legacy/\|/bandit/\|tinydir\.h\|\.clang-tidy'
}

# If $1 == xargs, set the variable OGDF_XARGS to the GNU xargs command (i.e.,
# xargs or gxargs). If $1 == find, do the analogous for OGDF_FIND.
use_gnu_tool () {
	var="$(echo "OGDF_$1" | tr a-z A-Z)"
	if "$1" --version 1> /dev/null 2>&1; then
		eval "export $var=$1"
	elif "g$1" --version 1> /dev/null 2>&1; then
		eval "export $var=g$1"
	else
		die "GNU $1 not found!"
	fi
}

# we sometimes need realpath, xargs, and find on Mac where the BSD version
# differs from the GNU one (you might have to install findutils for it to work)
test "$(realpath .)" = "$PWD" || die "The tool realpath does not work as expected!"
use_gnu_tool xargs
use_gnu_tool find

# allow injecting commands like 'set -x' to log executed statements
if [ -n "$OGDF_UTILS_PREQUEL" ]
then
	eval "$OGDF_UTILS_PREQUEL"
fi
