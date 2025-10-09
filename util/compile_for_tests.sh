#!/bin/bash
#
# Compiles the OGDF (using an in-source/out-of-source build).
#
# Author: Tilo Wiedera, Ivo Hedtke

. util/util-functions.sh || exit 123

usage () {
	cat >&2<<EOF
Usage: $0 <1st> <2nd> <3rd> <4th> <sourcedir> <builddir> [cmake arguments]

First argument:
   static    Static library build
   dynamic   Dynamic (shared) library build

Second argument:
   debug     Debug build
   release   Release build

Third argument:
   default_c Use the default compiler
   clang     Use clang++ (llvm-clang)
   gcc       Use g++ (gcc)

Fourth argument:
   default_s Use the default ILP solver
   gurobi    Use Gurobi as the ILP solver

Further arguments get passed to the CMake invocation.

Warning: Running this script may modify your build in <builddir>.
EOF
	exit 2
}

if [ "$#" -lt "6" ]
then
	usage
fi

set -o pipefail

if [ -z "$OGDF_KEEP_TMP" ]; then
	trap "rm -rf $tmp/CMakeFiles $tmp/libOGDF.a $tmp/libCOIN.a" EXIT
fi

#detect OS
unamestr=`uname`

# build and source paths
libtype=$1
buildtype=$2
compilertype=$3
ilpsolvertype=$4
tmp=`realpath $6`
sourcedir=`realpath $5`

mkdir -p $tmp
export CCACHE_BASEDIR="$sourcedir"

# CMake config according to the arguments
cmaketestargs=()
cmakeargs=()
cmakeargs+=("-DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=TRUE")
cmakeargs+=("-DOGDF_SEPARATE_TESTS=OFF" "-DOGDF_WARNING_ERRORS=ON")
case "$libtype" in
static)
	;;
dynamic)
	cmakeargs+=("-DBUILD_SHARED_LIBS=1")
	cmaketestargs+=("-DBUILD_SHARED_LIBS=1")
	;;
*)
	usage
esac

case "$buildtype" in
release)
	cmakeargs+=("-DCMAKE_BUILD_TYPE=Release")
	cmaketestargs+=("-DCMAKE_BUILD_TYPE=Release")
	;;
debug)
	cmakeargs+=("-DCMAKE_CXX_FLAGS_DEBUG='-g -O1'")
	cmakeargs+=("-DCMAKE_BUILD_TYPE=Debug")
	cmaketestargs+=("-DCMAKE_BUILD_TYPE=Debug")
	cmakeargs+=("-DOGDF_DEBUG_MODE=HEAVY")
	cmakeargs+=("-DOGDF_LEAK_CHECK=ON")
	cmakeargs+=("-DOGDF_MEMORY_MANAGER=MALLOC_TS")
	;;
*)
	usage
esac

case "$compilertype" in
gcc)
	cmakeargs+=("-DCMAKE_CXX_COMPILER='g++'")
	cmaketestargs+=("-DCMAKE_CXX_COMPILER='g++'")
	;;
clang)
	cmakeargs+=("-DCMAKE_CXX_COMPILER='clang++'")
	cmaketestargs+=("-DCMAKE_CXX_COMPILER='clang++'")
	;;
default_c|*)
esac

# Enable Gurobi if wanted.
if [ "$ilpsolvertype" = "gurobi" ]; then
  library=`$OGDF_FIND $GUROBI_HOME/lib/libgurobi*.so | head -n 1`
  if [ -z "$library" ]; then
    # For Mac and gurobi version >= 8, the library is not .so anymore but .dylib
    library=`$OGDF_FIND $GUROBI_HOME/lib/libgurobi*.dylib | head -n 1`
  fi
  cmakeargs+=("-DCOIN_SOLVER=GRB" "-DCOIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES='$GUROBI_HOME/include'" "-DCOIN_EXTERNAL_SOLVER_LIBRARIES='$library'")
fi

run_cmake() {
	echo cmake $@
	cmake "$@" -B "$tmp" -S "$sourcedir"
	ret=$?
	if [ $ret != 0 ]; then
		echo "CMake failed with exit code $ret"
		exit $ret
	fi
}

compile () {
	cmake --build $tmp --parallel "$cores" | grep -v 'Building CXX object'
	ret=$?
	if [ $ret != 0 ]; then
		echo "Make failed with exit code $ret"
		exit $ret
	fi
}

echo "::group::($(date -Iseconds)) Initial run of cmake"
run_cmake "${cmakeargs[@]}" "${@:7}"
echo "::endgroup::"

# build
echo "::group::($(date -Iseconds)) First compile with all custom macros set"
echo "running make using $cores parallel jobs"
ogdf_flags="$(cmake -LA "$tmp" | grep OGDF_EXTRA_CXX_FLAGS:STRING)"
run_cmake "-DOGDF_WARNING_ERRORS=ON" "-D$ogdf_flags $(./util/get_macro_defs.sh)"
compile
echo "::endgroup::"

echo "::group::($(date -Iseconds)) Now recompile without custom macros"
run_cmake "-DOGDF_WARNING_ERRORS=ON" "-D$ogdf_flags"
compile
echo "::endgroup::"

echo "::group::($(date -Iseconds)) Now recompile tests as separate tests"
run_cmake "-DOGDF_WARNING_ERRORS=ON" "-DOGDF_SEPARATE_TESTS=ON"
compile
echo "::endgroup::"

echo "::group::($(date -Iseconds)) Now compile simple dependent project"
old_tmp="$tmp"
old_sourcedir="$sourcedir"
tmp="$old_tmp/build-doc-examples-special"
sourcedir="$old_sourcedir/doc/examples/special"
run_cmake "${cmaketestargs[@]}" "-DOGDF_DIR=$old_tmp"
compile
"$tmp/check-build-mode"
ret=$?
if [ $ret != 0 ]; then
  echo "check-build-mode failed with exit code $ret"
  exit $ret
fi
tmp="$old_tmp"
sourcedir="$old_sourcedir"
echo "::endgroup::"
