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

trap "rm -rf $tmp/CMakeFiles $tmp/libOGDF.a $tmp/libCOIN.a" EXIT

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
export CCACHE_BASEDIR="$tmp"
export CCACHE_NOHASHDIR=1

# CMake config according to the arguments
cmakecommand="(cd "$tmp" && cmake -DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=TRUE "
cmakecommand+="-DOGDF_SEPARATE_TESTS=OFF -DOGDF_WARNING_ERRORS=ON "
case "$libtype" in
static)
	;;
dynamic)
	cmakecommand+="-DBUILD_SHARED_LIBS=1 "
	;;
*)
	usage
esac

case "$buildtype" in
release)
	cmakecommand+="-DCMAKE_BUILD_TYPE=Release "
	;;
debug)
	cmakecommand+="-DCMAKE_CXX_FLAGS_DEBUG='-g -O1' "
	cmakecommand+="-DCMAKE_BUILD_TYPE=Debug "
	cmakecommand+="-DOGDF_DEBUG_MODE=HEAVY "
	cmakecommand+="-DOGDF_LEAK_CHECK=ON "
	cmakecommand+="-DOGDF_MEMORY_MANAGER=MALLOC_TS "
	;;
*)
	usage
esac

case "$compilertype" in
gcc)
	cmakecommand+="-DCMAKE_CXX_COMPILER='g++' "
	;;
clang)
	cmakecommand+="-DCMAKE_CXX_COMPILER='clang++' "
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
  cmakecommand+="-DCOIN_SOLVER=GRB -DCOIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES=$GUROBI_HOME/include -DCOIN_EXTERNAL_SOLVER_LIBRARIES=$library "
fi

cmakecommand+="$sourcedir ${@:7})"

echo "::group::($(date -Iseconds)) Initial run of cmake"
echo $cmakecommand
eval $cmakecommand || exit 1
echo "::endgroup::"

run_cmake() {
  echo cmake $@
  cmake "$@" "$tmp"
}

compile () {
	make -C $tmp -j "$cores" build-all | grep -v 'Building CXX object'
}

# build
echo "::group::($(date -Iseconds)) First compile with all custom macros set"
echo "running make using $cores parallel jobs"
ogdf_flags="$(cmake -LA "$tmp" | grep OGDF_EXTRA_CXX_FLAGS:STRING)"
run_cmake "-DOGDF_WARNING_ERRORS=ON" "-D$ogdf_flags $(./util/get_macro_defs.sh)"
compile || exit 1
echo "::endgroup::"

echo "::group::($(date -Iseconds)) Now recompile without custom macros"
run_cmake "-DOGDF_WARNING_ERRORS=ON" "-D$ogdf_flags"
compile || exit 1
echo "::endgroup::"

echo "::group::($(date -Iseconds)) Now recompile tests as separate tests"
run_cmake "-DOGDF_WARNING_ERRORS=ON" "-DOGDF_SEPARATE_TESTS=ON"
compile || exit 1
echo "::endgroup::"

exit $?
