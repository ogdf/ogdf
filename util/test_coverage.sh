#!/bin/bash
#
# Analyzes the code test coverage.
#
# required tools:
#  - clang++
#  - llvm-profdata merge
#  - llvm-cov show|export|report
#
# Author: Simon D. Fink

. util/util-functions.sh || exit 123

set -e
mkdir -p build-coverage/profraw coverage

## generate compile_commands.json
opts="cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_COMPILER='clang++' "

## debug config
opts+="-DCMAKE_CXX_FLAGS_DEBUG='-g -O1' "
opts+="-DBUILD_SHARED_LIBS=ON "
opts+="-DCMAKE_BUILD_TYPE=Debug "
opts+="-DOGDF_DEBUG_MODE=HEAVY "
opts+="-DOGDF_USE_ASSERT_EXCEPTIONS=ON "
opts+="-DOGDF_MEMORY_MANAGER=POOL_TS "

## include CGAL
opts+="-DOGDF_INCLUDE_CGAL=ON "
opts+="-DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=TRUE "

## coverage
opts+="-DOGDF_SEPARATE_TESTS=ON "
opts+="-DOGDF_ARCH=x86-64 "
opts+="-DCMAKE_CXX_FLAGS='-fprofile-instr-generate -fcoverage-mapping -femit-all-decls -Wall -Wextra -fdiagnostics-show-option' "
# llvm-cov show for multiple objects requires emit-all-decls, which is incompatible with SSE (thus the old arch)
# https://github.com/llvm/llvm-project/issues/32849#issuecomment-2353071071

## cmd-line args
opts+="$@"

export CCACHE_BASEDIR="$(pwd)"

# Compile!
echo "::group::($(date -Iseconds)) Compile"
cd build-coverage
echo $opts
eval $opts
make -j "$cores" build-all
cp compile_commands.json ..
cd ..
echo "::endgroup::"

# Run tests and report coverage
export LLVM_PROFILE_FILE="$(realpath build-coverage/profraw)/%p.profraw"
util/run_examples.sh
util/perform_separate_tests.sh build-coverage
echo "::group::($(date -Iseconds)) Collect coverage"
tests=$(printf -- "-object %s " build-coverage/test/bin/test-*)
lib="build-coverage/libOGDF.so"
llvm-profdata merge -sparse build-coverage/profraw/*.profraw -o coverage/coverage.profdata
llvm-cov show --format=text $lib $tests -instr-profile=coverage/coverage.profdata > coverage/coverage.txt
# llvm-cov show --format=html $lib $tests -instr-profile=coverage/coverage.profdata > coverage/coverage.html
llvm-cov export $lib $tests -instr-profile=coverage/coverage.profdata > coverage/coverage.json
llvm-cov export --format=lcov $lib $tests -instr-profile=coverage/coverage.profdata > coverage/coverage.lcov
llvm-cov report $lib $tests -instr-profile=coverage/coverage.profdata > coverage/report.txt
echo "::endgroup::"
