#!/bin/bash
#
# Analyzes coverage.
#
# Author: Tilo Wiedera, Ivo Hedtke

. util/util-functions.sh || exit 123

make_tmpdir $0

regex=''

for f in $(git ls-files include/ogdf/*.h src/ogdf/*.cpp | check_filter)
do
	regex="$regex|$f"
done

regex=${regex:1}

library=`$OGDF_FIND $GUROBI_HOME/lib/libgurobi*.so | head -n 1`
if [ -z "$library" ]; then
  # For Mac and gurobi version >= 8, the library is not .so anymore but .dylib
  library=`$OGDF_FIND $GUROBI_HOME/lib/libgurobi*.dylib | head -n 1`
fi

(cd $tmp && cmake ..) || die "First CMake run failed"
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="--coverage -fPIC -Wall -Wextra -fdiagnostics-show-option -g -O0" \
      -DCMAKE_EXE_LINKER_FLAGS="--coverage" \
      -DOGDF_SEPARATE_TESTS=ON \
      -DOGDF_WARNING_ERRORS=ON \
      -DOGDF_LEAK_CHECK=ON \
      -DOGDF_MEMORY_MANAGER=MALLOC_TS \
      -DOGDF_INCLUDE_CGAL=ON \
      -DCOIN_SOLVER=GRB \
      -DCOIN_EXTERNAL_SOLVER_INCLUDE_DIRECTORIES=$GUROBI_HOME/include \
      -DCOIN_EXTERNAL_SOLVER_LIBRARIES="$library" $tmp || die "Second CMake run failed"
make -C $tmp -j "$cores" build-all || die "Make failed"
util/perform_separate_tests.sh $tmp || die "Tests failed"
echo "LAUNCHING GCOVR..."
gcovr -r . --object-directory=$tmp/CMakeFiles -s -o /dev/null --filter="$regex" > $tmp/gcovr-summary || die "GCOVR failed"

cat $tmp/gcovr-summary | awk '/^lines:/ { lines=$2 } /^branches:/ { branches=$2 } END { print "\nTOTAL COVERAGE IS " (lines+branches)/2 "%" }'

util/run_examples.sh
