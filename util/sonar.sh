#!/bin/bash
#
# Analyzes the code quality and sends a report to sonarqube.
#
# required tools:
#  - clang-tidy
#  - gcovr
#
# Author: Tilo Wiedera, Simon D. Fink

. util/util-functions.sh || exit 123

set -e
mkdir -p build-sonar sonar

git ls-files --full-name -- 'src/ogdf/*.cpp' 'test/src/*.cpp' | check_filter | sort | uniq > sonar/all-sources.txt
git ls-files --full-name -- 'include/ogdf/*.h' | check_filter | sort | uniq > sonar/all-headers.txt

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
opts+="-DCMAKE_CXX_FLAGS='--coverage -Wall -Wextra -fdiagnostics-show-option' "
opts+="-DCMAKE_EXE_LINKER_FLAGS='--coverage' "

## cmd-line args
opts+="$@"

# Compile!
echo "::group::($(date -Iseconds)) Compile"
cd build-sonar
echo $opts
eval $opts
make -j "$cores" build-all
cp compile_commands.json ..
cp compile_commands.json ../sonar
cat CMakeFiles/OGDF.dir/depend.make \
  | cut -d ":" -f 2 \
  | grep "^ ../include/ogdf/" \
  | $OGDF_XARGS realpath --relative-to=".." -q \
  | check_filter | sort | uniq \
  > ../sonar/used-headers.txt
cd ..
echo "::endgroup::"

# Check for unused headers
diff sonar/used-headers.txt sonar/all-headers.txt | sed -ne 's/^> //p' > sonar/unused-headers.txt

# Run tests and report coverage
util/run_examples.sh
util/perform_separate_tests.sh build-sonar
echo "::group::($(date -Iseconds)) Collect coverage"
type -a gcov # see also gcovr --gcov-executable; needs to match compiler (here: clang -> `llvm-cov-15 gcov`)
gcovr --object-directory=build-sonar/CMakeFiles --root=. --sonarqube=sonar/gcovr.xml
echo "::endgroup::"

# Run clang-tidy
echo "::group::($(date -Iseconds)) Run clang-tidy"
echo "Using clang-tidy in $(which clang-tidy) on $cores cores"
type -a clang-tidy
cat sonar/unused-headers.txt sonar/all-sources.txt \
  | $OGDF_XARGS -P "$cores" -L 10 \
      clang-tidy -p . --quiet \
  > sonar/clang-tidy.txt
echo "::endgroup::"

# TODO Run cppcheck
# cppcheck -j "$cores" --enable=style $compilerincludes --force --xml --language=c++ $sourcefiles $unusedheaders
