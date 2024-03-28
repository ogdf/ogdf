#!/bin/bash
#
# Analyzes the code quality.
#
# required tools:
#  - clang-tidy
#
# optional tools:
#  - clang-tidy-cache (see https://github.com/matus-chochlik/ctcache)
#
# Author: Simon D. Fink

. util/util-functions.sh || exit 123

set -e
mkdir -p build-static-analysis static-analysis

git ls-files --full-name -- 'src/ogdf/*.cpp' 'test/src/*.cpp' | check_filter | sort | uniq > static-analysis/all-sources.txt
echo "$(cat static-analysis/all-sources.txt | wc -l) source files"
git ls-files --full-name -- 'include/ogdf/*.h' | check_filter | sort | uniq > static-analysis/all-headers.txt
echo "$(cat static-analysis/all-headers.txt | wc -l) header files"

## generate compile_commands.json
opts="cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_COMPILER='clang++' "

## debug config
opts+="-DCMAKE_CXX_FLAGS_DEBUG='-g -O1' "
opts+="-DBUILD_SHARED_LIBS=ON "
opts+="-DCMAKE_BUILD_TYPE=Debug "
opts+="-DOGDF_DEBUG_MODE=HEAVY "
opts+="-DOGDF_USE_ASSERT_EXCEPTIONS=OFF "
opts+="-DOGDF_MEMORY_MANAGER=POOL_TS "

## include CGAL
opts+="-DOGDF_INCLUDE_CGAL=ON "
opts+="-DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=TRUE "

## cmd-line args
opts+="$@"

# Compile!
echo "::group::($(date -Iseconds)) Compile"
cd build-static-analysis
echo $opts
eval $opts
make -j "$cores" build-all
cp compile_commands.json ..
cp compile_commands.json ../static-analysis
cat CMakeFiles/OGDF.dir/depend.make \
  | cut -d ":" -f 2 \
  | grep "include/ogdf" \
  | $OGDF_XARGS echo -n \
  | $OGDF_XARGS realpath --relative-to=".." -q \
  | check_filter | sort -u \
  > ../static-analysis/used-headers.txt
cd ..
echo "::endgroup::"

# Check for unused headers
echo "$(cat static-analysis/used-headers.txt | wc -l) used header files"
if [ ! -s static-analysis/used-headers.txt ]; then
  echo "error: no used headers found"
  exit 1
fi
diff static-analysis/used-headers.txt static-analysis/all-headers.txt | sed -ne 's/^> //p' > static-analysis/unused-headers.txt
echo "$(cat static-analysis/unused-headers.txt | wc -l) unused header files"

# Run clang-tidy
echo "::group::($(date -Iseconds)) Run clang-tidy"
if command -v clang-tidy-cache > /dev/null 2>&1; then
  clang_tidy_command="$(which clang-tidy-cache) $(which clang-tidy)"
else
  clang_tidy_command="$(which clang-tidy)"
fi
echo "Using \"$clang_tidy_command\" on $cores cores"
type -a clang-tidy
cat static-analysis/unused-headers.txt static-analysis/all-sources.txt \
  | $OGDF_XARGS -P "$cores" -L 10 \
      $clang_tidy_command -p . --quiet \
  > static-analysis/clang-tidy.txt
echo "::endgroup::"

# TODO Run cppcheck
# cppcheck -j "$cores" --enable=style $compilerincludes --force --xml --language=c++ $sourcefiles $unusedheaders
