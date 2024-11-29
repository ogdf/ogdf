#!/bin/bash
#
# Fixes #includes using the include-what-you-use (iwyu) tool.
# Warning: this script will make changes to your ogdf directory.
#
# required tools:
#  - include-what-you-use
#  - iwyu_tool.py (from iwyu repo)
#  - fix_includes.py ( - " - )
#
# Author: Simon D. Fink

. util/util-functions.sh || exit 123
set -eo pipefail

readonly CLANG_VERSION=18
readonly DOCKER_IMAGE="ogdf/clang:$CLANG_VERSION"

usage () {
    echo
    echo "$0 [FILES]"
    echo
    echo "Fixes all import statements using include-what-you-use."
    echo "If include-what-you-use is not found, a docker image is used instead."
    echo
    echo "WARNING: This will modify your files!."
    echo
    echo "If [FILES] are given, run iwyu only on these files or directories."
    echo "Note that the list will not be deduplicated, iwyu will be run twice if a file/directory is given twice."
    echo "For header files in include/ogdf, we create a dummy cpp file in src/ogdf/iwyu-dummy for iwyu to process."
}

original_args="$@"
while getopts "n" opt; do
  case $opt in
    n) nosudo=yes ;;
    *) usage && exit 255 ;;
  esac
done
shift $(($OPTIND-1))


# Check whether the correct iwyu is available. If not, try to find and use a container runtime.
if ! (include-what-you-use --version | grep "include-what-you-use 0.22 .* clang version 18") > /dev/null 2>&1; then
  if [ -z "$DOCKER_RUN_CMD" ]; then
    # Use podman if available
    if command -v podman > /dev/null 2>&1; then
      DOCKER_RUN_CMD="podman run"
      # (changing the user ID will lead to unsafe git repo location warnings)

    else
      # If user is not in docker group, use sudo.
      prefix="sudo "
      user="$(whoami)"
      if [ -n "$nosudo" ] || \
        groups "$user" | tr " " "\n" | grep "^docker$" || \
        [[ "$OSTYPE" == "darwin"* ]]; then
        prefix=""
      fi

      DOCKER_RUN_CMD="${prefix}docker run --user $(id -u "$user")"
    fi
  fi

  echo "include-what-you-use not found in a suitable version. Using container (via \`$DOCKER_RUN_CMD\`) instead."

  # Run this script in a docker container.
  repo_dir="$(pwd)"
  $DOCKER_RUN_CMD --rm -ti -w "$repo_dir" \
    -v "$repo_dir":"$repo_dir":rw,z \
    "$DOCKER_IMAGE" \
    "$0" $original_args
  exit_code=$?
  exit $exit_code
fi


rm -rf build-iwyu src/ogdf/iwyu-dummy iwyu-stats
mkdir -p build-iwyu src/ogdf/iwyu-dummy iwyu-stats

# iwyu runs on cpp files and only indirectly on headers, so ensure we have a .cpp for each .h
echo "Create associated cpp files for all headers"
find include/ogdf -name '*.h' | while read hfile; do
    hfile="${hfile/include\//}"
    cfile="src/${hfile/ogdf\//ogdf\/iwyu-dummy\/}"
    cfile="${cfile/.h/.cpp}"
    mkdir -p $(dirname $cfile)
    echo -e "#include <${hfile}> // IWYU pragma: associated \n" > $cfile
done
rm -rf src/ogdf/iwyu-dummy/lib

SRC_DIRS=("$@")
if [ ${#SRC_DIRS[@]} -eq 0 ]; then
  SRC_DIRS=("src/ogdf" "doc/examples" "test/src")
fi
SRC_DIRS=( "${SRC_DIRS[@]//include\/ogdf/src\/ogdf\/iwyu-dummy}" )
SRC_DIRS=( "${SRC_DIRS[@]//.h/.cpp}" )
echo "Will process the following ${#SRC_DIRS[@]} files/directories:"
printf "\t%s\n" "${SRC_DIRS[@]}"

readarray -t SRC_DIRS < <(for d in "${SRC_DIRS[@]}"; do
  if ! realpath --canonicalize-existing "$d"; then
    echo "One of the given paths does not exist!" >&2
    usage >&2
    rm -rf src/ogdf/iwyu-dummy
    exit 255
  fi
  d="${d//src\/ogdf/src\/ogdf\/iwyu-dummy}"
  if [ -e "$d" ]; then
      realpath "$d"
  fi
done | sort -u | awk 'a=="" || $0 !~ "^"a {a=$0; print}') # deduplicated https://stackoverflow.com/a/48672202
echo "Resolved sources to ${#SRC_DIRS[@]} files/directories:"
printf "\t%s\n" "${SRC_DIRS[@]}"

readarray -t SRC_HDR_DIRS < <(for d in "${SRC_DIRS[@]}"; do
  ARR=("$d" "${d//src\/ogdf\/iwyu-dummy/include\/ogdf}" "${d//src\/ogdf/include\/ogdf}")
  for dd in "${ARR[@]}"; do
    if [ -e "$dd" ]; then
      echo "$dd"
    fi
  done
done | sort -u | awk 'a=="" || $0 !~ "^"a {a=$0; print}')
echo "Resolved sources+headers to ${#SRC_HDR_DIRS[@]} files/directories:"
printf "\t%s\n" "${SRC_HDR_DIRS[@]}"


replace_all () {
  # replaces all occurrences of $1 with $2 if they are not followed by a comment,
  # so appending eg `// IWYU pragma: keep` prevents replacement
  grep "$1" -r --files-with-matches "${SRC_HDR_DIRS[@]}" | $OGDF_XARGS --no-run-if-empty sed -i "s|$1$|$2|g" || true
  if [ -n "$2" ]; then
    # only keep the first occurrence of $2 in a file
    grep "$2" -r --files-with-matches "${SRC_HDR_DIRS[@]}" | $OGDF_XARGS --no-run-if-empty sed -z -i "s|$2$||2g" || true
  fi
}

run_all_replacements () {
  replace_all   "#include <ogdf/basic/internal/config.h>"        "#include <ogdf/basic/basic.h>"
  replace_all   "#include <ogdf/basic/internal/config_autogen.h" "#include <ogdf/basic/basic.h>"
  replace_all   "#include <ogdf/basic/Graph_d.h>"   "#include <ogdf/basic/Graph.h>"
  replace_all   "#include <ogdf/lib/abacus/.*>"     "#include <ogdf/external/abacus.h>"
  replace_all   "#include <coin/.*>"                "#include <ogdf/external/coin.h>"
  replace_all   "#include <bandit/.*>"              "#include <testing.h>"

  replace_all   "#include <emmintrin.h>"    ""
  replace_all   "#include <pmmintrin.h>"    ""
  replace_all   "#include <xmmintrin.h>"    ""
  replace_all   "#include <bits/chrono.h>"  "#include <chrono>"
  replace_all   "#include <bits/.*>"        ""
  replace_all   "#include <built-in>"       ""

  # https://cplusplus.com/reference/clibrary/
  replace_all   "#include <assert.h>"   "#include <cassert>"
  replace_all   "#include <ctype.h>"    "#include <cctype>"
  replace_all   "#include <errno.h>"    "#include <cerrno>"
  replace_all   "#include <float.h>"    "#include <cfloat>"
  replace_all   "#include <iso646.h>"   "#include <ciso646>"
  replace_all   "#include <limits.h>"   "#include <climits>"
  replace_all   "#include <locale.h>"   "#include <clocale>"
  replace_all   "#include <math.h>"     "#include <cmath>"
  replace_all   "#include <setjmp.h>"   "#include <csetjmp>"
  replace_all   "#include <signal.h>"   "#include <csignal>"
  replace_all   "#include <stdarg.h>"   "#include <cstdarg>"
  replace_all   "#include <stdbool.h>"  "#include <cstdbool>"
  replace_all   "#include <stddef.h>"   "#include <cstddef>"
  replace_all   "#include <stdint.h>"   "#include <cstdint>"
  replace_all   "#include <stdio.h>"    "#include <cstdio>"
  replace_all   "#include <stdlib.h>"   "#include <cstdlib>"
  replace_all   "#include <string.h>"   "#include <cstring>"
  replace_all   "#include <time.h>"     "#include <ctime>"
  replace_all   "#include <uchar.h>"    "#include <cuchar>"
  replace_all   "#include <wchar.h>"    "#include <cwchar>"
  replace_all   "#include <wctype.h>"   "#include <cwctype>"
}


# generate compile_commands.json
echo "::group::($(date -Iseconds)) Run cmake"
echo "First run of cmake"
cmake -S . -B build-iwyu -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_COMPILER='clang++' \
  -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG='-g -O1' \
  -DOGDF_DEBUG_MODE=HEAVY -DOGDF_USE_ASSERT_EXCEPTIONS=ON -DOGDF_MEMORY_MANAGER=POOL_TS -DOGDF_SEPARATE_TESTS=ON
# the optional CGAL includes tend to confuse iwyu, so disable for now
# -DOGDF_INCLUDE_CGAL=ON -DCGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE=TRUE

# set all custom macros
echo "Update cmake with all flags set"
ogdf_flags="$(cmake -LA build-iwyu | grep OGDF_EXTRA_CXX_FLAGS:STRING)"
cmake "-DOGDF_WARNING_ERRORS=ON" "-D$ogdf_flags $(./util/get_macro_defs.sh)" build-iwyu

# ensure all ogdf includes are found via `-isystem` instead of -I to make iwyu use angle brackets instead of quotes
sed -i "s|-I$(pwd)/build-iwyu/include|-isystem $(pwd)/build-iwyu/include|g" build-iwyu/compile_commands.json
sed -i "s|-I$(pwd)/include|-isystem $(pwd)/include|g" build-iwyu/compile_commands.json
sed -i "s|-I$(pwd)/test/include|-isystem $(pwd)/test/include|g" build-iwyu/compile_commands.json
sed -i "s|-isystem $(pwd)/include/coin||g" build-iwyu/compile_commands.json
echo "::endgroup::"

IWYU_OPTIONS="-Xiwyu --cxx17ns -Xiwyu --max_line_length=100 -Xiwyu --error=0"
FIX_OPTIONS="--blank_lines --nosafe_headers --reorder --ignore_re 'lib/|coin/|iwyu-dummy/' --separate_project_includes=ogdf"

# enable "why" comments
#IWYU_OPTIONS+=" -Xiwyu --update_comments"
#FIX_OPTIONS+=" --comments --update_comments"
# disable "why" comments
#IWYU_OPTIONS+=" -Xiwyu --no_comments"
FIX_OPTIONS+=" --nocomments --noupdate_comments"

# run iwyu on all files
echo "::group::($(date -Iseconds)) Run iwyu"
iwyu_ret=0
set -x
iwyu_tool.py -p build-iwyu/ -j $cores "${SRC_DIRS[@]}" -- $IWYU_OPTIONS | tee iwyu-stats/iwyu-out.txt || iwyu_ret=$?
set +x
echo "iwyu_tool exited with code $iwyu_ret"
echo "::endgroup::"

# collect some stats
echo "Collect stats"
pushd iwyu-stats
set +e
grep -Po '(?<=\()(.*)(?= has correct)' iwyu-out.txt > files-correct.txt
grep -Po '(?<=The full include-list for )(.*)(?=:)' iwyu-out.txt > files-to-fix.txt
cat files-correct.txt files-to-fix.txt | sort | uniq -c | egrep -v '^\s+1 ' > files-double.txt
find "${SRC_HDR_DIRS[@]}" -type f -regex ".*\.[ch]\(pp\)?$" -exec realpath \{\} \; | sort -u > files-all.txt
comm -3 <(cat files-all.txt) <(sort -u files-correct.txt files-to-fix.txt) | egrep -v 'ogdf/lib' > files-missing.txt
grep " error: " -A 2 -B 1 iwyu-out.txt > iwyu-error.txt
set -e
popd

# apply all iwyu fixes
echo "::group::($(date -Iseconds)) Run fix_includes"
fix_ret=0
set -x
cat iwyu-stats/iwyu-out.txt | fix_includes.py $FIX_OPTIONS || fix_ret=$?
set +x
echo "fix_includes exited with code $fix_ret"
echo "::endgroup::"
echo "::group::($(date -Iseconds)) Clean-up"
run_all_replacements
git restore src/ogdf/lib include/ogdf/lib src/ogdf/external include/ogdf/external
rm -rf src/ogdf/iwyu-dummy
util/style/test_clang_format.sh -f
echo "::endgroup::"

# generate patch for all changes
git diff > iwyu-stats/iwyu-fixes.patch
git status

if [ $iwyu_ret != 0 ]; then
  echo "Warning: iwyu_tool exited with code $iwyu_ret"
  exit $iwyu_ret
fi

if [ $fix_ret != 0 ]; then
  echo "Warning: fix_includes exited with code $fix_ret"
  exit $ret
fi

# To test:
# util/test_includes.sh
# pushd build-iwyu && rm -rf ../src/ogdf/iwyu-dummy && cmake .. && make -j 12 build-all; popd
# util/test_self-sufficiency.sh
# util/style/test_all.sh -f

# To commit:
# rm -rf build-iwyu src/ogdf/iwyu-dummy iwyu-stats
# git add src/ogdf/ include/ogdf/
# git commit -m "iwyu fixes"

# To check for weird-looking iwyu changes:
# git diff old new -U0 | \
# 	grep -Ev "^(@@|index|---|\+\+\+) |^[+-]\s*(#\s*include <.*>|namespace|class|enum|struct|template|\})|^[+-]\s*$" | \
# 	tac | awk ' a !~ "^diff --git" || $0 !~ "^diff --git" {a=$0; print; if($0 ~ "^diff --git"){print "\n"}}' | tac
