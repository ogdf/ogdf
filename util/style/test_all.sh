#!/bin/bash
#
# Run all style checks
#
# Author: Simon D. Fink

prog_args="$@"

run_test() {
  echo "Running $1"
  util/style/$1 "$prog_args"
  status=$?
  if [ $status -ne 0 ]; then
    echo "$1 failed with exit code $status, aborting."
    exit 1
  fi
}

run_test test_clang_format.sh
run_test test_eols.sh
run_test test_indent_comments.sh
run_test test_license.sh
run_test test_line_length.sh
run_test test_macros.sh
run_test test_no_enums.sh
run_test test_no_typedefs.sh
run_test test_number_newlines.sh
run_test test_tabs.sh
echo "All checks completed!"
