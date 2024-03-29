#!/bin/bash
#
# Verify the output of fix_includes.sh
#
# Author: Simon D. Fink

set -e -o pipefail

if [ ! -d iwyu-stats ]; then
    echo "Missing directory ./iwyu-stats/! Please call fix_includes.sh first."
    exit 1
fi

cd iwyu-stats

if [ -s iwyu-error.txt ]; then
    echo "iwyu reported errors:"
    cat iwyu-error.txt
    echo
    exit 2
else
    echo "No errors"
fi

if [ -s files-missing.txt ]; then
    echo "iwyu did not process the following files:"
    cat files-missing.txt
    echo
    exit 3
else
    echo "No missing files"
fi

if [ -s files-double.txt ]; then
    echo "Problem with detecting correct and to-fix files! The following files have both states:"
    cat files-double.txt
    echo
    exit 4
else
    echo "No doubles"
fi

if [ -s iwyu-fixes.patch ]; then
    echo "iwyu requested changes to includes! Please download the artifact and apply iwyu-fixes.patch."
    git status
    echo
    exit 5
else
    echo "No changes requested"
fi

echo "All good!"
