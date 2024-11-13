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
    echo "iwyu requested changes to includes!"
    if [ -z $CI ]; then
      echo "Please apply ./iwyu-stats/iwyu-fixes.patch and call fix_includes.sh again."
    else
      echo "Please download the 'iwyu-stats' artifact from the run Summary page and apply iwyu-fixes.patch."
      echo "  https://github.com/ogdf/ogdf/actions/runs/${GITHUB_RUN_ID}#artifacts"
    fi
    echo
    echo "Changes requested:"
    git -C .. apply --stat iwyu-stats/iwyu-fixes.patch
    echo
    echo "If any of the changes requested are actually incorrect, see the following link for ways of fixing this:"
    echo "  https://github.com/include-what-you-use/include-what-you-use/blob/clang_18/docs/IWYUPragmas.md"
    exit 5
else
    echo "No changes requested"
fi

echo "All good!"
