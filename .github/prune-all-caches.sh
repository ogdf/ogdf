#!/bin/bash

KEYS=(
  "ccache-self-sufficiency"
  "ccache-coverage"
  "clang-tidy-cache"
  "ccache-static-analysis"
  "ccache-build-linux-gcc:9-debug"
  "ccache-build-linux-gcc:9-release"
  "ccache-build-linux-gcc:13-debug"
  "ccache-build-linux-gcc:13-release"
  "ccache-build-linux-clang:15-debug"
  "ccache-build-linux-clang:15-release"
  "ccache-build-macos-14-debug"
  "ccache-build-macos-14-release"
  "ccache-build-macos-15-debug"
  "ccache-build-macos-15-release"
  "ccache-build-macos-15-intel-debug"
  "ccache-build-macos-15-intel-release"
  "ccache-build-windows-debug"
  "ccache-build-windows-release"
  "ccache-package-debian"
  "ccache-package-macos"
  "ccache-package-windows"
)

cd $(dirname "$0")

# use CLI args or default from env / working dir
REPO="${3:-${GITHUB_REPOSITORY:-$(gh repo view --json nameWithOwner --jq .nameWithOwner)}}"

function usage() {
  gh api "/repos/$REPO/actions/cache/usage" | jq --raw-output \
    '"Repo: " + (.full_name | tostring) +
    "\nCache Count: " + (.active_caches_count | tostring) +
    "\nCache Size: " + (.active_caches_size_in_bytes / 1024 / 1024 / 1024 * 100 | round | . / 100 | tostring) + "GB"'
  echo
}

usage

for ref in $(gh cache list --json ref --jq ".[].ref" | sort | uniq); do
  echo "Pruning ref $ref..."
  for key in ${KEYS[@]}; do
    ./prune-caches.sh "$key" "$ref" "$REPO"
  done
  echo
done

usage
