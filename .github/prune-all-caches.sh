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
  "ccache-build-macos-13-debug"
  "ccache-build-macos-13-release"
  "ccache-build-macos-14-debug"
  "ccache-build-macos-14-release"
  "ccache-build-windows-debug"
  "ccache-build-windows-release"
)

cd $(dirname "$0")

for ref in $(gh cache list --json ref --jq ".[].ref" | sort | uniq); do
  echo "Pruning ref $ref..."
  for key in ${KEYS[@]}; do
    ./prune-caches.sh "$key" "$ref"
  done
done
