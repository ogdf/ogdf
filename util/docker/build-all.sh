#!/bin/bash

set -e
set -x

# on ggc:9-bookworm, apt needs to be pointed to the newer system libstdc++ via LD_LIBRARY_PATH
./build.sh "$@" gcc:9    bookworm --env="ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/"
./build.sh "$@" gcc:13   bookworm # here our gcc 13 actually brings newer libstdc++ than system gcc 12
./build.sh "$@" clang:15 bookworm
./build.sh "$@" clang:18 bookworm
