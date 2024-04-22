#!/bin/bash

usage() {
  echo "Invocation: ./build.sh <ci-registry> <tag of base image>"
}

die() {
  echo "$@"
  exit -1
}

is_higher_version() {
  highest_version="$(printf "$1\n$2" | sort -V | tail -n 1)"
  [[ "$highest_version" == "$1" ]]
}

[[ -n "$1" ]] || die "$(usage)"
[[ -n "$2" ]] || die "$(usage)"

base="$(echo "$2" | cut -d":" -f1)"
version="$(echo "$2" | cut -d":" -f2)"

if [[ "$2" =~ ^{gcc,clang}:[0-9][0-9\.]*$ ]]; then
  echo "Warning: Image not of the form ^{gcc,clang}:[0-9][0-9\.]*$"
  echo "and hence untested. Proceed with caution."
fi

if [ -z "$DOCKER_CMD" ]; then
  # Use podman if available
  if command -v podman > /dev/null 2>&1; then
    DOCKER_CMD="podman"

  else
    # If user is not in docker group, use sudo.
    prefix="sudo "
    if groups "$(whoami)" | tr " " "\n" | grep "^docker$" || \
      [[ "$OSTYPE" == "darwin"* ]]; then
      prefix=""
    fi
    DOCKER_CMD="${prefix}docker"
  fi
fi

if [ -z "$DOCKER_BUILD_CMD" ]; then
  DOCKER_BUILD_CMD="$DOCKER_CMD build"
fi

if [ -z "$DOCKER_PUSH_CMD" ]; then
  DOCKER_PUSH_CMD="$DOCKER_CMD push"
fi

# Build clang image first if needed.
if [[ "$base" = "clang" ]]; then
  simple_version="$(echo $version | cut -d"." -f1)"
  $DOCKER_BUILD_CMD \
    --build-arg "llvmver"="$simple_version" \
    -t $2 clang/ || die "Failed to build original clang image."
fi

# Install CGAL only for GCC_VERSION >= 6.3 or CLANG_VERSION >= 10.0 (these are
# the minimum supported versions for CGAL 5).
cgal_install="false"
if { [[ "$base" = "gcc" ]]   && is_higher_version "$version" "6.3"; } || \
   { [[ "$base" = "clang" ]] && is_higher_version "$version" "10"; } ; then
  cgal_install="true"
fi

# Build and push the image.
image="$1"/"$2"
$DOCKER_BUILD_CMD \
  --build-arg "compiler"="$base" \
  --build-arg "version"="$version" \
  --build-arg "CGAL_INSTALL"="$cgal_install" \
  -t "$image" . || die "Failed to build image."
$DOCKER_PUSH_CMD "$image"
