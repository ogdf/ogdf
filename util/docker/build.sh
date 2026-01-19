#!/bin/bash

usage() {
  echo "Invocation: ./build.sh <ci-registry> <tag of base image> <debian release name> [ld_lib_path]"
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
registry="$1"
shift

[[ -n "$1" ]] || die "$(usage)"
if [[ "$1" =~ ^{gcc,clang}:[0-9][0-9\.]*$ ]]; then
  echo "Warning: Image not of the form ^{gcc,clang}:[0-9][0-9\.]*$"
  echo "and hence untested. Proceed with caution."
fi
base="$(echo "$1" | cut -d":" -f1)"
version="$(echo "$1" | cut -d":" -f2)"
shift

[[ -n "$1" ]] || die "$(usage)"
release="$1"
shift

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

# Install CGAL 5/6 only for sufficiently recent compiler versions
cgal_install="false"
if { [[ "$base" = "gcc" ]]   && is_higher_version "$version" "11.4"; } || \
   { [[ "$base" = "clang" ]] && is_higher_version "$version" "15.0.7"; } ; then
  cgal_install="6.0.2"
elif { [[ "$base" = "gcc" ]]   && is_higher_version "$version" "6.3"; } || \
   { [[ "$base" = "clang" ]] && is_higher_version "$version" "10"; } ; then
  cgal_install="5.6.3"
fi
echo "Will install CGAL: $cgal_install"

if [[ "$base" = "clang" ]]; then
  # build clang image (with given debian version) first, result will be cached
  $DOCKER_BUILD_CMD \
    --build-arg "llvmver"="$version" \
    --build-arg "release"="$release" \
    -t "$base:$version-$release" clang/ || die "Failed to build original clang image."
fi


# Build and push the image.
image="$registry/$base:$version"
$DOCKER_BUILD_CMD \
  --build-arg "compiler"="$base" \
  --build-arg "version"="$version-$release" \
  --build-arg "ld_lib_path"="$@" \
  --build-arg "CGAL_INSTALL"="$cgal_install" \
  -t "$image" . || die "Failed to build image."
$DOCKER_PUSH_CMD "$image"
