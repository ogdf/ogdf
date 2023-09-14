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

# Initialize Dockerfile.
echo "FROM $2" > Dockerfile
cat Dockerfile.template >> Dockerfile

base="$(echo "$2" | cut -d":" -f1)"
version="$(echo "$2" | cut -d":" -f2)"

if [[ "$2" =~ ^{gcc,clang}:[0-9][0-9\.]*$ ]]; then
  echo "Warning: Image not of the form ^{gcc,clang}:[0-9][0-9\.]*$"
  echo "and hence untested. Proceed with caution."
fi

# Build clang image first if needed.
if [[ "$base" = "clang" ]]; then
  simple_version="$(echo $version | cut -d"." -f1)"
  sudo docker build \
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

# Install doxygen only for CLANG_VERSION >= 11.0 (since doxygen 1.9.3 requires
# glibc 2.29 which is not installed in the other tested images).
doxygen_install="false"
if { [[ "$base" = "clang" ]] && is_higher_version "$version" "11"; } ; then
  doxygen_install="true"
fi

# Build and push the image.
image="$1"/"$2"
sudo docker build \
  --no-cache \
  --build-arg "CGAL_INSTALL"="$cgal_install" \
  --build-arg "DOXYGEN_INSTALL"="$doxygen_install" \
  -t "$image" . || die "Failed to build image."
sudo docker push "$image"
