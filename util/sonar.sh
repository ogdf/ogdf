#!/bin/bash
#
# Analyzes the code quality and sends a report to sonarqube.
#
# Requires the environment variable $SONAR_URL to be set correctly.
#
# required tools:
#  - sonar-runner
#  - cppcheck
#  - cmake
#  - gcc
#  - gcovr
#
# Author: Tilo Wiedera

. util/util-functions.sh || exit 123

make_tmpdir $0

mkdir -p log
trap "rm -rf sonar-project.properties;" EXIT
includes="include, test/include, $tmp/include"
compilerincludes="$(echo "-I$includes" | sed 's/, / -I/g')"

for dir in $(echo | gcc -Wp,-v -x c++ - -fsyntax-only 2>&1 | awk '/#include <...> search starts here:/{f=1;next} /End of search list./{f=0} f')
do
	includes="$includes, $dir"
done

version=$(git rev-parse --short HEAD)
filename="sonar-project.properties"
cp util/$filename.template $filename
echo "sonar.projectVersion=$version" >> $filename
echo "sonar.cxx.includeDirectories=$includes" >> $filename
echo "sonar.password=$SONAR_PASSWORD" >> $filename
echo "sonar.login=$SONAR_ACCESS_TOKEN" >> $filename
echo "sonar.host.url=$SONAR_URL" >> $filename

(cd $tmp && cmake -DCMAKE_BUILD_TYPE=Debug \
                  -DCMAKE_CXX_FLAGS="--coverage -Wall -Wextra -fdiagnostics-show-option" \
                  -DCMAKE_EXE_LINKER_FLAGS="--coverage" \
                  -DOGDF_SEPARATE_TESTS=ON \
                  -DOGDF_WARNING_ERRORS=ON ..) || die "Could not run CMake"
make -C $tmp -j "$cores" COIN 2>/dev/null || die "make COIN failed"
make -C $tmp -j "$cores" build-all 2>log/gcc.log || die "make build-all failed"
util/perform_separate_tests.sh $tmp || die "performing separate tests failed"
gcovr -x -r . --object-directory=$tmp/CMakeFiles >log/gcovr.xml || die "gcovr failed"

sourcefiles="$(git ls-files 'src/ogdf/*.cpp' 'test/src/*.cpp' | check_filter)"
unusedheaders="$(util/find_unused_headers.sh $tmp | check_filter)"
cppcheck -j "$cores" --enable=style $compilerincludes --force --xml --language=c++ $sourcefiles $unusedheaders 2>log/cppcheck.xml || die "cppcheck failed"

# point sonar-scanner to certificate so it can communicate with sonarqube server over https
sonar_keystore="$(echo /sonar-scanner*/jre/lib/security/cacerts)"
export SONAR_SCANNER_OPTS=-Djavax.net.ssl.trustStore="$sonar_keystore"
sonar-scanner || die "sonar-runner failed"
