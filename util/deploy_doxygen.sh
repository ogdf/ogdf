#!/bin/bash
#
# Copy Doc of master to server
#
# Requires the environment variables $DOXYGEN_DEPLOY_URL and $DOXYGEN_DEPLOY_KEY
# to be set correctly.
#
# Author: Ivo Hedtke, Tilo Wiedera, Jöran Schierbaum

. util/util-functions.sh || exit 123

make_tmpdir $0
source_dir=`pwd`

(cd $tmp && cmake -DOGDF_FULL_DOC=OFF $source_dir)
make -s -C $tmp doc

eval $(ssh-agent -s)
ssh-add <(echo "$DOXYGEN_DEPLOY_KEY")

branch=$1
[[ -n "${branch}" ]] || exit 123

echo "deploying latest documentation for ${branch} to ${DOXYGEN_DEPLOY_URL}"
directory="/var/www/html/doc/${branch}/"
cd doc/html && tar czf - * | ssh -o StrictHostKeyChecking=no $DOXYGEN_DEPLOY_URL "mkdir -p ${directory} && cd ${directory} && tar xzf -"
exit $?
