#!/bin/bash
#
# Do some steps to generate a release
#
# NO WARRANTY (WARNING: This script contains state-based rm -f and mv -f!)
#
# Author: Max Ilsen, Stephan Beyer

TARGET_DIR=../releases
TMP_DIR="$TARGET_DIR"/tmp
TMP_FILE="tmp-$$"

README=README.md
DOC_HEADER=doc/ogdf-header.html
DOC_FOOTER=doc/ogdf-footer.html
PORTING_DIR=doc/porting
PORTING_FILE=doc/porting.md
RELNOTES_DIR=doc/relnotes
RELNOTES_FILE=doc/relnotes.md
UNREL_PORTING_GUIDE="$PORTING_DIR"/unreleased.md
VERSION_H=include/ogdf/basic/internal/version.h
CMAKELISTS=CMakeLists.txt

GITHUB_REPO="https://github.com/ogdf/ogdf.git"
GITHUB_URL="git@github.com:ogdf/ogdf.git"

MAIN_BRANCH="master"
RELEASE_BRANCH="new-release"
LATEST_RELEASE_BRANCH="latest-release"

usage () {
  echo
  echo "Usage:   $0 <name-of-new-release>"
  echo "Example: $0 chestnut"
  echo
}

die() {
  echo "== An error occured: $*" >&2
  exit 1
}

test_say() {
  yes_text="$1"
  shift
  no_text="$1"
  shift
  if [ "$@" ]
  then
    echo -n "$yes_text"
  else
    echo -n "$no_text"
  fi
}

test_say_file() {
  test_say "GOOD: exists" "BAD: does not exist!" -f "$1"
}

check_or_init_repo () {
  if [ -d "$1" ]
  then
    echo "=== It is there! Doing a pull --rebase on master"
    pushd "$1" || die "Directory $1 does not exist"
    git stash
    git checkout master || die "git checkout master failed"
    git pull --rebase || die "git pull --rebase failed"
    popd
  elif [ -f "$1" ]
  then
    die "There is a file (not a directory) with name $1 ... exiting."
  else
    echo "=== It is not there. I'm cloning."
    git clone "$2" "$1" || die "git clone $2 failed"
  fi
}

separator() {
  echo
  echo "==============================================================================="
  echo
}

yessir() {
  echo ">> Type 'yes' to continue."
  read MAKE_RELEASE_YESSIR_REPLY
  test "$MAKE_RELEASE_YESSIR_REPLY" = "yes" || die "Not a 'yes'"
  separator
}

fixthat() {
  echo "Uh oh, we noticed something unexpected:"
  echo " * $1"
  echo
  echo "Can you fix that?"
  yessir
}

ask_to_commit() {
  git --no-pager diff
  echo
  echo "These are all the changes for the commit \"$1\"."
  echo "Please check whether they are complete. If they are not, you can make additional changes now."
  echo "I will commit these changes once you continue."
  echo ">> Type 'yes' to continue."
  read MAKE_RELEASE_YESSIR_REPLY
  test "$MAKE_RELEASE_YESSIR_REPLY" = "yes" || die "Not a 'yes'"
  git commit --all --message "$1" || die "git commit failed."
  echo
  separator
}

insert_before_last_paragraph() {
  tac "$1" | sed '0,/^$/{s/^$/'"$2"'\n/}' | tac > "$TMP_FILE"
  mv -f "$TMP_FILE" "$1"
}


# initial check: new release name makes sense
RELEASE_LCNAME="$(echo "$@" | tr '[:upper:]' '[:lower:]')"
if test -z "$RELEASE_LCNAME" || echo $RELEASE_LCNAME | grep --quiet -v '^[a-z]*$'
then
  usage
  die "Name of new release should be the name of a tree (a single word, no digits, no special characters)."
fi

# initial check: in root of working tree, no changes made
cd "$(git rev-parse --show-toplevel)" || die "Changing directory to root of working tree failed."

if ! git diff --exit-code >/dev/null 2>&1
then
  echo "You might want to commit, stash or reset your changes."
  die "Your working tree is dirty."
fi

# intial check: on a local branch that is up-to-date with the newest main branch
# in particular, this release script should be in its newest state
git fetch --quiet $GITHUB_URL $MAIN_BRANCH || die "git fetch failed."
test "$(git rev-parse FETCH_HEAD)" = "$(git rev-parse HEAD)" ||
  die "Please run this from a branch that is up-to-date with $MAIN_BRANCH from $GITHUB_URL."

# create and move to release branch if necessary
if ! test "$(git rev-parse --abbrev-ref HEAD)" = "$RELEASE_BRANCH"
then
  git checkout FETCH_HEAD
  git switch -c $RELEASE_BRANCH || die "Creating and switching to branch $RELEASE_BRANCH failed."
fi


cat <<EOF

################################
#   Creating an OGDF release   #
################################

This script leads you through the process of making and releasing a new
OGDF version and assists you with it.

If you read this, we are already on the branch $RELEASE_BRANCH that contains the
up-to-date version of $MAIN_BRANCH on $GITHUB_URL.
A good start!

Now check if the following settings are sane:

The following files and directories are needed for the script to work:
 * Current directory:             "`pwd`"
 * Directory to save release in:  "$TARGET_DIR" [`test_say "GOOD: exists" "GOOD: will be created" -d "$TARGET_DIR"`]
 * Temporary directory name:      "$TMP_DIR" [`test_say "BAD: it exists already" "GOOD: will be created" -d "$TMP_DIR"`]
 * Temporary file name:           "$TMP_FILE" [`test_say "BAD: it exists already" "GOOD: will be created" -f "$TMP_FILE"`]

The following files and directories will be adapted:
 * README file:                   "$README" [`test_say_file "$README"`]
 * Header of OGDF doc:            "$DOC_HEADER" [`test_say_file "$DOC_HEADER"`]
 * Footer of OGDF doc:            "$DOC_FOOTER" [`test_say_file "$DOC_FOOTER"`]
 * Directory with porting guides: "$PORTING_DIR" [`test_say "GOOD: exists" "BAD: does not exist" -d "$PORTING_DIR"`]
 * Porting guide overview:        "$PORTING_FILE" [`test_say_file "$PORTING_FILE"`]
 * Directory with release notes:  "$RELNOTES_DIR" [`test_say "GOOD: exists" "BAD: does not exist" -d "$RELNOTES_DIR"`]
 * Release note overview:         "$RELNOTES_FILE" [`test_say_file "$RELNOTES_FILE"`]
 * Unreleased porting guide:      "$UNREL_PORTING_GUIDE" [`test_say_file "$UNREL_PORTING_GUIDE"`]
 * Version header:                "$VERSION_H" [`test_say_file "$VERSION_H"`]
 * CMakeLists.txt with version:   "$CMAKELISTS" [`test_say_file "$CMAKELISTS"`]

Please also make sure that you have write access to the following repo:
  $GITHUB_REPO

This script will not push or upload anything automatically but you should
be able to do that.

I will now fetch the tags from $GITHUB_URL to compile information about the previous and new release.
Sounds good?
EOF
yessir

### Get names of last and new release
git fetch --quiet --tags $GITHUB_URL || "git fetch failed."
LAST_RELEASE_TAG="$(git tag --list '*-[0-9][0-9][0-9][0-9][0-9][0-9]' | sort -t"-" -k"2" | tail -n1)"
if test -z "$LAST_RELEASE_TAG"
then
  die "Cannot detect name of last release."
fi
LAST_RELEASE_LCNAME="$(echo "$LAST_RELEASE_TAG" | cut -d"-" -f1)"
LAST_RELEASE_NAME="${LAST_RELEASE_LCNAME^}"

RELEASE_NAME="${RELEASE_LCNAME^}"
VERSION_NUMBER="$(date +%Y.%m)"
TAGNAME="$RELEASE_LCNAME-$(echo "$VERSION_NUMBER" | tr -d ".")"

cat <<EOF
I got:
 * Name of the new release:            "$RELEASE_NAME"
 * Version number of the new release:  "$VERSION_NUMBER"
 * Git tag name of the new release:    "$TAGNAME"

 * Name of the last release:           "$LAST_RELEASE_NAME"
 * Git tag name of the last release:   "$LAST_RELEASE_TAG"

Is this fine?
EOF
yessir

while git rev-parse --abbrev-ref "$TAGNAME" >/dev/null 2>&1
do
  fixthat "Tag $TAGNAME already exists! Maybe git tag -d $TAGNAME can help?"
done


### Update copyright year
echo "Updating copyright year in $DOC_FOOTER and $README..."
echo

# careful: the replaced strings contain long hyphens, not just dashes
CURRENT_YEAR=$(date "+%Y")
sed -i 's/\(1999&ndash;\)[0-9]\{4\}/\1'"$CURRENT_YEAR"/g $DOC_FOOTER
sed -i 's/\(1999–\)[0-9]\{4\}/\1'"$CURRENT_YEAR"/g $README

if git diff --exit-code >/dev/null 2>&1
then
  grep 1999 $DOC_FOOTER $README
  echo
  echo "Copyright year is already up to date. Does that look right?"
  yessir
else
  ask_to_commit "Update latest year of copyright to $CURRENT_YEAR"
fi


##### Update porting guides
REL_PORTING_GUIDE="$PORTING_DIR/$RELEASE_LCNAME.md"
if test -f "$REL_PORTING_GUIDE"
then
  fixthat "$REL_PORTING_GUIDE already exists. Expecting release notes in $UNREL_PORTING_GUIDE."
fi

cat >"$REL_PORTING_GUIDE" <<EOF
[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » $RELEASE_NAME

# Porting from $LAST_RELEASE_NAME to $RELEASE_NAME
$(tail -n+4 $UNREL_PORTING_GUIDE)
EOF

cat >"$UNREL_PORTING_GUIDE" <<EOF
[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » Unreleased

# Porting from $RELEASE_NAME to current unreleased version

There are currently no breaking changes.
EOF

sed -i "/from $LAST_RELEASE_NAME/{s/the current unreleased version"'\(.*\)'"unreleased/$RELEASE_NAME\1$RELEASE_LCNAME/}" \
  "$PORTING_FILE"
NEW_PORTING_LINE="  * [from $RELEASE_NAME to the current unreleased version](porting\\/unreleased.md)"
insert_before_last_paragraph "$PORTING_FILE" "$NEW_PORTING_LINE"

git add --intent-to-add "$REL_PORTING_GUIDE"
ask_to_commit "Update porting guides"


##### Update release notes
NEW_REL_LINE="  * [$VERSION_NUMBER ($RELEASE_NAME)](relnotes\\/$RELEASE_LCNAME.md)"
insert_before_last_paragraph "$RELNOTES_FILE" "$NEW_REL_LINE"

# prepare release notes template
git shortlog --no-merges "$LAST_RELEASE_TAG"...HEAD >"$TMP_FILE"

RELNOTES_TEMPLATE="$RELNOTES_DIR/$RELEASE_LCNAME-template.md"
RELNOTES_NEW="$RELNOTES_DIR/$RELEASE_LCNAME.md"
cat >"$RELNOTES_TEMPLATE" <<EOF
[OGDF](../../README.md) » [Release Notes](../relnotes.md) » $RELEASE_NAME

# OGDF $RELEASE_NAME ($VERSION_NUMBER)

Released $(date -I).

This release... [you may add some mentionable notes about the release here
or remove this paragraph]

Noteworthy changes:
$(cat $TMP_FILE | sed -ne 's/^ \{6\}/ * /p')

This release contains contributions by:
$(cat $TMP_FILE | sed -ne 's/^\([^ ].*\) (.*$/ * \1/p')
EOF

rm -f "$TMP_FILE"

cat <<EOF
Please update the release notes in $RELNOTES_NEW
  * I already updated the release notes overview and made a template for you based
    on the git history, see $RELNOTES_TEMPLATE
  * Make sure the template file is gone after this step
EOF
yessir

while test -f "$RELNOTES_TEMPLATE"
do
  fixthat "The template file $RELNOTES_TEMPLATE still exists."
done
while ! test -f "$RELNOTES_NEW"
do
  fixthat "You should have created the release notes \"$RELNOTES_NEW\" but they do not exist."
done

git add --intent-to-add "$RELNOTES_NEW"
ask_to_commit "Add release notes for $RELEASE_NAME ($VERSION_NUMBER)"


##### Change version number in files and create tag
bump_version() {
  new_version="$1"
  shift
  release_name="$1"
  new_cmake_version="$(echo "$new_version" | sed 's/-dev/.01/')"
  sed -i 's/v\.[^<]*</v. '"$new_version ($release_name)"'</' "$DOC_HEADER"
  sed -i 's/\(project(.*VERSION \)".*")/\1"'"$new_cmake_version"'")/' "$CMAKELISTS"
  sed -i 's/\(.*OGDF_VERSION \)".*"/\1"'"$new_version"'"/' "$VERSION_H"
  ask_to_commit "Bump to version $new_version"
}

bump_version "$VERSION_NUMBER" "$RELEASE_NAME"
git tag
git tag $TAGNAME || die "git tag failed."
bump_version "$VERSION_NUMBER-dev" "$RELEASE_NAME"

### Create archive with OGDF source code
mkdir -p "$TARGET_DIR" || die "creating target directory (to save release in) failed"
TARGET_DIR="$(realpath "$TARGET_DIR")"
SRCZIP="$TARGET_DIR/ogdf.v${VERSION_NUMBER}.zip"
git archive --prefix=OGDF/ -9 -o "$SRCZIP" $TAGNAME || die "git archive failed"

### Set latest release branch to new release
git branch -f $LATEST_RELEASE_BRANCH $TAGNAME

cat <<EOF
Finished!

We created a tag with an updated version and tagged it with $TAGNAME.
The branch $LATEST_RELEASE_BRANCH was updated to $TAGNAME, and
an archive with the OGDF source code was created here:
  $SRCZIP
The version in the branch $RELEASE_BRANCH was further updated to a post-release development version.

To publish the release, do the following steps:

 1. Push to the Github repository and create a new pull request
     * git push -fu $GITHUB_URL $RELEASE_BRANCH
     * Create new pull reqest for $RELEASE_BRANCH: "Release $VERSION_NUMBER ($RELEASE_NAME)"

 WAIT UNTIL PIPELINE IS GREEN! If it is red, fix stuff and repeat the process.

 2. Add $RELEASE_NAME page in ogdf.net's Wordpress
     * Look at the last release page as a "template"
     * Upload file $SRCZIP
     * Copy the release notes (e.g. from the output of "grip $REL_PORTING_GUIDE")
     * Rewrite beginning as "Today we have released"...
     * Update links from the release notes, e.g., to the porting guide
     * Use category "release"
     * Publish
     * Use a permalink "$RELEASE_LCNAME" (probably done automatically after publish)
     * Update http://ogdf.net/team if necessary

 3. Release on Github
     * Merge the pull request (do not rebase it! otherwise the tag is lost and you have to recreate it)
     * git push $GITHUB_URL $TAGNAME
     * git push $GITHUB_URL $LATEST_RELEASE_BRANCH
     * Draft a new release on https://github.com/ogdf/ogdf (mimic the other releases)

 4. Additional things to consider
     * Update https://github.com/ogdf/ogdf-wheel
     * Update the internal research GitLab repo
     * Update vcpkg and conan packages
EOF
