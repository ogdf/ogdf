[OGDF](../README.md) » [Developer's Guide](dev-guide.md) » git, GitHub, and Workflow

# git, GitHub, and Workflow

Find a brief overview of our core toolchain below.

## git

We use git for version control.
Use git to receive a local copy of this repository.

```
$ git clone git@github.com:ogdf/ogdf.git
```

We recommend to use SSH keys to authenticate yourself to the git server.
Do not use plain `http`.
To use SSH key authentication you need to create a pair of keys and authorize your GitHub account with the public key.

We refer to the
[GitHub Documentation](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
for an example of creating an SSH key.
You must add your public key [here](https://github.com/settings/keys).


Use `git status` to receive a synopsis of the current state.
```
$ git status
On branch peter-issue-39-readme
Your branch is up-to-date with 'origin/peter-issue-39-readme'.
Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

        modified:   main.cpp

Untracked files:
  (use "git add <file>..." to include in what will be committed)

        test.cpp
        Makefile

no changes added to commit (use "git add" and/or "git commit -a")
```

### Resources

* our [workflow](https://guides.github.com/introduction/flow/) explained in 5 minutes
* simple [git Tutorial](https://www.atlassian.com/git/tutorials/)
* [Pro Git](http://git-scm.com/book/en/v1) comprehensive book about git

## GitHub and Workflow

The OGDF master branch reflects the latest stable changes and can​**not** be directly modified.
Every feature (or bugfix) will be developed on a separate git branch (so called feature branches).
If a feature is completed, it will be integrated into the master branch, using a pull request.
We will occasionally review these changes and the continuous integration service will test every pull requests for errors.

### Example

```
$ git clone git@github.com:ogdf/ogdf.git
$ cd ogdf
$ git switch -c peter-awesome-feature
...make some changes...
$ git add -p
$ git commit -m "Introduce awesome feature"
$ git push -u origin peter-awesome-feature
```

The above example clones the entire repository and
creates a new local branch using the `git switch -c` command.
Then, some changes are made to the files which are added to the next commit using `git add -p`.
The commit itself is created by `git commit`.
If you use git for the first time, make sure you have author and e-mail address set before committing.

Finally, the changes are pushed back to the remote repository using `git push`.
Since the branch `peter-awesome-feature` does not exists in the remote repository, we have to
specify the branch to be created using the `-u` flag.

After performing the above steps you will be able to see the branch
`peter-awesome-feature` in our
[GitHub web interface](https://github.com/ogdf/ogdf/pulls).
You could now create a pull requests for this branch.

Make sure to be in a clean state when branching (i.e., creating a new feature branch).
Usually you want to start at the latest commit of the master branch:
```
$ git fetch origin master
$ git switch -c susi-feature origin/master
```

### Conventions

Our [Coding Standards](standards.md#version-control) provide information
about commit granularity, commit messages and naming convention for branches.

We use [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to ensure
that the code adheres to a certain code style.
For this style to be consistent, developers need to format their changes with
the correct version of (git-)clang-format (the currently used version is given
by the header of the `.clang-format` file).
The script `util/test_clang_format.sh` can help to detect formatting errors and
fix them (when passed the `-f` option). If the correct clang-format version is
not installed locally, the script can still correctly format the code by using a
docker image.

### Pull Requests

Every pull request constitutes a proposed solution for an issue. A pull request usually targets the master branch
(i.e., proposes changes to the master branch) from a feature branch.
A GitHub user may issue a pull request using the web interface.
Only core developers can accept pull request.
Usually they will review your changes and make comments using the GitHub interface.
Pull requests will always be rejected if the continuous integration service detects any errors.

If you want your changes to be reviewed but not merged, mark your pull request as work in progress
by prefixing your pull request title by  `WIP:` or `[WIP]`.

### Git helpers

There are some helper scripts in the `util/git` directory.
You can add this directory to your PATH of executables to have them handy.

#### git-ogdf-move-to-new-branch

Use `git-ogdf-move-to-new-branch <new-feature-branch>` if you accidentally
created some commits on master and want to move your commits to a new feature branch.

#### git-ogdf-push-new-branch

Use `git-ogdf-push-new-branch` to push your newly generated feature branch.

#### git-ogdf-fix-commit

Use `git-ogdf-fix-commit <commit>` if you want to change the commit `<commit>`
in your feature branch, for example, after a review. The helper will checkout
`<commit>` such that you can apply changes directly. After that use `git add ...`
or `git add -p` to stage your changes. After making sure (`git status`) that there
are no unstaged changes, you can finish your session by invoking
`git-ogdf-fix-commit` (without any argument).

Note that in some cases of using `git-ogdf-fix-commit`, applying the subsequent
commits could result in conflicts. In that case you need to fix them, stage these
fixes (`git add ...`) and run `git rebase --continue`.

Also note that you may need to force-push your branch using `git push -f`
because the script may change the history of the branch.

#### git-ogdf-remove-merged-branches

Use `git-ogdf-remove-merged-branches` in order to remove all merged local
branches and prune remote branches that have been removed.

### Issue Tracker

We are using the [issue tracker](https://github.com/ogdf/ogdf/issues) integrated into GitHub to keep track of tasks and bugs.
