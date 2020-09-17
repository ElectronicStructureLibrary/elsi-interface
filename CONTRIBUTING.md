# Contributing to ELSI

Contributions to the ELSI project are very welcome.

## The first step

[Open an issue](https://gitlab.com/elsi_project/elsi_interface/-/issues) to
briefly describe the proposed changes to the code. This is to avoid the
situation where two or more developers work on something similar without knowing
each other.

## Suggested workflow

1. Fork the project on `gitlab.com`
2. Clone the forked repository
3. Create a new branch (with a meaningful name) in it
4. Commit any changes to the newly created branch
5. Push the branch to the forked repository
6. Submit a merge request to the ELSI repository

## Merge requests

* Put `WIP:` in the merge request title to indicate work in progress.
* An opened merge request can be updated by pushing new commits to the branch or
  amending existing commits. The branch to be merged must
  1. Be up to date with the upstream master branch
  2. Not contain any revert commits
  3. Not contain any merge commits
  4. Pass the continuous integration tests
* Always rebase the branch onto the upstream master, instead of merging the
  upstream master into the branch. This helps maintain a linear, readable, and
  therefore more useful git history.
* It is highly recommended to clean up the commits in a merge request by an
  interactive rebase, i.e., `git rebase -i upstream/master`.
* Please make sure that each merge request introduces only one logical change to
  the code.

## Commit messages

* Write meaningful commit messages in this format:
  1. A single "title" line of 50 characters or less summarizing the commit. The
     title should be written in the imperative mood, e.g., "Fix typo in doc"
     instead of "This commit fixed typo in doc" or "typo fixed". This aligns
     with the choice made by git itself.
  2. An empty line.
  3. An optional "body" containing multiple lines, each line containing 72
     characters or less.
* When applicable, refer to existing issues and/or merge requests.

## Coding style

* In general, avoid tabs, trailing whitespaces, and consecutive blank lines.
* Maximal line length is 80 characters. Use continuation for longer lines.
* Nested blocks are indented by 3 white spaces.
* Variable names follow the `lower_case_with_underscore` convention. Do not use
  upper case except for constants, e.g., `MPI_COMM_WORLD`.
* Subroutines that ar part of the public API should be documented with the
  [Doxygen](https://www.doxygen.nl/index.html) tool.
