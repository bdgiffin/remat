# Contributing Guide

Welcome! 
If you are new to **GitHub** or **git**, this guide will walk you through the process from opening an issue to merging your code. It also covers our unit‑testing expectations using **GoogleTest (gtest)**.

> **General Procedure**  
> - Crate an **issue** in GitHub.  
> - Create a **branch** for the issue.  
> - Make **small, focused commits**.  
> - Open a **Pull Request (PR)** that references the issue (e.g., `Closes #123`).  
> - **All tests must pass** in CI. Add **new tests** for new features/bug fixes.  
> - Address review comments and then **merge** when approved.

---

## Before You Start
1. Make sure you have a [GitHub account](https://github.com/join).
2. Install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your computer.

---

## Step 1: Create an Issue
Create an issue [from within GitHub](https://docs.github.com/en/issues/tracking-your-work-with-issues/using-issues/creating-an-issue).

---

## Step 2: Create a Branch
Create a branch to work on an issue following the [suggested procedure in GitHub](https://docs.github.com/en/issues/tracking-your-work-with-issues/using-issues/creating-a-branch-for-an-issue).

---

## Step 3: Pull the Branch Locally
If you haven’t cloned the repo yet:
```bash
git clone git@github.com:bdgiffin/remat.git
cd remat
```

Fetch and check out your new branch:
```bash
git fetch origin
git checkout -b <branch-name> origin/<branch-name>
```

Confirm you’re on the branch:
```bash
git status
```

---

## Step 4: Make Local Changes
- Work in **small, logical chunks**.
- Try to keep your commits focused (one feature/fix per commit).
- Update documentation and comments as needed.

---

## Step 5: Commit and Push
Stage and commit your changes:
```bash
git add <files>
git commit -m "Short one-line summary of changes."
```

Push your branch:
```bash
git push origin <branch-name>
```

> Using `Closes #<ISSUE-NUMBER>` (or `Fixes #<ISSUE-NUMBER>`) in your commit message or PR description will automatically close the issue when the PR is merged.

---

## Step 6: Open a Pull Request (PR)
1. On GitHub, you’ll see a banner for your pushed branch → click **Compare & pull request**.
2. Fill in the PR title and description:
   - **What** you changed and **why**.
   - Link the issue (e.g., `Closes #12`).
   - Mention any **follow‑ups** or limitations.
3. Submit the PR.

> A CI pipeline has been set up to automatically build and run tests on your PR. Keep an eye on those checks.

---

## Step 7: Request and Address Code Review
- A maintainer will review your PR and leave comments.
- To address feedback:
  1. Make changes locally.
  2. Commit and **push to the same branch**.
  3. The PR updates automatically.
- Mark conversations as **Resolved** when addressed.

---

## Step 8: Merge and Close
- When CI is green and the PR is **approved**, click **Merge** (we typically use **Squash and merge** to keep history tidy, unless otherwise stated).
- The linked issue should close automatically if you used `Closes #<ISSUE>` in the PR description or commit message.
- Done! Thank you for contributing.

---

## Quick Reference: Common Commands
```bash
# clone repository
git clone git@github.com:bdgiffin/remat.git
cd remat

# create/switch branches
git checkout -b <branch-name>    # new branch from current
git checkout main                # switch to main branch

# sync with remote main
git fetch origin
git checkout main
git pull origin main

# stage and commit
git add path/to/file.cpp
git commit -m "Implement feature X (Closes #123)"

# push branch
git push origin <branch-name>
```

---

# Testing Policy (gtest via BLT/CMake & CI)

We use **GoogleTest (gtest)** integrated through **BLT** and **CMake**. All tests run automatically in **CI** on every PR.

**Hard rules:**  
- **All existing tests must pass** before a PR can be merged.  
- **New features must include new unit tests** demonstrating correctness.  
- **Bug fixes must include a regression test** that fails before the fix and passes after.  
- Keep tests **fast, deterministic, and isolated** (no network, no global state, no filesystem side‑effects unless temporary).  

**Targets:**  
- Prefer **unit tests** over integration tests when possible. Add integration tests for cross‑module behavior.

---

## Writing Tests (GoogleTest)

Place any unit tests under the `tests` subdirectory. Refer to the existing tests contained in this directory as examples of how to structure new tests.

> Prefer `EXPECT_*` over `ASSERT_*` unless a failure should **abort** the test immediately.

---

## Adding Tests with BLT/CMake

Make sure that the `tests/CMakeLists.txt` file is modified to include any new tests. Follow a similar pattern to how the other tests are declared using the `blt_add_executable` and `blt_add_test` macros.

---

## Building & Running Tests Locally

From within the locally created CMake `build` directory, you can execute the unit test suite using:
```bash
make test
```

---

## Reproducibility & Determinism

- Seed random number generators explicitly (e.g., with a fixed seed) for deterministic tests.
- Keep floating‑point tolerances realistic; prefer `EXPECT_NEAR(x, y, 1e-12)` for doubles.

---

## PR Checklist (Pre‑Merge)

- [ ] Issue created and linked in PR (`Closes #...`).
- [ ] Branch up‑to‑date with `main` and clean rebase if necessary.
- [ ] Code builds locally.
- [ ] **All tests pass locally**.
- [ ] **New/changed code has tests** (unit and/or integration).
- [ ] Docs and examples updated if needed.
- [ ] CI checks are **green**.

> **Reminder:** PRs with failing CI or missing tests will not be merged.

---

## Code Review Guidelines (for reviewers)

- Be kind and constructive—this is a learning‑oriented project.
- Prefer **specific, actionable** feedback (what & why).
- Request tests where coverage is missing.
- Check for clear names, comments, and small, focused commits.
- Approve when the code is correct, tested, and maintainable.
