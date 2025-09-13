# Contributing Guide

Welcome! 
This repository is an academic project. We welcome contributions from students, researchers, and collaborators. Even if you are new to **GitHub** or **git**, this guide will walk you through the process from opening an issue to merging your code. It also covers our unit‑testing expectations using **GoogleTest (gtest)** integrated via **BLT/CMake** and run automatically in CI.

> **TL;DR**  
> - Always start with an **issue**.  
> - Create a **branch** named after the issue.  
> - Make **small, focused commits**.  
> - Open a **Pull Request (PR)** that references the issue (e.g., `Closes #123`).  
> - **All tests must pass** in CI. Add **new tests** for new features/bug fixes.  
> - Address review comments and then **merge** when approved.

---

## Before You Start
1. Make sure you have a [GitHub account](https://github.com/join).
2. Install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your computer.
3. (Optional) Install a Git GUI or editor integration like [GitHub Desktop](https://desktop.github.com/) or [VS Code](https://code.visualstudio.com/).

> If you are new to command‑line git, you can still follow this guide with GitHub Desktop or VS Code’s built‑in source control.

---

## Step 1: Create an Issue
1. Go to the **Issues** tab of the repository.
2. Click **New issue**.
3. Choose a template (if available), otherwise describe:
   - **What** you intend to change (bug, feature, docs).
   - **Why** it’s needed (motivation/use case).
   - **How** you plan to implement it (brief plan).
4. Submit the issue.

> Issues help us coordinate work and give you feedback **before** you write code.

---

## Step 2: Create a Branch
Create a branch on GitHub based on the default branch (`main` unless stated otherwise). Use this naming pattern:

```
issue-<ISSUE-NUMBER>-short-description
```

Examples:
```
issue-12-fix-readme-typo
issue-34-add-sparse-solver
```

---

## Step 3: Pull the Branch Locally
If you haven’t cloned the repo yet:
```bash
git clone https://github.com/<ORG_OR_USERNAME>/<REPO>.git
cd <REPO>
```

Fetch and check out your new branch:
```bash
git fetch origin
git checkout -b issue-12-fix-readme-typo origin/issue-12-fix-readme-typo
```

Confirm you’re on the branch:
```bash
git status
```

---

## Step 4: Make Local Changes
- Work in **small, logical chunks**.
- Keep your commits focused (one feature/fix per commit).
- Update documentation and comments as needed.

> Tip: If your change is large, split it into multiple PRs (e.g., “infrastructure”, then “feature”).

---

## Step 5: Commit and Push
Stage and commit your changes:
```bash
git add <files>
git commit -m "Short summary (Closes #<ISSUE-NUMBER>)"
```

Push your branch:
```bash
git push origin issue-12-fix-readme-typo
```

> Using `Closes #<ISSUE-NUMBER>` (or `Fixes #<ISSUE-NUMBER>`) in your commit or PR description will automatically close the issue when the PR is merged.

---

## Step 6: Open a Pull Request (PR)
1. On GitHub, you’ll see a banner for your pushed branch → click **Compare & pull request**.
2. Fill in the PR title and description:
   - **What** you changed and **why**.
   - Link the issue (e.g., `Closes #12`).
   - Mention any **follow‑ups** or limitations.
3. Submit the PR.

> Our CI will automatically build and run tests on your PR. Keep an eye on those checks.

---

## Step 7: Request and Address Code Review
- A maintainer (or classmate) will review your PR and leave comments.
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
git clone https://github.com/<ORG_OR_USERNAME>/<REPO>.git
cd <REPO>

# create/switch branches
git checkout -b issue-123-topic    # new branch from current
git checkout main                  # switch to main

# sync with remote main
git fetch origin
git checkout main
git pull origin main

# stage and commit
git add path/to/file.cpp
git commit -m "Implement feature X (Closes #123)"

# push branch
git push origin issue-123-topic
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
- For new/changed code, aim for **≥ 80% line coverage** (where coverage tools are enabled).  
- Prefer **unit tests** over integration tests when possible. Add integration tests for cross‑module behavior.

---

## Writing Tests (GoogleTest)

Place tests under:
```
tests/
  unit/
    test_<module>_<topic>.cpp
  integration/
    test_<feature>_integration.cpp
```

**Naming conventions**
- Test files: `test_<unit>.cpp`
- Test suites: `ClassOrModuleNameTest`
- Test cases: `DoesX`, `HandlesY`, `RejectsZ`

**Basic example**
```cpp
#include <gtest/gtest.h>
#include "mylib/vector_ops.hpp"

TEST(VectorOpsTest, DotProductBasic) {
    std::vector<double> a{1.0, 2.0, 3.0};
    std::vector<double> b{4.0, 5.0, 6.0};
    EXPECT_DOUBLE_EQ(dot(a,b), 32.0);
}
```

**Fixtures (recommended for setup/teardown)**
```cpp
class MatrixFixture : public ::testing::Test {
protected:
    void SetUp() override {
        // allocate or init common state
    }
    void TearDown() override {
        // cleanup
    }
};

TEST_F(MatrixFixture, InvertIdentity) {
    // EXPECT_* or ASSERT_* calls here
}
```

**Parameterized tests (for multiple inputs)**
```cpp
class NormParamTest : public ::testing::TestWithParam<std::tuple<double,double,double>> {};

TEST_P(NormParamTest, ComputesL2) {
    auto [x,y,expected] = GetParam();
    EXPECT_NEAR(l2_norm(x,y), expected, 1e-12);
}

INSTANTIATE_TEST_SUITE_P(
    BasicCases, NormParamTest,
    ::testing::Values(
        std::make_tuple(3.0, 4.0, 5.0),
        std::make_tuple(5.0,12.0,13.0)
    )
);
```

> Prefer `EXPECT_*` over `ASSERT_*` unless a failure should **abort** the test immediately.

---

## Adding Tests with BLT/CMake

If BLT is included via submodule or vendor directory, your `CMakeLists.txt` might look like this (simplified):

```cmake
# Top-level CMakeLists.txt excerpt
cmake_minimum_required(VERSION 3.18)
project(myproj LANGUAGES CXX)

# BLT setup (path may differ)
set(BLT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/blt")
include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

# Library / sources
add_library(mylib src/vector_ops.cpp)
target_include_directories(mylib PUBLIC include)

# GoogleTest via BLT
blt_add_executable(NAME test_vector_ops
                   SOURCES tests/unit/test_vector_ops.cpp
                   DEPENDS_ON mylib gtest)

# Register test with CTest via BLT
blt_add_test(NAME vector_ops_unit
             COMMAND test_vector_ops)
```

Notes:
- `DEPENDS_ON` ensures your test links against `mylib` and `gtest` (provided by BLT).
- You can create multiple test executables, one per file or logical group.
- Prefer **small** test executables for faster failure isolation.

---

## Building & Running Tests Locally

From the repo root:
```bash
# Configure (choose Debug for dev)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug

# Build everything (lib + tests)
cmake --build build -j

# Run all tests with CTest
ctest --test-dir build --output-on-failure
```

**Run a single test by name (regex):**
```bash
ctest --test-dir build -R vector_ops_unit --output-on-failure
```

**Filter at gtest level:**
```bash
./build/tests/unit/test_vector_ops --gtest_filter=VectorOpsTest.DotProductBasic
```

**Increase verbosity:**
```bash
CTEST_OUTPUT_ON_FAILURE=1 ctest --test-dir build -j
```

> If your platform needs special toolchains (e.g., Intel/Clang), consult the project’s `docs/BUILDING.md` if present.

---

## Test Data & Temporary Files

- Put small test data under `tests/data/`.  
- Use OS‑provided temp directories for outputs; clean up after tests.  
- Avoid writing into the source tree.  
- Keep test data small and under version control when feasible.

---

## Style, Linting, and Static Analysis (if enabled)

- Follow the project’s code style (see `.clang-format`, `.clang-tidy`, etc. if present).  
- CI may run format/lint checks. Fix issues locally and re‑push.  
- Keep public headers minimal and documented.

---

## Reproducibility & Determinism

- Seed random number generators explicitly (e.g., with a fixed seed) for deterministic tests.
- Do not rely on system time, filesystem order, or network calls.
- Keep floating‑point tolerances realistic; prefer `EXPECT_NEAR(x, y, 1e-12)` for doubles.

---

## PR Checklist (Pre‑Merge)

- [ ] Issue created and linked in PR (`Closes #…`).  
- [ ] Branch up‑to‑date with `main` and clean rebase if necessary.  
- [ ] Code builds locally (Debug/Release as relevant).  
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

---

## Workflow Diagram (High‑Level)

```
Issue → Branch → Commit → Push → Pull Request → Review/CI → Fixups → Approve → Merge → Close
```

Thanks for helping build this project!
