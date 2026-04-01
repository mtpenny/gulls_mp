# Release Instructions

### For Existing Users
**Nothing changes!** Your existing workflow, parameter files, and scripts work identically.

### For New Features
```bash
# Validate inputs before running
python scripts/validate_inputs.py your_file.prm

# Bump version and create release
python scripts/bump_version.py patch    # 2.0.0 -> 2.0.1
# Edit CHANGELOG.md with your changes
python scripts/bump_version.py release  # Commit, tag, and push automatically
```

### For Matt's Workflow
```bash
# 1. Make your changes
# 2. Bump version (creates blank changelog entry)
python scripts/bump_version.py patch

# 3. Edit CHANGELOG.md (describe your changes)
# 4. Create release (auto-generates RELEASE_NOTES.md from CHANGELOG.md)
python scripts/bump_version.py release
```

After merging this PR into the main repo (`dev` or `main` branch), you can make the first release simply using
```bash
python scripts/bump_version.py release
```

> Note. this creates a release on the version of the repo you are cloned from. If you run it on your personal fork, it will create the release on your fork, not the original repo.

The release notes have been created and are included in this PR. `conf.py` and `gulls.cpp` already have their versions marked with v2.0.0. Those two files are the only place the current version number exists in code and are handled by `bump_version.py`, should you want to change the current version. You can revert a version number using the `--revert` tag. For example:
```bash
# v2.0.0 -> v1.0.0
python scripts/bump_version.py major --revert
```
However, this will not automatically fix the versioning in the changelog and release notes.

**That's it!** The script handles:
- ✅ Commits changes (with smart unstaged change detection)
- ✅ Creates and pushes tags
- ✅ Triggers release workflow automatically
- ✅ Handles existing tags gracefully
- ✅ Version numbers in code

## Workflow Triggers

### CI Workflow (`.github/workflows/test.yml`)
- **Triggers**: Push to any branch, pull requests
- **What it does**: Builds Gulls, runs smoke tests, validates inputs
- **Manual trigger**: Go to Actions tab → "Test" → "Run workflow"

### Release Workflow (`.github/workflows/release.yml`)
- **Triggers**: Push tags matching `v*` (e.g., `v2.0.0`, `v2.0.1`)
  - Note. creation and pushing of these tags are handled by the command `python scripts/bump_version.py release`
- **What it does**: Builds, tests, creates GitHub release with source/binary archives
- **Release notes**: 
  - Auto-generates `RELEASE_NOTES.md` from `CHANGELOG.md` entries during `release` command
  - RELEASE_NOTES.md can be generated without preforming a release using `python scripts/bump_version.py release --dry-run`
  - Prompts before replacing existing, outdated `RELEASE_NOTES.md` with different version
  - If `RELEASE_NOTES.md` exists and has the correct version → Uses your content + appends smoke test plots
  - If no `RELEASE_NOTES.md` → Generates from changelog + smoke test plots
- **How to trigger**: 
  ```bash
  # Manual (old way)
  git tag v2.0.0
  git push origin v2.0.0
  
  # Automated (new way)
  python scripts/bump_version.py release
  ```

### Documentation Workflow (`.github/workflows/docs.yml`)
- **Triggers**: Push to `main` branch
- **What it does**: Builds Sphinx documentation
- **Note**: May be redundant with separate `gulls-microlensing.github.io` repo. 
- **Warning**: Workflow is untested, but very basic and uses standard workflow automations. It will probably work.

## Open Questions for Review

1. **Documentation hosting**: Currently builds docs in this repo, but `gulls-microlensing.github.io` exists separately. Should we:
   - Make `gulls-microlensing.github.io` a submodule of this repo?
   - Remove the docs workflow from this repo?
   - Keep both (redundant but safe)?

2. **Version strategy**: The version bumping script creates Git tags. Do you want to:
   - Use it for official releases?
   - Keep manual/no version management?

## Technical Changes (Changes to the source code)
- **Buffer size fixes** - Required for CI environment (long paths)
- **Version updates** - v2.0.0 with proper date (October 2025)
- **Stub implementations** - GSL fallbacks for CI (Numerical Recipes preferred for production)
- **Release warnings** - Clear notices about GSL fallbacks in binary releases

## Files Changed
- Added: CI workflows, documentation, validation scripts, version management
- Removed: Build artifacts (hundreds of files - makes diff look larger than it is)
- Modified: Buffer sizes, version numbers, added stubs

## [Example Automated Release](https://github.com/AmberLee2427/gulls_mp/releases)

---

| ![https://github.com/user-attachments/assets/1d336d37-b768-4d1d-bbbf-ffaf42a97128](https://github.com/user-attachments/assets/1d336d37-b768-4d1d-bbbf-ffaf42a97128) |
| :-: |

---