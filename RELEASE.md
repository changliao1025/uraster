# Release Workflow

This document describes the process for releasing new versions of uraster.

## Prerequisites

1. Ensure all tests pass locally and in CI
2. Update version numbers in relevant files
3. Update CHANGELOG.md with new version information
4. Ensure documentation is up to date

## Version Numbering

uraster follows [Semantic Versioning](https://semver.org/):
- **MAJOR** version for incompatible API changes
- **MINOR** version for backwards-compatible functionality additions
- **PATCH** version for backwards-compatible bug fixes

## Release Process

### 1. Prepare the Release

```bash
# Create a release branch
git checkout -b release/v0.2.0

# Update version in pyproject.toml
# Update version in uraster/__init__.py
# Update version in CITATION.cff
# Update CHANGELOG.md

# Commit changes
git add .
git commit -m "Prepare release v0.2.0"
```

### 2. Update Version Numbers

Update the version number in these files:
- `pyproject.toml` (version field)
- `uraster/__init__.py` (__version__ variable)
- `CITATION.cff` (version field)
- `docs/conf.py` (release variable)

### 3. Update CHANGELOG.md

Add a new section for the release with:
- Release date
- Summary of changes
- Breaking changes (if any)
- New features
- Bug fixes
- Dependencies updates

### 4. Test the Release

```bash
# Install in development mode
conda install -c conda-forge gdal geovista vtk=9.3.0 pyearth
conda develop .

# Run all tests
python -m pytest tests/ -v

# Test building the package
python -m build

# Test the built package
pip install dist/uraster-*.whl
python -c "import uraster; print(uraster.__version__)"
```

### 5. Create Pull Request

```bash
# Push the release branch
git push origin release/v0.2.0

# Create a pull request to main branch
# Ensure all CI checks pass
# Get approval from maintainers
```

### 6. Create GitHub Release

1. Go to GitHub repository
2. Click "Releases" → "Create a new release"
3. Create a new tag: `v0.2.0`
4. Release title: `uraster v0.2.0`
5. Copy relevant section from CHANGELOG.md to release description
6. Check "Set as the latest release" if appropriate
7. Click "Publish release"

### 7. Automated Deployment

The GitHub Actions workflow will automatically:
1. Build the package
2. Run tests
3. Deploy to PyPI (if release is published)
4. Build conda package (manual upload to conda-forge required)

### 8. Post-Release Tasks

1. **Conda-forge Update**: Create a PR to the [conda-forge feedstock](https://github.com/conda-forge/uraster-feedstock) with the new version
2. **Documentation**: Ensure documentation is built and deployed
3. **Announcement**: Announce the release on relevant channels
4. **Monitor**: Watch for any issues reported by users

## Conda-forge Release

After PyPI release, update the conda-forge recipe:

1. Fork the [uraster-feedstock](https://github.com/conda-forge/uraster-feedstock)
2. Update `recipe/meta.yaml`:
   - Update version number
   - Update source URL and SHA256 hash
   - Update dependencies if needed
3. Create a pull request
4. Wait for CI to pass and maintainer approval

## Hotfix Releases

For critical bug fixes:

```bash
# Create hotfix branch from main
git checkout main
git pull origin main
git checkout -b hotfix/v0.1.1

# Make minimal changes to fix the issue
# Update version to patch level (e.g., 0.1.0 → 0.1.1)
# Update CHANGELOG.md

# Follow normal release process
```

## Rollback Procedure

If a release has critical issues:

1. **Immediate**: Remove the problematic release from PyPI if possible
2. **Communication**: Notify users via GitHub issues and documentation
3. **Fix**: Create a hotfix release with the fix
4. **Post-mortem**: Document what went wrong and how to prevent it

## Release Checklist

- [ ] All tests pass locally and in CI
- [ ] Version numbers updated in all files
- [ ] CHANGELOG.md updated
- [ ] Documentation updated
- [ ] Release branch created and tested
- [ ] Pull request created and approved
- [ ] GitHub release created
- [ ] PyPI deployment successful
- [ ] Conda-forge PR created
- [ ] Release announced
- [ ] Post-release monitoring

## Contact

For questions about the release process, contact the maintainers:
- Chang Liao (changliao.climate@gmail.com)