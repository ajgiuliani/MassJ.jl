# Contributing to MassJ.jl

Thanks for considering a contribution — bug reports, fixes, new features,
and documentation improvements are all welcome.

## Reporting bugs

Open an issue with the
[bug report template](.github/ISSUE_TEMPLATE/bug_report.md). The most useful
report includes:

1. A minimal Julia snippet that reproduces the problem.
2. The exact error message + stacktrace.
3. Your `MassJ` version (`pkg> status MassJ`) and Julia version
   (`julia --version`).
4. If the bug involves a specific file, either attach a small public-domain
   example or paste the first few lines so we can see the structure.

## Proposing a feature

Open an issue with the
[feature request template](.github/ISSUE_TEMPLATE/feature_request.md) and
sketch:

- The use case (what are you trying to compute?).
- What you'd expect the API to look like (one or two snippets).
- Whether a related feature already exists in MassJ that we should extend
  instead of duplicating.

For larger changes, open the issue *before* writing code so the design can
be discussed.

## Submitting a pull request

MassJ follows GitHub Flow with short-lived feature branches off `master`:

1. Fork the repository.
2. Create a branch off `master`: `git checkout -b feature/short-name`.
3. Make focused commits — one logical change per commit. The body of the
   commit message should explain *why* the change is needed.
4. Add or update tests in `test/runtests.jl`. We aim for full coverage of
   public API and any non-trivial branch.
5. Run the test suite locally:
   ```
   julia> ]
   pkg> activate .
   pkg> test
   ```
   All tests must pass on Julia 1.10 and 1.12 (CI verifies this).
6. Update the documentation:
   - manual page in `docs/src/man/`,
   - API entries in `docs/src/reference.md`,
   - `CHANGELOG.md` under an "Unreleased" section.
7. Open the pull request against `master` with the
   [PR template](.github/PULL_REQUEST_TEMPLATE.md).

## Coding conventions

- **Public API**: anything exported from `src/MassJ.jl` is considered public
  and follows SemVer. Breaking changes are reserved for major-version bumps.
- **Naming**: snake_case for functions, CamelCase for types, UPPER_CASE for
  module-level constants.
- **Docstrings**: every exported function and type needs one. Include at
  least one runnable example for top-level functions.
- **Tests**: new code goes with new tests. If you're fixing a bug, add a
  regression test that fails before your fix.
- **No dependencies without discussion**: a new entry in `Project.toml` is a
  cost (load time, conflict surface). Justify it in the PR description.

## Architecture notes

The package is structured around multiple dispatch:

- `MScontainer` — abstract supertype for spectrum data (`MSscan`,
  `MSscans`, `Chromatogram`, `MSrun`, …).
- `MethodType` — algorithm selectors (`SG`, `TBPD`, `TopHat`, …) passed via
  the `method` keyword to processing functions.
- `FilterType` — data selectors (`Level`, `RT`, `Precursor`, …) accepted by
  `extract`/`average`/`chromatogram` and composed in a single pass via the
  predicate machinery in `src/predicates.jl`.

Each file format has its own reader (`src/mzml.jl`, `src/mzxml.jl`, …);
exports live in `src/export.jl`. The yields pipeline is in `src/yields.jl`.

## Release process (for maintainers)

1. Update `CHANGELOG.md` with the release section.
2. Bump `version` in `Project.toml` and `CITATION.cff`.
3. Tag the merge commit: `git tag -a vX.Y.Z -m "..."` then `git push --tags`.
4. Use `JuliaRegistrator` (or `@JuliaRegistrator register()` as a PR
   comment) to register the new release in the General registry.

## Code of conduct

By participating, you agree to abide by the
[Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
Report unacceptable behaviour to the package maintainer at the email in
`CITATION.cff`.
