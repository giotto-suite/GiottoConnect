#!/usr/bin/env bash
# Repo-specific preparation, run by the shared pkgdown workflow
# (giotto-suite/.github) from the package checkout pkgdown will build.
set -euo pipefail

# 1. Files that are not site content. The package keeps design notes under
# docs/, which is also pkgdown's output folder, and pkgdown publishes every
# Markdown file at the package root, so AGENTS.md would become a page.
rm -rf docs AGENTS.md
echo "pre-build: removed docs/ (design notes) and AGENTS.md"

# 2. Articles live with the package, in vignettes/articles/. The shared workflow
# replaces the checkout's vignettes/ with this branch's, which has none, so
# restore the package's own.
git checkout HEAD -- vignettes/
echo "pre-build: restored package articles ($(ls vignettes/articles/*.Rmd | wc -l | tr -d ' ') Rmd)"
