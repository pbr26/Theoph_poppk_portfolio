# Git & GitHub Pages — Setup and Per-Phase Workflow

Run these in the **Terminal** (or RStudio Terminal tab) from the project root.

## One-time initialisation

``` bash
cd "theoph_poppk_portfolio"
git init
git branch -M main
git add .gitignore README.md GIT_SETUP.md _quarto.yml *.Rproj .Rprofile index.qmd R/ reports/ renv.lock
git commit -m "chore(phase0): project setup — renv, quarto, theme, paths"

# Create an empty repo named theoph-poppk-portfolio on GitHub, then:
git remote add origin https://github.com/<your-username>/theoph-poppk-portfolio.git
git push -u origin main
```

## Per-phase workflow

After each phase is built and verified, snapshot and push:

``` bash
# (in R) renv::snapshot(prompt = FALSE)   # if dependencies changed
git add -A
git commit -m "feat(phaseN): <short summary>"
git push
```

Suggested commit messages:

| Phase | Message                                               |
|------:|-------------------------------------------------------|
|     1 | `feat(phase1): EDA script + report`                   |
|     2 | `feat(phase2): PKNCA noncompartmental analysis`       |
|     3 | `feat(phase3): base 1-cmt oral PopPK model`           |
|     4 | `feat(phase4): IIV on CL, V, Ka`                      |
|     5 | `feat(phase5): weight covariate model`                |
|     6 | `feat(phase6): GOF / VPC / CWRES diagnostics`         |
|     7 | `feat(phase7): clinical simulation`                   |
|     8 | `chore(phase8): publish site to GitHub Pages (docs/)` |

## GitHub Pages (Phase 8)

1.  Render the site: `quarto render` (writes HTML into `docs/`).
2.  Commit and push `docs/`.
3.  On GitHub: **Settings → Pages → Source = Deploy from a branch**, **Branch = `main`**, **Folder = `/docs`**, then **Save**.
4.  The site goes live at `https://<your-username>.github.io/theoph-poppk-portfolio/`.

> Note: `.gitignore` ignores `renv/library/` but **keeps `renv.lock`** — always commit the lockfile.
