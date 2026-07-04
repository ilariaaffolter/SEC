# Setup & Git/RStudio workflow — beginner guide

This guide assumes you are **not** a programmer. It walks through everything
from a blank computer to running the analysis and syncing with GitHub. Copy-paste
the commands exactly. Lines starting with `#` are comments you don't need to type.

---

## Part 1 — Install the tools (once per computer)

1. **R** — https://cran.r-project.org/ → download for your OS → install.
2. **RStudio Desktop** (free) — https://posit.co/download/rstudio-desktop/ → install.
3. **Git** — https://git-scm.com/downloads → install with default options.
   - *Windows only*: also install **Rtools** (matching your R version) from
     https://cran.r-project.org/bin/windows/Rtools/ — needed to build some packages.
   - *macOS only*: open the Terminal app and run `xcode-select --install`.

To check Git is visible to RStudio: open RStudio → *Tools ▸ Global Options ▸ Git/SVN*
→ the "Git executable" box should point to your git. If empty, browse to it
(`C:/Program Files/Git/bin/git.exe` on Windows, `/usr/bin/git` on macOS/Linux).

---

## Part 2 — Connect RStudio to GitHub (once)

You already have a GitHub repository for this project. You need to (a) tell Git
who you are and (b) let RStudio authenticate to GitHub.

### 2.1 Tell Git who you are
In RStudio, open the **Console** (bottom-left) and run, with your details:
```r
install.packages("usethis")      # a helper package (only needed once)
usethis::use_git_config(
  user.name  = "Your Name",
  user.email = "ilariaaffolter@gmail.com"   # the email on your GitHub account
)
```

### 2.2 Create a GitHub access token (your "password" for Git)
GitHub no longer accepts your normal password from tools. You create a token once:
```r
usethis::create_github_token()
```
This opens GitHub in your browser. Keep the defaults, scroll down, click
**Generate token**, and **copy** the long string (starts with `ghp_...`).
Then store it on your computer:
```r
install.packages("gitcreds")
gitcreds::gitcreds_set()
# paste the ghp_... token when asked
```
That's it — RStudio can now pull/push without asking for a password.

> Prefer clicking over typing? Install **GitHub Desktop**
> (https://desktop.github.com/), sign in, and use it instead of steps 2.2–4 for
> commit/pull/push. RStudio and GitHub Desktop can be used on the same folder.

---

## Part 3 — Get the project onto your computer (once)

In RStudio: **File ▸ New Project ▸ Version Control ▸ Git**, then:
- **Repository URL**: paste the repo's URL (the green *Code* button on GitHub →
  HTTPS → copy). It looks like `https://github.com/<user>/<repo>.git`.
- **Create project as subdirectory of**: pick a folder (e.g. your Documents).
- Click **Create Project**.

RStudio downloads the project and opens it. From now on, **always open the
project by double-clicking `SEC.Rproj`** (or *File ▸ Open Project*). This is what
makes `here()` and `renv` work correctly.

---

## Part 4 — Build the R environment (once per computer)

With the project open, in the **Console**:
```r
source("setup.R")
```
This installs every package the analysis needs and writes `renv.lock`. It can
take 20–40 minutes the first time. When it finishes, commit the lockfile (Part 6).

If someone hands **you** a project that already has a `renv.lock`, you don't run
`setup.R` — instead run:
```r
install.packages("renv")   # if you don't have it
renv::restore()            # installs the exact versions recorded in renv.lock
```

---

## Part 5 — Add your data and run the analysis

1. Put your input files in **`data/raw/`** (see `data/raw/README.md` for the
   expected file names).
2. Open `analysis/DiffAnalysis_yeast_QTL.Rmd`.
3. Click **Knit** (to produce the HTML report), or run it chunk by chunk with the
   green ▶ buttons.

Generated figures/tables/RData go into **`output/`** (these are not committed).

---

## Part 6 — The daily Git loop (pull → work → commit → push)

Use the **Git** tab (top-right pane in RStudio). The rhythm is:

1. **Pull first** (download others' latest changes): Git tab → **Pull** (the down
   arrow ⬇). Do this every time before you start working.
2. **Do your work** (edit the `.Rmd`, etc.). Save your files.
3. **Commit** (record a checkpoint):
   - Git tab → tick the boxes next to the files you changed (the "Staged" column).
   - Click **Commit**.
   - Write a short message describing what you did (e.g. "Add MW calibration step").
   - Click **Commit**.
4. **Push** (upload to GitHub): click **Push** (the up arrow ⬆).

That's the whole loop. Commit often (small, described steps); push when you want
your work saved on GitHub / shared.

> **What to commit:** your code (`.Rmd`, `R/`, `docs/`), `renv.lock`, and config
> files. **Not** committed automatically: the package library (`renv/library/`),
> generated `output/`, and raw data — that's intentional (see `.gitignore`).

### Added or updated a package?
If you `install.packages("something")` for the analysis, record it so others get
it too:
```r
renv::snapshot()   # updates renv.lock
```
then commit `renv.lock`.

---

## Part 7 — Handing the project over

The person receiving it should:
1. Install R, RStudio, Git (Part 1).
2. Clone the repo (Part 3) — or just open the folder you send them in RStudio via
   `SEC.Rproj`.
3. Run `renv::restore()` (Part 4) to get the exact package versions.
4. Put the data in `data/raw/` (you may need to send the data separately — it's
   not in the repo by default).
5. Knit `analysis/DiffAnalysis_yeast_QTL.Rmd`.

---

## Troubleshooting

- **A package fails to compile (e.g. `Rmpfr`, `igraph`)** — it needs a *system*
  library:
  - Ubuntu/Debian Linux: `sudo apt-get install libgmp-dev libmpfr-dev libxml2-dev libcurl4-openssl-dev libssl-dev`
  - macOS (with Homebrew): `brew install gmp mpfr`
  - Windows: make sure **Rtools** is installed (Part 1).
- **`here()` points to the wrong place** — you opened a file instead of the
  project. Close and reopen via `SEC.Rproj`.
- **Push rejected / "updates were rejected"** — someone pushed before you. Click
  **Pull**, resolve anything flagged, then **Push** again.
- **Line-ending warnings (LF/CRLF)** on Windows — harmless; you can silence them
  with `git config --global core.autocrlf true` in the RStudio *Terminal* tab.
- **Authentication failed when pushing** — your token expired or wasn't saved.
  Redo step 2.2 (`gitcreds::gitcreds_set()`).
