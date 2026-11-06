# Analysis environment

The analyses were run on Windows with **R 4.4.1**, Arial and Cairo support. `package_versions.tsv` records package versions; `sessionInfo.txt` and `runtime.tsv` record the runtime. The root `renv.lock` fixes the dependencies.

If the matching packages are installed, run the individual module scripts directly. To restore dependencies to a separate local library, run from the repository root:

```sh
Rscript environment/restore_environment.R --library r-library
```

Then make that library available to subsequent `Rscript` commands. In the same PowerShell terminal:

```powershell
$env:R_LIBS = (Resolve-Path "r-library").Path
```

`restore_environment.R`, `project_library.R` and `bootstrap_package.tsv` support this optional restoration. Downloads require access to CRAN/Bioconductor archives; compiling packages on Windows may require Rtools44. Full analysis execution was verified with matching installed packages; a separate clean online installation of every dependency and other operating systems have not been fully tested.
