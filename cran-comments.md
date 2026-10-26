## Resubmission of PathoFMPred 0.1.0

Thank you for the review. We have expanded “PLS-LDA” in DESCRIPTION, added methodological references there in the requested DOI format, documented the class and meaning of the value returned by `predict_titan()`, and removed implicit writes to the user's home directory, package directory and working directory. Functions that create reports or download model collections now require an explicit output path. Examples, tests and vignette code use temporary directories.

The source package contains no fitted model collection. The optional Prov-GigaPath collection is downloaded only when a user explicitly calls `fetch_pathofmpred_models()` with a destination directory.

The package was built and checked on macOS Sonoma, arm64, R 4.6.0. `R CMD check --as-cran` completed with status OK, with remote CRAN incoming checks disabled because the CRAN archive metadata endpoint was temporarily unavailable from this machine. The local incoming checks, examples, tests, vignette rebuild and manual checks all passed.
