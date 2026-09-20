#!/usr/bin/env Rscript

# Fail closed when the public source tree crosses a licensing or fitted-object
# boundary. This audit is deliberately independent of package loading so it
# runs before R CMD check and can inspect files that are not installed.

root <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(file.exists(file.path(root, "DESCRIPTION")))

description <- read.dcf(file.path(root, "DESCRIPTION"))
if (!identical(unname(description[1L, "License"]), "MIT + file LICENSE")) {
  stop("DESCRIPTION must declare MIT + file LICENSE")
}
for (path in c("LICENSE", "LICENSE.md", "inst/licenses/MODEL_ACCESS.md")) {
  if (!file.exists(file.path(root, path))) stop("Missing licensing file: ", path)
}

tracked <- system2("git", c("-C", root, "ls-files"), stdout = TRUE)
if (!length(tracked)) stop("Could not enumerate tracked public-release files")
relative_rds <- tracked[grepl("[.]rds$", tracked)]
allowed_embedded <- file.path(
  "inst", "models", "pathofmpred_gigassl_example.rds"
)
unexpected <- setdiff(relative_rds, allowed_embedded)
if (length(unexpected)) {
  stop(
    "Public source tree contains an unapproved fitted RDS artifact: ",
    paste(unexpected, collapse = ", ")
  )
}

if (any(grepl("titan", basename(relative_rds), ignore.case = TRUE))) {
  stop("Public source tree must not contain a TITAN fitted object")
}

builder <- paste(readLines(file.path(root, "tools", "build_package_data.R"),
                           warn = FALSE), collapse = "\n")
forbidden_builder_patterns <- c(
  "file[.]copy\\(model_source",
  "file[.]copy\\(addition_source",
  "inst.*models.*TITAN.*package_file"
)
for (pattern in forbidden_builder_patterns) {
  if (grepl(pattern, builder, perl = TRUE)) {
    stop("Public data builder contains a fitted-object copy path: ", pattern)
  }
}
if (!grepl("Public package tree contains fitted", builder, fixed = TRUE)) {
  stop("Public data builder does not fail closed on embedded fitted objects")
}

store <- paste(readLines(file.path(root, "R", "object_store.R"), warn = FALSE),
               collapse = "\n")
required_store_text <- c(
  "TITAN-derived fitted parameters are deliberately unavailable",
  "GigaSSL = \"3ccfb101a24cc1ac9c6564d42c3392a0a139ccddd052cefb326db47d1d54ae1c\"",
  "ProvGigaPath = \"712b39f122f5efddc8d182700d84f47314943870de42f5596aa4bec6daeb67c8\""
)
for (text in required_store_text) {
  if (!grepl(text, store, fixed = TRUE)) {
    stop("Public object-store safeguard is missing: ", text)
  }
}

message("Public release boundary audit passed")
