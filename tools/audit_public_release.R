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
relative_rds <- tracked[grepl("[.]rds$", tracked) &
                          file.exists(file.path(root, tracked))]
if (length(relative_rds)) {
  stop(
    "Public source tree contains a fitted RDS artifact: ",
    paste(relative_rds, collapse = ", ")
  )
}

# Git tracking alone is insufficient: an ignored or untracked object anywhere
# under inst/ would still be bundled by R CMD build. This includes prediction
# reference distributions, not only fitted coefficients.
physical_rds <- list.files(
  file.path(root, "inst"), pattern = "[.]rds$",
  recursive = TRUE, full.names = TRUE
)
physical_relative <- substring(physical_rds, nchar(root) + 2L)
if (length(physical_relative)) {
  stop("Public inst/ contains a fitted RDS object: ",
       paste(physical_relative, collapse = ", "))
}
if (file.exists(file.path(root, "inst", "extdata", "model_registry.csv"))) {
  stop("The research atlas registry must not be bundled")
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
if (grepl("releases/download/models-v", store, fixed = TRUE) ||
    grepl("PathoFMPred.model_base_url", store, fixed = TRUE)) {
  stop("Public package must not contain a default research-asset source")
}
for (text in c("expected_sha256", "missing(base_url)",
               "missing(expected_sha256)")) {
  if (!grepl(text, store, fixed = TRUE)) {
    stop("Public object-store safeguard is missing: ", text)
  }
}

message("Public release boundary audit passed")
