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

# Git tracking alone is insufficient: an ignored or untracked object anywhere
# under inst/ would still be bundled by R CMD build. This includes prediction
# reference distributions, not only fitted coefficients.
physical_rds <- list.files(
  file.path(root, "inst"), pattern = "[.]rds$",
  recursive = TRUE, full.names = TRUE
)
physical_relative <- substring(physical_rds, nchar(root) + 2L)
if (!setequal(physical_relative, allowed_embedded)) {
  stop("Public inst/ contains an unexpected RDS object: ",
       paste(setdiff(physical_relative, allowed_embedded), collapse = ", "))
}

example_path <- file.path(root, allowed_embedded)
if (!file.exists(example_path)) {
  stop("The minimal Giga-SSL example collection is missing")
}
example <- readRDS(example_path)
if (!identical(example$foundation_model, "GigaSSL") ||
    length(example$feature_names) != 512L ||
    length(example$models) != 2L ||
    nrow(example$registry) != 2L ||
    any(example$registry$foundation_model != "GigaSSL") ||
    any(!startsWith(names(example$models), "GigaSSL__")) ||
    any(vapply(example$models, function(model) {
      !identical(model$foundation_model, "GigaSSL")
    }, logical(1)))) {
  stop("The embedded public example does not contain exactly two Giga-SSL models")
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
  "GigaSSL = \"874b84cb040d7c955b96f95b29ee181c9056f5e40a3feda28171f942063d0745\"",
  "ProvGigaPath = \"72868d9fb7562f3eafaf7cc287cb4e909e1767587662be07625c985a028da832\""
)
for (text in required_store_text) {
  if (!grepl(text, store, fixed = TRUE)) {
    stop("Public object-store safeguard is missing: ", text)
  }
}

message("Public release boundary audit passed")
