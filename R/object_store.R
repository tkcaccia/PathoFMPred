.pathofm_cache <- new.env(parent = emptyenv())

.pathofm_slug <- function(foundation_model) {
  switch(
    foundation_model,
    TITAN = "titan",
    GigaSSL = "gigassl",
    ProvGigaPath = "provgigapath",
    tolower(gsub("[^A-Za-z0-9]+", "_", foundation_model))
  )
}

#' Validate a PathoFMPred model collection
#'
#' @param object An object created by [create_pathofmpred_object()] or a
#'   distributed PathoFMPred model collection.
#' @return `TRUE`, invisibly, when the object is valid.
#' @export
validate_pathofmpred_object <- function(object) {
  required <- c("format_version", "foundation_model", "feature_names",
                "registry", "models", "provenance")
  missing <- setdiff(required, names(object))
  if (length(missing)) {
    stop("Invalid PathoFMPred object; missing: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.list(object$models) || is.null(names(object$models)) ||
      any(!nzchar(names(object$models))) || anyDuplicated(names(object$models))) {
    stop("object$models must be a uniquely named list.", call. = FALSE)
  }
  if (!is.data.frame(object$registry) || !"model_id" %in% names(object$registry)) {
    stop("object$registry must contain model_id.", call. = FALSE)
  }
  absent <- setdiff(object$registry$model_id, names(object$models))
  if (length(absent)) {
    stop("Registry entries have no fitted object: ",
         paste(utils::head(absent, 5L), collapse = ", "),
         if (length(absent) > 5L) " ..." else "", call. = FALSE)
  }
  invisible(TRUE)
}

.pathofm_user_model_dir <- function() {
  override <- getOption("PathoFMPred.model_dir", NULL)
  if (!is.null(override)) return(normalizePath(override, mustWork = FALSE))
  tools::R_user_dir("PathoFMPred", which = "data")
}

.pathofm_object_candidates <- function(foundation_model) {
  slug <- .pathofm_slug(foundation_model)
  filename <- paste0("pathofmpred_", slug, ".rds")
  c(
    system.file("models", filename, package = "PathoFMPred"),
    file.path(.pathofm_user_model_dir(), filename)
  )
}

.pathofm_load_object <- function(foundation_model, refresh = FALSE) {
  key <- paste(foundation_model, .pathofm_user_model_dir(), sep = "::")
  if (!refresh && exists(key, envir = .pathofm_cache, inherits = FALSE)) {
    return(get(key, envir = .pathofm_cache, inherits = FALSE))
  }
  candidates <- .pathofm_object_candidates(foundation_model)
  path <- candidates[nzchar(candidates) & file.exists(candidates)][1L]
  if (!length(path) || is.na(path) || !file.exists(path)) {
    stop(
      "No fitted ", foundation_model, " PathoFMPred object is installed. ",
      if (identical(foundation_model, "TITAN")) {
        paste0(
          "TITAN-derived fitted parameters are not distributed in the public package. ",
          "Create an object with create_pathofmpred_object() or use the ",
          "access-controlled private package."
        )
      } else {
        paste0(
          "Obtain a collection under its applicable terms, then place ",
          "pathofmpred_", .pathofm_slug(foundation_model), ".rds in ",
          .pathofm_user_model_dir(), ". Alternatively, call ",
          "fetch_pathofmpred_models() with a source URL and SHA-256 ",
          "provided by the rights holder."
        )
      },
      call. = FALSE
    )
  }
  object <- readRDS(path)
  validate_pathofmpred_object(object)
  assign(key, object, envir = .pathofm_cache)
  object
}

.pathofm_read_model <- function(model_id, foundation_model) {
  object <- .pathofm_load_object(foundation_model)
  artifact <- object$models[[model_id]]
  if (is.null(artifact)) {
    stop(
      "Model ", model_id, " is listed in the atlas registry but is not present ",
      "in the installed ", foundation_model, " object. Supply an authorized ",
      "collection containing this model.",
      call. = FALSE
    )
  }
  artifact
}

#' Download user-authorized PathoFMPred model collections
#'
#' The package does not provide model files, download URLs, or checksums. A
#' download occurs only after the user supplies a source and its SHA-256 digest.
#' Users must confirm that they may access and use the collection under its
#' applicable terms. TITAN-derived objects are not supported by this helper;
#' authorized local TITAN objects can be used directly.
#'
#' @param foundation_model One or both of `GigaSSL` and `ProvGigaPath`.
#' @param destination Directory in which model collections will be stored.
#' @param base_url User-supplied HTTPS release directory, or a `file://` URL
#'   for an authorized local mirror. There is no package default.
#' @param expected_sha256 Named character vector containing the expected
#'   64-character SHA-256 digest for every requested representation. Obtain
#'   these digests independently from the authorized source.
#' @param overwrite Replace an existing file.
#' @param quiet Passed to [utils::download.file()].
#' @return Named paths to the downloaded model collections.
#' @export
fetch_pathofmpred_models <- function(
    foundation_model,
    destination = .pathofm_user_model_dir(),
    base_url,
    expected_sha256,
    overwrite = FALSE, quiet = FALSE) {
  foundation_model <- unique(match.arg(
    foundation_model, c("GigaSSL", "ProvGigaPath"), several.ok = TRUE
  ))
  if (missing(base_url) || length(base_url) != 1L || is.na(base_url) ||
      !grepl("^(https|file)://", base_url)) {
    stop("base_url must be an explicit HTTPS or file:// source supplied by the user.",
         call. = FALSE)
  }
  if (missing(expected_sha256) || is.null(names(expected_sha256)) ||
      !all(foundation_model %in% names(expected_sha256)) ||
      anyNA(expected_sha256[foundation_model]) ||
      any(!grepl("^[0-9a-fA-F]{64}$", expected_sha256[foundation_model]))) {
    stop("expected_sha256 must contain a valid named SHA-256 digest for every requested representation.",
         call. = FALSE)
  }
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(character(length(foundation_model)), foundation_model)
  for (model in foundation_model) {
    filename <- paste0("pathofmpred_", .pathofm_slug(model), ".rds")
    target <- file.path(destination, filename)
    if (file.exists(target) && !isTRUE(overwrite)) {
      observed <- digest::digest(file = target, algo = "sha256",
                                 serialize = FALSE)
      if (!identical(tolower(unname(observed)),
                     tolower(unname(expected_sha256[[model]])))) {
        stop(
          "SHA-256 verification failed for the existing ", model,
          " object. Set overwrite = TRUE to download a clean copy.",
          call. = FALSE
        )
      }
      object <- readRDS(target)
      validate_pathofmpred_object(object)
      paths[[model]] <- normalizePath(target)
      next
    }
    url <- paste0(sub("/$", "", base_url), "/", filename)
    temporary <- tempfile(pattern = paste0(filename, "-"), tmpdir = destination)
    on.exit(unlink(temporary), add = TRUE)
    status <- utils::download.file(url, temporary, mode = "wb", quiet = quiet)
    if (!identical(status, 0L) || !file.exists(temporary)) {
      stop("Download failed for ", model, " from ", url, ".", call. = FALSE)
    }
    observed <- digest::digest(file = temporary, algo = "sha256",
                               serialize = FALSE)
    if (!identical(tolower(unname(observed)),
                   tolower(unname(expected_sha256[[model]])))) {
      stop("SHA-256 verification failed for ", model, ".", call. = FALSE)
    }
    object <- readRDS(temporary)
    validate_pathofmpred_object(object)
    if (!identical(object$foundation_model, model)) {
      stop("Downloaded object identifies itself as ", object$foundation_model,
           ", not ", model, ".", call. = FALSE)
    }
    backup <- NULL
    if (file.exists(target) && isTRUE(overwrite)) {
      backup <- tempfile(pattern = paste0(filename, "-backup-"),
                         tmpdir = destination)
      if (!file.rename(target, backup)) {
        stop("Could not preserve the existing model before overwrite: ",
             target, ".", call. = FALSE)
      }
    }
    if (!file.rename(temporary, target)) {
      if (!is.null(backup) && !file.rename(backup, target)) {
        stop("Could not install the validated download or restore the existing model. The backup is at ",
             backup, ".", call. = FALSE)
      }
      stop("Could not move the validated download to ", target, ".", call. = FALSE)
    }
    if (!is.null(backup)) unlink(backup)
    paths[[model]] <- normalizePath(target)
  }
  rm(list = ls(envir = .pathofm_cache), envir = .pathofm_cache)
  invisible(paths)
}
