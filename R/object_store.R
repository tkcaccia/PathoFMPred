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
        if (identical(foundation_model, "ProvGigaPath")) {
          "Call fetch_pathofmpred_models(\"ProvGigaPath\", destination = <your chosen directory>) to download the optional verified research collection, or provide an authorized local object."
        } else {
          paste0(
            "Obtain an authorized collection, then place ",
            "pathofmpred_", .pathofm_slug(foundation_model), ".rds in ",
            .pathofm_user_model_dir(), ". Alternatively, call ",
            "fetch_pathofmpred_models() with an explicit source URL and SHA-256."
          )
        }
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
#' The source package contains no fitted models. For `ProvGigaPath`, the
#' default is a versioned GitHub release asset whose SHA-256 digest is pinned
#' in the package. It is downloaded only when this function is called. Other
#' representations require an explicit authorized URL and SHA-256 digest.
#' TITAN-derived objects are not downloadable through this helper.
#'
#' @param foundation_model `ProvGigaPath` or `GigaSSL`.
#' @param destination Explicit directory in which model collections will be
#'   stored. No directory is created unless the function is called with this
#'   argument.
#' @param base_url Optional HTTPS release directory, or a `file://` URL for
#'   an authorized local mirror. The pinned Prov-GigaPath release is the
#'   default when this is omitted.
#' @param expected_sha256 Optional named character vector containing a
#'   64-character SHA-256 digest for every requested representation. The
#'   pinned Prov-GigaPath digest is the default for its release.
#' @param overwrite Replace an existing file.
#' @param quiet Passed to [utils::download.file()].
#' @return Named paths to the downloaded model collections.
#' @export
fetch_pathofmpred_models <- function(
    foundation_model = "ProvGigaPath",
    destination = NULL,
    base_url = NULL,
    expected_sha256 = NULL,
    overwrite = FALSE, quiet = FALSE) {
  foundation_model <- unique(match.arg(
    foundation_model, c("GigaSSL", "ProvGigaPath"), several.ok = TRUE
  ))
  pinned <- list(
    base_url = paste0(
      "https://github.com/tkcaccia/PathoFMPred/releases/download/",
      "prov-v0.4.0-chol-20261026"
    ),
    sha256 = "c0918ca947d0fbde89faec2a51b06288d5f75b407372a53e8c89acabe7e59dd4"
  )
  if (is.null(base_url) && identical(foundation_model, "ProvGigaPath")) {
    base_url <- pinned$base_url
  }
  if (is.null(expected_sha256) && identical(foundation_model, "ProvGigaPath") &&
      identical(base_url, pinned$base_url)) {
    expected_sha256 <- c(ProvGigaPath = pinned$sha256)
  }
  if (is.null(base_url) || length(base_url) != 1L || is.na(base_url) ||
      !grepl("^(https|file)://", base_url)) {
    stop("base_url must be an HTTPS or file:// source; only ProvGigaPath has a pinned default.",
         call. = FALSE)
  }
  if (is.null(expected_sha256) || is.null(names(expected_sha256)) ||
      !all(foundation_model %in% names(expected_sha256)) ||
      anyNA(expected_sha256[foundation_model]) ||
      any(!grepl("^[0-9a-fA-F]{64}$", expected_sha256[foundation_model]))) {
    stop("expected_sha256 must contain a valid named SHA-256 digest for every requested representation.",
         call. = FALSE)
  }
  if (is.null(destination) || length(destination) != 1L ||
      !is.character(destination) || is.na(destination) ||
      !nzchar(destination)) {
    stop("destination must be an explicit directory selected by the caller.",
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
