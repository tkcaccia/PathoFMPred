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

.pathofm_load_object <- function(foundation_model, allow_example = TRUE,
                                  refresh = FALSE) {
  key <- paste0(foundation_model, "::", allow_example)
  if (!refresh && exists(key, envir = .pathofm_cache, inherits = FALSE)) {
    return(get(key, envir = .pathofm_cache, inherits = FALSE))
  }
  candidates <- .pathofm_object_candidates(foundation_model)
  path <- candidates[nzchar(candidates) & file.exists(candidates)][1L]
  if ((!length(path) || is.na(path)) && allow_example &&
      identical(foundation_model, "GigaSSL")) {
    example <- system.file("models", "pathofmpred_gigassl_example.rds",
                           package = "PathoFMPred")
    if (nzchar(example) && file.exists(example)) path <- example
  }
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
          "Call fetch_pathofmpred_models(\"", foundation_model,
          "\") explicitly after package installation."
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
      "in the installed ", foundation_model, " object. Install the full object ",
      "with fetch_pathofmpred_models() if it is publicly available.",
      call. = FALSE
    )
  }
  artifact
}

#' Download optional full PathoFMPred model collections
#'
#' Downloads occur only after an explicit user call. The public package does
#' not download model files during installation, loading, checking, or vignette
#' construction. TITAN-derived fitted parameters are deliberately unavailable
#' through this function.
#'
#' @param foundation_model One or both of `GigaSSL` and `ProvGigaPath`.
#' @param destination Directory in which model collections will be stored.
#' @param base_url Release directory containing the model files. Advanced users
#'   may point this at an institutional mirror.
#' @param overwrite Replace an existing file.
#' @param quiet Passed to [utils::download.file()].
#' @return Named paths to the downloaded model collections.
#' @export
fetch_pathofmpred_models <- function(
    foundation_model = c("GigaSSL", "ProvGigaPath"),
    destination = .pathofm_user_model_dir(),
    base_url = getOption(
      "PathoFMPred.model_base_url",
      "https://github.com/tkcaccia/PathoFMPred/releases/download/models-v1"
    ),
    overwrite = FALSE, quiet = FALSE) {
  foundation_model <- unique(match.arg(
    foundation_model, c("GigaSSL", "ProvGigaPath"), several.ok = TRUE
  ))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(character(length(foundation_model)), foundation_model)
  for (model in foundation_model) {
    filename <- paste0("pathofmpred_", .pathofm_slug(model), ".rds")
    target <- file.path(destination, filename)
    if (file.exists(target) && !isTRUE(overwrite)) {
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
    object <- readRDS(temporary)
    validate_pathofmpred_object(object)
    if (!identical(object$foundation_model, model)) {
      stop("Downloaded object identifies itself as ", object$foundation_model,
           ", not ", model, ".", call. = FALSE)
    }
    if (file.exists(target) && isTRUE(overwrite)) unlink(target)
    if (!file.rename(temporary, target)) {
      stop("Could not move the validated download to ", target, ".", call. = FALSE)
    }
    paths[[model]] <- normalizePath(target)
  }
  rm(list = ls(envir = .pathofm_cache), envir = .pathofm_cache)
  invisible(paths)
}
