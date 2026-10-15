#' Read registries from locally available model collections
#'
#' The public package contains no fitted-model registry. This function reads
#' metadata from collections that the user has installed under the applicable
#' access terms. An empty data frame is returned when none is installed.
#'
#' @return A data frame with one row per locally available fitted model.
#' @export
pathofm_model_registry <- function() {
  path <- system.file("extdata", "model_registry.csv", package = "PathoFMPred")
  if (nzchar(path)) {
    registry <- utils::read.csv(path, check.names = FALSE,
                                stringsAsFactors = FALSE)
    if (!"foundation_model" %in% names(registry)) {
      registry$foundation_model <- "TITAN"
    }
    return(registry)
  }
  candidates <- c("TITAN", "GigaSSL", "ProvGigaPath")
  registries <- lapply(candidates, function(model) {
    object <- try(.pathofm_load_object(model), silent = TRUE)
    if (inherits(object, "try-error")) return(NULL)
    registry <- object$registry
    if (!"foundation_model" %in% names(registry)) {
      registry$foundation_model <- model
    }
    registry
  })
  registries <- Filter(Negate(is.null), registries)
  if (!length(registries)) {
    return(data.frame(
      model_id = character(), foundation_model = character(),
      cancer_type = character(), endpoint = character(),
      outcome_type = character(), family = character(),
      model_evidence_tier = character(), stringsAsFactors = FALSE
    ))
  }
  columns <- unique(unlist(lapply(registries, names), use.names = FALSE))
  registries <- lapply(registries, function(registry) {
    for (column in setdiff(columns, names(registry))) registry[[column]] <- NA
    registry[, columns, drop = FALSE]
  })
  rownames(registries) <- NULL
  do.call(rbind, registries)
}

.pathofm_predictable_mask <- function(registry) {
  if ("matched_representation_effect_threshold_crossing" %in% names(registry)) {
    matched <- registry$matched_representation_effect_threshold_crossing
    keep <- matched %in% TRUE
    unmatched_titan <- is.na(matched) &
      registry$foundation_model %in% "TITAN"
    # TITAN objects outside the three-representation common cohort come from
    # the supporting permutation/FDR-qualified TITAN candidate registry.
    keep[unmatched_titan] <- TRUE
    return(keep)
  }
  if ("representation_effect_threshold_crossing" %in% names(registry)) {
    return(registry$representation_effect_threshold_crossing %in% TRUE)
  }
  rep(FALSE, nrow(registry))
}

#' List locally available foundation-model representations
#'
#' @return A character vector. Availability means that at least one fitted
#'   research model collection is installed; it does not imply external
#'   validation or permission to redistribute the collection.
#' @export
available_foundation_models <- function() {
  candidates <- c("TITAN", "GigaSSL", "ProvGigaPath")
  candidates[vapply(candidates, function(model) {
    !inherits(try(.pathofm_load_object(model), silent = TRUE), "try-error")
  }, logical(1))]
}

#' @rdname pathofm_model_registry
#' @export
titan_model_registry <- function() {
  registry <- pathofm_model_registry()
  registry[registry$foundation_model == "TITAN", , drop = FALSE]
}

#' List available cancer-specific models
#'
#' @param foundation_model One or more of `TITAN`, `GigaSSL`, or
#'   `ProvGigaPath`.
#' @param cancer Optional TCGA cancer code or vector of codes.
#' @param outcome_type Optional `continuous` or `binary` filter.
#' @param family Optional endpoint-family filter.
#' @param include_limited_evidence Include limited-sample models only when they
#'   crossed the representation-specific matched effect threshold. Models
#'   tested below the threshold never enter prediction or reporting. The
#'   default is `FALSE`.
#' @return A filtered model-registry data frame.
#' @export
available_pathofm_models <- function(foundation_model = NULL, cancer = NULL,
                                     outcome_type = NULL, family = NULL,
                                     include_limited_evidence = FALSE) {
  registry <- pathofm_model_registry()
  if (!is.null(foundation_model)) {
    registry <- registry[registry$foundation_model %in% foundation_model, , drop = FALSE]
  }
  if (nrow(registry)) {
    installed <- lapply(unique(registry$foundation_model), function(model) {
      object <- try(.pathofm_load_object(model), silent = TRUE)
      if (inherits(object, "try-error")) character() else names(object$models)
    })
    names(installed) <- unique(registry$foundation_model)
    locally_available <- vapply(seq_len(nrow(registry)), function(index) {
      registry$model_id[index] %in% installed[[registry$foundation_model[index]]]
    }, logical(1))
    registry <- registry[locally_available, , drop = FALSE]
  }
  registry <- registry[.pathofm_predictable_mask(registry), , drop = FALSE]
  if (!isTRUE(include_limited_evidence)) {
    limited <- if ("model_evidence_tier" %in% names(registry)) {
      grepl("^exploratory_limited", registry$model_evidence_tier)
    } else {
      rep(FALSE, nrow(registry))
    }
    registry <- registry[!limited, , drop = FALSE]
  }
  if (!is.null(cancer)) registry <- registry[registry$cancer_type %in% toupper(cancer), , drop = FALSE]
  if (!is.null(outcome_type)) registry <- registry[registry$outcome_type %in% outcome_type, , drop = FALSE]
  if (!is.null(family)) registry <- registry[registry$family %in% family, , drop = FALSE]
  rownames(registry) <- NULL
  registry
}

#' @rdname available_pathofm_models
#' @export
available_titan_models <- function(cancer = NULL, outcome_type = NULL,
                                   family = NULL,
                                   include_limited_evidence = FALSE) {
  available_pathofm_models("TITAN", cancer, outcome_type, family,
                           include_limited_evidence)
}

.pathofm_model_path <- function(model_id, foundation_model = "TITAN") {
  stop(
    ".pathofm_model_path() is obsolete. Models are stored in one validated ",
    "PathoFMPred object per representation; use .pathofm_read_model().",
    call. = FALSE
  )
}

.pathofm_reference <- function(foundation_model = NULL) {
  if (!is.null(foundation_model)) {
    object <- try(.pathofm_load_object(foundation_model), silent = TRUE)
    if (!inherits(object, "try-error") && !is.null(object$prediction_reference)) {
      return(object$prediction_reference)
    }
  }
  path <- system.file("extdata", "prediction_reference.rds", package = "PathoFMPred")
  if (!nzchar(path)) return(list())
  readRDS(path)
}
