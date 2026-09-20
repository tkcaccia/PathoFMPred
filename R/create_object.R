.pathofm_read_table <- function(x, argument) {
  if (is.data.frame(x) || is.matrix(x)) {
    return(as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    return(utils::read.csv(x, check.names = FALSE, stringsAsFactors = FALSE))
  }
  stop(argument, " must be a data frame, matrix, or existing CSV path.",
       call. = FALSE)
}

#' Define model-development settings
#'
#' @param components Candidate PLS component counts.
#' @param outer_folds,inner_folds Numbers of nested cross-validation folds.
#' @param repeats Number of independently seeded nested evaluations.
#' @param scaling,method,backend Arguments passed to fastPLS.
#' @param binary_selection Inner-CV objective for binary outcomes. AUROC is the
#'   default because the resulting score ranking is independent of prevalence-
#'   specific class priors.
#' @param continuous_selection Inner-CV objective for continuous outcomes.
#' @param seed Initial random seed.
#' @param minimum_continuous Minimum labelled patients for a continuous model.
#' @param minimum_binary_class Minimum patients in each binary class.
#' @param rsvd_oversample,rsvd_power Optional overrides. `NULL` preserves the
#'   installed fastPLS defaults.
#' @return A settings list for [create_pathofmpred_object()].
#' @export
pathofmpred_control <- function(
    components = 1:20, outer_folds = 5L, inner_folds = 5L, repeats = 1L,
    scaling = "centering", method = "plssvd", backend = NULL,
    binary_selection = "AUROC", continuous_selection = "auto", seed = 1L,
    minimum_continuous = 50L, minimum_binary_class = 20L,
    rsvd_oversample = NULL, rsvd_power = NULL) {
  components <- sort(unique(as.integer(components)))
  if (!length(components) || any(!is.finite(components)) || any(components < 1L)) {
    stop("components must contain positive integers.", call. = FALSE)
  }
  list(
    components = components,
    outer_folds = as.integer(outer_folds), inner_folds = as.integer(inner_folds),
    repeats = as.integer(repeats), scaling = scaling, method = method,
    backend = backend, binary_selection = binary_selection,
    continuous_selection = continuous_selection, seed = as.integer(seed),
    minimum_continuous = as.integer(minimum_continuous),
    minimum_binary_class = as.integer(minimum_binary_class),
    rsvd_oversample = rsvd_oversample, rsvd_power = rsvd_power
  )
}

.pathofm_fastpls_args <- function(control) {
  out <- list(scaling = control$scaling, method = control$method)
  if (!is.null(control$backend)) out$backend <- control$backend
  if (!is.null(control$rsvd_oversample)) {
    out$rsvd_oversample <- as.integer(control$rsvd_oversample)
  }
  if (!is.null(control$rsvd_power)) out$rsvd_power <- as.integer(control$rsvd_power)
  out
}

.pathofm_pool_rows <- function(x, ids, method) {
  split_rows <- split(seq_len(nrow(x)), ids)
  pooled <- t(vapply(split_rows, function(index) {
    if (method == "mean") colMeans(x[index, , drop = FALSE]) else
      apply(x[index, , drop = FALSE], 2L, stats::median)
  }, numeric(ncol(x))))
  colnames(pooled) <- colnames(x)
  pooled
}

.pathofm_auc <- function(y, score) {
  keep <- is.finite(score) & !is.na(y)
  y <- as.integer(as.character(y[keep])); score <- score[keep]
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) /
    (n1 * n0)
}

.pathofm_best_threshold <- function(y, score) {
  y <- as.integer(as.character(y)); keep <- is.finite(score) & !is.na(y)
  y <- y[keep]; score <- score[keep]
  candidates <- sort(unique(score))
  if (!length(candidates)) return(NA_real_)
  ba <- vapply(candidates, function(cut) {
    call <- as.integer(score >= cut)
    0.5 * (mean(call[y == 1L] == 1L) + mean(call[y == 0L] == 0L))
  }, numeric(1))
  candidates[which.max(ba)]
}

.pathofm_binary_scores <- function(fit, component) {
  index <- match(component, as.integer(dimnames(fit$lda_scores)[[3L]]))
  if (is.na(index)) index <- match(component, fit$ncomp)
  drop(fit$lda_scores[, 2L, index] - fit$lda_scores[, 1L, index])
}

.pathofm_fit_endpoint <- function(X, y, type, control, endpoint, feature_names,
                                  class_labels = NULL) {
  extra <- .pathofm_fastpls_args(control)
  selected <- integer(control$repeats)
  repeat_metrics <- vector("list", control$repeats)
  for (repeat_index in seq_len(control$repeats)) {
    common <- c(list(
      Xdata = X, Ydata = y, ncomp = control$components,
      kfold_outer = control$outer_folds, kfold_inner = control$inner_folds,
      seed = control$seed + repeat_index - 1L, perm.test = FALSE,
      return_splits = TRUE
    ), extra)
    if (type == "binary") {
      common$classifier <- "lda"; common$selection <- control$binary_selection
    } else {
      common$selection <- control$continuous_selection
    }
    cv <- do.call(fastPLS::pls.double.cv, common)
    selected[repeat_index] <- as.integer(stats::median(as.integer(cv$bcomp)))
    if (type == "binary") {
      repeat_metrics[[repeat_index]] <- data.frame(
        repeat_index = repeat_index, n = length(y), positive = sum(y == "1"),
        negative = sum(y == "0"), AUROC = as.numeric(cv$AUROC),
        balanced_accuracy = as.numeric(cv$balanced_accuracy),
        stringsAsFactors = FALSE
      )
    } else {
      prediction <- as.numeric(cv$Ypred)
      repeat_metrics[[repeat_index]] <- data.frame(
        repeat_index = repeat_index, n = length(y), Q2 = as.numeric(cv$Q2Y),
        RMSE = sqrt(mean((y - prediction)^2)),
        spearman = suppressWarnings(stats::cor(y, prediction, method = "spearman")),
        stringsAsFactors = FALSE
      )
    }
  }
  tune_args <- c(list(
    Xdata = X, Ydata = y, ncomp = control$components,
    kfold = control$inner_folds, seed = control$seed + 100000L,
    fit = FALSE, return_splits = TRUE
  ), extra)
  if (type == "binary") {
    tune_args$classifier <- "lda"; tune_args$selection <- control$binary_selection
  } else {
    tune_args$selection <- control$continuous_selection
  }
  tuning <- do.call(fastPLS::pls.single.cv, tune_args)
  ncomp <- as.integer(tuning$best_ncomp)
  fit_args <- c(list(
    Xtrain = X, Ytrain = y, ncomp = ncomp, fit = TRUE,
    return_loadings = TRUE, seed = control$seed + 200000L
  ), extra)
  if (type == "binary") fit_args$classifier <- "lda"
  fitted <- do.call(fastPLS::pls, fit_args)
  threshold <- NA_real_
  if (type == "binary") {
    threshold <- .pathofm_best_threshold(y, .pathofm_binary_scores(tuning, ncomp))
  }
  model_id <- paste0(
    gsub("[^A-Za-z0-9]+", "_", endpoint), "__",
    sprintf("%08x", sum(utf8ToInt(endpoint) * seq_along(utf8ToInt(endpoint))))
  )
  artifact <- list(
    model_id = model_id, endpoint = endpoint, outcome_type = type,
    model = fitted, operating_threshold = threshold, ncomp = ncomp,
    feature_names = feature_names, training_feature_mean = colMeans(X),
    training_feature_min = apply(X, 2L, min),
    training_feature_max = apply(X, 2L, max),
    training_n = length(y),
    training_positive = if (type == "binary") sum(y == "1") else NA_integer_,
    training_negative = if (type == "binary") sum(y == "0") else NA_integer_,
    class_labels = if (type == "binary") class_labels else NULL,
    nested_cv = do.call(rbind, repeat_metrics),
    selected_components_across_repeats = selected,
    fastPLS_version = as.character(utils::packageVersion("fastPLS")),
    fastPLS_remote_sha = as.character(utils::packageDescription("fastPLS")$RemoteSha)
  )
  list(artifact = artifact, registry = data.frame(
    model_id = model_id, endpoint = endpoint, outcome_type = type,
    n = length(y), positive = artifact$training_positive,
    negative = artifact$training_negative, ncomp = ncomp,
    negative_label = if (type == "binary") class_labels$negative else NA_character_,
    positive_label = if (type == "binary") class_labels$positive else NA_character_,
    stringsAsFactors = FALSE
  ))
}

.pathofm_positive_class <- function(values, endpoint, positive_class) {
  observed <- unique(as.character(values))
  if (length(observed) != 2L) {
    stop("Binary outcome ", endpoint, " must contain exactly two observed classes.",
         call. = FALSE)
  }
  selected <- NULL
  if (!is.null(positive_class)) {
    if (!is.null(names(positive_class)) && endpoint %in% names(positive_class)) {
      selected <- as.character(positive_class[[endpoint]])
    } else if (length(positive_class) == 1L &&
               (is.null(names(positive_class)) || !nzchar(names(positive_class)))) {
      selected <- as.character(positive_class)
    }
  }
  canonical <- sort(observed)
  if (is.null(selected) && identical(canonical, c("0", "1"))) selected <- "1"
  if (is.null(selected) && identical(canonical, c("FALSE", "TRUE"))) selected <- "TRUE"
  if (is.null(selected)) {
    stop(
      "positive_class must identify the event class for binary outcome ",
      endpoint, ". Use a named character vector when fitting several outcomes.",
      call. = FALSE
    )
  }
  if (length(selected) != 1L || is.na(selected) || !selected %in% observed) {
    stop("positive_class for ", endpoint, " must be one of: ",
         paste(observed, collapse = ", "), ".", call. = FALSE)
  }
  list(negative = setdiff(observed, selected)[[1L]], positive = selected)
}

#' Build a portable PathoFMPred model collection
#'
#' Feature and outcome tables are linked only through `id_column`. Outcome IDs
#' must be unique. Feature IDs may repeat, for example when a patient has more
#' than one slide, and are pooled before model development.
#'
#' @param feature_table Data frame, matrix, or CSV path containing identifiers
#'   and numeric foundation-model features.
#' @param outcome_table Data frame or CSV path containing one unique row per
#'   identifier and one or more outcomes.
#' @param id_column Name of the identifier column present in both tables.
#' @param aggregation How repeated feature rows are pooled: `mean` or `median`.
#' @param foundation_model Label stored in the resulting object.
#' @param feature_columns Optional numeric feature columns. By default all
#'   columns except `id_column` are used.
#' @param outcome_columns Optional outcomes to fit. By default all columns
#'   except `id_column` are considered.
#' @param outcome_types Optional named character vector containing `continuous`
#'   or `binary`. Otherwise numeric two-level outcomes are binary and other
#'   numeric outcomes are continuous.
#' @param positive_class Event-class label for binary outcomes. Supply either
#'   one value when fitting a single binary outcome or a named character vector
#'   keyed by endpoint. Zero/one and FALSE/TRUE outcomes default to one and TRUE;
#'   all other binary labels require an explicit value.
#' @param control Settings created by [pathofmpred_control()].
#' @param output_file Optional `.rds` destination.
#' @return A `PathoFMPredObject` containing fitted models, a registry, nested-CV
#'   summaries, and a linkage/missingness audit.
#' @export
create_pathofmpred_object <- function(
    feature_table, outcome_table, id_column,
    aggregation = c("mean", "median"), foundation_model = "custom",
    feature_columns = NULL, outcome_columns = NULL, outcome_types = NULL,
    positive_class = NULL, control = pathofmpred_control(), output_file = NULL) {
  aggregation <- match.arg(aggregation)
  features <- .pathofm_read_table(feature_table, "feature_table")
  outcomes <- .pathofm_read_table(outcome_table, "outcome_table")
  if (length(id_column) != 1L || !is.character(id_column) ||
      !id_column %in% names(features) || !id_column %in% names(outcomes)) {
    stop("id_column must name one column present in both tables.", call. = FALSE)
  }
  outcome_ids <- as.character(outcomes[[id_column]])
  feature_ids <- as.character(features[[id_column]])
  if (anyNA(outcome_ids) || any(!nzchar(outcome_ids)) || anyDuplicated(outcome_ids)) {
    stop("outcome_table IDs must be non-missing and unique.", call. = FALSE)
  }
  if (anyNA(feature_ids) || any(!nzchar(feature_ids))) {
    stop("feature_table IDs must be non-missing; repeated IDs are allowed.",
         call. = FALSE)
  }
  if (is.null(feature_columns)) feature_columns <- setdiff(names(features), id_column)
  if (!length(feature_columns) || any(!feature_columns %in% names(features))) {
    stop("feature_columns must select columns in feature_table.", call. = FALSE)
  }
  numeric_features <- vapply(features[feature_columns], is.numeric, logical(1))
  if (!all(numeric_features)) {
    stop("All selected feature columns must be numeric: ",
         paste(feature_columns[!numeric_features], collapse = ", "), ".",
         call. = FALSE)
  }
  Xraw <- as.matrix(features[feature_columns]); storage.mode(Xraw) <- "double"
  if (any(!is.finite(Xraw))) stop("Feature values must be finite.", call. = FALSE)
  pooled <- .pathofm_pool_rows(Xraw, feature_ids, aggregation)
  pooled_ids <- rownames(pooled)
  common_ids <- intersect(pooled_ids, outcome_ids)
  if (length(common_ids) <= 100L) {
    stop("Only ", length(common_ids), " matching unique IDs were found; more than 100 are required.",
         call. = FALSE)
  }
  pooled <- pooled[common_ids, , drop = FALSE]
  outcomes <- outcomes[match(common_ids, outcome_ids), , drop = FALSE]
  if (is.null(outcome_columns)) outcome_columns <- setdiff(names(outcomes), id_column)
  if (!length(outcome_columns) || any(!outcome_columns %in% names(outcomes))) {
    stop("outcome_columns must select columns in outcome_table.", call. = FALSE)
  }
  audit <- list(
    matching_ids = length(common_ids),
    feature_rows = nrow(features), unique_feature_ids = length(unique(feature_ids)),
    outcome_rows = nrow(outcomes), aggregation = aggregation,
    unmatched_feature_ids = setdiff(pooled_ids, outcome_ids),
    unmatched_outcome_ids = setdiff(outcome_ids, pooled_ids),
    outcome_missingness = data.frame(
      endpoint = outcome_columns,
      missing = vapply(outcomes[outcome_columns], function(x) sum(is.na(x)), integer(1)),
      observed = vapply(outcomes[outcome_columns], function(x) sum(!is.na(x)), integer(1)),
      stringsAsFactors = FALSE
    )
  )
  message(
    length(common_ids), " IDs matched; ", length(audit$unmatched_feature_ids),
    " feature IDs and ", length(audit$unmatched_outcome_ids),
    " outcome IDs were unmatched. See object$data_audit for IDs and endpoint missingness."
  )
  models <- list(); registry <- list(); skipped <- list()
  for (endpoint in outcome_columns) {
    raw <- outcomes[[endpoint]]; keep <- !is.na(raw)
    values <- raw[keep]; X <- pooled[keep, , drop = FALSE]
    declared <- if (!is.null(outcome_types) && endpoint %in% names(outcome_types)) {
      outcome_types[[endpoint]]
    } else if ((is.numeric(values) || is.logical(values) || is.factor(values) || is.character(values)) &&
               length(unique(values)) == 2L) "binary" else "continuous"
    if (!declared %in% c("continuous", "binary")) {
      stop("outcome_types for ", endpoint, " must be continuous or binary.", call. = FALSE)
    }
    class_labels <- NULL
    if (declared == "continuous") {
      if (!is.numeric(values) || length(values) < control$minimum_continuous) {
        skipped[[endpoint]] <- "continuous outcome was nonnumeric or below minimum_continuous"
        next
      }
      y <- as.numeric(values)
    } else {
      labels <- .pathofm_positive_class(values, endpoint, positive_class)
      y <- factor(ifelse(as.character(values) == labels$positive, 1L, 0L),
                  levels = c(0L, 1L))
      class_labels <- labels
      if (min(table(y)) < control$minimum_binary_class) {
        skipped[[endpoint]] <- "binary outcome had fewer than two classes or was below minimum_binary_class"
        next
      }
    }
    fitted <- .pathofm_fit_endpoint(
      X, y, declared, control, endpoint, colnames(pooled), class_labels
    )
    models[[fitted$artifact$model_id]] <- fitted$artifact
    registry[[length(registry) + 1L]] <- fitted$registry
  }
  if (!length(models)) stop("No outcome met the model-development eligibility rules.", call. = FALSE)
  object <- structure(list(
    format_version = "1.0.0", foundation_model = foundation_model,
    feature_names = colnames(pooled), registry = do.call(rbind, registry),
    models = models, data_audit = audit, skipped_outcomes = skipped,
    control = control,
    provenance = list(
      package = "PathoFMPred",
      package_version = as.character(utils::packageVersion("PathoFMPred")),
      fastPLS_version = as.character(utils::packageVersion("fastPLS")),
      contains_feature_rows = FALSE, contains_outcome_rows = FALSE,
      note = "Fitted parameters and audit metadata only; source feature and outcome rows are not stored."
    )
  ), class = c("PathoFMPredObject", "list"))
  validate_pathofmpred_object(object)
  if (!is.null(output_file)) saveRDS(object, output_file, compress = "xz")
  object
}

#' @export
print.PathoFMPredObject <- function(x, ...) {
  cat("PathoFMPred model collection\n")
  cat("  representation:", x$foundation_model, "\n")
  cat("  fitted outcomes:", length(x$models), "\n")
  cat("  feature dimensions:", length(x$feature_names), "\n")
  cat("  matched IDs:", x$data_audit$matching_ids %||% NA_integer_, "\n")
  invisible(x)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
