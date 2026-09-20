.validate_runtime <- function(artifact) {
  installed <- as.character(utils::packageVersion("fastPLS"))
  required <- as.character(artifact$fastPLS_version)
  if (nzchar(required) && !identical(installed, required)) {
    stop("Model ", artifact$model_id, " requires fastPLS ", required,
         "; installed version is ", installed, ".", call. = FALSE)
  }
  required_sha <- as.character(artifact$fastPLS_remote_sha)
  installed_sha <- as.character(utils::packageDescription("fastPLS")$RemoteSha)
  if (length(required_sha) && !is.na(required_sha) && nzchar(required_sha) &&
      (!length(installed_sha) || is.na(installed_sha) ||
       !identical(installed_sha, required_sha))) {
    stop(
      "Model ", artifact$model_id, " requires fastPLS Git revision ",
      required_sha, "; the installed package reports ",
      if (length(installed_sha) && !is.na(installed_sha) && nzchar(installed_sha))
        installed_sha else "no RemoteSha",
      ". Reinstall the recorded Git revision before inference.", call. = FALSE
    )
  }
}

.prepare_pathofm_input <- function(features, cancer, patient_id) {
  x <- as.data.frame(features, check.names = FALSE)
  if (!nrow(x)) stop("features must contain at least one row.", call. = FALSE)
  cancer <- toupper(as.character(cancer))
  if (length(cancer) == 1L) cancer <- rep(cancer, nrow(x))
  if (length(cancer) != nrow(x) || anyNA(cancer) || any(!nzchar(cancer))) {
    stop("cancer must contain one non-missing TCGA cancer code per feature row.", call. = FALSE)
  }
  if (is.null(patient_id)) {
    patient_id <- rownames(x)
    if (is.null(patient_id) || any(!nzchar(patient_id))) patient_id <- paste0("sample_", seq_len(nrow(x)))
  }
  patient_id <- as.character(patient_id)
  if (length(patient_id) != nrow(x) || anyNA(patient_id) || any(!nzchar(patient_id))) {
    stop("patient_id must contain one non-missing identifier per feature row.", call. = FALSE)
  }
  mixed <- tapply(cancer, patient_id, function(z) length(unique(z)))
  if (any(mixed > 1L)) stop("All slides for a patient_id must have the same cancer code.", call. = FALSE)
  list(features = x, cancer = cancer, patient_id = patient_id)
}

.pool_patient_features <- function(x, patient_id, cancer, feature_names,
                                   foundation_model) {
  missing <- setdiff(feature_names, names(x))
  if (length(missing)) {
    stop("Missing ", foundation_model, " feature columns: ",
         paste(utils::head(missing, 10L), collapse = ", "),
         if (length(missing) > 10L) " ..." else "", call. = FALSE)
  }
  X <- as.matrix(x[, feature_names, drop = FALSE])
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) stop("All foundation-model feature values must be finite numeric values.", call. = FALSE)
  key <- paste(cancer, patient_id, sep = "\r")
  ids <- unique(key)
  group <- match(key, ids)
  pooled <- rowsum(X, group = group, reorder = FALSE) / tabulate(group, nbins = length(ids))
  first <- match(ids, key)
  list(X = pooled, patient_id = patient_id[first], cancer = cancer[first],
       n_slides = tabulate(group, nbins = length(ids)))
}

.prediction_percentile <- function(model_id, value, reference,
                                   foundation_model = "TITAN") {
  ref <- reference[[paste(foundation_model, model_id, sep = "::")]]
  # Backward compatibility with the original TITAN-only reference object.
  if (is.null(ref) && identical(foundation_model, "TITAN")) {
    ref <- reference[[model_id]]
  }
  if (is.null(ref) || !length(ref)) return(rep(NA_real_, length(value)))
  vapply(value, function(z) {
    if (!is.finite(z)) return(NA_real_)
    100 * (findInterval(z, ref, all.inside = TRUE) - 0.5) / length(ref)
  }, numeric(1))
}

#' Apply supported models for one pathology foundation-model representation
#'
#' Multiple rows with the same `patient_id` are interpreted as multiple slides
#' and mean-pooled before prediction, matching model development.
#'
#' @param cancer One TCGA cancer code per row, or one code recycled to all rows.
#' @param foundation_model One of `TITAN`, `GigaSSL`, or `ProvGigaPath`.
#' @param features A matrix or data frame containing the correctly named
#'   features for `foundation_model`.
#' @param patient_id Optional identifier per row. Repeated identifiers are pooled.
#' @param outcome_type One or more of `continuous` and `binary`.
#' @param family Optional endpoint-family filter.
#' @param warn_ood Warn when over 5 percent of dimensions fall outside training ranges.
#' @param include_limited_evidence Explicitly opt in to limited-sample models
#'   that still crossed the representation-specific matched effect threshold.
#'   Models tested below the threshold never enter prediction or reporting.
#' @return A long data frame with one row per patient and applicable model.
#' @export
predict_pathofm <- function(cancer, features,
                            foundation_model = c("TITAN", "GigaSSL", "ProvGigaPath"),
                            patient_id = NULL,
                            outcome_type = c("continuous", "binary"),
                            family = NULL, warn_ood = TRUE,
                            include_limited_evidence = FALSE) {
  foundation_model <- match.arg(foundation_model)
  input <- .prepare_pathofm_input(features, cancer, patient_id)
  registry <- available_pathofm_models(
    foundation_model, unique(input$cancer), outcome_type, family,
    include_limited_evidence = include_limited_evidence
  )
  unsupported <- setdiff(unique(input$cancer), unique(registry$cancer_type))
  if (length(unsupported)) {
    all_models <- available_pathofm_models(
      foundation_model, unsupported, outcome_type, family,
      include_limited_evidence = TRUE
    )
    suffix <- if (nrow(all_models) && !isTRUE(include_limited_evidence)) {
      " Only exploratory/limited-evidence models match; set include_limited_evidence = TRUE for explicit research-only opt-in."
    } else ""
    stop("No predictable bundled models are available for: ",
         paste(unsupported, collapse = ", "), ".", suffix, call. = FALSE)
  }
  limited_count <- sum(grepl("^exploratory_limited",
                             registry$model_evidence_tier), na.rm = TRUE)
  if (limited_count && isTRUE(include_limited_evidence)) {
    warning(
      "Explicit opt-in includes ", limited_count,
      " exploratory/limited-evidence model(s). Inspect class-size or ",
      "continuous sample-size and repeat-stability metadata.", call. = FALSE
    )
  }
  first_artifact <- .pathofm_read_model(registry$model_id[1L], foundation_model)
  pooled <- .pool_patient_features(input$features, input$patient_id, input$cancer,
                                   first_artifact$feature_names, foundation_model)
  reference <- .pathofm_reference(foundation_model)
  output <- vector("list", nrow(registry))
  for (i in seq_len(nrow(registry))) {
    info <- registry[i, , drop = FALSE]
    info_value <- function(name, default = NA) {
      if (name %in% names(info)) info[[name]][[1L]] else default
    }
    artifact <- .pathofm_read_model(info$model_id, foundation_model)
    .validate_runtime(artifact)
    keep <- pooled$cancer == info$cancer_type
    X <- pooled$X[keep, , drop = FALSE]
    ood <- rep(NA_real_, nrow(X))
    if (!is.null(artifact$training_feature_min) && !is.null(artifact$training_feature_max)) {
      ood <- rowMeans(sweep(X, 2L, artifact$training_feature_min, `<`) |
                       sweep(X, 2L, artifact$training_feature_max, `>`))
      if (warn_ood && any(ood > 0.05)) {
        warning(info$model_id, ": ", sum(ood > 0.05),
                " patient(s) have >5% of features outside the TCGA training range.",
                call. = FALSE)
      }
    }
    pred <- stats::predict(artifact$model, X, raw_scores = identical(info$outcome_type, "binary"))
    if (identical(info$outcome_type, "binary")) {
      score <- drop(pred$LDA_scores[, 2, 1] - pred$LDA_scores[, 1, 1])
      value <- rep(NA_real_, length(score))
      operating_threshold <- artifact$operating_threshold
      if (is.null(operating_threshold) || !is.finite(operating_threshold)) {
        operating_threshold <- artifact$operating_threshold_median
      }
      if (is.null(operating_threshold) || !is.finite(operating_threshold)) {
        operating_threshold <- suppressWarnings(as.numeric(
          info_value("binary_operating_threshold")
        ))
      }
      if (is.null(operating_threshold) || !is.finite(operating_threshold)) {
        operating_threshold <- suppressWarnings(as.numeric(
          info_value("binary_operating_threshold_median")
        ))
      }
      if (length(operating_threshold) && is.finite(operating_threshold)) {
        cls <- as.character(as.integer(score >= operating_threshold))
        class_rule <- "locked training-derived score threshold"
      } else {
        cls <- as.character(pred$Ypred[, 1])
        operating_threshold <- NA_real_
        class_rule <- "fitted LDA empirical-prior fallback"
      }
      pct <- .prediction_percentile(info$model_id, score, reference,
                                    foundation_model)
      class_labels <- artifact$class_labels
      cls_label <- if (!is.null(class_labels) &&
                       all(c("negative", "positive") %in% names(class_labels))) {
        ifelse(cls == "1", class_labels$positive, class_labels$negative)
      } else {
        cls
      }
    } else {
      value <- if (length(dim(pred$Ypred)) == 3L) drop(pred$Ypred[, 1, 1]) else drop(pred$Ypred)
      score <- rep(NA_real_, length(value))
      cls <- rep(NA_character_, length(value))
      cls_label <- rep(NA_character_, length(value))
      operating_threshold <- NA_real_
      class_rule <- "not applicable"
      pct <- .prediction_percentile(info$model_id, value, reference,
                                    foundation_model)
    }
    primary_auc <- suppressWarnings(as.numeric(info_value("screen_auc")))
    if (!is.finite(primary_auc)) {
      primary_auc <- suppressWarnings(as.numeric(info_value("primary_auc")))
    }
    primary_pr_auc <- suppressWarnings(as.numeric(info_value("binary_pr_auc")))
    if (!is.finite(primary_pr_auc)) {
      primary_pr_auc <- suppressWarnings(as.numeric(info_value("primary_pr_auc")))
    }
    primary_prevalence <- suppressWarnings(as.numeric(info_value("binary_prevalence")))
    if (!is.finite(primary_prevalence)) {
      primary_prevalence <- suppressWarnings(as.numeric(
        info_value("binary_observed_tcga_prevalence")
      ))
    }
    if (!is.finite(primary_prevalence) && identical(info$outcome_type, "binary")) {
      primary_prevalence <- as.numeric(info$positive) / as.numeric(info$n)
    }
    no_skill_pr_auc <- suppressWarnings(as.numeric(
      info_value("binary_no_skill_pr_auc")
    ))
    if (!is.finite(no_skill_pr_auc)) no_skill_pr_auc <- primary_prevalence
    matched_crossing <- info_value(
      "matched_representation_effect_threshold_crossing"
    )
    predictability_basis <- if (!is.na(matched_crossing)) {
      if (identical(info$outcome_type, "continuous"))
        "matched patient-level Q2 >= 0.20" else
        "matched patient-level AUROC >= 0.60"
    } else if (identical(info$outcome_type, "continuous")) {
      "supporting TITAN screen Q2 >= 0.20 with permutation/FDR qualification"
    } else {
      paste(
        "supporting TITAN screen balanced accuracy >= 0.60 under the",
        "documented PLS-LDA rule with permutation/FDR qualification"
      )
    }
    warnings <- character()
    if (grepl("^exploratory_limited", info$model_evidence_tier)) {
      warnings <- c(
        warnings,
        paste0("Exploratory/limited evidence: ",
               info$limited_evidence_reason)
      )
    }
    if (!is.na(info$site_robustness_warning) &&
        nzchar(info$site_robustness_warning)) {
      warnings <- c(warnings, info$site_robustness_warning)
    }
    grouped_adequacy <- if (
      "matched_tss_grouped_fold_adequacy_flag" %in% names(info)
    ) info$matched_tss_grouped_fold_adequacy_flag[[1L]] else NA_character_
    if (!is.na(grouped_adequacy) && identical(grouped_adequacy, "sparse grouped folds")) {
      warnings <- c(
        warnings,
        paste0(
          "Grouped-fold adequacy warning: ",
          info$matched_tss_grouped_fold_adequacy_reasons[[1L]],
          ". This is sensitivity to grouping by a barcode-derived cohort-structure ",
          "variable, not site-level validation or evidence of generalisability."
        )
      )
    }
    decision_rule_status <- if (
      "decision_rule_sensitivity_status" %in% names(info)
    ) info$decision_rule_sensitivity_status[[1L]] else NA_character_
    if (identical(info$outcome_type, "binary") &&
        identical(class_rule, "fitted LDA empirical-prior fallback")) {
      warnings <- c(
        warnings,
        paste0(
          "This legacy binary object lacks a locked training-derived threshold; ",
          "the reported class uses the fitted empirical-prior LDA fallback. ",
          "Inspect decision-rule metadata and do not interpret the score rank ",
          "as a probability."
        )
      )
    }
    output[[i]] <- data.frame(
      patient_id = pooled$patient_id[keep], cancer = pooled$cancer[keep],
      foundation_model = foundation_model,
      n_slides = pooled$n_slides[keep], model_id = info$model_id,
      family = info$family, endpoint = info$endpoint, outcome_type = info$outcome_type,
      prediction = value, predicted_class = cls, lda_score = score,
      predicted_class_label = cls_label,
      operating_threshold = operating_threshold,
      binary_class_rule = class_rule,
      # reference_percentile is retained for backward compatibility. For a
      # binary endpoint this is a rank of an uncalibrated score, never a
      # probability; reference_rank is the preferred explicit field name.
      reference_percentile = pct, reference_rank = pct,
      rank_interpretation = if (identical(info$outcome_type, "binary"))
        "Reference rank, not probability" else
        "TCGA OOF prediction percentile (not probability)",
      rank_is_probability = if (identical(info$outcome_type, "binary"))
        FALSE else NA,
      calibration_status = if (identical(info$outcome_type, "binary"))
        "uncalibrated; no probability estimate" else
        "not applicable to continuous regression output",
      output_units = info$output_units,
      tier = info$tier,
      tier_origin = ifelse(
        is.na(info$tier_origin),
        "TITAN within-cancer permutation/FDR-qualified screening tier",
        info$tier_origin
      ),
      predictable_for_selected_cancer_and_representation = TRUE,
      predictability_rule = predictability_basis,
      representation_effect_threshold_crossing = TRUE,
      matched_representation_effect_threshold_crossing =
        matched_crossing,
      resource_selection_basis = ifelse(
        is.na(info$resource_selection_basis),
        "TITAN permutation/FDR-qualified candidate catalogue",
        info$resource_selection_basis
      ),
      representation_specific_qualification = ifelse(
        is.na(info$representation_specific_qualification),
        "TITAN within-cancer permutation/FDR-qualified candidate",
        info$representation_specific_qualification
      ),
      ncomp = info$ncomp, training_n = info$n,
      training_positive = info$positive, training_negative = info$negative,
      training_prevalence = if (identical(info$outcome_type, "binary"))
        as.numeric(info$positive) / as.numeric(info$n) else NA_real_,
      limited_class_size = if (identical(info$outcome_type, "binary"))
        min(as.numeric(info$positive), as.numeric(info$negative)) < 50 else FALSE,
      model_evidence_tier = info$model_evidence_tier,
      default_inference = info$default_inference,
      limited_evidence_reason = info$limited_evidence_reason,
      primary_screen_q2 = info$screen_q2,
      primary_screen_rmse = info$screen_rmse,
      primary_screen_spearman = info$screen_spearman,
      primary_screen_balanced_accuracy = info$screen_balanced_accuracy,
      primary_screen_auc = primary_auc,
      primary_screen_pr_auc = primary_pr_auc,
      primary_screen_prevalence = primary_prevalence,
      primary_no_skill_pr_auc = no_skill_pr_auc,
      primary_binary_metric = info_value(
        "primary_binary_metric",
        if (identical(foundation_model, "TITAN"))
          "balanced accuracy under empirical-training-prior PLS-LDA rule" else NA
      ),
      binary_crossing_rule = info_value(
        "binary_crossing_rule",
        if (identical(foundation_model, "TITAN"))
          "balanced accuracy >= 0.60 in the supporting TITAN screen" else NA
      ),
      primary_binary_decision_rule = info_value("primary_binary_decision_rule"),
      decision_rule_sensitivity_status = info_value("decision_rule_sensitivity_status"),
      equal_prior_balanced_accuracy = info$equal_prior_balanced_accuracy,
      optimized_balanced_accuracy = info$optimized_balanced_accuracy,
      equal_prior_crossing = info$equal_prior_crossing,
      optimized_crossing = info$optimized_crossing,
      median_equal_component = info$median_equal_component,
      median_optimized_component = info$median_optimized_component,
      median_optimized_threshold = info$median_optimized_threshold,
      primary_estimate_label = info$primary_estimate_label,
      repeated_estimate_label = info$repeated_estimate_label,
      repeated_minus_primary_primary_metric =
        info$repeated_minus_primary_primary_metric,
      repeated_q2 = info$repeated_q2,
      repeated_rmse = info$repeated_rmse,
      repeated_spearman = info$repeated_spearman,
      repeated_sensitivity = info$repeated_sensitivity,
      repeated_specificity = info$repeated_specificity,
      repeated_balanced_accuracy = info$repeated_balanced_accuracy,
      repeated_auc = info$repeated_auc,
      repeated_pr_auc = info$binary_pr_auc,
      ppv_tcga_prevalence = info$binary_ppv_tcga_prevalence,
      npv_tcga_prevalence = info$binary_npv_tcga_prevalence,
      binary_minimum_outer_test_positive =
        info$binary_minimum_outer_test_positive,
      binary_minimum_outer_test_negative =
        info$binary_minimum_outer_test_negative,
      binary_minimum_inner_training_positive =
        info$binary_minimum_inner_training_positive,
      binary_minimum_inner_training_negative =
        info$binary_minimum_inner_training_negative,
      binary_selected_components_median =
        info$binary_selected_components_median,
      binary_selected_components_q1 = info$binary_selected_components_q1,
      binary_selected_components_q3 = info$binary_selected_components_q3,
      binary_selected_components_minimum =
        info$binary_selected_components_minimum,
      binary_selected_components_maximum =
        info$binary_selected_components_maximum,
      binary_repeat_score_spearman = info$binary_repeat_score_spearman,
      binary_repeat_class_agreement = info$binary_repeat_class_agreement,
      binary_all_repeat_class_agreement =
        info$binary_all_repeat_class_agreement,
      continuous_evidence_category = info$continuous_evidence_category,
      continuous_repeat_q2_sd = info$continuous_repeat_q2_sd,
      continuous_prediction_repeat_spearman =
        info$continuous_prediction_repeat_spearman,
      continuous_selected_components_median =
        info$continuous_selected_components_median,
      continuous_selected_components_min = info$continuous_selected_components_min,
      continuous_selected_components_max = info$continuous_selected_components_max,
      continuous_outer_fits_at_ceiling = info$continuous_outer_fits_at_ceiling,
      continuous_outer_fits = info$continuous_outer_fits,
      continuous_outer_fit_ceiling_fraction =
        info$continuous_outer_fit_ceiling_fraction,
      site_grouped_metric_name = info$site_grouped_metric_name,
      site_grouped_metric = info$site_grouped_metric,
      site_performance_delta = info$site_performance_delta,
      site_grouped_n_sites = info$site_grouped_n_sites,
      site_retained_effect = info$site_retained_effect,
      site_near_chance_or_worse = info$site_near_chance_or_worse,
      site_robustness_status = info$site_robustness_status,
      site_robustness_warning = ifelse(
        is.na(info$site_robustness_warning), "", info$site_robustness_warning
      ),
      site_grouped_validation_scope = info$site_grouped_validation_scope,
      site_grouped_outer_folds = info$site_grouped_outer_folds,
      site_grouped_minimum_test_patients = info$site_grouped_minimum_test_patients,
      site_grouped_maximum_test_patients = info$site_grouped_maximum_test_patients,
      site_grouped_minimum_test_sites = info$site_grouped_minimum_test_sites,
      site_grouped_maximum_test_sites = info$site_grouped_maximum_test_sites,
      site_grouped_minimum_test_positive = info$site_grouped_minimum_test_positive,
      site_grouped_maximum_test_positive = info$site_grouped_maximum_test_positive,
      site_grouped_inner_site_separation = info$site_grouped_inner_site_separation,
      matched_tss_n_codes = info$matched_tss_n_codes,
      matched_tss_realized_outer_folds = info$matched_tss_realized_outer_folds,
      matched_tss_outer_test_n_min = info$matched_tss_outer_test_n_min,
      matched_tss_outer_test_n_max = info$matched_tss_outer_test_n_max,
      matched_tss_outer_training_positive_min =
        info$matched_tss_outer_training_positive_min,
      matched_tss_outer_training_negative_min =
        info$matched_tss_outer_training_negative_min,
      matched_tss_any_single_class_outer_test_fold =
        info$matched_tss_any_single_class_outer_test_fold,
      matched_tss_inner_folds_realized_min =
        info$matched_tss_inner_folds_realized_min,
      matched_tss_inner_training_positive_min =
        info$matched_tss_inner_training_positive_min,
      matched_tss_inner_training_negative_min =
        info$matched_tss_inner_training_negative_min,
      matched_tss_any_single_class_inner_training_fold =
        info$matched_tss_any_single_class_inner_training_fold,
      matched_tss_grouped_fold_adequacy_flag =
        info$matched_tss_grouped_fold_adequacy_flag,
      matched_tss_grouped_fold_adequacy_reasons =
        info$matched_tss_grouped_fold_adequacy_reasons,
      matched_tss_metric_construction = info$matched_tss_metric_construction,
      prediction_warning = paste(warnings, collapse = " "),
      ood_fraction = ood, external_validation = info$external_validation,
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  ans <- do.call(rbind, output)
  rownames(ans) <- NULL
  class(ans) <- c("pathofm_predictions", class(ans))
  ans
}

.read_pathofmpred_object <- function(object) {
  if (is.character(object) && length(object) == 1L && file.exists(object)) {
    object <- readRDS(object)
  }
  validate_pathofmpred_object(object)
  object
}

#' Apply a user-created PathoFMPred object
#'
#' This function applies an object returned by [create_pathofmpred_object()].
#' It is the supported inference route for locally constructed TITAN objects
#' and for future foundation-model representations that are not distributed by
#' the package. Repeated feature rows are pooled using the aggregation rule
#' recorded when the object was created, unless `aggregation` is supplied.
#'
#' @param object A `PathoFMPredObject` or the path to an `.rds` file containing
#'   one.
#' @param features A data frame or matrix containing the named features used to
#'   fit `object`. An identifier column may be selected with `id_column`.
#' @param patient_id Optional identifier per feature row. Use either
#'   `patient_id` or `id_column`, not both. Repeated identifiers are pooled.
#' @param id_column Optional identifier-column name in `features`.
#' @param outcome_type One or both of `continuous` and `binary`.
#' @param endpoints Optional endpoint names to apply.
#' @param aggregation Optional `mean` or `median` override. By default the
#'   object's model-development aggregation is used.
#' @param warn_ood Warn when more than 5 percent of dimensions are outside the
#'   fitted training range.
#' @return A long data frame with one row per patient and fitted endpoint.
#'   Binary outputs include the original event-class label. Scores are
#'   uncalibrated and are not probabilities.
#' @export
predict_pathofmpred_object <- function(
    object, features, patient_id = NULL, id_column = NULL,
    outcome_type = c("continuous", "binary"), endpoints = NULL,
    aggregation = NULL, warn_ood = TRUE) {
  object <- .read_pathofmpred_object(object)
  x <- as.data.frame(features, check.names = FALSE)
  if (!nrow(x)) stop("features must contain at least one row.", call. = FALSE)
  if (!is.null(id_column)) {
    if (!is.null(patient_id)) {
      stop("Use either patient_id or id_column, not both.", call. = FALSE)
    }
    if (length(id_column) != 1L || !is.character(id_column) ||
        !id_column %in% names(x)) {
      stop("id_column must name one column in features.", call. = FALSE)
    }
    patient_id <- as.character(x[[id_column]])
  }
  if (is.null(patient_id)) {
    patient_id <- rownames(x)
    if (is.null(patient_id) || any(!nzchar(patient_id))) {
      patient_id <- paste0("sample_", seq_len(nrow(x)))
    }
  }
  patient_id <- as.character(patient_id)
  if (length(patient_id) != nrow(x) || anyNA(patient_id) ||
      any(!nzchar(patient_id))) {
    stop("patient_id must contain one non-missing identifier per feature row.",
         call. = FALSE)
  }
  missing_features <- setdiff(object$feature_names, names(x))
  if (length(missing_features)) {
    stop("Missing feature columns: ",
         paste(utils::head(missing_features, 10L), collapse = ", "),
         if (length(missing_features) > 10L) " ..." else "", call. = FALSE)
  }
  X <- as.matrix(x[, object$feature_names, drop = FALSE])
  storage.mode(X) <- "double"
  if (any(!is.finite(X))) {
    stop("All foundation-model feature values must be finite numeric values.",
         call. = FALSE)
  }
  if (is.null(aggregation)) aggregation <- object$data_audit$aggregation %||% "mean"
  aggregation <- match.arg(aggregation, c("mean", "median"))
  pooled <- .pathofm_pool_rows(X, patient_id, aggregation)
  pooled_ids <- rownames(pooled)
  n_rows <- as.integer(table(factor(patient_id, levels = pooled_ids)))

  outcome_type <- match.arg(outcome_type, c("continuous", "binary"),
                            several.ok = TRUE)
  registry <- object$registry
  registry <- registry[registry$outcome_type %in% outcome_type, , drop = FALSE]
  if (!is.null(endpoints)) {
    registry <- registry[registry$endpoint %in% endpoints, , drop = FALSE]
    absent <- setdiff(endpoints, object$registry$endpoint)
    if (length(absent)) {
      stop("Endpoints are not present in object: ",
           paste(absent, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (!nrow(registry)) stop("No fitted endpoint matched the requested filters.",
                            call. = FALSE)

  output <- vector("list", nrow(registry))
  for (i in seq_len(nrow(registry))) {
    info <- registry[i, , drop = FALSE]
    artifact <- object$models[[info$model_id[[1L]]]]
    .validate_runtime(artifact)
    ood <- rep(NA_real_, nrow(pooled))
    if (!is.null(artifact$training_feature_min) &&
        !is.null(artifact$training_feature_max)) {
      ood <- rowMeans(
        sweep(pooled, 2L, artifact$training_feature_min, `<`) |
          sweep(pooled, 2L, artifact$training_feature_max, `>`)
      )
      if (isTRUE(warn_ood) && any(ood > 0.05)) {
        warning(artifact$model_id, ": ", sum(ood > 0.05),
                " patient(s) have >5% of features outside the training range.",
                call. = FALSE)
      }
    }
    binary <- identical(artifact$outcome_type, "binary")
    pred <- stats::predict(artifact$model, pooled, raw_scores = binary)
    if (binary) {
      score <- drop(pred$LDA_scores[, 2L, 1L] - pred$LDA_scores[, 1L, 1L])
      threshold <- artifact$operating_threshold
      if (is.null(threshold) || !length(threshold) || !is.finite(threshold)) {
        stop("Binary artifact ", artifact$model_id,
             " has no finite training-derived operating threshold.", call. = FALSE)
      }
      class01 <- as.character(as.integer(score >= threshold))
      labels <- artifact$class_labels
      class_label <- if (!is.null(labels) &&
                          all(c("negative", "positive") %in% names(labels))) {
        ifelse(class01 == "1", labels$positive, labels$negative)
      } else class01
      value <- rep(NA_real_, length(score))
    } else {
      value <- if (length(dim(pred$Ypred)) == 3L) {
        drop(pred$Ypred[, 1L, 1L])
      } else drop(pred$Ypred)
      score <- rep(NA_real_, length(value))
      threshold <- NA_real_
      class01 <- rep(NA_character_, length(value))
      class_label <- rep(NA_character_, length(value))
    }
    output[[i]] <- data.frame(
      patient_id = pooled_ids,
      foundation_model = object$foundation_model,
      n_feature_rows = n_rows,
      aggregation = aggregation,
      model_id = artifact$model_id,
      endpoint = artifact$endpoint,
      outcome_type = artifact$outcome_type,
      prediction = value,
      lda_score = score,
      predicted_class = class01,
      predicted_class_label = class_label,
      operating_threshold = threshold,
      score_interpretation = if (binary) {
        "Uncalibrated discriminant score, not probability"
      } else {
        "Prediction in the outcome scale used for model fitting"
      },
      training_n = artifact$training_n,
      training_positive = artifact$training_positive,
      training_negative = artifact$training_negative,
      ncomp = artifact$ncomp,
      ood_fraction = ood,
      external_validation = "not established",
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  answer <- do.call(rbind, output)
  rownames(answer) <- NULL
  class(answer) <- c("pathofm_object_predictions", class(answer))
  answer
}

#' Apply the legacy TITAN model interface
#'
#' Backward-compatible wrapper for [predict_pathofm()] with
#' `foundation_model = "TITAN"`.
#' @inheritParams predict_pathofm
#' @export
predict_titan <- function(cancer, features, patient_id = NULL,
                          outcome_type = c("continuous", "binary"),
                          family = NULL, warn_ood = TRUE,
                          include_limited_evidence = FALSE) {
  predict_pathofm(
    cancer = cancer, features = features, foundation_model = "TITAN",
    patient_id = patient_id, outcome_type = outcome_type, family = family,
    warn_ood = warn_ood,
    include_limited_evidence = include_limited_evidence
  )
}

#' Construct a synthetic feature vector for examples and smoke tests
#'
#' The vector is the across-model mean of stored training-feature means for the
#' requested cancer and is not a patient record.
#'
#' @param cancer A cancer code with bundled models.
#' @param include_limited_evidence Use limited-evidence models when constructing
#'   the synthetic across-model feature mean. The default is `FALSE`.
#' @param foundation_model One of `TITAN`, `GigaSSL`, or `ProvGigaPath`.
#' @return A one-row data frame containing the named features required by the
#'   selected representation.
#' @export
example_pathofm_features <- function(cancer = "BRCA", foundation_model = "TITAN",
                                     include_limited_evidence = FALSE) {
  models <- available_pathofm_models(
    foundation_model, cancer,
    include_limited_evidence = include_limited_evidence
  )
  if (!nrow(models)) stop("No bundled models are available for ", cancer, ".", call. = FALSE)
  means <- lapply(models$model_id, function(id) {
    .pathofm_read_model(id, foundation_model)$training_feature_mean
  })
  value <- Reduce(`+`, means) / length(means)
  out <- as.data.frame(as.list(value), check.names = FALSE)
  rownames(out) <- paste0("synthetic_", toupper(cancer))
  out
}

#' @rdname example_pathofm_features
#' @export
example_titan_features <- function(cancer = "BRCA",
                                   include_limited_evidence = FALSE) {
  example_pathofm_features(cancer, "TITAN", include_limited_evidence)
}
