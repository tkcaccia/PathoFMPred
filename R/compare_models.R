#' Compare representation-specific models for one cancer-endpoint task
#'
#' Lists every fitted reference object available for the same task and
#' distinguishes object availability from a representation-specific matched-
#' atlas effect-threshold crossing. The fitted-object inventory is only a partial
#' operationalization of the multi-representation atlas. The function does not
#' automatically select the largest internal performance estimate or combine
#' predictions: representation choice must be fixed before external evaluation,
#' and the representations are not interchangeable.
#'
#' @param cancer TCGA cancer code.
#' @param endpoint Exact endpoint label.
#' @param outcome_type Optional `continuous` or `binary` filter.
#' @param family Optional endpoint-family filter, useful when an endpoint label
#'   occurs in more than one family.
#' @param include_limited_evidence Include limited-sample models only when they
#'   crossed the representation-specific matched effect threshold.
#' @return A registry data frame with object-availability, matched-threshold,
#'   resource-scope and research-use selection fields repeated on every row.
#' @export
compare_pathofm_models <- function(cancer, endpoint, outcome_type = NULL,
                                   family = NULL,
                                   include_limited_evidence = FALSE) {
  all_models <- pathofm_model_registry()
  all_models <- all_models[
    all_models$cancer_type %in% toupper(cancer), , drop = FALSE
  ]
  if (!is.null(outcome_type)) {
    all_models <- all_models[
      all_models$outcome_type %in% outcome_type, , drop = FALSE
    ]
  }
  if (!is.null(family)) {
    all_models <- all_models[all_models$family %in% family, , drop = FALSE]
  }
  all_models <- all_models[all_models$endpoint %in% endpoint, , drop = FALSE]
  if (!nrow(all_models)) {
    stop(paste(
      "No fitted reference object matches this cancer-endpoint task.",
      "This does not imply that the task was absent or negative in the matched",
      "atlas; the fitted-object inventory does not cover every representation-specific crossing."
    ),
         call. = FALSE)
  }

  models <- all_models[.pathofm_predictable_mask(all_models), , drop = FALSE]
  if (!isTRUE(include_limited_evidence)) {
    keep <- !grepl("^exploratory_limited", models$model_evidence_tier)
    if (!any(keep)) {
      stop(paste(
        "Fitted object(s) exist for this task, but no larger-sample",
        "representation-specific model crossed its matched effect threshold.",
        "Re-run with include_limited_evidence = TRUE only to include a",
        "limited-sample model that still crossed the threshold."
      ), call. = FALSE)
    }
    models <- models[keep, , drop = FALSE]
  }
  if (!nrow(models)) {
    stop(
      paste(
        "Fitted object(s) exist for this task, but none crossed the",
        "representation-specific matched effect threshold."
      ),
      call. = FALSE
    )
  }

  crossing <- models$matched_representation_effect_threshold_crossing
  models$matched_threshold_status <- ifelse(
    is.na(crossing),
    "not in the matched common-cohort analysis",
    ifelse(crossing,
           "crossed the representation-specific matched effect threshold",
           "tested below the representation-specific matched effect threshold")
  )
  models$controlled_object_available <- TRUE
  models$object_availability_class <- paste0(
    length(unique(all_models$foundation_model)),
    " representation-specific fitted object(s) exist for this task; ",
    "object count is not a consensus count"
  )
  models$resource_completeness <- paste(
    "Partial operationalization: alternative-representation objects were fitted",
    "for the shared TITAN-qualified matched-eligible target universe, not for",
    "every Giga-SSL or Prov-GigaPath effect-threshold crossing."
  )
  models$selection_guidance <- paste(
    "Do not interpret object presence as threshold crossing, consensus or",
    "representation-specific permutation/FDR qualification.",
    "Do not choose the numerically largest internal TCGA estimate or average",
    "representation-specific predictions post hoc. Use the model matching the",
    "available embedding; for a new comparison, lock the representation before",
    "outcome inspection. Prioritize cross-representation retention, stronger",
    "sample-size stratum and retained sensitivity to grouping by the barcode-derived",
    "TCGA tissue-source-site code where available. Inspect the grouped-fold",
    "adequacy flag separately; grouped retention with sparse folds is not evidence",
    "of generalisability. External validation remains required."
  )
  models
}
