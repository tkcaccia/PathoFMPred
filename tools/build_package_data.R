args <- commandArgs(trailingOnly = TRUE)
source_root <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else
  normalizePath("../titan-prediction", mustWork = TRUE)
package_root <- if (length(args) >= 2L) normalizePath(args[2], mustWork = TRUE) else
  normalizePath(getwd(), mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

registry <- fread(file.path(source_root, "models", "model_registry.csv"))
registry[, foundation_model := "TITAN"]
for (field in c(
  "resource_selection_basis", "representation_specific_qualification",
  "titan_catalogue_tier", "tier_origin",
  "representation_effect_threshold_crossing"
)) {
  if (!field %in% names(registry)) registry[, (field) := NA]
}
required_site_fields <- c(
  "site_grouped_metric_name", "site_grouped_metric", "site_performance_delta",
  "site_grouped_n_sites", "site_retained_effect",
  "site_near_chance_or_worse", "site_robustness_status",
  "site_robustness_warning", "site_grouped_validation_scope",
  "site_grouped_inner_site_separation"
)
if (!all(required_site_fields %in% names(registry))) {
  stop("Analysis model registry is missing required site-robustness fields")
}
required_reliability_fields <- c(
  "model_evidence_tier", "default_inference", "limited_evidence_reason",
  "binary_pr_auc", "binary_ppv_tcga_prevalence",
  "binary_npv_tcga_prevalence", "binary_minimum_outer_test_positive",
  "binary_minimum_inner_training_positive",
  "binary_selected_components_median", "binary_repeat_score_spearman",
  "binary_repeat_class_agreement"
)
if (!all(required_reliability_fields %in% names(registry))) {
  stop("Analysis model registry is missing binary reliability fields")
}
required_decision_rule_fields <- c(
  "primary_binary_decision_rule", "decision_rule_sensitivity_status",
  "equal_prior_balanced_accuracy", "optimized_balanced_accuracy",
  "equal_prior_crossing", "optimized_crossing",
  "median_equal_component", "median_optimized_component",
  "median_optimized_threshold"
)
if (!all(required_decision_rule_fields %in% names(registry))) {
  stop("Analysis model registry is missing binary decision-rule sensitivity fields")
}
continuous <- fread(file.path(source_root, "results", "tables", "continuous_screen.csv"))
binary <- fread(file.path(source_root, "results", "tables", "binary_screen.csv"))

continuous_metrics <- continuous[, .(
  family, cancer_type = tumor_type, endpoint, screen_q2 = q2,
  screen_rmse = rmse, screen_spearman = spearman,
  screen_q_value = q_value, screen_q_value_global = q_value_global
)]
binary_metrics <- binary[, .(
  family, cancer_type = tumor_type, endpoint,
  screen_balanced_accuracy = balanced_accuracy,
  screen_q_value = q_value, screen_q_value_global = q_value_global
)]

registry <- merge(registry, continuous_metrics,
                  by = c("family", "cancer_type", "endpoint"), all.x = TRUE)
registry <- merge(registry, binary_metrics,
                  by = c("family", "cancer_type", "endpoint"), all.x = TRUE,
                  suffixes = c("", "_binary"))
registry[, screen_q_value := fcoalesce(screen_q_value, screen_q_value_binary)]
registry[, screen_q_value_global := fcoalesce(screen_q_value_global, screen_q_value_global_binary)]
registry[, c("screen_q_value_binary", "screen_q_value_global_binary") := NULL]

# Repeated nested-CV performance is research-use evidence metadata. It lets outputs show
# the evidence supporting each call next to the training denominator, without
# implying that an uncalibrated LDA score is a probability.
performance <- fread(file.path(source_root, "results", "tables",
                               "screen_positive_performance_summary.csv"))
performance <- performance[, .(
  family, cancer_type = tumor_type, endpoint,
  repeated_q2 = repeated_q2_mean,
  repeated_rmse = repeated_rmse_mean,
  repeated_spearman = repeated_spearman_mean,
  repeated_sensitivity = repeated_sensitivity_mean,
  repeated_specificity = repeated_specificity_mean,
  repeated_balanced_accuracy = repeated_balanced_accuracy_mean,
  repeated_auc = repeated_auc_mean,
  nested_partitions,
  primary_estimate_label,
  repeated_estimate_label,
  repeated_minus_primary_primary_metric = fifelse(
    outcome_type == "continuous",
    repeated_q2_mean - primary_screen_q2,
    repeated_balanced_accuracy_mean - primary_screen_balanced_accuracy
  )
)]
registry <- merge(registry, performance,
                  by = c("family", "cancer_type", "endpoint"), all.x = TRUE)

# Attach the continuous-outcome reliability audit used by the inference
# report. These fields are produced after fitted-model creation and therefore
# are not present in the base analysis registry unless they are joined here.
continuous_reliability <- fread(file.path(
  source_root, "results", "tables", "continuous_reliability_by_model.csv"
))
continuous_reliability <- continuous_reliability[, .(
  family, cancer_type = tumor_type, endpoint,
  continuous_evidence_category = evidence_category,
  continuous_repeat_q2_sd = repeated_q2_sd,
  continuous_prediction_repeat_spearman = prediction_repeat_spearman_mean,
  continuous_selected_components_median = selected_components_median,
  continuous_selected_components_min = selected_components_min,
  continuous_selected_components_max = selected_components_max,
  continuous_outer_fits_at_ceiling = outer_fits_at_ceiling,
  continuous_outer_fits = outer_fits,
  continuous_outer_fit_ceiling_fraction = outer_fit_ceiling_fraction
)]
registry <- merge(
  registry, continuous_reliability,
  by = c("family", "cancer_type", "endpoint"), all.x = TRUE
)
registry[, file := basename(file)]
registry[, redistribution_status := paste(
  "private TITAN-derived artifact; public redistribution prohibited unless",
  "the TITAN rights holder grants written permission"
)]
limited <- registry$model_evidence_tier %in% c(
  "exploratory_limited_evidence",
  "exploratory_limited_continuous_evidence"
)
if (anyNA(registry$default_inference) || !any(registry$default_inference) ||
    any(registry$default_inference[limited]) ||
    any(registry$default_inference[
      registry$model_evidence_tier == "matched_tested_below_effect_threshold"
    ])) {
  stop("Invalid default/limited-evidence model status")
}
setorder(registry, cancer_type, outcome_type, family, endpoint)
registry[, source_file := file]
registry[, package_file := paste0("m_", substr(sha256, 1L, 24L), ".rds")]
registry[, file := package_file]
registry[, c("source_file", "package_file") := NULL]

continuous_oof <- fread(file.path(source_root, "results", "predictions",
                                  "continuous_repeated_oof_predictions.csv.gz"))
continuous_oof <- continuous_oof[, .(reference_value = mean(predicted)),
                                 by = .(family, cancer_type = tumor_type, endpoint, patient)]
binary_oof <- fread(file.path(source_root, "results", "predictions",
                              "binary_repeated_oof_predictions.csv.gz"))
binary_oof <- binary_oof[, .(reference_value = mean(lda_score)),
                         by = .(family, cancer_type = tumor_type, endpoint, patient)]
reference <- rbind(continuous_oof, binary_oof, use.names = TRUE)
reference <- merge(reference, registry[, .(model_id, family, cancer_type, endpoint)],
                   by = c("family", "cancer_type", "endpoint"), all.x = FALSE)
reference <- reference[is.finite(reference_value)]
reference_list <- lapply(split(reference$reference_value,
                               paste("TITAN", reference$model_id, sep = "::")), sort)

addition_registry_path <- file.path(
  source_root, "models", "foundation_models", "model_registry_additions.csv"
)
if (file.exists(addition_registry_path)) {
  additions <- fread(addition_registry_path)
  additions[, source_file := file]
  additions[, package_file := paste0("m_", substr(sha256, 1L, 24L), ".rds")]
  registry <- rbindlist(list(registry, additions), use.names = TRUE, fill = TRUE)
  registry[foundation_model != "TITAN", file := package_file]
  registry[, c("source_file", "package_file") := NULL]
  addition_reference <- readRDS(file.path(
    source_root, "models", "foundation_models", "prediction_reference_additions.rds"
  ))
  reference_list[names(addition_reference)] <- addition_reference
}

setorder(registry, foundation_model, cancer_type, outcome_type, family, endpoint)

# This is the public-package builder. Full fitted collections must never be
# copied into the package tree. Giga-SSL and Prov-GigaPath collections are
# separate, user-invoked downloads; TITAN objects remain private. Fail closed
# if a maintainer accidentally leaves any full object in a representation
# subdirectory. The tiny root-level Giga-SSL fixture is intentionally retained.
for (model in c("TITAN", "GigaSSL", "ProvGigaPath")) {
  model_dir <- file.path(package_root, "inst", "models", model)
  existing_files <- list.files(
    model_dir, pattern = "[.]rds$", full.names = TRUE
  )
  if (length(existing_files)) {
    stop(
      "Public package tree contains fitted ", model,
      " object(s). Remove them and use the post-install download mechanism."
    )
  }
}

binary_rows <- registry$outcome_type == "binary"
if (anyNA(registry$primary_binary_decision_rule[binary_rows]) ||
    anyNA(registry$decision_rule_sensitivity_status[binary_rows]) ||
    anyNA(registry$equal_prior_balanced_accuracy[binary_rows]) ||
    anyNA(registry$optimized_balanced_accuracy[binary_rows])) {
  stop("One or more binary package models lack operating-rule sensitivity metadata")
}
alternative_rows <- registry$foundation_model != "TITAN"
if (anyNA(registry$representation_effect_threshold_crossing[alternative_rows]) ||
    anyNA(registry$tier_origin[alternative_rows]) ||
    anyNA(registry$resource_selection_basis[alternative_rows]) ||
    anyNA(registry$representation_specific_qualification[alternative_rows])) {
  stop("One or more alternative-representation models lack selection/qualification metadata")
}
fwrite(registry, file.path(package_root, "inst", "extdata", "model_registry.csv"))
saveRDS(reference_list, file.path(package_root, "inst", "extdata", "prediction_reference.rds"),
        compress = "xz")

cat("Registry rows:", nrow(registry), "\n")
cat("Reference models:", length(reference_list), "\n")
cat("Reference values:", nrow(reference), "\n")
