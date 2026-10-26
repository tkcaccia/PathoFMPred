test_that("the package starts without fitted models or an atlas registry", {
  old <- getOption("PathoFMPred.model_dir")
  on.exit(options(PathoFMPred.model_dir = old), add = TRUE)
  options(PathoFMPred.model_dir = tempfile("pathofmpred-empty-"))
  installed <- available_foundation_models()
  expect_length(installed, 0L)
  registry <- pathofm_model_registry()
  expect_equal(nrow(registry), 0L)
  expect_equal(nrow(available_pathofm_models("GigaSSL")), 0L)
})

test_that("atlas examples require an authorized locally installed collection", {
  old <- getOption("PathoFMPred.model_dir")
  on.exit(options(PathoFMPred.model_dir = old), add = TRUE)
  options(PathoFMPred.model_dir = tempfile("pathofmpred-empty-"))
  expect_error(example_pathofm_features("COAD", "GigaSSL"),
               "No predictable")
})

test_that("runtime validation requires the recorded fastPLS version", {
  validator <- getFromNamespace(".validate_runtime", "PathoFMPred")
  expect_error(
    validator(list(model_id = "runtime_fixture", fastPLS_version = "0.0.0")),
    "requires fastPLS 0.0.0"
  )
  expect_silent(validator(list(
    model_id = "runtime_fixture",
    fastPLS_version = as.character(utils::packageVersion("fastPLS")),
    fastPLS_repository = "another-installation-source",
    fastPLS_remote_sha = "another-revision"
  )))
})

test_that("builder enforces ID rules, pools repeats, and reports missingness", {
  set.seed(14)
  ids <- sprintf("patient_%03d", seq_len(112L))
  features <- data.frame(
    patient_id = rep(c(ids, "feature_only"), each = 2L),
    feature_1 = stats::rnorm(226L),
    feature_2 = stats::rnorm(226L),
    check.names = FALSE
  )
  signal <- tapply(features$feature_1, features$patient_id, mean)[ids]
  outcomes <- data.frame(
    patient_id = c(ids, "outcome_only"),
    continuous = c(as.numeric(signal) + stats::rnorm(112L), NA_real_),
    binary = c(as.integer(signal > 0), NA_integer_),
    check.names = FALSE
  )
  settings <- pathofmpred_control(
    components = 1:2, outer_folds = 3, inner_folds = 3,
    repeats = 1, method = "simpls", seed = 5
  )
  object <- NULL
  expect_message(
    object <- create_pathofmpred_object(
      features, outcomes, "patient_id", control = settings
    ),
    "112 IDs matched"
  )
  expect_s3_class(object, "PathoFMPredObject")
  expect_length(object$models, 2L)
  expect_equal(object$data_audit$matching_ids, 112L)
  expect_equal(object$data_audit$unmatched_feature_ids, "feature_only")
  expect_equal(object$data_audit$unmatched_outcome_ids, "outcome_only")
  expect_equal(object$data_audit$input_outcome_rows, 113L)
  expect_equal(object$data_audit$outcome_rows, 112L)
  expect_true(all(object$data_audit$input_outcome_missingness$missing == 1L))
  expect_true(all(object$data_audit$outcome_missingness$missing == 0L))
  expect_true(all(vapply(object$models, function(x) {
    !is.null(x$nested_cv) && inherits(x$model, "fastPLS")
  }, logical(1))))

  duplicate_outcome <- rbind(outcomes, outcomes[1L, , drop = FALSE])
  expect_error(
    create_pathofmpred_object(
      features, duplicate_outcome, "patient_id", control = settings
    ),
    "unique"
  )
  expect_error(
    create_pathofmpred_object(
      features[features$patient_id %in% ids[1:100], ],
      outcomes[outcomes$patient_id %in% ids[1:100], ],
      "patient_id", control = settings
    ),
    "more than 100"
  )
  duplicate_feature_name <- features
  names(duplicate_feature_name)[3L] <- "feature_1"
  expect_error(
    create_pathofmpred_object(
      duplicate_feature_name, outcomes, "patient_id", control = settings
    ),
    "unique column names"
  )
})

test_that("builder requires and preserves the binary event-class meaning", {
  set.seed(29)
  ids <- sprintf("patient_%03d", seq_len(112L))
  features <- data.frame(
    patient_id = ids,
    feature_1 = stats::rnorm(112L),
    feature_2 = stats::rnorm(112L)
  )
  outcomes <- data.frame(
    patient_id = ids,
    status = ifelse(features$feature_1 > 0, "mutated", "wild-type")
  )
  settings <- pathofmpred_control(
    components = 1:2, outer_folds = 3, inner_folds = 3,
    repeats = 1, method = "simpls", seed = 6
  )
  expect_error(
    create_pathofmpred_object(
      features, outcomes, "patient_id", outcome_types = c(status = "binary"),
      control = settings
    ),
    "positive_class"
  )
  object <- suppressMessages(create_pathofmpred_object(
    features, outcomes, "patient_id", outcome_types = c(status = "binary"),
    positive_class = c(status = "mutated"), control = settings
  ))
  artifact <- object$models[[1L]]
  expect_identical(artifact$class_labels$positive, "mutated")
  expect_identical(artifact$class_labels$negative, "wild-type")
  expect_identical(object$registry$positive_label, "mutated")

  new_features <- data.frame(
    patient_id = rep(c("new_a", "new_b"), each = 2L),
    feature_1 = c(-1.0, -0.8, 1.0, 0.8),
    feature_2 = c(0.1, 0.2, -0.1, -0.2)
  )
  prediction <- suppressWarnings(predict_pathofmpred_object(
    object, new_features, id_column = "patient_id"
  ))
  expect_s3_class(prediction, "pathofm_object_predictions")
  expect_equal(nrow(prediction), 2L)
  expect_equal(prediction$n_feature_rows, c(2L, 2L))
  expect_true(all(prediction$predicted_class_label %in%
                    c("mutated", "wild-type")))
  expect_true(all(grepl("not probability", prediction$score_interpretation,
                        fixed = TRUE)))

  vector_label_object <- object
  vector_label_object$models[[1L]]$class_labels <-
    c(negative = "wild-type", positive = "mutated")
  vector_label_prediction <- suppressWarnings(predict_pathofmpred_object(
    vector_label_object, new_features, id_column = "patient_id"
  ))
  expect_equal(vector_label_prediction$predicted_class_label,
               prediction$predicted_class_label)

  path <- tempfile(fileext = ".rds")
  saveRDS(object, path)
  prediction_from_file <- suppressWarnings(predict_pathofmpred_object(
    path, new_features, id_column = "patient_id", endpoints = "status"
  ))
  expect_equal(prediction_from_file$predicted_class_label,
               prediction$predicted_class_label)
})

test_that("user-supplied fetch excludes TITAN and verifies checksums", {
  expect_error(fetch_pathofmpred_models("TITAN"), "arg")
  expect_error(fetch_pathofmpred_models("ProvGigaPath"), "destination")
  mirror <- file.path(tempdir(), "pathofmpred-mirror")
  destination <- file.path(tempdir(), "pathofmpred-download")
  dir.create(mirror, showWarnings = FALSE)
  source <- file.path(mirror, "pathofmpred_gigassl.rds")
  saveRDS(list(
    format_version = "1", foundation_model = "GigaSSL",
    feature_names = c("feature_1", "feature_2"),
    registry = data.frame(
      model_id = "synthetic_fixture", foundation_model = "GigaSSL",
      cancer_type = "COAD", endpoint = "synthetic_endpoint",
      outcome_type = "binary", family = "synthetic",
      matched_representation_effect_threshold_crossing = TRUE
    ),
    models = list(synthetic_fixture = list(model_id = "synthetic_fixture")),
    provenance = list(source = "synthetic unit test")
  ), source)
  sha <- digest::digest(file = source, algo = "sha256", serialize = FALSE)
  expect_error(fetch_pathofmpred_models("GigaSSL"), "base_url")
  expect_error(fetch_pathofmpred_models(
    "GigaSSL", base_url = paste0("file://", mirror)
  ), "expected_sha256")
  path <- fetch_pathofmpred_models(
    "GigaSSL", destination = destination,
    base_url = paste0("file://", mirror),
    expected_sha256 = c(GigaSSL = sha), overwrite = TRUE, quiet = TRUE
  )
  expect_true(file.exists(path[["GigaSSL"]]))
  expect_true(validate_pathofmpred_object(readRDS(path[["GigaSSL"]])))
  old <- getOption("PathoFMPred.model_dir")
  on.exit(options(PathoFMPred.model_dir = old), add = TRUE)
  options(PathoFMPred.model_dir = destination)
  expect_identical(available_foundation_models(), "GigaSSL")
  expect_equal(nrow(available_pathofm_models("GigaSSL", "COAD")), 1L)

  writeBin(as.raw(1:8), path[["GigaSSL"]])
  expect_error(
    fetch_pathofmpred_models(
      "GigaSSL", destination = destination,
      base_url = paste0("file://", mirror), overwrite = FALSE, quiet = TRUE,
      expected_sha256 = c(GigaSSL = sha)
    ),
    "SHA-256 verification failed"
  )
})

test_that("report creation requires an explicit output path", {
  expect_error(pathofm_sample_report("COAD", data.frame(x = 1)),
               "output_file")
  expect_error(titan_sample_report("COAD", data.frame(x = 1)),
               "output_file")
})
