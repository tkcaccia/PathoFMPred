test_that("installed model collections and registry agree", {
  installed <- available_foundation_models()
  expect_true("GigaSSL" %in% installed)
  registry <- pathofm_model_registry()
  expect_gt(nrow(registry), 0L)
  for (model in installed) {
    available <- available_pathofm_models(
      model, include_limited_evidence = TRUE
    )
    expect_gt(nrow(available), 0L)
    expect_true(all(available$foundation_model == model))
  }
})

test_that("minimal GigaSSL example predicts one binary and one continuous target", {
  models <- available_pathofm_models(
    "GigaSSL", "COAD", include_limited_evidence = TRUE
  )
  expect_true(all(c("APC", "MANTIS score") %in% models$endpoint))
  x <- example_pathofm_features(
    "COAD", "GigaSSL", include_limited_evidence = TRUE
  )
  answer <- suppressWarnings(predict_pathofm(
    "COAD", x, "GigaSSL", patient_id = "software_fixture",
    include_limited_evidence = TRUE
  ))
  expect_gte(nrow(answer), 2L)
  expect_true(all(c("binary", "continuous") %in% answer$outcome_type))
  binary <- answer[answer$endpoint == "APC" &
                     answer$outcome_type == "binary", , drop = FALSE]
  expect_equal(nrow(binary), 1L)
  expect_true(binary$predicted_class %in% c("0", "1"))
  expect_true(is.finite(binary$lda_score))
  expect_false(binary$rank_is_probability)
  expect_match(binary$rank_interpretation, "not probability", fixed = TRUE)
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
})

test_that("public fetch excludes TITAN and validates downloaded objects", {
  expect_error(fetch_pathofmpred_models("TITAN"), "arg")
  source <- system.file(
    "models", "pathofmpred_gigassl_example.rds", package = "PathoFMPred"
  )
  mirror <- file.path(tempdir(), "pathofmpred-mirror")
  destination <- file.path(tempdir(), "pathofmpred-download")
  dir.create(mirror, showWarnings = FALSE)
  file.copy(source, file.path(mirror, "pathofmpred_gigassl.rds"),
            overwrite = TRUE)
  path <- fetch_pathofmpred_models(
    "GigaSSL", destination = destination,
    base_url = paste0("file://", mirror), overwrite = TRUE, quiet = TRUE,
    verify = FALSE
  )
  expect_true(file.exists(path[["GigaSSL"]]))
  expect_true(validate_pathofmpred_object(readRDS(path[["GigaSSL"]])))
})
