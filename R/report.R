#' Create a single-sample pathology-foundation-model research output
#'
#' @param cancer A single TCGA cancer code.
#' @param features One or more slide-feature rows for one patient.
#' @param foundation_model One of `TITAN`, `GigaSSL`, or `ProvGigaPath`.
#' @param patient_id Optional patient/sample identifier for all rows.
#' @param output_file Explicit file stem or `.html`/`.pdf` path. Required;
#'   no report is written unless the caller chooses a destination.
#' @param format One of `html`, `pdf`, or `both`.
#' @param clinical_context Optional named list with `shared_features`,
#'   `pathology`, and `source` character fields. This information
#'   is displayed as supplied and is not used by any prediction model.
#' @param include_limited_evidence Explicitly opt in to limited-sample models
#'   that still crossed the matched effect threshold for the selected cancer
#'   and representation. Tested-below-threshold objects never enter reports.
#' @param quiet Passed to [rmarkdown::render()].
#' @return Invisibly, the generated file path(s).
#' @export
pathofm_sample_report <- function(cancer, features,
                                  foundation_model = c("TITAN", "GigaSSL", "ProvGigaPath"),
                                  patient_id = NULL,
                                  output_file = NULL,
                                  format = c("html", "pdf", "both"),
                                  clinical_context = NULL,
                                  include_limited_evidence = FALSE,
                                  quiet = TRUE) {
  foundation_model <- match.arg(foundation_model)
  format <- match.arg(format)
  if (is.null(output_file) || length(output_file) != 1L ||
      !is.character(output_file) || is.na(output_file) ||
      !nzchar(output_file)) {
    stop("output_file must be an explicit path selected by the caller.",
         call. = FALSE)
  }
  if (length(unique(toupper(cancer))) != 1L) stop("An output must contain one cancer type.", call. = FALSE)
  if (is.null(patient_id)) patient_id <- rep("sample", nrow(as.data.frame(features)))
  if (length(unique(patient_id)) != 1L) stop("An output must contain one patient_id.", call. = FALSE)
  if (!is.null(clinical_context)) {
    if (!is.list(clinical_context)) stop("clinical_context must be NULL or a named list.", call. = FALSE)
    allowed <- c("shared_features", "pathology", "source")
    if (is.null(names(clinical_context)) || any(!names(clinical_context) %in% allowed)) {
      stop("clinical_context names must be drawn from: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
    clinical_context <- lapply(clinical_context, as.character)
  }
  pred <- predict_pathofm(
    cancer, features, foundation_model = foundation_model,
    patient_id = patient_id,
    include_limited_evidence = include_limited_evidence
  )
  template <- system.file("rmarkdown", "templates", "pathofm-report", "skeleton", "skeleton.Rmd", package = "PathoFMPred")
  if (!nzchar(template)) stop("Report template is unavailable.", call. = FALSE)
  render_dir <- tempfile("pathofm-report-", tmpdir = tempdir())
  dir.create(file.path(render_dir, "skeleton"), recursive = TRUE)
  dir.create(file.path(render_dir, "resources"))
  on.exit(unlink(render_dir, recursive = TRUE), add = TRUE)
  render_input <- file.path(render_dir, "skeleton", "report.Rmd")
  report_css <- file.path(dirname(dirname(template)), "resources", "report.css")
  if (!file.copy(template, render_input) ||
      !file.copy(report_css, file.path(render_dir, "resources", "report.css"))) {
    stop("Could not prepare the report template in a temporary directory.",
         call. = FALSE)
  }
  stem <- sub("\\.(html|pdf)$", "", output_file, ignore.case = TRUE)
  formats <- if (format == "both") c("html", "pdf") else format
  outputs <- character(length(formats))
  for (i in seq_along(formats)) {
    extension <- formats[i]
    out <- paste0(basename(stem), ".", extension)
    outputs[i] <- rmarkdown::render(
      render_input, output_format = if (extension == "html") "html_document" else "pdf_document",
      output_file = out, output_dir = dirname(stem),
      params = list(predictions = pred, sample_id = unique(patient_id),
                    cancer = unique(toupper(cancer)),
                    clinical_context = clinical_context),
      envir = new.env(parent = globalenv()), quiet = quiet
    )
  }
  invisible(normalizePath(outputs, mustWork = FALSE))
}

#' @rdname pathofm_sample_report
#' @export
titan_sample_report <- function(cancer, features, patient_id = NULL,
                                output_file = NULL,
                                format = c("html", "pdf", "both"),
                                clinical_context = NULL,
                                include_limited_evidence = FALSE,
                                quiet = TRUE) {
  pathofm_sample_report(
    cancer, features, foundation_model = "TITAN", patient_id = patient_id,
    output_file = output_file, format = format,
    clinical_context = clinical_context,
    include_limited_evidence = include_limited_evidence, quiet = quiet
  )
}
