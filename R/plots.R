.short_endpoint <- function(x, width = 32L) {
  vapply(x, function(z) if (nchar(z) <= width) z else paste0(substr(z, 1L, width - 3L), "..."), character(1))
}

.filter_predictable_predictions <- function(predictions) {
  if ("predictable_for_selected_cancer_and_representation" %in% names(predictions)) {
    keep <- predictions$predictable_for_selected_cancer_and_representation %in% TRUE
  } else if ("matched_representation_effect_threshold_crossing" %in% names(predictions)) {
    keep <- predictions$matched_representation_effect_threshold_crossing %in% TRUE
  } else if ("representation_effect_threshold_crossing" %in% names(predictions)) {
    keep <- predictions$representation_effect_threshold_crossing %in% TRUE
  } else {
    stop(
      "Prediction rows lack representation-specific predictability metadata. ",
      "Regenerate them with predict_pathofm().",
      call. = FALSE
    )
  }
  predictions[keep, , drop = FALSE]
}

.pathofm_colour <- function(predictions) {
  model <- if ("foundation_model" %in% names(predictions)) {
    unique(as.character(predictions$foundation_model))
  } else {
    "TITAN"
  }
  colours <- c(TITAN = "#2C7FB8", GigaSSL = "#E07A3F",
               ProvGigaPath = "#2A9D8F")
  unname(colours[if (length(model) == 1L && model %in% names(colours)) model else "TITAN"])
}

#' Radar plot of continuous prediction percentiles
#'
#' @param predictions Output from [predict_titan()] for one patient.
#' @param max_endpoints Maximum number of endpoints, prioritising extremes.
#' @return A ggplot object. Endpoint names and original model predictions are
#'   placed at the outer radar corners to avoid overlap with the profile.
#' @export
plot_titan_radar <- function(predictions, max_endpoints = 12L) {
  predictions <- .filter_predictable_predictions(predictions)
  d <- predictions[predictions$outcome_type == "continuous" &
                     is.finite(predictions$reference_percentile), , drop = FALSE]
  if (!nrow(d)) stop("No continuous predictions with reference percentiles.", call. = FALSE)
  if (length(unique(d$patient_id)) != 1L) stop("Supply predictions for one patient.", call. = FALSE)
  d$extreme <- abs(d$reference_percentile - 50)
  if (nrow(d) > max_endpoints) {
    d <- d[order(d$extreme, decreasing = TRUE), , drop = FALSE]
    d <- utils::head(d, max_endpoints)
  }
  semantic_order <- c(
    "Lymphocyte Infiltration Signature Score", "TIL Regional Fraction",
    "MANTIS score", "MSIsensor score", "Aneuploidy score",
    "Deleted arm count", "Aneuploidy Score", "Silent Mutation Rate",
    "Nonsilent Mutation Rate", "SNV Neoantigens"
  )
  d$semantic_rank <- match(d$endpoint, semantic_order)
  d$semantic_rank[is.na(d$semantic_rank)] <- length(semantic_order) +
    rank(d$endpoint[is.na(d$semantic_rank)], ties.method = "first")
  d <- d[order(d$semantic_rank), , drop = FALSE]
  display_label <- c(
    "Lymphocyte Infiltration Signature Score" = "Lymphocyte signature",
    "TIL Regional Fraction" = "TIL regional fraction",
    "MANTIS score" = "MANTIS score",
    "MSIsensor score" = "MSIsensor score",
    "Deleted arm count" = "Deleted arm count",
    "Silent Mutation Rate" = "Silent mutation rate",
    "Nonsilent Mutation Rate" = "Nonsilent mutation rate",
    "SNV Neoantigens" = "SNV neoantigens"
  )
  raw_label <- unname(display_label[d$endpoint])
  raw_label[is.na(raw_label)] <- d$endpoint[is.na(raw_label)]
  raw_label[d$endpoint == "Aneuploidy score" & d$family == "aneuploidy"] <-
    "Aneuploidy (Taylor)"
  raw_label[d$endpoint == "Aneuploidy Score" & d$family == "thorsson"] <-
    "Aneuploidy (immune atlas)"
  duplicated_label <- duplicated(raw_label) | duplicated(raw_label, fromLast = TRUE)
  raw_label[duplicated_label] <- paste0(raw_label[duplicated_label], " [",
                                       d$family[duplicated_label], "]")
  d$label <- .short_endpoint(raw_label, 28L)
  d$site_sensitive <- grepl("^site-sensitive", d$site_robustness_status)
  d$corner_label <- sprintf("%s\npred. %s", d$label,
                            formatC(d$prediction, digits = 4, format = "fg"))
  n <- nrow(d)
  d$theta <- pi / 2 - 2 * pi * (seq_len(n) - 1L) / n
  d$x <- cos(d$theta) * d$reference_percentile / 100
  d$y <- sin(d$theta) * d$reference_percentile / 100
  d$label_x <- cos(d$theta) * 1.22
  d$label_y <- sin(d$theta) * 1.22
  d$hjust <- ifelse(cos(d$theta) > 0.25, 0,
                    ifelse(cos(d$theta) < -0.25, 1, 0.5))
  d$vjust <- ifelse(abs(cos(d$theta)) < 0.25, 0.5,
                    ifelse(sin(d$theta) > 0.25, 1,
                           ifelse(sin(d$theta) < -0.25, 0, 0.5)))
  polygon <- rbind(d, d[1L, , drop = FALSE])
  spokes <- data.frame(x = 0, y = 0, xend = cos(d$theta), yend = sin(d$theta))
  levels <- c(0.25, 0.50, 0.75, 1.00)
  rings <- do.call(rbind, lapply(levels, function(level) {
    z <- data.frame(level = level, theta = d$theta)
    z <- rbind(z, z[1L, , drop = FALSE])
    z$x <- cos(z$theta) * z$level
    z$y <- sin(z$theta) * z$level
    z
  }))
  ring_label_theta <- 1.90
  ring_labels <- data.frame(
    x = cos(ring_label_theta) * levels,
    y = sin(ring_label_theta) * levels,
    label = paste0(levels * 100, "%")
  )
  profile_colour <- .pathofm_colour(d)
  ggplot2::ggplot() +
    ggplot2::geom_path(data = rings,
                       ggplot2::aes(x = x, y = y, group = level),
                       color = "#D7DEE8", linewidth = 0.45) +
    ggplot2::geom_segment(data = spokes,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          color = "#D7DEE8", linewidth = 0.4) +
    ggplot2::geom_polygon(data = polygon, ggplot2::aes(x = x, y = y),
                          fill = profile_colour, alpha = 0.16,
                          color = profile_colour, linewidth = 1.05) +
    ggplot2::geom_point(data = d, ggplot2::aes(x = x, y = y),
                        color = "#F05D5E", fill = "white", shape = 21,
                        size = 3, stroke = 1) +
    ggplot2::geom_text(data = d,
                       ggplot2::aes(x = label_x, y = label_y,
                                    label = corner_label,
                                    hjust = hjust, vjust = vjust),
                       size = 3.0, lineheight = 1.02, color = "#26364A") +
    ggplot2::geom_text(data = ring_labels,
                       ggplot2::aes(x = x, y = y, label = label),
                       hjust = 1, vjust = -0.15, size = 2.4,
                       color = "#718096") +
    ggplot2::coord_equal(xlim = c(-1.62, 1.62), ylim = c(-1.48, 1.48),
                         clip = "off") +
    ggplot2::labs(
      title = "Continuous internally derived TCGA estimates",
      subtitle = paste(
        "Only endpoints predictable for this cancer and representation are shown.",
        "Corner labels give original predictions."
      )
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15,
                                         color = "#172B4D", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(color = "#5A6B7D", hjust = 0.5),
      plot.margin = ggplot2::margin(18, 28, 18, 28)
    )
}

#' Plot binary PLS-LDA outputs
#'
#' @param predictions Output from [predict_titan()] for one patient.
#' @return A ggplot object. Scores are displayed as TCGA out-of-fold score
#'   ranks explicitly labelled "not probability".
#' @export
plot_titan_binary <- function(predictions) {
  predictions <- .filter_predictable_predictions(predictions)
  d <- predictions[predictions$outcome_type == "binary" &
                     is.finite(predictions$reference_rank), , drop = FALSE]
  if (!nrow(d)) stop("No binary predictions with score ranks (not probability).", call. = FALSE)
  if (length(unique(d$patient_id)) != 1L) stop("Supply predictions for one patient.", call. = FALSE)
  raw_label <- ifelse(duplicated(d$endpoint) | duplicated(d$endpoint, fromLast = TRUE),
                      paste0(d$endpoint, " [", d$family, "]"), d$endpoint)
  d$label <- .short_endpoint(raw_label, 38L)
  d$site_sensitive <- grepl("^site-sensitive", d$site_robustness_status)
  d$label <- ifelse(d$site_sensitive,
                    paste0("[CODE-GROUPING-SENSITIVE] ", d$label), d$label)
  d$label <- paste0(d$label, "  (n=", d$training_n, "; +",
                    d$training_positive, "/-", d$training_negative, ")")
  d$label <- factor(d$label, levels = rev(d$label[order(d$reference_rank)]))
  d$class_label <- ifelse(d$predicted_class == "1", "Predicted positive", "Predicted negative")
  ggplot2::ggplot(d, ggplot2::aes(y = label, x = reference_rank, color = class_label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 50, xend = reference_rank, yend = label), color = "#CBD5E0", linewidth = 1) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(values = c("Predicted negative" = "#4C78A8", "Predicted positive" = "#F05D5E")) +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
    ggplot2::labs(title = "Binary research-model calls", subtitle = "Reference rank, not probability; [CODE-GROUPING-SENSITIVE] models require prominent caution", x = "Reference rank, not probability", y = NULL, color = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold", size = 15, color = "#172B4D"),
                   plot.subtitle = ggplot2::element_text(color = "#5A6B7D"))
}

#' Radar plot of PathoFMPred continuous outputs
#'
#' Representation-neutral alias for [plot_titan_radar()]. The selected
#' foundation model is carried in the prediction rows returned by
#' [predict_pathofm()]; no representation is selected inside the plotting
#' function.
#'
#' @inheritParams plot_titan_radar
#' @return A ggplot object.
#' @export
plot_pathofm_radar <- function(predictions, max_endpoints = 12L) {
  plot_titan_radar(predictions, max_endpoints = max_endpoints)
}

#' Plot PathoFMPred binary outputs
#'
#' Representation-neutral alias for [plot_titan_binary()]. Binary score ranks
#' remain explicitly uncalibrated ranks, not probabilities.
#'
#' @param predictions Output from [predict_pathofm()] for one patient.
#' @return A ggplot object.
#' @export
plot_pathofm_binary <- function(predictions) {
  plot_titan_binary(predictions)
}
