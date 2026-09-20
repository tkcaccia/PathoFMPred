.endpoint_glossary <- function(predictions) {
  p <- unique(predictions[, c("family", "endpoint", "outcome_type",
                              "output_units"), drop = FALSE])
  key <- paste(p$family, p$endpoint, sep = "::")
  meaning <- rep(NA_character_, nrow(p))
  source <- rep(NA_character_, nrow(p))

  put <- function(family, endpoint, text, citation) {
    hit <- key == paste(family, endpoint, sep = "::")
    meaning[hit] <<- text
    source[hit] <<- citation
  }
  taylor <- "Taylor et al., Cancer Cell 2018, Table S2; doi:10.1016/j.ccell.2018.03.007."
  thorsson <- "Thorsson et al., Immunity 2018, PanImmune_MS; doi:10.1016/j.immuni.2018.03.023."
  put("aneuploidy", "Aneuploidy score",
      "Number of chromosome arms with arm-level copy-number gain or loss; larger values indicate greater chromosomal imbalance.", taylor)
  put("aneuploidy", "Deleted arm count",
      "Number of chromosome arms classified as deleted; larger values indicate a greater arm-level deletion burden.", taylor)
  put("aneuploidy", "Genome doubling",
      "Binary whole-genome-doubling call derived from the tumour copy-number profile.", taylor)
  put("thorsson", "Aneuploidy Score",
      "Aneuploidy burden distributed with the TCGA immune landscape resource; larger values indicate more arm-level copy-number alterations.", thorsson)
  put("thorsson", "Lymphocyte Infiltration Signature Score",
      "Gene-expression signature summarising lymphocyte-associated transcriptional infiltration; it is not a direct cell count.", thorsson)
  put("thorsson", "TIL Regional Fraction",
      "Estimated fraction of tumour image regions containing tumour-infiltrating lymphocytes in the TCGA immune landscape resource.", thorsson)
  put("thorsson", "Nonsilent Mutation Rate",
      "Rate of coding mutations that alter the encoded protein. PathoFMPred reports the model output on the analysis log(1+x) scale.", thorsson)
  put("thorsson", "Silent Mutation Rate",
      "Rate of coding mutations that do not change the encoded amino acid. PathoFMPred reports the model output on the analysis log(1+x) scale.", thorsson)
  put("thorsson", "SNV Neoantigens",
      "Predicted neoantigen burden arising from single-nucleotide variants. PathoFMPred reports the model output on the analysis log(1+x) scale.", thorsson)
  put("microsatellite_instability", "MANTIS score",
      "Genome-wide microsatellite-instability score from the MANTIS algorithm; larger values indicate stronger MSI evidence.",
      "Kautto et al., Oncotarget 2017; doi:10.18632/oncotarget.20432. TCGA values were obtained from cBioPortal PanCancer Atlas records.")
  put("microsatellite_instability", "MSIsensor score",
      "Microsatellite-instability score from MSIsensor; larger values indicate a larger fraction of unstable microsatellite loci.",
      "Niu et al., Bioinformatics 2014; doi:10.1093/bioinformatics/btt755. TCGA values were obtained from cBioPortal PanCancer Atlas records.")
  put("microsatellite_instability", "MSI-H (MANTIS >0.4)",
      "Binary call indicating whether the source MANTIS score exceeds 0.4; the displayed LDA score is not a probability.",
      "Bonneville et al., JCO Precision Oncology 2017; doi:10.1200/PO.17.00073.")
  put("microsatellite_instability_sensitivity", "MSI-H strict (MANTIS >0.6)",
      "Sensitivity-analysis call: MANTIS >0.6 is positive, <0.4 is negative, and intermediate source values were excluded from training.",
      "Threshold documentation distributed with the cBioPortal TCGA PanCancer Atlas clinical sample files.")

  driver <- p$family == "driver_mutation" & is.na(meaning)
  meaning[driver] <- paste0(
    "Binary prediction of a qualifying protein-altering PASS mutation in ",
    p$endpoint[driver], "; the displayed LDA score is not a probability."
  )
  source[driver] <- "MC3 consensus calls: Ellrott et al., Cell Systems 2018, doi:10.1016/j.cels.2018.03.002; cancer-gene eligibility: Bailey et al., Cell 2018, doi:10.1016/j.cell.2018.02.060."
  pathway <- p$family == "oncogenic_pathway" & is.na(meaning)
  meaning[pathway] <- paste0(
    "Binary prediction that the ", p$endpoint[pathway],
    " oncogenic pathway carries at least one qualifying genomic alteration; the displayed LDA score is not a probability."
  )
  source[pathway] <- "Sanchez-Vega et al., Cell 2018, Table S4; doi:10.1016/j.cell.2018.03.035."
  other <- is.na(meaning)
  meaning[other] <- "Model output for the named endpoint. Interpret together with its family, analysed scale and source publication."
  source[other] <- "See the PathoFMPred model registry and the companion analysis source manifest."

  data.frame(
    family = p$family, endpoint = p$endpoint,
    label = ifelse(duplicated(p$endpoint) | duplicated(p$endpoint, fromLast = TRUE),
                   paste0(p$endpoint, " [", p$family, "]"), p$endpoint),
    meaning = meaning,
    scale = ifelse(p$outcome_type == "binary",
                   "PLS-LDA class and uncalibrated discriminant score",
                   p$output_units),
    source = source,
    stringsAsFactors = FALSE
  )
}
