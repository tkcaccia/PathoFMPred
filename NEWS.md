# PathoFMPred 0.1.0

First CRAN submission candidate.

- Fit portable cancer-specific PLS regression and PLS-LDA collections from
  user-supplied pathology representation and outcome tables.
- Link tables by a required patient ID, pool repeated feature rows by mean or
  median, and report unmatched IDs and endpoint-specific missingness.
- Apply locally fitted or separately authorized model collections to new
  representation features. Research reports distinguish uncalibrated binary
  scores from probabilities and identify limited-evidence models.
- Provide a vignette, reference manual and tests for building and using
  collections from TITAN, Giga-SSL, Prov-GigaPath or future representations.
- Distribute source code only. No pretrained encoders, embeddings, fitted
  model collections, patient-level data or research atlas registry are bundled.
  Optional downloads require an explicit user-supplied source and SHA-256
  checksum and never run during installation or package checks.
