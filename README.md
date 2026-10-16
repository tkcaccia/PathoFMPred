# PathoFMPred

PathoFMPred builds and applies portable collections of cancer-specific PLS
regression and PLS-LDA models from pathology foundation-model slide representations.
It links feature and outcome tables by an explicit patient identifier, pools
repeated slide rows at patient level, performs nested cross-validation, and
fits a final model for each eligible outcome.

The source code is available under the MIT license. The CRAN-oriented source
package contains no fitted model collection, atlas registry, slide embedding,
or patient-level data. Fitted collections obtained separately are governed by
their own terms; the package MIT license does not cover them.

These models provide internally derived TCGA research estimates. They have not
been externally validated and are not intended for diagnosis or treatment.
Binary scores and reference ranks are not probabilities.

See `inst/licenses/MODEL_ACCESS.md` for the package boundary. Earlier GitHub
model releases have been withdrawn from public access while their
redistribution rights are reviewed. Existing collaborators retain access
through the separate private repository.

## Install

```r
install.packages("fastPLS", repos = "https://cloud.r-project.org")
remotes::install_github("tkcaccia/PathoFMPred", dependencies = TRUE)
```

PathoFMPred checks the fastPLS provenance stored with each fitted model. The
current reference collections use fastPLS 0.3 from CRAN. Older collections
retain their Git revision metadata and may trigger a provenance warning.

## Use an authorized model collection

```r
library(PathoFMPred)

# Place an authorized collection at the path returned by:
tools::R_user_dir("PathoFMPred", which = "data")

# Alternatively, only when an authorized source has supplied both values:
fetch_pathofmpred_models(
  foundation_model = "GigaSSL",
  base_url = "https://authorized.example.org/release",
  expected_sha256 = c(GigaSSL = "<64-character SHA-256 digest>")
)
```

The URL and digest above are placeholders, not a working release. The user
must obtain the collection and digest under the applicable terms. Nothing is
downloaded during installation, package loading, checking, or vignette
construction. The helper stores verified objects in the user's R data
directory.

## Build a model collection

The first three arguments are mandatory. Outcome IDs must be unique. Feature
IDs may repeat and are pooled by patient mean or median. The function stops
when 100 or fewer unique IDs match and reports unmatched IDs plus missingness
for every outcome.

```r
features <- read.csv("slide_features.csv", check.names = FALSE)
outcomes <- read.csv("patient_outcomes.csv", check.names = FALSE)

object <- create_pathofmpred_object(
  feature_table = features,
  outcome_table = outcomes,
  id_column = "patient_id",
  aggregation = "mean",
  foundation_model = "my_representation",
  control = pathofmpred_control(
    components = 1:20,
    outer_folds = 5,
    inner_folds = 5,
    repeats = 5,
    method = "plssvd",
    binary_selection = "AUROC"
  ),
  output_file = "pathofmpred_my_representation.rds"
)

new_features <- read.csv("new_slide_features.csv", check.names = FALSE)
prediction <- predict_pathofmpred_object(
  object = object,
  features = new_features,
  id_column = "patient_id"
)
```

`rsvd_oversample` and `rsvd_power` remain `NULL` unless explicitly supplied, so
the installed fastPLS defaults are used.

For a binary outcome encoded with labels other than `0`/`1` or
`FALSE`/`TRUE`, define the event class explicitly. This prevents an
alphabetical factor order from silently reversing the scientific meaning:

```r
object <- create_pathofmpred_object(
  features, outcomes, "patient_id",
  outcome_types = c(KEAP1 = "binary"),
  positive_class = c(KEAP1 = "mutated")
)
```

## Rebuild the three representation collections locally

The [package vignette](vignettes/PathoFMPred-workflow.Rmd) gives non-running
examples for TITAN, Giga-SSL and Prov-GigaPath, including the official
embedding sources, conversion script, patient-ID linkage and fitting calls.
You must acquire each input under its own terms and prepare a patient-level
outcome table from the cited molecular sources. The function
`create_pathofmpred_object()` accepts any future slide representation with the same
ID and numeric-feature rules. A newly fitted collection is not necessarily
identical to the atlas reference collection unless its input snapshot, target
definitions, eligibility rules and fitting settings also match.

## TITAN policy

The CRAN source package does not distribute `pathofmpred_titan.rds`. Users may
obtain the official TCGA slide representations from the gated
[MahmoodLab/TITAN](https://huggingface.co/MahmoodLab/TITAN) repository and build
their own object with `create_pathofmpred_object()`. The vignette documents the
feature conversion and molecular, immune, pathway, aneuploidy, fusion, and MSI
sources. The resulting local object can be applied with
`predict_pathofmpred_object()`. Locally built TITAN objects should not be
redistributed without permission from the upstream rights holder.

The access-controlled `tkcaccia/PathoFMPred-private` repository contains the
project's research reference objects. Its access terms are separate from the
CRAN package licence.

## Apply an authorized fitted collection

```r
library(PathoFMPred)

available_pathofm_models("GigaSSL", "COAD") # after installing the collection
x <- example_pathofm_features("COAD", "GigaSSL")
pred <- predict_pathofm(
  cancer = "COAD",
  features = x,
  foundation_model = "GigaSSL",
  patient_id = "software_example"
)
pred[, c("endpoint", "outcome_type", "prediction",
         "predicted_class", "lda_score")]
```

No model is available in a fresh CRAN installation. To run a fully local
example without any third-party model, use `create_pathofmpred_object()` with
your own feature and outcome tables, then `predict_pathofmpred_object()`.
See `vignette("PathoFMPred-workflow")` for the workflow and licensing boundary.
