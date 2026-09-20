# PathoFMPred

PathoFMPred builds and applies portable collections of cancer-specific PLS
regression and PLS-LDA models from pathology foundation-model representations.
It links feature and outcome tables by an explicit patient identifier, pools
repeated slide rows at patient level, performs nested cross-validation, and
fits a final model for each eligible outcome.

The source code is available under the MIT license and is licensed separately
from the fitted model collections, registry and reference data. The
public package contains only two small Giga-SSL demonstration models. It does
not contain TITAN-derived fitted parameters. Full Giga-SSL and Prov-GigaPath
collections are downloaded only after an explicit user call.

These models provide internally derived TCGA research estimates. They have not
been externally validated and are not intended for diagnosis or treatment.
Binary scores and reference ranks are not probabilities.

See `inst/licenses/MODEL_ACCESS.md` for the CC BY 4.0 terms applied to original
registry/result material and the separate representation-specific conditions
for fitted objects. TITAN-derived objects are expressly excluded from the MIT
and CC BY 4.0 grants.

## Install

```r
remotes::install_github("tkcaccia/fastPLS", upgrade = "never")
remotes::install_github("tkcaccia/PathoFMPred", dependencies = TRUE)
```

## Download optional public model collections

```r
library(PathoFMPred)
fetch_pathofmpred_models(c("GigaSSL", "ProvGigaPath"))
```

Nothing is downloaded during installation, package loading, checking, or
vignette construction. `fetch_pathofmpred_models()` stores validated objects in
the user's R data directory. The download location can be changed with its
`destination` or `base_url` arguments.

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

## TITAN policy

The public repository does not distribute `pathofmpred_titan.rds`. Users may
obtain the official TCGA representations from the gated
[MahmoodLab/TITAN](https://huggingface.co/MahmoodLab/TITAN) repository and build
their own object with `create_pathofmpred_object()`. The vignette documents the
feature conversion and molecular, immune, pathway, aneuploidy, fusion, and MSI
sources. The resulting local object can be applied with
`predict_pathofmpred_object()`. Locally built TITAN objects should not be
redistributed without permission from the upstream rights holder.

The access-controlled `tkcaccia/PathoFMPred-private` repository contains the
project's TITAN, Giga-SSL, and Prov-GigaPath reference objects.

## Minimal public example

```r
library(PathoFMPred)

available_pathofm_models("GigaSSL", "COAD")
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

See `vignette("PathoFMPred-workflow")` for the complete workflow and licensing
boundaries.
