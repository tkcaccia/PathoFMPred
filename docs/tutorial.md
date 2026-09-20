# PathoFMPred quick tutorial

## 1. Inspect the minimal public example

```r
library(PathoFMPred)
available_pathofm_models("GigaSSL", "COAD")
```

The public package bundles only two small Giga-SSL COAD models. They are test
fixtures for the software interface, not a complete model atlas.

## 2. Download optional full public collections

```r
fetch_pathofmpred_models(c("GigaSSL", "ProvGigaPath"))
```

This is an explicit post-install operation. TITAN models are not available
through the public downloader.

## 3. Predict with the representation you possess

```r
x <- read.csv("my_gigassl_features.csv", check.names = FALSE)
pred <- predict_pathofm(
  cancer = "COAD",
  features = x,
  foundation_model = "GigaSSL",
  patient_id = x$patient_id
)
```

Repeat a patient identifier for multiple slides. Prediction pools those rows by
patient mean before applying a fitted model.

## 4. Build an object for TITAN or a future representation

```r
features <- read.csv("foundation_model_features.csv", check.names = FALSE)
outcomes <- read.csv("patient_outcomes.csv", check.names = FALSE)

object <- create_pathofmpred_object(
  features,
  outcomes,
  "patient_id",
  aggregation = "mean",
  foundation_model = "TITAN",
  control = pathofmpred_control(
    components = 1:20,
    outer_folds = 5,
    inner_folds = 5,
    repeats = 5,
    method = "plssvd",
    binary_selection = "AUROC"
  ),
  output_file = "pathofmpred_titan.rds"
)
```

The outcome table must contain unique IDs. The feature table may contain
repeated IDs and can be pooled by mean or median. More than 100 IDs must match.
The returned `data_audit` records unmatched IDs and endpoint-specific
missingness. The fitted object stores no original feature or outcome rows.

Do not place `pathofmpred_titan.rds` in the public repository. Consult the
installed vignette for upstream TITAN access and outcome-source details.
