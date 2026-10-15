# PathoFMPred quick tutorial

## 1. Check for an authorized local collection

```r
library(PathoFMPred)
available_pathofm_models("GigaSSL", "COAD")
```

The CRAN-oriented source package bundles no fitted models. This listing is
empty until a user installs an authorized collection.

## 2. Install a collection supplied under its applicable terms

```r
fetch_pathofmpred_models(
  foundation_model = "GigaSSL",
  base_url = "https://authorized.example.org/release",
  expected_sha256 = c(GigaSSL = "<SHA-256 supplied by the source>")
)
```

The URL and digest are placeholders. The helper has no built-in model source
and runs only when called. TITAN models are not available through it.

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
