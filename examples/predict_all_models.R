library(PathoFMPred)

# This example requires an authorized Giga-SSL collection installed locally.
x <- example_pathofm_features("COAD", "GigaSSL")
predictions <- predict_pathofm(
  cancer = "COAD", features = x, foundation_model = "GigaSSL",
  patient_id = "example_patient"
)
head(predictions)

# Obtain each collection under its applicable terms and install it locally.
# No fitted object or model download source is included in the CRAN package.
