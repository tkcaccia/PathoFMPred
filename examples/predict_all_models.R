library(PathoFMPred)

# The public package includes a two-model Giga-SSL software fixture.
x <- example_pathofm_features("COAD", "GigaSSL")
predictions <- predict_pathofm(
  cancer = "COAD", features = x, foundation_model = "GigaSSL",
  patient_id = "example_patient"
)
head(predictions)

# Full optional Giga-SSL and Prov-GigaPath collections are installed only after
# the user explicitly requests them.
# fetch_pathofmpred_models(c("GigaSSL", "ProvGigaPath"))
