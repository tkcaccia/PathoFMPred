# Model access and licensing boundary

The CRAN-oriented PathoFMPred source package contains code and documentation
under the MIT license. It does not include pretrained encoders, slide
embeddings, fitted PathoFMPred model collections, prediction-reference
distributions, the research atlas registry, pathology images, or molecular
source data.

The package can fit a collection from feature and outcome tables supplied by
the user. It can also load a locally available collection. Its optional
`fetch_pathofmpred_models()` helper includes a pinned URL and SHA-256 checksum
for a separate Prov-GigaPath coefficient collection. The asset is downloaded
only after the user calls the function. No data are downloaded at installation,
loading, checking, or vignette-build time. Other representations require a
user-supplied authorized URL and checksum or local construction.

The separately released Prov-GigaPath collection is licensed CC BY 4.0 with
attribution to the PathoFMPred authors, the Prov-GigaPath authors and the
released TCGA embedding dataset. Prov-GigaPath source code is Apache 2.0;
the TCGA embedding dataset is marked CC BY 4.0. The collection contains fitted
coefficients and aggregate registry metadata, but no foundation-model weights,
slide embeddings, patient rows, or patient-level prediction-reference scores.

The MIT license for this package does not grant rights in any separately
obtained encoder, embedding, fitted model, clinical data, or other third-party
material. Users must verify the access, use, attribution, and redistribution
terms for each input and fitted object before using or sharing it. In
particular, this package does not grant permission to redistribute objects
trained on TITAN outputs. The upstream terms are at
<https://github.com/mahmoodlab/TITAN#license-and-terms-of-use>.

The example code is illustrative and uses user-provided or synthetic data.
Prediction from these models is for research only. The package does not claim
external validation, calibration, clinical utility, or approval as a medical
device.
