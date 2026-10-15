# Model access and licensing boundary

The CRAN-oriented PathoFMPred source package contains code and documentation
under the MIT license. It does not include pretrained encoders, slide
embeddings, fitted PathoFMPred model collections, prediction-reference
distributions, the research atlas registry, pathology images, or molecular
source data.

The package can fit a collection from feature and outcome tables supplied by
the user. It can also load a locally available collection. Its optional
`fetch_pathofmpred_models()` helper requires the user to provide an explicit
source URL and an independently obtained SHA-256 checksum. Neither model URLs
nor checksums are provided by the source package. No data are downloaded at
installation, loading, checking, or vignette-build time.

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
