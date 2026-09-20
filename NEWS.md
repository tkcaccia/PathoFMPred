# PathoFMPred 0.3.0

- Relicensed contributor-authored source code under MIT while explicitly
  separating fitted-model, registry, reference-data and upstream terms.
- Pinned the required fastPLS 0.3 Git revision for reproducible installation
  and continuous integration.
- Runtime checks now enforce the fastPLS version and enforce its Git revision
  when R records `RemoteSha`; source-archive installations without that
  optional metadata receive a single provenance warning instead of a false
  incompatibility error.
- Added explicit binary event-class mapping to prevent factor-order reversal.
- Existing downloaded model collections are now checksum-verified when
  `verify = TRUE`, not only when newly downloaded.
- Added `predict_pathofmpred_object()` so locally created or access-controlled
  model collections can be applied without being installed as bundled package
  data.

# PathoFMPred 0.2.1

- Corrected model selection so every cancer-specific report, prediction table,
  binary plot and radar plot includes only endpoints that crossed the matched
  threshold for the selected representation. Limited-evidence opt-in can add
  only smaller-sample threshold-crossing models. Tested-below-threshold fitted
  objects remain in the raw registry for audit and never enter inference.
- Corrected radar plots to use the representation-specific predictable endpoint
  set and to label each corner with the original model prediction.

# PathoFMPred 0.2.0

- Renamed the package and primary API to PathoFMPred and added explicit
  `foundation_model` selection for TITAN, Giga-SSL and Prov-GigaPath schemas and
  representation-specific fitted objects.
- Aligned the matched binary benchmark around one estimand. Pooled inner
  out-of-fold AUROC selects the component count, outer out-of-fold AUROC is the
  primary statistic, and AUROC at least 0.60 defines a descriptive crossing.
  A training-only balanced-accuracy threshold supplies class calls,
  sensitivity, specificity, PPV and NPV. Reports include PR-AUC, prevalence and
  its no-skill reference. Empirical-prior and equal-prior calls remain
  secondary operating-rule sensitivities, not calibrated probabilities.
- Added `compare_pathofm_models()` so users can inspect representation-specific
  evidence without automatic post hoc model selection or averaging.
- Added representation-neutral `plot_pathofm_radar()` and
  `plot_pathofm_binary()` helpers while retaining the TITAN-named plotting
  functions as backwards-compatible aliases.
- Clarified that source, registries, schemas and synthetic fixtures are distinct
  from access-controlled fitted coefficient objects and their upstream terms.

# PathoFMPred 0.1.0

- Reports and prediction outputs now distinguish initial nested-CV primary
  screening estimates from means across five additional nested-CV repeats.

- Binary models with fewer than 50 TCGA participants in either development
  class are now excluded from default inference and require explicit opt-in.
- Added held-out average precision, cohort-specific PPV/NPV, fold class counts,
  component-selection distributions, learning curves and repeat prediction
  stability metadata.
- Report banners, figures and tables now place "not probability" immediately
  beside every binary TCGA out-of-fold score-rank label and make the absence of
  external validation prominent.
- Removed treatment and outcome narratives from the optional report context;
  reports accept descriptive pathology/source context only.
- Replaced effect-sounding tier descriptions in reports and documentation with
  the neutral labels "documented screening tier A/B".

- Bundles 323 cancer-specific TCGA discovery models: 219 continuous PLS
  regressions and 104 binary PLS-LDA classifiers across 27 cancers.
- Adds cancer-vector dispatch and patient-level pooling through `predict_titan()`.
- Adds HTML/PDF single-sample reports with continuous radar profiles, binary
  calls, reference percentiles, provenance and out-of-distribution diagnostics.
- Includes SHA-256 model registry and repeated out-of-fold prediction reference
  distributions without patient identifiers or patient-level training rows.
