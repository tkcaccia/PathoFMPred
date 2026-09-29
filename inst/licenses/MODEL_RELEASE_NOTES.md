# PathoFMPred optional model collections v3

This release contains the full Giga-SSL and Prov-GigaPath PathoFMPred fitted
model collections for explicit post-install download:

* `pathofmpred_gigassl.rds`, 795 fitted models
  SHA-256 `874b84cb040d7c955b96f95b29ee181c9056f5e40a3feda28171f942063d0745`
* `pathofmpred_provgigapath.rds`, 795 fitted models
  SHA-256 `72868d9fb7562f3eafaf7cc287cb4e909e1767587662be07625c985a028da832`

Both collections were rebuilt with the CRAN fastPLS 0.3 source release and
the documented 1 to 20 component analysis settings.

These assets are distributed separately from the MIT-licensed PathoFMPred
source code. To the extent that the PathoFMPred authors hold rights in the
fitted coefficients and accompanying metadata, the assets are provided under
CC BY 4.0. This grant does not replace or broaden upstream terms. Users must
retain attribution to PathoFMPred, Giga-SSL or Prov-GigaPath as applicable,
TCGA, and the endpoint-source studies recorded in the model registry.

Giga-SSL source code is released under MIT:
<https://github.com/trislaz/gigassl>

Prov-GigaPath source code is released under Apache-2.0:
<https://github.com/prov-gigapath/prov-gigapath>

The `seandavis/tcga_provgigapath_embeddings` dataset used by this project is
marked CC BY 4.0:
<https://huggingface.co/datasets/seandavis/tcga_provgigapath_embeddings>

The assets contain fitted parameters and provenance, not patient-level feature
or outcome rows. They provide internally validated TCGA research estimates,
have not undergone independent external validation, and are not clinical
models.

No TITAN-derived fitted parameters are included. TITAN-derived objects remain
excluded from public redistribution pending written permission from the
upstream rights holder.
