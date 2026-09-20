# PathoFMPred optional model collections v2

This release contains the full Giga-SSL and Prov-GigaPath PathoFMPred fitted
model collections for explicit post-install download:

* `pathofmpred_gigassl.rds`
  SHA-256 `3ccfb101a24cc1ac9c6564d42c3392a0a139ccddd052cefb326db47d1d54ae1c`
* `pathofmpred_provgigapath.rds`
  SHA-256 `712b39f122f5efddc8d182700d84f47314943870de42f5596aa4bec6daeb67c8`

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
