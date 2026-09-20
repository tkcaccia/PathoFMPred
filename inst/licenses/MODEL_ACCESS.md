# Model-asset access, licensing and upstream terms

The PathoFMPred source code and documentation authored by the package
contributors are distributed under the MIT license. The MIT grant expressly
excludes fitted model objects, model registries and prediction-reference data,
foundation-model embeddings, pathology images, molecular source data and
other third-party material.

The model registry and original result summaries are made available under
CC BY 4.0, subject to attribution of PathoFMPred and the cited source studies.
Fitted model objects have the representation-specific conditions below.
The CC BY 4.0 legal code is available at
<https://creativecommons.org/licenses/by/4.0/legalcode>.

## Public repository

The public `tkcaccia/PathoFMPred` repository contains a minimal two-model
Giga-SSL fixture. Full Giga-SSL and Prov-GigaPath collections are optional,
explicit post-install downloads. To the extent that the PathoFMPred authors
hold rights in the fitted coefficients and accompanying metadata, these
downstream objects are provided under CC BY 4.0. This permission does not
replace or broaden any upstream terms. Users must retain attribution to
PathoFMPred, Giga-SSL or Prov-GigaPath as applicable, the TCGA source data and
the endpoint-source studies recorded in the registry.

Giga-SSL source code is released under MIT at
<https://github.com/trislaz/gigassl>. Prov-GigaPath source code is released
under Apache-2.0 at <https://github.com/prov-gigapath/prov-gigapath>, and the
`seandavis/tcga_provgigapath_embeddings` dataset used by this project is
marked CC BY 4.0 at
<https://huggingface.co/datasets/seandavis/tcga_provgigapath_embeddings>.
The upstream notices and citations remain applicable.

The public repository does not contain a TITAN fitted-model collection and its
download function does not offer one.

## Access-controlled repository

`tkcaccia/PathoFMPred-private` contains the project reference collections:

* `pathofmpred_titan.rds`
* `pathofmpred_gigassl.rds`
* `pathofmpred_provgigapath.rds`

Access to the private repository does not grant a right to redistribute the
TITAN object. Existing collaborators retain access for the research
collaboration. The Giga-SSL and Prov-GigaPath objects are governed by the
public-object terms above.

## TITAN

Users who build a TITAN object locally must obtain the input representations
from the upstream TITAN repository and accept its terms. A locally fitted
`pathofmpred_titan.rds` must not be added to the public repository unless the
TITAN rights holder and the authors' institution have confirmed that
redistribution is permitted.

The upstream TITAN terms describe models trained on TITAN outputs as
derivatives, restrict them to non-commercial academic research and prohibit
redistribution without permission. The TITAN object is therefore excluded
from the PathoFMPred MIT and CC BY 4.0 grants.
The current upstream terms are published at
<https://github.com/mahmoodlab/TITAN#license-and-terms-of-use>.

## Research status

The fitted models provide internal TCGA research estimates. They have not been
independently externally validated and are not clinical devices.
