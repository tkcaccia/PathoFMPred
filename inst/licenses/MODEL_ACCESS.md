# Model-object access and upstream terms

The PathoFMPred R source code is distributed under the license stated in the
package DESCRIPTION. This does not replace or broaden the terms that govern
foundation-model weights, released embeddings, source pathology images, or
downstream fitted objects.

## Public repository

The public `tkcaccia/PathoFMPred` repository contains a minimal two-model
Giga-SSL fixture. Full Giga-SSL and Prov-GigaPath collections are optional,
explicit post-install downloads. Users remain responsible for reviewing and
complying with the upstream terms and dataset licenses.

The public repository does not contain a TITAN fitted-model collection and its
download function does not offer one.

## Access-controlled repository

`tkcaccia/PathoFMPred-private` contains the project reference collections:

* `pathofmpred_titan.rds`
* `pathofmpred_gigassl.rds`
* `pathofmpred_provgigapath.rds`

Access to the private repository does not grant a right to redistribute any
object. Existing collaborators retain access for the research collaboration.

## TITAN

Users who build a TITAN object locally must obtain the input representations
from the upstream TITAN repository and accept its terms. A locally fitted
`pathofmpred_titan.rds` must not be added to the public repository unless the
TITAN rights holder and the authors' institution have confirmed that
redistribution is permitted.

## Research status

The fitted models provide internal TCGA research estimates. They have not been
independently externally validated and are not clinical devices.
