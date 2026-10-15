# PathoFMPred rights review and permission requests

Prepared 15 October 2026. These drafts are for the maintainer to send, not a
legal determination. Do not treat silence as permission. Save replies and any
licence documents with the project release record.

## What needs a written answer

The proposed CRAN 0.5.0 source archive contains MIT-licensed code and
documentation only. It does not contain the fitted Giga-SSL, Prov-GigaPath or
TITAN objects, the research atlas registry, slide embeddings, or patient rows.
Its optional downloader has no built-in source URL or digest. CRAN submission
therefore need not wait for upstream permission to ship those excluded assets,
but the maintainer must still confirm that all package contributors agreed to
the MIT release and that no restricted third-party code or data remain.

The former public `models-v1`, `models-v2` and `models-v3` releases were
converted to drafts on 15 October 2026, so their Giga-SSL and Prov-GigaPath
assets are not publicly accessible through those releases. The requests below
seek clarification before any future public model release or proposal to
bundle models with CRAN. The TITAN object remains excluded from public
release. CRAN policy treats components downloaded during package use as part
of the rights audit when the package specifies their source.

Ask the relevant party to identify the legal rights holder and to confirm in
writing whether it authorizes worldwide, durable public redistribution of
downstream fitted coefficients, metadata and a small demonstration fit; whether
the authorization covers GitHub release assets and CRAN mirrors; whether
commercial users are permitted; required attribution and notices; and whether
the proposed CC BY 4.0 asset licence is compatible with upstream terms. A
reply from an author who cannot bind the rights holder should be referred to
the appropriate institutional licensing office.

## 1. Giga-SSL

To: Tristan Lazard <trislaz.tl@gmail.com>

Optional copy: Thomas Walter <thomas.walter@mines-paristech.fr>

Current contact: <https://trislaz.github.io/>. The paper gives Thomas
Walter's address: <https://openaccess.thecvf.com/content/CVPR2023W/CVMI/papers/Lazard_Giga-SSL_Self-Supervised_Learning_for_Gigapixel_Images_CVPRW_2023_paper.pdf>.

Subject: Permission clarification for downstream Giga-SSL PLS model release

Dear Dr Lazard,

We have trained cancer-specific PLS models on the released Giga-SSL TCGA
slide representations for a research benchmark. The fitted objects contain
coefficients, feature order and model metadata, but no original slide
embeddings or patient-level molecular outcome rows. We previously made a
GitHub research release containing 795 Giga-SSL fitted models, but withdrew
it from public access pending a rights review. We may later wish to
distribute these objects through a GitHub release or a small example through
an MIT-licensed R package on CRAN. The fitted objects would be licensed
separately, if allowed.

Could you please confirm who holds the rights to the released Giga-SSL TCGA
embeddings and whether models fitted on those embeddings may be redistributed
publicly? In particular, does the permission cover GitHub releases and CRAN
mirrors, commercial as well as non-commercial downstream users, and release of
the fitted objects under CC BY 4.0 with attribution? Please specify any
required notices or restrictions, and whether your repository's MIT licence
applies to the released embedding data as well as its source code. If another
institution must approve, would you please direct us to the appropriate
rights holder or licensing office?

We can send the exact object schema and release link for review. We would be
grateful for a written response that we can retain with the project records.

Kind regards,
Stefano Cacciatore
International Centre for Genetic Engineering and Biotechnology
stefano.cacciatore@icgeb.org

## 2. Prov-GigaPath

To: Hoifung Poon <hoifung@microsoft.com>; Carlo Bifulco
<carlo.bifulco@providence.org>

Optional copy: Sheng Wang <swang@cs.washington.edu>

The Nature paper lists these corresponding-author addresses:
<https://pubmed.ncbi.nlm.nih.gov/38778098/>. The model page is
<https://huggingface.co/prov-gigapath/prov-gigapath>.

Subject: Clarification of redistribution rights for downstream Prov-GigaPath models

Dear Dr Poon and Dr Bifulco,

We trained cancer-specific PLS models using the publicly posted TCGA
Prov-GigaPath embeddings from
<https://huggingface.co/datasets/seandavis/tcga_provgigapath_embeddings>.
The fitted objects contain coefficients and metadata, not the original slide
embeddings or patient-level outcomes. An earlier GitHub research release
contained 795 such fitted models, but we withdrew it from public access
pending a rights review. Our planned CRAN R package is code-only while we
clarify the rights to the fitted objects.

The model page identifies Apache-2.0 but also requires acceptance of access
conditions and describes research-only intended use. Could you clarify whether
the fitted coefficients are considered derivatives or otherwise subject to the
model access terms? May they be publicly redistributed on GitHub, potentially
through CRAN mirrors, and made available to commercial as well as academic
users? Is CC BY 4.0 an acceptable separate licence for our downstream fitted
objects, and what attribution, notices or user acceptance would be required?
Please also tell us whether Microsoft, Providence, or another rights holder
must issue formal approval.

We can provide the object schema and release link. A written response that
identifies the rights holder and scope of permission would help us document
the release accurately.

Kind regards,
Stefano Cacciatore
International Centre for Genetic Engineering and Biotechnology
stefano.cacciatore@icgeb.org

## 3. TCGA Prov-GigaPath embeddings dataset

To: Sean Davis <seandavi@gmail.com>

The dataset page marks the embeddings CC BY 4.0:
<https://huggingface.co/datasets/seandavis/tcga_provgigapath_embeddings>.
The address is published in the author's project page:
<https://github.com/seandavi/RBiocBook/blob/main/about.qmd>.

Subject: Provenance and licence scope of TCGA Prov-GigaPath embeddings

Dear Dr Davis,

We used your `tcga_provgigapath_embeddings` dataset to fit downstream PLS
models. We see that the Hugging Face dataset is marked CC BY 4.0. Could you
confirm that this licence applies to the embedding values and metadata used
for model training, and let us know whether you are aware of any additional
Prov-GigaPath or TCGA terms affecting public redistribution of fitted
downstream coefficients? Please provide the preferred attribution and a
stable dataset revision or citation. We are asking the Prov-GigaPath rights
holders separately about their model terms.

Kind regards,
Stefano Cacciatore
stefano.cacciatore@icgeb.org

## 4. TITAN, only if public redistribution is later proposed

To: Faisal Mahmood <faisalmahmood@bwh.harvard.edu>

Optional copy: Tong Ding <tong_ding@g.harvard.edu>

The upstream contact list and terms are at
<https://huggingface.co/MahmoodLab/TITAN> and
<https://github.com/mahmoodlab/TITAN#license-and-terms-of-use>.

Subject: Request for written permission concerning PLS models trained on TITAN outputs

Dear Professor Mahmood,

Our research project trained cancer-specific PLS models on the released TCGA
TITAN slide representations. We understand that the TITAN terms describe
models trained on TITAN outputs as derivatives and restrict redistribution.
We have excluded all TITAN fitted objects from the proposed public CRAN
package. We seek written clarification before considering any broader release.

Would your institution permit us to share the fitted coefficient objects
publicly, or only with named academic collaborators under controlled access?
If permission is possible, please specify the approved channels, audience,
licence, attribution, onward-sharing restrictions and any required agreement
for each user. We will not treat this request as permission until your rights
holder confirms it in writing.

Kind regards,
Stefano Cacciatore
stefano.cacciatore@icgeb.org

## 5. Institutional sign-off and author agreement

To: Martina Viviani <viviani@icgeb.org>, ICGEB Technology Transfer Office

The ICGEB office and contact are listed at
<https://vod.icgeb.org/icgeb_courses/degludec/>.

Subject: Rights review before CRAN submission of PathoFMPred

Dear Ms Viviani,

We are preparing PathoFMPred 0.5.0 for CRAN under MIT. The proposed CRAN
archive is code-only: it excludes pretrained encoders, embeddings, fitted
research models, the atlas registry and patient-level data. It can fit models
from user-provided tables and load independently obtained collections. The
earlier GitHub model releases have been withdrawn from public access while
rights are reviewed. All three reference collections remain in an
access-controlled private repository.

Could the Technology Transfer Office review (1) ownership of the package code
and documentation, including contributions from authors at other institutions;
(2) whether all listed package authors should approve MIT distribution;
(3) whether any model release can be restored after upstream permission is
documented; (4) the proposed CC BY 4.0 terms for downstream
fitted objects; and (5) any TCGA, upstream encoder or institutional conditions
that affect the release? We would appreciate written advice or referral to
the appropriate legal office before CRAN submission.

Kind regards,
Stefano Cacciatore
stefano.cacciatore@icgeb.org

For CRAN, retain a short affirmative email from every person whose work is
credited in `Authors@R`, or documented licence evidence for contributed
material. Suggested wording: "I approve inclusion of my contribution in
PathoFMPred under the MIT licence and authorize Stefano Cacciatore to submit
the package to CRAN."
