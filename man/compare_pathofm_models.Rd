\name{compare_pathofm_models}
\alias{compare_pathofm_models}
\title{Compare representation-specific models for one cancer-endpoint task}
\usage{
compare_pathofm_models(cancer, endpoint, outcome_type = NULL, family = NULL,
  include_limited_evidence = FALSE)
}
\arguments{
\item{cancer}{TCGA cancer code.}
\item{endpoint}{Exact endpoint label.}
\item{outcome_type}{Optional continuous or binary filter.}
\item{family}{Optional endpoint-family filter.}
\item{include_limited_evidence}{Include limited-sample models only when they
crossed the representation-specific matched effect threshold.}
}
\value{A registry data frame with object-availability, matched-threshold,
resource-scope and selection-guidance fields.}
\description{
Lists controlled reference objects for the same task while distinguishing object
availability from representation-specific matched-atlas effect-threshold crossing.
The inventory is a partial operationalization of the multi-representation atlas.
The function does not automatically select the largest internal estimate or combine
representation-specific predictions.
}
