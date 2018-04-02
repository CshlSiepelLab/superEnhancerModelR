#' Deletion data from the wap locus
#'
#' The expression data from Shin et al. (2016) for their deletion of subsets of the Wap enhancers.
#' Expression data is already normalized to the wild type. See paper for more details
#'
#' @format A data frame with 5 variables.
#' \describe{
#'   \item{condition}{A string describing the experimental condition}
#'   \item{expression}{the normalized expression level}
#'   \item{E*}{The 0/1 (absence/presence) coding for each enhancer}
#' }
#' @source \url{http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3606.html}
"wap"

#' Deletion data from the alpha-globin locus
#'
#' The expression data from Hay et al. (2016) for their deletion of subsets of the alpha-globin enhancers.
#' Expression data is already normalized to the wild type. See paper for more details
#'
#' @format A data frame with 5 variables.
#' \describe{
#'   \item{condition}{A string describing the experimental condition}
#'   \item{expression}{the normalized expression level}
#'   \item{E*}{The 0/1 (absence/presence) coding for each enhancer}
#' }
#' @source \url{http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3605.html}
"alpha.globin"
