#' @title Improved q-values for discrete uniform and homogeneous tests
#'
#'
#' @description This package implements five different versions of the q-value multiple testing procedure
#' proposed by Storey and Tibshirani (2003). The q-value method
#' is based on the false discovery rate (FDR); different versions of the q-value method can be defined depending
#' on the particular estimator used for the proportion of true null hypotheses, \eqn{\pi_0}{π_0}, which is plugged in the FDR formula. The
#' first version of the q-value uses the \eqn{\pi_0}{π_0} estimator in Storey (2002), with tunning parameter \eqn{\lambda}{lambda} = 0.5; whereas
#' the second version uses the \eqn{\pi_0}{π_0} estimator in Storey and Tibshirani (2003), which is based on an automatic method to
#' select the tunning parameter
#' \eqn{\lambda}{lambda}. These two methods are only appropriate when the P-values follow a continuous uniform distribution
#' under the global null hypothesis. This package also provides three other versions of the q-value for
#' homogeneous discrete uniform P-values, which often appear in practice.
#' The first discrete version of the q-value uses the \eqn{\pi_0}{π_0} estimator
#' proposed in Liang (2016).
#' The second discrete q-value method uses the estimator of \eqn{\pi_0}{π_0} proposed in Chen et al. (2014),
#' when simplified for the special case of homogeneous discrete P-values.
#' The third discrete version of the q-value employs a standard procedure but applied on randomized P-values.
#' Once the estimated q-values are computed, the q-value method
#' rejects the null hypotheses whose q-values are less than or equal to the nominal FDR level. All the versions of
#' the q-value method explained above can be seen in Cousido-Rocha et al. (2019).
#'
#'
#' @details
#' \itemize{
#' \item{Package: ‘DiscreteQvalue’}
#' \item{Version: 1.0}
#' \item{Maintainer: Marta Cousido Rocha \email{martacousido@@uvigo.es}}
#' \item{License: GPL-2}
#' }
#'
#' @return
#' \itemize{
#' \item{‘\link{DQ}’}
#' }
#'
#' @author
#' \itemize{
#' \item{Marta Cousido Rocha.}
#' \item{José Carlos Soage González.}
#' \item{Jacobo de Uña Álvarez.}
#' \item{Sebastian Döhler.}
#' }
#'
#'
#' @section Acknowledgements:
#' This work has received financial support of the Call 2015 Grants for PhD contracts for training of
#' doctors of the Ministry of Economy and Competitiveness, co-financed by the European Social Fund
#' (Ref. BES-2015-074958). The authors acknowledge support from MTM2014-55966-P project,
#' Ministry of Economy and Competitiveness, and MTM2017-89422-P project, Ministry of Economy,
#' Industry and Competitiveness, State Research Agency, and Regional Development Fund, UE. The
#' authors also acknowledge the financial support provided by the SiDOR research group through the
#' grant Competitive Reference Group, 2016-2019 (ED431C 2016/040), funded by the “Consellería
#' de Cultura, Educación e Ordenación Universitaria. Xunta de Galicia”.
#'
#'
#' @references
#'
#' \itemize{
#' \item{Chen, X., R. W. Doerge, and J. F. Heyse (2014). Methodology Multiple testing with discrete data: proportion of true null hypotheses and two adaptive FDR procedures. arXiv:1410.4274v2.}
#'
#' \item{Cousido-Rocha, M., J. de Uña-Álvarez, and S. Döhler (2019). Multiple testing methods for homogeneous discrete uniform P-values. Preprint.}
#'
#' \item{Liang, K. (2016). False discovery rate estimation for large-scale homogeneous discrete p-values. Biometrics 72, 639-648.}
#'
#' \item{Storey, J. D. (2002). A direct approach to false discovery rates. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64 (3), 479-498.}
#'
#' \item{Storey, J. and R. Tibshirani (2003). Statistical significance for genomewide studies. Proceedings of National Academy of Science 100, 9440-9445.}
#' }
#'
#'
"_PACKAGE"
#> [1] "_PACKAGE"
