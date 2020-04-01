#' @title Improved q-values for discrete uniform and homogeneous tests
#' @aliases DQ
#' @description Performs the five versions of the q-value method considered in Cousido-Rocha et al. (2019). The q-value method is based on the false discovery rate (FDR), and the  versions differ in the  estimator of the proportion of true null hypotheses, \eqn{\pi_0}{π_0}, which is plugged in the FDR estimator. Specifically, we consider as possible estimators for \eqn{\pi_0}{π_0}: two usual estimators for continuous and possibly heterogeneous P-values; an estimator for discrete P-values defined in two steps: firstly the P-values are randomized, and then the usual \eqn{\pi_0}{π_0} estimator for continuous P-values is applied; and the estimators recently proposed for discrete P-values by Liang (2016) and Chen et al. (2014).
#'
#'
#' @param pv A vector of P-values.
#' @param ss Support of the discrete distribution of the P-values. Only required for “Liang”, “Chen” and “Rand” methods which are specifically proposed for discrete P-values. If the P-values are continuous the methods “ST” and “SS” do not need this argument, hence “ss=NULL” by default.
#' @param ss_inf Logical. Default is FALSE. A variable indicating whether the support of the discrete distribution of the P-values is finite or infinite. See details.
#' @param method The q-value method. By default the “Chen” method is computed.
#' @param lambda The value of the tuning parameter to estimate \eqn{\pi_0}{π_0} when the method is “ST”. See details.
#'
#' @details The function implements the five different versions of the q-value method in Cousido-Rocha et al. (2019).
#' Three versions are novel adaptions for the case of homogeneous discrete uniform P-values, whereas the other two
#' are classical versions of the q-value method for P-values which follow a continuous uniform distribution under the global null hypothesis.
#' The classical versions are the q-value method based on the \eqn{\pi_0}{π_0} estimator proposed in Storey (2002) with tunning parameter \eqn{\lambda}{lambda} = 0.5,
#' and the q-value procedure which uses the  \eqn{\pi_0}{π_0} estimator in Storey and Tibshirani (2003) who proposed an automatic method to estimate \eqn{\pi_0}{π_0}.
#' We refer to these methods as “SS” and “ST”, respectively. On the other hand the three adaptations of the q-value method
#' for homogeneous discrete uniform P-values are: “Liang” which considers the  \eqn{\pi_0}{π_0} estimator proposed in
#' Liang (2016); “Chen” which uses a simplification for homogeneous discrete P-values of the algorithm for
#' the estimation of \eqn{\pi_0}{π_0} proposed in Chen et al. (2014); and “Rand” which employs the standard q-value procedure but applied to
#' randomized P-values. For details of the different q-value versions, in
#' particular for the novel adaptations for homogeneous discrete uniform P-values, see Cousido-Rocha et al. (2019).
#'
#' As we mentioned above the novel adaptations of the q-value method are developed for homogeneous discrete uniform P-values.
#' Specifically, suppose that we test a large number of null hypothesis, p, and that the P-values \eqn{\{pv_1, \dots, pv_p\}}{{pv_1, ... , pv_p}}
#' are observations of the random variables  \eqn{PV_i, i = 1, \dots , p.}{PV_i, i = 1, ... , p.}
#' Homogeneous means that all the P-values share an identical support \eqn{S}{S} with \eqn{0<t_1 < \dots <t_s<t_{s+1} \equiv 1}{0 < t_1 < ··· < t_s < t_{s+1} ≡ 1}.
#' On the other hand, making an abuse of language, we say that the P-values follow a discrete uniform distribution if it holds
#' \eqn{Pr(PV_i \leq t) = t}{Pr(PV_i ≤ t) = t} for \eqn{t \in S, i = 1, \dots, p}{t in S, i = 1, ..., p}.
#' The classical discrete uniform distribution is a particular case.
#'
#' The argument \code{“lambda”} must be a sequence of values in \eqn{[0, 1)}{[0, 1)}, for details of this parameter see Storey (2002) or Storey and Tibshirani (2003). The latter paper recommends the default value “lambda=seq(0.05,0.95,0.05)”.
#'
#' The support of the discrete distribution of the P-values can be finite or infinte. Hence the parameter “ss inf” must be “FALSE” if the support is finite and “TRUE” if the support is infinite. See examples where a poisson setting is considered.
#'
#' Finally, it is relevant to mention that Cousido-Rocha et al. (2019) verified (via simulations) that the results of the different q-values methods for dependent P-values are very similar to the ones corresponding to the independent setting.
#'
#'
#' @return A list containing the following components:
#' \item{pi0}{An estimate of the proportion of null P-values.}
#' \item{q.values}{A vector of the estimated q-values (the main quantities of interest).}
#'
#'
#'
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{José Carlos Soage González}
#' \item{Jacobo de Uña-Álvarez}
#' \item{Sebastian Döhler}
#' }
#'
#'
#' @references
#'
#' \itemize{
#' \item{Chen, X., R. W. Doerge, and J. F. Heyse (2014). Methodology Multiple testing with discrete data: proportion of true null hypotheses and two adaptive FDR procedures. arXiv:1410.4274v2.}
#'
#' \item{Cousido-Rocha, M., J. de Uña-Álvarez, and S. Döhler (2019). Multiple testing methods for homogeneous discrete uniform P-values. Preprint.}
#'
#' \item{Gibbons, J. D. and S. Chakraborti (1992). Nonparametric Statistical Inference. Third Edition. Marcel Dekker, Inc, New York.}
#'
#' \item{Liang, K. (2016). False discovery rate estimation for large-scale homogeneous discrete p-values. Biometrics 72, 639-648.}
#'
#' \item{Storey, J. D. (2002). A direct approach to false discovery rates. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 64 (3), 479-498.}
#'
#' \item{Storey, J. and R. Tibshirani (2003). Statistical significance for genomewide studies. Proceedings of National Academy of Science 100, 9440-9445.}
#' }
#'
#'
#'
#' @examples
#'
#' # We consider a simple simulated data set to illustrate the use of the DQ function.
#' # We have simulated the following situation.
#' # We have two groups, for example, 5 patients with tumor 1 and 5 patients
#' # with tumor 2. For each patient 100 variables are measured, for example,
#' # gene expression levels. Besides, the distributions of 30 of the variables
#' # are different in the two groups, and the differences are in location.
#'
#'
#' # We consider a collection of densities {f1=N(0,1), f2=N(0,1/4), f3=N(2,1), f4=N(2,1/4)}. In
#' # the first group (tumor 1) the sample of each variable (gene) comes from one of
#' # the four densities with the same probability 1/4. On the other hand, in the second
#' # group the sample of each variable comes from the same density as in the first
#' # group except for 30 randomly selected variables for which the density changes
#' # producing location differences. Specifically, if the variable follows
#' # f1 in the first group, its density, in the second group, is f3  producing a
#' # change on its location parameter.The situation for the other cases is as follows:
#' # the density f2 (f3 or f4) in group 1 leads to density f4 (f1 or f2, respectively)
#' # in the second one.
#'
#' set.seed(123)
#' p <- 100
#' n = m = 5
#' inds <- sample(1:4, p, replace = TRUE)
#' X <- matrix(rep(0, n * p), ncol = n)
#'
#' for (j in 1:p){
#'   if (inds[j] == 1){
#'     X[j, ] <- rnorm(n)}
#'   if (inds[j] == 2){
#'     X[j, ] <- rnorm(n, sd = sqrt(1/4))
#'   }
#'   if (inds[j] == 3){
#'     X[j, ]<-rnorm(n, mean = 2)
#'   }
#'   if (inds[j] == 4){
#'     X[j, ]<-rnorm(n, mean = 2, sd = sqrt(1/4))
#'   }
#' }
#'
#' rho <- 0.3
#'
#' ind <- sample(1:p, rho * p)
#' li <- length(ind)
#' indsy <- inds
#'
#' for (l in 1:li){
#'   if (indsy[ind[l]] == 1){indsy[ind[l]] = 3} else{
#'     if (indsy[ind[l]] == 2){indsy[ind[l]] = 4} else {
#'       if (indsy[ind[l]] == 3){indsy[ind[l]] = 1}
#'       else{indsy[ind[l]] = 2}}}}
#'
#' Y <- matrix(rep(0, m * p), ncol = m)
#'
#' for (j in 1:p){
#'   if (indsy[j] == 1){
#'     Y[j, ] <- rnorm(m)}
#'   if (indsy[j] == 2){
#'     Y[j, ] <- rnorm(m, sd = sqrt(1/4))
#'   }
#'   if (indsy[j] == 3){
#'     Y[j, ] <- rnorm(m, mean = 2)
#'   }
#'   if (indsy[j] == 4){
#'     Y[j, ]<- rnorm(m, mean = 2, sd = sqrt(1/4))
#'   }
#' }
#'
#'
#' # We can see which are the variables with different distributions in the two data sets.
#'
#' dif <- which(inds != indsy)
#'
#' # Cross table for (X,Y) density indexes:
#'
#' table(inds,indsy)
#'
#' # Our interest is to identify which variables have a different distribution in the two groups.
#' # Hence, since the differences between the distributions are in location, we applied
#' # Wilcoxon-Mann-Whitney test to verify for each variable the equality of distribution
#' # in the two groups.
#'
#' \donttest{
#'
#' library(exactRankTests)
#' library(coin)
#'
#' # We compute the P-values
#'
#' p <- nrow(X)
#' pv <- 1:p
#' for (i in 1:p){
#'   pv[i] <- wilcox.exact(X[i, ], Y[i, ])$p.value
#' }
#'
#' # When the sample size is small, in this case n=m=5, the distribution of
#' # the Wilcoxon's statistic is calibrated using an exact permutation test. Hence,
#' # the P-values are homogeneous discrete uniform distributed with support points
#' # of such distribution:
#'
#' ss <- c(1, 2, 4, 7, 12, 19, 28, 39, 53,69, 87, 106, 126) / 126
#'
#' # When the number of P-values is large enough "ss" is equal to:
#' sort(unique(pv))
#'
#' # For details about the Wilcoxon-Mann-Whitney test and its exact distribution
#' # see Section 9.2 of Gibbons and Chakraborti (1992).
#'
#' # We apply Chen method:
#'
#' R <- DQ(pv, ss = ss, method = "Chen")
#'
#' # The estimate of the proportion of null P-values:
#'
#' R$pi0
#'
#' # Summary of the vector of the estimated q-values:
#'
#' summary(R$q.values)
#'
#' # How many null hypotheses are rejected?
#' alpha <- 0.05
#' sum(R$q.values < alpha)
#'
#' # Which variables correspond to such null hypotheses?
#'
#' which(R$q.values < alpha)
#'
#' # Classification table (Decision at nominal level alpha vs. reality):
#'
#' table(R$q.values < alpha,inds != indsy)
#'
#' # The conclusion from the previous table is that Chen method reports
#' # 21 true positives and 9 false negatives.
#'
#'
#' # We can also apply Liang and SS methods as follows.
#' RLiang <- DQ(pv, ss = ss, method = "Liang")
#' RSS <- DQ(pv, ss = ss, method = "SS")
#'
#' # The next graphic help us to see that Liang method (for discrete P-values)
#' # is more powerful than SS method (only suitable for continuous P-values).
#'
#' plot(RSS$q.values,RLiang$q.values)
#' abline(a = 0, b = 1, col = 2, lty = 2)
#'
#' # We consider a poisson setting to show a case where the support of the discrete
#' # distribution of the P-values is infinte.
#' # We generate 100 values of a poisson distribution with event rate 10.
#' # Then we compute the probability that each of the values come from a
#' # a poisson distribution with event rate 10. This set of probabilities
#' # are considered as our set of P-values.
#'
#' p <- 100
#' N <- rpois(p, lambda = 10)
#' pv <- 1:p
#' for(i in 1:p){
#'   pv[i] <- ppois(N[i], lambda = 10)
#' }
#'
#' # It is well know that the support of a poisson distribution is infinite
#' # and equal to the natural numbers. Hence to know the support of the P-values
#' # defined above, we compute for 1,..., 50 their corresponding P-value.
#' # We only considered 50 values because for large values than 50 the P-value is 1 again.
#'
#' nn_0 <- 50
#' ss <- 1:(nn_0 + 1)
#' for (i in 1:(nn_0 + 1)){
#'  ss[i] <- ppois(i - 1, lambda = 10)
#' }
#'
#' # We eliminate repeated values.
#' ss <- unique(ss)
#' # For Chen method the relevant support points are only the values below tau[100] = 0.5.
#' # We define the support ss as such values. Then, we can apply Chen method. Of
#' # course s_inf = TRUE.
#'
#' indicator <- which(ss <= 0.5)
#' ssi <- ss[indicator]
#' R <- DQ(pv, ss = ssi, ss_inf = "TRUE", method = "Chen")
#' # For Liang method the relevant support values are also the values below 0.5, hence
#' # ss defined above is suitable.
#' R <- DQ(pv, ss = ssi, ss_inf = "TRUE", method = "Liang")
#'
#' # For Rand method the relevant support values are those ones with match with
#' # the P-values, and also the largest support points smaller than each one of such P-values.
#' # Hence, ss only includes the relevant points, as we can see below.
#' pv_unique <- unique(pv)
#' p_u <- length(pv_unique)
#' ind <- 1:p_u
#' for(i in 1:p_u){
#'   ind[i] <- which(ss == pv_unique[i])
#' }
#' p_u == length(ind)
#' ind_minus <- ind - 1
#' ind_final <- unique(sort(c(ind, ind_minus)))
#' ss <- ss[ind_final]
#' # Now, we can apply Rand method.
#' R <- DQ(pv, ss = ss, ss_inf = "TRUE", method = "Rand")
#'
#'}
#'
#'
#' @importFrom stats predict quantile runif smooth.spline density dnorm qnorm
#' @export
DQ <- function(pv, ss = NULL, ss_inf = FALSE, method = c("ST", "SS", "Liang", "Chen", "Rand"), lambda = seq(0.05, 0.95, 0.05)){

  # cat("Call:", "\n")
  # print(match.call())

  match.arg(method)

  if (missing(method)) {
    method <- "ST"
    cat("'ST' method used by default\n")
  }

  if(method == "ST" || method == "SS"){

    ss <- NULL
  }

 if(ss_inf == FALSE) {

  # Liang method
  if(method == "Liang"){

    p.grid <- ss
    NN <- length(ss)
    p.count.mat <- get.tbl.full(pv, NN, ss)
    pi0_Liang <- as.numeric(est.pi0.disc(p.count.mat, p.grid)$pi0)
    pi0_Liang <- max(pi0_Liang, 0)
    pi0_Liang <- min(pi0_Liang, 1)
    q_values_Liang <- q.value_pi0(pv, pi0_Liang)

    return(list(pi0 = pi0_Liang, q.values = q_values_Liang))
  }


  # Chen method
  if(method == "Chen"){

    pi0_Chen <- pi0_G(pv, ss)
    pi0_Chen <- max(pi0_Chen, 0)
    pi0_Chen <- min(pi0_Chen, 1)
    q_values_chen <- q.value_pi0(pv, pi0_Chen)
    return(list(pi0 = pi0_Chen, q.values = q_values_chen))
  }


  # Randomized method
  if(method == "Rand"){

    pi0_Rd <- Rp(pv, ss)
    pi0_Rd <- max(pi0_Rd, 0)
    pi0_Rd <- min(pi0_Rd, 1)
    q_values_Rd <- q.value_pi0(pv, pi0_Rd)
    return(list(pi0 = pi0_Rd, q.values = q_values_Rd))
  }


  ### Q-value method for continuous p-values

  if(method == "ST"){

    l_lambda <- length(lambda)

    pi0<-1:l_lambda
    for(i in 1:l_lambda){
      pi0[i] <- sum(pv >= lambda[i]) / ((1 - lambda[i]) * length(pv))
    }

    spi0 <- smooth.spline(lambda, pi0, df = 3)
    pi0Smooth <- predict(spi0, x = lambda)$y
    pi0_ST <- min(pi0Smooth[l_lambda], 1)
    pi0_ST <- max(pi0_ST, 0)

    q.values_ST <- q.value_pi0(pv, pi0_ST)

    return(list(pi0 = pi0_ST, q.values = q.values_ST))
  }

  if(method == "SS"){

    lambda <- 0.5
    pi0_SS <- sum(pv >= lambda) / ((1 - lambda) * length(pv))
    pi0_SS <- min(pi0_SS, 1)
    pi0_SS <- max(pi0_SS, 0)
    q.values_SS <- q.value_pi0(pv, pi0_SS)
    return(list(pi0 = pi0_SS, q.values = q.values_SS))
  }

 } else {

   ##### Liang method
   if(method == "Liang"){
    p.grid<-ss
    NN <- length(ss)
    p.count.mat <- get.tbl.full(pv, NN, ss)
    pi0_Liang <- est.pi0.disc(p.count.mat, p.grid)$pi0 # DEVOLVE A FUNCIÓN pi0 método Liang
    q_values_Liang <- q.value_pi0(pv, pi0_Liang) # DEVOLVE A FUNCIÓN q.values método Liang
    return(list(pi0 = pi0_Liang, q.values = q_values_Liang))
   }


   ### Chen method (the same as finite support)
   if(method == "Chen"){
    pv <- unlist(pv)
    pi0_Chen <- pi0_G(pv,ss) #DEVOLVE A FUNCIÓN pi0 método Chen
    q_values_chen <- q.value_pi0(pv, pi0_Chen) #DEVOLVE A FUNCIÓN q.values método chen
    return(list(pi0 = pi0_Chen, q.values = q_values_chen))
   }

   ### Randomized method
   if(method == "Rand"){
    pi0_Rd <- Rp(pv, ss) #DEVOLVE A FUNCIÓN pi0 método Rand
    q_values_Rd <- q.value_pi0(pv, pi0_Rd) #DEVOLVE A FUNCIÓN q.values método Rand
    return(list(pi0 = pi0_Rd, q.values = q_values_Rd))
   }

 }



}
