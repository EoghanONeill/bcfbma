# bcfbma

<!-- badges: start -->
<!-- badges: end -->

The goal of bcfbma is to provide an implementation of Bayesian Causal Forests using Bayesian Model Averaging (BCF-BMA). This combines BART-BMA (Hernandez et al. 2018) and BCF (Hahn et al. 2017).

Hernández, B., Raftery, A. E., Pennington, S. R., & Parnell, A. C. (2018). Bayesian additive regression trees using Bayesian model averaging. Statistics and computing, 28(4), 869-890.

Hahn, P. R., Murray, J. S., & Carvalho, C. (2017). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects. arXiv preprint arXiv:1706.09523.


## Installation

``` r
# install.packages("devtools")
devtools::install_github("EoghanONeill/bcfbma")
```

## Example

This is a basic example which shows you how to solve a common problem:

```{r example}
library(bcfbma)
## basic example code

# data generating process
p = 3 #two control variables and one moderator
n = 250
#
set.seed(1)

x = matrix(rnorm(n*p), nrow=n)
# create targeted selection
q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
# generate treatment variable
#pi = rep(0.5,n)
#pi=runif(n)/4+0.5
pi= pnorm(q)
z = rbinom(n,1,pi)
# tau is the true (homogeneous) treatment effect
tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
# generate the response using q, tau and z
mu = (q + tau*z)
# set the noise level relative to the expected mean function of Y
sigma = diff(range(q + tau*pi))/8
# draw the response variabligma = diff(range(q + tae with additive error
#y = mu + rnorm(n)#sigma*rnorm(n)
y = mu + sigma*rnorm(n)
# If you didn't know pi, you would estimate it here
pihat = pi#pnorm(q)
class(pihat)
pihat = as.matrix(pihat)


bcfBMA_rfunc_example<-bcfBMA(x,y,z,pihat,
                             a_mu=3,a_tau=1.5,nu=3,sigquant=0.9,c=5000,
                             pen_mu=12,pen_tau=12,num_cp_mu=20,num_cp_tau=20,
                             x.test=matrix(0.0,0,0),test_z = numeric(),test_pihat = matrix(0.0,0,0),
                             ntree_control=5,ntree_moderate=5,
                             alpha_mu=0.95,alpha_tau=0.25,beta_mu=1,beta_tau=3,split_rule_node=0,
                             gridpoint=1,maxOWsize=100,num_splits_mu =25, num_splits_tau =25,
                             gridsize_mu=10, gridsize_tau=10, include_pi= "control",
                             zero_split = 1, only_max_num_trees = 1,mu_or_tau_each_round = 1,
                             min_num_obs_after_mu_split = 5, min_num_obs_after_tau_split = 5,
                             exact_residuals = 1,spike_tree = 0,
                             transform_resids = 0)


#tauhat_bcfbma1  <-  bcfBMA_rfunc_example$fitted.values_tau

tauhat_bcfbma_linalg<-preds_bcfbma_lin_alg(bcfBMA_rfunc_example,num_iter=2000,newdata=NULL)

plot(tau, tauhat_bcfbma_linalg); abline(0,1)


```
