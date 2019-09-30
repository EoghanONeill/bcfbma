
#' @title Bayesian Causal Forest Using Bayesian Model Averaging (BCF-BMA)
#' 
#' @description This is an implementation of Bayesian causal forests (Hahn et al. 2018) using Bayesian Model Averaging, following the approach used in BART-BMA (Hernandez et al. 2018). 
#' The outcome is modelled as mu(x)+Z*tau(x), where mu(x) and tau(x) are sums-of-trees and Z is a treatment indicator variable. Therefore tau(x) gives the Individual Treatment Effect (ITE) estimates.
#' @param x.train Training data covariate matrix excluding the treatment and propensity scores.
#' @param y.train Training data outcome vector.
#' @param z Training data treatment vector. This should be a binary vector, equal to one for treated individuals.
#' @param pihat Training data propensity score estimates. Each column should be a vector of propensity score estimates.
#' @param a_mu This is a parameter that influences the variance of the terminal node values for trees that are not interacted with the treatment.
#' @param a_tau This is a parameter that influences the variance of the terminal node values for trees that are interacted with the treatment.
#' @param nu This is a hyperparameter in the distribution of the variance of the error term. THe inverse of the variance is distributed as Gamma (nu/2, nu*lambda/2).
#' @param sigquant ??
#' @param c This determines the size of Occam's Window
#' @param pen_mu This is a parameter used by the Pruned Exact Linear Time Algorithm when finding changepoints for mu trees.
#' @param pen_tau This is a parameter used by the Pruned Exact Linear Time Algorithm when finding changepoints for tau trees.
#' @param num_cp_mu This is a number between 0 and 100 that determines the proportion of changepoints proposed by the changepoint detection algorithm to keep when growing mu trees.
#' @param num_cp_tau This is a number between 0 and 100 that determines the proportion of changepoints proposed by the changepoint detection algorithm to keep when growing tau trees.
#' @param x.test Test data covariate matrix excluding the treatment and propensity scores.
#' @param test_z Test data treatment vector. This should be a binary vector, equal to one for treated individuals.
#' @param test_pihat Test data propensity score estimates. Each column should be a vector of propensity score estimates.
#' @param ntree_control Number of trees that are not interacteds with the treatment.
#' @param ntree_moderate Number of trees that are interacted with the treatment.
#' @param alpha_mu Parameter in prior probability of tree node splitting.
#' @param alpha_tau Parameter in prior probability of tree node splitting.
#' @param beta_mu Parameter in prior probability of tree node splitting.
#' @param beta_tau Parameter in prior probability of tree node splitting.
#' @param split_rule_node Binary variable. If equals 1, then find a new set of potential splitting points via a changepoint algorithm after adding each split to a tree. If equals zero, use the same set of potential split points for all splits in a tree.
#' @param gridpoint Binary variable. If equals 1, then a grid search changepoint detection algorithm will be used. If equals 0, then the Pruned Exact Linear Time (PELT) changepoint detection algorithm will be used (Killick et al. 2012). 
#' @param maxOWsize Maximum number of models to keep in Occam's window
#' @param num_splits_mu Maximum number of splits in a mu tree
#' @param num_splits_tau Maximum number of splits in a tau tree
#' @param gridsize_mu This integer determines the size of the grid across which to search if gridpoint=1 when constructing mu trees.
#' @param gridsize_tau This integer determines the size of the grid across which to search if gridpoint=1 when constructing tau trees.
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control" or "both" are HIGHLY recommended with observational data.
#' @param zero_split Binary variable. If equals 1, then zero split trees can be included in a sum-of-trees model. If equals zero, then only trees with at least one split can be included in a sum-of-trees model.
#' @param only_max_num_trees Binary variable. If equals 1, then only sum-of-trees models containing the maximum number of trees, num_rounds, are selected. If equals 0, then sum-of-trees models containing less than num_rounds trees can be selected. The default is only_max_num_trees=1.
#' @param mu_or_tau_each_round Binary variable. If equals 1, then a mu tree or a tau tree is added in each round. If equals 0, then a mu tree is added, followed by a tau tree, followed by a mu tree, and so on.
#' @param separate_tree_numbers Binary variable (irrelevant if mu_or_tau_each_round not equal to 1). If equals 1, and mu_or_tau_each_round equals 1, then num_splits_mu and num_splits_tau are the maximum numbers of mu and tau trees in the model. If equals zero, and mu_or_tau_each_round equals 1, then num_splits_mu + num_splits_tau is the maximum total number of trees, but there are not separate limits to the number of mu and tau trees.
#' @param min_num_obs_for_mu_split This integer determines the minimum number of observations in a (parent) mu tree node for the algorithm to consider potential splits of the node.
#' @param min_num_obs_after_mu_split This integer determines the minimum number of observations in a (mu tree) child node resulting from a split in order for a split to occur. If the left or right child node has less than this number of observations, then the split can not occur.
#' @param min_num_obs_for_tau_split This integer determines the minimum number of treated observations in a (parent) tau tree node for the algorithm to consider potential splits of the node.
#' @param min_num_obs_after_tau_split This integer determines the minimum number of treated observations in a (tau tree) child node resulting from a split in order for a split to occur. If the left or right child node has less than this number of observations, then the split can not occur.
#' @param exact_residuals Binary variable. If equal to 1, then trees are added to sum-of-tree models within each round of the algorithm by detecting changepoints in the exact residuals. If equals zero, then changepoints are detected in residuals that are constructed from approximate predictions.
#' @param transform_resids Binary variable. If equal to 1, then a Horvitz-Thompson transformation is applied to residuals before finding chanegpoints for building tau trees.

#' @export 
#' @return The following objects are returned by bcfbma:
#' \item{fitted.values_outcome}{The vector of predictions of the outcome for all training observations.} 
#' \item{fitted.values_mu}{The vector of fiited values of mu(x) for all training observations.}
#' \item{fitted.values_tau}{The vector of fiited values of tau(x) for all training observations. These are the in-sample ITE estimates.}
#' \item{sumoftrees_mu}{This is a list of lists of matrices. The outer list corresponds to a list of sum-of-tree models, and each element of the outer list is a list of matrices describing the structure of the mu(x) trees within a sum-of-tree model. See details.} 
#' \item{sumoftrees_tau}{This is a list of lists of matrices. The outer list corresponds to a list of sum-of-tree models, and each element of the outer list is a list of matrices describing the structure of the tau(x) trees within a sum-of-tree model. See details.} 
#' \item{obs_to_termNodesMatrix_mu}{This is a list of lists of matrices. The outer list corresponds to a list of sum-of-tree models, and each element of the outer list is a list of matrices describing to which node each of the observations is allocated to at all depths of each mu(x) trees within a sum-of-tree model. See details.} 
#' \item{obs_to_termNodesMatrix_tau}{This is a list of lists of matrices. The outer list corresponds to a list of sum-of-tree models, and each element of the outer list is a list of matrices describing to which node each of the observations is allocated to at all depths of each tau(x) trees within a sum-of-tree model. See details.}
#' \item{bic}{This is a vector of BICs for each sum-of-tree model.}
#' \item{sum_residuals_mu}{A list of lists (if more than one mu(x) tree) of vectors of partial residuals for each tree in mu(x) in each model.}
#' \item{sum_residuals_tau}{A list of lists (if more than one tau(x) tree) of vectors of partial residuals for each tree in tau(x) in each model.}
#' \item{numvars}{This is the total number of variables in the input training data matrix, excluding the pihat matrix.} 
#' \item{call}{match.call returns a call in which all of the specified arguments are specified by their full names.} 
#' \item{y_minmax}{Range of the input training data outcome vector.} 
#' \item{response}{Input taining data outcome vector.}
#' \item{nrowTrain}{number of observations in the input training data.} 
#' \item{sigma}{sd(y.train)/(max(y.train)-min(y.train))} 
#' \item{a_mu}{input parameter} 
#' \item{a_tau}{input parameter} 
#' \item{nu}{input parameter} 
#' \item{lambda}{parameter determined by the inputs sigma, sigquant, and nu} 
#' \item{numPSmethods}{Number of columns of the input matrix pihat. This should be the number of different propoensity score estimates used.} 
#' \item{include_pi2}{Equals 0 if input parameter include_pi is control, 1 if moderate, 2 if both, 4 if none.}
#' \item{z}{Input training data treatment vector. This should be a binary vector, equal to one for treated individuals.}


bcfBMA<-function(x,...)UseMethod("bcfBMA")

bcfBMA.default<-function(x.train,y.train,z,pihat,
                         a_mu=3,a_tau=3,nu=3,sigquant=0.9,c=1000,
                          pen_mu=12,pen_tau=12,num_cp_mu=20,num_cp_tau=20,
                          x.test=matrix(0.0,0,0),test_z = numeric(),test_pihat = matrix(0.0,0,0),
                          ntree_control=5,ntree_moderate=5,
                          alpha_mu=0.95,alpha_tau=0.95,beta_mu=1,beta_tau=1,split_rule_node=0,
                          gridpoint=0,maxOWsize=100, num_splits_mu =5, num_splits_tau =5, gridsize_mu=10, gridsize_tau=10,
                          include_pi= "control", zero_split=1, only_max_num_trees=1, mu_or_tau_each_round=1,separate_tree_numbers=1,
                         min_num_obs_for_mu_split=2, min_num_obs_after_mu_split=2,
                         min_num_obs_for_tau_split=2, min_num_obs_after_tau_split=2,
                         exact_residuals=1, 
                         spike_tree=0,s_t_hyperprior=1, 
                         p_s_t_mu=0.5, a_s_t_mu=1,b_s_t_mu=3,
                         p_s_t_tau=0.5, a_s_t_tau=1,b_s_t_tau=3,
                         lambda_poisson_mu=10,lambda_poisson_tau=10,
                         transform_resids=0){
  binary=FALSE
  start_mean=0
  start_sd=1
  mu_mu=0
  mu_tau=0
  sigma_mu_mu=0; 
  sigma_mu_tau=0; 
  
  #if(length(y.train) > 3 + ncol(cbind(pihat, x.train)) ){
    #df = data.frame((x.train),y.train)
    #lmf = lm(y.train~z+as.matrix(cbind(pihat, x.train)))
    #sigest = summary(lmf)$sigma
  #}else{
    sigest = sd(y.train)
  #}
  
  sigma=sigest/(max(y.train)-min(y.train))
  
  #sigma=sd(y.train)/(max(y.train)-min(y.train))
  qchi = qchisq(1.0-sigquant,nu,1,0);
  lambda = (sigma*sigma*qchi)/nu;
  include_pi2=-1
  if(include_pi=="control") {
    include_pi2 = 0
  }
  if(include_pi=="moderate") {
    include_pi2 = 1
  }
  if(include_pi=="both") {
    include_pi2 = 2
  }
  if(include_pi=="none") {
    include_pi2 = 4
  }
  if(include_pi2==-1) stop('include_pi must be equal to control, moderate, both, or none.')

  
  #PROBABLY NEED TO DO SOMETHING SIMILAR FOR z and test_z with binary_z and binary_test_z
  if(is.factor(y.train)) {
    # if(length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
    binary = TRUE
    #  y.train = as.numeric(y.train)-1
    stop("Response must be a numeric vector")
  } else {
    if((length(unique(y.train)) == 2) & (max(y.train) == 1) & (min(y.train) == 0)) {
      cat('NOTE: assumming numeric response is binary\n')
      binary = TRUE
      #stop("Response must be a numeric vector")
    }
  }
  
  if(!all(z %in% c(0,1))) stop("z values must be 0 or 1")
  if(!all(test_z %in% c(0,1))) stop("test_z values must be 0 or 1")
  if(gridsize_mu<1) stop("gridsize_mu must be a positive integer")
  if(gridsize_tau<1) stop("gridsize_tau must be a positive integer")
  
  if(is.vector(x.train) | is.factor(x.train)| is.data.frame(x.train)) x.train = as.matrix(x.train)
  if(is.vector(x.test) | is.factor(x.test)| is.data.frame(x.test)) x.test = as.matrix(x.test)
  if(is.vector(pihat) | is.factor(pihat)| is.data.frame(pihat)) pihat = as.matrix(pihat)
  if(is.vector(test_pihat) | is.factor(test_pihat)| is.data.frame(test_pihat)) pihat = as.matrix(test_pihat)
  
  if(is.matrix(x.train)) {
    if(nrow(x.test)>0) {
      if(!is.matrix(x.test)) stop('x.test must be a matrix')
    } 
  }
  if(is.matrix(x.train)) {
    if(nrow(pihat)>0) {
      if(!is.matrix(pihat)) stop('x.test must be a matrix')
    } 
  }
  if(is.matrix(x.train)) {
    if(nrow(test_pihat)>0) {
      if(!is.matrix(test_pihat)) stop('x.test must be a matrix')
    } 
  }
  #check input arguments:
  # if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
  #if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
  if((!is.matrix(x.train))) stop("argument x.train must be a double matrix")
  if((!is.matrix(x.test)) ) stop("argument x.test must be a double matrix")
  if((!is.matrix(pihat)) ) stop("argument x.test must be a double matrix")
  if((!is.matrix(test_pihat)) ) stop("argument x.test must be a double matrix")
  
  #PROBABLY NEED TO DO SOMETHING SIMILAR FOR z and test_z with binary_z and binary_test_z
  if(!binary) {
    if((!is.vector(y.train))) stop("argument y.train must be a double vector")
  }

  if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
  if(nrow(pihat) != length(y.train)) stop("number of rows in pihat must equal length of y.train")
  if(((nrow(pihat) >0) && (nrow(x.test) >0))&&(nrow(x.test) != nrow(test_pihat))) stop("number of rows in test_pihat must equal number of rows in x.test ")
  
  if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
  if(((nrow(pihat) >0) && (nrow(test_pihat) >0)) && (ncol(pihat)!=ncol(test_pihat))) stop("input test_pihat must have the same number of columns as pihat")
  #if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
  #if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
  #if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
  if(c<1)stop("Value of Occam's Window has to be greater than 0."); 
  if(num_cp_mu<0 || num_cp_mu>100)stop("Value of num_cp_mu should be a value between 1 and 100."); 
  if(num_cp_tau<0 || num_cp_tau>100)stop("Value of num_cp_tau should be a value between 1 and 100."); 
  
  
  if(mu_or_tau_each_round==1){
  bcfBMA_call=BCF_BMA_sumLikelihood_add_mu_or_tau(spike_tree,s_t_hyperprior, 
                                                  p_s_t_mu, a_s_t_mu,b_s_t_mu,
                                                  p_s_t_tau, a_s_t_tau,b_s_t_tau,
                                                  lambda_poisson_mu,lambda_poisson_tau,
                                                  x.train,y.train,z,pihat,
                                                  a_mu,a_tau,mu_mu,mu_tau,
                                                  nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                                  pen_mu,pen_tau,num_cp_mu,num_cp_tau,
                                                  x.test,test_z,test_pihat,ntree_control,ntree_moderate,
                                                  alpha_mu,alpha_tau,beta_mu,beta_tau,
                                                  split_rule_node,gridpoint,maxOWsize,
                                                  num_splits_mu,num_splits_tau,gridsize_mu, gridsize_tau,
                                                  include_pi2,zero_split,only_max_num_trees,separate_tree_numbers,
                                                  min_num_obs_for_mu_split, min_num_obs_after_mu_split,
                                                  min_num_obs_for_tau_split, min_num_obs_after_tau_split,
                                                  exact_residuals,
                                                  transform_resids)
  }else{
  bcfBMA_call=BCF_BMA_sumLikelihood(spike_tree,s_t_hyperprior, 
                                    p_s_t_mu, a_s_t_mu,b_s_t_mu,
                                    p_s_t_tau, a_s_t_tau,b_s_t_tau,
                                    lambda_poisson_mu,lambda_poisson_tau,
                                    x.train,y.train,z,pihat,
                                    a_mu,a_tau,mu_mu,mu_tau,
                                    nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                    pen_mu,pen_tau,num_cp_mu,num_cp_tau,
                                    x.test,test_z,test_pihat,ntree_control,ntree_moderate,
                                    alpha_mu,alpha_tau,beta_mu,beta_tau,
                                    split_rule_node,gridpoint,maxOWsize,
                                    num_splits_mu,num_splits_tau,gridsize_mu, gridsize_tau,
                                    include_pi2,zero_split,only_max_num_trees,
                                    min_num_obs_for_mu_split, min_num_obs_after_mu_split,
                                    min_num_obs_for_tau_split, min_num_obs_after_tau_split,
                                    exact_residuals,
                                    transform_resids)
  }
  
  
  if(length(bcfBMA_call)==13){
    #length of bcfBMA_call is 11 if test data was included in the call
    names(bcfBMA_call)<-c("fitted.values_outcome","fitted.values_mu","fitted.values_tau",
                          "sumoftrees_mu","sumoftrees_tau",
                          "obs_to_termNodesMatrix_mu","obs_to_termNodesMatrix_tau",
                          "bic",
                          "test.preds_outcome","test.preds_mu","test.preds_tau","sum_residuals_mu","sum_residuals_tau")
    bcfBMA_call[[12]]<-bcfBMA_call[[12]]#[[length(bcfBMA_call[[12]])]]
    bcfBMA_call[[13]]<-bcfBMA_call[[13]]#[[length(bcfBMA_call[[13]])]]
    bcfBMA_call$test_data<-x.test
    bcfBMA_call$test_pihat <- test_pihat
    bcfBMA_call$test_z <- test_z
  }else{
    names(bcfBMA_call)<-c("fitted.values_outcome","fitted.values_mu","fitted.values_tau",
                          "sumoftrees_mu","sumoftrees_tau",
                          "obs_to_termNodesMatrix_mu","obs_to_termNodesMatrix_tau",
                          "bic","sum_residuals_mu","sum_residuals_tau")
    bcfBMA_call[[9]]<-bcfBMA_call[[9]]#[[length(bcfBMA_call[[9]])]]
    bcfBMA_call[[10]]<-bcfBMA_call[[10]]#[[length(bcfBMA_call[[10]])]]
  }
  
  bcfBMA_call$numvars<-ncol(x.train) #ncol(training)
  bcfBMA_call$call<-match.call()
  bcfBMA_call[[4]]<-bcfBMA_call[[4]]#[[length(bcfBMA_call[[4]])]]
  bcfBMA_call[[5]]<-bcfBMA_call[[5]]#[[length(bcfBMA_call[[5]])]]
  bcfBMA_call[[6]]<-bcfBMA_call[[6]]#[[length(bcfBMA_call[[6]])]]
  bcfBMA_call[[7]]<-bcfBMA_call[[7]]#[[length(bcfBMA_call[[7]])]]
  bcfBMA_call$y_minmax<-range(y.train)
  bcfBMA_call$response<-y.train
  bcfBMA_call$nrowTrain<-nrow(x.train)
  bcfBMA_call$sigma<-sigma
  bcfBMA_call$a_mu<-a_mu
  bcfBMA_call$a_tau<-a_tau
  bcfBMA_call$nu<-nu
  bcfBMA_call$lambda<-lambda
  bcfBMA_call$numPSmethods<-ncol(pihat)
  bcfBMA_call$include_pi2<-include_pi2
  bcfBMA_call$z<-z
  bcfBMA_call$num_propscores <- ncol(pihat)
  bcfBMA_call$nrowtest<-nrow(x.test)
  class(bcfBMA_call)<-"bcfBMA"
  bcfBMA_call
}