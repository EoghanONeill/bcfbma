
#' @title Bayesian Causal Forest Using Bayesian Model Averaging
#' 
#' @description Insert description
#' @param x.train Training data covariate matrix excluding the treatment and propensity scores.
#' @param y.train Training data outcome vector.
#' @param z Training data treatment vector. This should be a binary vector, equal to one for treated individuals.
#' @param pihat Training data propensity score estimates. Each column should be a vector of propensity score estimates.
#' @param a_mu This is a parameter that influences the variance of the terminal node values for trees that are not interacted with the treatment.
#' @param a_tau This is a parameter that influences the variance of the terminal node values for trees that are interacted with the treatment.
#' @param nu ??
#' @param sigquant ??
#' @param c ? This determines the size of Occam's Window
#' @param pen_mu ??
#' @param pen_tau ??
#' @param num_cp_mu ??
#' @param num_cp_tau ??
#' @param x.test Test data covariate matrix excluding the treatment and propensity scores.
#' @param test_z Test data treatment vector. This should be a binary vector, equal to one for treated individuals.
#' @param test_pihat Test data propensity score estimates. Each column should be a vector of propensity score estimates.
#' @param ntree_control Number of trees that are not interacteds with the treatment.
#' @param ntree_moderate Number of trees that are interacted with the treatment.
#' @param alpha_mu Parameter in prior probability of tree node splitting.
#' @param alpha_tau Parameter in prior probability of tree node splitting.
#' @param beta_mu Parameter in prior probability of tree node splitting.
#' @param beta_tau Parameter in prior probability of tree node splitting.
#' @param split_rule_node ??
#' @param gridpoint ??
#' @param maxOWsize Maximum number of models to keep in Occam's window
#' @export
#' @return Include lots of details here.


bcfBMA<-function(x,...)UseMethod("bcfBMA")

bcfBMA.default<-function(x.train,y.train,z,pihat,
                         a_mu=3,a_tau=3,nu=3,sigquant=0.9,c=1000,
                          pen_mu=12,pen_tau=12,num_cp_mu=20,num_cp_tau=20,
                          x.test=matrix(0.0,0,0),test_z = numeric(),test_pihat = matrix(0.0,0,0),
                          ntree_control=5,ntree_moderate=5,
                          alpha_mu=0.95,alpha_tau=0.95,beta_mu=1,beta_tau=1,split_rule_node=0,
                          gridpoint=0,maxOWsize=100){
  binary=FALSE
  start_mean=0
  start_sd=1
  mu_mu=0
  mu_tau=0
  sigma_mu_mu=0; 
  sigma_mu_tau=0; 
  sigma=sd(y.train)/(max(y.train)-min(y.train))
  qchi = qchisq(1.0-sigquant,nu,1,0);
  lambda = (sigma*sigma*qchi)/nu;
  
  
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
  
  
  if(is.vector(x.train) | is.factor(x.train)| is.data.frame(x.train)) x.train = as.matrix(x.train)
  if(is.vector(x.test) | is.factor(x.test)| is.data.frame(x.test)) x.test = as.matrix(x.test)
  if(is.vector(pihat) | is.factor(pihat)| is.data.frame(pihat)) pihat = as.matrix(pihat)
  if(is.vector(test_pihat) | is.factor(test_pihat)| is.data.frame(test_pihat)) pihat = as.matrix(test_pihat)
  
  if(is.matrix(x.train)) {
    if(nrow(x.test)) {
      if(!is.matrix(x.test)) stop('x.test must be a matrix')
    } 
  }
  if(is.matrix(x.train)) {
    if(nrow(pihat)) {
      if(!is.matrix(pihat)) stop('x.test must be a matrix')
    } 
  }
  if(is.matrix(x.train)) {
    if(nrow(test_pihat)) {
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
  
  bcfBMA_call=BCF_BMA_sumLikelihood(x.train,y.train,z,pihat,
                                    a_mu,a_tau,mu_mu,mu_tau,
                                    nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                    pen_mu,pen_tau,num_cp_mu,num_cp_tau,
                                    x.test,test_z,test_pihat,ntree_control,ntree_moderate,
                                    alpha_mu,beta_mu,alpha_tau,beta_tau,
                                    split_rule_node,gridpoint,maxOWsize)
  
  if(length(bcfBMA_call)==11){
    #length of bcfBMA_call is 11 if test data was included in the call
    names(bcfBMA_call)<-c("fitted.values_outcome","fitted.values_mu","fitted.values_tau",
                          "sumoftrees_mu","sumoftrees_tau",
                          "obs_to_termNodesMatrix_mu","obs_to_termNodesMatrix_tau",
                          "bic",
                          "test.preds_outcome","test.preds_mu","test.preds_tau")
    bcfBMA_call$test_data<-x.test
  }else{
    names(bcfBMA_call)<-c("fitted.values_outcome","fitted.values_mu","fitted.values_tau",
                          "sumoftrees_mu","sumoftrees_tau",
                          "obs_to_termNodesMatrix_mu","obs_to_termNodesMatrix_tau",
                          "bic")
  }
  
  bcfBMA_call$numvars<-ncol(x.train) #ncol(training)
  bcfBMA_call$call<-match.call()
  bcfBMA_call[[4]]<-bcfBMA_call[[4]][[length(bcfBMA_call[[4]])]]
  bcfBMA_call[[5]]<-bcfBMA_call[[5]][[length(bcfBMA_call[[5]])]]
  bcfBMA_call[[6]]<-bcfBMA_call[[6]][[length(bcfBMA_call[[6]])]]
  bcfBMA_call[[7]]<-bcfBMA_call[[7]][[length(bcfBMA_call[[7]])]]
  bcfBMA_call$y_minmax<-range(y.train)
  bcfBMA_call$response<-y.train
  bcfBMA_call$nrowTrain<-nrow(x.train)
  bcfBMA_call$sigma<-sigma
  bcfBMA_call$a_mu<-a_mu
  bcfBMA_call$a_tau<-a_tau
  bcfBMA_call$nu<-nu
  bcfBMA_call$lambda<-lambda
  bcfBMA_call$numPSmethods<-ncol(pihat)
  
  class(bcfBMA_call)<-"bcfBMA"
  bcfBMA_call
}