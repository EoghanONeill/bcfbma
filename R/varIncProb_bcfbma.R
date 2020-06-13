#' @title Variable inclusion probabilities as defined by Linero (2018)
#' 
#' @description This measure defines the posterior inclusion probability of a variable as the model-probability weighted sum of indicator variables for whether the variable was used in any splitting rules in any of the trees in the sum-of-tree model.
#' @param object A bartBMA object obtained using the barBMA function.
#' @export 
#' @return A list containing: Vectors of posterior inclusion probabilities for 1. mu (control trees), i.e. the prognostic effect and 2. tau (treatment effect moderator trees), i.e. the treatment effect. The variables are ordered in the same order that they occur in columns of the input covariate matrix used to obtain the input bartBMA object.
varIncProb_bcfbma <-function(object,...){
  #object will be bcfBMA object.
  
  num_to_add_mu = 0
  num_to_add_tau = 0
  
  if(object$include_pi2 == 0) {
    #include_pi2 = 0
    # control
    num_to_add_mu = num_to_add_mu + object$numPSmethods
  }
  if(object$include_pi2 == 1) {
    #include_pi2 = 1
    #moderate
    
    num_to_add_tau = num_to_add_tau + object$numPSmethods
    
  }
  if(object$include_pi2 == 2) {
    #include_pi2 = 2
    # both
    
    num_to_add_mu = num_to_add_mu + object$numPSmethods
    num_to_add_tau = num_to_add_tau + object$numPSmethods
    
  }
  if(object$include_pi2 == 4) {
    #include_pi2 = 4
    # none
    
  }
  
  imp_vars2_mu=get_weighted_var_imp(num_vars=object$numvars + num_to_add_mu ,
                                    BIC=object$bic,
                                    sum_trees=object$sumoftrees_mu)
  res_mu<-apply((imp_vars2_mu[[3]]>0)*imp_vars2_mu[[1]],2,sum)
  
  imp_vars2_tau=get_weighted_var_imp(num_vars=object$numvars + num_to_add_tau ,
                                     BIC=object$bic,
                                     sum_trees=object$sumoftrees_tau)
  
  res_tau<-apply((imp_vars2_tau[[3]]>0)*imp_vars2_tau[[1]],2,sum)
  
  ret<-list()
  length(ret)<-2
  ret[[1]] <- res_mu
  ret[[2]] <- res_tau
  
  class(ret)<-"varPIP_list.bcfBMA"  
  ret
  
  
  
}