#' @title Variable inclusion probabilities as defined by Linero (2018)
#' 
#' @description This measure defines the posterior inclusion probability of a variable as the model-probability weighted sum of indicator variables for whether the variable was used in any splitting rules in any of the trees in the sum-of-tree model.
#' @param object A bartBMA object obtained using the barBMA function.
#' @export 
#' @return A vector of posterior inclusion probabilities. The variables are ordered in the same order that they occur in columns of the input covariate matrix used to obtain the input bartBMA object.
varIncProb_bcfbma <-function(object,...){
  #object will be bartBMA object.
  imp_vars2_mu=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees_mu)
  res_mu<-apply((imp_vars2_mu[[3]]>0)*imp_vars2_mu[[1]],2,sum)
  
  imp_vars2_tau=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees_tau)
  res_tau<-apply((imp_vars2_tau[[3]]>0)*imp_vars2_tau[[1]],2,sum)
  
  ret<-list()
  length(ret)<-2
  ret[[1]] <- res_mu
  ret[[2]] <- res_tau
  
  class(ret)<-"varPIP_list.bcfBMA"  
  ret
  
  
  
}