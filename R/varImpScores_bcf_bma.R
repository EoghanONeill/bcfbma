#' @title Variable importances as defined by Hernandez et al. (2018)
#' 
#' @description This measure defines the importance of a variable as the model-probability weighted sum of the number of splits on the variable of interest, divided by the sum over all variables of such weighted counts of splits.
#' @param object A bartBMA object obtained using the barBMA function.
#' @export 
#' @return A vector of variable importances. The variables are ordered in the same order that they occur in columns of the input covariate matrix used to obtain the input bartBMA object.
#' @examples 
#' #set the seed
#' set.seed(100)
#' #simulate some data
#' N <- 100
#' p<- 100
#' epsilon <- rnorm(N)
#' xcov <- matrix(runif(N*p), nrow=N)
#' y <- sin(pi*xcov[,1]*xcov[,2]) + 20*(xcov[,3]-0.5)^2+10*xcov[,4]+5*xcov[,5]+epsilon
#' epsilontest <- rnorm(N)
#' xcovtest <- matrix(runif(N*p), nrow=N)
#' ytest <- sin(pi*xcovtest[,1]*xcovtest[,2]) + 20*(xcovtest[,3]-0.5)^2+10*xcovtest[,4]+
#'   5*xcovtest[,5]+epsilontest
#' 
#' #Train the object 
#' bart_bma_example <- bartBMA(x.train = xcov,y.train=y,x.test=xcovtest,zero_split = 1, 
#'                             only_max_num_trees = 1,split_rule_node = 0)
#' #Obtain the variable importances
#' varImpScores(bart_bma_example)

varImpScores_bcf_bma<-function(object){
  #object will be bartBMA object.
  imp_vars2_mu=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees_mu)
  res_mu <-apply(imp_vars2_mu[[4]],2,sum)
  #create varImpPlot command
  vIP_mu <-rep(NA,length(res_mu))
  total_weighted_var_counts_mu <-sum(res_mu)
  #now get variable inclusion probabilities
  vIP_mu<-res_mu/total_weighted_var_counts_mu
  class(vIP_mu)<-"varImpScores.bcfBMA"
  
  
  imp_vars2_tau=get_weighted_var_imp(num_vars=object$numvars,BIC=object$bic,sum_trees=object$sumoftrees_tau)
  res_tau <-apply(imp_vars2_tau[[4]],2,sum)
  #create varImpPlot command
  vIP_tau <-rep(NA,length(res_tau))
  total_weighted_var_counts_tau <-sum(res_tau)
  #now get variable inclusion probabilities
  vIP_tau<-res_tau/total_weighted_var_counts_tau
  class(vIP_tau)<-"varImpScores.bcfBMA"
  
  
  ret<-list()
  length(ret)<-2
  ret[[1]] <- vIP_mu
  ret[[2]] <- vIP_tau
  
  class(ret)<-"varimp_list.bcfBMA"  
  ret
  
}