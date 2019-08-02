#' @title Treatment Effect estimates for BCF-BMA output obtained from the posterior probability weighted averaged of the posterior means for each model
#' 
#' @description This function produces treatment effect predictions from BCF-BMA by obtaining the posterior probability weighted average of the posterior means for each model.
#' @param object bcfBMA object obtained from function bcfBMA
#' @param num_iter Total number of iterations of the Gibbs sampler (including burn-in).
#' @param burnin Number of burn-on iterations of the Gibbs sampler.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residua;s.
#' @export 
#' @return The output is a list of length one. The one element in this list is a vector of prediction intervals???

preds_bcfbma_lin_alg<-function(object,num_iter,newdata=NULL,trainingdata){
  #object will be bartBMA object.
  
  
  
  
  if(is.null(newdata) && length(object)==31){
    #if test data specified separately
    ret<-preds_bcfbma_lin_alg_outsamp(object$sumoftrees_mu,object$sumoftrees_tau,
                                                 object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                                 object$response,
                                                 object$bic,num_iter,object$nrowTrain,
                                                 object$a_mu,object$a_tau,object$sigma,0,0,
                                                 object$nu,object$lambda,
                                                 object$z,
                                                 object$test_data,
                                                 object$test_pihat, object$test_z, object$include_pi2,
                                                 object$num_propscores, object$nrowtest)
  }else{if(is.null(newdata) && length(object)==25){
    #else return Pred Ints for training data
    ret<-preds_bcfbma_lin_alg_insamp(object$sumoftrees_mu,
                                     object$sumoftrees_tau,
                                   object$obs_to_termNodesMatrix_mu,
                                   object$obs_to_termNodesMatrix_tau,
                                   object$response,
                                   object$bic,
                                   num_iter,
                                   object$nrowTrain,
                                   object$a_mu,
                                   object$a_tau,
                                   object$sigma,0,0,
                                   object$nu,
                                   object$lambda,
                                   object$z)
    
    
  }else{
    #if test data included in call to object
    ret<-preds_bcfbma_lin_alg_outsamp(object$sumoftrees_mu,object$sumoftrees_tau,
                                    object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                    object$response,
                                    object$bic,num_iter,object$nrowTrain,
                                    object$a_mu,object$a_tau,object$sigma,0,0,
                                    object$nu,object$lambda,
                                    object$z,
                                    object$test_data,
                                    object$test_pihat, object$test_z, object$include_pi2, 
                                    object$num_propscores, object$nrowtest)
  }}
  
  
  
  ret
  
}