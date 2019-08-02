#' @title Prediction intervals for BCF-BMA Treatment Effect estimates
#' 
#' @description This function produces prediction intervals for BCF-BMA treatment effects estimates.
#' @param object bartBMA object obtained from function bartBMA
#' @param num_iter Total number of iterations of the Gibbs sampler (including burn-in).
#' @param burnin Number of burn-on iterations of the Gibbs sampler.
#' @param l_quant Lower quartile of the prediction interval.
#' @param u_quant Upper quartile of the prediction interval.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residua;s.
#' @export 
#' @return The output is a list of length one. The one element in this list is a vector of prediction intervals???

pred_ints_parallel_LDL_bcf_TE<-function(object,num_iter,l_quant,u_quant,newdata=NULL,num_cores=1){
  if(l_quant>0.5 ||u_quant<0 ||u_quant>1){stop("Lower quantile must be lower than 0.5 and greater than 0")}
  if(u_quant<0.5 ||u_quant<0 ||u_quant>1){stop("Upper quantile must be greater than 0.5 and less than 1")}
  #object will be bartBMA object.
  
  
  if(is.null(newdata) && length(object)==31){
    #if test data specified separately
    ret<-pred_ints_lin_alg_fields_LDL_outsamp_bcf(object$sumoftrees_mu,object$sumoftrees_tau,
                                              object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                              object$response,
                                              object$bic,num_iter,object$nrowTrain,
                                              object$a_mu,object$a_tau,object$sigma,0,0,
                                              object$nu,object$lambda,
                                              object$z,
                                              object$test_data,
                                              object$test_pihat, object$test_z, object$include_pi2, 
                                              object$num_propscores, object$nrowtest,l_quant,u_quant,num_cores)
    print("Currently does not support positive semi-definite covariance matrices.")
    #Code should be rewritten so that mean_vars_lin_alg outputs the 
    
    
  }else{ if(is.null(newdata) && length(object)==25){
    #else return Pred Ints for training data
    ret<-pred_ints_lin_alg_fields_LDL_insamp_bcf(object$sumoftrees_mu,object$sumoftrees_tau,
                                             object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                             object$response,
                                             object$bic,num_iter,object$nrowTrain,
                                             object$a_mu,object$a_tau,object$sigma,0,0,
                                             object$nu,object$lambda,
                                             object$z,l_quant,u_quant,num_cores)
    
    
  }else{
    #if test data included in call to object
    ret<-pred_ints_lin_alg_fields_LDL_outsamp_bcf(object$sumoftrees_mu,object$sumoftrees_tau,
                                              object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                              object$response,
                                              object$bic,num_iter,object$nrowTrain,
                                              object$a_mu,object$a_tau,object$sigma,0,0,
                                              object$nu,object$lambda,
                                              object$z,
                                              object$test_data,
                                              object$test_pihat, object$test_z, object$include_pi2, 
                                              object$num_propscores, object$nrowtest,l_quant,u_quant,num_cores)
    
    
    
  }
  }
  
  
  class(ret)<-"pred_intervals.bcfBMA"  
  ret
  
}