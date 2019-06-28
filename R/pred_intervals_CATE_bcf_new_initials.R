#' @title Prediction intervals for BCF-BMA estimates of CATE, CATT, PATE, PATT
#' 
#' @description This function produces prediction intervals for BCF-BMA treatment effects estimates by post-hoc Gibbs-sampling from the full conditionals of the terminal node parameters and the variance of the error term. See Hernandez et al. (2018) Appendix D for details.
#' @param object bartBMA object obtained from function bartBMA
#' @param num_iter Total number of iterations of the Gibbs sampler (including burn-in).
#' @param burnin Number of burn-on iterations of the Gibbs sampler.
#' @param l_quant Lower quartile of the prediction interval.
#' @param u_quant Upper quartile of the prediction interval.
#' @param newdata Test data for which predictions are to be produced. Default = NULL. If NULL, then produces prediction intervals for training data if no test data was used in producing the bartBMA object, or produces prediction intervals for the original test data if test data was used in producing the bartBMA object.
#' @param update_resids Option for whether to update the partial residuals in the gibbs sampler. If equal to 1, updates partial residuals, if equal to zero, does not update partial residuals. The defaullt setting is to update the partial residua;s.
#' @export 
#' @return The output is a list of length one. The one element in this list is a vector of prediction intervals???

pred_intervals_CATE_bcf_new_initials <-function(object,num_iter,burnin,l_quant,u_quant,newdata=NULL,update_resids=1,
                                                trainingdata,pihatdata){
  if(l_quant>0.5 ||u_quant<0 ||u_quant>1){stop("Lower quantile must be lower than 0.5 and greater than 0")}
  if(u_quant<0.5 ||u_quant<0 ||u_quant>1){stop("Upper quantile must be greater than 0.5 and less than 1")}
  #object will be bartBMA object.
  
  
  scaled_train_y <- scale_response_bcf(min(object$response),max(object$response),-0.5,0.5,object$response)
  
  get_resids <- get_initial_resids(trainingdata,
                                   pihatdata,
                                   object$sumoftrees_mu,
                                   object$sumoftrees_tau,
                                   scaled_train_y,
                                   object$z)
  
  diff_inital_resids_mu <- get_resids[[1]]
  new_pred_list1_mu <- get_resids[[2]]
  diff_inital_resids_tau <- get_resids[[3]]
  new_pred_list1_tau <- get_resids[[4]]
  
  
  if(update_resids==0){
    if(is.null(newdata) && length(object)==31){
      #if test data specified separately
      gs_chains<-gibbs_sampler_no_update_new_inits(object$sumoftrees_mu,object$sumoftrees_tau,object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                         object$response,
                                         object$bic,num_iter, burnin,object$nrowTrain,
                                         object$a_mu,object$a_tau,object$sigma,0,0,
                                         object$nu,object$lambda,
                                         diff_inital_resids_mu,
                                         diff_inital_resids_tau,
                                         object$z,
                                         object$test_data,
                                         object$test_pihat, object$test_z, object$include_pi2, object$num_propscores, object$nrowtest,
                                         new_pred_list1_mu,
                                         new_pred_list1_tau)
    }else{ if(is.null(newdata) && length(object)==25){
      #else return Pred Ints for training data
      gs_chains<-gibbs_sampler_no_update2_new_inits(object$sumoftrees_mu,object$sumoftrees_tau,object$obs_to_termNodesMatrix_mu,
                                          object$obs_to_termNodesMatrix_tau,
                                          object$response,
                                          object$bic,num_iter, burnin,object$nrowTrain,
                                          object$a_mu,object$a_tau,object$sigma,0,0,
                                          object$nu,object$lambda,
                                          diff_inital_resids_mu,
                                          diff_inital_resids_tau,
                                          object$z,
                                          new_pred_list1_mu,
                                          new_pred_list1_tau)
      
    }else{
      #if test data included in call to object
      gs_chains<-gibbs_sampler_no_update_new_inits(object$sumoftrees_mu,object$sumoftrees_tau,object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                         object$response,
                                         object$bic,num_iter, burnin,object$nrowTrain,
                                         object$a_mu,object$a_tau,object$sigma,0,0,
                                         object$nu,object$lambda,
                                         diff_inital_resids_mu,
                                         diff_inital_resids_tau,
                                         object$z,
                                         object$test_data,
                                         object$test_pihat, object$test_z, object$include_pi2, object$num_propscores, 
                                         object$nrowtest,
                                         new_pred_list1_mu,
                                         new_pred_list1_tau)
    }
    }
  }
  else{
    if(is.null(newdata) && length(object)==31){
      #if test data specified separately
      gs_chains<-gibbs_sampler_new_inits(object$sumoftrees_mu,object$sumoftrees_tau,object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                               object$response,
                               object$bic,num_iter, burnin,object$nrowTrain,
                               object$a_mu,object$a_tau,object$sigma,0,0,
                               object$nu,object$lambda,
                               diff_inital_resids_mu,
                               diff_inital_resids_tau,
                               object$z,
                               object$test_data,
                               object$test_pihat, object$test_z, object$include_pi2, object$num_propscores, 
                               object$nrowtest,
                               new_pred_list1_mu,
                               new_pred_list1_tau)
    }else{ if(is.null(newdata) && length(object)==25){
      #else return Pred Ints for training data
      gs_chains<-gibbs_sampler2_new_inits(object$sumoftrees_mu,object$sumoftrees_tau,object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                                object$response,
                                object$bic,num_iter, burnin,object$nrowTrain,
                                object$a_mu,object$a_tau,object$sigma,0,0,
                                object$nu,object$lambda,
                                diff_inital_resids_mu,
                                diff_inital_resids_tau,
                                object$z,
                                new_pred_list1_mu,
                                new_pred_list1_tau)
      
    }else{
      #if test data included in call to object
      gs_chains<-gibbs_sampler_new_inits(object$sumoftrees_mu,object$sumoftrees_tau,object$obs_to_termNodesMatrix_mu,object$obs_to_termNodesMatrix_tau,
                               object$response,
                               object$bic,num_iter, burnin,object$nrowTrain,
                               object$a_mu,object$a_tau,object$sigma,0,0,
                               object$nu,object$lambda,
                               diff_inital_resids_mu,
                               diff_inital_resids_tau,
                               object$z,
                               object$test_data,
                               object$test_pihat, object$test_z, object$include_pi2, object$num_propscores,
                               object$nrowtest,
                               new_pred_list1_mu,
                               new_pred_list1_tau)
    }
    }
  }
  #y_posterior_sum_trees<-gs_chains[[4]]
  #y_orig_post_sum_trees<-gs_chains[[5]]
  #sigma_chains<-gs_chains[[3]]
  if(is.null(newdata) && length(object)==31){
    #y_posterior_sum_trees<-gs_chains#[[8]]
    y_orig_post_sum_trees<-gs_chains[[2]]#[[9]]
    sigma_chains<-gs_chains[[1]]#[[3]]
  }else if(is.null(newdata) && length(object)==25){
    #y_posterior_sum_trees<-gs_chains#[[6]]
    y_orig_post_sum_trees<-gs_chains[[2]]#[[7]]
    sigma_chains<-gs_chains[[1]]#[[3]]
    
  }else{
    #y_posterior_sum_trees<-gs_chains#[[8]]
    y_orig_post_sum_trees<-gs_chains[[2]]#[[9]]
    sigma_chains<-gs_chains[[1]]#[[3]]
  } 
  
  sum_of_tree_BIC<- -0.5*object$bic
  weights<-exp(sum_of_tree_BIC-(max(sum_of_tree_BIC)+log(sum(exp(sum_of_tree_BIC-max(sum_of_tree_BIC))))))
  #final_length<-num_iter-burnin
  num_its_to_sample<-round(weights*(num_iter-burnin))
  #final_sigma_chain<-numeric(0)
  
  #final_y_chain<-matrix(nrow=0,ncol=ncol(y_posterior_sum_trees[[1]]))
  final_yorig_chain<-matrix(nrow=0,ncol=ncol(y_orig_post_sum_trees[[1]]))
  
  for(i in 1:length(sigma_chains)){
    sample_its<-sample(burnin:num_iter,num_its_to_sample[i])
    #final_sigma_chain<-c(final_sigma_chain,sigma_chains[[i]][sample_its])
    #now do the same for predicted response updates
    #post_y_i<-y_posterior_sum_trees[[i]]
    post_yorig_i<-y_orig_post_sum_trees[[i]]
    #final_y_chain<-rbind(final_y_chain,post_y_i[sample_its,])
    final_yorig_chain<-rbind(final_yorig_chain,post_yorig_i[sample_its,])
    
  }
  PI<-apply(final_yorig_chain,2,function(x)quantile(x,probs=c(l_quant,0.5,u_quant)))
  meanpreds<-apply(final_yorig_chain,2,mean)
  
  z_train2 <- object$z
  CATEs_across_iterations <- apply(final_yorig_chain,1,mean)
  CATE_est <- mean(CATEs_across_iterations)
  CATE_PI<-as.matrix(quantile(CATEs_across_iterations,probs=c(l_quant,0.5,u_quant)))
  
  PATE_var<- mean(apply(final_yorig_chain,1:2,function(x){(x-CATE_est)^2}))
  PATE_PI <- as.matrix(c(CATE_est + qnorm(l_quant)*sqrt(PATE_var),CATE_est ,  CATE_est + qnorm(u_quant)*sqrt(PATE_var) ))
  
  CATTs_across_iterations <- apply(final_yorig_chain[,z_train2],1,mean)
  CATT_est <- mean(CATTs_across_iterations)
  CATT_PI<-as.matrix(quantile(CATTs_across_iterations,probs=c(l_quant,0.5,u_quant)))
  
  PATT_var<- mean(apply(final_yorig_chain[,z_train2],1:2,function(x){(x-CATT_est)^2}))
  PATT_PI <- as.matrix(c(CATT_est + qnorm(l_quant)*sqrt(PATT_var),CATT_est ,  CATT_est + qnorm(u_quant)*sqrt(PATT_var) ))
  
  
  
  CATNTs_across_iterations <- apply(final_yorig_chain[,1-z_train2],1,mean)
  CATNT_est <- mean(CATNTs_across_iterations)
  CATNT_PI<-as.matrix(quantile(CATNTs_across_iterations,probs=c(l_quant,0.5,u_quant)))
  
  PATNT_var<- mean(apply(final_yorig_chain[,1-z_train2],1:2,function(x){(x-CATNT_est)^2}))
  PATNT_PI <- as.matrix(c(CATNT_est + qnorm(l_quant)*sqrt(PATNT_var),CATNT_est ,  CATNT_est + qnorm(u_quant)*sqrt(PATNT_var) ))
  
  
  ret<-list()
  ret$ITE_PI<-PI
  ret$ITE_mean <- meanpreds
  ret$CATE_est <- CATE_est
  ret$CATE_PI <- CATE_PI
  ret$CATT_est <- CATT_est
  ret$CATT_PI <- CATT_PI
  ret$CATNT_est <- CATNT_est
  ret$CATNT_PI <- CATNT_PI
  ret$PATE_PI <- PATE_PI
  ret$PATT_PI <- PATT_PI
  ret$PATNT_PI <- PATNT_PI
  
  
  class(ret)<-"CATE_intervals.bcfBMA"  
  ret
}