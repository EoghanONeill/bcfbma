//#########################################################################################################################//
  
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector get_imp_vars(NumericVector split_vars,int num_col,NumericVector current_vars){
  
  // NumericVector vars_chosen=sort_unique(split_vars);
  // 
  // if(vars_chosen[0]==0){
  //   vars_chosen.erase(0);
  // }
  //if(vars_chosen.size()!=0){
    for(int i=0;i<split_vars.size();i++){      
      if(split_vars[i]!=0){      
        current_vars[split_vars[i]-1]+=1;
      }
      
    }
  //}
  return(current_vars);
}
//#######################################################################################################################//
  
  #include <Rcpp.h>
  using namespace Rcpp;
//' @title Obtain weighted variable importances
//' @export
// [[Rcpp::export]]

List get_weighted_var_imp(int num_vars,NumericVector BIC,List sum_trees){
  
  NumericMatrix vars_for_all_trees(sum_trees.size(),num_vars);
  NumericMatrix weighted_vars_for_all_trees(sum_trees.size(),num_vars);
  //NumericVector weighted_BIC=BIC/sum(BIC);
  
  
  NumericVector BICi=-0.5*BIC;
  double max_BIC=max(BICi);
  
  // weighted_BIC is actually the posterior model probability
  NumericVector weighted_BIC(BIC.size());
  
  for(int k=0;k<BIC.size();k++){
    
    //NumericVector BICi=-0.5*BIC_weights;
    //double max_BIC=max(BICi);
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));
    weighted_BIC[k]=weight;
    //int num_its_to_sample = round(weight*(num_iter));
    
  }
  
  
  
  
  for(int i=0;i<sum_trees.size();i++){
    NumericVector selected_variables(num_vars);
    //for each set of trees loop over individual trees
    SEXP s = sum_trees[i];
    
    if(is<List>(s)){
      List tree_set=sum_trees[i];
      for(int j=0;j<tree_set.size();j++){
        //for each tree in current list get the variables selected for given tree and add to row i of vars_for_all_trees
        NumericMatrix tree_data=tree_set[j];
        //have current tree get the split variables used
        NumericVector tree_vars=tree_data(_,2);
        selected_variables=get_imp_vars(tree_vars,num_vars,selected_variables);
        
        //NumericVector temp_counts_so_far = vars_for_all_trees(i,_);
        //vars_for_all_trees(i,_)= temp_counts_so_far + selected_variables;
        if(j==tree_set.size()-1){
          vars_for_all_trees(i,_)= selected_variables;
        }
        
      }
    }else{
      NumericMatrix tree_data=sum_trees[i];
      //get variables selected for current tree and add to row i or vars_for_all_trees
      NumericVector tree_vars=tree_data(_,2);
      selected_variables=get_imp_vars(tree_vars,num_vars,selected_variables);
      vars_for_all_trees(i,_)=selected_variables;
      
    }
    weighted_vars_for_all_trees(i,_)=vars_for_all_trees(i,_)*weighted_BIC[i];
  }
  
  //get BIC and weight the vars by their BIC/sum(BIC)
  
  List ret(4);
  ret[0]=weighted_BIC; // This is actually the posterior model probability
  ret[1]=BIC;
  ret[2]=vars_for_all_trees;
  ret[3]=weighted_vars_for_all_trees;
  
  return(ret);  
}