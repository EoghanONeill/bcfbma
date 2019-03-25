//This code will take the output of BCF-BMA (list of sums of trees) and will update the predicted values for the terminal node
//means and model variance

//first take one set of sum of trees:
//this will work for the training data only, this will be updated for test data as an optional parameter later on 
//(will also have to update external predict function for test predictions)
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_nodes_gs(NumericMatrix tree_table){
  NumericVector terminal_nodes;
  arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false); 
  
  arma::vec colmat=arma_tree.col(4);
  arma::uvec term_nodes=arma::find(colmat==-1);
  term_nodes=term_nodes+1;
  return(wrap(term_nodes));
}

//################################################################################################################################//
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]

NumericVector find_term_obs_gs(NumericMatrix tree_matrix_temp,double terminal_node){
  //NumericVector term_obs;
  arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); 
  arma::uvec term_obs;
  
  for(int j=0;j<tree_matrix_temp.ncol();j++){
    arma::vec colmat=arma_tree_mat.col(j);
    term_obs=arma::find(colmat==terminal_node);
    if(term_obs.size()>0){
      break;
    }
  }
  
  return(wrap(term_obs));
}
//################################################################################################################################//
// [[Rcpp::export]]
NumericVector calc_rowsums(NumericMatrix predictions){
  arma::mat M1(predictions.begin(), predictions.nrow(), predictions.ncol(), false);
  arma::colvec predicted_values=sum(M1,1);
  return(wrap(predicted_values));
}
//############################################################################################//
// [[Rcpp::export]]
NumericVector calculate_resids(NumericMatrix predictions,NumericVector response){
  NumericVector resids=response.size();
  NumericVector row_sums=calc_rowsums(predictions);
  resids=response - row_sums;
  return(resids);
}

//################################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List update_Gibbs_mean_var(NumericMatrix tree_table,NumericVector resids,double a,double sigma,double mu_mu,IntegerVector terminal_nodes,  List term_obs_tree){
  List update_params(2);
  
  NumericVector mu_ij(terminal_nodes.size());
  NumericVector Tj(terminal_nodes.size());
  NumericVector new_mean(terminal_nodes.size());
  NumericVector new_var(terminal_nodes.size());
  
  for(int k=0;k< terminal_nodes.size();k++){
    //update the node means    
    //find which observations are assigned to terminal node k in tree i
    
    IntegerVector term_obs=term_obs_tree[k];
    
    //get the number of observations in node k
    NumericVector temp_obs=resids[term_obs];
    mu_ij[k]=std::accumulate(temp_obs.begin(),temp_obs.end(), 0.0);
    Tj[k]=term_obs.size();
    
    new_mean[k]=(mu_ij[k]+a*mu_mu)/(Tj[k]+a);
    
    new_var[k]=(1/((1/(pow(sigma,2)))*(Tj[k]+a)));
    
    NumericVector temp;
    term_obs=temp;
    
  }  
  
  update_params[0]=new_mean;
  update_params[1]=new_var;
  
  return(update_params);
}
//################################################################################################################################//

// [[Rcpp::export]]
double update_sigma(double a1,double b,NumericVector resids,int n){
  NumericVector sq_resid=resids*resids;
  double ssr= std::accumulate(sq_resid.begin(),sq_resid.end(),0.0);
  double shape=(a1+n/2);
  double rate =((ssr/2)+(1/b));
  RNGScope scope;
  //double tau =shape/rate;
  
  double tau =R::rgamma(shape,1/rate);
  double sigma = sqrt((1/tau));
  return sigma;
}

//################################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_node_means(NumericMatrix sum_tree, NumericVector term_nodes){
  NumericVector means=sum_tree(_,5);
  arma::vec node_means = as<arma::vec>(means);
  arma::uvec arma_term_nodes=as<arma::uvec>(term_nodes);
  arma_term_nodes=arma_term_nodes-1;
  arma::vec final_means=node_means.elem(arma_term_nodes);
  return(wrap(final_means));
  
}
//################################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_tree_info(List overall_sum_trees,List overall_sum_mat,int num_obs){
  List overall_term_nodes_trees(overall_sum_trees.size());
  List overall_term_obs_trees(overall_sum_trees.size());
  List overall_predictions(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = overall_sum_trees[i];
    
    NumericVector test_preds_sum_tree;
    
    if(is<List>(s)){
      //if current set of trees contains more than one tree
      List sum_tree=overall_sum_trees[i];
      
      List sum_tree_mat=overall_sum_mat[i];
      
      //save all info in list of list format the same as the trees.
      
      List term_nodes_trees(sum_tree.size());
      List term_obs_trees(sum_tree.size());
      NumericMatrix predictions(num_obs,sum_tree.size());
      
      for(int k =0;k<sum_tree.size();k++){
        
        NumericMatrix tree_table=sum_tree[k];
        NumericMatrix tree_mat=sum_tree_mat[k];
        NumericVector term_nodes=find_term_nodes_gs(tree_table);
        term_nodes_trees[k]=term_nodes;
        List term_obs_tree(term_nodes.size());
        NumericVector term_preds(num_obs);
        
        for(int j=0;j<term_nodes.size();j++){
          double terminal_node= term_nodes[j]; 
          NumericVector term_obs=find_term_obs_gs(tree_mat,terminal_node);
          NumericVector node_means=find_node_means(tree_table,term_nodes);
          term_obs_tree[j]=term_obs;
          double node_mean=node_means[j];
          term_preds[term_obs]=node_mean; 
        }          
        term_obs_trees[k]=term_obs_tree;
        
        predictions(_,k)=term_preds;
      } 
      overall_term_nodes_trees[i]=term_nodes_trees;
      overall_term_obs_trees[i]= term_obs_trees;
      overall_predictions[i]=predictions;
    }else{
      NumericMatrix sum_tree=overall_sum_trees[i];
      NumericMatrix tree_mat=overall_sum_mat[i];
      NumericVector term_nodes=find_term_nodes_gs(sum_tree);
      NumericVector node_means=find_node_means(sum_tree,term_nodes);
      List term_obs_tree(term_nodes.size());
      overall_term_nodes_trees[i]=term_nodes;
      NumericVector predictions(num_obs);
      
      for(int j=0;j<term_nodes.size();j++){
        double terminal_node= term_nodes[j];
        double node_mean=node_means[j];
        NumericVector term_obs=find_term_obs_gs(tree_mat,terminal_node);
        term_obs_tree[j]=term_obs;
        predictions[term_obs]=node_mean;
      }
      overall_term_obs_trees[i]= term_obs_tree;
      overall_predictions[i]=predictions;
    }  
  }    
  List ret(3);
  ret[0]=overall_term_nodes_trees;
  ret[1]=overall_term_obs_trees;
  ret[2]=overall_predictions;
  return(ret);
}
//################################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_tree_info_tau(List overall_sum_trees,List overall_sum_mat,int num_obs, NumericVector z){
  arma::vec z_a=Rcpp::as<arma::vec>(z);
  arma::uvec treated_obs_a=arma::find(z_a==1);
  NumericVector treated_obs=wrap(treated_obs_a);
  
  List overall_term_nodes_trees(overall_sum_trees.size());
  List overall_term_obs_trees(overall_sum_trees.size());
  List overall_predictions(overall_sum_trees.size());
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = overall_sum_trees[i];
    
    NumericVector test_preds_sum_tree;
    
    if(is<List>(s)){
      //if current set of trees contains more than one tree
      List sum_tree=overall_sum_trees[i];
      
      List sum_tree_mat=overall_sum_mat[i];
      
      //save all info in list of list format the same as the trees.
      
      List term_nodes_trees(sum_tree.size());
      List term_obs_trees(sum_tree.size());
      NumericMatrix predictions(num_obs,sum_tree.size());
      
      for(int k =0;k<sum_tree.size();k++){
        
        NumericMatrix tree_table=sum_tree[k];
        NumericMatrix tree_mat=sum_tree_mat[k];
        NumericVector term_nodes=find_term_nodes_gs(tree_table);
        term_nodes_trees[k]=term_nodes;
        List term_obs_tree(term_nodes.size());
        NumericVector term_preds(num_obs);
        
        for(int j=0;j<term_nodes.size();j++){
          double terminal_node= term_nodes[j]; 
          NumericVector term_obs=find_term_obs_gs(tree_mat,terminal_node);
          NumericVector treated_term_obs=intersect(treated_obs,term_obs);
          
          NumericVector node_means=find_node_means(tree_table,term_nodes);
          term_obs_tree[j]=treated_term_obs;
          double node_mean=node_means[j];
          term_preds[term_obs]=node_mean; 
        }          
        term_obs_trees[k]=term_obs_tree;
        
        predictions(_,k)=term_preds;
      } 
      overall_term_nodes_trees[i]=term_nodes_trees;
      overall_term_obs_trees[i]= term_obs_trees;
      overall_predictions[i]=predictions;
    }else{
      NumericMatrix sum_tree=overall_sum_trees[i];
      NumericMatrix tree_mat=overall_sum_mat[i];
      NumericVector term_nodes=find_term_nodes_gs(sum_tree);
      NumericVector node_means=find_node_means(sum_tree,term_nodes);
      List term_obs_tree(term_nodes.size());
      overall_term_nodes_trees[i]=term_nodes;
      NumericVector predictions(num_obs);
      
      for(int j=0;j<term_nodes.size();j++){
        double terminal_node= term_nodes[j];
        double node_mean=node_means[j];
        NumericVector term_obs=find_term_obs_gs(tree_mat,terminal_node);
        NumericVector treated_term_obs=intersect(treated_obs,term_obs);
        
        term_obs_tree[j]=treated_term_obs;
        predictions[term_obs]=node_mean;
      }
      overall_term_obs_trees[i]= term_obs_tree;
      overall_predictions[i]=predictions;
    }  
  }    
  List ret(3);
  ret[0]=overall_term_nodes_trees;
  ret[1]=overall_term_obs_trees;
  ret[2]=overall_predictions;//NOTE THAT THIS INCLUDES THE PREDICTED TE VALUES FOR BOTH TREATED AND NON-TREATED OBSERVATIONS
  return(ret);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix remove_curr_col(NumericMatrix predy,int i){
  arma::mat M=Rcpp::as<arma::mat>(predy);
  M.shed_col(i);
  NumericMatrix s=as<NumericMatrix>(wrap(M));
  return(s);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector get_new_mean(IntegerVector terminal_nodes,List new_mean_var){
  NumericVector node_means;
  for(int k=0;k<terminal_nodes.size();k++){
    NumericVector sd=new_mean_var[1];
    NumericVector temp_mean=new_mean_var[0];
    //double new_mean=temp_mean[k];
    
    double new_mean= R::rnorm(temp_mean[k],sqrt(sd[k]));
    node_means.push_back(new_mean);
    
  }
  
  return(node_means);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List update_predictions_gs(NumericMatrix tree_table,NumericVector new_mean,NumericVector new_var,int n,
                           IntegerVector terminal_nodes,List term_obs_tree){
  
  List updated_preds(2);
  NumericVector new_preds(n);
  
  for(int k=0;k<terminal_nodes.size();k++){
    
    NumericVector term_obs=term_obs_tree[k];        
    //update the mean of the selected tree nodes:
    tree_table(terminal_nodes[k]-1,5)= new_mean[k];
    tree_table(terminal_nodes[k]-1,6)=sqrt(new_var[k]);
    double newmean=new_mean[k];
    //update residuals for next iteration
    new_preds[term_obs]=newmean;
  }
  
  updated_preds[0]=tree_table;
  updated_preds[1]=new_preds;
  
  return(updated_preds);
}
//################################################################################################################################//
// [[Rcpp::export]]
NumericVector scale_response_gs(double a,double b,double c,double d,NumericVector y){
  NumericVector y_scaled = -((-b*c+a*d)/(-a+b))+((-c+d)*y/(-a+b));
  return(y_scaled);
}
//###########################################################################################################################//
// [[Rcpp::export]]
NumericVector get_original_gs(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
  NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
  return(original_y);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]

NumericVector find_internal_nodes_gs(NumericMatrix treetable){
  NumericVector internal_nodes;
  
  for(int l=0;l<treetable.nrow();l++){    
    if(treetable(l,4)==1){
      internal_nodes.push_back(l+1);
    }
  }
  
  NumericVector internal_nodes_sort = clone(internal_nodes);
  std::sort(internal_nodes.begin(), internal_nodes.end());
  
  return(internal_nodes_sort);
}
//###########################################################################################################################//
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List get_tree_info_test_data(NumericMatrix test_data,NumericMatrix tree_data) {
  // Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  
  //test_data is a nxp matrix with the same variable names as the training data the model was built on
  
  //tree_data is the tree table with the tree information i.e. split points and split variables and terminal node mean values
  
  //term_node_means is a vector storing the terminal node mean values
  
  arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);
  arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);
  NumericVector internal_nodes=find_internal_nodes_gs(tree_data);
  NumericVector terminal_nodes=find_term_nodes_gs(tree_data);
  arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);
  NumericVector tree_predictions;
  
  //now for each internal node find the observations that belong to the terminal nodes
  
  NumericVector predictions(test_data.nrow());
  List term_obs(terminal_nodes.size());
  for(int i=0;i<terminal_nodes.size();i++){
    arma::mat subdata=testd;
    int curr_term=terminal_nodes[i];
    int row_index;
    int term_node=terminal_nodes[i];
    
    if(curr_term % 2==0){
      //term node is left daughter
      row_index=terminal_nodes[i];
    }else{
      //term node is right daughter
      row_index=terminal_nodes[i]-1;
    }
    
    //save the left and right node data into arma uvec
    
    arma::vec left_nodes=arma_tree.col(0);
    arma::vec right_nodes=arma_tree.col(1);
    arma::mat node_split_mat;    
    node_split_mat.set_size(0,3);
    
    while(row_index!=1){
      //for each terminal node work backwards and see if the parent node was a left or right node
      //append split info to a matrix 
      int rd=0;
      arma::uvec parent_node=arma::find(left_nodes == term_node);
      
      if(parent_node.size()==0){
        parent_node=arma::find(right_nodes == term_node);
        rd=1;
      }
      
      //want to cout parent node and append to node_split_mat
      
      node_split_mat.insert_rows(0,1);
      node_split_mat(0,0)=tree_data(parent_node[0],2);
      node_split_mat(0,1)=tree_data(parent_node[0],3);
      node_split_mat(0,2)=rd;     
      row_index=parent_node[0]+1;
      term_node=parent_node[0]+1;
    }
    
    //once we have the split info, loop through rows and find the subset indexes for that terminal node!
    //then fill in the predicted value for that tree
    //double prediction = tree_data(term_node,5);
    arma::uvec pred_indices;
    int split= node_split_mat(0,0)-1;
    arma::vec tempvec = testd.col(split);
    double temp_split = node_split_mat(0,1);
    
    if(node_split_mat(0,2)==0){
      pred_indices = arma::find(tempvec <= temp_split);
    }else{
      pred_indices = arma::find(tempvec > temp_split);      
    }
    
    arma::uvec temp_pred_indices;
    arma::vec data_subset = testd.col(split);
    data_subset=data_subset.elem(pred_indices);
    
    //now loop through each row of node_split_mat
    int n=node_split_mat.n_rows;
    
    for(int j=1;j<n;j++){
      int curr_sv=node_split_mat(j,0);
      double split_p = node_split_mat(j,1);
      
      data_subset = testd.col(curr_sv-1);
      data_subset=data_subset.elem(pred_indices);
      
      if(node_split_mat(j,2)==0){
        //split is to the left
        temp_pred_indices=arma::find(data_subset <= split_p);
      }else{
        //split is to the right
        temp_pred_indices=arma::find(data_subset > split_p);
      }
      pred_indices=pred_indices.elem(temp_pred_indices);
      
      if(pred_indices.size()==0){
        continue;
      }
      
    }
    
    double nodemean=tree_data(terminal_nodes[i]-1,5);
    IntegerVector predind=as<IntegerVector>(wrap(pred_indices));
    predictions[predind]= nodemean;
    term_obs[i]=predind;
  } 
  List ret(3);
  ret[0] = terminal_nodes;
  ret[1] = term_obs;
  ret[2] = predictions;
  return(ret);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_tree_info_testdata_overall(List overall_sum_trees,int num_obs,NumericMatrix test_data){
  
  List overall_term_nodes_trees(overall_sum_trees.size());
  List overall_term_obs_trees(overall_sum_trees.size());
  List overall_predictions(overall_sum_trees.size());
  
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees
    SEXP s = overall_sum_trees[i];
    
    NumericVector test_preds_sum_tree;
    if(is<List>(s)){
      //if current set of trees contains more than one tree...usually does!
      List sum_tree=overall_sum_trees[i];        
      //save all info in list of list format the same as the trees.       
      List term_nodes_trees(sum_tree.size());
      List term_obs_trees(sum_tree.size());
      NumericMatrix predictions(num_obs,sum_tree.size());
      
      for(int k=0;k<sum_tree.size();k++){          
        NumericMatrix tree_table=sum_tree[k];
        List tree_info=get_tree_info_test_data(test_data, tree_table) ;
        NumericVector term_nodes=tree_info[0];
        term_nodes_trees[k]=term_nodes;
        term_obs_trees[k]=tree_info[1];
        NumericVector term_preds=tree_info[2];
        predictions(_,k)=term_preds;
      } 
      overall_term_nodes_trees[i]=term_nodes_trees;
      overall_term_obs_trees[i]= term_obs_trees;
      overall_predictions[i]=predictions;
    }else{
      NumericMatrix sum_tree=overall_sum_trees[i];
      List tree_info=get_tree_info_test_data(test_data, sum_tree) ;
      overall_term_nodes_trees[i]=tree_info[0];
      List term_obs_tree=tree_info[1];
      NumericVector term_preds=tree_info[2];
      NumericVector predictions=term_preds;   
      overall_term_obs_trees[i]= term_obs_tree;
      overall_predictions[i]=predictions;
    }  
  }    
  List ret(3);
  ret[0]=overall_term_nodes_trees;
  ret[1]=overall_term_obs_trees;
  ret[2]=overall_predictions;
  return(ret);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::depends(RcppArmadillo)]]
//' @title Obtain draws from gibbs sampler
//' @export
// [[Rcpp::export]]

List gibbs_sampler(List overall_sum_trees,
                   List overall_sum_mat,
                   NumericVector y,
                   NumericVector BIC_weights,
                   int num_iter,int burnin,int num_obs,int num_test_obs,
                   double a,double sigma,double mu_mu,double nu,double lambda,
                   List resids,
                   NumericMatrix test_data){
  if(burnin>=num_iter){
    throw std::range_error("Number of iterations has to be greater than the number of burn-in samples");
  }
  //get tree info just once don't repeat for each iteration:
  //need terminal nodes, terminal node means and observations terminal nodes refer to
  List tree_info=get_tree_info(overall_sum_trees,overall_sum_mat,num_obs);
  List test_tree_info=get_tree_info_testdata_overall(overall_sum_trees,num_test_obs,test_data);
  List overall_term_test_obs_trees=test_tree_info[1];
  
  double a1=nu/2;
  double b=2/(nu*lambda);
  double sigma2=sigma;
  
  double sigma_init=sigma2;
  List overall_term_nodes_trees=tree_info[0];
  List overall_term_obs_trees=tree_info[1];
  List overall_predictions=tree_info[2];
  NumericVector y_scaled=scale_response_gs(min(y),max(y),-0.5,0.5,y); 
  List prediction_list(overall_sum_trees.size());
  List prediction_list_orig(overall_sum_trees.size());
  List prediction_test_list(overall_sum_trees.size());
  List predictive_dist_train_list(overall_sum_trees.size());
  List predictive_dist_train_list_orig(overall_sum_trees.size());
  List prediction_test_list_orig(overall_sum_trees.size());
  List overall_sum_trees1=clone(overall_sum_trees);
  List overall_sum_mat1=clone(overall_sum_mat); 
  List sigma_chains(overall_sum_trees.size());
  //int one_tree=0;
  for(int i=0;i<overall_sum_trees.size();i++){
    //for each set of trees loop over individual trees    
    NumericVector sigma_its(num_iter);
    sigma2=sigma_init;
    SEXP s = overall_sum_trees[i];
    NumericVector test_preds_sum_tree;
    NumericMatrix sum_predictions;
    NumericMatrix sum_test_predictions;
    
    if(is<List>(s)){
      //if current set of trees contains more than one tree
      List sum_tree=overall_sum_trees1[i];
      List sum_tree_mat=overall_sum_mat1[i];
      List sum_term_nodes=overall_term_nodes_trees[i];
      List sum_term_obs=overall_term_obs_trees[i];
      List sum_term_test_obs=overall_term_test_obs_trees[i];
      List sum_resids0=resids[i];
      List sum_resids=clone(sum_resids0);
      
      NumericMatrix tree_predictions=overall_predictions[i];
      sum_predictions=clone(tree_predictions);
      NumericMatrix post_predictions(num_iter,num_obs);
      NumericMatrix post_predictions_PI(num_iter,num_obs);
      NumericMatrix post_predictions_orig_PI(num_iter,num_obs);
      NumericMatrix post_predictions_orig(num_iter,num_obs);
      NumericMatrix post_test_predictions(num_iter,num_test_obs);
      NumericMatrix post_test_predictions_orig(num_iter,num_test_obs);
      NumericMatrix sum_new_predictions(sum_predictions.nrow(),sum_predictions.ncol());
      NumericMatrix sum_new_test_predictions(num_test_obs,sum_predictions.ncol());
      
      for(int j=0;j<num_iter;j++){
        for(int k =0;k<sum_tree.size();k++){
          NumericMatrix tree_table=sum_tree[k];
          IntegerMatrix tree_mat=sum_tree_mat[k];
          //find terminal node means and observations associated with them
          IntegerVector term_nodes=sum_term_nodes[k];
          List term_obs=sum_term_obs[k];
          List term_test_obs=sum_term_test_obs[k];
          NumericVector predictions=sum_resids[k];
          //current predictions are the residuals for sum of trees!
          
          //update the means and predictions for tree
          List new_node_mean_var=update_Gibbs_mean_var(tree_table,predictions,a,sigma2,mu_mu,term_nodes,term_obs);
          NumericVector new_node_mean=get_new_mean(term_nodes,new_node_mean_var);
          NumericVector new_node_var=new_node_mean_var[1];
          //update predictions by setting predicted value for term_obs[termnode]=new mean value!
          
          List updated_preds=update_predictions_gs(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_obs);         
          NumericVector temp_preds=updated_preds[1];
          
          //NOW UPDATE THE RESIDUALS FOR USE IN NEXT ITERATION
          //THE PLACING OF THIS SECTION OF CODE HERE IS IMPORTANT
          //MUST BE BEFORE sum_new_predictions(_,k) is updated so that the
          //previous round's predictions can be added (i.e. removed from the residual before the new predictions are taken away to create the new residual) 
          for(int q=0;q<sum_tree.size();q++){
            if(q!=k){
              if(j==0){
                NumericVector temp_resids1= sum_resids[q];
                sum_resids[q]=temp_resids1+tree_predictions(_,k)-temp_preds;
              }else{
                NumericVector temp_resids1= sum_resids[q];
                sum_resids[q]=temp_resids1+sum_new_predictions(_,k)-temp_preds;
              }
            }
          }
          sum_new_predictions(_,k)=temp_preds;
          
          
          List updated_test_preds=update_predictions_gs(tree_table,new_node_mean,new_node_var,num_test_obs,term_nodes,term_test_obs);         
          NumericVector temp_test_preds=updated_test_preds[1];
          sum_new_test_predictions(_,k)=temp_test_preds;
          
        }
        
        
        NumericVector pred_obs=calc_rowsums(sum_new_predictions);
        NumericVector full_resids = y_scaled-pred_obs;
        //get overall predictions for current iteration and current sum of trees
        sigma2= update_sigma(a1,b,full_resids,num_obs);
        sigma_its[j]=sigma2;
        
        post_predictions(j,_)=pred_obs;
        NumericVector original_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_obs);    
        post_predictions_orig(j,_)=original_y;
        
        NumericVector pred_test_obs=calc_rowsums(sum_new_test_predictions);
        post_test_predictions(j,_)=pred_test_obs;
        NumericVector original_test_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_test_obs);    
        post_test_predictions_orig(j,_)=original_test_y;
        //prediction interval for training 
        
        for(int g=0;g<pred_obs.size();g++){
          //post_predictions_PI(j,g)=pred_obs[g];
          post_predictions_PI(j,g)=R::rnorm(pred_obs[g],sigma2);
        }
        
        NumericVector original_y_PI=get_original_gs(min(y),max(y),-0.5,0.5,post_predictions_PI(j,_));    
        post_predictions_orig_PI(j,_)=original_y_PI;
        
      }
      prediction_list[i]=post_predictions;
      prediction_list_orig[i]=post_predictions_orig;
      prediction_test_list[i]=post_test_predictions;
      prediction_test_list_orig[i]=post_test_predictions_orig;
      
      //predictive intervals
      
      predictive_dist_train_list[i]=post_predictions_PI;
      predictive_dist_train_list_orig[i]= post_predictions_orig_PI;
      
      
    }else{
      
      NumericVector sigma_its(num_iter);
      NumericMatrix post_predictions(num_iter,num_obs);
      NumericMatrix post_predictions_PI(num_iter,num_obs);
      NumericMatrix post_predictions_orig_PI(num_iter,num_obs);
      NumericMatrix post_predictions_orig(num_iter,num_obs);
      NumericMatrix post_test_predictions(num_iter,num_obs);
      NumericMatrix post_test_predictions_orig(num_iter,num_obs);
      NumericMatrix sum_predictions(num_obs,1);
      NumericMatrix sum_test_predictions(num_test_obs,1);
      
      NumericVector predictions=resids[i];
      
      for(int j=0;j<num_iter;j++){
        NumericMatrix tree_table=overall_sum_trees1[i];
        IntegerMatrix tree_mat=overall_sum_mat1[i];
        IntegerVector term_nodes=overall_term_nodes_trees[i];
        List term_obs=overall_term_obs_trees[i];
        List term_test_obs=overall_term_test_obs_trees[i];
        //find terminal node means and observations associated with them
        
        //current predictions are the residuals for sum of trees
        
        //update the means and predictions for tree
        List new_node_mean_var=update_Gibbs_mean_var(tree_table,predictions,a,sigma2,mu_mu,term_nodes,term_obs);
        NumericVector new_node_mean=get_new_mean(term_nodes,new_node_mean_var);      
        NumericVector new_node_var=new_node_mean_var[1];
        List updated_preds=update_predictions_gs(tree_table,new_node_mean,new_node_var,num_obs,term_nodes,term_obs);         
        NumericVector temp_preds=updated_preds[1];
        sum_predictions(_,1)=temp_preds;
        //NOW UPDATE RESIDUALS
        predictions=y_scaled-temp_preds;
        //get updated predictions for the test data
        List updated_test_preds=update_predictions_gs(tree_table,new_node_mean,new_node_var,num_test_obs,term_nodes,term_test_obs);         
        NumericVector temp_test_preds=updated_test_preds[1];
        //sum_predictions(_,1)=temp_preds;
        sum_test_predictions(_,1)=temp_test_preds;
        
        //get overall predictions for current iteration and current sum of trees
        
        
        //NumericVector pred_obs=calc_rowsums(sum_predictions);
        NumericVector pred_obs=temp_preds;
        
        //NumericVector pred_test_obs=calc_rowsums(sum_test_predictions);
        NumericVector pred_test_obs=temp_test_preds;
        
        post_predictions(j,_)=pred_obs;
        post_test_predictions(j,_)=pred_test_obs;
        NumericVector full_resids = y_scaled-pred_obs;
        sigma2= update_sigma(a1,b,full_resids,num_obs);
        sigma_its[j]=sigma2;
        NumericVector original_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_obs);   
        NumericVector original_test_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_test_obs);   
        post_predictions_orig(j,_)=original_y;
        post_test_predictions_orig(j,_)=original_test_y;
        
        //prediction interval for training 
        for(int g=0;g<pred_obs.size();g++){
          //post_predictions_PI(j,g)=pred_obs[g];
          post_predictions_PI(j,g)=R::rnorm(pred_obs[g],sigma2);
        }
        
        NumericVector original_y_PI=get_original_gs(min(y),max(y),-0.5,0.5,post_predictions_PI(j,_));    
        post_predictions_orig_PI(j,_)=original_y_PI;
        
      }
      
      prediction_list[i]=post_predictions;
      prediction_list_orig[i]=post_predictions_orig;
      prediction_test_list[i]=post_test_predictions;
      prediction_test_list_orig[i]=post_test_predictions_orig;
      
      //predictive intervals
      
      predictive_dist_train_list[i]=post_predictions_PI;
      predictive_dist_train_list_orig[i]= post_predictions_orig_PI;
    }
    sigma_chains[i]=sigma_its;
    
  }
  
  
  List ret(7);
  NumericVector test2=sigma_chains[0];
  ret[0]= prediction_list;
  ret[1]= prediction_list_orig;
  ret[2]=sigma_chains;
  ret[3]=prediction_test_list;
  ret[4]=prediction_test_list_orig;
  ret[5]=predictive_dist_train_list;
  ret[6]=predictive_dist_train_list_orig;
  return(ret); 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List gibbs_sampler2(List overall_sum_trees_mu,List overall_sum_trees_tau,
                    List overall_sum_mat_mu,List overall_sum_mat_tau,
                    NumericVector y,NumericVector BIC_weights,
                    int num_iter,int burnin,int num_obs,double a_mu,double a_tau,double sigma,
                    double mu_mu_mu,double mu_mu_tau,double nu,double lambda,
                    List resids_mu,List resids_tau, NumericVector z){
  if(burnin>=num_iter){
    throw std::range_error("Number of iterations has to be greater than the number of burn-in samples");
  }
  
  //get tree info just once don't repeat for each iteration:
  //need terminal nodes, terminal node means and observations terminal nodes refer to
  List tree_info_mu=get_tree_info(overall_sum_trees_mu,overall_sum_mat_mu,num_obs);
  List tree_info_tau=get_tree_info_tau(overall_sum_trees_tau,overall_sum_mat_tau,num_obs,z);
  //List test_tree_info_mu=get_tree_info_testdata_overall(overall_sum_trees_mu,num_test_obs,test_data);
  //List overall_term_test_obs_trees_mu=test_tree_info_mu[1];
  //List test_tree_info_tau=get_tree_info_testdata_overall_tau(overall_sum_trees_tau,num_test_obs,test_data);
  //List overall_term_test_obs_trees_tau=test_tree_info_tau[1];
  
  double a1=nu/2;
  double b=2/(nu*lambda);
  double sigma2=sigma;
  double sigma_init=sigma2;
  
  List overall_term_nodes_trees_mu=tree_info_mu[0];
  List overall_term_obs_trees_mu=tree_info_mu[1];
  List overall_predictions_mu=tree_info_mu[2];
  
  List overall_term_nodes_trees_tau=tree_info_tau[0];
  List overall_term_obs_trees_tau=tree_info_tau[1];
  List overall_predictions_tau=tree_info_tau[2];
  
  
  NumericVector y_scaled=scale_response_gs(min(y),max(y),-0.5,0.5,y); 
  List prediction_list(overall_sum_trees_mu.size());
  List prediction_list_orig(overall_sum_trees_mu.size());
  //List prediction_test_list(overall_sum_trees_mu.size());
  List predictive_dist_train_list(overall_sum_trees_mu.size());
  List predictive_dist_train_list_orig(overall_sum_trees_mu.size());
  
  List TE_list(overall_sum_trees_mu.size());
  List TE_list_orig(overall_sum_trees_mu.size());
  List TE_test_list(overall_sum_trees_mu.size());
  List TE_test_list_orig(overall_sum_trees_mu.size());
  
  //List prediction_test_list_orig(overall_sum_trees_mu.size());
  List overall_sum_trees1_mu=clone(overall_sum_trees_mu);
  List overall_sum_mat1_mu=clone(overall_sum_trees_mu);
  List overall_sum_trees1_tau=clone(overall_sum_trees_mu);
  List overall_sum_mat1_tau=clone(overall_sum_trees_mu);
  
  List sigma_chains(overall_sum_trees_mu.size());
  //int one_tree=0;
  Rcout << "Line 846. \n";
  for(int i=0;i<overall_sum_trees_mu.size();i++){
    //for each set of trees loop over individual trees    
    NumericVector sigma_its(num_iter);
    sigma2=sigma_init;
    SEXP s_mu = overall_sum_trees_mu[i];
    SEXP s_tau = overall_sum_trees_tau[i];
    
    //NumericVector test_preds_sum_tree;
    NumericMatrix sum_predictions_mu;
    NumericMatrix sum_predictions_tau;
    
    //NumericMatrix sum_test_predictions_mu;
    //NumericMatrix sum_test_predictions_tau;
    Rcout << "Line 860. Round  i= " << i << " .\n";
    if(is<List>(s_mu)){
    if(is<List>(s_tau)){
      Rcout << "Line 863. \n";
      //if current set of trees contains more than one tree
      List sum_tree_mu=overall_sum_trees1_mu[i];
      Rcout << "Line 866. \n";
      List sum_tree_mat_mu=overall_sum_mat1_mu[i];
      List sum_term_nodes_mu=overall_term_nodes_trees_mu[i];
      List sum_term_obs_mu=overall_term_obs_trees_mu[i];
      //List sum_term_test_obs_mu=overall_term_test_obs_trees_mu[i];
      List sum_tree_tau=overall_sum_trees1_tau[i];
      List sum_tree_mat_tau=overall_sum_mat1_tau[i];
      List sum_term_nodes_tau=overall_term_nodes_trees_tau[i];
      List sum_term_obs_tau=overall_term_obs_trees_tau[i];
      //List sum_term_test_obs_tau=overall_term_test_obs_trees_tau[i];
      Rcout << "Line 875. \n";
      List sum_resids0_mu=resids_mu[i];         //NOTE: IF s_mu or s_tau NOT LISTS, THEN resids_mu[i] and resids_tau[i] might not be lists... need to check this
      List sum_resids_mu=clone(sum_resids0_mu);
      List sum_resids0_tau=resids_tau[i];
      List sum_resids_tau=clone(sum_resids0_tau);
      Rcout << "Line 879. \n";
      NumericMatrix tree_predictions_mu=overall_predictions_mu[i];
      sum_predictions_mu=clone(tree_predictions_mu);
      NumericMatrix tree_predictions_tau=overall_predictions_tau[i];
      sum_predictions_tau=clone(tree_predictions_tau);
      
      NumericMatrix post_predictions(num_iter,num_obs);
      NumericMatrix post_predictions_PI(num_iter,num_obs);
      NumericMatrix post_predictions_orig_PI(num_iter,num_obs);
      NumericMatrix post_predictions_orig(num_iter,num_obs);
      //NumericMatrix post_test_predictions(num_iter,num_test_obs);
      //NumericMatrix post_test_predictions_orig(num_iter,num_test_obs);
      
      NumericMatrix post_TEs(num_iter,num_obs);
      //NumericMatrix post_test_TEs(num_iter,num_obs);
      NumericMatrix post_TEs_orig(num_iter,num_obs);
      //NumericMatrix post_test_TEs_orig(num_iter,num_obs);
      
      NumericMatrix sum_new_predictions_mu(sum_predictions_mu.nrow(),sum_predictions_mu.ncol());
      //NumericMatrix sum_new_test_predictions_mu(sum_predictions_mu.nrow(),sum_predictions_mu.ncol());
      NumericMatrix sum_new_predictions_tau(sum_predictions_tau.nrow(),sum_predictions_tau.ncol());
      //NumericMatrix sum_new_test_predictions_tau(sum_predictions_tau.nrow(),sum_predictions_tau.ncol());
      
      NumericMatrix sum_new_TEs(sum_predictions_tau.nrow(),sum_predictions_tau.ncol());
      //NumericMatrix sum_new_test_TEs(sum_predictions_tau.nrow(),sum_predictions_tau.ncol());
      
      Rcout << "Line 904. \n";
      for(int j=0;j<num_iter;j++){
        // NOW LOOP OVER MU TREES IN SUM OF TREE MODEL i
        for(int k =0;k<sum_tree_mu.size();k++){
          NumericMatrix tree_table_mu=sum_tree_mu[k];
          IntegerMatrix tree_mat_mu=sum_tree_mat_mu[k];
          //find terminal node means and observations associated with them
          IntegerVector term_nodes_mu=sum_term_nodes_mu[k];
          List term_obs_mu=sum_term_obs_mu[k];
          //  List term_test_obs=sum_term_test_obs[k];
          NumericVector predictions_mu=sum_resids_mu[k];
          //current predictions are the residuals for sum of trees!
          
          //update the means and predictions for tree
          
          List new_node_mean_var=update_Gibbs_mean_var(tree_table_mu,predictions_mu,a_mu,sigma2,mu_mu_mu,term_nodes_mu,term_obs_mu);
          NumericVector new_node_mean=get_new_mean(term_nodes_mu,new_node_mean_var);
          NumericVector new_node_var=new_node_mean_var[1];
          //update predictions by setting predicted value for term_obs[termnode]=new mean value!
          
          List updated_preds_mu=update_predictions_gs(tree_table_mu,new_node_mean,new_node_var,num_obs,term_nodes_mu,term_obs_mu);         
          NumericVector temp_preds_mu=updated_preds_mu[1];
          //NOW UPDATE THE RESIDUALS FOR USE IN NEXT ITERATION
          //THE PLACING OF THIS SECTION OF CODE HERE IS IMPORTANT
          //MUST BE BEFORE sum_new_predictions(_,k) is updated so that the
          //previous round's predictions can be added (i.e. removed from the residual before the new predictions are taken away to create the new residual) 
          Rcout << "Line 930. \n";
          
          for(int l=0;l<sum_tree_mu.size();l++){
            if(l!=k){
              if(j==0){
                NumericVector temp_resids1_mu= sum_resids_mu[l];
                sum_resids_mu[l]=temp_resids1_mu+sum_predictions_mu(_,k)-temp_preds_mu;
              }else{
                NumericVector temp_resids1_mu= sum_resids_mu[l];
                sum_resids_mu[l]=temp_resids1_mu+sum_new_predictions_mu(_,k)-temp_preds_mu;
              }
            }
          }
          for(int l=0;l<sum_tree_tau.size();l++){
              if(j==0){
                NumericVector temp_resids1_tau= sum_resids_tau[l];
                sum_resids_tau[l]=temp_resids1_tau+sum_predictions_mu(_,k)-temp_preds_mu;
              }else{
                NumericVector temp_resids1_tau= sum_resids_tau[l];
                sum_resids_tau[l]=temp_resids1_tau+sum_new_predictions_mu(_,k)-temp_preds_mu;
              }
            }
          sum_new_predictions_mu(_,k)=temp_preds_mu;
          
          Rcout << "Line 954. \n";
          //List updated_test_preds_mu=update_predictions_gs(tree_table_mu,new_node_mean,new_node_var,num_obs,term_nodes_mu,term_test_obs_mu);         
          //NumericVector temp_test_preds_mu=updated_test_preds_mu[1];
          //sum_new_test_predictions_mu(_,k)=temp_test_preds_mu;
          
          
        }
        // NOW LOOP OVER TAU TREES IN SUM OF TREE MODEL i
        for(int k =0;k<sum_tree_tau.size();k++){
          Rcout << "Line 963. \n";
          NumericMatrix tree_table_tau=sum_tree_tau[k];
          Rcout << "Line 968. \n";
          IntegerMatrix tree_mat_tau=sum_tree_mat_tau[k];
          //find terminal node means and observations associated with them
          IntegerVector term_nodes_tau=sum_term_nodes_tau[k];
          List term_obs_tau=sum_term_obs_tau[k];
          //  List term_test_obs_tau=sum_term_test_obs_tau[k];
          NumericVector predictions_tau=sum_resids_tau[k];
          //current predictions are the residuals for sum of trees!
          
          //update the means and predictions for tree
          Rcout << "Line 978. \n";
          List new_node_mean_var=update_Gibbs_mean_var(tree_table_tau,predictions_tau,a_tau,sigma2,mu_mu_tau,term_nodes_tau,term_obs_tau);
          NumericVector new_node_mean=get_new_mean(term_nodes_tau,new_node_mean_var);
          NumericVector new_node_var=new_node_mean_var[1];
          //update predictions by setting predicted value for term_obs[termnode]=new mean value!
          Rcout << "Line 983. \n";
          List updated_preds_tau=update_predictions_gs(tree_table_tau,new_node_mean,new_node_var,num_obs,term_nodes_tau,term_obs_tau);         
          Rcout << "Line 985. \n";
          NumericVector temp_preds_tau=updated_preds_tau[1];
          Rcout << "Line 987. \n";
          NumericVector temp_preds_tau_z=temp_preds_tau*z;
          Rcout << "Line 989. \n";
          //NOW UPDATE THE RESIDUALS FOR USE IN NEXT ITERATION
          //THE PLACING OF THIS SECTION OF CODE HERE IS IMPORTANT
          //MUST BE BEFORE sum_new_predictions(_,k) is updated so that the
          //previous round's predictions can be added (i.e. removed from the residual before the new predictions are taken away to create the new residual) 
          for(int l=0;l<sum_tree_mu.size();l++){
              if(j==0){
                NumericVector temp_resids1_mu= sum_resids_mu[l];
                sum_resids_mu[l]=temp_resids1_mu+z*sum_predictions_tau(_,k)-temp_preds_tau_z;
              }else{
                NumericVector temp_resids1_mu= sum_resids_mu[l];
                sum_resids_mu[l]=temp_resids1_mu+sum_new_predictions_tau(_,k)-temp_preds_tau_z;
              }
            }
          
          for(int l=0;l<sum_tree_tau.size();l++){
            if(l!=k){
              if(j==0){
                NumericVector temp_resids1_tau= sum_resids_tau[l];
                sum_resids_tau[l]=temp_resids1_tau+z*sum_predictions_tau(_,k)-temp_preds_tau_z;
              }else{
                NumericVector temp_resids1_tau= sum_resids_tau[l];
                sum_resids_tau[l]=temp_resids1_tau+sum_new_predictions_tau(_,k)-temp_preds_tau_z;
              }
            }
          }
          Rcout << "Line 1015. \n";
          sum_new_predictions_tau(_,k)=temp_preds_tau_z;
          sum_new_TEs(_,k)=temp_preds_tau;
          Rcout << "Line 1018. \n";
          
          //List updated_test_preds_tau=update_predictions_gs(tree_table_tau,new_node_mean,new_node_var,num_obs,term_nodes_tau,term_test_obs_tau);         
          //NumericVector temp_test_preds_tau=updated_test_preds_tau[1];
          //sum_new_test_predictions_tau(_,k)=z*temp_test_preds_tau;
          //sum_new_test_TEs(_,k)=temp_test_preds_tau;
          
          
        }
        Rcout << "Line 1027. \n";
        
        NumericVector pred_obs_mu=calc_rowsums(sum_new_predictions_mu);
        NumericVector pred_obs_tau=calc_rowsums(sum_new_predictions_tau);
        NumericVector pred_TEs=calc_rowsums(sum_new_TEs);
        //NumericVector pred_test_TEs=calc_rowsums(sum_new_test_TEs);
        Rcout << "Line 1033. \n";
        
        NumericVector pred_obs=pred_obs_mu+pred_obs_tau;
        NumericVector full_resids = y_scaled-pred_obs;
        //get overall predictions for current iteration and current sum of trees
        sigma2= update_sigma(a1,b,full_resids,num_obs);
        sigma_its[j]=sigma2;  
        
        post_predictions(j,_)=pred_obs;
        post_TEs(j,_)=pred_TEs;
        //post_test_TEs(j,_)=pred_test_TEs;
        
        NumericVector original_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_obs);    
        post_predictions_orig(j,_)=original_y;
        
        NumericVector original_TEs=get_original_gs(min(y),max(y),-0.5,0.5,pred_TEs);    
        post_TEs_orig(j,_)=original_TEs;
        //NumericVector original_test_TEs=get_original_gs(min(y),max(y),-0.5,0.5,pred_test_TEs);    
        //post_test_TEs_orig(j,_)=original_TEs;
        
        //NumericVector pred_test_obs_mu=calc_rowsums(sum_new_test_predictions_mu);
        //NumericVector pred_test_obs_tau=calc_rowsums(sum_new_test_predictions_tau);
        //NumericVector pred_test_obs=pred_test_obs_mu+pred_test_obs_tau;
        
        //post_test_predictions(j,_)=pred_test_obs;
        //NumericVector original_test_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_test_obs);    
        //post_test_predictions_orig(j,_)=original_test_y;
        //prediction interval for training 
        Rcout << "Line 1055. \n";
        
        for(int g=0;g<pred_obs.size();g++){
          post_predictions_PI(j,g)=R::rnorm(pred_obs[g],sigma2);
        }
        NumericVector original_y_PI=get_original_gs(min(y),max(y),-0.5,0.5,post_predictions_PI(j,_));    
        post_predictions_orig_PI(j,_)=original_y_PI;
        
      }
      Rcout << "Line 1064. \n";
      
      prediction_list[i]=post_predictions;
      prediction_list_orig[i]=post_predictions_orig;
      //prediction_test_list[i]=post_test_predictions;
      //prediction_test_list_orig[i]=post_test_predictions_orig;
      
      //predictive intervals
      predictive_dist_train_list[i]=post_predictions_PI;
      predictive_dist_train_list_orig[i]= post_predictions_orig_PI;
      
      TE_list[i]=post_TEs;
      TE_list_orig[i]=post_TEs_orig;
      //TE_test_list[i]=post_test_TEs;
      //TE_test_list_orig[i]=post_test_TEs_orig;
      
      
      Rcout << "Line 1081. \n";
      
    }else{
      
      throw std::range_error("MU IS A LIST AND TAU IS A MATRIX");
      
    }
    }else{//WHEN NOT A LIST OF MU TREES, i.e. there is just one mu tree
      if(is<List>(s_tau)){
        throw std::range_error("MU IS A MATRIX AND TAU IS A LIST");
        
      }else{
        Rcout << "Line 1093. \n";
        
        NumericVector sigma_its(num_iter);
        NumericMatrix post_predictions(num_iter,num_obs);
        NumericMatrix post_predictions_PI(num_iter,num_obs);
        NumericMatrix post_predictions_orig_PI(num_iter,num_obs);
        NumericMatrix post_predictions_orig(num_iter,num_obs);
        //NumericMatrix post_test_predictions(num_iter,num_obs);
        //NumericMatrix post_test_predictions_orig(num_iter,num_obs);
        NumericMatrix post_TEs(num_iter,num_obs);
        //NumericMatrix post_test_TEs(num_iter,num_obs);
        NumericMatrix post_TEs_orig(num_iter,num_obs);
        //NumericMatrix post_test_TEs_orig(num_iter,num_obs);
        
        NumericMatrix sum_predictions_mu(num_obs,1);
        //NumericMatrix sum_test_predictions_mu(num_test_obs,1);
        NumericMatrix sum_predictions_tau(num_obs,1);
        //NumericMatrix sum_test_predictions_taau(num_test_obs,1);
        
        NumericMatrix sum_new_predictions_mu(num_obs,1);
        //NumericMatrix sum_new_test_predictions_mu(num_test_obs,1);
        NumericMatrix sum_new_predictions_tau(num_obs,1);
        //NumericMatrix sum_new_test_predictions_tau(num_test_obs,1);
        
        NumericMatrix sum_new_TEs(num_obs,1);
        //NumericMatrix sum_new_test_TEs(num_test_obs,1);
        
        NumericVector predictions_mu=resids_mu[i];
        NumericVector predictions_tau=resids_tau[i];
        NumericMatrix tree_predictions_mu=overall_predictions_mu[i];
        NumericMatrix tree_predictions_tau=overall_predictions_tau[i];
        Rcout << "Line 1124. \n";
        
        for(int j=0;j<num_iter;j++){
          Rcout << "Line 1127. \n";
          
          NumericMatrix tree_table_mu=overall_sum_trees1_mu[i];
          IntegerMatrix tree_mat_mu=overall_sum_mat1_mu[i];
          IntegerVector term_nodes_mu=overall_term_nodes_trees_mu[i];
          List term_obs_mu=overall_term_obs_trees_mu[i];
          //  List term_test_obs_mu=overall_term_test_obs_trees_mu[i];
          //find terminal node means and observations associated with them
          
          //current predictions are the residuals for sum of trees
          
          //update the means and predictions for tree
          List new_node_mean_var_mu=update_Gibbs_mean_var(tree_table_mu,predictions_mu,a_mu,sigma2,mu_mu_mu,term_nodes_mu,term_obs_mu);
          NumericVector new_node_mean_mu=get_new_mean(term_nodes_mu,new_node_mean_var_mu);      
          NumericVector new_node_var_mu=new_node_mean_var_mu[1];
          List updated_preds_mu=update_predictions_gs(tree_table_mu,new_node_mean_mu,new_node_var_mu,num_obs,term_nodes_mu,term_obs_mu);         
          NumericVector temp_preds_mu=updated_preds_mu[1];
          sum_predictions_mu(_,1)=temp_preds_mu;
          //NOW UPDATE PARTIAL RESIDUALS
          //MU IS UNAFFECTED, BUT TAU SHOULD BE AFFECTED
          if(j==0){
            NumericVector temp_resids1_tau= predictions_tau;
            predictions_tau=temp_resids1_tau+sum_predictions_mu(_,1)-temp_preds_mu;
          }else{
            NumericVector temp_resids1_tau= predictions_tau;
            predictions_tau=temp_resids1_tau+sum_new_predictions_mu(_,1)-temp_preds_mu;
          }
          sum_new_predictions_mu(_,1)=temp_preds_mu;
          Rcout << "Line 1155. \n";
          
          //get updated predictions for the test data
          //List updated_test_preds_mu=update_predictions_gs(tree_table_mu,new_node_mean_mu,new_node_var_mu,num_test_obs,term_nodes_mu,term_test_obs_mu);         
          //NumericVector temp_test_preds_mu=updated_test_preds_mu[1];
          //sum_test_predictions_mu(_,1)=z*temp_test_preds_mu;
          //sum_new_test_TEs(_,1)=temp_test_preds_tau;
          
          
          NumericMatrix tree_table_tau=overall_sum_trees1_tau[i];
          IntegerMatrix tree_mat_tau=overall_sum_mat1_tau[i];
          IntegerVector term_nodes_tau=overall_term_nodes_trees_tau[i];
          List term_obs_tau=overall_term_obs_trees_tau[i];
          //  List term_test_obs_tau=overall_term_test_obs_trees_tau[i];
          //find terminal node means and observations associated with them
          
          //current predictions are the residuals for sum of trees
          
          //update the means and predictions for tree
          List new_node_mean_var_tau=update_Gibbs_mean_var(tree_table_tau,predictions_tau,a_tau,sigma2,mu_mu_tau,term_nodes_tau,term_obs_tau);
          NumericVector new_node_mean_tau=get_new_mean(term_nodes_tau,new_node_mean_var_tau);      
          NumericVector new_node_var_tau=new_node_mean_var_tau[1];
          List updated_preds_tau=update_predictions_gs(tree_table_tau,new_node_mean_tau,new_node_var_tau,num_obs,term_nodes_tau,term_obs_tau);         
          NumericVector temp_preds_tau=updated_preds_tau[1];
          NumericVector temp_preds_tau_z=z*temp_preds_tau;
          
          sum_predictions_tau(_,1)=temp_preds_tau;
          //NOW UPDATE PARTIAL RESIDUALS
          //TAU IS UNAFFECTED, BUT MU SHOULD BE AFFECTED
          if(j==0){
            NumericVector temp_resids1_mu= predictions_mu;
            predictions_mu=temp_resids1_mu+z*sum_predictions_tau(_,1)-temp_preds_tau_z;
          }else{
            NumericVector temp_resids1_mu= predictions_mu;
            predictions_mu=temp_resids1_mu+sum_new_predictions_tau(_,1)-temp_preds_tau_z;
          }
          sum_new_predictions_tau(_,1)=temp_preds_tau_z;
          sum_new_TEs(_,1)=temp_preds_tau;
          
          //get updated predictions for the test data
          //List updated_test_preds_tau=update_predictions_gs(tree_table_tau,new_node_mean_tau,new_node_var_tau,num_test_obs,term_nodes_tau,term_test_obs_tau);         
          //NumericVector temp_test_preds_tau=updated_test_preds_tau[1];
          //sum_test_predictions_tau(_,1)=temp_test_preds_tau;
          
          
          Rcout << "Line 1200. \n";
          
          
          //get overall predictions for current iteration and current sum of trees
          NumericVector pred_obs_mu=temp_preds_mu;
          NumericVector pred_obs_tau=temp_preds_tau;
          NumericVector pred_obs=pred_obs_mu+pred_obs_tau;
          //NumericVector pred_test_obs_mu=temp_test_preds_mu;
          //NumericVector pred_test_obs_tau=temp_test_preds_tau;
          //NumericVector pred_test_obs=pred_test_obs_mu+pred_test_obs_tau;
          
          post_predictions(j,_)=pred_obs;
          post_TEs(j,_)=temp_preds_tau;
          
          //post_test_predictions(j,_)=pred_test_obs;
          //post_test_TEs(j,_)=temp_test_preds_tau;
          
          NumericVector full_resids=y_scaled-pred_obs;
          sigma2= update_sigma(a1,b,full_resids,num_obs);
          sigma_its[j]=sigma2;
          
          
          NumericVector original_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_obs);   
          //NumericVector original_test_y=get_original_gs(min(y),max(y),-0.5,0.5,pred_test_obs);   
          post_predictions_orig(j,_)=original_y;
          //post_test_predictions_orig(j,_)=original_test_y;
          
          NumericVector original_TEs=get_original_gs(min(y),max(y),-0.5,0.5,temp_preds_tau);    
          post_TEs_orig(j,_)=original_TEs;
          //NumericVector original_test_TEs=get_original_gs(min(y),max(y),-0.5,0.5,temp_test_preds_tau);    
          //post_test_TEs_orig(j,_)=original_TEs;
          
          
          //prediction interval for training 
          for(int g=0;g<pred_obs.size();g++){
            post_predictions_PI(j,g)=R::rnorm(pred_obs[g],sigma2);
          }
          
          NumericVector original_y_PI=get_original_gs(min(y),max(y),-0.5,0.5,post_predictions_PI(j,_));    
          post_predictions_orig_PI(j,_)=original_y_PI;
          
        }
        
        Rcout << "Line 1243. \n";
        
        prediction_list[i]=post_predictions;
        prediction_list_orig[i]=post_predictions_orig;
        //prediction_test_list[i]=post_test_predictions;
        //prediction_test_list_orig[i]=post_test_predictions_orig;
        
        //predictive intervals
        predictive_dist_train_list[i]=post_predictions_PI;
        predictive_dist_train_list_orig[i]= post_predictions_orig_PI;
        
        TE_list[i]=post_TEs;
        TE_list_orig[i]=post_TEs_orig;
        //TE_test_list[i]=post_test_TEs;
        //TE_test_list_orig[i]=post_test_TEs_orig;
        Rcout << "Line 1258. \n";
        
      }
    }
    
    Rcout << "Line 1263. \n";
    
    sigma_chains[i]=sigma_its;
    
  }
  Rcout << "Line 1268. \n";
  
  List ret(9);
  NumericVector test2=sigma_chains[0];
  ret[0]= prediction_list;
  ret[1]= prediction_list_orig;
  ret[2]=sigma_chains;
  //ret[3]=prediction_test_list;
  //ret[4]=prediction_test_list_orig;
  ret[3]=predictive_dist_train_list;
  ret[4]=predictive_dist_train_list_orig;
  ret[5]=TE_list;
  ret[6]=TE_list_orig;
  ret[7]=TE_test_list;
  ret[8]=TE_test_list_orig;
  return(ret); 
}