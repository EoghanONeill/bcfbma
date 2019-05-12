
//######################################################################################################################//

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix add_rows_bcf(NumericMatrix prior_tree_table_temp,int grow_node){  // function add_rows_bcf takes matrix and integer as arguments
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);  // makes a arma matrix M from matrix argument. Not sure if as<>() (R to C++?) is necessary? Not sure if "Rcpp::" because using rcpp namesapce.
  M(grow_node-1,5)=0;		// Assign value 0 to matrix M row grow_node-1, column 5
  M(grow_node-1,6)=0;		// Assign value 0 to matrix M row grow_node-1, column 6
  M(grow_node-1,0)=grow_node+1;		// Assign value grow_node+1 to matrix M row grow_node-1, column 0 ?? Matrix indices appear to start from 0.
  M(grow_node-1,1)=grow_node+2;		// Assign value grow_node+2 to matrix M row grow_node-1, column 1 ??
  M.insert_rows(grow_node,2);			// Insert two (second argument) rows in M at (not sure if after or before) row number grow_node. Rows are set to zero by default.
  M(grow_node,4)=-1;					// Insert value -1 to matrix M row number grow_node, column 4
  M(grow_node+1,4)=-1;				// Insert value -1 to matrix M row number grow_node+1, column 4
  NumericMatrix t=as<NumericMatrix>(wrap(M));		// Now convert M to a Numeric martrix (from arma mat) and call it t. ?? First use wrap to convert to R object then use as to convert back to C++?
  IntegerVector rname=seq_len(t.nrow());			// Create a integer vector (1, 2, ..., number_rows_in_t) ??. seq_len (part of rccp) "creates an integer sugar expression whose ith element expnds to i". 
  
  List dimnms = // two vec. with static names		// creates list. I'm not sure what this line does. List is an rcpp type
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  t.attr("dimnames") = dimnms;		// R objects have attributes. This line sets the attribute of t called "dimnames" to be equal to the list dimnms. 
  
  return(t);		// return the new matrix
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix addcol_bcf(NumericMatrix prior_tree_matrix_temp,int grow_node,NumericVector ld_obs,NumericVector rd_obs){		// create a function called addcol_bcf. take a matrix, integer, and two vectors as arguments.
  int ncol=prior_tree_matrix_temp.ncol();		// integer named ncol is equal to the number of columns in the input matrix
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);	// input matrix converted to a arma mat M
  M.insert_cols(ncol,1);		// insert one column of zeros after the final column of M (note: prviously columns numbered 0 to ncol-1, so added col number ncol).
  for(int i =0;i<ld_obs.size();i++){		// for loop, length of input vectot ld_obs
    try{
      if(ld_obs[i]>prior_tree_matrix_temp.nrow()){		// create an error if the value of entry i in input vector ld_obs is greater than the number of rows in the input matrix
        throw std::range_error("can't add col because ld row index is out of range");
      }
    }catch(...){
      ::Rf_error("there is a problem adding col to mat don't know why");
    }
    M(ld_obs[i],ncol)=grow_node+1;		// Enter value grow_node+1 in matrix M row ld_obs[i], column ncol (i.e. the last column)
  }
  for(int i =0;i<rd_obs.size();i++){		// for loop, length of input vectot rd_obs
    try{
      if(rd_obs[i]>prior_tree_matrix_temp.nrow()){	// create an error if the value of entry i in input vector rd_obs is greater than the number of rows in the input matrix
        throw std::range_error("can't add col because rd row index is out of range");
      }
    }catch(...){
      ::Rf_error("there is a problem adding rd col to mat");
    }    
    M(rd_obs[i],ncol)=grow_node+2;		// Enter value grow_node+2 in matrix M row rd_obs[i], column ncol (i.e. the new last column)
  }
  return(wrap(M));		// return the matrix
} 

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter_to_end_tree_bcf(int grow_node,NumericMatrix prior_tree_table_temp,double left_daughter){	// create a function with integer, matrix, double input
  int nrow=prior_tree_table_temp.nrow();						// number of rows in input matrix
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);		// create a matrix M equal to the input matrix
  M(grow_node-1,5)=0;											// Assign value 0 to matrix M row grow_node-1, column 5
  M(grow_node-1,6)=0;											// Assign value 0 to matrix M row grow_node-1, column 6
  M.insert_rows(nrow,2);										// create two rows of zeros after final row of M 
  M(grow_node-1,0)=left_daughter;								// Assign value left_daughter to matrix M, row number grow_node-1, and column number 0.
  M(grow_node-1,1)=left_daughter+1;							// Assign value left_daughter+1 to matrix M, row number grow_node-1, and column number 1.
  M(left_daughter-1,4)=-1;									// Assign value -1 to Matrix M row left_daughter-1, column 4
  M(left_daughter,4)=-1;										// Assign value -1 to Matrix M row left_daughter, column 4
  
  NumericMatrix s=as<NumericMatrix>(wrap(M));					// convert M to numeric matrix called s
  IntegerVector rname=seq_len(s.nrow());						// rname is vector 1,2,3,..., up to the number of rows in the matrix s
  
  List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  s.attr("dimnames") = dimnms;								// add attributes "dimnames" to the R object s.
  
  return(s);	// return the matrix
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter_to_end_mat_bcf(double d,NumericMatrix prior_tree_matrix_temp,double left_daughter,NumericVector ld_obs,NumericVector rd_obs){
  int ncol_mat=prior_tree_matrix_temp.ncol();		// number of columns in input matrix
  arma::mat N=Rcpp::as<arma::mat>(prior_tree_matrix_temp);	// saves input matrix as arma mat N
  arma::vec colmat=N.col(d);			// extracts column number d+1 (why is d a double rather than an integer?) ASSUMING CHANGES DON'T PROPAGATE TO N
  NumericVector colmat2=wrap(colmat);		// convert from arma vec to numeric vector 
  
  if(d+1==ncol_mat){		// if d+1 equls the number of columns
    //update prior_tree_matrix
    //insert extra column for the new split node
    
    N.insert_cols(ncol_mat,1);
    int nrow_mat=prior_tree_matrix_temp.nrow();
    NumericVector colmatzero(nrow_mat);
    colmatzero[ld_obs]=left_daughter;
    colmatzero[rd_obs]=left_daughter+1;
    colmat=Rcpp::as<arma::vec>(colmatzero);
    N.col(d+1)=colmat;			// put this new column in as column d+2 (index starts at 0 so there are now d+2 columns) of N. this replaces the last column of zeros (note the if condition).
    
  }else{
    colmat2[ld_obs]=left_daughter;		// change ld_obs th element to left_daughter. (or vector of obs indexed by ld_obs)
    colmat2[rd_obs]=left_daughter+1;	// change rd_obs th element to left_daughter+1. (or vector of obs indexed by rd_obs)
    colmat=Rcpp::as<arma::vec>(colmat2);	// convert from numeric vector to arma vec and write over colmat
    N.col(d)=colmat;  // write over column d+1 (index starts at 0 so there are d+1 columns). (note that in this case d+1 is not equal to the number of columns in N).
  }
  
  return(wrap(N));		// convert N to R object and return the matrix
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector remove_zero_bcf(NumericVector nodes_at_depth){
  arma::vec nodes_at_depth2=Rcpp::as<arma::vec>(nodes_at_depth);		// converts to arma vec
  arma::vec ret=nodes_at_depth2.elem(arma::find(nodes_at_depth2!=0));		// Defines a arma vec without all the zero values [forst finds indices, then takes vector of elements defined by indices]
  return(wrap(ret));		// converts to R  object and returns
}

//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector order_intvec_bcf(IntegerVector x) {		// input is a vector of integers.  Output is vector of numbers with 1 indicating the position of the highest input, 2 the second highest, and so on.
  IntegerVector sorted = clone(x).sort();		// sorts ascending?? 
  std::reverse(sorted.begin(), sorted.end());		// reverses sorted. Could have obtained this in one line? std::sort(clone(x).begin(), clone(x).end(), std::greater<>())
  
  return match(sorted, x); //  match is the Rcpp sugar version of the R function match, which returns a vector of the positions of the first matches of the first argument in the second.
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector get_gnp_bcf(NumericVector nodes_at_depth,int grow_node){		
  arma::uvec grow_node_pos=arma::find(as<arma::vec>(nodes_at_depth)==grow_node);	// uvec is vector of unsigned integers for index. Set equal to the index number for member(members?) of nodes_at_depth equal to grow_node
  
  return(wrap(grow_node_pos));  // converts from C++ object ro R object and returns the vector of indices
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_term_nodes_bcf(NumericMatrix tree_table){
  arma::mat arma_tree(tree_table.begin(),tree_table.nrow(), tree_table.ncol(), false);	// define a arma mat called arma_tree. Initialize using the dimensions and point to beginning of the data. No explicit memory allocation is needed.
  arma::vec colmat=arma_tree.col(4);		// Defines colmat as the fifth column of arma_tree (index starts at 0)
  arma::uvec term_nodes=arma::find(colmat==-1);		// vector of indices giving positions of values in colmat that are equal to -1
  term_nodes=term_nodes+1;		// Add one to all elements of term_nodes
  
  return(wrap(term_nodes));		// convert term_nodes to an R object and return
} 

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_term_obs_bcf(NumericMatrix tree_matrix_temp,double terminal_node){
  arma::mat arma_tree_mat(tree_matrix_temp.begin(),tree_matrix_temp.nrow(), tree_matrix_temp.ncol(), false); // read tree_matrix_temp as an arma mat call it arma_tree_mat
  arma::uvec term_obs;	// create an unsigned integer vector called term_obs (for indices?)
  
  for(int j=0;j<tree_matrix_temp.ncol();j++){	// loop over number of columns
    arma::vec colmat=arma_tree_mat.col(j);	// let colmat equal the j+1th column of arma_tree_mat. Note index starts at 0.
    term_obs=arma::find(colmat==terminal_node);	// vector of indices giving positions of colmat elements that are equal to terminal_node
    if(term_obs.size()>0){	// If the vector of indices is of non-empty size. i.e. if any elements of the j+1th column are equal to terminal_node
      break;		// breaks out of loop? Therefore term_obs will be the first nonempty set of indices (if any). i.e. indices for the leftmost column that has elements equal to terminal_node
    }
  }
  
  return(term_obs);	// Returns the vector of indices. (could be empty?). Not converted to r object by wrap()?
}

//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double likelihood_function_mu_bcf(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
  double tree_log_lik;		// create a double (not initialized)
  NumericVector terminal_nodes=find_term_nodes_bcf(treetable_temp);		// find term nodes function defined line 168.  gives index of values of treetable_temp that are term nodes (indices from 1 to length of vector). Why not integer vector?
  double b=terminal_nodes.size();		// length of term nodes vector
  IntegerVector n(b);		// integer vector called n of length b. (to be filled with numbers of obs in each terminal node?)
  double term1=0;		// initialize a double 
  double term2=0;		// initialize a double
  arma::uvec term_obs;		// create an unsigned integer vector term_obs (not initialized)
  int ni;		// create an integer called ni (not initialized)
  arma::vec y_temp2=as<arma::vec>(y_temp);		// convert y_temp input to an arma vec called y_temp2
  for(int i=0;i< b;i++){		// loop over the length of the terminal nodes vector
    term_obs=find_term_obs_bcf(obs_to_nodes_temp,terminal_nodes[i]);	// use the find_term_obs_bcf function defined above. Returns a vector of indices
    arma::vec y_k=y_temp2.elem(term_obs);	// extracts the elements of y_temp corresponding to the term_obs indices (presumably the observations for individuals falling in the relevant node)
    ni=term_obs.size();		// length of term_obs. Presumably the number of individuals falling in the relevant node
    n[i]=ni;		// save ni in vector n. Presumably n gives numbers of obs in each terminal node
    double ybar=0;		// initialize a double called ybar
    if(y_k.size()!=0){ 	// If there is a nonzero number of y observations (for the relevant terminal node?)
      ybar=mean(y_k); 	// Then let ybar be the mean of the observations (in the node)
    }
    term1+=log(ni+a);		// iteratively add the log of ni+a (get the log of the products of ni+a for all i)
    arma::vec y_k_sq=pow(y_k,2);	// let y_k_sq be the squares of all the obsrevations in y_k (observations in the terminal node)
    double sum_yksq=sum(y_k_sq);	// let sum_yksq be the sum of these squares
    double b2=pow(ni*ybar +a*mu,2)/(ni+a);		//  b2 is a square of an expression all divided by (ni+a). 
    term2+=(sum_yksq+a*pow(mu,2)-b2+nu*lambda);	// iteratively add function of parameters and sum of squares
  }
  tree_log_lik=(b/2)*log(a)-0.5*term1-((y_temp.size()+nu)/2)*log(term2);	// function of paramters, number of observations, and term1 and term2 obtained in for loop
  
  // for(int i=0;i<b;i++){		// loop over the length of the terminal nodes vector
  //   if(n[i]<=5){	// If the number of obs in term node i is less than 5
  //     tree_log_lik=tree_log_lik-(100000);	// Then take 100000 from the tree_log_lik. (Should this 100000 be entered as a constant above?). Presumably it is an arbitrary large number
  //   }
  //   else{
  //     tree_log_lik=tree_log_lik;	// otherwise keep the tree_log_lik the same. Over the loop, loglik is very negative if any term node has less than 5 obs
  //   }
  // }
  
  return(tree_log_lik);		// returns the tree_log_lik
}  


//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double likelihood_function_bcf(NumericVector y_temp,NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
  double tree_log_lik;		// create a double (not initialized)
  NumericVector terminal_nodes=find_term_nodes_bcf(treetable_temp);		// find term nodes function defined line 168.  gives index of values of treetable_temp that are term nodes (indices from 1 to length of vector). Why not integer vector?
  double b=terminal_nodes.size();		// length of term nodes vector
  IntegerVector n(b);		// integer vector called n of length b. (to be filled with numbers of obs in each terminal node?)
  double term1=0;		// initialize a double 
  double term2=0;		// initialize a double
  arma::uvec term_obs;		// create an unsigned integer vector term_obs (not initialized)
  int ni;		// create an integer called ni (not initialized)
  arma::vec y_temp2=as<arma::vec>(y_temp);		// convert y_temp input to an arma vec called y_temp2
  for(int i=0;i< b;i++){		// loop over the length of the terminal nodes vector
    term_obs=find_term_obs_bcf(obs_to_nodes_temp,terminal_nodes[i]);	// use the find_term_obs_bcf function defined above. Returns a vector of indices
    arma::vec y_k=y_temp2.elem(term_obs);	// extracts the elements of y_temp corresponding to the term_obs indices (presumably the observations for individuals falling in the relevant node)
    ni=term_obs.size();		// length of term_obs. Presumably the number of individuals falling in the relevant node
    n[i]=ni;		// save ni in vector n. Presumably n gives numbers of obs in each terminal node
    double ybar=0;		// initialize a double called ybar
    if(y_k.size()!=0){ 	// If there is a nonzero number of y observations (for the relevant terminal node?)
      ybar=mean(y_k); 	// Then ley ybar be the mean of the observations (in the node)
    }
    term1+=log(ni+a);		// iteratively add the log of ni+a (get the log of the products of ni+a for all i)
    arma::vec y_k_sq=pow(y_k,2);	// let y_k_sq be the squares of all the obsrevations in y_k (observations in the terminal node)
    double sum_yksq=sum(y_k_sq);	// let sum_yksq be the sum of these squares
    double b2=pow(ni*ybar +a*mu,2)/(ni+a);		//  b2 is a square of an expression all divided by (ni+a). 
    term2+=(sum_yksq+a*pow(mu,2)-b2+nu*lambda);	// iteratively add function of parameters and sum of squares
  }
  tree_log_lik=(b/2)*log(a)-0.5*term1-((y_temp.size()+nu)/2)*log(term2);	// function of paramters, number of observations, and term1 and term2 obtained in for loop
  
  // for(int i=0;i<b;i++){		// loop over the length of the terminal nodes vector
  //   if(n[i]<=5){	// If the number of obs in term node i is less than 5
  //     tree_log_lik=tree_log_lik-(100000);	// Then take 100000 from the tree_log_lik. (Should this 100000 be entered as a constant above?). Presumably it is an arbitrary large number
  //   }
  //   else{
  //     tree_log_lik=tree_log_lik;	// otherwise keep the tree_log_lik the same. Over the loop, loglik is very negative if any term node has less than 5 obs
  //   }
  // }
  
  return(tree_log_lik);		// returns the tree_log_lik
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_internal_nodes_bcf(NumericMatrix treetable){
  
  arma::mat arma_tree(treetable.begin(),treetable.nrow(), treetable.ncol(), false); 	// saves the treetable as an arma mat called arma_tree
  arma::vec colmat=arma_tree.col(4);		// defines an arma vec called colmat equal to the fifth column of the tree table
  arma::uvec term_nodes=arma::find(colmat==1);		// index vector of the positions of the elements of colmat that are equal to 1 (presumably 1 deotes an internal node)
  term_nodes=term_nodes+1;		// add one to all the elements of the index vector
  
  return(term_nodes);		// return the index vector (not converted to R object by wrap())
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double find_prev_nonterm_bcf(arma::uvec find_nonterm,NumericVector prev){
  
  double ret=0;		// initialize a double called ret. Why is it not an integer?
  int z=prev.size();		// define z equal to the length of the input vector prev
  for(int j=0;j<z;j++){		// loop the length of prev
    arma::uvec term_equal = arma::find(find_nonterm==prev[j]);		// term_equal is the vector of indices of all elements of find_term that are equal to prev [j]
    ret+=term_equal.size(); // iteratively add the numbers of elements of find_nonterm equal to each element of prev[j]
  }  
  
  return(ret);	// return the total number of elements of find_nonterm equal to elements of prev. (Includes double counting if possible).
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::uvec find_nodes_to_update_bcf(arma::uvec all_ld,double left_daughter){
  arma::uvec gr_ld=arma::find(all_ld>=left_daughter);  		// gr_ld is a vector of indices of all the elements of all_ld that are >= left_daughter
  return(gr_ld);		// returns the vector of indices
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_tree_to_middle_bcf(NumericVector node_to_update,NumericMatrix prior_tree_table_temp,int grow_node,double left_daughter){
  for(int i=0;i<node_to_update.size();i++){		// loop the length of node_to_update
    if((prior_tree_table_temp(node_to_update[i],0)!=0) && (prior_tree_table_temp(node_to_update[i],1)!=0)){
      prior_tree_table_temp(node_to_update[i],0)+=2;		// Add 2 to 1st column row number node_to_update[i] (if nonzero?)
      prior_tree_table_temp(node_to_update[i],1)+=2;		// Add 2 to 2nd column row number node_to_update[i] (if nonzero?)
    }
  }
  
  prior_tree_table_temp(grow_node-1,5)=0;		// set grow_node^th row, 6th column entry to zero
  prior_tree_table_temp(grow_node-1,6)=0;		// set grow_node^th row, 7th column entry to zero (why isn't left_daughter an int?)
  
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);		// let arma mat M equal prior_tree_table_temp
  M.insert_rows(left_daughter-1,2);		// Insert two rows at left_daughter^th row
  M(left_daughter-1,4)=-1;		// Let entry in left_daughter^th row, 5th column be -1. (denote terminal node?)
  M(left_daughter,4)=-1;      // Let entry in left_daughter+1^th row, 5th column be -1. (denote terminal node?)
  
  M(grow_node-1,0)=left_daughter;		// Let entry in grow_node^th row, first column of M be left_daughter (??)
  M(grow_node-1,1)=left_daughter+1;		// Let entry in grow_node^th row, second column of M be left_daughter+1 (??)
  NumericMatrix t=as<NumericMatrix>(wrap(M));		// now convert to NumericMatrix and call it t
  IntegerVector rname=seq_len(t.nrow());		// rname is a vector 1,2,3,..., up to the number of rows in the matrix
  
  List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  t.attr("dimnames") = dimnms;	// give t the attribute dimnames (list defined above)
  
  return(t);		// return the matrix t
}

//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix update_grow_obs_bcf(NumericMatrix prior_tree_matrix_temp,double grow_node,double left_daughter,double d,NumericVector ld_obs,NumericVector rd_obs){
  arma::mat prior_tree_matrix_temp2(prior_tree_matrix_temp.begin(),prior_tree_matrix_temp.nrow(),prior_tree_matrix_temp.ncol(),false);	// make an arma mat from the input matrix
  arma::vec ptm2=prior_tree_matrix_temp2.col(d);	// d+1^th column of the matrix. Why isn't d an integer?
  NumericVector ptm=wrap(ptm2);					// convert from arma mat to numeric vector
  ptm[ld_obs]=left_daughter;						// set elements of ptm indexed by ld_obs equal to left_daughter
  ptm[rd_obs]=left_daughter+1;					// set elements of ptm indexed by rd_obs equal to left_daughter+1
  prior_tree_matrix_temp(_,d)=ptm;				// set d+1^th column of input matrix equal to ptm
  
  return(prior_tree_matrix_temp);		// return the matrix
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix find_obs_to_update_grow_bcf(NumericMatrix prior_tree_matrix_temp,double left_daughter,double d,NumericVector ld_obs,NumericVector rd_obs){
  
  int rows=prior_tree_matrix_temp.nrow();		// number of rows in input matrix
  int cols=prior_tree_matrix_temp.ncol();		// number of columns in input matrix
  int elements=rows*cols;						// number of elements in input matrix
  std::vector<double> rows_obs(elements);		// create a vector of length elements (not initialized??, doesn't matter count++ within if statement, then resize_bcf)
  std::vector<double> cols_obs(elements);		// create a vectpr of length elements (not initialized??)
  int count=0;								// initialize integer count equal to 0.
  for(int i=0;i<prior_tree_matrix_temp.nrow();i++){			// loop of length equal to number of rows
    for(int j=0;j<prior_tree_matrix_temp.ncol();j++){		// loop of length equal to number of columns
      if(prior_tree_matrix_temp(i,j)>=left_daughter){		// if element in row i column j is greater than ot equal to left_daughter
        rows_obs[count]=i;								// set row_obs[count] equal to row number.
        cols_obs[count]=j;								// set col_obs[count] equal to column number. 
        count++;										// add 1 to the count. NOTE: COUNT++ IS WITHIN THE IF STATEMENT
      }
    }
  }
  rows_obs.resize(count);		// resize_bcf rows_obs so that it only contains its first count elements
  cols_obs.resize(count);		// resize_bcf cols_obs so that it only contains its first count elements
  
  if(rows_obs.size()!=0){			// if rows_obs is of nonzero length, i.e. a nonzero number of elements are greater than left_daughter
    for(int k=0;k< count;k++){		// loop of length count (equals length of rows_obs and cols_obs)
      if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])<left_daughter){			// if element given by k^th indices in rows_obs and cols_obs is less than left_daughter
        prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=prior_tree_matrix_temp(rows_obs[k],cols_obs[k]);	// sets equal to itself. But this if statement should never be true by construction of rows_obs and cols_obs?
      }else if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])==0){	
        prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=0;	// if already equals zero, then this line is unnecessary? Make empty if statement
      }else{   
        int temp=prior_tree_matrix_temp(rows_obs[k],cols_obs[k])+2;		// add two to relevant element of matrix. Why isn't this a numeric (part of numeric matrix)
        prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=temp;			// this and the above line could be one line?
      }
    }
  }
  //update prior_tree_matrix
  //insert extra column for the new split node 
  if(prior_tree_matrix_temp.ncol()<=d+1){								// If the number of columns is greater than d+1
    arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);		// M equals arma mat copy of the matrix
    M.insert_cols(prior_tree_matrix_temp.ncol(),1);  				// insert one column after last column
    NumericMatrix prior_tree_matrix=as<NumericMatrix>(wrap(M));		// convert back to NumericMatrix and save as prior_tree_matrix
  }
  
  arma::mat prior_tree_matrix_temp2(prior_tree_matrix_temp.begin(),prior_tree_matrix_temp.nrow(),prior_tree_matrix_temp.ncol(),false);	// copy as arma mat with new name
  arma::vec ptm2=prior_tree_matrix_temp2.col(d+1);		// create vector ptm2 equal to d+2^th column. But are there always this many columns? 
  NumericVector ptm=wrap(ptm2);	// convert to numeric vector. Name ptm.
  ptm[ld_obs]=left_daughter;		// Set all elements indexed by ld_obs equal to left_daughter
  ptm[rd_obs]=left_daughter+1;	// Set all elements indexed by rd_obs equal to left_daughter+1
  prior_tree_matrix_temp(_,d+1)=ptm;  // set the d+2^th column equal to ptm
  
  return(prior_tree_matrix_temp);	// return the matrix
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_daughter_obs_bcf(arma::mat& xmat,NumericVector obs_to_update,int split_var,double split_point){		// Reference to xmat? What is the purpose of "&"?
  
  arma::uvec ld_obs_ind(obs_to_update.size());		// create a vector of unsigned integers called ld_obd_ind (presumably indices). Not initialized?
  arma::uvec rd_obs_ind(obs_to_update.size());		// create a vector of unsigned integers called rd_obd_ind (presumably indices). Not initialized?
  
  arma::colvec sv_col;								// create a arma colvec called sv_col (not initialized?)
  List daughter_obs(2);								// create a List called daughter_obs with two elements. (List is a rcpp type)
  
  sv_col=xmat.col(split_var-1);						// let sv_col be equal to the split_var^th column of xmat
  
  arma::uvec obs_to_update_arma=as<arma::uvec>(obs_to_update);				// Let obs_to_update_arma be a vector of unsige integers equal to obs_to_update. (why isn't obs_to_update an IntegerVector?)
  ld_obs_ind = arma::find(sv_col.elem(obs_to_update_arma)<=split_point);		// ld_obs_ind is the indices of sv_col.elem(obs_to_update_arma)that are <= split_point. LHS is elements of sv_col given by indices in obs_to_update_arma.
  rd_obs_ind = arma::find(sv_col.elem(obs_to_update_arma)>split_point);		// rd_obs_ind is the indices of sv_col.elem(obs_to_update_arma)that are > split_point.
  
  NumericVector ld_ind2(as<NumericVector>(wrap(ld_obs_ind)));		// Let ld_ind2 (converted from arma uvec ld_obs_ind). Why not IntegerVector?
  NumericVector rd_ind2(as<NumericVector>(wrap(rd_obs_ind)));		// Let rd_ind2 (converted from arma uvec rd_obs_ind). Why not IntegerVector?
  
  NumericVector ld_obs=obs_to_update[ld_ind2];		// ld_obs is the vector of elements of obs_to_update given by the indices in ld_ind2 [obs_to_update might itself be a set of indices?]
  NumericVector rd_obs=obs_to_update[rd_ind2];		// rd_obs is the vector of elements of obs_to_update given by the indices in rd_ind2
  
  daughter_obs[0]=ld_obs;			// Fist element of list is the vector ld_obs
  daughter_obs[1]=rd_obs;			// Second element of list is the vector rd_obs
  
  return(daughter_obs);		// Return the list
  
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector find_term_cols_bcf(NumericMatrix tree_matrix_temp,int terminal_node){
  
  arma::mat tree_matrix_temp2(tree_matrix_temp.begin(),tree_matrix_temp.nrow(),tree_matrix_temp.ncol(),false);	// convert input matrix to arma mat called tree_matrix_temp2
  int count=0;												// Initialize count at 0
  std::vector<double> term_cols(tree_matrix_temp.ncol());		// term_cols is a vector of length equal to the number of columns in input matrix. Why not IntegerVector>
  
  for(int j=0;j<tree_matrix_temp.ncol();j++){					// loop of length equal to the number of columns in the input matrix
    
    arma::vec tempcol=tree_matrix_temp2.col(j);				// arma vec tempcol equals j^th column of matrix
    arma::uvec term_nodes=find(tempcol==terminal_node);		// unsigned integer vec term_nodes equals indices of elements of tempcol that are equal to terminal_node (input value).
    
    if(term_nodes.size()>0){		// If term_nodes has nonzero length, i.e. some elements of tempcol are equal to terminal_node
      term_cols[count]=j;		// Then let the count+1^th (starting from zero) element of term_cols be j (the number minus 1 of the column in the loop) 
      count++;		// only increase count when the condition is true. term_cols therefore contains indexes (starting at 0) all the columns with at least some elements equal to terminal_nodes
    }
  }
  
  term_cols.resize(count);	// reduce size to nonempty entries. Note that only all of term_cols filled in if all columns have some elements equal to terminal_nodes.
  
  return(wrap(term_cols));		// convert to R object (NumericVector) and return.
  
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector get_grow_obs_bcf(arma::mat& xmat,NumericVector grow_obs,int split_var){
  
  arma::vec sv_col=xmat.col(split_var-1);				// arma vec sv_col is split_var^th element of input matrix
  arma::uvec grow_obs2(as<arma::uvec>(grow_obs));		// convert input vector to an unsigned integer vector called grow_obs2. Why isn't the input vector an IntegerVector?
  arma::vec get_min=sv_col.elem(grow_obs2);			// Obtain the elements of the split_var^th column of the input matrix indexed by grow_obs
  
  return(wrap(get_min));		// convert the vector to an R object and return
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector get_grow_obs_in_z_bcf(arma::vec& z_ar,NumericVector grow_obs){
  
  arma::uvec grow_obs2(as<arma::uvec>(grow_obs));		// convert input vector to an unsigned integer vector called grow_obs2. Why isn't the input vector an IntegerVector?
  arma::vec get_min=z_ar.elem(grow_obs2);			// Obtain the elements of the split_var^th column of the input matrix indexed by grow_obs
  
  return(wrap(get_min));		// convert the vector to an R object and return
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List grow_tree_bcf(arma::mat& xmat,NumericVector y,NumericMatrix prior_tree_matrix,int grow_node,NumericMatrix prior_tree_table,int splitvar,double splitpoint,
               NumericVector terminal_nodes,NumericVector grow_obs,double d,NumericVector get_min,arma::mat& data_curr_node)
{
  
  NumericMatrix prior_tree_matrix_temp=clone(prior_tree_matrix);		// create a copy
  NumericMatrix prior_tree_table_temp=clone(prior_tree_table);		// create a copy
  double yy=xmat.n_cols;												// set yy equal to the number of columns of xmat
  IntegerVector xx=seq_len(yy);										// create an integer vector of {1,2,3,...,yy}
  prior_tree_table_temp(grow_node-1,3)=splitpoint;					// For matrix prior_tree_table_temp, set the grow_node^th row, fourth column entry to splitpoint
  prior_tree_table_temp(grow_node-1,2)=splitvar;						// For matrix prior_tree_table_temp, set the grow_node^th row, third column entry to splitvar 
  prior_tree_table_temp(grow_node-1,4)=1;								// For matrix prior_tree_table_temp, set the grow_node^th row, fifth column entry to 1
  //get data subset for left and right daughter nodes
  List daughter_obs=get_daughter_obs_bcf(xmat,grow_obs,splitvar,splitpoint);		// apply function from line 397 to get the list of two vectors ld_obs and rd_obs
  NumericVector ld_obs=daughter_obs[0];										// first element of list
  NumericVector rd_obs=daughter_obs[1];										// second element of list
  
  if(prior_tree_table_temp.nrow()==grow_node){											// If the number of rows of prior_tree_table_temp equals grow_node
    prior_tree_table_temp=add_rows_bcf(prior_tree_table_temp,grow_node);					// add two rows using the function defined on line 22. (add to tree table)
    prior_tree_matrix_temp=addcol_bcf(prior_tree_matrix_temp,grow_node,ld_obs,rd_obs);		// add a column using function defined on line 47. (add to tree matrix)
  }else{
    //if grow node is in the middle of the tree
    NumericVector nodes_d;															// create a numeric vector called nodes_d
    nodes_d=prior_tree_matrix_temp(_,d);											// set nodes_d equal to the d+1^th column of prior_tree_matrix_temp (d is an input double. Perhaps d should be an integer)
    NumericVector nodes_at_depth=sort_unique(nodes_d);								// obtains unique elemnts of nodes_d, and then sort in ascending order
    NumericVector nodes_at_depth1=remove_zero_bcf(nodes_at_depth);						// then remove zeros and save as nodes_at_depth1
    NumericVector gn_pos=get_gnp_bcf(nodes_at_depth1, grow_node);						// get_gnp_bcf defined on line 157. gn_pos is vector of index numbers of elements of nodes_at_depth1 that are equal to grow_node
    arma::uvec prev_uvec= arma::find(as<arma::vec>(nodes_at_depth1)<grow_node);		// vector of indices of elements of nodes_at_depth1 that are less than grow_node
    arma::vec nontermvec=Rcpp::as<arma::vec>(nodes_at_depth1);						// Convert nodes_at_depth1 to an arma vec called nontermvec
    nontermvec=nontermvec.elem(prev_uvec);								// redefine nontermvec such that it only contains its elements indexed by prev_uvec
    NumericVector prev= as<NumericVector>(wrap(nontermvec));			// convert to NumericVector and call prev
    double prev_nonterm=0;												// create a double called prev_nonterm. Initialize to zero
    if(prev.size()!=0){														// If prev has nonzero size
      if(prior_tree_table.ncol()<5) throw std::range_error("Line 531");
      arma::uvec find_nonterm=find_internal_nodes_bcf(prior_tree_table);		// Use function defined on line 241. Returns index vector of internal nodes. Note not applied to temporary version of table
      //should only find internal nodes at the current depth
      prev_nonterm=find_prev_nonterm_bcf(find_nonterm,prev);					// Function defined on line 256. Total number of elements of find_nonterm equal to elements of prev.
    }
    double left_daughter=grow_node +2*(prev_nonterm)+(nodes_at_depth1.size()-gn_pos[0]);	// define left_daughter
    NumericVector ptt=prior_tree_table(_,1);												// ptt is equal to second column of prior_tree_table. (Note: not temp version of prior_tree_table)
    arma::uvec node_to_update=find_nodes_to_update_bcf(as<arma::uvec>(ptt),left_daughter);		// node_to_update is vector of indices of all elements of ptt that are >= left_daughter. Function defined on line 273.
    //increase the node number of nodes after the grow node by two (because 2 daughter nodes are now appended to grow node)
    //do this for all observations except those that already belong to a terminal node (a value of 0)
    if(node_to_update.size()==0){																				// If node_to_update is of length zero
      if(prior_tree_matrix_temp.ncol()>d+1){																	// If prior_tree_matrix_temp has more than d+1 columns
        prior_tree_table_temp=set_daughter_to_end_tree_bcf(grow_node,prior_tree_table_temp,left_daughter);		// Use function defined on line 79 to add rows and change some values of prior_tree_table_temp
        prior_tree_matrix_temp=update_grow_obs_bcf(prior_tree_matrix_temp,grow_node,left_daughter,d+1,ld_obs,rd_obs);	// Set elements of d+2^th column indexed by ld_obs equal to left_daughter and elements indexed by rd_obs equal to left_daughter+1
      }else{
        //if the daughter node number already exists in the tree and existing node numbers have to be updated
        //daughter nodes need to be added to the end of the table not in the center of it
        prior_tree_table_temp=set_daughter_to_end_tree_bcf(grow_node,prior_tree_table_temp,left_daughter);		// Use function defined on line 79 to add rows and change some values of prior_tree_table_temp
        prior_tree_matrix_temp=set_daughter_to_end_mat_bcf(d,prior_tree_matrix_temp,left_daughter,ld_obs,rd_obs);	// See function defined on line 107. In last column set the elements indexed by ld_obs to left_daughter and those indexe by rd_obs to left_daughter+1. Should left_daughter be an integer?
      }
    }else{
      //if the daughter node number already exists in the tree and existing node numbers have to be updated
      prior_tree_table_temp=set_tree_to_middle_bcf(wrap(node_to_update),prior_tree_table_temp,grow_node,left_daughter);	// line 283 defines the function. Adds two tows and changes first and second column.
      prior_tree_matrix_temp=find_obs_to_update_grow_bcf(prior_tree_matrix_temp,left_daughter,d,ld_obs,rd_obs);		// line 330 defines the function.
    }
  }
  
  List ret(2);					// list containing two matrices
  ret[0]=prior_tree_table_temp;	// new tree table
  ret[1]=prior_tree_matrix_temp;	// new tree matrix
  
  return(ret);					// return the list
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericMatrix set_daughter_bcf(int left_daughter,int right_daughter,IntegerVector ld_obs,IntegerVector rd_obs,NumericMatrix tree_matrix_temp,double term_cols){
  
  arma::mat tree_matrix_temp2(tree_matrix_temp.begin(),tree_matrix_temp.nrow(),tree_matrix_temp.ncol(),false);	// red the input matrix as an arma mat 
  arma::vec arma_col=tree_matrix_temp2.col(term_cols+1);		// arma_col is term_cols+2^th column of tree_matrix_temp2
  NumericVector col(as<NumericVector>(wrap(arma_col)));		// col is arma_col converted to a Numeric_Vector
  col[ld_obs]=left_daughter;									// set the elements of col indexed by ld_obs equal to left_daughter
  col[rd_obs]=right_daughter;									// set the elements of col indexed by rd_obs equal to right_daughter
  tree_matrix_temp(_,term_cols+1)=col;						// Set the term_cols+2^th column of tree_matrix_temp equal to col
  
  return(tree_matrix_temp);									// return the matrix
}

//######################################################################################################################//

// [[Rcpp::export]]

IntegerVector order__bcf(NumericVector x) {	// gives vector of position of largest value, then position of second largest value, and so on.
  NumericVector sorted = clone(x).sort();		// sorted is x in ascending order
  std::reverse(sorted.begin(), sorted.end()); // reverse so that it is in descending order. Could use one line with std::sort(clone(x).begin(), clone(x).end(), std::greater<>())
  
  return match(sorted, x);	//  match is the Rcpp sugar version of the R function match, which returns a vector of the positions of the first matches of the first argument in the second.
}
//######################################################################################################################//

// [[Rcpp::export]]

IntegerVector orderforOW__bcf(NumericVector x) {	// gives vector of position of smallest value, then position of second smallest value, and so on.
  NumericVector sorted = clone(x).sort();		// sorted is x in ascending order
  //std::reverse(sorted.begin(), sorted.end()); // reverse so that it is in descending order. Could use one line with std::sort(clone(x).begin(), clone(x).end(), std::greater<>())
  
  return match(sorted, x);	//  match is the Rcpp sugar version of the R function match, which returns a vector of the positions of the first matches of the first argument in the second.
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double get_tree_prior_bcf(NumericMatrix tree_table,NumericMatrix tree_matrix,double alpha,double beta){
  double propsplit=1;
  IntegerVector d;
  IntegerVector d2;
  int col=tree_matrix.ncol();
  std::vector<int> int_nodes_index(100*col); // Why  100* ? 
  std::vector<int> term_nodes_index(100*col); //
  
  int index_count=0;  
  int index_count2=0;
  arma::uvec internal_nodes_prop=find_internal_nodes_bcf(tree_table);
  NumericVector terminal_nodes_prop_wrapped=find_term_nodes_bcf(tree_table);
  arma::uvec terminal_nodes_prop=as<arma::uvec>(terminal_nodes_prop_wrapped);
  arma::mat tree_matrix2(tree_matrix.begin(),tree_matrix.nrow(),tree_matrix.ncol(),false);
  int count=internal_nodes_prop.size();
  int count_term_nodes=terminal_nodes_prop.size();
  
  if(count==0) propsplit=1-alpha;
  
  for(int k=0;k<count;k++){ 
    for(int j=0;j<tree_matrix.ncol();j++){
      arma::vec armacol=tree_matrix2.col(j);
      arma::uvec found=find(armacol==internal_nodes_prop[k]);      
      if(found.size()>0){        
        int_nodes_index[index_count]=j;
        index_count++;
        break;
      }        
    }
    int_nodes_index.resize(index_count);
    if(int_nodes_index.size()!=0){      
      d=unique(as<IntegerVector>(wrap(int_nodes_index)));
      double d1=d[0];
      propsplit*=alpha*pow((d1+1),-beta) ;  
    }
    std::vector<int> temp(col);
    int_nodes_index=temp;
    index_count=0;
  } 
  
  if(count!=0){  
    for(int k=0;k<count_term_nodes;k++){ 
      for(int j=0;j<tree_matrix.ncol();j++){
        arma::vec armacol=tree_matrix2.col(j);
        arma::uvec found=find(armacol==terminal_nodes_prop[k]);      
        if(found.size()>0){        
          term_nodes_index[index_count2]=j;
          index_count2++;
          break;
        }        
      }
      term_nodes_index.resize(index_count2);
      if(term_nodes_index.size()!=0){      
        d2=unique(as<IntegerVector>(wrap(term_nodes_index)));
        double d1=d2[0];
        propsplit*=1-alpha*pow((d1+1),-beta) ;  
      }
      std::vector<int> temp2(col);
      term_nodes_index=temp2;
      index_count2=0;
    }
  }
  
  return(propsplit);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericMatrix start_tree_bcf(double start_mean,double start_sd){
  
  NumericMatrix treemat(1,7);											// crete matrix treemat with 1 row, 7 columns (not initialized).
  double rand=R::rnorm(start_mean,start_sd);							// rand is a draw from a standard normal distribution. An alternative to rnorm would be randn() from armadillo.
  NumericVector testrow = NumericVector::create(0,0,0,0,-1,rand,0);	// define the vector testrow
  for(int k=0;k<1;k++){												// loop of length 1? What is the purpose of this loop? Why not use 0 instead of k in the inner loop?
    for(int j=0;j<7;j++){											// loop of length 7 (number of columns)
      treemat(k,j)=testrow[j];									// set values of treemat (just one row). But why not do this directly instead of first defining testrow?
    }
  }
  List dimnms = // two vec. with static names
    List::create(CharacterVector::create("1"),
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  treemat.attr("dimnames") = dimnms;		// give attributes to the R object treemat
  
  return(treemat);						// return the matrix
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericMatrix start_matrix_bcf(int n){
  NumericMatrix mat(n,1);					// Matrix with n rows, 1 column
  std::fill(mat.begin(), mat.end(), 1);	// set all matrix elements to 1
  return(mat);							// return the matrix
}

//######################################################################################################################//

// [[Rcpp::export]]
List evaluate_model_occams_window_bcf(NumericVector tree_lik,double lowest_BIC,double c,List tree_list,List tree_mat_list,IntegerVector tree_parent){
  IntegerVector sorted_lik_index=order__bcf(tree_lik);		// Function order__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on.
  std::vector<double> to_be_removed(tree_lik.size());		// create a vector of length equal to that of tree_lik (uninitialized). Why isn't this an integer vector?
  int s=0;												// set s equal to 0
  
  // check tree is in Occam's window if it isn't then delete it from list                    
  while((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC) > c){	// while element of tree_lik -lowestBIC >c. (starts with largest, and then s increases, giving next largest element of tree_lik, and so on
    if(s==(tree_lik.size()-1)){								// finish if at last element of tree_lik. Should this be tree_lik.size() instead of tree_lik.size()-1? Want to be able to remove the last tree?
      break;												// break out of while loop?
    }
    //delete tree from tree list
    if((tree_lik[sorted_lik_index[s]-1])-(lowest_BIC)>c){	// same condition as for while loop above. Presuably tree_lik actually includes the BICs
      //set indicator for index of trees to be removed
      to_be_removed[s]=sorted_lik_index[s]-1;				// include the index number sorted_lik_index[s]-1 in to_be_removed. (Records models to be removed).
      s+=1;												// increase s by 1
    }
  }
  
  to_be_removed.resize(s);											// resize_bcf to_be_removed to exlude elements that aren't filled in
  IntegerVector remove_order_index(to_be_removed.size());				// remove_order_index is of same length as to_be_removed
  //delete elements from the higest index down 
  remove_order_index=order__bcf(wrap(to_be_removed));						// Gives vector of position of largest element, then position of second largest argument, and so on. Also convert to R object.
  
  for(int j=0;j<s;j++){												// loop of size s. either s<=(tree_lik.size()-1 (only equal if all trees in occam's window)
    tree_list.erase(to_be_removed[remove_order_index[j]-1]);		// removes element with index given by to_be_removed[remove_order_index[j]-1]. Note -1 because order__bcf gives values from 1 to s.
    tree_mat_list.erase(to_be_removed[remove_order_index[j]-1]);	// removes element with index given by to_be_removed[remove_order_index[j]-1]. Note -1 because order__bcf gives values from 1 to s.
    tree_lik.erase(to_be_removed[remove_order_index[j]-1]);			// removes element with index given by to_be_removed[remove_order_index[j]-1]. Note -1 because order__bcf gives values from 1 to s.
    tree_parent.erase(to_be_removed[remove_order_index[j]-1]);		// removes element with index given by to_be_removed[remove_order_index[j]-1]. Note -1 because order__bcf gives values from 1 to s.
  }
  List ret(4);				// create a list of length 4
  ret(0)=tree_lik;			// first element of list
  ret(1)=tree_list;			// second element of list
  ret(2)=tree_mat_list;		// third element of list
  ret(3)=tree_parent;			// fourth element of list
  
  return(ret);				// return the list
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector get_testdata_term_obs_bcf(NumericMatrix test_data,NumericMatrix tree_data,NumericVector term_node_means) {
  //Function to make predictions from test data, given a single tree and the terminal node predictions, this function will be called
  //for each tree accepted in Occam's Window.
  arma::mat arma_tree(tree_data.begin(), tree_data.nrow(), tree_data.ncol(), false);		// read matrix tree data, call it arma_tree
  arma::mat testd(test_data.begin(), test_data.nrow(), test_data.ncol(), false);			// read matrix test_data, call it testd
  NumericVector terminal_nodes=find_term_nodes_bcf(tree_data);								// find term nodes function defined line 168. Gives index of values of tree_data that are term nodes (indices from 1 to length of vector). Why not integer vector?
  arma::vec arma_terminal_nodes=Rcpp::as<arma::vec>(terminal_nodes);						// convert terminal_nodes to arma_vec. Could use uvec
  NumericVector tree_predictions;															// this vector is not used anywhere. delete line?
  //for each internal node find the observations that belong to the terminal nodes
  NumericVector predictions(test_data.nrow());											// Create a vector of length equal to the number of rows of test_data. Not initialized.
  if(terminal_nodes.size()==1){
    double nodemean=tree_data(terminal_nodes[0]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    predictions=rep(nodemean,test_data.nrow());
  }
else{
  for(int i=0;i<terminal_nodes.size();i++){						// loop of same length as terminal_nodes
    arma::mat subdata=testd;									// arma mat equal to testd (test data)
    int curr_term=terminal_nodes[i];							// curr_term is i^th element of terminal_nodes
    int row_index;												// create row_index, not initialized
    int term_node=terminal_nodes[i];							// term_node is i^th element of terminal_nodes
    if(curr_term % 2==0){										// If curr_term is even (even index value of current terminal node)
      row_index=terminal_nodes[i];								// row_ndex is i^th element of terminal_nodes 
    }else{
      row_index=terminal_nodes[i]-1;							// Otherwise row_ndex is i^th element of terminal_nodes MINUS ONE
    }
    //save the left and right node data into arma uvec
    arma::vec left_nodes=arma_tree.col(0);				// left_nodes is first column of arma_tree (tree_data)
    arma::vec right_nodes=arma_tree.col(1);				// right_nodes is second column of arma_tree (tree_data)
    arma::mat node_split_mat;							// create matrix node_split_mat (uninitialized)
    node_split_mat.set_size(0,3);						// node_split_mat has 0 rows, 3 columns?
    while(row_index!=1){									// While row_index is not equal to 1. 
      //for each terminal node work backwards and see if the parent node was a left or right node
      //append split info to a matrix 
      int rd=0;													// integer variable rd initialized as 0
      arma::uvec parent_node=arma::find(left_nodes == term_node);	// parent_node is vector of indices of elements of left_nodes (first col of tree_data) equal to i^th element of termnal_nodes
      if(parent_node.size()==0){							// If no elements of left_nodes are equal to i^th element of termina_nodes
        parent_node=arma::find(right_nodes == term_node);	// Indices of elements of right_nodes (second column of tree_data) that are equal to i^th element of terminal_nodes
        rd=1;												// change value of rd to 1.
      }

      node_split_mat.insert_rows(0,1);					// insert the first row of node_split_mat. Elements are 0 by deffault. (0 is the row number, 1 is the number of rows to insert)
      node_split_mat(0,0)=tree_data(parent_node[0],2);	// 1st row 1st column entry is first index is from tree_data(row= first element of parent_node, 3rd column)
      node_split_mat(0,1)=tree_data(parent_node[0],3);	// 1st row 2nd column entry is first index is from tree_data(row= first element of parent_node, 4th column)
      node_split_mat(0,2)=rd;								// 1st row 3rd column entry is rd=1.
      row_index=parent_node[0] +1;						// row_index is reset to first element of parent_node +1.
      term_node=parent_node[0]+1;							// term_node is reset to first element of parent_node +1.
    }
    //fill in the predicted value for tree

    arma::uvec pred_indices;								// create vector of unsigned integers. Not initialized.
    int split= node_split_mat(0,0)-1;						// first row, first column of node_split_mat tree_data(parent_node[0],2) for last time row_index is nonzero in while loop above. Assuming row_index can't start equal to 1, because then node_split_mat would have no rows.
    arma::vec tempvec = testd.col(split);					// column of test_data indexed by split
    double temp_split = node_split_mat(0,1);				// first row, second column of node_split_mat tree_data(parent_node[0],3) for last time row_index is nonzero in while loop above. Assuming row_index can't start equal to 1, because then node_split_mat would have no rows.
    if(node_split_mat(0,2)==0){								// If rd=0 (and row_index didn't start at 1), i.e. if at least one element of left_nodes is equal to i^th element of terminal_nodes.
      pred_indices = arma::find(tempvec <= temp_split);		// Indices for elements of testd.col(split) that are <= node_split_mat(0,1)
    }else{													// If rd!=0, i.e. if no element of left_nodes is equal to i^th element of terminal_nodes.
      pred_indices = arma::find(tempvec > temp_split);	// Indices for elements of testd.col(split) that are > node_split_mat(0,1) 
    }
    arma::uvec temp_pred_indices;						// vector of unsigned integers
    arma::vec data_subset = testd.col(split);			// column of test_data indexed by split. Why not use tempvec? Equal and already defined.
    data_subset=data_subset.elem(pred_indices);			// let data_subset be those elements of testd.col(split) <= or > temp_split depeding on whether rd=0.
    int n=node_split_mat.n_rows;						// n is number of rows of node_split_mat
    for(int j=1;j<n;j++){								// loop of length equal to number of rows of node_split_mat MINUS 1. Note no j=0.
      int curr_sv=node_split_mat(j,0);				// curr_sv is element in j+1^th row, 1st column of node_split_mat
      double split_p = node_split_mat(j,1);			// split_p is element in j+1^th row, 2nd column of node_split_mat
      data_subset = testd.col(curr_sv-1);				// data_subset is curr_sv^th column of test data
      data_subset=data_subset.elem(pred_indices);		// data_subset is elements of curr_sv^th column of test data indexed by pred_indices
      if(node_split_mat(j,2)==0){									// If j+1^th column, 3rd row of node_split_mat (rd) is zero. i.e. split is to the left.
        //split is to the left
        temp_pred_indices=arma::find(data_subset <= split_p);	// indices of data_subset elements <= split_p [if split is to left, find test data that goes left?]
      }else{
        //split is to the right
        temp_pred_indices=arma::find(data_subset > split_p);	// indices of data_subset elements > split_p [if split is to right, find test data that goes right?]
      }
      pred_indices=pred_indices.elem(temp_pred_indices);			// keep elements of pred_indices indexed by temp_pred_indices (presumably keep data falling to correct side of split_p and pred_indices)
      if(pred_indices.size()==0){									// If pred_indices is now empty
        continue;													// Then skip to the next iteration of the for-loop.
      }
    }
    double nodemean=tree_data(terminal_nodes[i]-1,5);				// let nodemean equal tree_data row terminal_nodes[i]^th row , 6th column. The minus 1 is because terminal nodes consists of indices starting at 1, but need indices to start at 0.
    IntegerVector predind=as<IntegerVector>(wrap(pred_indices));	// convert pred_indices from a uvec to the IntegerVector predind
    predictions[predind]= nodemean;			// Let the elements of predictions indexed by predind be equal to nodemean

    
  } 
}
  return(predictions);						// return the vector of predictions
}
//######################################################################################################################//

// [[Rcpp::export]]
List resize_bcf(const List& x, int n ){
  List y(n) ;								// create list with n elements
  for( int i=0; i<n; i++) y[i] = x[i] ;	// set all n elements to first n elements of x
  
  return y ;								// return the list
}
//######################################################################################################################//

// [[Rcpp::export]]
List resize_bigger_bcf( const List& x, int n ){
  int oldsize = x.size() ;					// let oldsize be length of x
  List y(n) ;									// create list with n elements
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;	// set first oldsize elements of y equal to those of x. (this allows y to be longer than x)
  return y ;									// return the list
}
//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat J_bcf(NumericMatrix treetable_temp,NumericMatrix obs_to_nodes_temp,NumericVector tree_term_nodes){
  //this function will make a binary nxb matrix where each column assigns observations to terminal nodes
  
  //create J_bcf matrix with correct dimensions and fill with zeros
  arma::mat Jmat(obs_to_nodes_temp.nrow(), tree_term_nodes.size());		// create matrix with number of rows equal to that of obs_to_nodes_temp, and number of columns equal to length of tree_term_nodes
  Jmat.zeros();
  
  //for each terminal node get the observations associated with it and set column
  for(int i=0;i<tree_term_nodes.size();i++){						// loop of same length as tree_term_nodes
    double tn=tree_term_nodes[i];									// let tn equal the i-1^th value in tree_term_nodes
    arma::uvec term_obs=find_term_obs_bcf(obs_to_nodes_temp,tn);		// function find_term_obs_bcf. Gives indices of elements equal to tn (for leftmost column of obs_to_nodes_temp that has elements equal to tn).
    //assign term_obs to the correct index of J_bcf
    NumericVector term_obs2=as<NumericVector>(wrap(term_obs));		// convert term_obs to Numeric Vector. Could use IntegerVector?
    NumericVector obs_col(obs_to_nodes_temp.nrow());				// create a vector obs_col of length equal to number of rows of obs_to_nodes_temp. Other values automatically 0.
    obs_col[term_obs2]=1;											// elements of obs_col indexed by term_obs2 set equal to 1. Other values automatically 0.
    arma::vec colmat=Rcpp::as<arma::vec>(obs_col);					// convert obs_col to vec called colmat 
    Jmat.col(i)= colmat;											// Set i-1^th column of Jmat equal to comat
  }
  return(Jmat);			// return the arma mat Jmat
}

//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector mu_vector_bcf(List sum_treetable,int n){
  NumericVector mu_vec;										// create a NumericVector called mu_vec. Not initialized
  
  for(int j=0;j<sum_treetable.size();j++){						// loop of length equal to that of sum_treetable list
    NumericMatrix curr_tree=sum_treetable[j];					// create a list equal to j-1^th element of sum_treetable
    NumericVector tree_term_nodes=find_term_nodes_bcf(curr_tree);	// WHAT IS THE PURPOSE OF THIS LINE? find_term_nodes_bcf on line . Gives indices of curr_tree elements that correspond to term_nodes.
    NumericVector term_means1=curr_tree(_,5);					// term_means1 equals 6th column of curr_tree.
    NumericVector term_means=remove_zero_bcf(term_means1);			// Function defined on line 137. Removes zero from vector.
    
    for(int i=0;i<term_means.size();i++){						// loop of length equal to that of term_means (after zeros removed)
      mu_vec.push_back(term_means[i]);							// Append term_means[i] to mu_vec. End up with vectors for each j stacked into one mu vector (j indexes the trees, and i appears to index terminal nodes within trees).
    }
  }
  return(mu_vec);												// return the vector
}
//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat W_bcf(List sum_treetable ,List sum_obs_to_nodes,int n){
  //this will take in a list of obs to node matrices for each tree in the sum make the J_bcf matrix assigning observations to terminal nodes
  //J_bcf is an nxb_j matrix. It will then iteratively append the J_bcf matrix to itself to get the overall W_bcf matrix which has dimensions nxsumb_j
  //Rcout << "Size of tree table list = " << sum_treetable.size() <<" \n";
  //Rcout << "Size of tree mat list = " << sum_obs_to_nodes.size() <<" \n";
  //create empty matrix to which we will append individual J_bcf matrices
  arma::mat W_bcf(n,0);															// create a matrix with n(number of obs) rows and zero columns
  int upsilon=0;															// initialize the variable upsilon equal to zero
  for(int j=0;j<sum_obs_to_nodes.size();j++){								// loop of length equal to that of sum_obs_to_nodes (presumably equals the number of trees in the sum)
    
    NumericMatrix curr_tree=sum_treetable[j];								// curr tree is j-1^th element of the list sum_treetable
    //Rcout << "After Line 883.\n";
    NumericMatrix curr_obs_nodes=sum_obs_to_nodes[j];						// curr_pbs_nodes is j-1^th element of sum_obs_to_nodes
    NumericVector tree_term_nodes=find_term_nodes_bcf(curr_tree);				// tree_term_nodes gives indices of curr_tree elements that correspond to term_nodes.
    
    int b_j=tree_term_nodes.size();											// b_j is number of terminal nodes in tree j
    //will make J_bcf as we go in BART-BMA no need to create it again here....
    arma::mat Jmat=J_bcf(curr_tree,curr_obs_nodes,tree_term_nodes);				// create the J_bcf matrix where each column assigns observatins to terminal nodes. Function defined on line 786
    //Rcout << "Get to 890.\n";
    //Rcout << "Number of rows of J_mat "<< Jmat.n_rows <<".\n";
    //Rcout << "Number of rows of curr_obs_nodes "<< curr_obs_nodes.nrow() <<".\n";
    //Rcout << "Number of rows of curr_tree "<< curr_tree.nrow() <<".\n";
    
    //Rcout << "Number of cols of W_bcf "<< W_bcf.n_rows <<".\n";
    
    W_bcf.insert_cols(upsilon,Jmat);											// Add matrix J_bcf to right of W_bcf matrix
    //Rcout << "Get to 893.\n";
    upsilon+=b_j;															// update point to add matrix to indez of last column +1
    //Rcout << "Get to 900.\n";
    
  }
  
  return(W_bcf);				// return the matrix
  // ret[1]=mu_vec;
}
//###########################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat W_tauround1_bcf(NumericMatrix treetable , NumericMatrix obs_to_nodes,int n){
  //this will take in a list of obs to node matrices for each tree in the sum make the J_bcf matrix assigning observations to terminal nodes
  //J_bcf is an nxb_j matrix. It will then iteratively append the J_bcf matrix to itself to get the overall W_bcf matrix which has dimensions nxsumb_j
  //Rcout << "Size of tree table list = " << sum_treetable.size() <<" \n";
  //Rcout << "Size of tree mat list = " << sum_obs_to_nodes.size() <<" \n";
  //create empty matrix to which we will append individual J_bcf matrices
  arma::mat W_bcf(n,0);															// create a matrix with n(number of obs) rows and zero columns
  int upsilon=0;															// initialize the variable upsilon equal to zero

    NumericMatrix curr_tree=treetable;								// curr tree is j-1^th element of the list sum_treetable
    //Rcout << "After Line 883.\n";
    NumericMatrix curr_obs_nodes=obs_to_nodes;						// curr_pbs_nodes is j-1^th element of sum_obs_to_nodes
    NumericVector tree_term_nodes=find_term_nodes_bcf(curr_tree);				// tree_term_nodes gives indices of curr_tree elements that correspond to term_nodes.
    
    int b_j=tree_term_nodes.size();											// b_j is number of terminal nodes in tree j
    //will make J_bcf as we go in BART-BMA no need to create it again here....
    
    //Rcout << "number of rows of curr_obs_nodes = " << curr_obs_nodes.nrow() << ".\n";
    arma::mat Jmat=J_bcf(curr_tree,curr_obs_nodes,tree_term_nodes);				// create the J_bcf matrix where each column assigns observatins to terminal nodes. Function defined on line 786
    //Rcout << "Get to 890.\n";
    //Rcout << "Number of rows of J_mat "<< Jmat.n_rows <<".\n";
    //Rcout << "Number of rows of curr_obs_nodes "<< curr_obs_nodes.nrow() <<".\n";
    //Rcout << "Number of rows of curr_tree "<< curr_tree.nrow() <<".\n";
    
    //Rcout << "Number of rows of W_bcf "<< W_bcf.n_rows <<".\n";
    //Rcout << "Number of cols of W_bcf "<< W_bcf.n_cols <<".\n";
    
    W_bcf.insert_cols(upsilon,Jmat);											// Add matrix J_bcf to right of W_bcf matrix
    //Rcout << "Get to 965.\n";
    upsilon+=b_j;															// update point to add matrix to indez of last column +1
    //Rcout << "Get to 900.\n";
    
  
  
  return(W_bcf);				// return the matrix
  // ret[1]=mu_vec;
}
//############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sumtree_likelihood_function_bcf_bcf(NumericVector y_temp,List sum_treetable_mu ,List sum_treetable_tau ,
                                       List sum_obs_to_nodes_mu,List sum_obs_to_nodes_tau,
                                       int n,double a_mu,double a_tau,double nu,double lambda, NumericVector z){
  
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec	
  //make W_bcf and mu matrices for the sum of trees
  
  arma::mat Wmat_mu=W_bcf(sum_treetable_mu,sum_obs_to_nodes_mu,n);	// Create W_bcf matrix. Function defined on line 829
  //Rcout << "Wmat_mu works.\n";
  arma::mat Wmat_tau=W_bcf(sum_treetable_tau,sum_obs_to_nodes_tau,n);	// Create W_bcf matrix. Function defined on line 829
  //Rcout << "Wmat_tau works.\n";
  
  double b_mu=Wmat_mu.n_cols;
  double b_tau=Wmat_tau.n_cols;
  Wmat_tau.each_col()%=z_ar;
  arma::mat Wmat = join_rows(Wmat_mu,Wmat_tau);
  
  double b=Wmat.n_cols;									// b is number of columns of W_bcf matrix (omega in the paper)
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);			// convert input y_temp to arma vec called yvec
  arma::mat y(n,1);										// create a matrix with n (number of observations) rows and 1 column
  y.col(0)=yvec;										// set first column of y equal to yvec
  //get exponent
  double expon=(n+nu+b)/2;								// set the expoenent (equation 5 in the paper)
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;								// yty = y transpose y (sum of squares)
  
  //get t(y)inv(psi)J_bcf
  arma::mat ytW=y.t()*Wmat;								// y transpose W_bcf
  //get t(J_bcf)inv(psi)J_bcf  
  arma::mat WtW=Wmat.t()*Wmat;							// W_bcf transpose W_bcf
  //get jpsij +aI
  arma::mat aI(b,b);									// create b by b matrix called aI. NOT INIIALIZED. 
  aI=aI.eye();										// a times b by b identity matrix. The .eye() turns aI into an identity matrix.
  arma::vec a_vec_mu = a_mu*arma::ones<arma::vec>(b_mu);
  arma::vec a_vec_tau = a_tau*arma::ones<arma::vec>(b_tau);
  arma::vec a_vec(b);
  a_vec.head(b_mu) = a_vec_mu;
  a_vec.tail(b_tau) = a_vec_tau;
  aI.diag() = a_vec;
  
  arma::mat sec_term=WtW+aI;							//
  arma::mat sec_term_inv=sec_term.i();					// matrix inverse expression in middle of eq 5 in the paper. The .i() obtains the matrix inverse.
  //get t(J_bcf)inv(psi)y
  arma::mat third_term=Wmat.t()*y;						// W_bcf transpose Y
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;			// matrix expression in middle of equation 5
  arma::mat rel=-expon*log(nu*lambda - mvm +yty);		// log of all of equation 5 (i.e. the log of the marginal likelihood of the sum of tree model)
  double rel2=as<double>(wrap(rel));					// convert to double (from 1 by 1 arma mat)
  return(rel2);											// return the log marginal likelihood
}    
//############################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sumtree_likelihood_tau_round1_bcf(NumericVector y_temp,NumericMatrix treetable_tau ,
                                           NumericMatrix obs_to_nodes_tau,
                                           int n,double a_mu,double a_tau,double nu,double lambda, NumericVector z){
  
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec	
  //make W_bcf and mu matrices for the sum of trees
  
  arma::mat Wmat_tau=W_tauround1_bcf(treetable_tau,obs_to_nodes_tau,n);	// Create W_bcf matrix. Function defined on line 829
  //Rcout << "Wmat_tau works.\n";
  
  double b_tau=Wmat_tau.n_cols;
  Wmat_tau.each_col()%=z_ar;
  arma::mat Wmat = Wmat_tau;
  
  double b=Wmat.n_cols;									// b is number of columns of W_bcf matrix (omega in the paper)
  arma::vec yvec=Rcpp::as<arma::vec>(y_temp);			// convert input y_temp to arma vec called yvec
  arma::mat y(n,1);										// create a matrix with n (number of observations) rows and 1 column
  y.col(0)=yvec;										// set first column of y equal to yvec
  //get exponent
  double expon=(n+nu+b)/2;								// set the expoenent (equation 5 in the paper)
  //get y^Tpsi^{-1}y
  // arma::mat psi_inv=psi.i();
  arma::mat yty=y.t()*y;								// yty = y transpose y (sum of squares)
  
  //get t(y)inv(psi)J_bcf
  arma::mat ytW=y.t()*Wmat;								// y transpose W_bcf
  //get t(J_bcf)inv(psi)J_bcf  
  arma::mat WtW=Wmat.t()*Wmat;							// W_bcf transpose W_bcf
  //get jpsij +aI
  arma::mat aI(b,b);									// create b by b matrix called aI. NOT INIIALIZED. 
  aI=aI.eye();										// a times b by b identity matrix. The .eye() turns aI into an identity matrix.
  //arma::vec a_vec_mu = a_mu*arma::ones<arma::vec>(b_mu);
  arma::vec a_vec_tau = a_tau*arma::ones<arma::vec>(b_tau);
  arma::vec a_vec(b);
  //a_vec.head(b_mu) = a_vec_mu;
  a_vec.tail(b_tau) = a_vec_tau;
  aI.diag() = a_vec;
  
  arma::mat sec_term=WtW+aI;							//
  arma::mat sec_term_inv=sec_term.i();					// matrix inverse expression in middle of eq 5 in the paper. The .i() obtains the matrix inverse.
  //get t(J_bcf)inv(psi)y
  arma::mat third_term=Wmat.t()*y;						// W_bcf transpose Y
  //get m^TV^{-1}m
  arma::mat mvm= ytW*sec_term_inv*third_term;			// matrix expression in middle of equation 5
  arma::mat rel=-expon*log(nu*lambda - mvm +yty);		// log of all of equation 5 (i.e. the log of the marginal likelihood of the sum of tree model)
  double rel2=as<double>(wrap(rel));					// convert to double (from 1 by 1 arma mat)
  return(rel2);											// return the log marginal likelihood
}    
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_mu_bcf(NumericVector resids,arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,
                           double a,double mu,double nu,double lambda,double c,double lowest_BIC,
                           int parent,NumericMatrix cp_mat,double alpha,double beta,int maxOWsize//,int first_round
                             ){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;													// create integer variable (not initialized)
  NumericMatrix treetable_c=treetable;							// copy the input matrix, call it treetable_c
  NumericMatrix treemat_c=tree_mat;								// copy the input matrix, call it treemat_c
  
  NumericVector terminal_nodes=find_term_nodes_bcf(treetable_c);		// terminal_nodes gives indices of treetable_c elements that correspond to term_nodes.
  IntegerVector change_node1;										// cerates an integer vector. THIS VECTOR IS NOT USED IN THIS FUNCTION
  int list_size=1000;												// create an integer initialied equal to 1000. Why 1000?
  std::vector<double> tree_lik(list_size);						// create a vector tree_lik of length 1000.
  List proposal_tree;												// create a list called proposal_tree
  List ret(9);													// create a list ret of length 9
  bool no_tree_err=0;												// create a bool variable no_tree_err initialized equal to 0 (FALSE)
  List likeliest_tree;											// create a list likeliest_tree
  List tree_list(list_size);										// create a list, tree_list, of length 1000
  List tree_mat_list(list_size);									// create a list, tree_mat_list, of length 1000
  int count=0;													// create a variabke count, initialize equal to 0.
  std::vector<int> tree_parent(list_size);						// create a vector, tree_parent, of length 1000
  int best_sv;													// create a variable best_sv. Not initialized
  double best_sp;													// create a variable best_sp. Not initialized
  double tree_prior=0;											// create a variable tree_prior. Initialized equal to 0.
  List changetree;												// create a list changetree
  double BIC;														// create a variable BIC
  int p;															// create a variable p
  List eval_model;												// create a list eval_model
  NumericVector int_nodes;										// create a vector int_nodes
  arma::colvec curr_col=data.col(0);										// Let the arma colvec, curr_col, equal the 1st column of the input matrix data
  arma::uvec grow_obs=find_term_obs_bcf(treemat_c,terminal_nodes[0]);			// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[0] (for leftmost column of treemat_c that has elements equal to terminal_nodes[0]).
  NumericVector d1=unique(find_term_cols_bcf(treemat_c,terminal_nodes[0]));	// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[0]. Unique funtion removes duplicated of columns. Unique also sorts descending. Why not IntegerVector
  arma::mat data_curr_node=data.rows(grow_obs);				// matrix consisting of first grow_obs rows of data. 
  double d=d1[0];															// index of rightmost column of treemat_c with at least one element equal to terminal_nodes[0]
  NumericVector get_min=get_grow_obs_bcf(data,wrap(grow_obs),cp_mat(0,0)+1);	// obtain the elements of the cp_mat(0,0)+1^th column of data that are indexed by grow_obs
  double lik;																// create a variable called lik. Not initialized.
  
  for(int l=0;l<terminal_nodes.size();l++){										//	vector of length equal to that of terminal_nodes
    //loop over each terminal node												//
    grow_obs=find_term_obs_bcf(treemat_c,terminal_nodes[l]);						// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[l] (letter l) (for leftmost column of treemat_c that has elements equal to terminal_nodes[l] (letter l)).
    //depth of tree at current terminal node									//
    d1=unique(find_term_cols_bcf(treemat_c,terminal_nodes[l]));						// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[l] (letter l). Unique funtion removes duplcated of columns. Unique also sorts descending.
    data_curr_node=data.rows(grow_obs);								// matrix consisting of first grow_obs rows of data. Note grow_obs changed on line 926, therefore not duplicating line 919.
    d=d1[0];																	// index of rightmost column of treemat_c with at least one element equal to terminal_nodes[l] (letter l)
    int w=cp_mat.nrow();														// w is number of rows of cp_mat
    if(data_curr_node.n_rows<=2){												// if data_curr_node has 2 rows or less.
      throw std::range_error("not enough obs in node to grow any further");	// throw an error message. Not enough observations.
      //continue;
    }
    for(int k=0;k<w;k++){														// loop of length w, the number of rows of cp_mat
      split_var=cp_mat(k,0)+1;												// split_var is k+1^th row, 1st column, of cp_mat, +1
      arma::colvec curr_cols=data.col(split_var-1);							// curr_cols is the split_var^tgh column of data
      get_min=get_grow_obs_bcf(data,wrap(grow_obs),split_var);					// obtain the elements of the split_var^th column of data that are indexed by grow_obs 
      
      //Removing unnecessary lines
      //if(get_min.size()<=2){													// If get_min has 2 or less observations. (too few variables to split on?)
      //  throw std::range_error("obs in this terminal node are too small");	//
      //}
      
      double split_point=cp_mat(k,1);											// variable split_point equals element in k+1^th row 2nd column of cp_mat
      arma::vec curr_cols2=data_curr_node.col(split_var-1);					// curr_cols2 is split_var^th column of data_curr_node
      //arma::vec get_min_a=Rcpp::as<arma::vec>(get_min);		// converts to arma vec
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));	// ld_prop is elements of curr_cols2 <= split_point
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2 > split_point));		// rd_prop is elements of curr_cols2 > split_point
      
      if(ld_prop.size()<=2 || rd_prop.size()<=2){									// if 2 or less observations in either ld_prop or rd_prop
        continue;																// skip to next iteration of the loop
      }
      proposal_tree=grow_tree_bcf(data,resids,treemat_c,terminal_nodes[l],treetable_c,split_var,split_point,terminal_nodes,wrap(grow_obs),d,get_min,data_curr_node); // elaborate function on line 469. creates list of 2 elements, tree matrix and tree table. Appears to grow node indexed by terminal_nodes[l] (letter l) (?).
      
      
      //NumericMatrix test =proposal_tree[0];											// first element is tree table
      //NumericMatrix test1 =proposal_tree[1];											// second element is tree matrix
      
      //if(test1.ncol()==3){															// If tree matrix has 3 columns
      //  NumericVector u1=unique(test1(_,0));										// set u1 equal to the unique (ordered descending) elements of 1st column of tree matrix
      //  NumericVector u2=unique(test1(_,1));										// set u2 equal to the unique (ordered descending) elements of 2nd column of tree matrix
      //  NumericVector u3=unique(test1(_,2));										// set u3 equal to the unique (ordered descending) elements of 3rd column of tree matrix
      //}
      
      
      
      
      //if(first_round==1){																	// If input value first_round equals 1 (number one)
        lik=likelihood_function_bcf(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);	// set lik equal to tree likelihood defined on line 201.
        
      //}else{
        //have a sum of trees
      //  lik=likelihood_function_bcf(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);	// Same as line in if-statement above. What is the purpose of the if-statement. set lik equal to tree likelihood defined on line 201.
      //}
      
      
      //NumericMatrix temptestingtabcols = proposal_tree[0];
      //if(temptestingtabcols.ncol()<5) throw std::range_error("Line 1021");
      
      
      tree_prior=get_tree_prior_bcf(proposal_tree[0],proposal_tree[1],alpha,beta);	// defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
      int_nodes=find_term_nodes_bcf(proposal_tree[0]);							// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector?
      p=int_nodes.size();														// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
      BIC=-2*(lik+log(tree_prior))+p*log(data.n_rows);						// data.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      //BIC=-2*(lik)+p*log(data.n_rows);						// data.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      if(BIC<lowest_BIC){														// if statement for updating lowest BIC...etc.
        lowest_BIC=BIC;														// update (input variable) lowest_BIC
        best_sv=split_var;													// set a value for, or update best_sv
        best_sp=split_point;												// set a value for, or update split_point
        likeliest_tree=proposal_tree;										// set a value for, or update likeliest_tree
        tree_list[count]=proposal_tree[0];									// add an element to the list of tree tables
        tree_mat_list[count]=proposal_tree[1];								// add an element to the list of tree matrices
        tree_lik[count]=BIC;												// add an element to the vector of tree liklihoods (BICs)
        tree_parent[count]=parent;											// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
        count++;															// increase the count
        if(count==(tree_list.size()-1)){									// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
          list_size=list_size*2;											// multiply list size by 2
          tree_list=resize_bigger_bcf(tree_list,list_size);					// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);			// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_lik.resize(list_size);										// increase tree_lik to twice its previous size
          tree_parent.resize(list_size);									// increase tree_parent to twice its previous size
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){											// If in Occam's window (but not the new minimum BIC)
          if(is<NumericMatrix>(proposal_tree[0])){						// If proposal_tree[0] is a Numeric Matrix, the do nothing
          }else{															// If proposal_tree[0] is NOT a NumericMatrix
            throw std::range_error("proposal tree not a matrix");			// Then throw an error
          }
          tree_list[count]=proposal_tree[0];								// add an element to the list of tree tables
          tree_mat_list[count]=proposal_tree[1];							// add an element to the list of tree matrices
          tree_lik[count]=BIC;											// add an element to the vector of tree liklihoods (BICs)
          tree_parent[count]=parent;										// add an element to the vector tree_parent. NOTE: All elements will be equal to the inout value of parent which is not changed in this function.
          count++;														// increase the count
          if(count==(tree_list.size()-1)){								// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
            list_size=list_size*2;										// multiply list size by 2
            tree_list=resize_bigger_bcf(tree_list,list_size);				// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);		// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_lik.resize(list_size);									// increase tree_lik to twice its previous size
            tree_parent.resize(list_size);								// increase tree_parent to twice its previous size
          }
        }
      }
    }  
  }
  tree_list=resize_bcf(tree_list,count);				// remove the values in tree_list that aren't filled in
  tree_mat_list=resize_bcf(tree_mat_list,count);		// remove the values in tree_mat_list that aren't filled in
  tree_lik.resize(count);							// remove the values in tree_lik that aren't filled in
  tree_parent.resize(count);						// remove the values in tree_parent that aren't filled in
  if(count>0){									// If these lists are nonempty
    eval_model=evaluate_model_occams_window_bcf(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    NumericVector testlik =eval_model[0];		// testlik is tree_lik after removing models outside Occam's window
    List testtree =eval_model[1];				// testtree is tree_list after removing models outside Occam's window
    List testmat =eval_model[2];				// testmat is tree_mat_list after removing models outside Occam's window
    IntegerVector testpar =eval_model[3];		// testpar is tree_parent after removing models outside Occam's window
    
    if(testlik.size()>0){									// If a nonzero number of models remain in Occam's window
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){						// maxOWsize is an input variable
        IntegerVector owindices=orderforOW__bcf(testlik);		// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on.
        owindices=owindices-1;							// Presumably the match function in orderforOW__bcf gives indices beginning at 1, and therefore 1 must be taken away from all index values.
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);				// create vector temp_olik of size maxOWsize
        List temp_otrees(maxOWsize);					// create List temp_otrees of size maxOWsize
        List temp_omat(maxOWsize);						// create List temp_omat of size maxOWsize
        IntegerVector temp_oparent(maxOWsize);			// create IntegerVector temp_oparent of size maxOWsize
        for(int t=0;t<maxOWsize;t++){					// loop of length maxOWsize
          temp_olik[t]=testlik[owindices[t]];			// temp_olik is maxOWsize largest BIC elements of testlik ordered by BIC (descending?)
          temp_otrees[t]=testtree[owindices[t]];		// temp_otrees is maxOWsize largest BIC elements of testtree ordered by BIC (descending?)
          temp_omat[t]=testmat[owindices[t]];			// temp_omat is maxOWsize largest BIC elements of testmat ordered by BIC (descending?)
          temp_oparent[t]=testpar[owindices[t]];		// temp_oparent is maxOWsize largest BIC elements of testpar ordered by BIC (descending?)
        }
        testlik=temp_olik;			// reset testlik equal to temp_olik
        testtree=temp_otrees;		// reset testtree equal to temp_otrees
        testmat=temp_omat;			// reset testmat equal to temp_omat
        testpar=temp_oparent;		// reset testpar equal to temp_oparent
      }
      ret[0]=lowest_BIC;				// first element of output list. Lowest BIC
      ret[1]=best_sv;					// second element of output list. Best splitting variable
      ret[2]=best_sp;					// third element of ouput list. Best splitting point
      ret[3]=likeliest_tree;			// fourth element of output list. List containing treee table and tree matrix of lowest BIC tree
      ret[4]=testtree;				// fifth element of output list. List of tree tables
      ret[5]=testlik;					// sixth element of output list. Vector of BICs
      ret[6]=testmat;					// seventh element of output list. List of tree matrices
      ret[7]=testpar;					// eighth element of output list. Vector with all elements equal to the input value of parent?
      ret[8]=no_tree_err;				// ninth element of output list. Boolean equal to false
      
      return (ret);					// return the list ret
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;					// Boolean equal to true.
      List gr(1);						// Create a list, gr, of length one.
      gr[0]=no_tree_err;				// First element of gr is the Boolean no_tree_err.
      return(gr);						// Return the list gr.
    }
  }else{
    no_tree_err=1;						// Boolean equal to true.
    List gr(1);							// Create a list, gr, of length one.
    gr[0]=no_tree_err;					// First element of gr is the Boolean no_tree_err.
    return(gr);							// Return the list gr.
  }
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_bcf(NumericVector resids,arma::mat& data,NumericMatrix treetable,NumericMatrix tree_mat,double a,double mu,double nu,double lambda,double c,double lowest_BIC,int parent,NumericMatrix cp_mat,double alpha,double beta,int maxOWsize,int first_round){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;													// create integer variable (not initialized)
  NumericMatrix treetable_c=treetable;							// copy the input matrix, call it treetable_c
  NumericMatrix treemat_c=tree_mat;								// copy the input matrix, call it treemat_c
  
  NumericVector terminal_nodes=find_term_nodes_bcf(treetable_c);		// terminal_nodes gives indices of treetable_c elements that correspond to term_nodes.
  IntegerVector change_node1;										// cerates an integer vector. THIS VECTOR IS NOT USED IN THIS FUNCTION
  int list_size=1000;												// create an integer initialied equal to 1000. Why 1000?
  std::vector<double> tree_lik(list_size);						// create a vector tree_lik of length 1000.
  List proposal_tree;												// create a list called proposal_tree
  List ret(9);													// create a list ret of length 9
  bool no_tree_err=0;												// create a bool variable no_tree_err initialized equal to 0 (FALSE)
  List likeliest_tree;											// create a list likeliest_tree
  List tree_list(list_size);										// create a list, tree_list, of length 1000
  List tree_mat_list(list_size);									// create a list, tree_mat_list, of length 1000
  int count=0;													// create a variabke count, initialize equal to 0.
  std::vector<int> tree_parent(list_size);						// create a vector, tree_parent, of length 1000
  int best_sv;													// create a variable best_sv. Not initialized
  double best_sp;													// create a variable best_sp. Not initialized
  double tree_prior=0;											// create a variable tree_prior. Initialized equal to 0.
  List changetree;												// create a list changetree
  double BIC;														// create a variable BIC
  int p;															// create a variable p
  List eval_model;												// create a list eval_model
  NumericVector int_nodes;										// create a vector int_nodes
  arma::colvec curr_col=data.col(0);										// Let the arma colvec, curr_col, equal the 1st column of the input matrix data
  arma::uvec grow_obs=find_term_obs_bcf(treemat_c,terminal_nodes[0]);			// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[0] (for leftmost column of treemat_c that has elements equal to terminal_nodes[0]).
  NumericVector d1=unique(find_term_cols_bcf(treemat_c,terminal_nodes[0]));	// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[0]. Unique funtion removes duplicated of columns. Unique also sorts descending. Why not IntegerVector
  arma::mat data_curr_node=data.rows(grow_obs);				// matrix consisting of first grow_obs rows of data.
  double d=d1[0];															// index of rightmost column of treemat_c with at least one element equal to terminal_nodes[0]
  NumericVector get_min=get_grow_obs_bcf(data,wrap(grow_obs),cp_mat(0,0)+1);	// obtain the elements of the cp_mat(0,0)+1^th column of data that are indexed by grow_obs
  double lik;																// create a variable called lik. Not initialized.
  
  for(int l=0;l<terminal_nodes.size();l++){										//	vector of length equal to that of terminal_nodes
    //loop over each terminal node												//
    grow_obs=find_term_obs_bcf(treemat_c,terminal_nodes[l]);						// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[l] (letter l) (for leftmost column of treemat_c that has elements equal to terminal_nodes[l] (letter l)).
    //depth of tree at current terminal node									//
    d1=unique(find_term_cols_bcf(treemat_c,terminal_nodes[l]));						// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[l] (letter l). Unique funtion removes duplcated of columns. Unique also sorts descending.
    data_curr_node=data.rows(grow_obs);								// matrix consisting of first grow_obs rows of data. Note grow_obs changed on line 926, therefore not duplicating line 919.
    d=d1[0];																	// index of rightmost column of treemat_c with at least one element equal to terminal_nodes[l] (letter l)
    int w=cp_mat.nrow();														// w is number of rows of cp_mat
    if(data_curr_node.n_rows<=2){												// if data_curr_node has 2 rows or less.
      throw std::range_error("not enough obs in node to grow any further Line 1169");	// throw an error message. Not enough observations.
      //continue;
    }
    for(int k=0;k<w;k++){														// loop of length w, the number of rows of cp_mat
      split_var=cp_mat(k,0)+1;												// split_var is k+1^th row, 1st column, of cp_mat, +1
      arma::colvec curr_cols=data.col(split_var-1);							// curr_cols is the split_var^tgh column of data
      get_min=get_grow_obs_bcf(data,wrap(grow_obs),split_var);					// obtain the elements of the split_var^th column of data that are indexed by grow_obs 
      
      //Removing unnecessary lines
      //if(get_min.size()<=2){													// If get_min has 2 or less observations. (too few variables to split on?)
      //  throw std::range_error("obs in this terminal node are too small");	//
      //}
      
      double split_point=cp_mat(k,1);											// variable split_point equals element in k+1^th row 2nd column of cp_mat
      arma::vec curr_cols2=data_curr_node.col(split_var-1);					// curr_cols2 is split_var^th column of data_curr_node
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));	// ld_prop is elements of curr_cols2 <= split_point
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));		// rd_prop is elements of curr_cols2 > split_point
      
      if(ld_prop.size()<=2 || rd_prop.size()<=2){									// if 2 or less observations in either ld_prop or rd_prop
        continue;																// skip to next iteration of the loop
      }
      proposal_tree=grow_tree_bcf(data,resids,treemat_c,terminal_nodes[l],treetable_c,split_var,split_point,terminal_nodes,wrap(grow_obs),d,get_min,data_curr_node); // elaborate function on line 469. creates list of 2 elements, tree matrix and tree table. Appears to grow node indexed by terminal_nodes[l] (letter l) (?).
      
      
      
      //NumericMatrix test =proposal_tree[0];											// first element is tree table
      //NumericMatrix test1 =proposal_tree[1];											// second element is tree matrix
      
      //if(test1.ncol()==3){															// If tree matrix has 3 columns
      //  NumericVector u1=unique(test1(_,0));										// set u1 equal to the unique (ordered descending) elements of 1st column of tree matrix
      //  NumericVector u2=unique(test1(_,1));										// set u2 equal to the unique (ordered descending) elements of 2nd column of tree matrix
      //  NumericVector u3=unique(test1(_,2));										// set u3 equal to the unique (ordered descending) elements of 3rd column of tree matrix
      //}
      
      
      
      
      if(first_round==1){																	// If input value first_round equals 1 (number one)
        lik=likelihood_function_bcf(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);	// set lik equal to tree likelihood defined on line 201.
        
      }else{
        //have a sum of trees
        lik=likelihood_function_bcf(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);	// Same as line in if-statement above. What is the purpose of the if-statement. set lik equal to tree likelihood defined on line 201.
      }
      NumericMatrix temptestingtabcols = proposal_tree[0];
      if(temptestingtabcols.ncol()<5) throw std::range_error("Line 1208");
      tree_prior=get_tree_prior_bcf(proposal_tree[0],proposal_tree[1],alpha,beta);	// defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
      int_nodes=find_term_nodes_bcf(proposal_tree[0]);							// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector?
      p=int_nodes.size();														// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
      BIC=-2*(lik+log(tree_prior))+p*log(data.n_rows);						// data.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      //BIC=-2*(lik)+p*log(data.n_rows);						// data.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      if(BIC<lowest_BIC){														// if statement for updating lowest BIC...etc.
        lowest_BIC=BIC;														// update (input variable) lowest_BIC
        best_sv=split_var;													// set a value for, or update best_sv
        best_sp=split_point;												// set a value for, or update split_point
        likeliest_tree=proposal_tree;										// set a value for, or update likeliest_tree
        tree_list[count]=proposal_tree[0];									// add an element to the list of tree tables
        tree_mat_list[count]=proposal_tree[1];								// add an element to the list of tree matrices
        tree_lik[count]=BIC;												// add an element to the vector of tree liklihoods (BICs)
        tree_parent[count]=parent;											// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
        count++;															// increase the count
        if(count==(tree_list.size()-1)){									// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
          list_size=list_size*2;											// multiply list size by 2
          tree_list=resize_bigger_bcf(tree_list,list_size);					// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);			// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_lik.resize(list_size);										// increase tree_lik to twice its previous size
          tree_parent.resize(list_size);									// increase tree_parent to twice its previous size
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){											// If in Occam's window (but not the new minimum BIC)
          if(is<NumericMatrix>(proposal_tree[0])){						// If proposal_tree[0] is a Numeric Matrix, the do nothing
          }else{															// If proposal_tree[0] is NOT a NumericMatrix
            throw std::range_error("proposal tree not a matrix");			// Then throw an error
          }
          tree_list[count]=proposal_tree[0];								// add an element to the list of tree tables
          tree_mat_list[count]=proposal_tree[1];							// add an element to the list of tree matrices
          tree_lik[count]=BIC;											// add an element to the vector of tree liklihoods (BICs)
          tree_parent[count]=parent;										// add an element to the vector tree_parent. NOTE: All elements will be equal to the inout value of parent which is not changed in this function.
          count++;														// increase the count
          if(count==(tree_list.size()-1)){								// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
            list_size=list_size*2;										// multiply list size by 2
            tree_list=resize_bigger_bcf(tree_list,list_size);				// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);		// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_lik.resize(list_size);									// increase tree_lik to twice its previous size
            tree_parent.resize(list_size);								// increase tree_parent to twice its previous size
          }
        }
      }
    }  
  }
  tree_list=resize_bcf(tree_list,count);				// remove the values in tree_list that aren't filled in
  tree_mat_list=resize_bcf(tree_mat_list,count);		// remove the values in tree_mat_list that aren't filled in
  tree_lik.resize(count);							// remove the values in tree_lik that aren't filled in
  tree_parent.resize(count);						// remove the values in tree_parent that aren't filled in
  if(count>0){									// If these lists are nonempty
    eval_model=evaluate_model_occams_window_bcf(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    NumericVector testlik =eval_model[0];		// testlik is tree_lik after removing models outside Occam's window
    List testtree =eval_model[1];				// testtree is tree_list after removing models outside Occam's window
    List testmat =eval_model[2];				// testmat is tree_mat_list after removing models outside Occam's window
    IntegerVector testpar =eval_model[3];		// testpar is tree_parent after removing models outside Occam's window
    
    if(testlik.size()>0){									// If a nonzero number of models remain in Occam's window
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){						// maxOWsize is an input variable
        IntegerVector owindices=orderforOW__bcf(testlik);		// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on.
        owindices=owindices-1;							// Presumably the match function in orderforOW__bcf gives indices beginning at 1, and therefore 1 must be taken away from all index values.
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);				// create vector temp_olik of size maxOWsize
        List temp_otrees(maxOWsize);					// create List temp_otrees of size maxOWsize
        List temp_omat(maxOWsize);						// create List temp_omat of size maxOWsize
        IntegerVector temp_oparent(maxOWsize);			// create IntegerVector temp_oparent of size maxOWsize
        for(int t=0;t<maxOWsize;t++){					// loop of length maxOWsize
          temp_olik[t]=testlik[owindices[t]];			// temp_olik is maxOWsize largest BIC elements of testlik ordered by BIC (descending?)
          temp_otrees[t]=testtree[owindices[t]];		// temp_otrees is maxOWsize largest BIC elements of testtree ordered by BIC (descending?)
          temp_omat[t]=testmat[owindices[t]];			// temp_omat is maxOWsize largest BIC elements of testmat ordered by BIC (descending?)
          temp_oparent[t]=testpar[owindices[t]];		// temp_oparent is maxOWsize largest BIC elements of testpar ordered by BIC (descending?)
        }
        testlik=temp_olik;			// reset testlik equal to temp_olik
        testtree=temp_otrees;		// reset testtree equal to temp_otrees
        testmat=temp_omat;			// reset testmat equal to temp_omat
        testpar=temp_oparent;		// reset testpar equal to temp_oparent
      }
      ret[0]=lowest_BIC;				// first element of output list. Lowest BIC
      ret[1]=best_sv;					// second element of output list. Best splitting variable
      ret[2]=best_sp;					// third element of ouput list. Best splitting point
      ret[3]=likeliest_tree;			// fourth element of output list. List containing treee table and tree matrix of lowest BIC tree
      ret[4]=testtree;				// fifth element of output list. List of tree tables
      ret[5]=testlik;					// sixth element of output list. Vector of BICs
      ret[6]=testmat;					// seventh element of output list. List of tree matrices
      ret[7]=testpar;					// eighth element of output list. Vector with all elements equal to the input value of parent?
      ret[8]=no_tree_err;				// ninth element of output list. Boolean equal to false
      
      return (ret);					// return the list ret
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;					// Boolean equal to true.
      List gr(1);						// Create a list, gr, of length one.
      gr[0]=no_tree_err;				// First element of gr is the Boolean no_tree_err.
      return(gr);						// Return the list gr.
    }
  }else{
    no_tree_err=1;						// Boolean equal to true.
    List gr(1);							// Create a list, gr, of length one.
    gr[0]=no_tree_err;					// First element of gr is the Boolean no_tree_err.
    return(gr);							// Return the list gr.
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_tau_bcf(NumericVector resids,arma::mat& x_moderate_a,
                            NumericMatrix tree_table_tau,NumericMatrix tree_mat_tau,
                                   double a_mu,double a_tau,double mu_tau,double nu,double lambda,
                                   double c,
                                   double lowest_BIC,
                                   int parent,
                                   NumericMatrix cp_mat,
                                   double alpha_mu,double beta_mu,double alpha_tau,double beta_tau,
                                   int maxOWsize,//int first_round,
                                   NumericVector z){
  
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;																// create an integer split_var. Not initialized
  NumericMatrix tree_table_tau_c=tree_table_tau;										// create a matrix treetable_c equal to the input matrix.
  NumericMatrix tree_mat_tau_c=tree_mat_tau;											// create a matric tree_mat_tau_c equal to the input matrix.
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  
  NumericVector terminal_nodes=find_term_nodes_bcf(tree_table_tau_c);					// find term nodes function defined line 168. Gives index of values of treetable_c that are term nodes (indices from 1 to length of vector). Why not integer vector?
  //Rcout << "terminal_nodes[0] equals " << terminal_nodes[0] << ".\n" ;
  //Rcout << "length terminal_nodes equals " << terminal_nodes.size() << ".\n" ;
  
  IntegerVector change_node1;													// create an IntegerVector. Initially all values are 0, but size not given.
  int list_size=1000;															// create a list of size 1000. Why 1000.
  std::vector<double> tree_lik(list_size);									// create a vector of size 1000
  List proposal_tree;															// create an empty list
  List ret(9);																// create a list with 9 elements
  bool no_tree_err=0;															// create a boolean initialized equal to 0 (false).
  List likeliest_tree;														// create an empty list
  List tree_list(list_size);													// create a list of size 1000
  List tree_mat_list(list_size);												// create a list of size 1000
  int count=0;																// create an integer variable, initialized equal to 0.
  std::vector<int> tree_parent(list_size);									// create an integer vector of size 1000
  int best_sv;																// create an integer vector. Not initialized.
  double best_sp;																// create a double variable
  double tree_prior=1;														// create a double variable, initialized equal to 0.
  List changetree;															// create an empty list.
  double BIC;																	// create a double variable. Not initialized.
  int p;																		// create an integer variable. Not initialized.
  //int p_other_mu=0;
  //int p_other_tau=0;
  List eval_model;															// create an empty list.
  NumericVector int_nodes;													// create a NumericVector. Initially all values are 0, but size not given.
  //NumericVector other_int_nodes_mu;
  //NumericVector other_int_nodes_tau;
  arma::colvec curr_col=x_moderate_a.col(0);											// let curr_col be an arma colvec equal to the first column of the inut arma mat x_moderate_a.
  arma::uvec grow_obs=find_term_obs_bcf(tree_mat_tau_c,terminal_nodes[0]);				// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[0] (for leftmost column of tree_mat_tau_c that has elements equal to terminal_nodes[0]).
  //Rcout << "length of grow_obs equals " << grow_obs.n_elem<< ".\n";
  NumericVector d1=unique(find_term_cols_bcf(tree_mat_tau_c,terminal_nodes[0]));		// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[0]. Unique funtion removes duplicated of columns. Unique also sorts descending. Why not IntegerVector
  arma::mat data_curr_node=x_moderate_a.rows(grow_obs);					// matrix consisting of first grow_obs rows of x_moderate_a.
  double d=d1[0];																// index of rightmost column of tree_mat_tau_c with at least one element equal to terminal_nodes[0]
  NumericVector get_min=get_grow_obs_bcf(x_moderate_a,wrap(grow_obs),cp_mat(0,0)+1);		// obtain the elements of the cp_mat(0,0)+1^th column of x_moderate_a that are indexed by grow_obs
  double lik;																	// create a variable called lik. Not initialized.
  
  for(int l=0;l<terminal_nodes.size();l++){										//	vector of length equal to that of terminal_nodes
    //loop over each terminal node
    grow_obs=find_term_obs_bcf(tree_mat_tau_c,terminal_nodes[l]);						// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[l] (letter l) (for leftmost column of tree_mat_tau_c that has elements equal to terminal_nodes[l] (letter l)).
    //depth of tree at current terminal node
    d1=unique(find_term_cols_bcf(tree_mat_tau_c,terminal_nodes[l]));						// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[l] (letter l). Unique funtion removes duplcated of columns. Unique also sorts descending.
    data_curr_node=x_moderate_a.rows(grow_obs);								// matrix consisting of first grow_obs rows of x_moderate_a. Note grow_obs changed on line 926, therefore not duplicating line 919.
    d=d1[0];																	// index of rightmost column of tree_mat_tau_c with at least one element equal to terminal_nodes[l] (letter l)
    int w=cp_mat.nrow();														// w is number of rows of cp_mat
    if(data_curr_node.n_rows<=2){												// if data_curr_node has 2 rows or less.
      //Rcout << "terminal_nodes[l] equals " << terminal_nodes[l]<< ".\n" ;
      //Rcout << "num obs in grow_obs" << grow_obs.n_elem<< ".\n" ;
      //Rcout << " num obs in data_curr_node" << data_curr_node.n_rows<< ".\n" ;
      //Rcout << " iteration number " << l<< ".\n" ;
      throw std::range_error("not enough obs in node to grow any further Line 1360");	// throw an error message. Not enough observations.
      //continue;
    }
    //Rcout << " iteration number " << l<< ".\n" ;
    
    for(int k=0;k<w;k++){														// loop of length w, the number of rows of cp_mat
      //p_other_mu=0;
      //Rcout << "inner iteration number " << k<< ".\n" ;
      split_var=cp_mat(k,0)+1;												// split_var is k+1^th row, 1st column of cp_mat +1
      arma::colvec curr_cols=x_moderate_a.col(split_var-1);							// curr_cols is the split_var^th column of x_moderate_a
      get_min=get_grow_obs_bcf(x_moderate_a,wrap(grow_obs),split_var);					// obtain the elements of the split_var^th column of x_moderate_a that are indexed by grow_obs 
      NumericVector z_growvec = get_grow_obs_in_z_bcf(z_ar, wrap(grow_obs));
      
      //Removing unnecessary lines
      //if(get_min.size()<=2){													// If get_min has 2 or less observations. (too few variables to split on?)
      //  throw std::range_error("obs in this terminal node are too small");
      //}
      
      double split_point=cp_mat(k,1);											// variable split_point equals element in k+1^th row 2nd column of cp_mat
      arma::vec curr_cols2=data_curr_node.col(split_var-1);					// curr_cols2 is split_var^th column of data_curr_node
      arma::vec get_min_a=Rcpp::as<arma::vec>(get_min);		// converts to arma vec
      arma::vec z_growvec_a=Rcpp::as<arma::vec>(z_growvec);		// converts to arma vec
      arma::vec get_min_a_treated = get_min_a.elem(arma::find(z_growvec_a==1));	
      
      arma::vec ld_prop=get_min_a_treated.elem(arma::find(get_min_a_treated <= split_point));	// ld_prop is elements of curr_cols2 <= split_point
      arma::vec rd_prop=get_min_a_treated.elem(arma::find(get_min_a_treated > split_point));		// rd_prop is elements of curr_cols2 > split_point
      
      if(ld_prop.size()<=2 || rd_prop.size()<=2){									// if 2 or less observations in either ld_prop or rd_prop
        continue;																// skip to next iteration of the loop
      }
      
      //Rcout << "error after new code .\n" ;
      
      
      proposal_tree=grow_tree_bcf(x_moderate_a,resids,tree_mat_tau_c,terminal_nodes[l],tree_table_tau_c,split_var,split_point,terminal_nodes,wrap(grow_obs),d,get_min,data_curr_node);	// elaborate function on line 469. creates list of 2 elements, tree matrix and tree table. Appears to grow node indexed by terminal_nodes[l] (letter l) (?).
      
      
      
      //NumericMatrix test =proposal_tree[0];										// first element is tree table
      //NumericMatrix test1 =proposal_tree[1];										// second element is tree matrix
      
      //if(test1.ncol()==3){														// If tree matrix has 3 columns
      //  NumericVector u1=unique(test1(_,0));									// set u1 equal to the unique (ordered descending) elements of 1st column of tree matrix
      //  NumericVector u2=unique(test1(_,1));									// set u2 equal to the unique (ordered descending) elements of 2nd column of tree matrix
      //  NumericVector u3=unique(test1(_,2));									// set u3 equal to the unique (ordered descending) elements of 3rd column of tree matrix
      //}
      
      
      
      //if(first_round==1){																// If input value first_round equals 1 (number one)
        lik=sumtree_likelihood_tau_round1_bcf(resids,proposal_tree[0],proposal_tree[1],resids.size(),a_mu,a_tau,nu,lambda,z);
        
      //}else{
       // throw std::range_error("get_best_split_tau_bcf should only be used in the first round");	// throw an error.
      //}
      //at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      //NumericMatrix temptestingtabcols = proposal_tree[0];
      //if(temptestingtabcols.ncol()<5) throw std::range_error("Line 1021");
      
      
      tree_prior=get_tree_prior_bcf(proposal_tree[0],proposal_tree[1],alpha_tau,beta_tau);	// defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
      int_nodes=find_term_nodes_bcf(proposal_tree[0]);							// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector?
      p=int_nodes.size();														// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
      BIC=-2*(lik+log(tree_prior))+p*log(x_moderate_a.n_rows);	
      
      //BIC=-2*(lik)+(p_other_mu+p)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      if(BIC<lowest_BIC){											// if statement for updating lowest BIC...etc.
        lowest_BIC=BIC;											// update (input variable) lowest_BIC
        best_sv=split_var;										// set a value for, or update best_sv
        best_sp=split_point;									// set a value for, or update split_point
        likeliest_tree=proposal_tree;							// set a value for, or update likeliest_tree
        tree_list[count]=proposal_tree[0];						// add an element to the list of tree tables
        tree_mat_list[count]=proposal_tree[1];					// add an element to the list of tree matrices
        tree_lik[count]=BIC;									// add an element to the vector of tree liklihoods (BICs)
        tree_parent[count]=parent;								// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
        count++;												// increase the count
        if(count==(tree_list.size()-1)){						// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
          list_size=list_size*2;									// multiply list size by 2
          tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
          tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){										// If in Occam's window (but not the new minimum BIC)
          if(is<NumericMatrix>(proposal_tree[0])){					// If proposal_tree[0] is a Numeric Matrix, do nothing
            //std::cout<<"its a matrix "<<"\n";						// if this is used, it should say "it's a matrix".
          }else{														// If proposal_tree[0] is NOT a NumericMatrix
            throw std::range_error("proposal tree not a matrix");	// Then throw an error
          }
          tree_list[count]=proposal_tree[0];							// add an element to the list of tree tables
          tree_mat_list[count]=proposal_tree[1];						// add an element to the list of tree matrices
          tree_lik[count]=BIC;										// add an element to the vector of tree liklihoods (BICs)
          tree_parent[count]=parent;									// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
          count++;													// increase the count
          if(count==(tree_list.size()-1)){							// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
            list_size=list_size*2;									// multiply list size by 2
            tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
            tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
          }
        }
      }
    }  
  }
  
  tree_list=resize_bcf(tree_list,count);								// remove the values in tree_list that aren't filled in
  tree_mat_list=resize_bcf(tree_mat_list,count);						// remove the values in tree_mat_list that aren't filled in
  tree_lik.resize(count);											// remove the values in tree_lik that aren't filled in
  tree_parent.resize(count);										// remove the values in tree_parent that aren't filled in
  if(count>0){													// If these lists are nonempty
    //Rcout << "count >0. error after 1532 .\n" ;
    eval_model=evaluate_model_occams_window_bcf(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    NumericVector testlik =eval_model[0];						// testlik is tree_lik after removing models outside Occam's window
    List testtree =eval_model[1];    							// testtree is tree_list after removing models outside Occam's window
    List testmat =eval_model[2]; 								// testmat is tree_mat_list after removing models outside Occam's window
    IntegerVector testpar =eval_model[3];						// testpar is tree_parent after removing models outside Occam's window
    
    if(testlik.size()>0){										// If a nonzero number of models remain in Occam's window
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){							// maxOWsize is an input variable
        IntegerVector owindices=orderforOW__bcf(testlik);			// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on.
        owindices=owindices-1;								// Presumably the match function in orderforOW__bcf gives indices beginning at 1, and therefore 1 must be taken away from all index values.
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);					// create vector temp_olik of size maxOWsize
        List temp_otrees(maxOWsize);						// create List temp_otrees of size maxOWsize
        List temp_omat(maxOWsize);							// create List temp_omat of size maxOWsize
        IntegerVector temp_oparent(maxOWsize);				// create IntegerVector temp_oparent of size maxOWsize
        for(int t=0;t<maxOWsize;t++){						// loop of length maxOWsize
          temp_olik[t]=testlik[owindices[t]];				// temp_olik is maxOWsize largest BIC elements of testlik ordered by BIC (descending?)
          temp_otrees[t]=testtree[owindices[t]];			// temp_otrees is maxOWsize largest BIC elements of testtree ordered by BIC (descending?)
          temp_omat[t]=testmat[owindices[t]];				// temp_omat is maxOWsize largest BIC elements of testmat ordered by BIC (descending?)
          temp_oparent[t]=testpar[owindices[t]];			// temp_oparent is maxOWsize largest BIC elements of testpar ordered by BIC (descending?)
        }
        testlik=temp_olik;			// reset testlik equal to temp_olik
        testtree=temp_otrees;		// reset testtree equal to temp_otrees
        testmat=temp_omat;			// reset testmat equal to temp_omat
        testpar=temp_oparent;		// reset testpar equal to temp_oparent
      }
      ret[0]=lowest_BIC;				// first element of output list. Lowest BIC
      ret[1]=best_sv;					// second element of output list. Best splitting variable
      ret[2]=best_sp;					// third element of ouput list. Best splitting point
      ret[3]=likeliest_tree;			// fourth element of output list. List containing treee table and tree matrix of lowest BIC tree
      ret[4]=testtree;				// fifth element of output list. List of tree tables
      ret[5]=testlik;					// sixth element of output list. List of BICs
      ret[6]=testmat;					// seventh element of output list. List of tree matrices
      ret[7]=testpar;					// eighth element of output list. Vector with all elements equal to the input value of parent?
      ret[8]=no_tree_err;				// ninth element of output list. Boolean equal to false
      
      return (ret);					// return the list ret
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;					// Boolean equal to true.
      List gr(1);						// Create a list, gr, of length one.
      gr[0]=no_tree_err;				// First element of gr is the Boolean no_tree_err.
      return(gr);						// Return the list gr.
    }
  }else{
    no_tree_err=1;						// Boolean equal to true.
    List gr(1);							// Create a list, gr, of length one.
    gr[0]=no_tree_err;					// First element of gr is the Boolean no_tree_err.
    //Rcout << "count =0. error after 1582 .\n" ;
    
    return(gr);							// Return the list gr.
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_tau_round1_bcf(NumericVector resids,arma::mat& x_moderate_a,NumericMatrix tree_table_tau,NumericMatrix tree_mat_tau,
                               double a_mu,double a_tau,double mu_mu,double mu_tau,double nu,double lambda,double c,
                               double lowest_BIC,int parent,NumericMatrix cp_mat,
                               double alpha_mu,double beta_mu,double alpha_tau,double beta_tau,
                               int maxOWsize,//int first_round,
                               List prev_sum_trees_mu,List prev_sum_trees_mat_mu,
                               NumericVector y_scaled,IntegerVector parent2,int i, NumericVector z){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;																// create an integer split_var. Not initialized
  NumericMatrix tree_table_tau_c=tree_table_tau;										// create a matrix treetable_c equal to the input matrix.
  NumericMatrix tree_mat_tau_c=tree_mat_tau;											// create a matric tree_mat_tau_c equal to the input matrix.
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  
  NumericVector terminal_nodes=find_term_nodes_bcf(tree_table_tau_c);					// find term nodes function defined line 168. Gives index of values of treetable_c that are term nodes (indices from 1 to length of vector). Why not integer vector?
  //Rcout << "terminal_nodes[0] equals " << terminal_nodes[0] << ".\n" ;
  //Rcout << "length terminal_nodes equals " << terminal_nodes.size() << ".\n" ;
  
  IntegerVector change_node1;													// create an IntegerVector. Initially all values are 0, but size not given.
  int list_size=1000;															// create a list of size 1000. Why 1000.
  std::vector<double> tree_lik(list_size);									// create a vector of size 1000
  List proposal_tree;															// create an empty list
  List ret(9);																// create a list with 9 elements
  bool no_tree_err=0;															// create a boolean initialized equal to 0 (false).
  List likeliest_tree;														// create an empty list
  List tree_list(list_size);													// create a list of size 1000
  List tree_mat_list(list_size);												// create a list of size 1000
  int count=0;																// create an integer variable, initialized equal to 0.
  std::vector<int> tree_parent(list_size);									// create an integer vector of size 1000
  int best_sv;																// create an integer vector. Not initialized.
  double best_sp;																// create a double variable
  double tree_prior=1;														// create a double variable, initialized equal to 0.
  List changetree;															// create an empty list.
  double BIC;																	// create a double variable. Not initialized.
  int p;																		// create an integer variable. Not initialized.
  int p_other_mu=0;
  //int p_other_tau=0;
  List eval_model;															// create an empty list.
  NumericVector int_nodes;													// create a NumericVector. Initially all values are 0, but size not given.
  NumericVector other_int_nodes_mu;
  NumericVector other_int_nodes_tau;
  arma::colvec curr_col=x_moderate_a.col(0);											// let curr_col be an arma colvec equal to the first column of the inut arma mat x_moderate_a.
  arma::uvec grow_obs=find_term_obs_bcf(tree_mat_tau_c,terminal_nodes[0]);				// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[0] (for leftmost column of tree_mat_tau_c that has elements equal to terminal_nodes[0]).
  //Rcout << "length of grow_obs equals " << grow_obs.n_elem<< ".\n";
  NumericVector d1=unique(find_term_cols_bcf(tree_mat_tau_c,terminal_nodes[0]));		// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[0]. Unique funtion removes duplicated of columns. Unique also sorts descending. Why not IntegerVector
  arma::mat data_curr_node=x_moderate_a.rows(grow_obs);					// matrix consisting of first grow_obs rows of x_moderate_a.
  double d=d1[0];																// index of rightmost column of tree_mat_tau_c with at least one element equal to terminal_nodes[0]
  NumericVector get_min=get_grow_obs_bcf(x_moderate_a,wrap(grow_obs),cp_mat(0,0)+1);		// obtain the elements of the cp_mat(0,0)+1^th column of x_moderate_a that are indexed by grow_obs
  double lik;																	// create a variable called lik. Not initialized.
  
  for(int l=0;l<terminal_nodes.size();l++){										//	vector of length equal to that of terminal_nodes
    //loop over each terminal node
    grow_obs=find_term_obs_bcf(tree_mat_tau_c,terminal_nodes[l]);						// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[l] (letter l) (for leftmost column of tree_mat_tau_c that has elements equal to terminal_nodes[l] (letter l)).
    //depth of tree at current terminal node
    d1=unique(find_term_cols_bcf(tree_mat_tau_c,terminal_nodes[l]));						// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[l] (letter l). Unique funtion removes duplcated of columns. Unique also sorts descending.
    data_curr_node=x_moderate_a.rows(grow_obs);								// matrix consisting of first grow_obs rows of x_moderate_a. Note grow_obs changed on line 926, therefore not duplicating line 919.
    d=d1[0];																	// index of rightmost column of tree_mat_tau_c with at least one element equal to terminal_nodes[l] (letter l)
    int w=cp_mat.nrow();														// w is number of rows of cp_mat
    if(data_curr_node.n_rows<=2){												// if data_curr_node has 2 rows or less.
      //Rcout << "terminal_nodes[l] equals " << terminal_nodes[l]<< ".\n" ;
      //Rcout << "num obs in grow_obs" << grow_obs.n_elem<< ".\n" ;
      //Rcout << " num obs in data_curr_node" << data_curr_node.n_rows<< ".\n" ;
      //Rcout << " iteration number " << l<< ".\n" ;
      throw std::range_error("not enough obs in node to grow any further Line 1360");	// throw an error message. Not enough observations.
      //continue;
    }
    //Rcout << " iteration number " << l<< ".\n" ;
    
    for(int k=0;k<w;k++){														// loop of length w, the number of rows of cp_mat
      p_other_mu=0;
      //Rcout << "inner iteration number " << k<< ".\n" ;
      split_var=cp_mat(k,0)+1;												// split_var is k+1^th row, 1st column of cp_mat +1
      arma::colvec curr_cols=x_moderate_a.col(split_var-1);							// curr_cols is the split_var^th column of x_moderate_a
      get_min=get_grow_obs_bcf(x_moderate_a,wrap(grow_obs),split_var);					// obtain the elements of the split_var^th column of x_moderate_a that are indexed by grow_obs 
      NumericVector z_growvec = get_grow_obs_in_z_bcf(z_ar, wrap(grow_obs));
      
      //Removing unnecessary lines
      //if(get_min.size()<=2){													// If get_min has 2 or less observations. (too few variables to split on?)
      //  throw std::range_error("obs in this terminal node are too small");
      //}
      
      double split_point=cp_mat(k,1);											// variable split_point equals element in k+1^th row 2nd column of cp_mat
      arma::vec curr_cols2=data_curr_node.col(split_var-1);					// curr_cols2 is split_var^th column of data_curr_node
      arma::vec get_min_a=Rcpp::as<arma::vec>(get_min);		// converts to arma vec
      arma::vec z_growvec_a=Rcpp::as<arma::vec>(z_growvec);		// converts to arma vec
      arma::vec get_min_a_treated = get_min_a.elem(arma::find(z_growvec_a==1));	
      
      arma::vec ld_prop=get_min_a_treated.elem(arma::find(get_min_a_treated <= split_point));	// ld_prop is elements of curr_cols2 <= split_point
      arma::vec rd_prop=get_min_a_treated.elem(arma::find(get_min_a_treated > split_point));		// rd_prop is elements of curr_cols2 > split_point
      
      if(ld_prop.size()<=2 || rd_prop.size()<=2){									// if 2 or less observations in either ld_prop or rd_prop
        continue;																// skip to next iteration of the loop
      }
      
      //Rcout << "error after new code .\n" ;
      
      
      proposal_tree=grow_tree_bcf(x_moderate_a,resids,tree_mat_tau_c,terminal_nodes[l],tree_table_tau_c,split_var,split_point,terminal_nodes,wrap(grow_obs),d,get_min,data_curr_node);	// elaborate function on line 469. creates list of 2 elements, tree matrix and tree table. Appears to grow node indexed by terminal_nodes[l] (letter l) (?).
      
      
      //NumericMatrix test =proposal_tree[0];										// first element is tree table
      //NumericMatrix test1 =proposal_tree[1];										// second element is tree matrix
      
      //if(test1.ncol()==3){														// If tree matrix has 3 columns
      //  NumericVector u1=unique(test1(_,0));									// set u1 equal to the unique (ordered descending) elements of 1st column of tree matrix
      //  NumericVector u2=unique(test1(_,1));									// set u2 equal to the unique (ordered descending) elements of 2nd column of tree matrix
      //  NumericVector u3=unique(test1(_,2));									// set u3 equal to the unique (ordered descending) elements of 3rd column of tree matrix
      //}
      
      
      
      
      //if(first_round==1){																// If input value first_round equals 1 (number one)
      
        //lik=likelihood_function_bcf(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);	// set lik equal to tree likelihood defined on line 201.
      
        SEXP s = prev_sum_trees_mu[parent2[i]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
        if(is<List>(s)){												// is s is a list
          List prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          
          List st_tau(1);													// create list, st, of length 2.
          List st_mat_tau(1);												// create lisr, st_mat of length 2.
          st_tau[0]=proposal_tree[0];												// let the first elemetn of st be sum_trees2.
          st_mat_tau[0]=proposal_tree[1];										// let the first element of st_mat be sum_trees_mat2.
          
          //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
          //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
          lik=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2,st_tau,prev_sum_trees_mat_mu2,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          //????????????????? check over this
          // This can all probably be made more efficient by taking account of the fact that there should be just one mu tree
          for(int t=0;t<prev_sum_trees_mu2.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=prev_sum_trees_mu2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=prev_sum_trees_mat_mu2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1412");
            tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not well defined
          NumericMatrix temptestingtabcols = proposal_tree[0];
          if(temptestingtabcols.ncol()<5) throw std::range_error("Line 1417");
          tree_prior*=get_tree_prior_bcf(proposal_tree[0],proposal_tree[1],alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }else{															// if s is not a list
          NumericMatrix prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
          NumericMatrix prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          other_int_nodes_mu = find_term_nodes_bcf(prev_sum_trees_mu2);
          p_other_mu=other_int_nodes_mu.size();
          List st_mu(1);													// create list, st, of length 2.
          List st_mat_mu(1);												// create lisr, st_mat of length 2.
          st_mu[0]=prev_sum_trees_mu2;												// let the first elemetn of st be sum_trees2.
          st_mat_mu[0]=prev_sum_trees_mat_mu2;										// let the first element of st_mat be sum_trees_mat2.
          // return(st);
          List st_tau(1);													// create list, st, of length 2.
          List st_mat_tau(1);												// create lisr, st_mat of length 2.
          st_tau[0]=proposal_tree[0];												// let the first elemetn of st be sum_trees2.
          st_mat_tau[0]=proposal_tree[1];										// let the first element of st_mat be sum_trees_mat2.
          //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
          //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
          lik=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            //ONLY ONE TREE IN mu(x) function
            //NumericMatrix tree=prev_sum_trees_mu2;									// let tree equal (t+1)^th element of st
            //NumericMatrix mat=prev_sum_trees_mat_mu2;								// let mat equal (t+1)^th element of st_mat
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(prev_sum_trees_mu2.ncol()<5) throw std::range_error("Line 1438");
            tree_prior*=get_tree_prior_bcf(prev_sum_trees_mu2,prev_sum_trees_mat_mu2,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not well defined
          NumericMatrix temptestingtabcols = proposal_tree[0];
          if(temptestingtabcols.ncol()<5) throw std::range_error("Line 1443");
          tree_prior*=get_tree_prior_bcf(proposal_tree[0],proposal_tree[1],alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }  
      //}else{
      //  throw std::range_error("get_best_split_tau_round1_bcf should only be used in the first round");	// throw an error.
      //}
      //at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      int_nodes=find_term_nodes_bcf(proposal_tree[0]);				// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector? 
      p=int_nodes.size();											// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
      BIC=-2*(lik+log(tree_prior))+(p_other_mu+p)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      //BIC=-2*(lik)+(p_other_mu+p)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      if(BIC<lowest_BIC){											// if statement for updating lowest BIC...etc.
        lowest_BIC=BIC;											// update (input variable) lowest_BIC
        best_sv=split_var;										// set a value for, or update best_sv
        best_sp=split_point;									// set a value for, or update split_point
        likeliest_tree=proposal_tree;							// set a value for, or update likeliest_tree
        tree_list[count]=proposal_tree[0];						// add an element to the list of tree tables
        tree_mat_list[count]=proposal_tree[1];					// add an element to the list of tree matrices
        tree_lik[count]=BIC;									// add an element to the vector of tree liklihoods (BICs)
        tree_parent[count]=parent;								// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
        count++;												// increase the count
        if(count==(tree_list.size()-1)){						// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
          list_size=list_size*2;									// multiply list size by 2
          tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
          tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){										// If in Occam's window (but not the new minimum BIC)
          if(is<NumericMatrix>(proposal_tree[0])){					// If proposal_tree[0] is a Numeric Matrix, do nothing
            //std::cout<<"its a matrix "<<"\n";						// if this is used, it should say "it's a matrix".
          }else{														// If proposal_tree[0] is NOT a NumericMatrix
            throw std::range_error("proposal tree not a matrix");	// Then throw an error
          }
          tree_list[count]=proposal_tree[0];							// add an element to the list of tree tables
          tree_mat_list[count]=proposal_tree[1];						// add an element to the list of tree matrices
          tree_lik[count]=BIC;										// add an element to the vector of tree liklihoods (BICs)
          tree_parent[count]=parent;									// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
          count++;													// increase the count
          if(count==(tree_list.size()-1)){							// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
            list_size=list_size*2;									// multiply list size by 2
            tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
            tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
          }
        }
      }
    }  
  }
  
  tree_list=resize_bcf(tree_list,count);								// remove the values in tree_list that aren't filled in
  tree_mat_list=resize_bcf(tree_mat_list,count);						// remove the values in tree_mat_list that aren't filled in
  tree_lik.resize(count);											// remove the values in tree_lik that aren't filled in
  tree_parent.resize(count);										// remove the values in tree_parent that aren't filled in
  if(count>0){													// If these lists are nonempty
    //Rcout << "count >0. error after 1532 .\n" ;
    eval_model=evaluate_model_occams_window_bcf(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    NumericVector testlik =eval_model[0];						// testlik is tree_lik after removing models outside Occam's window
    List testtree =eval_model[1];    							// testtree is tree_list after removing models outside Occam's window
    List testmat =eval_model[2]; 								// testmat is tree_mat_list after removing models outside Occam's window
    IntegerVector testpar =eval_model[3];						// testpar is tree_parent after removing models outside Occam's window
    
    if(testlik.size()>0){										// If a nonzero number of models remain in Occam's window
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){							// maxOWsize is an input variable
        IntegerVector owindices=orderforOW__bcf(testlik);			// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on.
        owindices=owindices-1;								// Presumably the match function in orderforOW__bcf gives indices beginning at 1, and therefore 1 must be taken away from all index values.
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);					// create vector temp_olik of size maxOWsize
        List temp_otrees(maxOWsize);						// create List temp_otrees of size maxOWsize
        List temp_omat(maxOWsize);							// create List temp_omat of size maxOWsize
        IntegerVector temp_oparent(maxOWsize);				// create IntegerVector temp_oparent of size maxOWsize
        for(int t=0;t<maxOWsize;t++){						// loop of length maxOWsize
          temp_olik[t]=testlik[owindices[t]];				// temp_olik is maxOWsize largest BIC elements of testlik ordered by BIC (descending?)
          temp_otrees[t]=testtree[owindices[t]];			// temp_otrees is maxOWsize largest BIC elements of testtree ordered by BIC (descending?)
          temp_omat[t]=testmat[owindices[t]];				// temp_omat is maxOWsize largest BIC elements of testmat ordered by BIC (descending?)
          temp_oparent[t]=testpar[owindices[t]];			// temp_oparent is maxOWsize largest BIC elements of testpar ordered by BIC (descending?)
        }
        testlik=temp_olik;			// reset testlik equal to temp_olik
        testtree=temp_otrees;		// reset testtree equal to temp_otrees
        testmat=temp_omat;			// reset testmat equal to temp_omat
        testpar=temp_oparent;		// reset testpar equal to temp_oparent
      }
      ret[0]=lowest_BIC;				// first element of output list. Lowest BIC
      ret[1]=best_sv;					// second element of output list. Best splitting variable
      ret[2]=best_sp;					// third element of ouput list. Best splitting point
      ret[3]=likeliest_tree;			// fourth element of output list. List containing treee table and tree matrix of lowest BIC tree
      ret[4]=testtree;				// fifth element of output list. List of tree tables
      ret[5]=testlik;					// sixth element of output list. List of BICs
      ret[6]=testmat;					// seventh element of output list. List of tree matrices
      ret[7]=testpar;					// eighth element of output list. Vector with all elements equal to the input value of parent?
      ret[8]=no_tree_err;				// ninth element of output list. Boolean equal to false
      
      return (ret);					// return the list ret
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;					// Boolean equal to true.
      List gr(1);						// Create a list, gr, of length one.
      gr[0]=no_tree_err;				// First element of gr is the Boolean no_tree_err.
      return(gr);						// Return the list gr.
    }
  }else{
    no_tree_err=1;						// Boolean equal to true.
    List gr(1);							// Create a list, gr, of length one.
    gr[0]=no_tree_err;					// First element of gr is the Boolean no_tree_err.
    //Rcout << "count =0. error after 1582 .\n" ;
    
    return(gr);							// Return the list gr.
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_sum_tau_bcf(NumericVector resids,arma::mat& x_moderate_a,NumericMatrix tree_table_tau,NumericMatrix tree_mat_tau,
                            double a_mu,double a_tau,double mu_mu,double mu_tau,double nu,double lambda,double c,
                            double lowest_BIC,int parent,NumericMatrix cp_mat,
                            double alpha_mu,double beta_mu,double alpha_tau,double beta_tau,
                            int maxOWsize,//int first_round,
                            List prev_sum_trees_mu,List prev_sum_trees_tau,List prev_sum_trees_mat_mu,List prev_sum_trees_mat_tau,
                            NumericVector y_scaled,IntegerVector parent2,int i,NumericVector z){
  //this function will search through all predictive split points and return those within Occam's Window.
  int split_var;																// create an integer split_var. Not initialized
  NumericMatrix tree_table_tau_c=tree_table_tau;										// create a matrix treetable_c equal to the input matrix.
  NumericMatrix tree_mat_tau_c=tree_mat_tau;											// create a matric tree_mat_tau_c equal to the input matrix.
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  
  NumericVector terminal_nodes=find_term_nodes_bcf(tree_table_tau_c);					// find term nodes function defined line 168. Gives index of values of treetable_c that are term nodes (indices from 1 to length of vector). Why not integer vector?
  IntegerVector change_node1;													// create an IntegerVector. Initially all values are 0, but size not given.
  int list_size=1000;															// create a list of size 1000. Why 1000.
  std::vector<double> tree_lik(list_size);									// create a vector of size 1000
  List proposal_tree;															// create an empty list
  List ret(9);																// create a list with 9 elements
  bool no_tree_err=0;															// create a boolean initialized equal to 0 (false).
  List likeliest_tree;														// create an empty list
  List tree_list(list_size);													// create a list of size 1000
  List tree_mat_list(list_size);												// create a list of size 1000
  int count=0;																// create an integer variable, initialized equal to 0.
  std::vector<int> tree_parent(list_size);									// create an integer vector of size 1000
  int best_sv;																// create an integer vector. Not initialized.
  double best_sp;																// create a double variable
  double tree_prior=1;														// create a double variable, initialized equal to 0.
  List changetree;															// create an empty list.
  double BIC;																	// create a double variable. Not initialized.
  //int p;																		// create an integer variable. Not initialized.
  int p_other_mu=0;
  int p_other_tau=0;
  List eval_model;															// create an empty list.
  NumericVector int_nodes;													// create a NumericVector. Initially all values are 0, but size not given.
  NumericVector other_int_nodes_mu;
  NumericVector other_int_nodes_tau;
  arma::colvec curr_col=x_moderate_a.col(0);											// let curr_col be an arma colvec equal to the first column of the inut arma mat x_moderate_a.
  arma::uvec grow_obs=find_term_obs_bcf(tree_mat_tau_c,terminal_nodes[0]);				// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[0] (for leftmost column of tree_mat_tau_c that has elements equal to terminal_nodes[0]).
  NumericVector d1=unique(find_term_cols_bcf(tree_mat_tau_c,terminal_nodes[0]));		// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[0]. Unique funtion removes duplicated of columns. Unique also sorts descending. Why not IntegerVector
  arma::mat data_curr_node=x_moderate_a.rows(grow_obs);					// matrix consisting of first grow_obs rows of x_moderate_a.
  double d=d1[0];																// index of rightmost column of tree_mat_tau_c with at least one element equal to terminal_nodes[0]
  NumericVector get_min=get_grow_obs_bcf(x_moderate_a,wrap(grow_obs),cp_mat(0,0)+1);		// obtain the elements of the cp_mat(0,0)+1^th column of x_moderate_a that are indexed by grow_obs
  double lik;																	// create a variable called lik. Not initialized.
  
  for(int l=0;l<terminal_nodes.size();l++){										//	vector of length equal to that of terminal_nodes
    //loop over each terminal node
    grow_obs=find_term_obs_bcf(tree_mat_tau_c,terminal_nodes[l]);						// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[l] (letter l) (for leftmost column of tree_mat_tau_c that has elements equal to terminal_nodes[l] (letter l)).
    //depth of tree at current terminal node
    d1=unique(find_term_cols_bcf(tree_mat_tau_c,terminal_nodes[l]));						// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[l] (letter l). Unique funtion removes duplcated of columns. Unique also sorts descending.
    data_curr_node=x_moderate_a.rows(grow_obs);								// matrix consisting of first grow_obs rows of x_moderate_a. Note grow_obs changed on line 926, therefore not duplicating line 919.
    d=d1[0];																	// index of rightmost column of tree_mat_tau_c with at least one element equal to terminal_nodes[l] (letter l)
    int w=cp_mat.nrow();														// w is number of rows of cp_mat
    if(data_curr_node.n_rows<=2){												// if data_curr_node has 2 rows or less.
      throw std::range_error("not enough obs in node to grow any further Line 1602");	// throw an error message. Not enough observations.
      //continue;
    }
    for(int k=0;k<w;k++){														// loop of length w, the number of rows of cp_mat
      p_other_mu=0;
      p_other_tau=0;
      split_var=cp_mat(k,0)+1;												// split_var is k+1^th row, 1st column of cp_mat +1
      arma::colvec curr_cols=x_moderate_a.col(split_var-1);							// curr_cols is the split_var^th column of x_moderate_a
      get_min=get_grow_obs_bcf(x_moderate_a,wrap(grow_obs),split_var);					// obtain the elements of the split_var^th column of x_moderate_a that are indexed by grow_obs 
      NumericVector z_growvec = get_grow_obs_in_z_bcf(z_ar, wrap(grow_obs));
      
      //Removing unnecessary lines
      //if(get_min.size()<=2){													// If get_min has 2 or less observations. (too few variables to split on?)
      //  throw std::range_error("obs in this terminal node are too small");
      //}
      
      double split_point=cp_mat(k,1);											// variable split_point equals element in k+1^th row 2nd column of cp_mat
      arma::vec curr_cols2=data_curr_node.col(split_var-1);					// curr_cols2 is split_var^th column of data_curr_node
      arma::vec get_min_a=Rcpp::as<arma::vec>(get_min);		// converts to arma vec
      arma::vec z_growvec_a=Rcpp::as<arma::vec>(z_growvec);		// converts to arma vec
      arma::vec get_min_a_treated = get_min_a.elem(arma::find(z_growvec_a==1));	
      
      arma::vec ld_prop=get_min_a_treated.elem(arma::find(get_min_a_treated <= split_point));	// ld_prop is elements of curr_cols2 <= split_point
      arma::vec rd_prop=get_min_a_treated.elem(arma::find(get_min_a_treated > split_point));		// rd_prop is elements of curr_cols2 > split_point
      
      //arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));	// ld_prop is elements of curr_cols2 <= split_point
      //arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));		// rd_prop is elements of curr_cols2 > split_point
      
      if(ld_prop.size()<=2 || rd_prop.size()<=2){									// if 2 or less observations in either ld_prop or rd_prop
        continue;																// skip to next iteration of the loop
      }
      proposal_tree=grow_tree_bcf(x_moderate_a,resids,tree_mat_tau_c,terminal_nodes[l],tree_table_tau_c,split_var,split_point,terminal_nodes,wrap(grow_obs),d,get_min,data_curr_node);	// elaborate function on line 469. creates list of 2 elements, tree matrix and tree table. Appears to grow node indexed by terminal_nodes[l] (letter l) (?).
      
      
      //NumericMatrix test =proposal_tree[0];										// first element is tree table
      //NumericMatrix test1 =proposal_tree[1];										// second element is tree matrix
      
      //if(test1.ncol()==3){														// If tree matrix has 3 columns
      //  NumericVector u1=unique(test1(_,0));									// set u1 equal to the unique (ordered descending) elements of 1st column of tree matrix
      //  NumericVector u2=unique(test1(_,1));									// set u2 equal to the unique (ordered descending) elements of 2nd column of tree matrix
      //  NumericVector u3=unique(test1(_,2));									// set u3 equal to the unique (ordered descending) elements of 3rd column of tree matrix
      //}
      
      
      
      // if(first_round==1){	
      //   throw std::range_error("get_best_split_sum_tau_bcf should not be used in the first round");	// throw an error.
      //   // If input value first_round equals 1 (number one)
      //   //lik=likelihood_function_bcf(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);	// set lik equal to tree likelihood defined on line 201.
      //   
      // }else{
        SEXP s_mu = prev_sum_trees_mu[parent2[i]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
        SEXP s_tau = prev_sum_trees_tau[parent2[i]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
        if(is<List>(s_mu)){												// is s is a list
          if(is<List>(s_tau)){												// is s is a list
            //Rcout << "\n RELEVANT LINE 1703.\n";
            
            List prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            //Rcout << "Length of prev_sum_trees_tau = " << prev_sum_trees_tau.size() <<".\n";
            
            //Rcout << "Length of mu mat list = " << prev_sum_trees_mat_mu2.size() <<".\n";
            //Rcout << "Length of tau mat list = " << sum_trees_mat_tau2.size() <<".\n";
            
            //NumericMatrix testexamplemat5 = prev_sum_trees_mat_mu2[0];
            //Rcout << "number of rows of mat into likelihood = " << testexamplemat5.nrow() <<".\n";
            
            sum_trees_tau2.push_back(proposal_tree[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            sum_trees_mat_tau2.push_back(proposal_tree[1]);	
            //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2,sum_trees_tau2,prev_sum_trees_mat_mu2,sum_trees_mat_tau2,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            //????????????????? check over this
            // This can all probably be made more efficient by taking account of the fact that there should be just one mu tree
            for(int t=0;t<prev_sum_trees_mu2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=prev_sum_trees_mu2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=prev_sum_trees_mat_mu2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1660");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<sum_trees_tau2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=sum_trees_tau2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=sum_trees_mat_tau2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1667");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            
          }else{
            List prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            NumericMatrix sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            NumericMatrix sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_tau(2);													// create list, st, of length 2.
            List st_mat_tau(2);												// create lisr, st_mat of length 2.
            st_tau[0]=sum_trees_tau2;												// let the first elemetn of st be sum_trees2.
            st_tau[1]=proposal_tree[0];										// let the second element of st be proposal_tree[0] (tree table).
            st_mat_tau[0]=sum_trees_mat_tau2;										// let the first element of st_mat be sum_trees_mat2.
            st_mat_tau[1]=proposal_tree[1];		
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2,st_tau,prev_sum_trees_mat_mu2,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            //????????????????? check over this
            // This can all probably be made more efficient by taking account of the fact that there should be just one mu tree
            for(int t=0;t<prev_sum_trees_mu2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=prev_sum_trees_mu2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=prev_sum_trees_mat_mu2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1689");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1695");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
          }
          
        }else{															// if s is not a list
          if(is<List>(s_tau)){												// is s is a list
            
            NumericMatrix prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
            NumericMatrix prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_mu(1);													// create list, st, of length 2.
            List st_mat_mu(1);												// create lisr, st_mat of length 2.
            st_mu[0]=prev_sum_trees_mu2;												// let the first elemetn of st be sum_trees2.
            st_mat_mu[0]=prev_sum_trees_mat_mu2;										// let the first element of st_mat be sum_trees_mat2.
            // return(st);
            List sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            
            sum_trees_tau2.push_back(proposal_tree[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            sum_trees_mat_tau2.push_back(proposal_tree[1]);	
            //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,sum_trees_tau2,st_mat_mu,sum_trees_mat_tau2,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
              NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1723");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<sum_trees_tau2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=sum_trees_tau2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=sum_trees_mat_tau2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1730");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
          }else{
            NumericMatrix prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
            NumericMatrix prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_mu(1);													// create list, st, of length 2.
            List st_mat_mu(1);												// create lisr, st_mat of length 2.
            st_mu[0]=prev_sum_trees_mu2;												// let the first elemetn of st be sum_trees2.
            st_mat_mu[0]=prev_sum_trees_mat_mu2;										// let the first element of st_mat be sum_trees_mat2.
            NumericMatrix sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            NumericMatrix sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_tau(2);													// create list, st, of length 2.
            List st_mat_tau(2);												// create lisr, st_mat of length 2.
            st_tau[0]=sum_trees_tau2;												// let the first elemetn of st be sum_trees2.
            st_tau[1]=proposal_tree[0];										// let the second element of st be proposal_tree[0] (tree table).
            st_mat_tau[0]=sum_trees_mat_tau2;										// let the first element of st_mat be sum_trees_mat2.
            st_mat_tau[1]=proposal_tree[1];		
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
              NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1753");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1760");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
          }
        }  
      //}
      //at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      //int_nodes=find_term_nodes_bcf(proposal_tree[0]);				// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector? 
      //p=int_nodes.size();											// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
      BIC=-2*(lik+log(tree_prior))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      
      //BIC=-2*(lik)+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      // Rcout << "p_other_mu+p_other_tau=" << p_other_mu+p_other_tau << ".\n";
      // Rcout << "lik=" << lik << ".\n";
      // Rcout << "BIC=" << BIC << ".\n";
      // Rcout << "lowest_BIC=" << lowest_BIC << ".\n";
      // Rcout << "LINE 1885 tree_prior=" << tree_prior << ".\n";
      
      if(BIC<lowest_BIC){											// if statement for updating lowest BIC...etc.
        lowest_BIC=BIC;											// update (input variable) lowest_BIC
        best_sv=split_var;										// set a value for, or update best_sv
        best_sp=split_point;									// set a value for, or update split_point
        likeliest_tree=proposal_tree;							// set a value for, or update likeliest_tree
        tree_list[count]=proposal_tree[0];						// add an element to the list of tree tables
        tree_mat_list[count]=proposal_tree[1];					// add an element to the list of tree matrices
        tree_lik[count]=BIC;									// add an element to the vector of tree liklihoods (BICs)
        tree_parent[count]=parent;								// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
        count++;												// increase the count
        if(count==(tree_list.size()-1)){						// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
          list_size=list_size*2;									// multiply list size by 2
          tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
          tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){										// If in Occam's window (but not the new minimum BIC)
          if(is<NumericMatrix>(proposal_tree[0])){					// If proposal_tree[0] is a Numeric Matrix, do nothing
            //std::cout<<"its a matrix "<<"\n";						// if this is used, it should say "it's a matrix".
          }else{														// If proposal_tree[0] is NOT a NumericMatrix
            throw std::range_error("proposal tree not a matrix");	// Then throw an error
          }
          tree_list[count]=proposal_tree[0];							// add an element to the list of tree tables
          tree_mat_list[count]=proposal_tree[1];						// add an element to the list of tree matrices
          tree_lik[count]=BIC;										// add an element to the vector of tree liklihoods (BICs)
          tree_parent[count]=parent;									// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
          count++;													// increase the count
          if(count==(tree_list.size()-1)){							// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
            list_size=list_size*2;									// multiply list size by 2
            tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
            tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
          }
        }
      }
    }  
  }
  tree_list=resize_bcf(tree_list,count);								// remove the values in tree_list that aren't filled in
  tree_mat_list=resize_bcf(tree_mat_list,count);						// remove the values in tree_mat_list that aren't filled in
  tree_lik.resize(count);											// remove the values in tree_lik that aren't filled in
  tree_parent.resize(count);										// remove the values in tree_parent that aren't filled in
  if(count>0){													// If these lists are nonempty
    eval_model=evaluate_model_occams_window_bcf(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    NumericVector testlik =eval_model[0];						// testlik is tree_lik after removing models outside Occam's window
    List testtree =eval_model[1];    							// testtree is tree_list after removing models outside Occam's window
    List testmat =eval_model[2]; 								// testmat is tree_mat_list after removing models outside Occam's window
    IntegerVector testpar =eval_model[3];						// testpar is tree_parent after removing models outside Occam's window
    
    if(testlik.size()>0){										// If a nonzero number of models remain in Occam's window
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){							// maxOWsize is an input variable
        IntegerVector owindices=orderforOW__bcf(testlik);			// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on.
        owindices=owindices-1;								// Presumably the match function in orderforOW__bcf gives indices beginning at 1, and therefore 1 must be taken away from all index values.
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);					// create vector temp_olik of size maxOWsize
        List temp_otrees(maxOWsize);						// create List temp_otrees of size maxOWsize
        List temp_omat(maxOWsize);							// create List temp_omat of size maxOWsize
        IntegerVector temp_oparent(maxOWsize);				// create IntegerVector temp_oparent of size maxOWsize
        for(int t=0;t<maxOWsize;t++){						// loop of length maxOWsize
          temp_olik[t]=testlik[owindices[t]];				// temp_olik is maxOWsize largest BIC elements of testlik ordered by BIC (descending?)
          temp_otrees[t]=testtree[owindices[t]];			// temp_otrees is maxOWsize largest BIC elements of testtree ordered by BIC (descending?)
          temp_omat[t]=testmat[owindices[t]];				// temp_omat is maxOWsize largest BIC elements of testmat ordered by BIC (descending?)
          temp_oparent[t]=testpar[owindices[t]];			// temp_oparent is maxOWsize largest BIC elements of testpar ordered by BIC (descending?)
        }
        testlik=temp_olik;			// reset testlik equal to temp_olik
        testtree=temp_otrees;		// reset testtree equal to temp_otrees
        testmat=temp_omat;			// reset testmat equal to temp_omat
        testpar=temp_oparent;		// reset testpar equal to temp_oparent
      }
      ret[0]=lowest_BIC;				// first element of output list. Lowest BIC
      ret[1]=best_sv;					// second element of output list. Best splitting variable
      ret[2]=best_sp;					// third element of ouput list. Best splitting point
      ret[3]=likeliest_tree;			// fourth element of output list. List containing treee table and tree matrix of lowest BIC tree
      ret[4]=testtree;				// fifth element of output list. List of tree tables
      ret[5]=testlik;					// sixth element of output list. List of BICs
      ret[6]=testmat;					// seventh element of output list. List of tree matrices
      ret[7]=testpar;					// eighth element of output list. Vector with all elements equal to the input value of parent?
      ret[8]=no_tree_err;				// ninth element of output list. Boolean equal to false
      
      return (ret);					// return the list ret
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;					// Boolean equal to true.
      List gr(1);						// Create a list, gr, of length one.
      gr[0]=no_tree_err;				// First element of gr is the Boolean no_tree_err.
      return(gr);						// Return the list gr.
    }
  }else{
    no_tree_err=1;						// Boolean equal to true.
    List gr(1);							// Create a list, gr, of length one.
    gr[0]=no_tree_err;					// First element of gr is the Boolean no_tree_err.
    return(gr);							// Return the list gr.
  }
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List get_best_split_sum_mu_bcf(NumericVector resids,arma::mat& x_control_a,NumericMatrix tree_table_mu,NumericMatrix tree_mat_mu,
                           double a_mu,double a_tau,double mu_mu,double mu_tau,double nu,double lambda,double c,
                           double lowest_BIC,int parent,NumericMatrix cp_mat,
                           double alpha_mu,double beta_mu,double alpha_tau,double beta_tau,
                           int maxOWsize,//int first_round,
                           List prev_sum_trees_mu,List prev_sum_trees_tau,List prev_sum_trees_mat_mu,List prev_sum_trees_mat_tau,
                           NumericVector y_scaled,IntegerVector parent2,int i,NumericVector z){
  //this function will search through all predictive split points and return those within Occam's Window.
  List example_tree_tab = prev_sum_trees_mu[parent2[i]];
  List example_tree_mat = prev_sum_trees_mu[parent2[i]];
  
  // Rcout << "Length of input tree table list used within best split sum = " << example_tree_tab.size() << ".\n";
  // Rcout << "Length of input tree mat list used within best split sum= " << example_tree_mat.size() << ".\n";
  int split_var;																// create an integer split_var. Not initialized
  NumericMatrix tree_table_mu_c=tree_table_mu;										// create a matrix treetable_c equal to the input matrix.
  NumericMatrix tree_mat_mu_c=tree_mat_mu;											// create a matric tree_mat_mu_c equal to the input matrix.
  
  NumericVector terminal_nodes=find_term_nodes_bcf(tree_table_mu_c);					// find term nodes function defined line 168. Gives index of values of treetable_c that are term nodes (indices from 1 to length of vector). Why not integer vector?
  IntegerVector change_node1;													// create an IntegerVector. Initially all values are 0, but size not given.
  int list_size=1000;															// create a list of size 1000. Why 1000.
  std::vector<double> tree_lik(list_size);									// create a vector of size 1000
  List proposal_tree;															// create an empty list
  List ret(9);																// create a list with 9 elements
  bool no_tree_err=0;															// create a boolean initialized equal to 0 (false).
  List likeliest_tree;														// create an empty list
  List tree_list(list_size);													// create a list of size 1000
  List tree_mat_list(list_size);												// create a list of size 1000
  int count=0;																// create an integer variable, initialized equal to 0.
  std::vector<int> tree_parent(list_size);									// create an integer vector of size 1000
  int best_sv;																// create an integer vector. Not initialized.
  double best_sp;																// create a double variable
  double tree_prior=1;														// create a double avriable, initialized equal to 0.
  List changetree;															// create an empty list.
  double BIC;																	// create a double variable. Not initialized.
  //int p;																		// create an integer variable. Not initialized.
  int p_other_mu=0;
  int p_other_tau=0;
  List eval_model;															// create an empty list.
  NumericVector int_nodes;													// create a NumericVector. Initially all values are 0, but size not given.
  NumericVector other_int_nodes_mu;
  NumericVector other_int_nodes_tau;
  arma::colvec curr_col=x_control_a.col(0);											// let curr_col be an arma colvec equal to the first column of the inut arma mat x_control_a.
  arma::uvec grow_obs=find_term_obs_bcf(tree_mat_mu_c,terminal_nodes[0]);				// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[0] (for leftmost column of tree_mat_mu_c that has elements equal to terminal_nodes[0]).
  NumericVector d1=unique(find_term_cols_bcf(tree_mat_mu_c,terminal_nodes[0]));		// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[0]. Unique funtion removes duplicated of columns. Unique also sorts descending. Why not IntegerVector
  arma::mat data_curr_node=x_control_a.rows(grow_obs);					// matrix consisting of first grow_obs rows of x_control_a.
  double d=d1[0];																// index of rightmost column of tree_mat_mu_c with at least one element equal to terminal_nodes[0]
  NumericVector get_min=get_grow_obs_bcf(x_control_a,wrap(grow_obs),cp_mat(0,0)+1);		// obtain the elements of the cp_mat(0,0)+1^th column of x_control_a that are indexed by grow_obs
  double lik;																	// create a variable called lik. Not initialized.
  // Rcout << "get to line 2383.\n";
  for(int l=0;l<terminal_nodes.size();l++){										//	vector of length equal to that of terminal_nodes
    // Rcout << "get to line 2383. l= " << l << ". \n";
    
    //loop over each terminal node
    grow_obs=find_term_obs_bcf(tree_mat_mu_c,terminal_nodes[l]);						// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[l] (letter l) (for leftmost column of tree_mat_mu_c that has elements equal to terminal_nodes[l] (letter l)).
    //depth of tree at current terminal node
    d1=unique(find_term_cols_bcf(tree_mat_mu_c,terminal_nodes[l]));						// d1 is a vector of indexes of (starting at 0) all the columns with at least some elements equal to terminal_nodes[l] (letter l). Unique funtion removes duplcated of columns. Unique also sorts descending.
    data_curr_node=x_control_a.rows(grow_obs);								// matrix consisting of first grow_obs rows of x_control_a. Note grow_obs changed on line 926, therefore not duplicating line 919.
    d=d1[0];																	// index of rightmost column of tree_mat_mu_c with at least one element equal to terminal_nodes[l] (letter l)
    int w=cp_mat.nrow();														// w is number of rows of cp_mat
    if(data_curr_node.n_rows<=2){												// if data_curr_node has 2 rows or less.
      throw std::range_error("not enough obs in node to grow any further Line 1919");	// throw an error message. Not enough observations.
      //continue;
    }
    // Rcout << "get to line 2398. l= " << l << ". \n";
    
    for(int k=0;k<w;k++){														// loop of length w, the number of rows of cp_mat
      // Rcout << "get to line 2401. l= " << l << ". k= " << k << ". \n";
      
      p_other_mu=0;
      p_other_tau=0;
      split_var=cp_mat(k,0)+1;												// split_var is k+1^th row, 1st column of cp_mat +1
      arma::colvec curr_cols=x_control_a.col(split_var-1);							// curr_cols is the split_var^th column of x_control_a
      get_min=get_grow_obs_bcf(x_control_a,wrap(grow_obs),split_var);					// obtain the elements of the split_var^th column of x_control_a that are indexed by grow_obs 
      
      //Removing unnecessary lines
      //if(get_min.size()<=2){													// If get_min has 2 or less observations. (too few variables to split on?)
      //  throw std::range_error("obs in this terminal node are too small");
      //}
      
      double split_point=cp_mat(k,1);											// variable split_point equals element in k+1^th row 2nd column of cp_mat
      arma::vec curr_cols2=data_curr_node.col(split_var-1);					// curr_cols2 is split_var^th column of data_curr_node
      
      arma::vec ld_prop=curr_cols2.elem(arma::find(curr_cols2 <= split_point));	// ld_prop is elements of curr_cols2 <= split_point
      arma::vec rd_prop=curr_cols2.elem(arma::find(curr_cols2> split_point));		// rd_prop is elements of curr_cols2 > split_point
      
      if(ld_prop.size()<=2 || rd_prop.size()<=2){									// if 2 or less observations in either ld_prop or rd_prop
        continue;																// skip to next iteration of the loop
      }
      proposal_tree=grow_tree_bcf(x_control_a,resids,tree_mat_mu_c,terminal_nodes[l],tree_table_mu_c,split_var,split_point,terminal_nodes,wrap(grow_obs),d,get_min,data_curr_node);	// elaborate function on line 469. creates list of 2 elements, tree matrix and tree table. Appears to grow node indexed by terminal_nodes[l] (letter l) (?).
      
      
      
      //NumericMatrix test =proposal_tree[0];										// first element is tree table
      //NumericMatrix test1 =proposal_tree[1];										// second element is tree matrix
      
      //if(test1.ncol()==3){														// If tree matrix has 3 columns
      //  NumericVector u1=unique(test1(_,0));									// set u1 equal to the unique (ordered descending) elements of 1st column of tree matrix
      //  NumericVector u2=unique(test1(_,1));									// set u2 equal to the unique (ordered descending) elements of 2nd column of tree matrix
      //  NumericVector u3=unique(test1(_,2));									// set u3 equal to the unique (ordered descending) elements of 3rd column of tree matrix
      //}
      
      // Rcout << "get to line 2431. l= " << l << ". k= " << k << ". \n";
      
      //if(first_round==1){	
      //  throw std::range_error("get_best_split_sum_tau_bcf should not be used in the first round");	// throw an error.
        // If input value first_round equals 1 (number one)
        //lik=likelihood_function_bcf(resids,proposal_tree[0],proposal_tree[1],a,mu,nu,lambda);	// set lik equal to tree likelihood defined on line 201.
        
      //}else{
        
        SEXP s_mu = prev_sum_trees_mu[parent2[i]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
        SEXP s_tau = prev_sum_trees_tau[parent2[i]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
        
        if(is<List>(s_mu)){												// is s is a list
          if(is<List>(s_tau)){												// is s is a list
            //Rcout << "RELEVANT LINE 2027.\n";
            // Rcout << "get to line 2446. l= " << l << ". k= " << k << ". \n";
            
            List prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            
            prev_sum_trees_mu2.push_back(proposal_tree[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            prev_sum_trees_mat_mu2.push_back(proposal_tree[1]);	
            //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2,sum_trees_tau2,prev_sum_trees_mat_mu2,sum_trees_mat_tau2,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            //????????????????? check over this
            // This can all probably be made more efficient by taking account of the fact that there should be just one mu tree
            // Rcout << "get to line 2460. l= " << l << ". k= " << k << ". \n";
            
            for(int t=0;t<prev_sum_trees_mu2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=prev_sum_trees_mu2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=prev_sum_trees_mat_mu2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1977");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<sum_trees_tau2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=sum_trees_tau2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=sum_trees_mat_tau2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 1984");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            // Rcout << "get to line 2480. l= " << l << ". k= " << k << ". \n";
            
          }else{
            //Rcout << "\n RELEVANT LINE 2057.\n";
            
            List prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            prev_sum_trees_mu2.push_back(proposal_tree[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            prev_sum_trees_mat_mu2.push_back(proposal_tree[1]);	
            NumericMatrix sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            NumericMatrix sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_tau(1);													// create list, st, of length 2.
            List st_mat_tau(1);												// create lisr, st_mat of length 2.
            st_tau[0]=sum_trees_tau2;												// let the first elemetn of st be sum_trees2.
            st_mat_tau[0]=sum_trees_mat_tau2;										// let the first element of st_mat be sum_trees_mat2.
            //Rcout << "After Line 2047.\n";
            //Rcout << "Size of tree table list before sumlik = " << prev_sum_trees_mu2.size() <<" \n";
            //Rcout << "Size of tree mat list before sumlik = " << prev_sum_trees_mat_mu2.size() <<" \n";
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2,st_tau,prev_sum_trees_mat_mu2,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            
            //????????????????? check over this
            // This can all probably be made more efficient by taking account of the fact that there should be just one mu tree
            for(int t=0;t<prev_sum_trees_mu2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=prev_sum_trees_mu2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=prev_sum_trees_mat_mu2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 2006");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 2013");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
          }
          
        }else{															// if s is not a list
          if(is<List>(s_tau)){												// is s is a list
            //Rcout << "RELEVANT LINE 2093.\n";
            
            NumericMatrix prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
            NumericMatrix prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_mu(2);													// create list, st, of length 2.
            List st_mat_mu(2);												// create lisr, st_mat of length 2.
            st_mu[0]=prev_sum_trees_mu2;												// let the first elemetn of st be sum_trees2.
            st_mu[1]=proposal_tree[0];												// let the first elemetn of st be sum_trees2.
            st_mat_mu[0]=prev_sum_trees_mat_mu2;										// let the first element of st_mat be sum_trees_mat2.
            st_mat_mu[1]=proposal_tree[1];										// let the first element of st_mat be sum_trees_mat2.
            // return(st);
            List sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            List sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            
            //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
            //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,sum_trees_tau2,st_mat_mu,sum_trees_mat_tau2,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
              NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 2040");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<sum_trees_tau2.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=sum_trees_tau2[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=sum_trees_mat_tau2[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 2047");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
          }else{
            //Rcout << "RELEVANT LINE 2124.\n";
            
            NumericMatrix prev_sum_trees_mu2=prev_sum_trees_mu[parent2[i]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
            NumericMatrix prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent2[i]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_mu(2);													// create list, st, of length 2.
            List st_mat_mu(2);												// create lisr, st_mat of length 2.
            st_mu[0]=prev_sum_trees_mu2;												// let the first elemetn of st be sum_trees2.
            st_mu[1]=proposal_tree[0];										// let the second element of st be proposal_tree[0] (tree table).
            st_mat_mu[0]=prev_sum_trees_mat_mu2;										// let the first element of st_mat be sum_trees_mat2.
            st_mat_mu[1]=proposal_tree[1];		
            NumericMatrix sum_trees_tau2=prev_sum_trees_tau[parent2[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
            NumericMatrix sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent2[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
            List st_tau(1);													// create list, st, of length 2.
            List st_mat_tau(1);												// create lisr, st_mat of length 2.
            st_tau[0]=sum_trees_tau2;												// let the first elemetn of st be sum_trees2.
            st_mat_tau[0]=sum_trees_mat_tau2;										// let the first element of st_mat be sum_trees_mat2.
            lik=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
            for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
              NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
              other_int_nodes_mu = find_term_nodes_bcf(tree);
              p_other_mu+=other_int_nodes_mu.size();
              NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 2070");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
            for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
              NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
              other_int_nodes_tau = find_term_nodes_bcf(tree);
              p_other_tau+=other_int_nodes_tau.size();
              NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
              //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
              if(tree.ncol()<5) throw std::range_error("Line 2077");
              tree_prior*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            }
          }
        }  
      //} function not used in first round
      
      //at the moment tree prior is only for current tree need to get it for entire sum of tree list.
      
      //int_nodes=find_term_nodes_bcf(proposal_tree[0]);				// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector? 
      //p=int_nodes.size();											// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
      BIC=-2*(lik+log(tree_prior))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      //BIC=-2*(lik)+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
      // Rcout << "p_other_mu+p_other_tau=" << p_other_mu+p_other_tau << ".\n";
      // Rcout << "lik=" << lik << ".\n";
      // Rcout << "BIC=" << BIC << ".\n";
      // Rcout << "lowest_BIC=" << lowest_BIC << ".\n";
      // Rcout << "LINE 2247 tree_prior=" << tree_prior << ".\n";
      
      if(BIC<lowest_BIC){											// if statement for updating lowest BIC...etc.
        lowest_BIC=BIC;											// update (input variable) lowest_BIC
        best_sv=split_var;										// set a value for, or update best_sv
        best_sp=split_point;									// set a value for, or update split_point
        likeliest_tree=proposal_tree;							// set a value for, or update likeliest_tree
        tree_list[count]=proposal_tree[0];						// add an element to the list of tree tables
        tree_mat_list[count]=proposal_tree[1];					// add an element to the list of tree matrices
        tree_lik[count]=BIC;									// add an element to the vector of tree liklihoods (BICs)
        tree_parent[count]=parent;								// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
        count++;												// increase the count
        if(count==(tree_list.size()-1)){						// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
          list_size=list_size*2;									// multiply list size by 2
          tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
          tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
          tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
        }
      }else{
        if((BIC)-(lowest_BIC)<=c){										// If in Occam's window (but not the new minimum BIC)
          if(is<NumericMatrix>(proposal_tree[0])){					// If proposal_tree[0] is a Numeric Matrix, do nothing
            //std::cout<<"its a matrix "<<"\n";						// if this is used, it should say "it's a matrix".
          }else{														// If proposal_tree[0] is NOT a NumericMatrix
            throw std::range_error("proposal tree not a matrix");	// Then throw an error
          }
          tree_list[count]=proposal_tree[0];							// add an element to the list of tree tables
          tree_mat_list[count]=proposal_tree[1];						// add an element to the list of tree matrices
          tree_lik[count]=BIC;										// add an element to the vector of tree liklihoods (BICs)
          tree_parent[count]=parent;									// add an element to the vector tree_parent. NOTE: All elements will be equal to the input value of parent which is not changed in this function.
          count++;													// increase the count
          if(count==(tree_list.size()-1)){							// Increase list size if the length of the various lists and vectors were not set to be long enough (i.e. 1000 was too small)
            list_size=list_size*2;									// multiply list size by 2
            tree_list=resize_bigger_bcf(tree_list,list_size);			// increase tree_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_mat_list=resize_bigger_bcf(tree_mat_list,list_size);	// increase tree_mat_list to twice its previous size. resize_bigger_bcf defined on line 776
            tree_lik.resize(list_size);								// increase tree_lik to twice its previous size
            tree_parent.resize(list_size);							// increase tree_parent to twice its previous size
          }
        }
      }
    }  
  }
  // Rcout << "get to line 2638.\n";
  
  tree_list=resize_bcf(tree_list,count);								// remove the values in tree_list that aren't filled in
  tree_mat_list=resize_bcf(tree_mat_list,count);						// remove the values in tree_mat_list that aren't filled in
  tree_lik.resize(count);											// remove the values in tree_lik that aren't filled in
  tree_parent.resize(count);										// remove the values in tree_parent that aren't filled in
  if(count>0){													// If these lists are nonempty
    eval_model=evaluate_model_occams_window_bcf(wrap(tree_lik),lowest_BIC,log(c),wrap(tree_list),wrap(tree_mat_list),wrap(tree_parent));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    NumericVector testlik =eval_model[0];						// testlik is tree_lik after removing models outside Occam's window
    List testtree =eval_model[1];    							// testtree is tree_list after removing models outside Occam's window
    List testmat =eval_model[2]; 								// testmat is tree_mat_list after removing models outside Occam's window
    IntegerVector testpar =eval_model[3];						// testpar is tree_parent after removing models outside Occam's window
    
    if(testlik.size()>0){										// If a nonzero number of models remain in Occam's window
      //check if number of trees to be returned is greater than maxOWsize if so only return the best maxOWsize models
      if(testlik.size()>maxOWsize){							// maxOWsize is an input variable
        IntegerVector owindices=orderforOW__bcf(testlik);			// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on.
        owindices=owindices-1;								// Presumably the match function in orderforOW__bcf gives indices beginning at 1, and therefore 1 must be taken away from all index values.
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);					// create vector temp_olik of size maxOWsize
        List temp_otrees(maxOWsize);						// create List temp_otrees of size maxOWsize
        List temp_omat(maxOWsize);							// create List temp_omat of size maxOWsize
        IntegerVector temp_oparent(maxOWsize);				// create IntegerVector temp_oparent of size maxOWsize
        for(int t=0;t<maxOWsize;t++){						// loop of length maxOWsize
          temp_olik[t]=testlik[owindices[t]];				// temp_olik is maxOWsize largest BIC elements of testlik ordered by BIC (descending?)
          temp_otrees[t]=testtree[owindices[t]];			// temp_otrees is maxOWsize largest BIC elements of testtree ordered by BIC (descending?)
          temp_omat[t]=testmat[owindices[t]];				// temp_omat is maxOWsize largest BIC elements of testmat ordered by BIC (descending?)
          temp_oparent[t]=testpar[owindices[t]];			// temp_oparent is maxOWsize largest BIC elements of testpar ordered by BIC (descending?)
        }
        testlik=temp_olik;			// reset testlik equal to temp_olik
        testtree=temp_otrees;		// reset testtree equal to temp_otrees
        testmat=temp_omat;			// reset testmat equal to temp_omat
        testpar=temp_oparent;		// reset testpar equal to temp_oparent
      }
      ret[0]=lowest_BIC;				// first element of output list. Lowest BIC
      ret[1]=best_sv;					// second element of output list. Best splitting variable
      ret[2]=best_sp;					// third element of ouput list. Best splitting point
      ret[3]=likeliest_tree;			// fourth element of output list. List containing treee table and tree matrix of lowest BIC tree
      ret[4]=testtree;				// fifth element of output list. List of tree tables
      ret[5]=testlik;					// sixth element of output list. List of BICs
      ret[6]=testmat;					// seventh element of output list. List of tree matrices
      ret[7]=testpar;					// eighth element of output list. Vector with all elements equal to the input value of parent?
      ret[8]=no_tree_err;				// ninth element of output list. Boolean equal to false
      
      return (ret);					// return the list ret
    }else{
      //if no trees are found within Occam's window function will return an error to main
      no_tree_err=1;					// Boolean equal to true.
      List gr(1);						// Create a list, gr, of length one.
      gr[0]=no_tree_err;				// First element of gr is the Boolean no_tree_err.
      return(gr);						// Return the list gr.
    }
  }else{
    no_tree_err=1;						// Boolean equal to true.
    List gr(1);							// Create a list, gr, of length one.
    gr[0]=no_tree_err;					// First element of gr is the Boolean no_tree_err.
    return(gr);							// Return the list gr.
  }
}

//######################################################################################################################//

// [[Rcpp::export]]
NumericVector update_mean_var_bcf(NumericMatrix tree_table,NumericMatrix tree_matrix,NumericVector resids,double a){
  List update_params(1);															// create a list of length 1
  NumericVector terminal_nodes;													// create a Numeric vector. All values initially 0.
  arma::uvec term_obs;															// create an unsigned integer vector. Not initialized or given a length.
  terminal_nodes= find_term_nodes_bcf(tree_table);									// find term nodes function defined line 168. Gives index of values of treetable_temp that are term nodes (indices from 1 to length of vector). Why not integer vector?
  NumericVector Tj(terminal_nodes.size());										// Tj is a vector of length equal to length of terminal_nodes
  NumericVector new_mean(terminal_nodes.size());									// new_mean is a vector of length equal to length of terminal_nodes
  arma::vec armaresids=as<arma::vec>(resids);										// copy input vector resids to an arma vec
  
  for(int k=0;k< terminal_nodes.size();k++){    									// for-loop of length equal to that of terminal_nodes.
    term_obs=find_term_obs_bcf(tree_matrix,terminal_nodes[k]);						// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[k] (for leftmost column of tree_matrix that has elements equal to terminal_nodes[k]). 
    //get the number of observations in node k
    Tj[k]=term_obs.size();														// Tj[k] is length of term_obs. Number of obs in k^th termonal node? Tj should be an IntegerVEctor
    NumericVector  get_mean(term_obs.size());									// create vector get_mean of length equal to number of termnal observations
    for(int i=0;i<Tj[k];i++){													// for-loop of length equal to number of terminal observations
      get_mean[i]=resids[term_obs[i]];										// let get_mean[i] be the element of input vector resids indexed by term_obs[i].
    } 
    double sum_resids=std::accumulate(get_mean.begin(),get_mean.end(),0.0);		// sum_resides is sum of values in get_mean
    new_mean[k]=sum_resids/(Tj[k]+a);											// New mean for terminal node k. new_mean[k] is sum of residuals divided by (number of resids +a)
    arma::uvec temp;															// create an uninitialized unsigned integer vector. Length not given.
    term_obs=temp;																// reset term_obs equal to temp. (for next iteration)
  }  
  
  return(new_mean);		// return vector of new means for each terminal node.
}
//######################################################################################################################//

// [[Rcpp::export]]
List update_predictions_bcf(NumericMatrix tree_table,NumericMatrix tree_matrix,NumericVector new_mean,int n){
  
  List updated_preds(2);												// create list of length 2.
  NumericVector new_preds(n);											// create vector of length n.
  NumericVector terminal_nodes;										// create vector. No initial length.
  arma::uvec term_obs;												// create unsigned integer vector. No initial length.
  terminal_nodes=find_term_nodes_bcf(tree_table);							// find term nodes function defined line 168. Gives index of values of treetable_temp that are term nodes (indices from 1 to length of vector). Why not integer vector?
  
  for(int k=0;k<terminal_nodes.size();k++){							// for-loop of length equal to number of terminal nodes
    term_obs=find_term_obs_bcf(tree_matrix,terminal_nodes[k]);        	// function find_term_obs_bcf. Gives indices of elements equal to terminal_nodes[k] (for leftmost column of tree_matrix that has elements equal to terminal_nodes[k]). 
    //update the terminal node mean of the selected tree nodes:
    tree_table(terminal_nodes[k]-1,5)= new_mean[k];					// reset tree_table terminal_nodes[k]^th row, 6^th column equal to new_mean[k] (new_mean is an input vector)
    IntegerVector term_obs2=wrap(term_obs);							// create R object from term_obs
    new_preds[term_obs2]=new_mean[k];								// reset the predictions for observations in node k to the new mean for node k
  }
  updated_preds[0]=tree_table;		// first element of updated_preds is the updated tree table.
  updated_preds[1]=new_preds;			// second element of updated_preds is the updated vector of predictions for all observations.
  
  return(updated_preds);				// return list updated_preds.
}
//######################################################################################################################//

using namespace Rcpp;
using namespace std;

const double flagval = __DBL_MIN__; 
inline double flag(double a, bool b) { return b ? a : flagval; }

// [[Rcpp::export]]
NumericVector subsetter_bcf(NumericVector a, LogicalVector b) { 		// This funtion can be found at http://gallery.rcpp.org/articles/stl-transform-for-subsetting/index.html
  NumericVector a1=clone(a);										// copy input vector
  transform(a1.begin(), a1.end(), b.begin(), a1.begin(), flag);	// Mark values of a1 for which b is false. Applies operation flag sequentially to the elements of a1 and b (two arguments). Stores the output in a1.
  NumericVector res = NumericVector(sum(b));						// res is a vector equal to the number of true values
  remove_copy(a1.begin(), a1.end(), res.begin(), flagval);		// put all the values of a1 corresponding to true values of b into res
  
  return res;    		// return a vector of the true values of a
}
//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector order_inc__bcf(NumericVector x) {
  NumericVector sorted = clone(x).sort();		// sorts ascending. Smallest to largest.
  return match(sorted, x);					// returns indices of position of smallest, followed by position of second smallest, and so on.
}

//######################################################################################################################//

// [[Rcpp::export]]
List min_which2_bcf(NumericVector array,int n,double minout,int whichout){
  // Function to find minimum of an array with n elements that is put in min
  minout=array[0];				// let minout equal first element of input vector "array". Input is unnecessary? Reset before the input value is used.
  whichout=0;						// reset input variable whichout to 0. Input is unnecessary? Reset before the input value is used.
  
  for(int i=1;i<n;i++){			// for-loop of length n-1.
    if(array[i]< minout){		// if i+1^th element of array is less than minout (input variable).
      minout= array[i];		// reset minout equal to the i+1^th element of array. (Equals minimum of all elements of array).
      whichout=i;				// set whichout equal to i. Indexes the position in array of the minumum of all elements of array.
    }
  }
  List ret(2);					// list with two elements.
  ret[0]=minout;					// minimum of all elements of array.
  ret[1]=whichout;				// Indexes the position in array of the minumum of all elements of array.
  
  return(ret);					// return the list.
}
//######################################################################################################################//

#include <Rmath.h>
// [[Rcpp::export]]
double mll_meanvar2_bcf(double x, double x2, int n){
  double sigsq=(x2-((x*x)/n))/n;				// sets value of sigsq.
  if(sigsq<=0){sigsq=0.00000000001;}			// if sigsq <=0, reset it to very small positive number.
  
  return(n*(log(2*M_PI)+log(sigsq)+1)); /* M_PI is in Rmath.h  */
}
//######################################################################################################################//

// [[Rcpp::export]]
IntegerVector PELT_meanvar_norm2_bcf(NumericVector resp,double pen)
  // 0 by default, nonzero indicates error in code 
{
  int n=resp.size();								// n is length of input vector resp.
  NumericVector y2=cumsum(pow(resp,2));			// y2 is vector of cumulative sums of squares of elements of resp.
  y2.push_front(0);								// add a 0 to the front of vector y2. First element.
  NumericVector y=cumsum(resp);					// y is the cumulative sum of elements of resp.
  y.push_front(0);								// add a 0 to the front of vector y. First element.
  IntegerVector cptsout(n,0);    					// cptsout is a vector of n zeroes.
  IntegerVector lastchangecpts(2*(n+1));			// lastchangecpts is a vector of length 2*(n+1). (initially zeroes)
  NumericVector lastchangelike(n+1);				// lastchangelike is a vector of length n+1. (initially zeroes)
  IntegerVector checklist(n+1);					// checklist is a vector of length n+1. (initially zeroes)
  int nchecklist;									// create variable nchecklist.
  double minout;									// create variable minout.
  NumericVector tmplike(n+1);						// tmplike is a vector with n+1 elements.
  IntegerVector tmpt(n+1);						// tmpt is a vector with n+1 elements.
  int tstar,i,whichout,nchecktmp;					// create a set of integers.
  double mll_meanvar();
  void min_which();
  lastchangelike[0]= -pen;						// first element of lastchangelike is negative of input pen.
  lastchangecpts[0]=0; lastchangecpts[n]=0;		// set first and n+1^th element of lastchangecpts to 0.
  double x=y[1];									// let x be second element of y.
  double x2=y2[1];								// let x2 be second element of y2.
  lastchangelike[1]=mll_meanvar2_bcf(x,x2,1);			// apply function defined on line 1382.
  lastchangecpts[1]=0; lastchangecpts[n+1]=1;		// set second element to 0. Set n+2^th element to 1
  lastchangelike[2]=mll_meanvar2_bcf(y[2],y2[2],2);	// apply function defined on line 1382.
  lastchangecpts[2]=0; lastchangecpts[n+2]=2;		// set third element to 0. Set n+3^th element to 2
  lastchangelike[3]=mll_meanvar2_bcf(y[3],y2[3],3);	// apply function defined on line 1382.
  lastchangecpts[3]=0; lastchangecpts[n+3]=3;		// set fourth element to 0. Set n+4^th element to 3
  
  minout=lastchangelike[checklist[0]] + mll_meanvar2_bcf(x,x2,0)+pen;		// minoyt is first element of lastchangelike  + mll_meanvar2_bcf(x,x2,0)+pen
  whichout=0;				// whichout set to 0
  
  nchecklist=2;			// nchecklist set to 2
  checklist[0]=0;			// first element of checklist set to 0
  checklist[1]=2;			// second element of checklist set to 2
  
  for(tstar=4;tstar<(n+1);tstar++){					// loop of length n-3. n is  length of input resp
    R_CheckUserInterrupt(); // checks if user has interrupted the R session and quits if true 
    
    for(i=0;i<nchecklist;i++){						// loop of length nchecklist
      tmplike[i]=lastchangelike[checklist[i]] + mll_meanvar2_bcf(y[tstar]- y[checklist[i]],y2[tstar]-y2[checklist[i]],tstar-checklist[i])+pen;
    }
    List mw=min_which2_bcf(tmplike,nchecklist,minout,whichout); // updates minout and whichout with min and which element 
    NumericVector tempmin=mw[0];
    minout=tempmin[0];
    lastchangelike[tstar]=minout;
    whichout=mw[1];
    lastchangecpts[tstar]=checklist[whichout]; 
    lastchangecpts[n+tstar]=tstar;
    
    // Update checklist for next iteration, first element is next tau 
    nchecktmp=0;
    for(i=0;i<nchecklist;i++){
      if(tmplike[i]<= (lastchangelike[tstar]+pen)){
        checklist[nchecktmp]=checklist[i];
        nchecktmp+=1;
      }
    }
    checklist[nchecktmp]=tstar-1;  // atleast 2 obs per seg
    nchecktmp+=1;
    nchecklist=nchecktmp;
  } // end taustar
  
  // put final set of changepoints together
  int ncpts=0;
  int last=n;
  while(last!=0){
    cptsout[ncpts]=lastchangecpts[n+last];
    last=lastchangecpts[last];
    ncpts+=1;
  }
  
  IntegerVector cptsoutret=cptsout[cptsout>0];
  std::sort(cptsoutret.begin(), cptsoutret.end());
  
  return(cptsoutret);
}

//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double SS_bcf(arma::vec x, arma::vec y, double split){            
  double meanLeft, meanRight;									// create two double variables
  int n = x.n_rows;											// n is length of input vector x
  arma::vec meanAll(n);										// meanAll is a vector of length n
  arma::vec xLeft = y.elem(arma::find(x <= split));			// xLeft is vector of elements of y that correspond to x<= split.
  arma::vec xRight = y.elem(arma::find(x > split));			// xRight is vector of elements of y than correspond to x>split
  meanLeft = mean(xLeft);										// mean of observations sent to left of split.
  meanRight = mean(xRight);									// mean of observations sent to right of split.
  
  for(int j=0;j<n;j++) {										// for-loop of length of n.
    
    if(x[j]<=split){										// if observatino is to left of split
      meanAll[j] = meanLeft;								// set meanAll (the prediction) to meanLeft
    }else{													// if observation is to right of split
      meanAll[j] = meanRight;								// set menAll (the prediction) to meanRight
    }  
  }
  arma::vec test_resids=y-meanAll;							// test_resids is difference between outcome and prediction
  double tSS=as<double>(wrap(trans(y-meanAll)*(y-meanAll)));	// tSS is sum of squared residuals
  
  return tSS;			// retuen the double (sum of squared residuals).
}
//######################################################################################################################//
// [[Rcpp::export]]

List gridCP_bcf(arma::vec x, arma::vec y, int gridSize = 10) {		// grid of change points??
  
  NumericVector out(gridSize-2);										// vector of length gridSize-2
  NumericVector cp_strength(gridSize-2);								// vector of length gridSize-2
  arma::vec no_split(gridSize-2);										// arma vec of length gridSize-2
  double xGrid, gridStep = (max(x)-min(x))/((double)gridSize-1.0);	// double xGrid uninitialized. Double gridstep is range of x divided by gridSize-1? [size of steps between potential spli points?]
  xGrid = min(x);														// set xGrid equal to minimum of x
  List currSS(2);														// list with 2 elements
  
  for(int i=1;i<(gridSize-1);i++) {									// loop of length gridSize-1. Note loop starts at i=1, not i=0. (loop of potential split points).
    xGrid += gridStep;												// add gridStep to xGrid in each iteration.
    arma::vec ld_size= y.elem(arma::find(x <= xGrid));				// vector of elements of y for observations sent left.
    arma::vec rd_size=y.elem(arma::find(x > xGrid));    			// vector of elements of y for observations sent right.
    
    if(ld_size.size()>2 && rd_size.size()>2)						// if more than 2 observations sent left AND more than 2 observations sent right.
    {
      out[i-1]=xGrid;												// i^th element of out set equal to xGrid, which will be i*gridStep for iteration i
      double testSS=SS_bcf(x,y,xGrid);								// sum of squared residusals for i^th potential splitting point. Function SS_bcf Defined on line 1475.
      cp_strength[i-1] = testSS;									// set i^th element of cp_strength equal to SSR for i^th potential splitting point.
      no_split[i-1]=0;											// set i^th element of no_split equal to 0 (to indicate that potential split i is viable in that at least 2 observations went left and two went right)
    }else{															// if less than 3 observations sent left or less than 3 observations sent right
      no_split[i-1]=1;											// set i^th element of no_split equal to 1 (to indicate that potential split i is not viable in that the resulting left and/or right node has less than the minimum number of observations)
    }
  }
  arma::uvec to_remove=find(no_split ==1);							// indexes of split points to be removed (because not enough observations sent left or right).
  IntegerVector remove_order_index=as<IntegerVector>(wrap(order__bcf(as<NumericVector> (wrap(to_remove)))));	// apply function order. gives index of larest element, then index of second largest elemtne etc. But to_remove is already a vector of indexes? So just reverse and converts to IntegerVector?
  
  if(to_remove.size()>0){											// if at least one potential split must be removed (because not enough ons sent left or right)
    for(int k=0;k<remove_order_index.size();k++){				// loop of length equal to the numebr of potential splits to be removed
      out.erase(to_remove[remove_order_index[k]-1]);			// remove one of the potential splits from out
      cp_strength.erase(to_remove[remove_order_index[k]-1]);	// remove one of the potential splits from cp_strength
    }
  }
  if(out.size()>0){				// if at least some viable splitting points (i.e. some potential splits send above the minimum bout left and right).
    currSS[0]= out;				// first element of list is the vector of viable splitting values of x
    currSS[1]= cp_strength;		// secpnd element of list is the RSS from the resulting split.
    
    return currSS;				// return the list currSS
  }else{							// if there are no potential splits.
    List ret(2);				// create a list of length 2.
    ret[0]=0;					// set first element of list equal to 0.
    ret[1]=0;					// set second element of list equal to 0.
    
    return ret;					// return the list.
  }
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_gridpoint_cpmat_mu_bcf(NumericMatrix data,NumericVector resp,int gridsize,int num_cp){
  int err1=0;											// create variable err1. Initialize equal to 0.
  int numvar = data.ncol();							// numvar is number of columns in the input matrix data.
  int numrows=gridsize*numvar-1;						// numrows is gridsize*numvar-1. ??Number of poential splitting points (per variable) by number of variables minus 1.
  int numcols=3;										// create variable numcols. Initialize equal to 0.
  arma::mat change_points(numrows,numcols);			// create arma mat changepoints. Dimensions numrows by numcols.
  int row_count=-1;									// create variable row_count. Initialize equal to -1.
  IntegerVector sorted_var;							// create IntegerVector sorted_var. Length not specified.
  
  for(int i=0;i<data.ncol();i++){						// for-loop of length equal to number of columns in the input data (number of variables).
    NumericVector coli=data(_,i);					// coli is i+1^th column of data.	(data is the X matrix)
    arma::vec coli_ar=as<arma::vec>(coli);			// convert coli to  arma vec.
    arma::vec t=Rcpp::as<arma::vec>(resp);			// t is input resp converted to arma vec. (resp is the residuals?)
    List ans1 = gridCP_bcf(coli,t,gridsize);			// apply gridCP_bcf (defined on line 1500) to column i+1. returns list of two elements.
    NumericVector ans2=ans1[0];						// vector of viable splitting values of the variable in the i+1^th column. 
    NumericVector cp_strength=ans1[1];				// vector of residual sums of squares from the splits.
    
    if(ans1.size()!=1 && ans2[0]!=0){				// if ans1 has more than 1 element AND there are at least some viable splitting points. [how can gridCP_bcf output have anything other than 2 elements? is the first condition necessary?]
      for(int j=0;j<ans2.size();j++){				// for-loop of length equal to number of viable potential splitting points.
        row_count+=1;									// add 1 to row_count in each iteration. Why not use j instead of row_count?
        change_points(row_count,0)=i;					// let the row_count+1^th row, first column entry of change_points equal the variable number (number of column in data)
        change_points(row_count,1)=ans2[j];				// let the row_count+1^th row, second column entry of change_points equal the j+1^th viable splitting value of the variable.
        change_points(row_count,2)=cp_strength[j];		// let the row_count+1^th row, third column entry of change_points equal the resisudual sum of squares from the split corresponding to the row
      }   
    }
  }
  if(row_count+1!=(int) change_points.n_rows){						// if row_count+1 equals numbreof rows of change_points (note change_points initially set to have numbr of rows equal to a number >= max numebr of potential splits. There could be less rows because of non-viable splits).
    change_points.shed_rows(row_count+1,change_points.n_rows-1);	// remove all extra rows that were no filled in by the for-loops above.
  }
  arma::vec te=change_points.col(2);									// set te equal to third column of change_points. The residuals SS_bcf for all potential splits.
  NumericVector col_to_order=as<NumericVector>(wrap(te));				// convert residual SS_bcf vector to NumericVector.
  IntegerVector ordered_dev=order_inc__bcf(col_to_order);					// gives indices of position of smallest, followed by position of second smallest element, and so on/
  ordered_dev=ordered_dev-1;											// take one away from all indices (so that they are between 0 and number of splits -1).
  change_points.shed_col(2);											// removes third column.
  int cp=change_points.n_rows;										// number of rows of change_points. Number of potential splits.
  double cp_prop=(double)num_cp/(double)100;							// cp_prop is the input value num_cp divided by 100. [presumably the fraction of possible splits to keep]
  int num_cp2=round(cp*(cp_prop));									// Nultiply the number of possible splits by cp_prop. Gives number of potential splits to keep? Rounded to nearest integer.
  num_cp=round(num_cp2);												// re-set num_cp to the number of potential splits kept. Why round again?
  
  if(num_cp==0 && cp!=0){												// This condition only occurs if original input value of num_cp is 0?
    num_cp=cp;														// set num_cp equal to cp (number of potential splits).
  }
  if(cp<num_cp){														// This condition only occurs if original value of num_cp >=100?
    num_cp=change_points.n_rows;									// set num_cp equal to the nnumber of potential splitting variables given in change_points.
  }
  if(num_cp==0){														// This occurs if cp=0? no rows filled in for change_points?
    err1=1;															// re-set the integer err1 equal to 1. This could be a bool instead of an int? Records an error?
  }
  if(err1==0){														// if err1=0, i.e. num_cp !=0.
    arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);					// indices of smallest residual SS_bcf, then second smallest, and so on.
    t=t.subvec(0,num_cp-1);											// keep indices for the best num_cp potential splits.
    change_points=change_points.rows(t);							// Keep rows of matrix change_points that correspond to best potential splits.
  }
  List ret(2);					// list of length 2
  ret[0]=wrap(change_points);		// first element is a matrix (2 columns). Rows correspond to different potential splits. First column gives column index of splitting variables in the data matrix. Second column gives splitting point values.
  ret[1]=err1;					// second element is 0 if no error. 1 if there is an error? (i.e. if no potential splits are returned?).
  
  return(ret);					// return the list
}
//######################################################################################################################//

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_gridpoint_cpmat_tau_bcf(NumericMatrix data,NumericVector resp,int gridsize,int num_cp, NumericVector z){
  
  arma::mat x_moderate_a(data.begin(), data.nrow(), data.ncol(), false);	
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  arma::vec resp_a=Rcpp::as<arma::vec>(resp);		// converts to arma vec
  
  // create matrix of data for treated households only
  arma::mat x_mod_treat_a = x_moderate_a.rows(arma::find(z_ar==1));	
  NumericMatrix x_mod_treat = wrap(x_mod_treat_a);	
  
  // create vector of residuals for treated households only
  arma::vec resp_treat_a = resp_a.elem(arma::find(z_ar==1));	
  NumericVector resp_treat = wrap(resp_treat_a);
  
  
  int err1=0;											// create variable err1. Initialize equal to 0.
  int numvar = x_mod_treat.ncol();							// numvar is number of columns in the input matrix x_mod_treat.
  int numrows=gridsize*numvar-1;						// numrows is gridsize*numvar-1. ??Number of poential splitting points (per variable) by number of variables minus 1.
  int numcols=3;										// create variable numcols. Initialize equal to 0.
  arma::mat change_points(numrows,numcols);			// create arma mat changepoints. Dimensions numrows by numcols.
  int row_count=-1;									// create variable row_count. Initialize equal to -1.
  IntegerVector sorted_var;							// create IntegerVector sorted_var. Length not specified.
  
  for(int i=0;i<x_mod_treat.ncol();i++){						// for-loop of length equal to number of columns in the input x_mod_treat (number of variables).
    NumericVector coli=x_mod_treat(_,i);					// coli is i+1^th column of x_mod_treat.	(x_mod_treat is the X matrix)
    arma::vec coli_ar=as<arma::vec>(coli);			// convert coli to  arma vec.
    arma::vec t=Rcpp::as<arma::vec>(resp_treat);			// t is input resp_treat converted to arma vec. (resp_treat is the residuals?)
    List ans1 = gridCP_bcf(coli,t,gridsize);			// apply gridCP_bcf (defined on line 1500) to column i+1. returns list of two elements.
    NumericVector ans2=ans1[0];						// vector of viable splitting values of the variable in the i+1^th column. 
    NumericVector cp_strength=ans1[1];				// vector of residual sums of squares from the splits.
    
    if(ans1.size()!=1 && ans2[0]!=0){				// if ans1 has more than 1 element AND there are at least some viable splitting points. [how can gridCP_bcf output have anything other than 2 elements? is the first condition necessary?]
      for(int j=0;j<ans2.size();j++){				// for-loop of length equal to number of viable potential splitting points.
        row_count+=1;									// add 1 to row_count in each iteration. Why not use j instead of row_count?
        change_points(row_count,0)=i;					// let the row_count+1^th row, first column entry of change_points equal the variable number (number of column in x_mod_treat)
        change_points(row_count,1)=ans2[j];				// let the row_count+1^th row, second column entry of change_points equal the j+1^th viable splitting value of the variable.
        change_points(row_count,2)=cp_strength[j];		// let the row_count+1^th row, third column entry of change_points equal the resisudual sum of squares from the split corresponding to the row
      }   
    }
  }
  if(row_count+1!=(int) change_points.n_rows){						// if row_count+1 equals numbreof rows of change_points (note change_points initially set to have numbr of rows equal to a number >= max numebr of potential splits. There could be less rows because of non-viable splits).
    change_points.shed_rows(row_count+1,change_points.n_rows-1);	// remove all extra rows that were no filled in by the for-loops above.
  }
  arma::vec te=change_points.col(2);									// set te equal to third column of change_points. The residuals SS_bcf for all potential splits.
  NumericVector col_to_order=as<NumericVector>(wrap(te));				// convert residual SS_bcf vector to NumericVector.
  IntegerVector ordered_dev=order_inc__bcf(col_to_order);					// gives indices of position of smallest, followed by position of second smallest element, and so on/
  ordered_dev=ordered_dev-1;											// take one away from all indices (so that they are between 0 and number of splits -1).
  change_points.shed_col(2);											// removes third column.
  int cp=change_points.n_rows;										// number of rows of change_points. Number of potential splits.
  double cp_prop=(double)num_cp/(double)100;							// cp_prop is the input value num_cp divided by 100. [presumably the fraction of possible splits to keep]
  int num_cp2=round(cp*(cp_prop));									// Nultiply the number of possible splits by cp_prop. Gives number of potential splits to keep? Rounded to nearest integer.
  num_cp=round(num_cp2);												// re-set num_cp to the number of potential splits kept. Why round again?
  
  if(num_cp==0 && cp!=0){												// This condition only occurs if original input value of num_cp is 0?
    num_cp=cp;														// set num_cp equal to cp (number of potential splits).
  }
  if(cp<num_cp){														// This condition only occurs if original value of num_cp >=100?
    num_cp=change_points.n_rows;									// set num_cp equal to the nnumber of potential splitting variables given in change_points.
  }
  if(num_cp==0){														// This occurs if cp=0? no rows filled in for change_points?
    err1=1;															// re-set the integer err1 equal to 1. This could be a bool instead of an int? Records an error?
  }
  if(err1==0){														// if err1=0, i.e. num_cp !=0.
    arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);					// indices of smallest residual SS_bcf, then second smallest, and so on.
    t=t.subvec(0,num_cp-1);											// keep indices for the best num_cp potential splits.
    change_points=change_points.rows(t);							// Keep rows of matrix change_points that correspond to best potential splits.
  }
  List ret(2);					// list of length 2
  ret[0]=wrap(change_points);		// first element is a matrix (2 columns). Rows correspond to different potential splits. First column gives column index of splitting variables in the x_mod_treat matrix. Second column gives splitting point values.
  ret[1]=err1;					// second element is 0 if no error. 1 if there is an error? (i.e. if no potential splits are returned?).
  
  return(ret);					// return the list
}
//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_pelt_cpmat_mu_bcf(NumericMatrix data,NumericVector resp,double pen,int num_cp){
  int err1=0;									// create variable err1. Initialize equal to 0.
  int n=data.nrow();							// n i number of rows in input data.
  arma::mat change_points;					// create an arma mat change_points. No initial dimensions.
  change_points.set_size(10000000,3);			// set dimensions of change_pints to 10 million rows, 3 columns. Do not initialize the values.
  int row_count=-1;							// create variable row_count. Initialize equal to -1.
  IntegerVector sorted_var;					// create an IntegerVector sorted_var. Length not set.
  
  for(int i=0;i<data.ncol();i++){					// Loop of length equal to number of columns in data. (Number of variables).
    NumericVector coli=data(_,i);				// coli is i+1^th column of data. (vector of observations of i+1^th covariate)
    arma::vec coli_ar=as<arma::vec>(coli);		// create an arma vec copy of coli
    sorted_var=order_inc__bcf(coli);				// gives indices of position of smallest, followed by position of second smallest element, and so on. 
    sorted_var=sorted_var-1;					// takes one away from indices so that they are between 0 and number of observations minus 1.
    IntegerVector ans = PELT_meanvar_norm2_bcf(resp[sorted_var],pen);			// PELT_meanvar_norm2_bcf defined on line 1391. Not sure what this does.
    ans=ans-1;									// Presumably ans is a set of indices for changepoints and 1 is taken away so that the indices have smallest value 0.
    arma::vec resp_ar=as<arma::vec>(resp);		// copy input resp to an arma vec.
    NumericVector allMean(n);					// vector of length equal to number of observations in data.
    
    if(ans.size()!=1 && ans[0]!=n){							// if ans is not of length 1 and the first element of ans is not equal to the number of observations.
      double meanRight;									// create variable meanRight. Not initialized.
      double meanLeft;									// create variable meanLeft. Not initialized.
      IntegerVector pelt_sp=sorted_var[ans];				// elements of sorted_var indexed by ans. ??Indices of splitting points chosen by PELT algorithm? (for variable i+1)
      NumericVector split_values=coli[pelt_sp];			// values of splitting points chosen by PELT algorithm? (for variable i+1)
      NumericVector unique_sp=unique(split_values);		// keep unique split values. unique sorts in descending order.
      NumericVector all_sp=coli[pelt_sp];					// values of splitting points chosen by PELT algorithm? (for variable i+1). This seems unnecessary. Same as split_values. Neither is changed in later lines within the loop.
      
      if(unique_sp.size()<all_sp.size()){					// If some split values were not unique
        IntegerVector index=match(unique_sp,all_sp);	// positions of first matches of elements of unique in all_sp.
        arma::vec indexarma=as<arma::vec>(index);		// convert index to arma vec. Why not arma uvec?
        IntegerVector t4=ifelse(is_na(index),1,0);		// t4 is a vector, equal to 1 if the corresponding value of index is NA. Otherwise equal to 0. How can an NA be returned?
        arma::vec t4arma=as<arma::vec>(t4);				// copy t4 as arma vec.
        arma::uvec t42=find(t4arma==0);					// find indices of elements of t4arma that are not equal to 0. (indices of elemetns that are not NA.
        IntegerVector index2(t42.size());				// create variable index2 of length equal to the number of elements that are not NA.
        
        arma::vec index3=indexarma.elem(t42);			// index3 is the not NA elements of index.
        
        index2=wrap(index3);							// convert index3 to an R object.
        index2=index2-1;								// take one away from the indices
        ans=index2;										// set ans equal to index2
      }
      
      for(int j=0;j<ans.size()-1;j++){								// loop of length equal to that of ans MINUS ONE. (kept splitting points minus 1?)
        arma::vec xLeft;											// create arma vec xLeft
        arma::vec xRight;											// create arma vec xRight
        double splitpoint;											// create variable splitpoint
        if(unique_sp.size()<all_sp.size()){							// If some split values were not unique
          splitpoint = all_sp[j];									// value of j+1^th split point (chosen by PELT)
          xLeft = resp_ar.elem(arma::find(coli_ar<=splitpoint));	// observations (in input resp) sent left by the j+1^th splitpoint.
          xRight = resp_ar.elem(arma::find(coli_ar>splitpoint));	// observations (in input resp) sent right by the j+1^th splitpoint.
        }else{														// If all split values were unique
          splitpoint = coli[sorted_var[ans[j]]];					// value of j+1^th split point (chosen by PELT)
          xLeft = resp_ar.elem(arma::find(coli_ar<=splitpoint));	// observations (in input resp) sent left by the j+1^th splitpoint.
          xRight = resp_ar.elem(arma::find(coli_ar>splitpoint));	// observations (in input resp) sent right by the j+1^th splitpoint.
        }
        if(xLeft.size()>2 && xRight.size()>2)						// If more than 2 obs sent left AND more than 2 obs sent right.
        {
          row_count +=1;											// increment row_count by 1 for each iteration (note: keeps increasing across splitting points, and different potetial splitting variables)
          meanLeft = mean(xLeft);									// find the mean for observations sent left.
          meanRight = mean(xRight);								// find the mean for observations sent right.
          arma::uvec left_ind=find(coli_ar<=splitpoint);			// indices for covariate values below the splitting point
          arma::uvec right_ind=find(coli_ar>splitpoint);			// indices for covariate values above the splitting point
          NumericVector left_ind2=wrap(left_ind);					// convert indices of obs left of split to NumericVector. Why not IntegerVector?
          NumericVector right_ind2=wrap(right_ind);				// convert indices of obs right of split to NumericVector. Why not IntegerVector?
          allMean[left_ind2]=meanLeft;							// set allMean for obs sent left to mean of outcomes for obs sent left. (left node predictions).
          allMean[right_ind2]=meanRight;							// set allMean for obs sent right to mean of outcomes for obs sent right. (right node predictions).
          arma::vec allMean_ar=as<arma::vec>(allMean);			// create an arma vec copy of allMean.
          
          double cp_strength;																// create a double called cp_strength.
          cp_strength=as<double>(wrap(trans(resp_ar-allMean_ar)*(resp_ar-allMean_ar)));	// cp_strength is sum of squared residuals for the proposed splitting point (change-point).
          change_points(row_count,0)=i;													// let first column of change_points index the variable (the column in data for the variable, beginning at column 0).
          change_points(row_count,1)=coli[sorted_var[ans[j]]];							// the second column of change_points is the covariate value for the splitpoint.
          change_points(row_count,2)=cp_strength;											// the third column of change_points is the sum of squared residuals for the predictions resulting from the split.
        }
      }          
    }
  }
  change_points.shed_rows(row_count+1,change_points.n_rows-1);	//remove rows of change_points that were not filled in.
  arma::vec te=change_points.col(2);								// te is an arma vec copy of the third column of change_points. This is the sums of squared errors.
  NumericVector col_to_order=as<NumericVector>(wrap(te));			// convert te to a NumericVector called col_to_order.
  IntegerVector ordered_dev=order_inc__bcf(col_to_order);				// returns indices of position of smallest SS_bcf residuals, followed by position of second smallest, and so on.
  ordered_dev=ordered_dev-1;										// takes one away from the indices (so that they are between 0 and the number of potential splits minus 1).
  change_points.shed_col(2);										// remove the third column from change_points.
  int cp=change_points.n_rows;									// cp is the number of rows in change_points. the number of selected potential splitting points across all variables.
  double cp_prop=(double)num_cp/(double)100;						// cp_prop is input value num_cp divided by 100. Perhaps this is the proportion of potential splits to keep?
  int num_cp2=round(cp*(cp_prop));								// number of potential splits to keep? (rounded to nearest integer)
  num_cp=round(num_cp2);											// (rounding in this line is probably unnecessary).
  if(num_cp==0 && cp!=0){											// num_cp equals zero, but the number of potential splits is nonzero (i.e. the original entered value of num_cp is 0... or cp and the original input value of num_cp are very small and nonzero such that the product cp*cp_prop is a small positive number, which is rounded down to zero, making num_cp then be set equal to 0)
    num_cp=cp;													// set num_cp equal to cp. The number of potential splits to keep equals the total number of potential splits.
  }
  if(cp<num_cp){													// if cp is less than num_cp. (?num_cp has original value greater than 100?)
    num_cp=change_points.n_rows;								// keep all potential splits. Why not use cp here instead of change_points.n_rows
  }
  if(num_cp==0){													// This occurs if cp is 0(given the resetting in the if statement above), i.e. no viable potential splitting points returned.
    err1=1;														// re-set err1 equal to 1. (record an error?)
  }
  if(err1==0){													// If err1 equals 0 (no error?)
    arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);				// convert indices of positions of splits (ordered descending by SSR) to arma uvec.
    t=t.subvec(0,num_cp-1);										// keep the indices of the best num_cp splits.
    change_points=change_points.rows(t);						// keep the rows in change_points for only the best num_cp splits.
  }
  List ret(2);					// create a list of length 2.
  ret[0]=wrap(change_points);		// convert the matrix to an R object. First column indexes splitting variable (gives column number in data matrix, starting at zero). Second column gives covariate value for split point.
  ret[1]=err1;					// Integer equals 1 to record an error, 0 otherwise.
  
  return(ret);					// return the list.
}

//######################################################################################################################//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List make_pelt_cpmat_tau_bcf(NumericMatrix data,NumericVector resp,double pen,int num_cp, NumericVector z){
  
  arma::mat x_moderate_a(data.begin(), data.nrow(), data.ncol(), false);	
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  arma::vec resp_a=Rcpp::as<arma::vec>(resp);		// converts to arma vec
  
  // create matrix of data for treated households only
  arma::mat x_mod_treat_a = x_moderate_a.rows(arma::find(z_ar==1));	
  NumericMatrix x_mod_treat = wrap(x_mod_treat_a);	
  
  // create vector of residuals for treated households only
  arma::vec resp_treat_a = resp_a.elem(arma::find(z_ar==1));	
  NumericVector resp_treat = wrap(resp_treat_a);
  
  int err1=0;									// create variable err1. Initialize equal to 0.
  int n=x_mod_treat.nrow();							// n i number of rows in input x_mod_treat.
  arma::mat change_points;					// create an arma mat change_points. No initial dimensions.
  change_points.set_size(10000000,3);			// set dimensions of change_pints to 10 million rows, 3 columns. Do not initialize the values.
  int row_count=-1;							// create variable row_count. Initialize equal to -1.
  IntegerVector sorted_var;					// create an IntegerVector sorted_var. Length not set.
  
  for(int i=0;i<x_mod_treat.ncol();i++){					// Loop of length equal to number of columns in x_mod_treat. (Number of variables).
    NumericVector coli=x_mod_treat(_,i);				// coli is i+1^th column of x_mod_treat. (vector of observations of i+1^th covariate)
    arma::vec coli_ar=as<arma::vec>(coli);		// create an arma vec copy of coli
    sorted_var=order_inc__bcf(coli);				// gives indices of position of smallest, followed by position of second smallest element, and so on. 
    sorted_var=sorted_var-1;					// takes one away from indices so that they are between 0 and number of observations minus 1.
    IntegerVector ans = PELT_meanvar_norm2_bcf(resp_treat[sorted_var],pen);			// PELT_meanvar_norm2_bcf defined on line 1391. Not sure what this does.
    ans=ans-1;									// Presumably ans is a set of indices for changepoints and 1 is taken away so that the indices have smallest value 0.
    arma::vec resp_treat_ar=as<arma::vec>(resp_treat);		// copy input resp_treat to an arma vec.
    NumericVector allMean(n);					// vector of length equal to number of observations in x_mod_treat.
    
    if(ans.size()!=1 && ans[0]!=n){							// if ans is not of length 1 and the first element of ans is not equal to the number of observations.
      double meanRight;									// create variable meanRight. Not initialized.
      double meanLeft;									// create variable meanLeft. Not initialized.
      IntegerVector pelt_sp=sorted_var[ans];				// elements of sorted_var indexed by ans. ??Indices of splitting points chosen by PELT algorithm? (for variable i+1)
      NumericVector split_values=coli[pelt_sp];			// values of splitting points chosen by PELT algorithm? (for variable i+1)
      NumericVector unique_sp=unique(split_values);		// keep unique split values. unique sorts in descending order.
      NumericVector all_sp=coli[pelt_sp];					// values of splitting points chosen by PELT algorithm? (for variable i+1). This seems unnecessary. Same as split_values. Neither is changed in later lines within the loop.
      
      if(unique_sp.size()<all_sp.size()){					// If some split values were not unique
        IntegerVector index=match(unique_sp,all_sp);	// positions of first matches of elements of unique in all_sp.
        arma::vec indexarma=as<arma::vec>(index);		// convert index to arma vec. Why not arma uvec?
        IntegerVector t4=ifelse(is_na(index),1,0);		// t4 is a vector, equal to 1 if the corresponding value of index is NA. Otherwise equal to 0. How can an NA be returned?
        arma::vec t4arma=as<arma::vec>(t4);				// copy t4 as arma vec.
        arma::uvec t42=find(t4arma==0);					// find indices of elements of t4arma that are not equal to 0. (indices of elemetns that are not NA.
        IntegerVector index2(t42.size());				// create variable index2 of length equal to the number of elements that are not NA.
        
        arma::vec index3=indexarma.elem(t42);			// index3 is the not NA elements of index.
        
        index2=wrap(index3);							// convert index3 to an R object.
        index2=index2-1;								// take one away from the indices
        ans=index2;										// set ans equal to index2
      }
      
      for(int j=0;j<ans.size()-1;j++){								// loop of length equal to that of ans MINUS ONE. (kept splitting points minus 1?)
        arma::vec xLeft;											// create arma vec xLeft
        arma::vec xRight;											// create arma vec xRight
        double splitpoint;											// create variable splitpoint
        if(unique_sp.size()<all_sp.size()){							// If some split values were not unique
          splitpoint = all_sp[j];									// value of j+1^th split point (chosen by PELT)
          xLeft = resp_treat_ar.elem(arma::find(coli_ar<=splitpoint));	// observations (in input resp_treat) sent left by the j+1^th splitpoint.
          xRight = resp_treat_ar.elem(arma::find(coli_ar>splitpoint));	// observations (in input resp_treat) sent right by the j+1^th splitpoint.
        }else{														// If all split values were unique
          splitpoint = coli[sorted_var[ans[j]]];					// value of j+1^th split point (chosen by PELT)
          xLeft = resp_treat_ar.elem(arma::find(coli_ar<=splitpoint));	// observations (in input resp_treat) sent left by the j+1^th splitpoint.
          xRight = resp_treat_ar.elem(arma::find(coli_ar>splitpoint));	// observations (in input resp_treat) sent right by the j+1^th splitpoint.
        }
        if(xLeft.size()>2 && xRight.size()>2)						// If more than 2 obs sent left AND more than 2 obs sent right.
        {
          row_count +=1;											// increment row_count by 1 for each iteration (note: keeps increasing across splitting points, and different potetial splitting variables)
          meanLeft = mean(xLeft);									// find the mean for observations sent left.
          meanRight = mean(xRight);								// find the mean for observations sent right.
          arma::uvec left_ind=find(coli_ar<=splitpoint);			// indices for covariate values below the splitting point
          arma::uvec right_ind=find(coli_ar>splitpoint);			// indices for covariate values above the splitting point
          NumericVector left_ind2=wrap(left_ind);					// convert indices of obs left of split to NumericVector. Why not IntegerVector?
          NumericVector right_ind2=wrap(right_ind);				// convert indices of obs right of split to NumericVector. Why not IntegerVector?
          allMean[left_ind2]=meanLeft;							// set allMean for obs sent left to mean of outcomes for obs sent left. (left node predictions).
          allMean[right_ind2]=meanRight;							// set allMean for obs sent right to mean of outcomes for obs sent right. (right node predictions).
          arma::vec allMean_ar=as<arma::vec>(allMean);			// create an arma vec copy of allMean.
          
          double cp_strength;																// create a double called cp_strength.
          cp_strength=as<double>(wrap(trans(resp_treat_ar-allMean_ar)*(resp_treat_ar-allMean_ar)));	// cp_strength is sum of squared residuals for the proposed splitting point (change-point).
          change_points(row_count,0)=i;													// let first column of change_points index the variable (the column in x_mod_treat for the variable, beginning at column 0).
          change_points(row_count,1)=coli[sorted_var[ans[j]]];							// the second column of change_points is the covariate value for the splitpoint.
          change_points(row_count,2)=cp_strength;											// the third column of change_points is the sum of squared residuals for the predictions resulting from the split.
        }
      }          
    }
  }
  change_points.shed_rows(row_count+1,change_points.n_rows-1);	//remove rows of change_points that were not filled in.
  arma::vec te=change_points.col(2);								// te is an arma vec copy of the third column of change_points. This is the sums of squared errors.
  NumericVector col_to_order=as<NumericVector>(wrap(te));			// convert te to a NumericVector called col_to_order.
  IntegerVector ordered_dev=order_inc__bcf(col_to_order);				// returns indices of position of smallest SS_bcf residuals, followed by position of second smallest, and so on.
  ordered_dev=ordered_dev-1;										// takes one away from the indices (so that they are between 0 and the number of potential splits minus 1).
  change_points.shed_col(2);										// remove the third column from change_points.
  int cp=change_points.n_rows;									// cp is the number of rows in change_points. the number of selected potential splitting points across all variables.
  double cp_prop=(double)num_cp/(double)100;						// cp_prop is input value num_cp divided by 100. Perhaps this is the proportion of potential splits to keep?
  int num_cp2=round(cp*(cp_prop));								// number of potential splits to keep? (rounded to nearest integer)
  num_cp=round(num_cp2);											// (rounding in this line is probably unnecessary).
  if(num_cp==0 && cp!=0){											// num_cp equals zero, but the number of potential splits is nonzero (i.e. the original entered value of num_cp is 0... or cp and the original input value of num_cp are very small and nonzero such that the product cp*cp_prop is a small positive number, which is rounded down to zero, making num_cp then be set equal to 0)
    num_cp=cp;													// set num_cp equal to cp. The number of potential splits to keep equals the total number of potential splits.
  }
  if(cp<num_cp){													// if cp is less than num_cp. (?num_cp has original value greater than 100?)
    num_cp=change_points.n_rows;								// keep all potential splits. Why not use cp here instead of change_points.n_rows
  }
  if(num_cp==0){													// This occurs if cp is 0(given the resetting in the if statement above), i.e. no viable potential splitting points returned.
    err1=1;														// re-set err1 equal to 1. (record an error?)
  }
  if(err1==0){													// If err1 equals 0 (no error?)
    arma::uvec t=Rcpp::as<arma::uvec>(ordered_dev);				// convert indices of positions of splits (ordered descending by SSR) to arma uvec.
    t=t.subvec(0,num_cp-1);										// keep the indices of the best num_cp splits.
    change_points=change_points.rows(t);						// keep the rows in change_points for only the best num_cp splits.
  }
  List ret(2);					// create a list of length 2.
  ret[0]=wrap(change_points);		// convert the matrix to an R object. First column indexes splitting variable (gives column number in x_mod_treat matrix, starting at zero). Second column gives covariate value for split point.
  ret[1]=err1;					// Integer equals 1 to record an error, 0 otherwise.
  
  return(ret);					// return the list.
}


//###################################################################################//

// [[Rcpp::export]]

List get_best_trees_mu_bcf(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                       double a_mu,double a_tau,double mu_mu,double mu_tau,
                       double nu,double lambda,double c,double sigma_mu_mu,double sigma_mu_tau,
                       List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                       double lowest_BIC,//int first_round,
                       IntegerVector parent,
                       List resids_cp_mat_mu,IntegerVector err_list,NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                       double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,bool is_test_data,
                       double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                       bool split_rule_node,bool gridpoint,int maxOWsize, int num_splits_mu,int num_splits_tau,int gridsize_mu, bool zero_split
){
  List eval_model;										// create a list. Length not given.
  NumericVector lik_list;									// create a vector. Length not givem
  List best_subset;										// create a list. length not givn
  int overall_size=1000;									// create a variable overall_size. Initialized as 1000.
  List overall_trees(overall_size);						// create a list of length 1000.
  NumericVector overall_lik2;								// create a vector. length not given.
  IntegerVector overall_parent2;							// create a vector. length not give,
  List overall_mat(overall_size);							// create a list of length 1000.
  int overall_count=0;									// create a variable overall_count, Initialized equal to 0.
  std::vector<int> overall_parent(overall_size);			// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);			// create a vector of length 1000.
  NumericVector test_preds;								// create a vector. length not given.
  List cp_mat_list=resids_cp_mat_mu;
  
  
  
  if(zero_split==1){
    
  //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  
  // Rcout << "get to likelihood function mu round one. \n";
  overall_trees[0]= tree_table_mu[0];
  overall_mat[0]= tree_mat_mu[0];
  overall_parent[0]=-1;
  overall_parent2[0]=-1;
  double lik_temp=likelihood_function_bcf(resids(_,0),tree_table_mu[0],tree_mat_mu[0],a_mu,mu_mu,nu,lambda);
  double tree_prior_temp=get_tree_prior_bcf(tree_table_mu[0],tree_mat_mu[0],alpha_mu,beta_mu);
  double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(x_control_a.n_rows);
  overall_lik[0]= lowest_BIC_temp;
  overall_count=1;  
  
  }
  
  // Rcout << "get past likelihood function mu round one. \n";
  
  
  
  
  
  
    
    for(int j=0;j<num_splits_mu;j++){									// for-loop of length 5.
      int lsize=1000;										// reate a variable lsize. Initialized as 1000.
      List table_subset_curr_round(lsize);				// create a list of length 1000.
      std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000.
      List mat_subset_curr_round(lsize);					// create a list of length 1000.
      std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
      int count=0;										// create a vector count. Initialize as 0.
      for(int i=0;i<tree_table_mu.size();i++){				// create a length equal to the length of the input list tree_table_mu
        parent=-1;									// reset input vector parent to have one element, equal to -1.
        //NumericMatrix temp_list=cp_mat_list[0];		// indexes columns of splitting variables for potential splitting points. let temp_list equal the first element of the input list cp_mat_list
        best_subset=get_best_split_mu_bcf(resids(_,0),x_control_a,tree_table_mu[i],tree_mat_mu[i],
                                      a_mu,mu_mu,nu,lambda,log(c),lowest_BIC,
                                      parent[0],cp_mat_list[0],
                                                           alpha_mu,beta_mu,maxOWsize//,first_round
                                            ); // defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
        if(best_subset.size()==1){						// if get_best_split_bcf resulted in a (?no tree?) error.
          continue;									// skip to the next iteration of the iner for-loop.
        }
        List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_bcf
        List temp_mat=best_subset[6];					// list of tree matrices
        lik_list=best_subset[5];						// vector of BICs
        IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?)
        if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
          throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
        }
        if(lik_list.size()==0){							// if vector of BICs has no elements
          throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
        }
        if(min(lik_list)<lowest_BIC){					// If the minimum of the list of BICs is less than the input value of lowest_BIC.
          lowest_BIC=min(lik_list);					// reset lowest_BIC equal to the minimum of the list of BICs
        }
        for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
          table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
          lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
          mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
          parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
          count++;										// increment the count by 1.
          
          if(count==(lsize-1)){														// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
            lsize=lsize*2;															// Double the length
            table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
            mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
            lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
            parent_curr_round.resize(lsize);										// double length of parent_curr_round
          }
        }
      }
      table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
      mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
      lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
      parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
      
      if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
        break;															// break out of the for-loop,
      }
      for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
        overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
        overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
        overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
        overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
        overall_count++;												// increment overall_count by one.
        
        if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
          overall_size=overall_size*2;								// double the size
          overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
          overall_lik.resize(overall_size);							// double the size of overall_lik
          overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
          overall_parent.resize(overall_size);						// double the size of overall_parent
        }
      }
      overall_trees=resize_bcf(overall_trees,overall_count);			// remove remaining spaces that are not filled in
      overall_lik.resize(overall_count);							// remove remaining spaces that are not filled in
      overall_mat=resize_bcf(overall_mat,overall_count);				// remove remaining spaces that are not filled in
      overall_parent.resize(overall_count);						// remove remaining spaces that are not filled in
      eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_mu_list, tree_lik, tree_parent
      overall_lik2=eval_model[0];									// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
      overall_trees=eval_model[1];								// list of tree tables in Occam's window
      overall_mat=eval_model[2];									// list of tree matrices in Occa's window
      overall_count=overall_trees.size();							// re-set overall_count to number of models in Occam's window
      overall_parent2=eval_model[3];								// tree parent vector for all models in Occam's window.
      //add in check to see if OW accepted more than the top maxOW models...
      if(overall_lik2.size()>maxOWsize){							// If more than maxOWsize models kept in Occam's window
        IntegerVector owindices=orderforOW__bcf(overall_lik2);			// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
        owindices=owindices-1;									// take one away from indices so that they are in the correct range.
        //get the top maxOWsize indices to keep in OW
        NumericVector temp_olik(maxOWsize);						// create vector of length maxOWsize
        List temp_otrees(maxOWsize);							// create list of length maxOWsize
        List temp_omat(maxOWsize);								// create list of length maxOWsize
        IntegerVector temp_oparent(maxOWsize);					// create vector of length maxOWsize
        
        //now only select those elements
        for(int t=0;t<maxOWsize;t++){							// for-loop of length equal to size of Occam's window
          temp_olik[t]=overall_lik2[owindices[t]];			// keep the BICs (likelihoods?) for the top maxOWsize models
          temp_otrees[t]=overall_trees[owindices[t]];			// keep the tree tables for the top maxOWsize models
          temp_omat[t]= overall_mat[owindices[t]];			// keep the tree matrices for the top maxOWsize models
          temp_oparent[t]=overall_parent2[owindices[t]];		// keep the tree parent vectors for the top maxOWsize models
        }
        
        overall_lik2=temp_olik;									// reset overall_lik2 to keep only top maxOWsize models
        overall_trees=temp_otrees;								// reset overall_trees to keep only top maxOWsize models
        overall_mat=temp_omat;									// reset overall_mat to keep only top maxOWsize models
        overall_count=overall_trees.size();						// reset overall_count to keep only top maxOWsize models
        overall_parent2=temp_oparent;							// reset overall_parent2 to keep only top maxOWsize models
      }
      
      tree_table_mu=table_subset_curr_round;								// reset tree_table_mu equal to the list of tables produced in the current iteration of the outer loop.
      IntegerVector temp1(table_subset_curr_round.size(),1);			// create a vector of ones of length equal to that of the list of tree tables for the current round.
      err_list=temp1;													// set err_list equal to temp1
      if(overall_trees.size()<overall_size-1){						// if length of overall_trees is less than overall_size-1
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// increse size of overall_trees
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// increse size of overall_mat
        overall_lik.resize(overall_size);							// increse size of overall_lik
        overall_parent.resize(overall_size);						// increse size of overall_parent
      }else{															// if length of overall_trees is greater than or equal to overall_size-1
        overall_size=2*overall_size;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the length of overall_trees
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the length of overall_mat
        overall_lik.resize(overall_size);							// double the length of overall_lik
        overall_parent.resize(overall_size);						// double the length of overall_parent
      }
      tree_mat_mu=mat_subset_curr_round;									// reset tree_mat_mu to list of treematrices obtained in current round of for-loop)
      parent=parent_curr_round;										// reset parent to parent vector obtained in current round of for-loop.
      
      if(split_rule_node==1){											// If input boolean split_rule_node = 1 (true)
        NumericVector temp_preds;									// create vector. Length not given.
        List updated_curr_preds;									// create list. Length not given.
        NumericVector new_mean;										// create vector. Length not given.
        lowest_BIC=min(overall_lik2);								// update lowest_BIC
        NumericMatrix curr_resids(resids.nrow(),resids.ncol());		// curr_resids is a matrix with dimensions of input matrix resids
        
        for(int k=0;k<table_subset_curr_round.size();k++){			// for-loop of length equal to number of models proposed in the current iteration
          NumericVector terminal_nodes;							// create a vector called terminal_nodes.
          
          if(parent_curr_round[k]==-1){							// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
            new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a_mu);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
          }else{													// if the k+1^th element of parent_curr_round does not equal -1 
            new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
          }  
          
          terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);		// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
          updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_control_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
          NumericVector test_res;									// create vector called test_res. Length not given.
          
          if(parent_curr_round[k]==-1){							// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
            test_res=resids(_,0);								// let test_res equal the first column of resids
          }else{													// if the k+1^th element of parent_curr_round does not equal -1 
            test_res=resids(_,parent_curr_round[k]);			// let test_res equal the parent_curr_round[k]+1^th column of resids
          }
          
          NumericVector curr_test_res=updated_curr_preds[1];		// Let curr_test_res be the vector of updated predictions
          
          if(parent_curr_round[k]==-1){							// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
            curr_resids(_,0)=test_res-curr_test_res;			// let first column of curr_resids be ?the updated residuals?
          }else{													// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
            curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
          }
        }
        List temp(0);																	// list of length zero.
        
        cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
        
        for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
          List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
          if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
            cp_mat_list1=make_pelt_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),pen_mu,num_cp_mu);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
          }else{																			// If gridpoint equals 1 (or non-zero).
            cp_mat_list1=make_gridpoint_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),gridsize_mu,num_cp_mu);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
          }
          
          cp_mat_list.push_back(cp_mat_list1[0]);		// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
        }
        
      }
    }
    overall_trees=resize_bcf(overall_trees,overall_count);							// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);								// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);											// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);										// remove remaining spaces that are not filled in
  NumericVector temp_preds;													// create a vector
  List updated_preds;															// create a list
  NumericVector new_mean;														// create a vector
  NumericMatrix overall_test_preds(x_control_test.nrow(),overall_trees.size()); 	// create a matrix with no. of rows equal to that of input x_control_test, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_control_a.n_rows,overall_trees.size());					// create a matrix with no. of rows equal to that of x_control_a (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);												// reset lowest_BIC to minimum of overall_lik2
  
  for(int k=0;k<overall_trees.size();k++){																// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;																		// create a vector
    
    if(overall_parent2[k]==-1){																			// If k+1^th element of overall_parent2 equals -1
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,0),a_mu);						// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
    }else{    
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
    }  
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);													// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_control_a.n_rows);				// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_control_test,overall_trees[k],new_mean);				// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    temp_preds=updated_preds[1];																		// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																		// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;													// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);					// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);																// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);	// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_test_values=sum(M2,1);															// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}
//###################################################################################//

// [[Rcpp::export]]

List get_best_trees_mu_bcf_2(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                           double a_mu,double a_tau,double mu_mu,double mu_tau,
                           double nu,double lambda,double c,double sigma_mu_mu,double sigma_mu_tau,
                           List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                           double lowest_BIC,//int first_round,
                           IntegerVector parent,
                           List resids_cp_mat_mu,IntegerVector err_list,NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                           double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,bool is_test_data,
                           double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                           bool split_rule_node,bool gridpoint,int maxOWsize, int num_splits_mu,int num_splits_tau,int gridsize_mu, bool zero_split,
                           IntegerVector no_more_mu_trees
){
  List eval_model;										// create a list. Length not given.
  NumericVector lik_list;									// create a vector. Length not givem
  List best_subset;										// create a list. length not givn
  int overall_size=1000;									// create a variable overall_size. Initialized as 1000.
  List overall_trees(overall_size);						// create a list of length 1000.
  NumericVector overall_lik2;								// create a vector. length not given.
  IntegerVector overall_parent2;							// create a vector. length not give,
  List overall_mat(overall_size);							// create a list of length 1000.
  int overall_count=0;									// create a variable overall_count, Initialized equal to 0.
  std::vector<int> overall_parent(overall_size);			// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);			// create a vector of length 1000.
  NumericVector test_preds;								// create a vector. length not given.
  List cp_mat_list=resids_cp_mat_mu;
  
  
  
  if(zero_split==1 && no_more_mu_trees[0] !=1){
    
    //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
    
    // Rcout << "get to likelihood function mu round one. \n";
    overall_trees[0]= tree_table_mu[0];
    overall_mat[0]= tree_mat_mu[0];
    overall_parent[0]=-1;
    overall_parent2[0]=-1;
    double lik_temp=likelihood_function_bcf(resids(_,0),tree_table_mu[0],tree_mat_mu[0],a_mu,mu_mu,nu,lambda);
    double tree_prior_temp=get_tree_prior_bcf(tree_table_mu[0],tree_mat_mu[0],alpha_mu,beta_mu);
    double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(x_control_a.n_rows);
    overall_lik[0]= lowest_BIC_temp;
    overall_count=1;  
    
  }
  
  // Rcout << "get past likelihood function mu round one. \n";
  
  
  
  
  
  
  
  for(int j=0;j<num_splits_mu;j++){									// for-loop of length 5.
    int lsize=1000;										// reate a variable lsize. Initialized as 1000.
    List table_subset_curr_round(lsize);				// create a list of length 1000.
    std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000.
    List mat_subset_curr_round(lsize);					// create a list of length 1000.
    std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
    int count=0;										// create a vector count. Initialize as 0.
    for(int i=0;i<tree_table_mu.size();i++){				// create a length equal to the length of the input list tree_table_mu
      parent=-1;									// reset input vector parent to have one element, equal to -1.
      //NumericMatrix temp_list=cp_mat_list[0];		// indexes columns of splitting variables for potential splitting points. let temp_list equal the first element of the input list cp_mat_list
      best_subset=get_best_split_mu_bcf(resids(_,0),x_control_a,tree_table_mu[i],tree_mat_mu[i],
                                        a_mu,mu_mu,nu,lambda,log(c),lowest_BIC,
                                        parent[0],cp_mat_list[0],
                                                             alpha_mu,beta_mu,maxOWsize//,first_round
                                          ); // defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
      if(best_subset.size()==1){						// if get_best_split_bcf resulted in a (?no tree?) error.
        continue;									// skip to the next iteration of the iner for-loop.
      }
      List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_bcf
      List temp_mat=best_subset[6];					// list of tree matrices
      lik_list=best_subset[5];						// vector of BICs
      IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?)
      if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
        throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
      }
      if(lik_list.size()==0){							// if vector of BICs has no elements
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      if(min(lik_list)<lowest_BIC){					// If the minimum of the list of BICs is less than the input value of lowest_BIC.
        lowest_BIC=min(lik_list);					// reset lowest_BIC equal to the minimum of the list of BICs
      }
      for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
        table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
        lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
        mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
        parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
        count++;										// increment the count by 1.
        
        if(count==(lsize-1)){														// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
          lsize=lsize*2;															// Double the length
          table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
          mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
          lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
          parent_curr_round.resize(lsize);										// double length of parent_curr_round
        }
      }
    }
    table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
    mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
    lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
    parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0
      break;															// break out of the for-loop,
    }
    for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
      overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
      overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
      overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
      overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
      overall_count++;												// increment overall_count by one.
      
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    overall_trees=resize_bcf(overall_trees,overall_count);			// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);							// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);				// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);						// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_mu_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];									// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];								// list of tree tables in Occam's window
    overall_mat=eval_model[2];									// list of tree matrices in Occa's window
    overall_count=overall_trees.size();							// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];								// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){							// If more than maxOWsize models kept in Occam's window
      IntegerVector owindices=orderforOW__bcf(overall_lik2);			// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;									// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);						// create vector of length maxOWsize
      List temp_otrees(maxOWsize);							// create list of length maxOWsize
      List temp_omat(maxOWsize);								// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);					// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){							// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];			// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];			// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];			// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];		// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;									// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;								// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;									// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();						// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;							// reset overall_parent2 to keep only top maxOWsize models
    }
    
    tree_table_mu=table_subset_curr_round;								// reset tree_table_mu equal to the list of tables produced in the current iteration of the outer loop.
    IntegerVector temp1(table_subset_curr_round.size(),1);			// create a vector of ones of length equal to that of the list of tree tables for the current round.
    err_list=temp1;													// set err_list equal to temp1
    if(overall_trees.size()<overall_size-1){						// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// increse size of overall_mat
      overall_lik.resize(overall_size);							// increse size of overall_lik
      overall_parent.resize(overall_size);						// increse size of overall_parent
    }else{															// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;								// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the length of overall_mat
      overall_lik.resize(overall_size);							// double the length of overall_lik
      overall_parent.resize(overall_size);						// double the length of overall_parent
    }
    tree_mat_mu=mat_subset_curr_round;									// reset tree_mat_mu to list of treematrices obtained in current round of for-loop)
    parent=parent_curr_round;										// reset parent to parent vector obtained in current round of for-loop.
    
    if(split_rule_node==1){											// If input boolean split_rule_node = 1 (true)
      NumericVector temp_preds;									// create vector. Length not given.
      List updated_curr_preds;									// create list. Length not given.
      NumericVector new_mean;										// create vector. Length not given.
      lowest_BIC=min(overall_lik2);								// update lowest_BIC
      NumericMatrix curr_resids(resids.nrow(),resids.ncol());		// curr_resids is a matrix with dimensions of input matrix resids
      
      for(int k=0;k<table_subset_curr_round.size();k++){			// for-loop of length equal to number of models proposed in the current iteration
        NumericVector terminal_nodes;							// create a vector called terminal_nodes.
        
        if(parent_curr_round[k]==-1){							// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a_mu);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
        }else{													// if the k+1^th element of parent_curr_round does not equal -1 
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
        }  
        
        terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);		// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
        updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_control_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
        NumericVector test_res;									// create vector called test_res. Length not given.
        
        if(parent_curr_round[k]==-1){							// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          test_res=resids(_,0);								// let test_res equal the first column of resids
        }else{													// if the k+1^th element of parent_curr_round does not equal -1 
          test_res=resids(_,parent_curr_round[k]);			// let test_res equal the parent_curr_round[k]+1^th column of resids
        }
        
        NumericVector curr_test_res=updated_curr_preds[1];		// Let curr_test_res be the vector of updated predictions
        
        if(parent_curr_round[k]==-1){							// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          curr_resids(_,0)=test_res-curr_test_res;			// let first column of curr_resids be ?the updated residuals?
        }else{													// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
          curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
        }
      }
      List temp(0);																	// list of length zero.
      
      cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
      
      for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
        List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
        if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
          cp_mat_list1=make_pelt_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),pen_mu,num_cp_mu);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																			// If gridpoint equals 1 (or non-zero).
          cp_mat_list1=make_gridpoint_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),gridsize_mu,num_cp_mu);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        
        cp_mat_list.push_back(cp_mat_list1[0]);		// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
      }
      
    }
  }
  overall_trees=resize_bcf(overall_trees,overall_count);							// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);								// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);											// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);										// remove remaining spaces that are not filled in
  NumericVector temp_preds;													// create a vector
  List updated_preds;															// create a list
  NumericVector new_mean;														// create a vector
  NumericMatrix overall_test_preds(x_control_test.nrow(),overall_trees.size()); 	// create a matrix with no. of rows equal to that of input x_control_test, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_control_a.n_rows,overall_trees.size());					// create a matrix with no. of rows equal to that of x_control_a (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);												// reset lowest_BIC to minimum of overall_lik2
  
  for(int k=0;k<overall_trees.size();k++){																// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;																		// create a vector
    
    if(overall_parent2[k]==-1){																			// If k+1^th element of overall_parent2 equals -1
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,0),a_mu);						// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
    }else{    
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
    }  
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);													// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_control_a.n_rows);				// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_control_test,overall_trees[k],new_mean);				// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    temp_preds=updated_preds[1];																		// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																		// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;													// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);					// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);																// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);	// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_test_values=sum(M2,1);															// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}

//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_mu_bcf(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                           double a_mu,double a_tau,double mu_mu,double mu_tau,
                           double nu,double lambda,double c,double sigma_mu_mu,double sigma_mu_tau,
                           List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                           double lowest_BIC,//int first_round,
                           IntegerVector parent,List resids_cp_mat_mu,IntegerVector err_list,
                           NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                           double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,
                           bool is_test_data,double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                           bool split_rule_node,bool gridpoint,int maxOWsize,
                           List prev_sum_trees_mu,List prev_sum_trees_tau,List prev_sum_trees_mat_mu,List prev_sum_trees_mat_tau,
                           NumericVector y_scaled,int num_splits_mu,int num_splits_tau,int gridsize_mu,bool zero_split
){
  List example_tree_tab2 = prev_sum_trees_mu[parent[0]];
  List example_tree_mat2 = prev_sum_trees_mu[parent[0]];
  
  // Rcout << "Length of input tree table list get_best_trees_sum_mu_bcf = " << example_tree_tab2.size() << ".\n";
  // Rcout << "Length of input tree mat list get_best_trees_sum_mu_bcf= " << example_tree_mat2.size() << ".\n";
  // Rcout << "Get to Line 3625 in get_best_trees_sum_mu_bcf.\n";
  List eval_model;										// create a list
  NumericVector lik_list;								// create a vector
  List best_subset;									// create a list
  int overall_size=1000;								// create variable overall_size. Initialize equal to 1000.
  List overall_trees(overall_size);					// create a list of length 1000.
  NumericVector overall_lik2;							// create a vector.
  IntegerVector overall_parent2;						// create a vector.
  List overall_mat(overall_size);						// create a list of length 1000.
  int overall_count=0;								// create a variable initialized equal to 0.
  std::vector<int> overall_parent(overall_size);		// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);		// create a vector of length 1000.
  NumericVector test_preds;							// create a vector.
  List cp_mat_list=resids_cp_mat_mu;
  
  
  
  
  // Rcout << "Get to Line 3643 in get_best_trees_sum_mu_bcf.\n";
  
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  if(zero_split==1){
  for(int q=0; q<parent.size();q++){
  
  SEXP s_mu = prev_sum_trees_mu[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
    SEXP s_tau = prev_sum_trees_tau[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
    
    if(is<List>(s_mu)){												// is s is a list
      if(is<List>(s_tau)){												// is s is a list
        //Rcout << "RELEVANT LINE 2027.\n";
        //Rcout << "Get to Line 3310 in get_best_trees_sum_mu_bcf.\n";
        List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        
        prev_sum_trees_mu2_temp.push_back(tree_table_mu[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
        prev_sum_trees_mat_mu2_temp.push_back(tree_mat_mu[0]);	


        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,sum_trees_tau2_temp,prev_sum_trees_mat_mu2_temp,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood

        double tree_prior_temp=1;
        int p_other_mu=0;
        int p_other_tau=0;
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        // Rcout << "Get to Line 3673 in get_best_trees_sum_mu_bcf.\n";
          for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1977");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1984");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        
        // Rcout << "Get to Line 3693 in get_best_trees_sum_mu_bcf.\n";
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          
        overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
        overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
        overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
        
        overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
        
        overall_count++;												// increment overall_count by one.
        //Rcout << "Get to Line 3357 in get_best_trees_sum_mu_bcf.\n";
        
      }else{
        //Rcout << "\n RELEVANT LINE 2057.\n";
        // Rcout << "Get to Line 3707 in get_best_trees_sum_mu_bcf.\n";
        List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        prev_sum_trees_mu2_temp.push_back(tree_table_mu[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
        prev_sum_trees_mat_mu2_temp.push_back(tree_mat_mu[0]);	
        NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_tau(1);													// create list, st, of length 2.
        List st_mat_tau(1);												// create lisr, st_mat of length 2.
        st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
        st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.

        
        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,st_tau,prev_sum_trees_mat_mu2_temp,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood

        double tree_prior_temp=1;
        int p_other_mu=0;
        int p_other_tau=0;
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        // Rcout << "Get to Line 3727 in get_best_trees_sum_mu_bcf.\n";
          for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 2006");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 2013");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          //Rcout << "Get to Line 3401 in get_best_trees_sum_mu_bcf.\n";
          overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
          overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
          overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
          
          overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
          
          overall_count++;												// increment overall_count by one.
          //Rcout << "Get to Line 3409 in get_best_trees_sum_mu_bcf.\n";
        
      }
      
    }else{															// if s is not a list
      if(is<List>(s_tau)){												// is s is a list
        //Rcout << "RELEVANT LINE 2093.\n";
        // Rcout << "Get to Line 3762 in get_best_trees_sum_mu_bcf.\n";
        NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
        NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_mu(2);													// create list, st, of length 2.
        List st_mat_mu(2);												// create lisr, st_mat of length 2.
        st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
        st_mu[1]=tree_table_mu[0];												// let the first elemetn of st be sum_trees2.
        st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
        st_mat_mu[1]=tree_mat_mu[0];										// let the first element of st_mat be sum_trees_mat2.
        // return(st);
        List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        

        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,sum_trees_tau2_temp,st_mat_mu,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
        double tree_prior_temp=1;
        
        int p_other_mu=0;
        int p_other_tau=0;
        
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        //Rcout << "Get to Line 3438 in get_best_trees_sum_mu_bcf.\n";
        for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
          NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 2040");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 2047");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
        //Rcout << "Get to Line 3458 in get_best_trees_sum_mu_bcf.\n";
        overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
        overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
        overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
        
        overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
        
        overall_count++;												// increment overall_count by one.
        //Rcout << "Get to Line 3466 in get_best_trees_sum_mu_bcf.\n";
      }else{
        //Rcout << "RELEVANT LINE 2124.\n";
        //Rcout << "Get to Line 3469 in get_best_trees_sum_mu_bcf.\n";
        NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
        NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_mu(2);													// create list, st, of length 2.
        List st_mat_mu(2);												// create lisr, st_mat of length 2.
        st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
        st_mu[1]=tree_table_mu[0];										// let the second element of st be proposal_tree[0] (tree table).
        st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
        st_mat_mu[1]=tree_mat_mu[0];		
        NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_tau(1);													// create list, st, of length 2.
        List st_mat_tau(1);												// create lisr, st_mat of length 2.
        st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
        st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.
        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
        double tree_prior_temp=1;
        
        int p_other_mu=0;
        int p_other_tau=0;
        
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        
        //Rcout << "Get to Line 3493 in get_best_trees_sum_mu_bcf.\n";
        for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
          NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 2070");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 2077");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
        //Rcout << "Get to Line 3513 in get_best_trees_sum_mu_bcf.\n";
        overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
        overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
        overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
        
        overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
        //Rcout << "Get to Line 3519 in get_best_trees_sum_mu_bcf.\n";
        overall_count++;
      }
    }
    if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
      overall_size=overall_size*2;								// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
      overall_lik.resize(overall_size);							// double the size of overall_lik
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
      overall_parent.resize(overall_size);						// double the size of overall_parent
    }
  
  }
  // Rcout << "Get to Line 3878 in get_best_trees_sum_mu_bcf.\n";
  overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
  eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_mu_list, tree_lik, tree_parent
  overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
  overall_trees=eval_model[1];										// list of tree tables in Occam's window
  overall_mat=eval_model[2];											// list of tree matrices in Occam's window
  overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
  overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
  //Rcout << "Get to Line 3543 in get_best_trees_sum_mu_bcf.\n";
  //add in check to see if OW accepted more than the top maxOW models...
  if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
    //find the maxOWsize best models and continue with those!
    IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
    owindices=owindices-1;											// take one away from indices so that they are in the correct range.
    //get the top maxOWsize indices to keep in OW
    NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
    List temp_otrees(maxOWsize);									// create list of length maxOWsize
    List temp_omat(maxOWsize);										// create list of length maxOWsize
    IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
    
    //now only select those elements
    for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
      temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
      temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
      temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
      temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
    }
    
    overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
    overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
    overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
    overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
    overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
  }
  //Rcout << "Get to Line 3569 in get_best_trees_sum_mu_bcf.\n";
  if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
    overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
    overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
    overall_lik.resize(overall_size);								// increse size of overall_lik
    overall_parent.resize(overall_size);							// increse size of overall_parent
  }else{																// if length of overall_trees is greater than or equal to overall_size-1
    overall_size=2*overall_size;									// double the size
    overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
    overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
    overall_lik.resize(overall_size);								// double the length of overall_lik
    overall_parent.resize(overall_size);							// double the length of overall_parent
  }
  
  } //end if zero split=1 code
  
  // Rcout << "Get to Line 3931 in get_best_trees_sum_mu_bcf.\n";

  
  for(int j=0;j<num_splits_mu;j++){									// create a for-loop of length 5.
    // Rcout << "Get to Line 3935 in get_best_trees_sum_mu_bcf.\n";
    int lsize=1000;										// create a variable initialized equal to 1000.
    List table_subset_curr_round(lsize);				// create a list of length 1000
    std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000
    List mat_subset_curr_round(lsize);					// create a list of length 1000.
    std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
    int count=0;										// create a varialble. Initialized equal to 0.
    for(int i=0;i<tree_table_mu.size();i++){				// for loop of length equal to that of the list tree_table_mu (input or updated in previous iteration)
      //if(first_round==1){								// if input first_round equals 1.
      //  throw std::range_error("get_best_trees_sum_mu_bcf should not be used in the first round");	// throw an error.
        //parent=-1;									// reset input parent to -1
        //NumericMatrix temp_list=cp_mat_list[0];		// create matrix temp_list equal to first element of cp_mat_list
        //best_subset=get_best_split_bcf(resids(_,0),D1,tree_table_mu[i],tree_mat_mu[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,first_round);	// defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
      //}else{											// if input first_round is not equal to 1.
      
        if(err_list[i]==0){												// if i+1^th element of input vector err_list equals zero (no error?)
          NumericMatrix test_tree=tree_table_mu[i];						// create test_tree equal to i+1^th element of input list tree_table_mu
          NumericMatrix test_treemat=tree_mat_mu[i];						// create test_treemat equal to i+1^th element of input list tree_mat_mu
          NumericMatrix test_cpmat= cp_mat_list[parent[i]];			// create a matrix test_cpmat equal to the parent[i]+1^th element of cp_mat_list
          //need to append current tree_table_mu[i] to its parent sum_of_trees   
          // Rcout << "Get to Line 3954 in get_best_trees_sum_mu_bcf. i = " << i << "\n";
          // Rcout << "parent = "<< parent << "  .\n";
          // Rcout << "prev_sum_trees_mu.size() = " << prev_sum_trees_mu.size() << ".\n";
          // Rcout << "prev_sum_trees_tau.size() = " << prev_sum_trees_tau.size() << ".\n";
          // Rcout << "tree_table_mu.size() = " << tree_table_mu.size() << ".\n";
          
          List example_tree_tab = prev_sum_trees_mu[parent[i]];
          List example_tree_mat = prev_sum_trees_mu[parent[i]];
          // Rcout << " i value = " << i << ".\n";
          // Rcout << "Length of input tree table list used to get best split sum = " << example_tree_tab.size() << ".\n";
          // Rcout << "Length of input tree mat list used to get best split sum= " << example_tree_mat.size() << ".\n";
          
          best_subset=get_best_split_sum_mu_bcf(resids(_,parent[i]),x_control_a,tree_table_mu[i],tree_mat_mu[i],
                                            a_mu,a_tau,mu_mu,mu_tau,nu,lambda,log(c),
                                            lowest_BIC,parent[i],cp_mat_list[parent[i]],
                                            alpha_mu,beta_mu,alpha_tau,beta_tau,
                                            maxOWsize,//first_round,
                                            prev_sum_trees_mu,prev_sum_trees_tau,prev_sum_trees_mat_mu,prev_sum_trees_mat_tau,
                                            y_scaled,parent,i,z);	// defined on line 1074. Returns the BIC, best splitting variable (column number), value of covariate for splitting variable, list including tree table and tree matrix, list of tree tables, list of BICs, list of tree matrices, tree parent vector.
          
          // Rcout << "Get to Line 3968 in get_best_trees_sum_mu_bcf.\n";
          
          // return(best_subset);
        }else if(err_list[i]==1){		// if i+1^th element of input vector err_list equals one (no error?).
          continue;						// skip to next iteration of for-loop.
        }else{							// if i+1^th element of input vector err_list is not 0 or 1 ???.
          List ret_list(6);				// list of length 6.
          ret_list[0]=9999;				// first element of list equals 9999.
          ret_list[1]=err_list[i];		// second element is i+1^th element of input vector err_list. (some number not equal to 0 or 1).
          ret_list[2]=i;					// third element of list is number indexing the element of tree_table_mu (index starts from 0).
          ret_list[3]=j;					// fourth element is between 0 and 4. Index of outer for-loop.
          ret_list[4]=tree_table_mu;			// fifth element of list is the list of tree tables.
          ret_list[5]=err_list;			// sixth element of list is vector err_list.
          return(ret_list);				// return the list.
          throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");	// throw an error.
        }
      //} tget_best_trees_sum_mu_bcf is never used in the first round
      
      
      if(best_subset.size()==1){	// if get_best_split_bcf resulted in a (?no ?sum-of? tree?) error.
        continue;				// skip to the next iteration of the inner for-loop.
      }
      // Rcout << "Get to Line 3986 in get_best_trees_sum_mu_bcf.\n";
      List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_sum
      List temp_mat=best_subset[6];					// list of tree matrices
      lik_list=best_subset[5];						// vector of BICs
      IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?)
      if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
        throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
      }
      if(lik_list.size()==0){							// if vector of BICs has no elements
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){						// If the minimum of the list of BICs is less than the input value of lowest_BIC.
        lowest_BIC=min(lik_list);						// reset lowest_BIC equal to the minimum of the list of BICs
      }
      for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
        table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
        lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
        mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
        parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
        count++;										// increment the count by 1.
        // Rcout << "Get to Line 4007 in get_best_trees_sum_mu_bcf.\n";
        if(count==(lsize-1)){							// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
          lsize=lsize*2;								// Double the length
          table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
          mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
          lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
          parent_curr_round.resize(lsize);										// double length of parent_curr_round
        }
      }
    }
    table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
    mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
    lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
    parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
    // Rcout << "Get to Line 4021 in get_best_trees_sum_mu_bcf.\n";
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0 
      break;															// break out of the for-loop,
    }
    for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
      overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
      overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
      overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
      overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
      overall_count++;												// increment overall_count by one.
      
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    // Rcout << "Get to Line 4040 in get_best_trees_sum_mu_bcf.\n";
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_mu_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      // Rcout << "Get to Line 4061 in get_best_trees_sum_mu_bcf.\n";
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      // Rcout << "Get to Line 4069 in get_best_trees_sum_mu_bcf.\n";
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    // Rcout << "Get to Line 4076 in get_best_trees_sum_mu_bcf.\n";
    tree_table_mu=table_subset_curr_round;									// reset tree_table_mu equal to the list of tables produced in the current iteration of the outer loop.
    IntegerVector temp1(table_subset_curr_round.size(),1);				// create a vector of ones of length equal to that of the list of tree tables for the current round.
    err_list=temp1;														// set err_list equal to temp1
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    tree_mat_mu=mat_subset_curr_round;										// reset tree_mat_mu to list of treematrices obtained in current round of for-loop)
    parent=parent_curr_round;											// reset parent to parent vector obtained in current round of for-loop.
    // Rcout << "Get to Line 4094 in get_best_trees_sum_mu_bcf.\n";
    if(split_rule_node==1){												// If input boolean split_rule_node = 1 (true)
      NumericVector temp_preds;										// create vector. Length not given.
      List updated_curr_preds;										// create list. Length not given.
      NumericVector new_mean;											// create vector. Length not given.
      lowest_BIC=min(overall_lik2);									// update lowest_BIC
      NumericMatrix curr_resids(resids.nrow(),resids.ncol());			// curr_resids is a matrix with dimensions of input matrix resids
      // Rcout << "Get to Line 4101 in get_best_trees_sum_mu_bcf.\n";
      for(int k=0;k<table_subset_curr_round.size();k++){				// for-loop of length equal to number of models proposed in the current iteration
        NumericVector terminal_nodes;								// create a vector called terminal_nodes.
        
        if(parent_curr_round[k]==-1){								// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a_mu);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
        }else{														// if the k+1^th element of parent_curr_round does not equal -1 
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
        }  
        // Rcout << "Get to Line 4110 in get_best_trees_sum_mu_bcf.\n";
        terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);			// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
        updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_control_a.n_rows);	// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
        NumericVector test_res;												// create vector called test_res. Length not given.
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          test_res=resids(_,0);											// let test_res equal the first column of resids
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 
          test_res=resids(_,parent_curr_round[k]);						// let test_res equal the parent_curr_round[k]+1^th column of resids
        }
        
        NumericVector curr_test_res=updated_curr_preds[1];					// Let curr_test_res be the vector of updated predictions
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          curr_resids(_,0)=test_res-curr_test_res;						// let first column of curr_resids be ?the updated residuals?
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
          curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
        }
      }
      List temp(0);																	// list of length zero.
      // Rcout << "Get to Line 4130 in get_best_trees_sum_mu_bcf.\n";
      cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
      
      for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
        List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
        if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
          cp_mat_list1=make_pelt_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),pen_mu,num_cp_mu);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																			// If gridpoint equals 1 (or non-zero).
          cp_mat_list1=make_gridpoint_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),gridsize_mu,num_cp_mu);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        
        cp_mat_list.push_back(cp_mat_list1[0]);											// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
      }
      
    }
  }
  // Rcout << "Get to Line 4146 in get_best_trees_sum_mu_bcf.\n";
  overall_trees=resize_bcf(overall_trees,overall_count);		// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);			// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);					// remove remaining spaces that are not filled in
  NumericVector temp_preds;								// create a vector
  List updated_preds;										// create a list
  NumericVector new_mean;									// create a vector
  NumericMatrix overall_test_preds(x_control_test.nrow(),overall_trees.size());	// create a matrix with no. of rows equal to that of input test_data, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_control_a.n_rows,overall_trees.size());		// create a matrix with no. of rows equal to that of D1 (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);							// reset lowest_BIC to minimum of overall_lik2
  // Rcout << "Get to Line 4157 in get_best_trees_sum_mu_bcf.\n";
  for(int k=0;k<overall_trees.size();k++){			// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;					// create a vector
    
    if(overall_parent2[k]==-1){						// If k+1^th element of overall_parent2 equals -1
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,0),a_mu);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
    }else{    
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
    }  
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);		// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_control_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    //Rcout << "Get to Line 3824 in get_best_trees_sum_mu_bcf.\n";
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_control_test,overall_trees[k],new_mean);		// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    //Rcout << "Get to Line 3825 in get_best_trees_sum_mu_bcf.\n";
    temp_preds=updated_preds[1];																// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;											// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
    //Rcout << "Get to Line 3829 in get_best_trees_sum_mu_bcf.\n";
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);			// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);														// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  // Rcout << "Get to Line 4179 in get_best_trees_sum_mu_bcf.\n";
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);		// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  // Rcout << "Get to Line 4181 in get_best_trees_sum_mu_bcf.\n";
  arma::colvec predicted_test_values=sum(M2,1);		// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_mu_bcf_2(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                               double a_mu,double a_tau,double mu_mu,double mu_tau,
                               double nu,double lambda,double c,double sigma_mu_mu,double sigma_mu_tau,
                               List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                               double lowest_BIC,//int first_round,
                               IntegerVector parent,List resids_cp_mat_mu,IntegerVector err_list,
                               NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                               double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,
                               bool is_test_data,double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                               bool split_rule_node,bool gridpoint,int maxOWsize,
                               List prev_sum_trees_mu,List prev_sum_trees_tau,List prev_sum_trees_mat_mu,List prev_sum_trees_mat_tau,
                               NumericVector y_scaled,int num_splits_mu,int num_splits_tau,int gridsize_mu,bool zero_split,
                               IntegerVector no_more_mu_trees
){
  List example_tree_tab2 = prev_sum_trees_mu[parent[0]];
  List example_tree_mat2 = prev_sum_trees_mu[parent[0]];
  
  // Rcout << "Length of input tree table list get_best_trees_sum_mu_bcf = " << example_tree_tab2.size() << ".\n";
  // Rcout << "Length of input tree mat list get_best_trees_sum_mu_bcf= " << example_tree_mat2.size() << ".\n";
  // Rcout << "Get to Line 3625 in get_best_trees_sum_mu_bcf.\n";
  List eval_model;										// create a list
  NumericVector lik_list;								// create a vector
  List best_subset;									// create a list
  int overall_size=1000;								// create variable overall_size. Initialize equal to 1000.
  List overall_trees(overall_size);					// create a list of length 1000.
  NumericVector overall_lik2;							// create a vector.
  IntegerVector overall_parent2;						// create a vector.
  List overall_mat(overall_size);						// create a list of length 1000.
  int overall_count=0;								// create a variable initialized equal to 0.
  std::vector<int> overall_parent(overall_size);		// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);		// create a vector of length 1000.
  NumericVector test_preds;							// create a vector.
  List cp_mat_list=resids_cp_mat_mu;
  
  
  
  
   // Rcout << "Get to Line 4525 in get_best_trees_sum_mu_bcf.\n";
  
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  if(zero_split==1){
    if( is_true(all(no_more_mu_trees==1)) ){
      
    }else{
      
    for(int q=0; q<parent.size();q++){
      // Rcout << "Line 4531.\n";
      if(no_more_mu_trees[q]!=1){
      // Rcout << "Line 4533.\n";
      SEXP s_mu = prev_sum_trees_mu[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
      SEXP s_tau = prev_sum_trees_tau[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
      
      if(is<List>(s_mu)){												// is s is a list
        if(is<List>(s_tau)){												// is s is a list
          //Rcout << "RELEVANT LINE 2027.\n";
          //Rcout << "Get to Line 3310 in get_best_trees_sum_mu_bcf.\n";
          List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          
          prev_sum_trees_mu2_temp.push_back(tree_table_mu[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
          prev_sum_trees_mat_mu2_temp.push_back(tree_mat_mu[0]);	
          
          
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,sum_trees_tau2_temp,prev_sum_trees_mat_mu2_temp,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          
          double tree_prior_temp=1;
          int p_other_mu=0;
          int p_other_tau=0;
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          // Rcout << "Get to Line 3673 in get_best_trees_sum_mu_bcf.\n";
          for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1977");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1984");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          
          // Rcout << "Get to Line 3693 in get_best_trees_sum_mu_bcf.\n";
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          
          overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
          overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
          overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
          
          overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
          
          overall_count++;												// increment overall_count by one.
          //Rcout << "Get to Line 3357 in get_best_trees_sum_mu_bcf.\n";
          
        }else{
          //Rcout << "\n RELEVANT LINE 2057.\n";
          // Rcout << "Get to Line 3707 in get_best_trees_sum_mu_bcf.\n";
          List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          prev_sum_trees_mu2_temp.push_back(tree_table_mu[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
          prev_sum_trees_mat_mu2_temp.push_back(tree_mat_mu[0]);	
          NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_tau(1);													// create list, st, of length 2.
          List st_mat_tau(1);												// create lisr, st_mat of length 2.
          st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
          st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.
          
          
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,st_tau,prev_sum_trees_mat_mu2_temp,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          
          double tree_prior_temp=1;
          int p_other_mu=0;
          int p_other_tau=0;
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          // Rcout << "Get to Line 3727 in get_best_trees_sum_mu_bcf.\n";
          for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 2006");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 2013");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          //Rcout << "Get to Line 3401 in get_best_trees_sum_mu_bcf.\n";
          overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
          overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
          overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
          
          overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
          
          overall_count++;												// increment overall_count by one.
          //Rcout << "Get to Line 3409 in get_best_trees_sum_mu_bcf.\n";
          
        }
        
      }else{															// if s is not a list
        if(is<List>(s_tau)){												// is s is a list
          //Rcout << "RELEVANT LINE 2093.\n";
          // Rcout << "Get to Line 3762 in get_best_trees_sum_mu_bcf.\n";
          NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
          NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_mu(2);													// create list, st, of length 2.
          List st_mat_mu(2);												// create lisr, st_mat of length 2.
          st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
          st_mu[1]=tree_table_mu[0];												// let the first elemetn of st be sum_trees2.
          st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
          st_mat_mu[1]=tree_mat_mu[0];										// let the first element of st_mat be sum_trees_mat2.
          // return(st);
          List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          
          
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,sum_trees_tau2_temp,st_mat_mu,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          double tree_prior_temp=1;
          
          int p_other_mu=0;
          int p_other_tau=0;
          
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          //Rcout << "Get to Line 3438 in get_best_trees_sum_mu_bcf.\n";
          for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
            NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 2040");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 2047");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          //Rcout << "Get to Line 3458 in get_best_trees_sum_mu_bcf.\n";
          overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
          overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
          overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
          
          overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
          
          overall_count++;												// increment overall_count by one.
          //Rcout << "Get to Line 3466 in get_best_trees_sum_mu_bcf.\n";
        }else{
          //Rcout << "RELEVANT LINE 2124.\n";
          //Rcout << "Get to Line 3469 in get_best_trees_sum_mu_bcf.\n";
          NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
          NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_mu(2);													// create list, st, of length 2.
          List st_mat_mu(2);												// create lisr, st_mat of length 2.
          st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
          st_mu[1]=tree_table_mu[0];										// let the second element of st be proposal_tree[0] (tree table).
          st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
          st_mat_mu[1]=tree_mat_mu[0];		
          NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_tau(1);													// create list, st, of length 2.
          List st_mat_tau(1);												// create lisr, st_mat of length 2.
          st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
          st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          double tree_prior_temp=1;
          
          int p_other_mu=0;
          int p_other_tau=0;
          
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          
          //Rcout << "Get to Line 3493 in get_best_trees_sum_mu_bcf.\n";
          for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
            NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 2070");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 2077");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_control_a.n_rows);			// x_control_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          //Rcout << "Get to Line 3513 in get_best_trees_sum_mu_bcf.\n";
          overall_trees[overall_count]=tree_table_mu[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
          overall_mat[overall_count]=tree_mat_mu[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
          overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
          
          overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
          //Rcout << "Get to Line 3519 in get_best_trees_sum_mu_bcf.\n";
          overall_count++;
        }
      }
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    }
    // Rcout << "Line 4762.\n";
    
     // Rcout << "Get to Line 4764 in get_best_trees_sum_mu_bcf.\n";
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_mu_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //Rcout << "Get to Line 3543 in get_best_trees_sum_mu_bcf.\n";
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    //Rcout << "Get to Line 3569 in get_best_trees_sum_mu_bcf.\n";
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    // Rcout << "Line 4814.\n";
    }
  } //end if zero split=1 code
  
   // Rcout << "Get to Line 4818 in get_best_trees_sum_mu_bcf.\n";
  
  
  for(int j=0;j<num_splits_mu;j++){									// create a for-loop of length 5.
    // Rcout << "Get to Line 3935 in get_best_trees_sum_mu_bcf.\n";
    int lsize=1000;										// create a variable initialized equal to 1000.
    List table_subset_curr_round(lsize);				// create a list of length 1000
    std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000
    List mat_subset_curr_round(lsize);					// create a list of length 1000.
    std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
    int count=0;										// create a varialble. Initialized equal to 0.
    for(int i=0;i<tree_table_mu.size();i++){				// for loop of length equal to that of the list tree_table_mu (input or updated in previous iteration)
      //if(first_round==1){								// if input first_round equals 1.
      //  throw std::range_error("get_best_trees_sum_mu_bcf should not be used in the first round");	// throw an error.
        //parent=-1;									// reset input parent to -1
        //NumericMatrix temp_list=cp_mat_list[0];		// create matrix temp_list equal to first element of cp_mat_list
        //best_subset=get_best_split_bcf(resids(_,0),D1,tree_table_mu[i],tree_mat_mu[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,first_round);	// defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
      //}else{											// if input first_round is not equal to 1.
        if(err_list[i]==0){												// if i+1^th element of input vector err_list equals zero (no error?)
          NumericMatrix test_tree=tree_table_mu[i];						// create test_tree equal to i+1^th element of input list tree_table_mu
          NumericMatrix test_treemat=tree_mat_mu[i];						// create test_treemat equal to i+1^th element of input list tree_mat_mu
          NumericMatrix test_cpmat= cp_mat_list[parent[i]];			// create a matrix test_cpmat equal to the parent[i]+1^th element of cp_mat_list
          //need to append current tree_table_mu[i] to its parent sum_of_trees   
          // Rcout << "Get to Line 3954 in get_best_trees_sum_mu_bcf. i = " << i << "\n";
          // Rcout << "parent = "<< parent << "  .\n";
          // Rcout << "prev_sum_trees_mu.size() = " << prev_sum_trees_mu.size() << ".\n";
          // Rcout << "prev_sum_trees_tau.size() = " << prev_sum_trees_tau.size() << ".\n";
          // Rcout << "tree_table_mu.size() = " << tree_table_mu.size() << ".\n";
          
          List example_tree_tab = prev_sum_trees_mu[parent[i]];
          List example_tree_mat = prev_sum_trees_mu[parent[i]];
          // Rcout << " i value = " << i << ".\n";
          // Rcout << "Length of input tree table list used to get best split sum = " << example_tree_tab.size() << ".\n";
          // Rcout << "Length of input tree mat list used to get best split sum= " << example_tree_mat.size() << ".\n";
          
          best_subset=get_best_split_sum_mu_bcf(resids(_,parent[i]),x_control_a,tree_table_mu[i],tree_mat_mu[i],
                                                a_mu,a_tau,mu_mu,mu_tau,nu,lambda,log(c),
                                                lowest_BIC,parent[i],cp_mat_list[parent[i]],
                                                alpha_mu,beta_mu,alpha_tau,beta_tau,
                                                maxOWsize,//first_round,
                                                prev_sum_trees_mu,prev_sum_trees_tau,prev_sum_trees_mat_mu,prev_sum_trees_mat_tau,
                                                y_scaled,parent,i,z);	// defined on line 1074. Returns the BIC, best splitting variable (column number), value of covariate for splitting variable, list including tree table and tree matrix, list of tree tables, list of BICs, list of tree matrices, tree parent vector.
          
          // Rcout << "Get to Line 3968 in get_best_trees_sum_mu_bcf.\n";
          
          // return(best_subset);
        }else if(err_list[i]==1){		// if i+1^th element of input vector err_list equals one (no error?).
          continue;						// skip to next iteration of for-loop.
        }else{							// if i+1^th element of input vector err_list is not 0 or 1 ???.
          List ret_list(6);				// list of length 6.
          ret_list[0]=9999;				// first element of list equals 9999.
          ret_list[1]=err_list[i];		// second element is i+1^th element of input vector err_list. (some number not equal to 0 or 1).
          ret_list[2]=i;					// third element of list is number indexing the element of tree_table_mu (index starts from 0).
          ret_list[3]=j;					// fourth element is between 0 and 4. Index of outer for-loop.
          ret_list[4]=tree_table_mu;			// fifth element of list is the list of tree tables.
          ret_list[5]=err_list;			// sixth element of list is vector err_list.
          return(ret_list);				// return the list.
          throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");	// throw an error.
        }
      //}   get_best_split_sum_mu_bcf_2 is never used in the first round 
      
      
      if(best_subset.size()==1){	// if get_best_split_bcf resulted in a (?no ?sum-of? tree?) error.
        continue;				// skip to the next iteration of the inner for-loop.
      }
      // Rcout << "Get to Line 3986 in get_best_trees_sum_mu_bcf.\n";
      List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_sum
      List temp_mat=best_subset[6];					// list of tree matrices
      lik_list=best_subset[5];						// vector of BICs
      IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?)
      if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
        throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
      }
      if(lik_list.size()==0){							// if vector of BICs has no elements
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){						// If the minimum of the list of BICs is less than the input value of lowest_BIC.
        lowest_BIC=min(lik_list);						// reset lowest_BIC equal to the minimum of the list of BICs
      }
      for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
        table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
        lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
        mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
        parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
        count++;										// increment the count by 1.
        // Rcout << "Get to Line 4007 in get_best_trees_sum_mu_bcf.\n";
        if(count==(lsize-1)){							// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
          lsize=lsize*2;								// Double the length
          table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
          mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
          lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
          parent_curr_round.resize(lsize);										// double length of parent_curr_round
        }
      }
    }
    table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
    mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
    lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
    parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
    // Rcout << "Get to Line 4021 in get_best_trees_sum_mu_bcf.\n";
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0 
      break;															// break out of the for-loop,
    }
    for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
      overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
      overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
      overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
      overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
      overall_count++;												// increment overall_count by one.
      
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    // Rcout << "Get to Line 4040 in get_best_trees_sum_mu_bcf.\n";
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_mu_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      // Rcout << "Get to Line 4061 in get_best_trees_sum_mu_bcf.\n";
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      // Rcout << "Get to Line 4069 in get_best_trees_sum_mu_bcf.\n";
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    // Rcout << "Get to Line 4076 in get_best_trees_sum_mu_bcf.\n";
    tree_table_mu=table_subset_curr_round;									// reset tree_table_mu equal to the list of tables produced in the current iteration of the outer loop.
    IntegerVector temp1(table_subset_curr_round.size(),1);				// create a vector of ones of length equal to that of the list of tree tables for the current round.
    err_list=temp1;														// set err_list equal to temp1
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    tree_mat_mu=mat_subset_curr_round;										// reset tree_mat_mu to list of treematrices obtained in current round of for-loop)
    parent=parent_curr_round;											// reset parent to parent vector obtained in current round of for-loop.
    // Rcout << "Get to Line 4094 in get_best_trees_sum_mu_bcf.\n";
    if(split_rule_node==1){												// If input boolean split_rule_node = 1 (true)
      NumericVector temp_preds;										// create vector. Length not given.
      List updated_curr_preds;										// create list. Length not given.
      NumericVector new_mean;											// create vector. Length not given.
      lowest_BIC=min(overall_lik2);									// update lowest_BIC
      NumericMatrix curr_resids(resids.nrow(),resids.ncol());			// curr_resids is a matrix with dimensions of input matrix resids
      // Rcout << "Get to Line 4101 in get_best_trees_sum_mu_bcf.\n";
      for(int k=0;k<table_subset_curr_round.size();k++){				// for-loop of length equal to number of models proposed in the current iteration
        NumericVector terminal_nodes;								// create a vector called terminal_nodes.
        
        if(parent_curr_round[k]==-1){								// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,0),a_mu);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
        }else{														// if the k+1^th element of parent_curr_round does not equal -1 
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],resids(_,parent_curr_round[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
        }  
        // Rcout << "Get to Line 4110 in get_best_trees_sum_mu_bcf.\n";
        terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);			// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
        updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_control_a.n_rows);	// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
        NumericVector test_res;												// create vector called test_res. Length not given.
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          test_res=resids(_,0);											// let test_res equal the first column of resids
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 
          test_res=resids(_,parent_curr_round[k]);						// let test_res equal the parent_curr_round[k]+1^th column of resids
        }
        
        NumericVector curr_test_res=updated_curr_preds[1];					// Let curr_test_res be the vector of updated predictions
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          curr_resids(_,0)=test_res-curr_test_res;						// let first column of curr_resids be ?the updated residuals?
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
          curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
        }
      }
      List temp(0);																	// list of length zero.
      // Rcout << "Get to Line 4130 in get_best_trees_sum_mu_bcf.\n";
      cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
      
      for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
        List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
        if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
          cp_mat_list1=make_pelt_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),pen_mu,num_cp_mu);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																			// If gridpoint equals 1 (or non-zero).
          cp_mat_list1=make_gridpoint_cpmat_mu_bcf(wrap(x_control_a),curr_resids(_,f),gridsize_mu,num_cp_mu);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        
        cp_mat_list.push_back(cp_mat_list1[0]);											// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
      }
      
    }
  }
  // Rcout << "Get to Line 4146 in get_best_trees_sum_mu_bcf.\n";
  overall_trees=resize_bcf(overall_trees,overall_count);		// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);			// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);					// remove remaining spaces that are not filled in
  NumericVector temp_preds;								// create a vector
  List updated_preds;										// create a list
  NumericVector new_mean;									// create a vector
  NumericMatrix overall_test_preds(x_control_test.nrow(),overall_trees.size());	// create a matrix with no. of rows equal to that of input test_data, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_control_a.n_rows,overall_trees.size());		// create a matrix with no. of rows equal to that of D1 (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);							// reset lowest_BIC to minimum of overall_lik2
  // Rcout << "Get to Line 4157 in get_best_trees_sum_mu_bcf.\n";
  for(int k=0;k<overall_trees.size();k++){			// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;					// create a vector
    
    if(overall_parent2[k]==-1){						// If k+1^th element of overall_parent2 equals -1
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,0),a_mu);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
    }else{    
      new_mean=update_mean_var_bcf(overall_trees[k],overall_mat[k],resids(_,overall_parent2[k]),a_mu);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
    }  
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);		// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_control_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    //Rcout << "Get to Line 3824 in get_best_trees_sum_mu_bcf.\n";
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_control_test,overall_trees[k],new_mean);		// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    //Rcout << "Get to Line 3825 in get_best_trees_sum_mu_bcf.\n";
    temp_preds=updated_preds[1];																// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;											// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
    //Rcout << "Get to Line 3829 in get_best_trees_sum_mu_bcf.\n";
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);			// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);														// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  // Rcout << "Get to Line 4179 in get_best_trees_sum_mu_bcf.\n";
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);		// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  // Rcout << "Get to Line 4181 in get_best_trees_sum_mu_bcf.\n";
  arma::colvec predicted_test_values=sum(M2,1);		// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_tau_bcf(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                                       double a_mu,double a_tau,double mu_mu,double mu_tau,double nu,double lambda,double c,
                                       double sigma_mu_mu,double sigma_mu_tau,List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                                       double lowest_BIC,//int first_round,
                                       IntegerVector parent,List resids_cp_mat_tau,IntegerVector err_list,
                                       NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                                       double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,bool is_test_data,
                                       double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                                       bool split_rule_node,bool gridpoint,int maxOWsize,
                                       int num_splits_mu,int num_splits_tau,int gridsize_tau,bool zero_split,IntegerVector no_more_tau_trees
){
  List eval_model;										// create a list
  NumericVector lik_list;								// create a vector
  List best_subset;									// create a list
  int overall_size=1000;								// create variable overall_size. Initialize equal to 1000.
  List overall_trees(overall_size);					// create a list of length 1000.
  NumericVector overall_lik2;							// create a vector.
  IntegerVector overall_parent2;						// create a vector.
  List overall_mat(overall_size);						// create a list of length 1000.
  int overall_count=0;								// create a variable initialized equal to 0.
  std::vector<int> overall_parent(overall_size);		// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);		// create a vector of length 1000.
  NumericVector test_preds;							// create a vector.
  List cp_mat_list=resids_cp_mat_tau;
  
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  arma::mat temp_all_resids_a(resids.begin(), resids.nrow(), resids.ncol(), false);	
  arma::mat temp_treated_resids_a = temp_all_resids_a.rows(arma::find(z_ar==1));	
  NumericMatrix temp_treated_resids = wrap(temp_treated_resids_a);	
  

  if(zero_split==1 && no_more_tau_trees[0] != 1){
    
    //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
    
     //Rcout << "get to likelihood function tau round one. \n";
    overall_trees[0]= tree_table_tau[0];
    overall_mat[0]= tree_mat_tau[0];
    overall_parent[0]=-1;
    overall_parent2[0]=-1;
    
    double lik_temp=sumtree_likelihood_tau_round1_bcf(resids(_,0),tree_table_tau[0],tree_mat_tau[0],resids(_,0).size(),a_mu,a_tau,nu,lambda,z);
    double tree_prior_temp=get_tree_prior_bcf(tree_table_tau[0],tree_mat_tau[0],alpha_tau,beta_tau);
    double lowest_BIC_temp=-2*(lik_temp+log(tree_prior_temp))+1*log(x_control_a.n_rows);
    overall_lik[0]= lowest_BIC_temp;
    overall_count=1;  
    
  }
  
  //Rcout << "Line 4263.\n";
  
  
  for(int j=0;j<num_splits_tau;j++){									// create a for-loop of length 5.
    int lsize=1000;										// create a variable initialized equal to 1000.
    List table_subset_curr_round(lsize);				// create a list of length 1000
    std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000
    List mat_subset_curr_round(lsize);					// create a list of length 1000.
    std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
    int count=0;										// create a varialble. Initialized equal to 0.
    for(int i=0;i<tree_table_tau.size();i++){				// for loop of length equal tothat of the list tree_table (input or updated in previous iteration)
      parent=-1;									// reset input vector parent to have one element, equal to -1.
      best_subset=get_best_split_tau_bcf(resids(_,0),x_moderate_a,
                                         tree_table_tau[i],tree_mat_tau[i],
                                         a_mu, a_tau,mu_tau,nu,lambda,
                                         log(c),
                                         lowest_BIC,
                                         parent[0],
                                         cp_mat_list[0],
                                         alpha_mu,beta_mu,alpha_tau,beta_tau,
                                         maxOWsize,//first_round,
                                         z); // defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
      if(best_subset.size()==1){						// if get_best_split_bcf resulted in a (?no tree?) error.
        continue;									// skip to the next iteration of the iner for-loop.
      }
      List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_bcf
      List temp_mat=best_subset[6];					// list of tree matrices
      lik_list=best_subset[5];						// vector of BICs
      IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?)
      //Rcout << "Line 4293.\n";
      
      if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
        throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
      }
      if(lik_list.size()==0){							// if vector of BICs has no elements
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      if(min(lik_list)<lowest_BIC){					// If the minimum of the list of BICs is less than the input value of lowest_BIC.
        lowest_BIC=min(lik_list);					// reset lowest_BIC equal to the minimum of the list of BICs
      }
      for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
        table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
        lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
        mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
        parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
        count++;										// increment the count by 1.
        
        if(count==(lsize-1)){														// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
          lsize=lsize*2;															// Double the length
          table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
          mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
          lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
          parent_curr_round.resize(lsize);										// double length of parent_curr_round
        }
      }
    }
    //Rcout << "Line 4320.\n";
    
    table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
    mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
    lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
    parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0 
      break;															// break out of the for-loop,
    }
    for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
      overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
      overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
      overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
      overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
      overall_count++;												// increment overall_count by one.
      
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    //Rcout << "Line 4380.\n";
    
    tree_table_tau=table_subset_curr_round;									// reset tree_table_tau equal to the list of tables produced in the current iteration of the outer loop.
    IntegerVector temp1(table_subset_curr_round.size(),1);				// create a vector of ones of length equal to that of the list of tree tables for the current round.
    err_list=temp1;														// set err_list equal to temp1
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    tree_mat_tau=mat_subset_curr_round;										// reset tree_mat to list of treematrices obtained in current round of for-loop)
    parent=parent_curr_round;											// reset parent to parent vector obtained in current round of for-loop.
    
    if(split_rule_node==1){												// If input boolean split_rule_node = 1 (true)
      NumericVector temp_preds;										// create vector. Length not given.
      List updated_curr_preds;										// create list. Length not given.
      NumericVector new_mean;											// create vector. Length not given.
      lowest_BIC=min(overall_lik2);									// update lowest_BIC
      NumericMatrix curr_resids(resids.nrow(),resids.ncol());			// curr_resids is a matrix with dimensions of input matrix resids
      
      for(int k=0;k<table_subset_curr_round.size();k++){				// for-loop of length equal to number of models proposed in the current iteration
        NumericVector terminal_nodes;								// create a vector called terminal_nodes.
        
        NumericMatrix temp_all_mat = mat_subset_curr_round[k];
        
        arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
        arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
        NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);
        
        if(parent_curr_round[k]==-1){								// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
        }else{														// if the k+1^th element of parent_curr_round does not equal -1 
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,parent_curr_round[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
        }  
        
        terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);			// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
        updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_moderate_a.n_rows);	// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
        NumericVector test_res;												// create vector called test_res. Length not given.
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          test_res=resids(_,0);											// let test_res equal the first column of resids
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 
          test_res=resids(_,parent_curr_round[k]);						// let test_res equal the parent_curr_round[k]+1^th column of resids
        }
        //Rcout << "Line 4431.\n";
        
        NumericVector temp_preds = updated_curr_preds[1];
        NumericVector curr_test_res=z*temp_preds;				//might be unnecessary to multiply by z because will probably subset correctly when this vector is used	// Let curr_test_res be the vector of updated predictions
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          curr_resids(_,0)=test_res-curr_test_res;						// let first column of curr_resids be ?the updated residuals?
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
          curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
        }
      }
      List temp(0);																	// list of length zero.
      
      cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
      
      for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
        List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
        if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
          cp_mat_list1=make_pelt_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),pen_tau,num_cp_tau,z);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																			// If gridpoint equals 1 (or non-zero).
          cp_mat_list1=make_gridpoint_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),gridsize_tau,num_cp_tau,z);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        
        cp_mat_list.push_back(cp_mat_list1[0]);											// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
      }
      
    }
  }
  overall_trees=resize_bcf(overall_trees,overall_count);		// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);			// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);					// remove remaining spaces that are not filled in
  NumericVector temp_preds;								// create a vector
  List updated_preds;										// create a list
  NumericVector new_mean;									// create a vector
  NumericMatrix overall_test_preds(x_moderate_test.nrow(),overall_trees.size());	// create a matrix with no. of rows equal to that of input x_moderate_test, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_moderate_a.n_rows,overall_trees.size());		// create a matrix with no. of rows equal to that of x_moderate_a (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);							// reset lowest_BIC to minimum of overall_lik2
  //Rcout << "Line 4469.\n";
  
  
  for(int k=0;k<overall_trees.size();k++){			// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;					// create a vector
    
    NumericMatrix temp_all_mat = overall_mat[k];
    //Rcout << "Line 4476.\n";
    
    arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
    arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
    NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);	
    //Rcout << "Line 4480. k = " << k << ".\n";
    //Rcout << "overall_parent2[k] = " << overall_parent2[k] << ".\n";
    if(overall_parent2[k]==-1){						// If k+1^th element of overall_parent2 equals -1
      //Rcout << "Line 4484.\n";
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
      //Rcout << "Line 4486.\n";
      
    }else{    
      //Rcout << "Line 4489.\n";
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,overall_parent2[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
      //Rcout << "Line 4491.\n";
    }  
    //Rcout << "Line 4486.\n";
    
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);		// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_moderate_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_moderate_test,overall_trees[k],new_mean);		// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    temp_preds=updated_preds[1];																// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;											// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
  }
  //Rcout << "Line 4494.\n";
  
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);			// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);														// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);		// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_test_values=sum(M2,1);		// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_tau_round1_bcf(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                                   double a_mu,double a_tau,double mu_mu,double mu_tau,double nu,double lambda,double c,
                                   double sigma_mu_mu,double sigma_mu_tau,List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                                   double lowest_BIC,//int first_round,
                                   IntegerVector parent,List resids_cp_mat_tau,IntegerVector err_list,
                                   NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                                   double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,bool is_test_data,
                                   double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                                   bool split_rule_node,bool gridpoint,int maxOWsize,
                                   List prev_sum_trees_mu,
                                   List prev_sum_trees_mat_mu,
                                   NumericVector y_scaled,int num_splits_mu,int num_splits_tau,int gridsize_tau,bool zero_split
){
  List eval_model;										// create a list
  NumericVector lik_list;								// create a vector
  List best_subset;									// create a list
  int overall_size=1000;								// create variable overall_size. Initialize equal to 1000.
  List overall_trees(overall_size);					// create a list of length 1000.
  NumericVector overall_lik2;							// create a vector.
  IntegerVector overall_parent2;						// create a vector.
  List overall_mat(overall_size);						// create a list of length 1000.
  int overall_count=0;								// create a variable initialized equal to 0.
  std::vector<int> overall_parent(overall_size);		// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);		// create a vector of length 1000.
  NumericVector test_preds;							// create a vector.
  List cp_mat_list=resids_cp_mat_tau;
  
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  arma::mat temp_all_resids_a(resids.begin(), resids.nrow(), resids.ncol(), false);	
  arma::mat temp_treated_resids_a = temp_all_resids_a.rows(arma::find(z_ar==1));	
  NumericMatrix temp_treated_resids = wrap(temp_treated_resids_a);	
  
  
  if(zero_split==1){
    
  for(int q=0; q<parent.size();q++){
    
  
  SEXP s = prev_sum_trees_mu[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
  if(is<List>(s)){												// is s is a list

    List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
    List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
    
    List st_tau(1);													// create list, st, of length 2.
    List st_mat_tau(1);												// create lisr, st_mat of length 2.
    st_tau[0]=tree_table_tau[0];												// let the first elemetn of st be sum_trees2.
    st_mat_tau[0]=tree_mat_tau[0];										// let the first element of st_mat be sum_trees_mat2.
    
    double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,st_tau,prev_sum_trees_mat_mu2_temp,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
    double tree_prior_temp=1;
    int p_other_mu=0;
    NumericVector other_int_nodes_mu;
    
    
      for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
      NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
      other_int_nodes_mu = find_term_nodes_bcf(tree);
      p_other_mu+=other_int_nodes_mu.size();
      NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
      //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
      if(tree.ncol()<5) throw std::range_error("Line 1412");
      tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
    }
      tree_prior_temp*=get_tree_prior_bcf(tree_table_tau[0],tree_mat_tau[0],alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
 
 NumericVector int_nodes_tau=find_term_nodes_bcf(tree_table_tau[0]);				// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector? 
 int p_other_tau=int_nodes_tau.size();											// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
 double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);	
 
 overall_trees[overall_count]=tree_table_tau[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
 overall_mat[overall_count]=tree_mat_tau[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
 overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)

  overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
 overall_count++;												// increment overall_count by one.
 
 // Rcout << "get to end ifelse function tau round one. \n";
 
  }else{															// if s is not a list

    NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
    NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
    int p_other_mu=0;
    NumericVector other_int_nodes_mu;

    other_int_nodes_mu = find_term_nodes_bcf(prev_sum_trees_mu2_temp);
    p_other_mu=other_int_nodes_mu.size();
    List st_mu(1);													// create list, st, of length 2.
    List st_mat_mu(1);												// create lisr, st_mat of length 2.
    st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
    st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
    // return(st);
    List st_tau(1);													// create list, st, of length 2.
    List st_mat_tau(1);												// create lisr, st_mat of length 2.
    st_tau[0]=tree_table_tau[0];												// let the first elemetn of st be sum_trees2.
    st_mat_tau[0]=tree_mat_tau[0];										// let the first element of st_mat be sum_trees_mat2.

    double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
    double tree_prior_temp=1;

    tree_prior_temp*=get_tree_prior_bcf(prev_sum_trees_mu2_temp,prev_sum_trees_mat_mu2_temp,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)

    tree_prior_temp*=get_tree_prior_bcf(tree_table_tau[0],tree_mat_tau[0],alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)

  NumericVector int_nodes_tau=find_term_nodes_bcf(tree_table_tau[0]);				// find term nodes function defined line 168. Gives index of values of proposal_tree[0] that are term nodes (indices from 1 to length of vector). Why not integer vector? 
  int p_other_tau=int_nodes_tau.size();											// p is length of int_nodes. Number of terminal nodes is used as numbr of parameters/ (B in equation 7 of the paper)
  double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);	

  overall_trees[overall_count]=tree_table_tau[0];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
  overall_mat[overall_count]=tree_mat_tau[0];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
  overall_parent[overall_count]=parent[q];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)

  overall_lik[overall_count]=BIC;			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
  overall_count++;												// increment overall_count by one.
  

  }
  
  if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
    overall_size=overall_size*2;								// double the size
    overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
    overall_lik.resize(overall_size);							// double the size of overall_lik
    overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
    overall_parent.resize(overall_size);						// double the size of overall_parent
  }
  
  
}

  overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
  eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
  overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
  overall_trees=eval_model[1];										// list of tree tables in Occam's window
  overall_mat=eval_model[2];											// list of tree matrices in Occam's window
  overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
  overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
  //add in check to see if OW accepted more than the top maxOW models...
  if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
    //find the maxOWsize best models and continue with those!
    IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
    owindices=owindices-1;											// take one away from indices so that they are in the correct range.
    //get the top maxOWsize indices to keep in OW
    NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
    List temp_otrees(maxOWsize);									// create list of length maxOWsize
    List temp_omat(maxOWsize);										// create list of length maxOWsize
    IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
    
    //now only select those elements
    for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
      temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
      temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
      temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
      temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
    }
    
    overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
    overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
    overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
    overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
    overall_parent2=temp_oparent;	
  
  }
  if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
    overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
    overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
    overall_lik.resize(overall_size);								// increse size of overall_lik
    overall_parent.resize(overall_size);							// increse size of overall_parent
  }else{																// if length of overall_trees is greater than or equal to overall_size-1
    overall_size=2*overall_size;									// double the size
    overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
    overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
    overall_lik.resize(overall_size);								// double the length of overall_lik
    overall_parent.resize(overall_size);							// double the length of overall_parent
  }
  
  
  } // end if zero_split=1 code
  
  
  
  for(int j=0;j<num_splits_tau;j++){									// create a for-loop of length 5.
    int lsize=1000;										// create a variable initialized equal to 1000.
    List table_subset_curr_round(lsize);				// create a list of length 1000
    std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000
    List mat_subset_curr_round(lsize);					// create a list of length 1000.
    std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
    int count=0;										// create a varialble. Initialized equal to 0.
    for(int i=0;i<tree_table_tau.size();i++){				// for loop of length equal to that of the list tree_table (input or updated in previous iteration)
      //if(first_round==1){								// if input first_round equals 1.
        //parent=-1;									// reset input parent to -1
        NumericMatrix temp_list=cp_mat_list[0];		// create matrix temp_list equal to first element of cp_mat_list
        // while the tau tree matrices are new, they are being appended to a sum of tree model, and the correct corresponding mu matrix must be in each model
        
        //NEED TO FIND OUT WHAT TO DO WITH parent index and cp_mat_list, is it different for j==0?
        // Even for j==0, cp_mat_list has more than one element because there are different residuals from the different single mu tree models
        // therefore it seems to make sense to use parent[i] where possible... hopefully this is correct
        best_subset=get_best_split_tau_round1_bcf(resids(_,parent[i]),x_moderate_a,tree_table_tau[i],tree_mat_tau[i],
                                              a_mu,a_tau,mu_mu,mu_tau,nu,lambda,log(c),
                                              lowest_BIC,parent[i],cp_mat_list[parent[i]],
                                              alpha_mu,beta_mu,alpha_tau,beta_tau,
                                              maxOWsize,//first_round,
                                              prev_sum_trees_mu,prev_sum_trees_mat_mu,
                                              y_scaled,parent,i,z);
        
        //something between the following two functions should be used
        //best_subset=get_best_split_bcf(resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,first_round);	// defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
        //best_subset=get_best_split_sum(resids(_,parent[i]),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[i],cp_mat_list[parent[i]],alpha,beta,maxOWsize,first_round,prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i);
      // }else{		// if input first_round is not equal to 1.
      //   throw std::range_error("get_best_trees_sum_tau_round1_bcf should only be used in the first round");	// throw an error.
      //   
      //   //if(err_list[i]==0){												// if i+1^th element of input vector err_list equals zero (no error?)
      //   //	NumericMatrix test_tree=tree_table[i];						// create test_tree equal to i+1^th element of input list tree_table
      //   //	NumericMatrix test_treemat=tree_mat[i];						// create test_treemat equal to i+1^th element of input list tree_mat
      //   //	NumericMatrix test_cpmat= cp_mat_list[parent[i]];			// create a matrix test_cpmat equal to the parent[i]+1^th element of cp_mat_list
      //   //need to append current tree_table[i] to its parent sum_of_trees   
      //   //	best_subset=get_best_split_sum(resids(_,parent[i]),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[i],cp_mat_list[parent[i]],alpha,beta,maxOWsize,first_round,prev_sum_trees,prev_sum_trees_mat,y_scaled,parent,i);	// defined on line 1074. Returns the BIC, best splitting variable (column number), value of covariate for splitting variable, list including tree table and tree matrix, list of tree tables, list of BICs, list of tree matrices, tree parent vector.
      //   // return(best_subset);
      //   //	}else if(err_list[i]==1){		// if i+1^th element of input vector err_list equals one (no error?).
      //   //	continue;						// skip to next iteration of for-loop.
      //   //	}else{							// if i+1^th element of input vector err_list is not 0 or 1 ???.
      //   //	List ret_list(6);				// list of length 6.
      //   //	ret_list[0]=9999;				// first element of list equals 9999.
      //   //	ret_list[1]=err_list[i];		// second element is i+1^th element of input vector err_list. (some number not equal to 0 or 1).
      //   //	ret_list[2]=i;					// third element of list is number indexing the element of tree_table (index starts from 0).
      //   //	ret_list[3]=j;					// fourth element is between 0 and 4. Index of outer for-loop.
      //   //	ret_list[4]=tree_table;			// fifth element of list is the list of tree tables.
      //   //	ret_list[5]=err_list;			// sixth element of list is vector err_list.
      //   //	return(ret_list);				// return the list.
      //   //	throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");	// throw an error.
      //   //}
      // }    
      
      
      if(best_subset.size()==1){	// if get_best_split_bcf resulted in a (?no ?sum-of? tree?) error.
        //Rcout << "There was a no tree error. \n";
        continue;				// skip to the next iteration of the inner for-loop.
      }
      List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_sum
      List temp_mat=best_subset[6];					// list of tree matrices
      lik_list=best_subset[5];						// vector of BICs
      IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?...not -1 if past first round)
      if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
        throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
      }
      if(lik_list.size()==0){							// if vector of BICs has no elements
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){						// If the minimum of the list of BICs is less than the input value of lowest_BIC.
        lowest_BIC=min(lik_list);						// reset lowest_BIC equal to the minimum of the list of BICs
      }
      for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
        table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
        lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
        mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
        parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
        count++;										// increment the count by 1.
        
        if(count==(lsize-1)){							// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
          lsize=lsize*2;								// Double the length
          table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
          mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
          lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
          parent_curr_round.resize(lsize);										// double length of parent_curr_round
        }
      }
    }
    table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
    mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
    lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
    parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0 
      break;															// break out of the for-loop,
    }
    for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
      overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
      overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
      overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
      overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
      overall_count++;												// increment overall_count by one.
      
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    
    tree_table_tau=table_subset_curr_round;									// reset tree_table_tau equal to the list of tables produced in the current iteration of the outer loop.
    IntegerVector temp1(table_subset_curr_round.size(),1);				// create a vector of ones of length equal to that of the list of tree tables for the current round.
    err_list=temp1;														// set err_list equal to temp1
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    tree_mat_tau=mat_subset_curr_round;										// reset tree_mat to list of treematrices obtained in current round of for-loop)
    parent=parent_curr_round;											// reset parent to parent vector obtained in current round of for-loop.
    
    if(split_rule_node==1){												// If input boolean split_rule_node = 1 (true)
      NumericVector temp_preds;										// create vector. Length not given.
      List updated_curr_preds;										// create list. Length not given.
      NumericVector new_mean;											// create vector. Length not given.
      lowest_BIC=min(overall_lik2);									// update lowest_BIC
      NumericMatrix curr_resids(resids.nrow(),resids.ncol());			// curr_resids is a matrix with dimensions of input matrix resids
      
      for(int k=0;k<table_subset_curr_round.size();k++){				// for-loop of length equal to number of models proposed in the current iteration
        NumericVector terminal_nodes;								// create a vector called terminal_nodes.
        
        NumericMatrix temp_all_mat = mat_subset_curr_round[k];
        
        arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
        arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
        NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);
        
        if(parent_curr_round[k]==-1){								// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
        }else{														// if the k+1^th element of parent_curr_round does not equal -1 
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,parent_curr_round[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
        }  
        
        terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);			// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
        updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_moderate_a.n_rows);	// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
        NumericVector test_res;												// create vector called test_res. Length not given.
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          test_res=resids(_,0);											// let test_res equal the first column of resids
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 
          test_res=resids(_,parent_curr_round[k]);						// let test_res equal the parent_curr_round[k]+1^th column of resids
        }
        
        
        NumericVector temp_preds = updated_curr_preds[1];
        NumericVector curr_test_res=z*temp_preds;				//might be unnecessary to multiply by z because will probably subset correctly when this vector is used	// Let curr_test_res be the vector of updated predictions
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          curr_resids(_,0)=test_res-curr_test_res;						// let first column of curr_resids be ?the updated residuals?
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
          curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
        }
      }
      List temp(0);																	// list of length zero.
      
      cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
      
      for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
        List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
        if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
          cp_mat_list1=make_pelt_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),pen_tau,num_cp_tau,z);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																			// If gridpoint equals 1 (or non-zero).
          cp_mat_list1=make_gridpoint_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),gridsize_tau,num_cp_tau,z);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        
        cp_mat_list.push_back(cp_mat_list1[0]);											// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
      }
      
    }
  }
  overall_trees=resize_bcf(overall_trees,overall_count);		// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);			// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);					// remove remaining spaces that are not filled in
  NumericVector temp_preds;								// create a vector
  List updated_preds;										// create a list
  NumericVector new_mean;									// create a vector
  NumericMatrix overall_test_preds(x_moderate_test.nrow(),overall_trees.size());	// create a matrix with no. of rows equal to that of input x_moderate_test, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_moderate_a.n_rows,overall_trees.size());		// create a matrix with no. of rows equal to that of x_moderate_a (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);							// reset lowest_BIC to minimum of overall_lik2
  
  
  for(int k=0;k<overall_trees.size();k++){			// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;					// create a vector
    
    NumericMatrix temp_all_mat = overall_mat[k];
    
    arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
    arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
    NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);	
    
    if(overall_parent2[k]==-1){						// If k+1^th element of overall_parent2 equals -1
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
    }else{    
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,overall_parent2[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
    }  
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);		// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_moderate_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_moderate_test,overall_trees[k],new_mean);		// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    temp_preds=updated_preds[1];																// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;											// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);			// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);														// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);		// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_test_values=sum(M2,1);		// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_tau_bcf(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                            double a_mu,double a_tau,double mu_mu,double mu_tau,
                            double nu,double lambda,double c,double sigma_mu_mu,double sigma_mu_tau,
                            List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                            double lowest_BIC,//int first_round,
                            IntegerVector parent,
                            List resids_cp_mat_tau,IntegerVector err_list,
                            NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                            double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,
                            bool is_test_data,
                            double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                            bool split_rule_node,bool gridpoint,int maxOWsize,
                            List prev_sum_trees_mu,List prev_sum_trees_tau,List prev_sum_trees_mat_mu,List prev_sum_trees_mat_tau,
                            NumericVector y_scaled,
                            int num_splits_mu,int num_splits_tau, int gridsize_tau,bool zero_split
){
  List eval_model;										// create a list
  NumericVector lik_list;								// create a vector
  List best_subset;									// create a list
  int overall_size=1000;								// create variable overall_size. Initialize equal to 1000.
  List overall_trees(overall_size);					// create a list of length 1000.
  NumericVector overall_lik2;							// create a vector.
  IntegerVector overall_parent2;						// create a vector.
  List overall_mat(overall_size);						// create a list of length 1000.
  int overall_count=0;								// create a variable initialized equal to 0.
  std::vector<int> overall_parent(overall_size);		// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);		// create a vector of length 1000.
  NumericVector test_preds;							// create a vector.
  
  List cp_mat_list = resids_cp_mat_tau;
  
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  arma::mat temp_all_resids_a(resids.begin(), resids.nrow(), resids.ncol(), false);	
  arma::mat temp_treated_resids_a = temp_all_resids_a.rows(arma::find(z_ar==1));	
  NumericMatrix temp_treated_resids = wrap(temp_treated_resids_a);	
  
  
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  
  
  ////BEGGINING OF ZERO SPLIT TREE CODE
  ///////////////////////////////////////////////////////////////////////
  
  if(zero_split==1){
    
  for(int q=0; q<parent.size();q++){
    
    
    SEXP s_mu = prev_sum_trees_mu[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
    SEXP s_tau = prev_sum_trees_tau[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
    if(is<List>(s_mu)){												// is s is a list
      if(is<List>(s_tau)){												// is s is a list
        //Rcout << "\n RELEVANT LINE 1703.\n";
        
        List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]


        sum_trees_tau2_temp.push_back(tree_table_tau[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
        sum_trees_mat_tau2_temp.push_back(tree_mat_tau[0]);	

        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,sum_trees_tau2_temp,prev_sum_trees_mat_mu2_temp,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood

        double tree_prior_temp=1;
        
        int p_other_mu=0;
        int p_other_tau=0;
        
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        
        
        
        for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1660");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          // Rcout << "LINE 4378 t =" << t << ".\n";
          // Rcout << "LINE 4378 tree_prior_temp =" << tree_prior_temp << ".\n";
          
        }
        for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1667");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          // Rcout << "LINE 4390 t =" << t << ".\n";
          // Rcout << "LINE 4390 tree_prior_temp =" << tree_prior_temp << ".\n";
        }
        
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
        // Rcout << "LINE 4395 BIC =" << BIC << ".\n";
        // Rcout << "LINE 4395 lik_temp =" << lik_temp << ".\n";
        // Rcout << "LINE 4395 tree_prior_temp =" << tree_prior_temp << ".\n";
        
        overall_trees[overall_count]=tree_table_tau[0];		
        overall_mat[overall_count]=tree_mat_tau[0];			
        overall_parent[overall_count]=parent[q];				
        
        overall_lik[overall_count]=BIC;			
        
        overall_count++;												// increment overall_count by one.
        
      }else{
        List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_tau(2);													// create list, st, of length 2.
        List st_mat_tau(2);												// create lisr, st_mat of length 2.
        st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
        st_tau[1]=tree_table_tau[0];										// let the second element of st be proposal_tree[0] (tree table).
        st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.
        st_mat_tau[1]=tree_mat_tau[0];		
        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,st_tau,prev_sum_trees_mat_mu2_temp,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
        
        double tree_prior_temp=1;
        
        int p_other_mu=0;
        int p_other_tau=0;
        
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        
        
        for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1689");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1695");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
        // Rcout << "LINE 4440 BIC =" << BIC << ".\n";
        
        overall_trees[overall_count]=tree_table_tau[0];		
        overall_mat[overall_count]=tree_mat_tau[0];			
        overall_parent[overall_count]=parent[q];				
        
        overall_lik[overall_count]=BIC;			
        
        overall_count++;												// increment overall_count by one.
        
      }
      
    }else{															// if s is not a list
      if(is<List>(s_tau)){												// is s is a list
        
        NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
        NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_mu(1);													// create list, st, of length 2.
        List st_mat_mu(1);												// create lisr, st_mat of length 2.
        st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
        st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
        // return(st);
        List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        
        sum_trees_tau2_temp.push_back(tree_table_tau[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
        sum_trees_mat_tau2_temp.push_back(tree_mat_tau[0]);	
        //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
        //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,sum_trees_tau2_temp,st_mat_mu,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
        
        double tree_prior_temp=1;
        
        int p_other_mu=0;
        int p_other_tau=0;
        
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        
        
        for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
          NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1723");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1730");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
        // Rcout << "LINE 4500 BIC =" << BIC << ".\n";
        
        overall_trees[overall_count]=tree_table_tau[0];		
        overall_mat[overall_count]=tree_mat_tau[0];			
        overall_parent[overall_count]=parent[q];				
        
        overall_lik[overall_count]=BIC;			
        
        overall_count++;												// increment overall_count by one.
        
      }else{
        NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
        NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_mu(1);													// create list, st, of length 2.
        List st_mat_mu(1);												// create lisr, st_mat of length 2.
        st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
        st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
        NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
        NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
        List st_tau(2);													// create list, st, of length 2.
        List st_mat_tau(2);												// create lisr, st_mat of length 2.
        st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
        st_tau[1]=tree_table_tau[0];										// let the second element of st be proposal_tree[0] (tree table).
        st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.
        st_mat_tau[1]=tree_mat_tau[0];		
        double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
        double tree_prior_temp=1;
        
        int p_other_mu=0;
        int p_other_tau=0;
        
        NumericVector other_int_nodes_mu;
        NumericVector other_int_nodes_tau;
        
        for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
          NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
          other_int_nodes_mu = find_term_nodes_bcf(tree);
          p_other_mu+=other_int_nodes_mu.size();
          NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1753");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
          NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
          other_int_nodes_tau = find_term_nodes_bcf(tree);
          p_other_tau+=other_int_nodes_tau.size();
          NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
          //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
          if(tree.ncol()<5) throw std::range_error("Line 1760");
          tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
        }
        double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
        // Rcout << "LINE 4553 BIC =" << BIC << ".\n";
        
        overall_trees[overall_count]=tree_table_tau[0];		
        overall_mat[overall_count]=tree_mat_tau[0];			
        overall_parent[overall_count]=parent[q];				
        
        overall_lik[overall_count]=BIC;			
        
        overall_count++;												// increment overall_count by one.
        
      }
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
  }
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    
  
  
////END OF ZERO SPLIT TREE CODE
///////////////////////////////////////////////////////////////////////
  } // end of if zero_split=1 code
  
  
  
  
  for(int j=0;j<num_splits_tau;j++){									// create a for-loop of length 5.
    int lsize=1000;										// create a variable initialized equal to 1000.
    List table_subset_curr_round(lsize);				// create a list of length 1000
    std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000
    List mat_subset_curr_round(lsize);					// create a list of length 1000.
    std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
    int count=0;										// create a varialble. Initialized equal to 0.
    for(int i=0;i<tree_table_tau.size();i++){				// for loop of length equal to that of the list tree_table (input or updated in previous iteration)
      // if(first_round==1){								// if input first_round equals 1.
      //   throw std::range_error("get_best_trees_sum_tau_bcf should not be used in the first round");	// throw an error.
      //   //parent=-1;									// reset input parent to -1
      //   //NumericMatrix temp_list=cp_mat_list[0];		// create matrix temp_list equal to first element of cp_mat_list
      //   //best_subset=get_best_split_bcf(resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,first_round);	// defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
      // }else{											// if input first_round is not equal to 1.
        if(err_list[i]==0){												// if i+1^th element of input vector err_list equals zero (no error?)
          //NumericMatrix test_tree=tree_table[i];						// create test_tree equal to i+1^th element of input list tree_table
          //NumericMatrix test_treemat=tree_mat[i];						// create test_treemat equal to i+1^th element of input list tree_mat
          //NumericMatrix test_cpmat= cp_mat_list[parent[i]];			// create a matrix test_cpmat equal to the parent[i]+1^th element of cp_mat_list
          //need to append current tree_table[i] to its parent sum_of_trees   
          
          List prev_sum_trees_mu2=prev_sum_trees_mu[parent[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List sum_trees_tau2=prev_sum_trees_tau[parent[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          //Rcout << "Length of mat list = " << prev_sum_trees_mat_mu2.size() <<".\n";
          //NumericMatrix testexamplemat5 = prev_sum_trees_mat_mu2[0];
          //Rcout << "number of rows of mat into likelihood = " << testexamplemat5.nrow() <<".\n";
          
          best_subset=get_best_split_sum_tau_bcf(resids(_,parent[i]),x_moderate_a,tree_table_tau[i],tree_mat_tau[i],
                                             a_mu,a_tau,mu_mu,mu_tau,nu,lambda,log(c),
                                             lowest_BIC,parent[i],cp_mat_list[parent[i]],
                                             alpha_mu,beta_mu,alpha_tau,beta_tau,
                                             maxOWsize,//first_round,
                                             prev_sum_trees_mu,prev_sum_trees_tau,prev_sum_trees_mat_mu,prev_sum_trees_mat_tau,
                                             y_scaled,parent,i,z);	// defined on line 1074. Returns the BIC, best splitting variable (column number), value of covariate for splitting variable, list including tree table and tree matrix, list of tree tables, list of BICs, list of tree matrices, tree parent vector.
          // return(best_subset);
        }else if(err_list[i]==1){		// if i+1^th element of input vector err_list equals one (no error?).
          continue;						// skip to next iteration of for-loop.
        }else{							// if i+1^th element of input vector err_list is not 0 or 1 ???.
          List ret_list(6);				// list of length 6.
          ret_list[0]=9999;				// first element of list equals 9999.
          ret_list[1]=err_list[i];		// second element is i+1^th element of input vector err_list. (some number not equal to 0 or 1).
          ret_list[2]=i;					// third element of list is number indexing the element of tree_table (index starts from 0).
          ret_list[3]=j;					// fourth element is between 0 and 4. Index of outer for-loop.
          ret_list[4]=tree_table_tau;			// fifth element of list is the list of tree tables.
          ret_list[5]=err_list;			// sixth element of list is vector err_list.
          return(ret_list);				// return the list.
          throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");	// throw an error.
        }
      //}    
      if(best_subset.size()==1){	// if get_best_split_bcf resulted in a (?no ?sum-of? tree?) error.
        continue;				// skip to the next iteration of the inner for-loop.
      }
      List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_sum
      List temp_mat=best_subset[6];					// list of tree matrices
      lik_list=best_subset[5];						// vector of BICs
      IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?)
      if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
        throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
      }
      if(lik_list.size()==0){							// if vector of BICs has no elements
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){						// If the minimum of the list of BICs is less than the input value of lowest_BIC.
        lowest_BIC=min(lik_list);						// reset lowest_BIC equal to the minimum of the list of BICs
      }
      for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
        table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
        lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
        mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
        parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
        count++;										// increment the count by 1.
        
        if(count==(lsize-1)){							// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
          lsize=lsize*2;								// Double the length
          table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
          mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
          lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
          parent_curr_round.resize(lsize);										// double length of parent_curr_round
        }
      }
    }
    table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
    mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
    lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
    parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0 
      break;															// break out of the for-loop,
    }
    for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
      overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
      overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
      overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
      overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
      overall_count++;												// increment overall_count by one.
      
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    
    tree_table_tau=table_subset_curr_round;									// reset tree_table equal to the list of tables produced in the current iteration of the outer loop.
    IntegerVector temp1(table_subset_curr_round.size(),1);				// create a vector of ones of length equal to that of the list of tree tables for the current round.
    err_list=temp1;														// set err_list equal to temp1
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    tree_mat_tau=mat_subset_curr_round;										// reset tree_mat to list of treematrices obtained in current round of for-loop)
    parent=parent_curr_round;											// reset parent to parent vector obtained in current round of for-loop.
    
    if(split_rule_node==1){												// If input boolean split_rule_node = 1 (true)
      NumericVector temp_preds;										// create vector. Length not given.
      List updated_curr_preds;										// create list. Length not given.
      NumericVector new_mean;											// create vector. Length not given.
      lowest_BIC=min(overall_lik2);									// update lowest_BIC
      NumericMatrix curr_resids(resids.nrow(),resids.ncol());			// curr_resids is a matrix with dimensions of input matrix resids
      
      for(int k=0;k<table_subset_curr_round.size();k++){				// for-loop of length equal to number of models proposed in the current iteration
        NumericVector terminal_nodes;								// create a vector called terminal_nodes.
        
        NumericMatrix temp_all_mat = mat_subset_curr_round[k];
        
        arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
        arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
        NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);
        
        if(parent_curr_round[k]==-1){								// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
        }else{														// if the k+1^th element of parent_curr_round does not equal -1 
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,parent_curr_round[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
        }  
        
        terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);			// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
        updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_moderate_a.n_rows);	// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
        NumericVector test_res;												// create vector called test_res. Length not given.
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          test_res=resids(_,0);											// let test_res equal the first column of resids
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 
          test_res=resids(_,parent_curr_round[k]);						// let test_res equal the parent_curr_round[k]+1^th column of resids
        }
        
        NumericVector temppreds = updated_curr_preds[1];
        NumericVector curr_test_res=z*temppreds;				//might be unnecessary to multiply by z because will probably subset correctly when this vector is used	// Let curr_test_res be the vector of updated predictions
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          curr_resids(_,0)=test_res-curr_test_res;						// let first column of curr_resids be ?the updated residuals?
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
          curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
        }
      }
      List temp(0);																	// list of length zero.
      
      cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
      
      for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
        List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
        if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
          cp_mat_list1=make_pelt_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),pen_tau,num_cp_tau,z);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																			// If gridpoint equals 1 (or non-zero).
          cp_mat_list1=make_gridpoint_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),gridsize_tau,num_cp_tau,z);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        
        cp_mat_list.push_back(cp_mat_list1[0]);											// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
      }
      
    }
  }
  overall_trees=resize_bcf(overall_trees,overall_count);		// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);			// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);					// remove remaining spaces that are not filled in
  NumericVector temp_preds;								// create a vector
  List updated_preds;										// create a list
  NumericVector new_mean;									// create a vector
  NumericMatrix overall_test_preds(x_moderate_test.nrow(),overall_trees.size());	// create a matrix with no. of rows equal to that of input x_moderate_test, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_moderate_a.n_rows,overall_trees.size());		// create a matrix with no. of rows equal to that of x_moderate_a (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);							// reset lowest_BIC to minimum of overall_lik2
  
  
  for(int k=0;k<overall_trees.size();k++){			// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;					// create a vector
    
    
    NumericMatrix temp_all_mat = overall_mat[k];
    
    arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
    arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
    NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);	
    
    if(overall_parent2[k]==-1){						// If k+1^th element of overall_parent2 equals -1
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
    }else{    
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,overall_parent2[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
    }  
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);		// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_moderate_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_moderate_test,overall_trees[k],new_mean);		// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    temp_preds=updated_preds[1];																// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;											// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);			// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);														// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);		// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_test_values=sum(M2,1);		// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}
//######################################################################################################################//
// [[Rcpp::export]]

List get_best_trees_sum_tau_bcf_2(arma::mat& x_control_a,arma::mat& x_moderate_a,NumericVector z,NumericMatrix resids,
                                double a_mu,double a_tau,double mu_mu,double mu_tau,
                                double nu,double lambda,double c,double sigma_mu_mu,double sigma_mu_tau,
                                List tree_table_mu,List tree_mat_mu,List tree_table_tau,List tree_mat_tau,
                                double lowest_BIC,//int first_round,
                                IntegerVector parent,
                                List resids_cp_mat_tau,IntegerVector err_list,
                                NumericMatrix x_control_test,NumericMatrix x_moderate_test,NumericVector test_z,
                                double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,
                                bool is_test_data,
                                double pen_mu,int num_cp_mu,double pen_tau,int num_cp_tau,
                                bool split_rule_node,bool gridpoint,int maxOWsize,
                                List prev_sum_trees_mu,List prev_sum_trees_tau,List prev_sum_trees_mat_mu,List prev_sum_trees_mat_tau,
                                NumericVector y_scaled,
                                int num_splits_mu,int num_splits_tau, int gridsize_tau,bool zero_split,IntegerVector no_more_tau_trees
){
  List eval_model;										// create a list
  NumericVector lik_list;								// create a vector
  List best_subset;									// create a list
  int overall_size=1000;								// create variable overall_size. Initialize equal to 1000.
  List overall_trees(overall_size);					// create a list of length 1000.
  NumericVector overall_lik2;							// create a vector.
  IntegerVector overall_parent2;						// create a vector.
  List overall_mat(overall_size);						// create a list of length 1000.
  int overall_count=0;								// create a variable initialized equal to 0.
  std::vector<int> overall_parent(overall_size);		// create a vector of length 1000.
  std::vector<double> overall_lik(overall_size);		// create a vector of length 1000.
  NumericVector test_preds;							// create a vector.
  
  List cp_mat_list = resids_cp_mat_tau;
  
  arma::vec z_ar=Rcpp::as<arma::vec>(z);		// converts to arma vec
  arma::mat temp_all_resids_a(resids.begin(), resids.nrow(), resids.ncol(), false);	
  arma::mat temp_treated_resids_a = temp_all_resids_a.rows(arma::find(z_ar==1));	
  NumericMatrix temp_treated_resids = wrap(temp_treated_resids_a);	
  
  // Rcout << "Line 6501 no_more_tau_trees = " << no_more_tau_trees << ".\n";
  // Rcout << "parent = " << parent << ".\n";
  
  //////   //COMMENTED OUT ATTEMPT TO FIX NO ZERO SPLIT TREES BUG
  
  
  ////BEGGINING OF ZERO SPLIT TREE CODE
  ///////////////////////////////////////////////////////////////////////
  
  if(zero_split==1){
    if( is_true(all(no_more_tau_trees==1)) ){
    }
    else{
    for(int q=0; q<parent.size();q++){
      
      // Rcout << "Line 6513. q = " << q <<".\n";
      
      // Rcout << "no_more_tau_trees[q] = " << no_more_tau_trees[q] <<".\n";
      
      
      if(no_more_tau_trees[q]!=1){
        // Rcout << "Line 6519. q = " << q <<".\n";
        
      SEXP s_mu = prev_sum_trees_mu[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
      SEXP s_tau = prev_sum_trees_tau[parent[q]];									// s is a pointer to S expression type equal to the element of the sumtrees input list indexed by the i^th element of the input Integer vetor parent2. i is also an input integer.
      if(is<List>(s_mu)){												// is s is a list
        if(is<List>(s_tau)){												// is s is a list
          //Rcout << "\n RELEVANT LINE 1703.\n";
          
          List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          
          
          sum_trees_tau2_temp.push_back(tree_table_tau[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
          sum_trees_mat_tau2_temp.push_back(tree_mat_tau[0]);	
          
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,sum_trees_tau2_temp,prev_sum_trees_mat_mu2_temp,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          
          double tree_prior_temp=1;
          
          int p_other_mu=0;
          int p_other_tau=0;
          
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          
          
          
          for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1660");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            // Rcout << "LINE 4378 t =" << t << ".\n";
            // Rcout << "LINE 4378 tree_prior_temp =" << tree_prior_temp << ".\n";
            
          }
          for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1667");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
            // Rcout << "LINE 4390 t =" << t << ".\n";
            // Rcout << "LINE 4390 tree_prior_temp =" << tree_prior_temp << ".\n";
          }
          
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          // Rcout << "LINE 4395 BIC =" << BIC << ".\n";
          // Rcout << "LINE 4395 lik_temp =" << lik_temp << ".\n";
          // Rcout << "LINE 4395 tree_prior_temp =" << tree_prior_temp << ".\n";
          
          overall_trees[overall_count]=tree_table_tau[0];		
          overall_mat[overall_count]=tree_mat_tau[0];			
          overall_parent[overall_count]=parent[q];				
          
          overall_lik[overall_count]=BIC;			
          
          overall_count++;												// increment overall_count by one.
          
        }else{
          List prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_tau(2);													// create list, st, of length 2.
          List st_mat_tau(2);												// create lisr, st_mat of length 2.
          st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
          st_tau[1]=tree_table_tau[0];										// let the second element of st be proposal_tree[0] (tree table).
          st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.
          st_mat_tau[1]=tree_mat_tau[0];		
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,prev_sum_trees_mu2_temp,st_tau,prev_sum_trees_mat_mu2_temp,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          
          double tree_prior_temp=1;
          
          int p_other_mu=0;
          int p_other_tau=0;
          
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          
          
          for(int t=0;t<prev_sum_trees_mu2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=prev_sum_trees_mu2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=prev_sum_trees_mat_mu2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1689");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1695");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          // Rcout << "LINE 4440 BIC =" << BIC << ".\n";
          
          overall_trees[overall_count]=tree_table_tau[0];		
          overall_mat[overall_count]=tree_mat_tau[0];			
          overall_parent[overall_count]=parent[q];				
          
          overall_lik[overall_count]=BIC;			
          
          overall_count++;												// increment overall_count by one.
          
        }
        
      }else{															// if s is not a list
        if(is<List>(s_tau)){												// is s is a list
          
          NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
          NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_mu(1);													// create list, st, of length 2.
          List st_mat_mu(1);												// create lisr, st_mat of length 2.
          st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
          st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
          // return(st);
          List sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          
          sum_trees_tau2_temp.push_back(tree_table_tau[0]);						// append the treetable proposal_tree[0] to the end of the list sum_trees2
          sum_trees_mat_tau2_temp.push_back(tree_mat_tau[0]);	
          //NumericMatrix prop_table_tau_temp = proposal_tree[0];						// append the treetable proposal_tree[0] to the end of the list sum_trees2
          //NumericMatrix prop_mat_tau_temp =proposal_tree[1];					// append the tree matreix proposal_tree[1] to the end of the list sum_trees_mat2
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,sum_trees_tau2_temp,st_mat_mu,sum_trees_mat_tau2_temp,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          
          double tree_prior_temp=1;
          
          int p_other_mu=0;
          int p_other_tau=0;
          
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          
          
          for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
            NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1723");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          for(int t=0;t<sum_trees_tau2_temp.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=sum_trees_tau2_temp[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=sum_trees_mat_tau2_temp[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1730");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          // Rcout << "LINE 4500 BIC =" << BIC << ".\n";
          
          overall_trees[overall_count]=tree_table_tau[0];		
          overall_mat[overall_count]=tree_mat_tau[0];			
          overall_parent[overall_count]=parent[q];				
          
          overall_lik[overall_count]=BIC;			
          
          overall_count++;												// increment overall_count by one.
          
        }else{
          NumericMatrix prev_sum_trees_mu2_temp=prev_sum_trees_mu[parent[q]];				// sum_trees2 is the element of the input list sum_trees indexed by parent2[i] 
          NumericMatrix prev_sum_trees_mat_mu2_temp=prev_sum_trees_mat_mu[parent[q]];		// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_mu(1);													// create list, st, of length 2.
          List st_mat_mu(1);												// create lisr, st_mat of length 2.
          st_mu[0]=prev_sum_trees_mu2_temp;												// let the first elemetn of st be sum_trees2.
          st_mat_mu[0]=prev_sum_trees_mat_mu2_temp;										// let the first element of st_mat be sum_trees_mat2.
          NumericMatrix sum_trees_tau2_temp=prev_sum_trees_tau[parent[q]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          NumericMatrix sum_trees_mat_tau2_temp=prev_sum_trees_mat_tau[parent[q]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List st_tau(2);													// create list, st, of length 2.
          List st_mat_tau(2);												// create lisr, st_mat of length 2.
          st_tau[0]=sum_trees_tau2_temp;												// let the first elemetn of st be sum_trees2.
          st_tau[1]=tree_table_tau[0];										// let the second element of st be proposal_tree[0] (tree table).
          st_mat_tau[0]=sum_trees_mat_tau2_temp;										// let the first element of st_mat be sum_trees_mat2.
          st_mat_tau[1]=tree_mat_tau[0];		
          double lik_temp=sumtree_likelihood_function_bcf_bcf(y_scaled,st_mu,st_tau,st_mat_mu,st_mat_tau,y_scaled.size(),a_mu,a_tau,nu,lambda,z);  // Defined on line 855. Returns the lof marginal likelihood
          double tree_prior_temp=1;
          
          int p_other_mu=0;
          int p_other_tau=0;
          
          NumericVector other_int_nodes_mu;
          NumericVector other_int_nodes_tau;
          
          for(int t=0;t<st_mu.size();t++){									// for-loop of length equal to that of st (which should be length 2)
            NumericMatrix tree=st_mu[t];									// let tree equal (t+1)^th element of st
            other_int_nodes_mu = find_term_nodes_bcf(tree);
            p_other_mu+=other_int_nodes_mu.size();
            NumericMatrix mat=st_mat_mu[t];								// let mat equal (t+1)^th element of st_mat
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1753");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_mu,beta_mu);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          for(int t=0;t<st_tau.size();t++){							// for-loop of length equal to that of sum_trees2
            NumericMatrix tree=st_tau[t];							// tree equals (t+1)^th element of the list sum_trees2 (tree table)
            other_int_nodes_tau = find_term_nodes_bcf(tree);
            p_other_tau+=other_int_nodes_tau.size();
            NumericMatrix mat=st_mat_tau[t];						// mat equals (t+1)^th element of the list sum_trees_mat2 (tree matrix)
            //THIS SHOULD PROBABLY BE CHANGED TO *= , and actually the prior is still probably not properly defined
            if(tree.ncol()<5) throw std::range_error("Line 1760");
            tree_prior_temp*=get_tree_prior_bcf(tree,mat,alpha_tau,beta_tau);			// iteratively add to tree_prior. get_tree_prior_bcf defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees?)
          }
          double BIC=-2*(lik_temp+log(tree_prior_temp))+(p_other_mu+p_other_tau)*log(x_moderate_a.n_rows);			// x_moderate_a.nrows is number of obs. Not sure why tree_prior is included here. Only need likelihood?
          // Rcout << "LINE 4553 BIC =" << BIC << ".\n";
          
          overall_trees[overall_count]=tree_table_tau[0];		
          overall_mat[overall_count]=tree_mat_tau[0];			
          overall_parent[overall_count]=parent[q];				
          
          overall_lik[overall_count]=BIC;			
          
          overall_count++;												// increment overall_count by one.
          
        }
        if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
          overall_size=overall_size*2;								// double the size
          overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
          overall_lik.resize(overall_size);							// double the size of overall_lik
          overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
          overall_parent.resize(overall_size);						// double the size of overall_parent
        }
      }
      }
    }
    // Rcout << "Line 6760.\n";
    
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    // Rcout << "Line 6767.\n";
    
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    
    // Rcout << "Line 6810.\n";
    
    
    ////END OF ZERO SPLIT TREE CODE
    ///////////////////////////////////////////////////////////////////////
  } // end if all no_more_trees
  } // end of if zero_split=1 code
  
  
  // Rcout << "Line 6818.\n";
  
  
  for(int j=0;j<num_splits_tau;j++){									// create a for-loop of length 5.
    int lsize=1000;										// create a variable initialized equal to 1000.
    List table_subset_curr_round(lsize);				// create a list of length 1000
    std::vector<double> lik_subset_curr_round(lsize);	// create a vector of length 1000
    List mat_subset_curr_round(lsize);					// create a list of length 1000.
    std::vector<int> parent_curr_round(lsize);			// create a vector of length 1000.
    int count=0;										// create a varialble. Initialized equal to 0.
    for(int i=0;i<tree_table_tau.size();i++){				// for loop of length equal to that of the list tree_table (input or updated in previous iteration)
      // if(first_round==1){								// if input first_round equals 1.
      //   throw std::range_error("get_best_trees_sum_tau_bcf should not be used in the first round");	// throw an error.
      //   //parent=-1;									// reset input parent to -1
      //   //NumericMatrix temp_list=cp_mat_list[0];		// create matrix temp_list equal to first element of cp_mat_list
      //   //best_subset=get_best_split_bcf(resids(_,0),D1,tree_table[i],tree_mat[i],a,mu,nu,lambda,log(c),lowest_BIC,parent[0],cp_mat_list[0],alpha,beta,maxOWsize,first_round);	// defined on line 890, Returns list including BICS, tree tables, tree matrices, splitting variables, splitting points, and so on.
      // }else{											// if input first_round is not equal to 1.
        if(err_list[i]==0){												// if i+1^th element of input vector err_list equals zero (no error?)
          //NumericMatrix test_tree=tree_table[i];						// create test_tree equal to i+1^th element of input list tree_table
          //NumericMatrix test_treemat=tree_mat[i];						// create test_treemat equal to i+1^th element of input list tree_mat
          //NumericMatrix test_cpmat= cp_mat_list[parent[i]];			// create a matrix test_cpmat equal to the parent[i]+1^th element of cp_mat_list
          //need to append current tree_table[i] to its parent sum_of_trees   
          
          List prev_sum_trees_mu2=prev_sum_trees_mu[parent[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List prev_sum_trees_mat_mu2=prev_sum_trees_mat_mu[parent[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          List sum_trees_tau2=prev_sum_trees_tau[parent[i]];						// sum_trees2 is the element of the input list sum_trees indexed by parent2[i]
          List sum_trees_mat_tau2=prev_sum_trees_mat_tau[parent[i]];				// sum_trees_mat2 is the element of the input list sum_trees_mat indexed by parent2[i]
          //Rcout << "Length of mat list = " << prev_sum_trees_mat_mu2.size() <<".\n";
          //NumericMatrix testexamplemat5 = prev_sum_trees_mat_mu2[0];
          //Rcout << "number of rows of mat into likelihood = " << testexamplemat5.nrow() <<".\n";
          
          best_subset=get_best_split_sum_tau_bcf(resids(_,parent[i]),x_moderate_a,tree_table_tau[i],tree_mat_tau[i],
                                                 a_mu,a_tau,mu_mu,mu_tau,nu,lambda,log(c),
                                                 lowest_BIC,parent[i],cp_mat_list[parent[i]],
                                                 alpha_mu,beta_mu,alpha_tau,beta_tau,
                                                 maxOWsize,//first_round,
                                                 prev_sum_trees_mu,prev_sum_trees_tau,prev_sum_trees_mat_mu,prev_sum_trees_mat_tau,
                                                 y_scaled,parent,i,z);	// defined on line 1074. Returns the BIC, best splitting variable (column number), value of covariate for splitting variable, list including tree table and tree matrix, list of tree tables, list of BICs, list of tree matrices, tree parent vector.
          // return(best_subset);
        }else if(err_list[i]==1){		// if i+1^th element of input vector err_list equals one (no error?).
          continue;						// skip to next iteration of for-loop.
        }else{							// if i+1^th element of input vector err_list is not 0 or 1 ???.
          List ret_list(6);				// list of length 6.
          ret_list[0]=9999;				// first element of list equals 9999.
          ret_list[1]=err_list[i];		// second element is i+1^th element of input vector err_list. (some number not equal to 0 or 1).
          ret_list[2]=i;					// third element of list is number indexing the element of tree_table (index starts from 0).
          ret_list[3]=j;					// fourth element is between 0 and 4. Index of outer for-loop.
          ret_list[4]=tree_table_tau;			// fifth element of list is the list of tree tables.
          ret_list[5]=err_list;			// sixth element of list is vector err_list.
          return(ret_list);				// return the list.
          throw std::range_error("err_list[i] is neither 0 nor 1...something is wrong here!");	// throw an error.
        }
      //}
      
      
      if(best_subset.size()==1){	// if get_best_split_bcf resulted in a (?no ?sum-of? tree?) error.
        continue;				// skip to the next iteration of the inner for-loop.
      }
      
      List temp_trees=best_subset[4];					// list of tree tables returned by get_best_split_sum
      List temp_mat=best_subset[6];					// list of tree matrices
      lik_list=best_subset[5];						// vector of BICs
      IntegerVector temp_parent=best_subset[7];		// ?vector with all elements equal to the input value of parent? (?-1?)
      if(temp_parent.size()!= temp_trees.size()){		// if length of temp_parent is not equal to the length of the list of tree tables.
        throw std::range_error("there should be a parent for each tree!!!");	// throw an error.
      }
      if(lik_list.size()==0){							// if vector of BICs has no elements
        throw std::range_error("like list size is 0 need to find out why should have broke from list before now!");
      }
      
      if(min(lik_list)<lowest_BIC){						// If the minimum of the list of BICs is less than the input value of lowest_BIC.
        lowest_BIC=min(lik_list);						// reset lowest_BIC equal to the minimum of the list of BICs
      }
      for(int k=0;k<temp_trees.size();k++){				// for-loop of length equal to the length of the temporary list of tree tables
        table_subset_curr_round[count]=temp_trees[k];	// across all iterations, add all elements of the temporary list of tree tables to the list table_subset_curr_round
        lik_subset_curr_round[count]=lik_list[k];		// add all BICs to the vector lik_subset_curr_round
        mat_subset_curr_round[count]=temp_mat[k];		// add all tree matrices to the list mat_subset_curr_round
        parent_curr_round[count]=temp_parent[k];		// add all elements of temp_parent to parent_curr_round
        count++;										// increment the count by 1.
        
        if(count==(lsize-1)){							// If count equals lsize-1, i.e. the length of the lists and vectors needs to be increased. (if inital value of 1000, or previous reset value is not enough)
          lsize=lsize*2;								// Double the length
          table_subset_curr_round=resize_bigger_bcf(table_subset_curr_round,lsize);	// double length of table_subset_curr_round
          mat_subset_curr_round=resize_bigger_bcf(mat_subset_curr_round,lsize);		// double length of mat_subset_curr_round
          lik_subset_curr_round.resize(lsize);									// double length of lik_subset_curr_round
          parent_curr_round.resize(lsize);										// double length of parent_curr_round
        }
      }
    }
    table_subset_curr_round=resize_bcf(table_subset_curr_round,count);		// remove remaining spaces that are not filled in
    mat_subset_curr_round=resize_bcf(mat_subset_curr_round,count);			// remove remaining spaces that are not filled in
    lik_subset_curr_round.resize(count);								// remove remaining spaces that are not filled in
    parent_curr_round.resize(count);									// remove remaining spaces that are not filled in
    
    if(table_subset_curr_round.size()==0){								// If length of table_subset_curr_round is 0 
      break;															// break out of the for-loop,
    }
    for(int k=0;k<table_subset_curr_round.size();k++){					// for-loop of length table_subset_curr_round (potential splits/addittions to trees?)
      overall_trees[overall_count]=table_subset_curr_round[k];		// include all elements of table_subset_curr_round in overall_trees (across iterations of this loop)
      overall_lik[overall_count]=lik_subset_curr_round[k];			// include all elements of lik_subset_curr_round in overall_lik (across iterations of this loop)
      overall_mat[overall_count]=mat_subset_curr_round[k];			// include all elements of mat_subset_curr_round in overall_mat (across iterations of this loop)
      overall_parent[overall_count]=parent_curr_round[k];				// include all elements of parent_curr_round in overall_parent (across iterations of this loop)
      overall_count++;												// increment overall_count by one.
      
      if(overall_count==(overall_size-1)){							// if the length of the overall lists and vectors to be filled is not large enough
        overall_size=overall_size*2;								// double the size
        overall_trees=resize_bigger_bcf(overall_trees,overall_size);	// double the size of overall_trees
        overall_lik.resize(overall_size);							// double the size of overall_lik
        overall_mat=resize_bigger_bcf(overall_mat,overall_size);		// double the size of overall_mat
        overall_parent.resize(overall_size);						// double the size of overall_parent
      }
    }
    overall_trees=resize_bcf(overall_trees,overall_count);					// remove remaining spaces that are not filled in
    overall_lik.resize(overall_count);									// remove remaining spaces that are not filled in
    overall_mat=resize_bcf(overall_mat,overall_count);						// remove remaining spaces that are not filled in
    overall_parent.resize(overall_count);								// remove remaining spaces that are not filled in
    eval_model=evaluate_model_occams_window_bcf(as<NumericVector>(wrap(overall_lik)),lowest_BIC,log(c),overall_trees,overall_mat,as<IntegerVector>(wrap(overall_parent)));	// removes models (trees?) outside Occam's window and returns a list of four elements: tree_list, tree_mat_list, tree_lik, tree_parent
    overall_lik2=eval_model[0];											// vector of all BICS for models in Occam's window. (maybe likelihoods rather than BICs?)
    overall_trees=eval_model[1];										// list of tree tables in Occam's window
    overall_mat=eval_model[2];											// list of tree matrices in Occam's window
    overall_count=overall_trees.size();									// re-set overall_count to number of models in Occam's window
    overall_parent2=eval_model[3];										// tree parent vector for all models in Occam's window.
    //add in check to see if OW accepted more than the top maxOW models...
    if(overall_lik2.size()>maxOWsize){									// If more than maxOWsize models kept in Occam's window
      //find the maxOWsize best models and continue with those!
      IntegerVector owindices=orderforOW__bcf(overall_lik2);					// Function orderforOW__bcf defined on line 555. Gives vector of position of largest element, then position of second largest argument, and so on. 
      owindices=owindices-1;											// take one away from indices so that they are in the correct range.
      //get the top maxOWsize indices to keep in OW
      NumericVector temp_olik(maxOWsize);								// create vector of length maxOWsize
      List temp_otrees(maxOWsize);									// create list of length maxOWsize
      List temp_omat(maxOWsize);										// create list of length maxOWsize
      IntegerVector temp_oparent(maxOWsize);							// create vector of length maxOWsize
      
      //now only select those elements
      for(int t=0;t<maxOWsize;t++){									// for-loop of length equal to size of Occam's window
        temp_olik[t]=overall_lik2[owindices[t]];					// keep the BICs (likelihoods?) for the top maxOWsize models
        temp_otrees[t]=overall_trees[owindices[t]];					// keep the tree tables for the top maxOWsize models
        temp_omat[t]= overall_mat[owindices[t]];					// keep the tree matrices for the top maxOWsize models
        temp_oparent[t]=overall_parent2[owindices[t]];				// keep the tree parent vectors for the top maxOWsize models
      }
      
      overall_lik2=temp_olik;											// reset overall_lik2 to keep only top maxOWsize models
      overall_trees=temp_otrees;										// reset overall_trees to keep only top maxOWsize models
      overall_mat=temp_omat;											// reset overall_mat to keep only top maxOWsize models
      overall_count=overall_trees.size();								// reset overall_count to keep only top maxOWsize models
      overall_parent2=temp_oparent;									// reset overall_parent2 to keep only top maxOWsize models
    }
    
    tree_table_tau=table_subset_curr_round;									// reset tree_table equal to the list of tables produced in the current iteration of the outer loop.
    IntegerVector temp1(table_subset_curr_round.size(),1);				// create a vector of ones of length equal to that of the list of tree tables for the current round.
    err_list=temp1;														// set err_list equal to temp1
    if(overall_trees.size()<overall_size-1){							// if length of overall_trees is less than overall_size-1
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// increse size of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// increse size of overall_mat
      overall_lik.resize(overall_size);								// increse size of overall_lik
      overall_parent.resize(overall_size);							// increse size of overall_parent
    }else{																// if length of overall_trees is greater than or equal to overall_size-1
      overall_size=2*overall_size;									// double the size
      overall_trees=resize_bigger_bcf(overall_trees,overall_size);		// double the length of overall_trees
      overall_mat=resize_bigger_bcf(overall_mat,overall_size);			// double the length of overall_mat
      overall_lik.resize(overall_size);								// double the length of overall_lik
      overall_parent.resize(overall_size);							// double the length of overall_parent
    }
    tree_mat_tau=mat_subset_curr_round;										// reset tree_mat to list of treematrices obtained in current round of for-loop)
    parent=parent_curr_round;											// reset parent to parent vector obtained in current round of for-loop.
    
    if(split_rule_node==1){												// If input boolean split_rule_node = 1 (true)
      NumericVector temp_preds;										// create vector. Length not given.
      List updated_curr_preds;										// create list. Length not given.
      NumericVector new_mean;											// create vector. Length not given.
      lowest_BIC=min(overall_lik2);									// update lowest_BIC
      NumericMatrix curr_resids(resids.nrow(),resids.ncol());			// curr_resids is a matrix with dimensions of input matrix resids
      
      for(int k=0;k<table_subset_curr_round.size();k++){				// for-loop of length equal to number of models proposed in the current iteration
        NumericVector terminal_nodes;								// create a vector called terminal_nodes.
        
        NumericMatrix temp_all_mat = mat_subset_curr_round[k];
        
        arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
        arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
        NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);
        
        if(parent_curr_round[k]==-1){								// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
        }else{														// if the k+1^th element of parent_curr_round does not equal -1 
          new_mean=update_mean_var_bcf(table_subset_curr_round[k],temp_treated_mat,temp_treated_resids(_,parent_curr_round[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that parent_curr_round[k]+1^th column of resids is used as an input vector.
        }  
        
        terminal_nodes=find_term_nodes_bcf(table_subset_curr_round[k]);			// find term nodes function defined line 168.  gives index of values of table_subset_curr_round[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
        updated_curr_preds=update_predictions_bcf(table_subset_curr_round[k],mat_subset_curr_round[k],new_mean,x_moderate_a.n_rows);	// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
        NumericVector test_res;												// create vector called test_res. Length not given.
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          test_res=resids(_,0);											// let test_res equal the first column of resids
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 
          test_res=resids(_,parent_curr_round[k]);						// let test_res equal the parent_curr_round[k]+1^th column of resids
        }
        
        NumericVector temppreds = updated_curr_preds[1];
        NumericVector curr_test_res=z*temppreds;				//might be unnecessary to multiply by z because will probably subset correctly when this vector is used	// Let curr_test_res be the vector of updated predictions
        
        if(parent_curr_round[k]==-1){										// if the k+1^th element of parent_curr_round equals -1 (don't know how it could take a different value)
          curr_resids(_,0)=test_res-curr_test_res;						// let first column of curr_resids be ?the updated residuals?
        }else{																// if the k+1^th element of parent_curr_round does not equal -1 (don't know how it could take a different value)
          curr_resids(_,parent_curr_round[k])=test_res-curr_test_res;		// let parent_curr_round[k]^th column of curr_resids be ?the updated residuals?
        }
      }
      List temp(0);																	// list of length zero.
      
      cp_mat_list=temp;																// reset cp_mat_list to be a list of length zero
      
      for(int f=0;f<curr_resids.ncol();f++){											// loop of length equal to number of columns of curr_resids (equals number of columns of resids at start of loop).
        List cp_mat_list1;															// create a list cp_mat_list1. Length noy given.
        if(gridpoint==0){																// If input Boolean gridpoint equals 0 (false). (i.e. indicator is 1 for grid-search , or 0 for PELT)
          cp_mat_list1=make_pelt_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),pen_tau,num_cp_tau,z);			// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																			// If gridpoint equals 1 (or non-zero).
          cp_mat_list1=make_gridpoint_cpmat_tau_bcf(wrap(x_moderate_a),curr_resids(_,f),gridsize_tau,num_cp_tau,z);	// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        
        cp_mat_list.push_back(cp_mat_list1[0]);											// Add the first element of cp_mat_list1 to the end of the list cp_mat_list. cp_mat_list1[0] is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
      }
      
    }
  }
  // Rcout << "Line 7039.\n";
  
  overall_trees=resize_bcf(overall_trees,overall_count);		// remove remaining spaces that are not filled in
  overall_mat=resize_bcf(overall_mat,overall_count);			// remove remaining spaces that are not filled in
  overall_lik.resize(overall_count);						// remove remaining spaces that are not filled in
  overall_parent.resize(overall_count);					// remove remaining spaces that are not filled in
  NumericVector temp_preds;								// create a vector
  List updated_preds;										// create a list
  NumericVector new_mean;									// create a vector
  NumericMatrix overall_test_preds(x_moderate_test.nrow(),overall_trees.size());	// create a matrix with no. of rows equal to that of input x_moderate_test, and no. of columns equal to the number of elements of overall_trees
  NumericMatrix overallpreds(x_moderate_a.n_rows,overall_trees.size());		// create a matrix with no. of rows equal to that of x_moderate_a (an input data matrix), and no. of columns equal to the number of elements of overall_trees
  lowest_BIC=min(overall_lik2);							// reset lowest_BIC to minimum of overall_lik2
  
  
  for(int k=0;k<overall_trees.size();k++){			// loop of length equal to that of list overall_trees
    NumericVector terminal_nodes;					// create a vector
    
    
    NumericMatrix temp_all_mat = overall_mat[k];
    
    arma::mat temp_all_mat_a(temp_all_mat.begin(), temp_all_mat.nrow(), temp_all_mat.ncol(), false);	
    arma::mat temp_treated_mat_a = temp_all_mat_a.rows(arma::find(z_ar==1));	
    NumericMatrix temp_treated_mat = wrap(temp_treated_mat_a);	
    
    if(overall_parent2[k]==-1){						// If k+1^th element of overall_parent2 equals -1
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,0),a_tau);	// return vector of new means for each terminal node. Defined on line 1285. Note that first column of resids is used as an input vector.
    }else{    
      new_mean=update_mean_var_bcf(overall_trees[k],temp_treated_mat,temp_treated_resids(_,overall_parent2[k]),a_tau);		// return vector of new means for each terminal node. Defined on line 1285. Note that overall_parent2[k]+1^th column of resids is used as an input vector.
    }  
    terminal_nodes=find_term_nodes_bcf(overall_trees[k]);		// find term nodes function defined line 168. Gives index of values of overall_trees[k] that are term nodes (indices from 1 to length of vector). Why not integer vector?
    updated_preds=update_predictions_bcf(overall_trees[k],overall_mat[k],new_mean,x_moderate_a.n_rows);		// Defined in line 1313. Returns list. First element is updated tree table (6th column contains new predictions). Second element is vector of updated predictions.
    //get the predicted values for the test data.
    if(is_test_data) test_preds=get_testdata_term_obs_bcf(x_moderate_test,overall_trees[k],new_mean);		// If input boolean is_test_data is true, return the vector of predictions for tree overall_trees[k].
    temp_preds=updated_preds[1];																// (re)set temp_preds equal to the vector of updated predictions
    overallpreds(_,k)=temp_preds;																// Let the k+1^th column of overallpreds be temp_preds
    if(is_test_data)overall_test_preds(_,k)=test_preds;											// If input boolean is_test_data is true, let the k+1^th column of overallpreds be the vector of predictions for tree overall_trees[k].
  }
  arma::mat M1(overallpreds.begin(), overallpreds.nrow(), overallpreds.ncol(), false);			// create arma mat copy of overallpreds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_values=sum(M1,1);														// vector of the sum of all elements in each row of overallpreds. (vetor of final predictions)?
  arma::mat M2(overall_test_preds.begin(), overall_test_preds.nrow(), overall_test_preds.ncol(), false);		// create arma mat copy of overall_test_preds [each column appears to be the vector of predictions for a particular tree (in the sum of trees)]
  arma::colvec predicted_test_values=sum(M2,1);		// vector of the sum of all elements in each row of overall_test_preds. (vetcor of final predictions?
  List ret(7);				// list of length 8.
  ret[0]=overall_lik2;		// first element of list is a vector of BICs?
  ret[1]=overall_trees;		// list of tree tables.
  ret[2]=overall_mat;			// list of tree matrices.
  //ret[3]=predicted_values;	// vector of (in-sample) predicted values.
  ret[3]=overall_parent2;		// vector of tree parent numbers (not sure what possible values this can take. Could every element be -1?)
  ret[4]=wrap(M1);			// (in-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  ret[5]=lowest_BIC;			// lowest BIC among trees
  ret[6]=wrap(M2);			// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspons to predictions from different (single) trees.
  return(ret);				// return the list
}
//######################################################################################################################//
// [[Rcpp::export]]
NumericVector scale_response_bcf(double a,double b,double c,double d,NumericVector y){
  NumericVector y_scaled = -((-b*c+a*d)/(-a+b))+((-c+d)*y/(-a+b));
  
  return(y_scaled);
}
//######################################################################################################################//

// [[Rcpp::export]]
NumericVector get_original_bcf(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
  NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
  
  return(original_y); // reverse scaling of predictions of scaled variable (??)
}
//###########################################################################################################################//

//' @title Obtain BCFBMA predictions, trees, BICs etc. to be called by R functions
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List BCF_BMA_sumLikelihood(NumericMatrix data,NumericVector y, NumericVector z, NumericMatrix pihat,
                           double a_mu,double a_tau,double mu_mu,double mu_tau,double nu,double lambda,double c,
                           double sigma_mu_mu,double sigma_mu_tau,double pen_mu,double pen_tau,int num_cp_mu,int num_cp_tau,
                           NumericMatrix test_data,NumericVector test_z,NumericMatrix test_pihat,
                           int ntree_control,int ntree_moderate,
                           double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,bool split_rule_node,bool gridpoint,int maxOWsize,
                           int num_splits_mu,int num_splits_tau,int gridsize_mu, int gridsize_tau, int include_pi2, bool zero_split, bool only_max_num_trees){
  bool is_test_data=0;					// create bool is_test_data. Initialize equal to 0.
  if(test_data.nrow()>0){					// If test data has non-zero number of rows.
    is_test_data=1;						// set is_test_data equal to 1.
  }
  if(y.size() !=data.nrow()){				// If the length of input vector y is not equal to the nunber of rows in the input data (covariates)
    if(y.size()<data.nrow()){			// If the length of y is less than the number of rows in data
      throw std::range_error("Response length is smaller than the number of observations in the data"); 
    }else{								// If the length of y is greater than the number of rows in data
      throw std::range_error("Response length is greater than the number of observations in the data"); 
    }
  }
  if(z.size() !=data.nrow()){				// If the length of input vector z is not equal to the nunber of rows in the input data (covariates)
    if(z.size()<data.nrow()){			// If the length of z is less than the number of rows in data
      throw std::range_error("Treatment indicator vector length is smaller than the number of observations in the data"); 
    }else{								// If the length of z is greater than the number of rows in data
      throw std::range_error("Treatment indicator vector length is greater than the number of observations in the data"); 
    }
  }
  if(pihat.nrow() !=data.nrow()){				// If the nunber of rows in the input matrix pihat is not equal to the nunber of rows in the input data (covariates)
    if(pihat.nrow()<data.nrow()){			// If the nunber of rows in the input matrix pihat is less than the number of rows in data
      throw std::range_error("The nunber of rows in the input matrix pihat is smaller than the number of observations in the data"); 
    }else{								// If the nunber of rows in the input matrix pihat is greater than the number of rows in data
      throw std::range_error("The nunber of rows in the input matrix pihat is greater than the number of observations in the data"); 
    }
  }
  //check test data has the same number of variables as training data
  if(test_data.nrow()>0 && (data.ncol() != test_data.ncol())){	// If the number of rows in the test data is >0 AND the number of columns (variables) is not equal to that of data (the training data)
    throw std::range_error("Test data and training data must have the same number of variables. BART BMA assumes variables are in the same order."); 
  }
  if(test_z.size() != test_data.nrow()){	// If the number of rows in the test data covariate matrix is not equal to that of the test data treatment indicator variable
    throw std::range_error("Test data covariates and test data treatment indicator variable must have the same number of observations."); 
  }
  if(test_data.nrow() != test_pihat.nrow()){	// If the number of rows in the test data covariate matrix is not equal to that of the test data propensity score estimates matrix
    throw std::range_error("Test data covariates and test data propensity score estimates must have the same number of observations."); 
  }
  if(test_pihat.nrow()>0 && (pihat.ncol() != test_pihat.ncol())){	// If the number of rows in the test data propensity score estimates is >0 AND the number of columns (variables) is not equal to that of the training data propensity score estimates
    throw std::range_error("Test data propensity score estimates and training data propensity score estimates must have the same number of columns. BART BMA assumes variables are in the same order."); 
  }
  //check value of c is greater than 1!
  //	if(c<1){
  //		throw std::range_error("Value of Occam's Window has to be greater than 0."); 
  //	}
  if(num_cp_mu<0 || num_cp_mu>100){		// If input num_cp_mu is <0 or >100
    throw std::range_error("Value of num_cp_mu should be a value between 1 and 100."); 
  }
  if(num_cp_tau<0 || num_cp_tau>100){		// If input num_cp_tau is <0 or >100
    throw std::range_error("Value of num_cp_tau should be a value between 1 and 100."); 
  }
  // Now add propensity score estimates matrix as new leftmost column of data matrix. Call the resulting matrix x_control (to be consistent with terminology used by bcf package).
  arma::mat D1(data.begin(), data.nrow(), data.ncol(), false);				// copy the covariate data matrix into an arma mat
  arma::mat pihat_a(pihat.begin(), pihat.nrow(), pihat.ncol(), false);				// copy the pihat matrix into an arma mat
  arma::mat x_control_a=D1;				// create a copy of data arma mat called x_control_a
  if((include_pi2==0) | (include_pi2==2) ){
  if(pihat.nrow()>0 ){
  x_control_a.insert_cols(0,pihat_a);		// add propensity scores as new leftmost columns of x_control_a
  }
  }
  //Rcout << "Number of columns of matrix" << x_control_a.n_cols << ".\n";
  NumericMatrix x_control=wrap(x_control_a);	// convert x_control_a to a NumericMatrix called x_control
    // Name the matrix without the estimated propensity scores x_moderate.[CAN REMOVE THE DUPLICATION AND ADD x_control, x_moderate, and include_pi as input parameters later]
    //NumericMatrix x_moderate = data;	// x_moderate matrix is the covariate data without the propensity scores
  arma::mat x_moderate_a=D1;			// create arma mat copy of x_moderate.
  if((include_pi2==1)| (include_pi2==2) ){
    if(pihat.nrow()>0 ){
      x_moderate_a.insert_cols(0,pihat_a);		// add propensity scores as new leftmost columns of x_control_a
    }
  }
  NumericMatrix x_moderate=wrap(x_moderate_a);	// convert x_control_a to a NumericMatrix called x_control
  //Rcout << "Get to Line 5001  "  << ".\n";
  // Add test propensity scores to test data matrix
  arma::mat T1(test_data.begin(), test_data.nrow(), test_data.ncol(), false);				// copy the covariate test_data matrix into an arma mat
  arma::mat pihat_a_test(test_pihat.begin(), test_pihat.nrow(), test_pihat.ncol(), false);				// copy the test_pihat matrix into an arma mat
  arma::mat x_control_test_a=T1;				// create a copy of test_data arma mat called x_control_test_a
  if((include_pi2==0)| (include_pi2==2) ){
  if(test_pihat.nrow()>0 ){
  x_control_test_a.insert_cols(0,pihat_a_test);		// add propensity scores as new leftmost columns of x_control_test_a
  }
  }
  NumericMatrix x_control_test=wrap(x_control_test_a);	// convert x_control_test_a to a NumericMatrix called x_control_test
    // Name the matrix without the estimated propensity scores x_moderate_test.[CAN REMOVE THE DUPLICATION AND ADD x_control_test, x_moderate_test, and include_pi as input parameters later]
    //NumericMatrix x_moderate_test = test_data;	// x_moderate_test matrix is the covariate test_data without the propensity scores
  arma::mat x_moderate_test_a=T1;			// create arma mat copy of x_moderate_test.
  if((include_pi2==1)| (include_pi2==2) ){
    if(test_pihat.nrow()>0 ){
      x_moderate_test_a.insert_cols(0,pihat_a_test);		// add propensity scores as new leftmost columns of x_control_a
    }
  }
  NumericMatrix x_moderate_test=wrap(x_moderate_test_a);	// convert x_control_test_a to a NumericMatrix called x_control_test
  
  //Rcout << "Get to Line 5022  "  << ".\n";
  
  // NOT SURE IF SEPARATE INITIAL TREE MATRIX REQUIRED FOR mu(x) and tau(x)
  // BUT STILL DESIRABLE TO END UP WITH SEPARATE LISTS AND MATRICES FOR mu(x) and tau(x) trees
  //
  
  
  //	NumericMatrix treetable=start_tree_bcf(mu,sigma_mu);								// create matrix treetable (defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
  //	NumericMatrix treemat=start_matrix_bcf(data.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  //initialize the tree table and matrix
  
  NumericMatrix treetable_mu=start_tree_bcf(mu_mu,sigma_mu_mu);								// create matrix treetable (defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
  NumericMatrix treemat_mu=start_matrix_bcf(x_control.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  
  //	MIGHT NEED JUST NUMBER OF TREATED OBSERVATIONS for treemat_tau
  NumericMatrix treetable_tau=start_tree_bcf(mu_tau,sigma_mu_tau);								// create matrix treetable (defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
  
  // 	
  
  // create matrix of data for treated households only (not sure if this will be used)
  //	arma::mat x_mod_treat_a = x_moderate_a.rows(arma::find(z_ar==1));	
  //	NumericMatrix x_mod_treat = wrap(x_mod_treat_a);	
  //	NumericMatrix treemat_tau=start_matrix_bcf(x_mod_treat.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  
  // FOR NOW TRY TO PUT ALL OBSERVATIONS IN THE TAU TREE MATRIX OUTPUT
  NumericMatrix treemat_tau=start_matrix_bcf(x_moderate.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  
  
  //	MAYBE INCLUDE treetable_tau and treemat_tau HERE
  
  NumericVector y_scaled=scale_response_bcf(min(y),max(y),-0.5,0.5,y);		// re-scale the outcome variable
  double n=D1.n_rows;					// number of rows of the data matrix (number of training observations)
  //	CHANGE ABOVE LINE IF REMOVE data and D1
  
  // BELOW FUNCTIONS AND LINES ONLY USED DIRECTLY HERE FOR OBTAINING lowest_BIC
  //	//	double lik=likelihood_function_bcf(y_scaled,treetable,treemat,a,mu,nu,lambda);
  //	double lik=likelihood_function_bcf(y_scaled,treetable,treemat,a,mu,nu,lambda);	// defined on line 201. Returns log likelihood of initial model with no trees?
  //	double tree_prior=get_tree_prior_bcf(treetable,treemat,alpha,beta);				// defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees? but the tree is empty  at this stage?)
  //	double lowest_BIC=-2*(lik+log(tree_prior))+1*log(n);						// BIC actually not the BIC because add 2*log(prior). Closer to the Bayes factor? (Not noted in paper).
  
  //	THIS IS PROBABLY WRONG: SHOULD START WITH MODEL Y=a+bZ, WHERE THERE IS A mu(x) stump AND tau(x) STUMP.]
  // 	currently: will add tau stump after first mu(x) tree grown
  double lik=likelihood_function_mu_bcf(y_scaled,treetable_mu,treemat_mu,a_mu,mu_mu,nu,lambda);	// defined on line 201. Returns log likelihood of initial model with no trees?
  if(treetable_mu.ncol()<5) throw std::range_error("Line 4081");
  double tree_prior=get_tree_prior_bcf(treetable_mu,treemat_mu,alpha_mu,beta_mu);				// defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees? but the tree is empty  at this stage?)
  double lowest_BIC=-2*(lik+log(tree_prior))+1*log(n);						// BIC actually not the BIC because add 2*log(prior). Closer to the Bayes factor? (Not noted in paper).

  // create initial tree table lists for mu(x) and tau(x)
  // Need to check if doing this correctly for tau(x)
  List tree_table_mu;								// create tree table.
  List tree_mat_mu;									// create tree matrix.
  tree_table_mu.push_back(treetable_mu);				// append treetable to list tree_table. (adding the model with no trees as the first element of the list).
  tree_mat_mu.push_back(treemat_mu);					// append treemat to list tree_mat. (adding the model with no trees as the first element of the list).
  
  List tree_table_tau;								// create tree table.
  List tree_mat_tau;									// create tree matrix.
  // these two lines might be unnecessary
  tree_table_tau.push_back(treetable_tau);				// append treetable to list tree_table. (adding the model with no trees as the first element of the list).
  tree_mat_tau.push_back(treemat_tau);					// append treemat to list tree_mat. (adding the model with no trees as the first element of the list).
  
  
  //	Not sure if should create separate CART_BMA_mu and CART_BMA_tau lists or just write over when moving from mu trees to tau trees
  // Probably unnecessary to have separate BART_BMA, but can edit this later
  List CART_BMA_mu;									// create list.
  List CART_BMA_tau;									// create list.
  
  
  arma::mat r;									// create matrix r.
  arma::colvec yarma=clone(y_scaled);				// create arma colvec copy, yarma, of scaled outcome vector.
  r.insert_cols(0,yarma);							// let the first column of r be yarma (the scaled outcome vector).
  NumericMatrix resids=wrap(r);					// create NumericMatrix copy of r called resids.
  
  
  // FOR NOW, TRYING TO LEAVE ALL SUBSETTING TO TREATED RESIDUALS TO FUNCTIONS (e.g. get_best_trees_sum_tau_bcf)	
  // create vector of treated outcome observations only?
  //	arma::mat yarma_treated = yarma.elem(arma::find(z_ar==1));	
  //	create vector of treated resids only?
  //	perhaps unnecessary if first should remove prediction from first mu tree
  //	arma::mat r_treated;									// create matrix r.
  //	r_treated.insert_cols(0,yarma_treated);					// let the first column of r be yarma (the scaled outcome vector).
  //	NumericMatrix resids=wrap(r_treated);					// create NumericMatrix copy of r called resids.
  
  
  
  //int first_round;								// create a variable first_round. Not initialized.
  
  //	NOT SURE IF SHOULD KEEP ALL OF THESE
  //	List overall_trees(max(ntree_control,ntree_moderate));					// create a list of length equal to the input value num_rounds.
  //	List overall_mat;								// create a list.
  List overall_lik;								// create a list.
  
  List overall_trees_mu(ntree_control);					// create a list of length equal to the input value num_rounds.
  List overall_mat_mu;								// create a list.
  //	List overall_lik_mu;								// create a list.
  
  List overall_trees_tau(ntree_moderate);					// create a list of length equal to the input value num_rounds.
  List overall_mat_tau;								// create a list.
  //	List overall_lik_tau;								// create a list.
  
  
  NumericMatrix prev_round_preds_outcome;  				// create a matrix.
  NumericMatrix prev_round_preds_mu;  				// create a matrix.
  NumericMatrix prev_round_preds_tau;  				// create a matrix.
  
  NumericVector prev_round_BIC;					// create a vector.
  NumericVector prev_round_BIC2;					// create a vector.
  arma::mat prev_round_preds2_outcome;					// create an arma mat.
  arma::mat prev_round_preds2_mu;					// create an arma mat.
  arma::mat prev_round_preds2_tau;					// create an arma mat.
  NumericMatrix prev_round_test_preds_outcome;			// create a matrix. (default values are zeroes?)
  NumericMatrix prev_round_test_preds_mu;			// create a matrix. (default values are zeroes?)
  NumericMatrix prev_round_test_preds_tau;			// create a matrix. (default values are zeroes?)
  
  arma::mat prev_round_test_preds2_outcome;				// create an arma mat
  arma::mat prev_round_test_preds2_mu;				// create an arma mat
  arma::mat prev_round_test_preds2_tau;				// create an arma mat
  
  arma::mat overall_overall_sum_test_preds_outcome;		// create an arma mat
  arma::mat overall_overall_sum_test_preds_mu;		// create an arma mat
  arma::mat overall_overall_sum_test_preds_tau;		// create an arma mat
  
  arma::colvec predicted_test_values_outcome;				// create an arma colvec
  arma::colvec predicted_test_values_mu;				// create an arma colvec
  arma::colvec predicted_test_values_tau;				// create an arma colvec
  
  
  //	Not sure if should keep these, or separately have lists for mu and tau trees, resids, mat	
  //	List prev_sum_trees;							// create a list
  // List prev_sum_tree_resids;						// create a list
  //	List prev_sum_trees_mat;						// create a list
  //	List cp_mat_list;								// create a list
  
  List prev_sum_trees_mu;								// create a list
  List prev_sum_tree_resids_mu;						// create a list
  List prev_sum_trees_mat_mu;							// create a list
  List cp_mat_list_mu;								// create a list
  
  List prev_sum_trees_tau;							// create a list
  List prev_sum_tree_resids_tau;						// create a list
  List prev_sum_trees_mat_tau;						// create a list
  List cp_mat_list_tau;								// create a list	
  
  
  
  int oo_size=300;								// create a variable. Initialized equal to 300.
  List overall_overall_sum_tree_resids_mu(oo_size);	// create a list of length 300.
  List overall_overall_sum_tree_resids_tau(oo_size);	// create a list of length 300.
  
  List overall_overall_sum_BIC(oo_size);			// create a list of length 300.
  int oo_count=0;									// create a variable, Initialized equal to 0.
  
  List overall_overall_sum_trees_mu(oo_size);			// create a list of length 300.
  List overall_overall_sum_trees_mat_mu(oo_size);		// create a list of length 300.
  
  List overall_overall_sum_trees_tau(oo_size);		// create a list of length 300.
  List overall_overall_sum_trees_mat_tau(oo_size);	// create a list of length 300.
  
  arma::mat overall_overall_sum_preds_outcome;			// create an arma mat.
  arma::mat overall_overall_sum_preds_mu;			// create an arma mat.
  arma::mat overall_overall_sum_preds_tau;			// create an arma mat.
  
  //	parent indexes a whole sum-of-tree model including mu(x) and tau(x) trees
  IntegerVector prev_par;							// create a vector
  
  arma::colvec predicted_values_outcome;					// create an arma colvec
  arma::colvec predicted_values_mu;					// create an arma colvec
  arma::colvec predicted_values_tau;					// create an arma colvec
  
  
  for(int j=0;j<max(ntree_control,ntree_moderate);j++){					// create a for-loop of length equal to the input value num_rounds.
    // Rcout << "Beginning of loop number = " << j << ".\n";
    // Rcout << "ntree_control = " << ntree_control << ".\n";
    // Rcout << "ntree_moderate = " << ntree_moderate << ".\n";
    
    if(j<ntree_control){
      if(j>0){
        if(only_max_num_trees==1){
          lowest_BIC=100000;
        }
        }
      int overall_size=300;						// create a variable overall_size. Initialized equal to 0.
      List overall_sum_tree_resids_mu(overall_size);	// create a list of length 300
      List overall_sum_tree_resids_tau(overall_size);	// create a list of length 300
      
      List overall_sum_trees_mu(overall_size);		// create a list of length 300
      List overall_sum_trees_mat_mu(overall_size);	// create a list of length 300
      
      List overall_sum_trees_tau(overall_size);		// create a list of length 300
      List overall_sum_trees_mat_tau(overall_size);	// create a list of length 300
      
      // Rcout << "Get to Line 5212 in loop j = " << j << ".\n";
      int overall_count=0;						// set overall_count equal to 0.
      //		parent indexes the whole models to which the a tree can be appended
      IntegerVector parent;						// create vector.
      NumericVector curr_round_lik;				// create vector	// To be filled with BICs for whole models suggested after a tree appended
      List curr_round_trees_mu;						// create list.
      List curr_round_trees_tau;						// create list.
      
      List curr_round_mat_mu;						// create list.
      List curr_round_mat_tau;						// create list.
      
      NumericVector curr_BIC;						// create vector.
      //		next line probably shouldn't need to be duplicated for curr_round_parent_mu and curr_round_parent_tau
      IntegerVector curr_round_parent;			// create vector.
      NumericVector overall_sum_BIC;				// create vector.
      
      arma::mat overall_sum_preds_outcome;				// create matrix.
      arma::mat overall_sum_preds_mu;				// create matrix.
      arma::mat overall_sum_preds_tau;				// create matrix.
      
      arma::mat overall_sum_test_preds_outcome;			// create matrix.
      arma::mat overall_sum_test_preds_mu;			// create matrix.
      arma::mat overall_sum_test_preds_tau;			// create matrix.
      //Rcout << "Get to Line 5235 in loop j = " << j << ".\n";
      if(j==0){									// If in the first round of the for-loop.
        parent.push_back(0);					// append a 0 to the end of the parent vector. (first and only element of parent vector so far).
        //first_round=1;							// set the variable first_round equal to 1.
      }else{										// If not in the first round of the for-loop.
        //first_round=0;							// set the variable first_round equal to 0.
      }
      //		The _mu	in resids_cp_mat_mu is probably unnecessary, but including it to remove ambiguity.
      //		Replace with just resids_cp_mat for memory efficiency after code is all working.
      //		similarly, cp_mat_list_mu is used (perhaps unnecessarily) here instead of cp_mat_list
      List resids_cp_mat_mu(resids.ncol());												// create a list of length equal to the number of colmns of resids
      int resids_count=0;																// create a variable resids_count. Initialize equal to zero.
      std::vector<int> err_list(resids.ncol());										// create a vector err_list of length equal to the number of colmns of resids
      //get best splits
      for(int f=0;f<resids.ncol();f++){												// for-loop of length equal to the unmber of columns of resids
        if(gridpoint==0){															// If input gridpoint equals 0. i.e. the PELT method will be used.
          cp_mat_list_mu=make_pelt_cpmat_mu_bcf(x_control,resids(_,f),pen_mu,num_cp_mu);				// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																		// If input gridpoint equals 1. i.e. the PELT method will be used.
          cp_mat_list_mu=make_gridpoint_cpmat_mu_bcf(x_control,resids(_,f),gridsize_mu,num_cp_mu);			// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        resids_cp_mat_mu[resids_count]=cp_mat_list_mu[0];									// let the resids_count+1^th element be a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. 
        err_list[resids_count]=cp_mat_list_mu[1];										// let the resids_count+1^th element be a list that records error(s?)
        resids_count++;																// increment resids_count by 1.
      }
      resids_cp_mat_mu=resize_bcf(resids_cp_mat_mu,resids_count);								// remove elements of resids_cp_mat_mu that are not fileld in. How is  it possible that elements are not filled in??
      err_list.resize(resids_count);													// remove elements of err_list that are not fileld in. How is  it possible that elements are not filled in??
      parent=seq_len(tree_table_mu.size())-1;											// set parent equal to a vector 0,1,2,3,..., up to the length of tree_table minus one.
      if(is_true(all(as<IntegerVector>(wrap(err_list))==1))){							// if all elements of err_list equal 1.
        if(j==0){																	// If in the first round of the for-loop.
          throw std::range_error("No split points could be found to grow trees");
        }else{																		// If not in the first round of the for-loop.
          throw std::range_error("No Mu trees can be grown for the number of iterations desired, as no splits were found.Please try fewer iterations.");
        }
      } 
      //Rcout << "Get to Line 5269 in loop j = " << j << ".\n";
      //get current set of trees.
      if(j==0){						// If in the first round of the for-loop.
        CART_BMA_mu=get_best_trees_mu_bcf(x_control_a, x_moderate_a,z,resids,
                                      a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                      tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                      lowest_BIC,//first_round,
                                      parent,resids_cp_mat_mu,as<IntegerVector>(wrap(err_list)),
                                      x_control_test,x_moderate_test,test_z,
                                      alpha_mu,alpha_tau,beta_mu,beta_tau,
                                      is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                      split_rule_node,gridpoint,maxOWsize,num_splits_mu,num_splits_tau,gridsize_mu,zero_split);
        
      }else{							// If not in the first round of the for-loop.
        //if j >0 then sum of trees become a list so need to read in list and get likelihood for each split point and terminal node

        CART_BMA_mu=get_best_trees_sum_mu_bcf(x_control_a, x_moderate_a,z,resids,
                                          a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                          tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                          lowest_BIC,//first_round,
                                          parent,resids_cp_mat_mu,as<IntegerVector>(wrap(err_list)),
                                          x_control_test,x_moderate_test,test_z,
                                          alpha_mu,alpha_tau,beta_mu,beta_tau,
                                          is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                          split_rule_node,gridpoint,maxOWsize,
                                          prev_sum_trees_mu,prev_sum_trees_tau,prev_sum_trees_mat_mu,
                                          prev_sum_trees_mat_tau,y_scaled,num_splits_mu,num_splits_tau,gridsize_mu,zero_split);	// function defined on line 1953.
      }
       //Rcout << "Get to after get best trees in mu round in loop j = " << j << ".\n";
      
      curr_round_lik=CART_BMA_mu[0];							// vector of BICs (for whole sum-of tree-models after suggested trees added). Should be ordered ascending
      curr_round_trees_mu=CART_BMA_mu[1];						// list of tree tables
      //curr_round_trees_tau=CART_BMA_mu[2];
      curr_round_mat_mu=CART_BMA_mu[2];							// list of tree matrices
      //curr_round_mat_tau=CART_BMA_mu[4];							// list of tree matrices
      curr_round_parent=CART_BMA_mu[3];						// vector of tree parent numbers
      NumericMatrix curr_round_preds_mu=CART_BMA_mu[4];			// (in-sample single tree predictions for trees to add) matrix rows correspond to different units/individuals, columns corresponds to predictions from different (single or sums-of?) trees.
      curr_BIC=CART_BMA_mu[5];								// lowest BIC among trees?
      NumericMatrix curr_round_test_preds_mu=CART_BMA_mu[6];	// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspond to predictions from different (single or sums-of?) trees.
      if(curr_round_lik.size()==0) {						// If number of sum of tree models is zero?
        // Rcout << "curr_round_lik.size()==0 BREAK in mu round in loop j = " << j << ".\n";
        //REMOVE THIS ERROR IF WANT TO ALLOW LESS THAN MAX NUMBER OF TREES
        //throw std::range_error("No mu trees chosen in round");
        
        break;											// break out of for-loop
      } 
      //Rcout << "Get to Line 5313 in loop j = " << j << ".\n";
      if(curr_BIC[0]<lowest_BIC){							// If the lowest BIC obtained by get_best_trees_sum is less than the currently saved lowest value
        lowest_BIC=curr_BIC[0];							// reset lowest_BIC to the new lowest value
      }
      tree_table_mu=List();									// reset tree_table to an empty list
      tree_mat_mu=List();									// reset tree_mat to an empty list.
      tree_table_tau=List();									// reset tree_table to an empty list
      tree_mat_tau=List();									// reset tree_mat to an empty list.
      int lsize=curr_round_lik.size();					// create a variable equal to the number of sum of tree models returned by get_best_trees_sum
      tree_table_mu=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_mat_mu=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_table_tau=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_mat_tau=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_outcome(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_mu(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_tau(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      
      NumericMatrix temp_test_preds_outcome(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_test_preds_mu(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_test_preds_tau(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      
      NumericMatrix temp_resids(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericVector temp_parent(curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericVector temp_BIC(curr_round_lik.size());		// create a vector of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      List temp_sum_trees_mu(lsize);							// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_trees_tau(lsize);							// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      List temp_sum_tree_resids(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_tree_resids_mu(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_tree_resids_tau(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      List temp_sum_trees_mat_mu(lsize);						// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_trees_mat_tau(lsize);						// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
       //Rcout << "Get to Line 5347 in loop j = " << j << ".\n";
      
      int count=0; 										// create a count variable. Initialized equal to zero.
      for(int k=0;k<curr_round_lik.size();k++){			// create a for-loop of length equal to the number of sum of tree models returned by get_best_trees_sum
        tree_table_mu[count]=start_tree_bcf(mu_mu,sigma_mu_mu);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
        tree_mat_mu[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
        tree_table_tau[count]=start_tree_bcf(mu_tau,sigma_mu_tau);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
        tree_mat_tau[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
        if(j==0){										// if in first round of for-loop
          temp_preds_outcome(_,k)=curr_round_preds_mu(_,k);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
          temp_preds_mu(_,k)=curr_round_preds_mu(_,k);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
          NumericVector zerovec1(n,0);
          temp_preds_tau(_,k)=zerovec1;
          if(is_test_data==1){       
            //Rcout << "Get to Line 5362 in loop j = " << j << "k= "<< k << ".\n";
            temp_test_preds_outcome(_,k)=curr_round_test_preds_mu(_,k);
            temp_test_preds_mu(_,k)=curr_round_test_preds_mu(_,k);
            NumericVector zerovectest(test_data.nrow(), 0.0);
            temp_test_preds_tau(_,k)=zerovectest;
          }
          // If there is test data, let the k+1^th column of temp_test_preds be the k+1^th column of curr_round_test_preds. These are the out-of-sample predictions of from the k+1^th model.
          temp_resids(_,k)=y_scaled-temp_preds_outcome(_,k);	// Let the k+1^th column of temp_resids be the outcome minus the predictons from the k+1^th model 
          temp_parent[k]=-1;							// Let the K=1^th element of temp_parent be -1.
          temp_BIC[k]=curr_round_lik[k];				// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
          temp_sum_trees_mu[count]=curr_round_trees_mu[k];	// Add the tree table to the list temp_sum_trees
          //temp_sum_trees_tau[count]=curr_round_trees_tau[k];	// Add the tree table to the list temp_sum_trees
          
          temp_sum_trees_mat_mu[count]=curr_round_mat_mu[k];	// Add the tree matrix to temp_sum_trees_mat
          //temp_sum_trees_mat_tau[count]=curr_round_mat_tau[k];	// Add the tree matrix to temp_sum_trees_mat
          //NOT SURE ABOUT THE FOLLOWING LINE
          temp_sum_tree_resids[count]=resids(_,0);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
          temp_sum_tree_resids_mu[count]=resids(_,0);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
          // Rcout << "LENGTH OF resids(_,0)= " << resids(_,0).size() << ".\n";
        }else{											// If not in the first round of the for-loop.
          NumericVector curr_temp_pred_outcome=curr_round_preds_mu(_,k) + 
            prev_round_preds_mu(_,curr_round_parent[k])+ 
            z*prev_round_preds_tau(_,curr_round_parent[k]);	// curr_temp_pred is the sum of the current round predictions and the previous round predictions?? Each round is for one tree? and explain more of the residuals in each round?
          NumericVector curr_temp_pred_mu=curr_round_preds_mu(_,k) + 
            prev_round_preds_mu(_,curr_round_parent[k]);					
          NumericVector curr_temp_pred_tau=prev_round_preds_tau(_,curr_round_parent[k]);	
          
          NumericVector curr_temp_test_pred_outcome;			// create a vector
          NumericVector curr_temp_test_pred_mu;			// create a vector
          NumericVector curr_temp_test_pred_tau;			// create a vector
          
          if(is_test_data==1) {						// If there is test data.
            //Rcout << "Get to Line 5393 in loop j = " << j << "k= "<< k << ".\n";
            curr_temp_test_pred_outcome=curr_round_test_preds_mu(_,k) + prev_round_test_preds_mu(_,curr_round_parent[k])+test_z*prev_round_test_preds_tau(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            curr_temp_test_pred_mu=curr_round_test_preds_mu(_,k) + prev_round_test_preds_mu(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            curr_temp_test_pred_tau=prev_round_test_preds_tau(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            
            temp_test_preds_outcome(_,k) = curr_temp_test_pred_outcome;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            temp_test_preds_mu(_,k) = curr_temp_test_pred_mu;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            temp_test_preds_tau(_,k) = curr_temp_test_pred_tau;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            //Rcout << "Get to Line 5401 in loop j = " << j << "k= "<< k << ".\n";
          }
          temp_BIC[k]=curr_round_lik[k];									// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
          temp_preds_outcome(_,k)=curr_temp_pred_outcome;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          temp_preds_mu(_,k)=curr_temp_pred_mu;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          temp_preds_tau(_,k)=curr_temp_pred_tau;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          
          temp_resids(_,k)=y_scaled-curr_temp_pred_outcome;						// Let the k+1^th column of temp_resids be the new residuals after appending the k+1^th
          temp_parent[k] = k;												// Let k+1^th element of temp_parent equal k.
          temp_sum_trees_mu[count]=curr_round_trees_mu[k];						// Add the tree table to the list temp_sum_trees.
          //temp_sum_trees_tau[count]=curr_round_trees_tau[k];						// Add the tree table to the list temp_sum_trees.
          
          temp_sum_trees_mat_mu[count]=curr_round_mat_mu[k];					// Add the tree matrix to temp_sum_trees_mat.
          //temp_sum_trees_mat_tau[count]=curr_round_mat_tau[k];					// Add the tree matrix to temp_sum_trees_mat.
          temp_sum_tree_resids[count]=resids(_,curr_round_parent[k]);		// Add the curr_round_parent[k]+1^th column of resids to temp_sum_tree_resids. resids contains predictions from previous rounds?
          temp_sum_tree_resids_mu[count]=resids(_,curr_round_parent[k]);		// Add the curr_round_parent[k]+1^th column of resids to temp_sum_tree_resids. resids contains predictions from previous rounds?
          
        }
        count++;													// Increment the count by 1.(Note this is within the innermost for-loop).
      }  
      if(curr_round_lik.size()==0){									// if no new trees outputted in current round (by get_best_trees_sum? Why not throw this error earlier, at line 2350)
        throw std::range_error("No trees chosen in last round");
      }
      // Rcout << "Get to line 5422 in loop j = " << j << ".\n";
      for(int k=0;k<curr_round_lik.size();k++){	// create a for-loop of length equal to the number of sum of tree models returned by get_best_trees_sum
        int size_mat=300;						// create a variable, Initializd equal to 300.
        List sum_of_trees_mu(size_mat);			// create a list of length 300.
        List sum_of_tree_resids_mu(size_mat);		// create a list of length 300.
        List sum_of_tree_resids_tau(size_mat);		// create a list of length 300.
        List sum_of_trees_mat_mu(size_mat);		// create a list of length 300.
        int count=0;							// create a variable count equal to 0. (Count was already define, so could remove "int" at start of this line and just reset cound to 0).
        
        // Rcout << "Get to line 5431 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        List sum_of_trees_tau(size_mat);			// create a list of length 300.
        List sum_of_trees_mat_tau(size_mat);		// create a list of length 300.
        //if (j==1) List sum_of_trees_tau(1);			// create a list of length 300.
        //if (j==1) List sum_of_trees_mat_tau(1);		// create a list of length 300.
        //if (j>1) List sum_of_trees_tau = prev_sum_trees_tau[curr_round_parent[k]];
        //if (j>1) List sum_of_trees_mat_tau = prev_sum_trees_mat_tau[curr_round_parent[k]];
        // Rcout << "Get to line 5438 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        
        if(curr_round_parent[k]==-1){			// If the k+1^th element of curr_round_parent is -1, do nothing. (-1 is a terminal node?)
        }else{									// If the k+1^th element of curr_round_parent is not equal to -1.
          // NEED TO THINK MORE ABOUT j==1 CASE. Also need the prev_sum_trees to be correct
          if(j==1){							// If in the SECOND round of the outter for-loop??
            
            sum_of_trees_tau = resize_bcf(sum_of_trees_tau,1);			// create a list of length 300.
            sum_of_trees_mat_tau = resize_bcf(sum_of_trees_mat_tau,1);		// create a list of length 300.
            
            // Rcout << "LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
            // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
            
            // Rcout << "Get to 5449 in mu round in loop j = " << j << ".\n";
            List other_tree_mulist=prev_sum_trees_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            NumericMatrix other_tree_mu=other_tree_mulist[0];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            
            //NumericMatrix other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            // Rcout << "Get to 5454 in mu round in loop j = " << j << ".\n";
            
            NumericMatrix other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            
            //NumericVector other_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];	// create vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
            
            //MIGHT NEED TO CHANGE THIS LINE
            List other_resids_mulist=prev_sum_tree_resids_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            NumericVector other_resids_mu=other_resids_mulist[0];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            // Rcout << "Line 5465 LENGTH OF other_resids_mu = " << other_resids_mu.size() << ".\n";
            NumericVector other_resids_tau=prev_sum_tree_resids_tau[curr_round_parent[k]];	// create vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
            
            sum_of_trees_mu[count]= other_tree_mu;										// add other_tree to the list sum_of_trees
            sum_of_trees_tau[count]= other_tree_tau;										// add other_tree to the list sum_of_trees
            sum_of_tree_resids_mu[count]=other_resids_mu-curr_round_preds_mu(_,k);									// add other_resids to the list sum_of_tree_resids
            sum_of_tree_resids_tau[count]=other_resids_tau-curr_round_preds_mu(_,k);									// add other_resids to the list sum_of_tree_resids
            
            
            List other_tree_mu_mat_list=prev_sum_trees_mat_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            NumericMatrix other_mat_mu=other_tree_mu_mat_list[0];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            NumericMatrix other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            sum_of_trees_mat_mu[count]=other_mat_mu;										// create a add other_mat to the list sum_of_trees_mat
            sum_of_trees_mat_tau[count]=other_mat_tau;										// create a add other_mat to the list sum_of_trees_mat
            count++;																// increment the count variable.
            // Rcout << "LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
            // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
            
            if(count==(size_mat-1)){												// If list size is too small
              size_mat=size_mat*2;												// double the size
              sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);					// double the length of sum_of_trees
              //sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);					// double the length of sum_of_trees
              sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,size_mat);		// double the length of sum_of_tree_resids
              sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,size_mat);		// double the length of sum_of_tree_resids
              sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);			// double the length of sum_of_trees_mat
              //sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);			// double the length of sum_of_trees_mat
            }
            // Rcout << "on  line 5490 LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
            // Rcout << "on  line 5491 LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
            // Rcout << "loop number" << j << "\n,";
          }else{
            List other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            List other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            List other_tree_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
            List other_tree_resids_tau=prev_sum_tree_resids_tau[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
            List other_mat_mu=prev_sum_trees_mat_mu[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            //List other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            for(int f=0;f<other_tree_mu.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
              if(is<NumericMatrix>(other_tree_mu[f])){								// if f+1^th element of other_tree is a matrix, do nothing
              }else{																// if f+1^th element of other_tree is not a matrix
                throw std::range_error("Line 5505 tree is not a numeric matrix!");		// throw an error
              }
              NumericMatrix treetoadd_mu=other_tree_mu[f];								// create matrix treetoadd equal to f+1^th element of other_tree
              if(is<NumericVector>(other_tree_resids_mu[f])){						// if f+1^th element of other_tree_resids is a NumericVector, do nothing
              }else{																// if f+1^th element of other_tree_resids is not a NumericVector
              	throw std::range_error("other resids not a numeric matrix!");	// throw an error
              }
              
              NumericVector resids_prevroundtemp_mu=other_tree_resids_mu[f];
              NumericVector residstoadd_mu=resids_prevroundtemp_mu-curr_round_preds_mu(_,k);						// create vector residstoadd equal to f+1^th element of other_tree_resids
              if(is<NumericMatrix>(other_mat_mu[f])){								// if f+1^th element of other_mat is a NumericMatrix, do nothing
                
              }else{																// if f+1^th element of other_mat is not a NumericMatrix
                throw std::range_error(" other mat not a numeric matrix!");		// throw an error
              }
              NumericMatrix mattoadd_mu=other_mat_mu[f];								// create matrix mattoadd equal to f+1^th element of other_mat
              
              sum_of_trees_mu[count]=treetoadd_mu;										// add treetoadd to sum_of_trees
              sum_of_tree_resids_mu[count]=residstoadd_mu;								// add residstoadd to sum_of_tree_resids
              sum_of_trees_mat_mu[count]=mattoadd_mu;									// add mattoadd to sum_of_trees_mat
              count++;															// inremet the count variable.
              
              if(count==(size_mat-1)){											// If list size is too small
                size_mat=size_mat*2;											// double the size
                sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);				// double the length of sum_of_trees
                sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,size_mat);	// double the length of sum_of_tree_resids
                sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);		// double the length of sum_of_trees_mat
              }
            }
            
            if(other_tree_tau.size()>size_mat){											// If list size is too small
              sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,other_tree_tau.size());	// increase the length of sum_of_tree_resids
            }
            for(int f=0;f<other_tree_tau.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
              NumericVector resids_prevroundtemp_tau=other_tree_resids_tau[f];
              NumericVector residstoadd_tau=resids_prevroundtemp_tau-curr_round_preds_mu(_,k);						// create vector residstoadd equal to f+1^th element of other_tree_resids

              sum_of_tree_resids_tau[f]=residstoadd_tau;								// add residstoadd to sum_of_tree_resids
            }
            
            
            
            
            List temp1_prev_sum_tree = prev_sum_trees_tau[curr_round_parent[k]];
            List temp1_prev_sum_tree_mat = prev_sum_trees_tau[curr_round_parent[k]];
            
            sum_of_trees_tau = resize_bcf(sum_of_trees_tau,temp1_prev_sum_tree.size());			// create a list of length 300.
            sum_of_trees_mat_tau = resize_bcf(sum_of_trees_mat_tau,temp1_prev_sum_tree_mat.size());		// create a list of length 300.
            sum_of_trees_tau = prev_sum_trees_tau[curr_round_parent[k]];
            sum_of_trees_mat_tau = prev_sum_trees_mat_tau[curr_round_parent[k]];
            // Rcout << "on  line 5553 LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
            
          }
          // Rcout << "on  line 5556 LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
          // Rcout << "on  line 5557 LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
          // Rcout << "loop number" << j << "\n,";
          
          sum_of_trees_mu[count]=temp_sum_trees_mu[k];										// add k+1^th element of temp_sum_trees to sum_of_trees
          sum_of_tree_resids_mu[count]=temp_sum_tree_resids_mu[k];							// add k+1^th element of temp_sum_tree_resids to sum_of_tree_resids
          sum_of_trees_mat_mu[count]=temp_sum_trees_mat_mu[k];								// add k+1^th element of temp_sum_trees_mat to sum_of_trees_mat
          count++;																	// increment the count variable
          
          // Rcout << "LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
          // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
          
          
          if(count==(size_mat-1)){													// If list size is too small
            size_mat=size_mat*2;													// double the size
            sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);						// double the length of sum_of_trees
            if(j==0) sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);						// double the length of sum_of_trees
            
            sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);				// double the length of sum_of_trees_mat_mu
            if(j==0) sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);				// double the length of sum_of_trees_mat_tau
            sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,size_mat);			// double the length of sum_of_trees_mat
          }
        }
        // Rcout << "Get to line 5581 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        sum_of_trees_mu=resize_bcf(sum_of_trees_mu,count);										// remove spaces that are not filled in.
        if(j==0) sum_of_trees_tau=resize_bcf(sum_of_trees_tau,count);
        									// remove spaces that are not filled in.
        sum_of_trees_mat_mu=resize_bcf(sum_of_trees_mat_mu,count);								// remove spaces that are not filled in.
        if(j==0) sum_of_trees_mat_tau=resize_bcf(sum_of_trees_mat_tau,count);								// remove spaces that are not filled in.
        sum_of_tree_resids_mu=resize_bcf(sum_of_tree_resids_mu,count);							// remove spaces that are not filled in.
        // Rcout << "Get to line 5587 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        if(j>0){
          List other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          sum_of_tree_resids_tau=resize_bcf(sum_of_tree_resids_tau,other_tree_tau.size());							// remove spaces that are not filled in.
        }
        // Rcout << "Get to line 5591 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        if(curr_round_parent[k]!=-1){															// If the k+1^th element of curr_round_parent is -1, (-1 is a terminal node?)
          overall_sum_trees_mu[overall_count]=sum_of_trees_mu;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_tau[overall_count]=sum_of_trees_tau;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
          // Rcout <<"ADD ELEMENT TO LIST OF SUM TREE LISTS OF LENGTH"<< sum_of_trees_tau.size() << ".\n";
          overall_sum_tree_resids_mu[overall_count]=sum_of_tree_resids_mu;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
          overall_sum_tree_resids_tau[overall_count]=sum_of_tree_resids_tau;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_mat_mu[overall_count]=sum_of_trees_mat_mu;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_mat_tau[overall_count]=sum_of_trees_mat_tau;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
          overall_sum_BIC=temp_BIC;															// Let overall_sum_BIC equal temp_BIC, the vector of BICs.
          overall_sum_preds_outcome=Rcpp::as<arma::mat>(temp_preds_outcome);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          overall_sum_preds_mu=Rcpp::as<arma::mat>(temp_preds_mu);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          overall_sum_preds_tau=Rcpp::as<arma::mat>(temp_preds_tau);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          if(is_test_data==1){	// If there is test data, overall_sum_test_preds equal temp_test_preds, the matrix of out-of-sample predictions (columns correpond to different models?)
            //Rcout << "Get to Line 5609 in loop j = " << j << ". k= "<< k << " .\n";
            overall_sum_test_preds_outcome=Rcpp::as<arma::mat>(temp_test_preds_outcome);
            overall_sum_test_preds_mu=Rcpp::as<arma::mat>(temp_test_preds_mu);
            overall_sum_test_preds_tau=Rcpp::as<arma::mat>(temp_test_preds_tau);
          }					
          overall_count++;																	// increment overall_count
          if(overall_count==(overall_size-1)){												// If overall_size is too small
            overall_size=overall_size*2;													// double the size
            overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
            overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
            overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
            overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
            overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
            overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
          }
        }  
      }
      // Rcout << "Get to line 5619 in loop j = " << j << ".\n";
      if(only_max_num_trees==0){
        
      //check if there were any trees from the previous round that didn't have daughter trees grown.
      //create vector to count number of possible parents for previous round
      if(j>0){																					// If not the first round of the outter loop
        IntegerVector prev_par_no_child=match(prev_par,curr_round_parent);						// create a vector equal to indices of the positions of the (the first matches of the) elements of prev_par in curr_round_parent.
        if(any(is_na(prev_par_no_child))){														// any of the vector of matches are NA (no match?)
          IntegerVector t4=ifelse(is_na(prev_par_no_child),1,0);								// create a vector t4 equal to 1 for the NA values, 0 otherwise.
          for(int h=0;h<prev_par_no_child.size();h++){										// for-loop of length equal to that of prev_par_no_child
            if(t4[h]==1){																	// If h+1^th element of vector of matches is NA
              if(prev_round_BIC2[h]-lowest_BIC<=log(c)){									// If the h+1^th model (from the previous round?) is in Occam's window
                SEXP s_mu = prev_sum_trees_mu[h];												// create a pointer to S expression type equal to the h+1^th element of prev_sum_trees (a tree table or list of tree tables from the previous round?)
                SEXP s_tau = prev_sum_trees_tau[h];												// create a pointer to S expression type equal to the h+1^th element of prev_sum_trees (a tree table or list of tree tables from the previous round?)
                //SEXP s_resid_mu = prev_sum_tree_resids_mu[h]; 
                SEXP s_resid_tau = prev_sum_tree_resids_tau[h];
                
                // if(is<List>(s_resid_mu)){	
                //   if(is<List>(s_resid_tau)){
                //     Rcout << "mu resid list and tau resid list. mu round j= " << j << ".\n"; 
                //   }else{
                //     Rcout << "mu resid list and tau resid vector?. mu round j= " << j << ".\n"; 
                //     
                //   }
                // }else{
                //   if(is<List>(s_resid_tau)){
                //     Rcout << "mu resid vector and tau resid list. mu round j= " << j << ".\n"; 
                //     
                //   }else{
                //     Rcout << "mu resid vector and tau resid vector. mu round j= " << j << ".\n"; 
                //     
                //   }
                // }
                
                if(is<List>(s_mu)){														// If prev_sum_trees[h] is a list

                  if(is<List>(s_tau)){
                    List tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    List tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                    List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                    NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }
                    List treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    List treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 5673 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }else{
                    List tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    NumericMatrix tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                      NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }                    
                    List treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    NumericMatrix treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 5708 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }


                }else{																// If prev_sum_trees[h] is NOT a list
                  if(is<List>(s_tau)){


                    NumericMatrix tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    List tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                      NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }                    
                    NumericMatrix treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    List treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 5749 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }else{
                    NumericMatrix tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    NumericMatrix tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                      NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }                    
                    NumericMatrix treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    NumericMatrix treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 5784 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }

                }
              }
            }
          }
        }
      }
      }
      prev_round_preds_outcome=temp_preds_outcome;															// let prev_round_preds equal to the matrix of predictions.
      prev_round_preds_mu=temp_preds_mu;															// let prev_round_preds equal to the matrix of predictions.
      prev_round_preds_tau=temp_preds_tau;															// let prev_round_preds equal to the matrix of predictions.
      if(is_test_data==1){
        //Rcout << "Get to Line 5802 in loop j = " << j  << ".\n";
        prev_round_test_preds_outcome=temp_test_preds_outcome;								// if there is test data, let prev_round_test_preds equal the test data predictions
        prev_round_test_preds_mu=temp_test_preds_mu;								// if there is test data, let prev_round_test_preds equal the test data predictions
        prev_round_test_preds_tau=temp_test_preds_tau;								// if there is test data, let prev_round_test_preds equal the test data predictions
      }
      prev_round_BIC=temp_BIC;																// let prev_round_BIC equal the vector of BICs
      prev_round_BIC2=temp_BIC;																// let prev_round_BIC2 equal the vector of BICs
      prev_round_preds2_outcome=Rcpp::as<arma::mat>(temp_preds_outcome);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      prev_round_preds2_mu=Rcpp::as<arma::mat>(temp_preds_mu);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      prev_round_preds2_tau=Rcpp::as<arma::mat>(temp_preds_tau);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      if(is_test_data==1){
        //Rcout << "Get to Line 5813 in loop j = " << j  << ".\n";
        prev_round_test_preds2_outcome=Rcpp::as<arma::mat>(temp_test_preds_outcome);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
        prev_round_test_preds2_mu=Rcpp::as<arma::mat>(temp_test_preds_mu);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
        prev_round_test_preds2_tau=Rcpp::as<arma::mat>(temp_test_preds_tau);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
      }
      resids=temp_resids;																		// let resids equal the matrix of residuals
      parent=temp_parent;																		// let parent equal the parent vector
      overall_sum_trees_mu=resize_bcf(overall_sum_trees_mu,overall_count);								// remove spaces that are not filled in.
      overall_sum_trees_tau=resize_bcf(overall_sum_trees_tau,overall_count);								// remove spaces that are not filled in.
      overall_sum_tree_resids_mu=resize_bcf(overall_sum_tree_resids_mu,overall_count);					// remove spaces that are not filled in.
      overall_sum_tree_resids_tau=resize_bcf(overall_sum_tree_resids_tau,overall_count);					// remove spaces that are not filled in.
      overall_sum_trees_mat_mu=resize_bcf(overall_sum_trees_mat_mu,overall_count);						// remove spaces that are not filled in.
      overall_sum_trees_mat_tau=resize_bcf(overall_sum_trees_mat_tau,overall_count);						// remove spaces that are not filled in.
      
      
      if(j==0){																		// if in the first round of the outer for-loop (j==0)
        prev_sum_trees_mu=temp_sum_trees_mu;														// let prev_sum_trees equal the list of tree tables (from the current round)
        //ADDING TO MU TREES, therefore nothing yet added to tau trees
        List prev_sum_trees_tau(prev_sum_trees_mu.size());
        prev_sum_trees_mat_mu=temp_sum_trees_mat_mu;												// let prev_sum_trees_mat equal the list of tree matrice
        List prev_sum_trees_mat_tau(prev_sum_trees_mu.size());
        
        for(int p=0;p<curr_round_lik.size();p++){
          prev_sum_trees_tau[p] = start_tree_bcf(0,0); // maybe should be empty list rather than list of empty trees?
          prev_sum_trees_mat_tau[p] = start_matrix_bcf(n); // maybe should be empty list? .. model shouldn't even have stub tree
        }
        prev_sum_tree_resids_mu=temp_sum_tree_resids_mu;											// let prev_sum_tree_resids equal the list of residual vectors (from the current round)
        //prev_sum_tree_resids_tau=;											// let prev_sum_tree_resids equal the list of residual vectors (from the current round)
        
        //NumericMatrix test=prev_sum_trees[0];												// create a matrix test equal to the first element of temp_sum_trees (first obtained tree table)
        overall_sum_trees_mu=resize_bcf(temp_sum_trees_mu,temp_sum_trees_mu.size());						// remove spaces that are not filled in.
        overall_sum_trees_tau=resize_bcf(prev_sum_trees_tau,temp_sum_trees_mu.size());// INTENTIONALLY USING SIZE OF MU						// remove spaces that are not filled in.
        overall_sum_trees_mu=temp_sum_trees_mu;													// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
        overall_sum_trees_tau=prev_sum_trees_tau; // not sure about this													// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
        overall_sum_tree_resids_mu=resize_bcf(temp_sum_tree_resids_mu,temp_sum_tree_resids_mu.size());	// remove spaces that are not filled in.
        overall_sum_tree_resids_mu=temp_sum_tree_resids_mu;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)
        overall_sum_tree_resids_tau=resize_bcf(prev_sum_tree_resids_tau,prev_sum_tree_resids_tau.size());	// remove spaces that are not filled in.
        overall_sum_tree_resids_tau=prev_sum_tree_resids_tau;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)
        overall_sum_trees_mat_mu=temp_sum_trees_mat_mu;											// let overall_sum_trees_mat equal the RESIZED list of tree matrice
        overall_sum_trees_mat_tau=prev_sum_trees_mat_tau;											// let overall_sum_trees_mat equal the RESIZED list of tree matrice
        overall_sum_trees_mat_mu=resize_bcf(temp_sum_trees_mat_mu,temp_sum_trees_mat_mu.size());				// remove spaces that are not filled in.
        overall_sum_trees_mat_tau=resize_bcf(prev_sum_trees_mat_tau,temp_sum_trees_mat_mu.size());				// remove spaces that are not filled in.
        overall_sum_BIC=temp_BIC;															// let overall_sum_BIC equal the vector of BICs
        overall_sum_preds_outcome= Rcpp::as<arma::mat>(temp_preds_outcome);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
        overall_sum_preds_mu= Rcpp::as<arma::mat>(temp_preds_mu);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
        overall_sum_preds_tau= Rcpp::as<arma::mat>(temp_preds_tau);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
        if(is_test_data==1){ // if there is test data, let overall_sum_test_preds be an arma mat copy of the test data predictions
          //Rcout << "Get to Line 5860 in loop j = " << j  << ".\n";
          overall_sum_test_preds_outcome= Rcpp::as<arma::mat>(temp_test_preds_outcome);
          overall_sum_test_preds_mu= Rcpp::as<arma::mat>(temp_test_preds_mu);
          overall_sum_test_preds_tau= Rcpp::as<arma::mat>(temp_test_preds_tau);
        }
      }else{																					// if not in the first round of the outer for-loop. (i.e. j>0)
        
        prev_sum_trees_mu=overall_sum_trees_mu;													// let prev_sum_trees equal the list of tree tables (up to the current round??)
        prev_sum_trees_tau=overall_sum_trees_tau;													// let prev_sum_trees equal the list of tree tables (up to the current round??)
        prev_sum_tree_resids_mu=overall_sum_tree_resids_mu;										// let prev_sum_tree_resids equal the list of residual vectors (up tp the current round)
        prev_sum_tree_resids_tau=overall_sum_tree_resids_tau;										// let prev_sum_tree_resids equal the list of residual vectors (up tp the current round)
        prev_sum_trees_mat_mu=overall_sum_trees_mat_mu;											// let prev_sum_trees_mat equal the list of tree matrice (up to the current round??)
        prev_sum_trees_mat_tau=overall_sum_trees_mat_tau;											// let prev_sum_trees_mat equal the list of tree matrice (up to the current round??)
        
        
        prev_round_BIC2=overall_sum_BIC;													// let prev_round_BIC2 equal the vector of BICs (up to the current round??)
        prev_round_preds2_outcome=overall_sum_preds_outcome;												// let prev_round_preds2 equal the predictions matrix
        prev_round_preds2_mu=overall_sum_preds_mu;												// let prev_round_preds2 equal the predictions matrix
        prev_round_preds2_tau=overall_sum_preds_tau;												// let prev_round_preds2 equal the predictions matrix
        if(is_test_data==1){
          //Rcout << "Get to Line 5880 in loop j = " << j  << ".\n";
          prev_round_test_preds2_outcome=overall_sum_test_preds_outcome;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
          prev_round_test_preds2_mu=overall_sum_test_preds_mu;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
          prev_round_test_preds2_tau=overall_sum_test_preds_tau;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
        }
      }
      overall_overall_sum_trees_mu[oo_count]=overall_sum_trees_mu;									// Add the current outer loop's table list to overall_overall_sum_trees (list of lists?? OR list of lists of lists??)
      overall_overall_sum_trees_tau[oo_count]=overall_sum_trees_tau;									// Add the current outer loop's table list to overall_overall_sum_trees (list of lists?? OR list of lists of lists??)
      overall_overall_sum_tree_resids_mu[oo_count]=overall_sum_tree_resids_mu;						// Add the current outer loop's list of residual vectors to the list overall_overall_sum_tree_resids (list of lists)
      overall_overall_sum_tree_resids_tau[oo_count]=overall_sum_tree_resids_tau;						// Add the current outer loop's list of residual vectors to the list overall_overall_sum_tree_resids (list of lists)
      overall_overall_sum_trees_mat_mu[oo_count]=overall_sum_trees_mat_mu;							// Add the current outer loop's list of matrices to the list overall_overall_sum_trees_mat (list of lists)
      overall_overall_sum_trees_mat_tau[oo_count]=overall_sum_trees_mat_tau;							// Add the current outer loop's list of matrices to the list overall_overall_sum_trees_mat (list of lists)
      overall_overall_sum_BIC[oo_count]=overall_sum_BIC;										// Add the vector of BICs to the list overall_overall_sum_BIC. (list of vectors?)
      oo_count ++;																					// increment the count 
      if(oo_count==(oo_size-1)){																		// If lists are not large enough
        oo_size=oo_size*2;																			// double the size.
        overall_overall_sum_trees_mu=resize_bigger_bcf(overall_overall_sum_trees_mu,oo_size);					// double the length of overall_overall_sum_trees
        overall_overall_sum_trees_tau=resize_bigger_bcf(overall_overall_sum_trees_tau,oo_size);					// double the length of overall_overall_sum_trees
        overall_overall_sum_tree_resids_mu=resize_bigger_bcf(overall_overall_sum_tree_resids_mu,oo_size);		// double the length of overall_overall_sum_tree_resids
        overall_overall_sum_tree_resids_tau=resize_bigger_bcf(overall_overall_sum_tree_resids_tau,oo_size);		// double the length of overall_overall_sum_tree_resids
        overall_overall_sum_trees_mat_mu=resize_bigger_bcf(overall_overall_sum_trees_mat_mu,oo_size);			// double the length of overall_overall_sum_trees_mat
        overall_overall_sum_trees_mat_tau=resize_bigger_bcf(overall_overall_sum_trees_mat_tau,oo_size);			// double the length of overall_overall_sum_trees_mat
        overall_overall_sum_BIC=resize_bigger_bcf(overall_overall_sum_BIC,oo_size);						// double the length of overall_overall_sum_BIC
      }    
      overall_overall_sum_preds_outcome=overall_sum_preds_outcome;											// let overall_overall_sum_preds equal the prediction matrix.
      overall_overall_sum_preds_mu=overall_sum_preds_mu;											// let overall_overall_sum_preds equal the prediction matrix.
      overall_overall_sum_preds_tau=overall_sum_preds_tau;											// let overall_overall_sum_preds equal the prediction matrix.
      
      if(is_test_data==1){
        //Rcout << "Get to Line 5909 in loop j = " << j  << ".\n";
        overall_overall_sum_test_preds_outcome=overall_sum_test_preds_outcome;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
        overall_overall_sum_test_preds_mu=overall_sum_test_preds_mu;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
        overall_overall_sum_test_preds_tau=overall_sum_test_preds_tau;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
      }
      //overall_trees_mu[j]=curr_round_trees_mu;														// let the j+1^th element of the list overall_trees be the list of trees obtained in the current round of the outer loop (list of lists of tree table matrices)
      //overall_trees_tau[j]=curr_round_trees_tau;														// let the j+1^th element of the list overall_trees be the list of trees obtained in the current round of the outer loop (list of lists of tree table matrices)
      //overall_mat_mu.push_back(curr_round_mat_mu);													// append the list of current round of tree matrices to the list overall_mat (list of lists of matrices)
      //overall_mat_tau.push_back(curr_round_mat_tau);													// append the list of current round of tree matrices to the list overall_mat (list of lists of matrices)
      overall_lik.push_back(curr_round_lik);													// append the current round's vector of BICs to the list overall_lik (list of vectors of BICs).
      prev_par=seq_len(overall_sum_trees_mu.size())-1;											// let prev_par equal a sequence 0,1,2,3,..., up to length of overall_sum_trees minus one
      
      
      
    }
    //	END OF mu(x) TREE CODE
    
    // START OF tau(x) TREE CODE	
    
    if(j<ntree_moderate){
      if(j>0){
        if(only_max_num_trees==1){
          lowest_BIC=100000;
        }
        }
      int overall_size=300;						// create a variable overall_size. Initialized equal to 0.
      List overall_sum_tree_resids_mu(overall_size);	// create a list of length 300
      List overall_sum_tree_resids_tau(overall_size);	// create a list of length 300
      
      List overall_sum_trees_mu(overall_size);		// create a list of length 300
      List overall_sum_trees_mat_mu(overall_size);	// create a list of length 300
      
      List overall_sum_trees_tau(overall_size);		// create a list of length 300
      List overall_sum_trees_mat_tau(overall_size);	// create a list of length 300
      
      
      int overall_count=0;						// set overall_count equal to 0.
      //		parent indexes the whole models to which the a tree can be appended
      IntegerVector parent;						// create vector.
      NumericVector curr_round_lik;				// create vector	// To be filled with BICs for whole models suggested after a tree appended
      List curr_round_trees_mu;						// create list.
      List curr_round_trees_tau;						// create list.
      
      List curr_round_mat_mu;						// create list.
      List curr_round_mat_tau;						// create list.
      
      NumericVector curr_BIC;						// create vector.
      //		next line probably shouldn't need to be duplicated for curr_round_parent_mu and curr_round_parent_tau
      IntegerVector curr_round_parent;			// create vector.
      NumericVector overall_sum_BIC;				// create vector.
      
      arma::mat overall_sum_preds_outcome;				// create matrix.
      arma::mat overall_sum_preds_mu;				// create matrix.
      arma::mat overall_sum_preds_tau;				// create matrix.
      
      arma::mat overall_sum_test_preds_outcome;			// create matrix.
      arma::mat overall_sum_test_preds_mu;			// create matrix.
      arma::mat overall_sum_test_preds_tau;			// create matrix.
      
      if(j==0){									// If in the first round of the for-loop.
        parent.push_back(0);					// append a 0 to the end of the parent vector. (first and only element of parent vector so far).
        //first_round=1;							// set the variable first_round equal to 1.
      }else{										// If not in the first round of the for-loop.
        //first_round=0;							// set the variable first_round equal to 0.
      }
      //		The _tau in resids_cp_mat_tau is probably unnecessary, but including it to remove ambiguity.
      //		Replace with just resids_cp_mat for memory efficiency after code is all working.
      //		similarly, cp_mat_list_tau is used (perhaps unnecessarily) here instead of cp_mat_list
      List resids_cp_mat_tau(resids.ncol());												// create a list of length equal to the number of colmns of resids
      int resids_count=0;																// create a variable resids_count. Initialize equal to zero.
      std::vector<int> err_list(resids.ncol());										// create a vector err_list of length equal to the number of colmns of resids
      //get best splits
      for(int f=0;f<resids.ncol();f++){												// for-loop of length equal to the unmber of columns of resids
        if(gridpoint==0){															// If input gridpoint equals 0. i.e. the PELT method will be used.
          cp_mat_list_tau=make_pelt_cpmat_tau_bcf(x_moderate,resids(_,f),gridsize_tau,num_cp_tau,z);				// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																		// If input gridpoint equals 1. i.e. the PELT method will be used.
          cp_mat_list_tau=make_gridpoint_cpmat_tau_bcf(x_moderate,resids(_,f),gridsize_tau,num_cp_tau,z);			// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        resids_cp_mat_tau[resids_count]=cp_mat_list_tau[0];									// let the resids_count+1^th element be a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. 
        err_list[resids_count]=cp_mat_list_tau[1];										// let the resids_count+1^th element be a list that records error(s?)
        resids_count++;																// increment resids_count by 1.
      }
      resids_cp_mat_tau=resize_bcf(resids_cp_mat_tau,resids_count);								// remove elements of resids_cp_mat_tau that are not fileld in. How is  it possible that elements are not filled in??
      err_list.resize(resids_count);													// remove elements of err_list that are not fileld in. How is  it possible that elements are not filled in??
      parent=seq_len(tree_table_tau.size())-1;											// set parent equal to a vector 0,1,2,3,..., up to the length of tree_table minus one.
      if(is_true(all(as<IntegerVector>(wrap(err_list))==1))){							// if all elements of err_list equal 1.
        if(j==0){																	// If in the first round of the for-loop.
          throw std::range_error("No split points could be found to grow trees");
        }else{																		// If not in the first round of the for-loop.
          throw std::range_error("No Tau trees can be grown for the number of iterations desired, as no splits were found.Please try fewer iterations.");
        }
      }
      // Rcout << "Get to 5984 in tau round in loop j = " << j << ".\n";
      
      //get current set of trees.
      if(j==0){						// If in the first round of the for-loop.
        //CART_BMA_tau=get_best_trees_tau(x_control_a, x_moderate_a,z,resids,
        //a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
        //tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
        //lowest_BIC,first_round,parent,resids_cp_mat_tau,as<IntegerVector>(wrap(err_list)),
        //x_control_test,x_moderate_test,test_z,
        //alpha_mu,alpha_tau,beta_mu,beta_tau,
        //is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
        //split_rule_node,gridpoint,maxOWsize);
        CART_BMA_tau=get_best_trees_sum_tau_round1_bcf(x_control_a, x_moderate_a,z,resids,
                                                   a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                                   tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                                   lowest_BIC,//first_round,
                                                   parent,resids_cp_mat_tau,as<IntegerVector>(wrap(err_list)),
                                                   x_control_test,x_moderate_test,test_z,
                                                   alpha_mu,alpha_tau,beta_mu,beta_tau,
                                                   is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                                   split_rule_node,gridpoint,maxOWsize,
                                                   prev_sum_trees_mu,
                                                   //prev_sum_trees_tau,
                                                   prev_sum_trees_mat_mu,
                                                   //prev_sum_trees_mat_tau,
                                                   y_scaled,num_splits_mu,num_splits_tau,gridsize_tau,zero_split);	// function defined on line 1953.
      }else{							// If not in the first round of the for-loop.
        //if j >0 then sum of trees become a list so need to read in list and get likelihood for each split point and terminal node
        CART_BMA_tau=get_best_trees_sum_tau_bcf(x_control_a, x_moderate_a,z,resids,
                                            a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                            tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                            lowest_BIC,//first_round,
                                            parent,resids_cp_mat_tau,as<IntegerVector>(wrap(err_list)),
                                            x_control_test,x_moderate_test,test_z,
                                            alpha_mu,alpha_tau,beta_mu,beta_tau,
                                            is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                            split_rule_node,gridpoint,maxOWsize,
                                            prev_sum_trees_mu,prev_sum_trees_tau,
                                            prev_sum_trees_mat_mu,prev_sum_trees_mat_tau,y_scaled,
                                            num_splits_mu,num_splits_tau,gridsize_tau,zero_split);	// function defined on line 1953.
      }
      // Rcout << "Get to 6023 in tau round in loop j = " << j << ".\n";
      
      curr_round_lik=CART_BMA_tau[0];							// vector of BICs (for whole sum-of tree-models after suggested trees added). Should be ordered ascending
      //curr_round_trees_mu=CART_BMA_tau[1];						// list of tree tables
      curr_round_trees_tau=CART_BMA_tau[1];					
      //curr_round_mat_mu=CART_BMA_tau[3];							// list of tree matrices
      curr_round_mat_tau=CART_BMA_tau[2];							// list of tree matrices
      curr_round_parent=CART_BMA_tau[3];						// vector of tree parent numbers
      NumericMatrix curr_round_preds_tau=CART_BMA_tau[4];			// (in-sample single tree predictions for trees to add) matrix rows correspond to different units/individuals, columns corresponds to predictions from different (single or sums-of?) trees.
      curr_BIC=CART_BMA_tau[5];								// lowest BIC among trees?
      NumericMatrix curr_round_test_preds_tau=CART_BMA_tau[6];	// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspond to predictions from different (single or sums-of?) trees.
      if(curr_round_lik.size()==0) {						// If number of sum of tree models is zero?
        // Rcout << "curr_round_lik.size()==0 BREAK in tau round in loop j = " << j << ".\n";
        //REMOVE THIS ERROR IF WANT TO ALLOW LESS THAN MAX NUMBER OF TREES
        //throw std::range_error("No tau trees chosen in round");
        break;											// break out of for-loop
      } 
      
      if(curr_BIC[0]<lowest_BIC){							// If the lowest BIC obtained by get_best_trees_sum is less than the currently saved lowest value
        lowest_BIC=curr_BIC[0];							// reset lowest_BIC to the new lowest value
      }
      tree_table_mu=List();									// reset tree_table to an empty list
      tree_mat_mu=List();									// reset tree_mat to an empty list.
      tree_table_tau=List();									// reset tree_table to an empty list
      tree_mat_tau=List();									// reset tree_mat to an empty list.
      int lsize=curr_round_lik.size();					// create a variable equal to the number of sum of tree models returned by get_best_trees_sum
      tree_table_mu=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_mat_mu=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_table_tau=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_mat_tau=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_outcome(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_mu(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_tau(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      
      NumericMatrix temp_test_preds_outcome(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_test_preds_mu(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_test_preds_tau(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      
      NumericMatrix temp_resids(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericVector temp_parent(curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericVector temp_BIC(curr_round_lik.size());		// create a vector of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      List temp_sum_trees_mu(lsize);							// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_trees_tau(lsize);							// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      List temp_sum_tree_resids(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_tree_resids_mu(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_tree_resids_tau(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      List temp_sum_trees_mat_mu(lsize);						// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_trees_mat_tau(lsize);						// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum

      int count=0; 										// create a count variable. Initialized equal to zero.
      for(int k=0;k<curr_round_lik.size();k++){			// create a for-loop of length equal to the number of sum of tree models returned by get_best_trees_sum
        tree_table_mu[count]=start_tree_bcf(mu_mu,sigma_mu_mu);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
        tree_mat_mu[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
        tree_table_tau[count]=start_tree_bcf(mu_tau,sigma_mu_tau);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
        tree_mat_tau[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
        if(j==0){										// if in first round of for-loop
          temp_preds_outcome(_,k)=prev_round_preds_mu(_,curr_round_parent[k])+z*curr_round_preds_tau(_,k);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
          temp_preds_mu(_,k)=prev_round_preds_mu(_,curr_round_parent[k]);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
          temp_preds_tau(_,k)=curr_round_preds_tau(_,k);
          if(is_test_data==1){
            //Rcout << "Get to Line 6103 in loop j = " << j  << ".\n";
            temp_test_preds_outcome(_,k)=prev_round_test_preds_mu(_,curr_round_parent[k])+z*curr_round_test_preds_tau(_,k);
            temp_test_preds_mu(_,k)=prev_round_test_preds_mu(_,curr_round_parent[k]);
            temp_test_preds_tau(_,k)=curr_round_test_preds_tau(_,k);
          }
          // If there is test data, let the k+1^th column of temp_test_preds be the k+1^th column of curr_round_test_preds. These are the out-of-sample predictions of from the k+1^th model.
          temp_resids(_,k)=y_scaled-temp_preds_outcome(_,k);	// Let the k+1^th column of temp_resids be the outcome minus the predictons from the k+1^th model 
          temp_parent[k]=-1;							// Let the K=1^th element of temp_parent be -1.
          temp_BIC[k]=curr_round_lik[k];				// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
          //temp_sum_trees_mu[count]=curr_round_trees_mu[k];	// Add the tree table to the list temp_sum_trees
          temp_sum_trees_tau[count]=curr_round_trees_tau[k];	// Add the tree table to the list temp_sum_trees
          
          //temp_sum_trees_mat_mu[count]=curr_round_mat_mu[k];	// Add the tree matrix to temp_sum_trees_mat
          temp_sum_trees_mat_tau[count]=curr_round_mat_tau[k];	// Add the tree matrix to temp_sum_trees_mat
          //NOT SURE ABOUT THE FOLLOWING LINE
          //temp_sum_tree_resids[count]=resids(_,0);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
          //temp_sum_tree_resids_mu[count]=y_scaled-z*curr_round_preds_tau(_,k);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
          temp_sum_tree_resids_tau[count]=resids(_,curr_round_parent[k]);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
          
          
        }else{											// If not in the first round of the for-loop.
          NumericVector curr_temp_pred_outcome=prev_round_preds_mu(_,curr_round_parent[k]) + 
            z*prev_round_preds_tau(_,curr_round_parent[k]) + 
            z*curr_round_preds_tau(_,k);	// curr_temp_pred is the sum of the current round predictions and the previous round predictions?? Each round is for one tree? and explain more of the residuals in each round?
          NumericVector curr_temp_pred_mu=prev_round_preds_mu(_,curr_round_parent[k]);					
          NumericVector curr_temp_pred_tau=prev_round_preds_tau(_,curr_round_parent[k])+curr_round_preds_tau(_,k);	
          
          NumericVector curr_temp_test_pred_outcome;			// create a vector
          NumericVector curr_temp_test_pred_mu;			// create a vector
          NumericVector curr_temp_test_pred_tau;			// create a vector
          
          if(is_test_data==1) {						// If there is test data.
            //Rcout << "Get to Line 6135 in loop j = " << j  << ".\n";
            curr_temp_test_pred_outcome=prev_round_test_preds_mu(_,curr_round_parent[k])+
              test_z*prev_round_test_preds_tau(_,curr_round_parent[k])+
              test_z*curr_round_test_preds_tau(_,k);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            curr_temp_test_pred_mu=prev_round_test_preds_mu(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            curr_temp_test_pred_tau=prev_round_test_preds_tau(_,curr_round_parent[k])+
              curr_round_test_preds_tau(_,k);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            
            temp_test_preds_outcome(_,k) = curr_temp_test_pred_outcome;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            temp_test_preds_mu(_,k) = curr_temp_test_pred_mu;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            temp_test_preds_tau(_,k) = curr_temp_test_pred_tau;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            
          }
          temp_BIC[k]=curr_round_lik[k];									// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
          // Rcout << "TEST LINE 6052 temp_BIC[k]=  "<< temp_BIC[k] << " .\n";
          
          
          temp_preds_outcome(_,k)=curr_temp_pred_outcome;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          temp_preds_mu(_,k)=curr_temp_pred_mu;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          temp_preds_tau(_,k)=curr_temp_pred_tau;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          
          temp_resids(_,k)=y_scaled-curr_temp_pred_outcome;						// Let the k+1^th column of temp_resids be the new residuals after appending the k+1^th
          temp_parent[k] = k;												// Let k+1^th element of temp_parent equal k.
          //temp_sum_trees_mu[count]=curr_round_trees_mu[k];						// Add the tree table to the list temp_sum_trees.
          temp_sum_trees_tau[count]=curr_round_trees_tau[k];						// Add the tree table to the list temp_sum_trees.
          
          //temp_sum_trees_mat[count]=curr_round_mat_mu[k];					// Add the tree matrix to temp_sum_trees_mat.
          temp_sum_trees_mat_tau[count]=curr_round_mat_tau[k];					// Add the tree matrix to temp_sum_trees_mat.
          temp_sum_tree_resids_tau[count]=resids(_,curr_round_parent[k]);		// Add the curr_round_parent[k]+1^th column of resids to temp_sum_tree_resids. resids contains predictions from previous rounds?
        }
        count++;													// Increment the count by 1.(Note this is within the innermost for-loop).
      }  
      if(curr_round_lik.size()==0){									// if no new trees outputted in current round (by get_best_trees_sum? Why not throw this error earlier, at line 2350)
        throw std::range_error("No trees chosen in last round");
      }
      // Rcout << "GET TO LINE 6151 .\n";
      
      for(int k=0;k<curr_round_lik.size();k++){	// create a for-loop of length equal to the number of sum of tree models returned by get_best_trees_sum
        int size_mat=300;						// create a variable, Initializd equal to 300.
        List sum_of_trees_mu(size_mat);			// create a list of length 300.
        List sum_of_trees_tau(size_mat);			// create a list of length 300.
        List sum_of_tree_resids_mu(size_mat);		// create a list of length 300.
        List sum_of_tree_resids_tau(size_mat);		// create a list of length 300.
        List sum_of_trees_mat_mu(size_mat);		// create a list of length 300.
        List sum_of_trees_mat_tau(size_mat);		// create a list of length 300.
        int count=0;							// create a variable count equal to 0. (Count was already define, so could remove "int" at start of this line and just reset cound to 0).
        
        if(curr_round_parent[k]==-1){			// If the k+1^th element of curr_round_parent is -1, do nothing. (-1 is a terminal node?)
        }else{									// If the k+1^th element of curr_round_parent is not equal to -1.
          // NEED TO THINK MORE ABOUT j==1 CASE. Also need the prev_sum_trees to be correct
          
          if(j==0){
            sum_of_trees_mu = resize_bcf(sum_of_trees_mu,1);			// create a list of length 300.
            sum_of_trees_mat_mu = resize_bcf(sum_of_trees_mat_mu,1);		// create a list of length 300.
            
            // Rcout << "LENGTH OF LIST SUM_OF_TREES_MU = " << sum_of_trees_mu.size() << ".\n";
            // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_MU = " << sum_of_trees_mat_mu.size() << ".\n";
            
            NumericMatrix other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            //NumericMatrix other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            NumericVector other_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];	// create vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
            // Rcout << "LENGTH OF LIST other_resids_mu = " << other_resids_mu.size() << ".\n";
            sum_of_trees_mu[count]= other_tree_mu;										// add other_tree to the list sum_of_trees
            //sum_of_trees_tau[count]= other_tree_tau;										// add other_tree to the list sum_of_trees
            sum_of_tree_resids_mu[count]=other_resids_mu-z*curr_round_preds_tau(_,k);									// add other_resids to the list sum_of_tree_resids
            NumericMatrix other_mat_mu=prev_sum_trees_mat_mu[curr_round_parent[k]];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            //NumericMatrix other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            sum_of_trees_mat_mu[count]=other_mat_mu;										// create a add other_mat to the list sum_of_trees_mat
            //sum_of_trees_mat_tau[count]=other_mat_tau;										// create a add other_mat to the list sum_of_trees_mat
            //count++;																// increment the count variable.
            
            //if(count==(size_mat-1)){												// If list size is too small
            //  size_mat=size_mat*2;												// double the size
            //  sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);					// double the length of sum_of_trees
            //  sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);					// double the length of sum_of_trees
            //  sum_of_tree_resids=resize_bigger_bcf(sum_of_tree_resids,size_mat);		// double the length of sum_of_tree_resids
            //  sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);			// double the length of sum_of_trees_mat
            //  sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);			// double the length of sum_of_trees_mat
            //}
          
            }else{
              List other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
              List other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
              List other_tree_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
              List other_tree_resids_tau=prev_sum_tree_resids_tau[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
              //List other_mat_mu=prev_sum_trees_mat_mu[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
              List other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
              for(int f=0;f<other_tree_tau.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
                if(is<NumericMatrix>(other_tree_tau[f])){								// if f+1^th element of other_tree is a matrix, do nothing
                }else{																// if f+1^th element of other_tree is not a matrix
                  throw std::range_error("Line 6211 tree is not a numeric matrix!");		// throw an error
                }
                NumericMatrix treetoadd_tau=other_tree_tau[f];								// create matrix treetoadd equal to f+1^th element of other_tree
                if(is<NumericVector>(other_tree_resids_tau[f])){						// if f+1^th element of other_tree_resids is a NumericVector, do nothing
                }else{																// if f+1^th element of other_tree_resids is not a NumericVector
                	throw std::range_error("other resids not a numeric matrix!");	// throw an error
                }
                NumericVector resids_prevroundtemp_tau=other_tree_resids_tau[f];						// create vector residstoadd equal to f+1^th element of other_tree_resids
                NumericVector residstoadd_tau=resids_prevroundtemp_tau-z*curr_round_preds_tau(_,k);
                
                if(is<NumericMatrix>(other_mat_tau[f])){								// if f+1^th element of other_mat is a NumericMatrix, do nothing
                  
                }else{																// if f+1^th element of other_mat is not a NumericMatrix
                  throw std::range_error(" other mat not a numeric matrix!");		// throw an error
                }
                NumericMatrix mattoadd_tau=other_mat_tau[f];								// create matrix mattoadd equal to f+1^th element of other_mat
                
                sum_of_trees_tau[count]=treetoadd_tau;										// add treetoadd to sum_of_trees
                sum_of_tree_resids_tau[count]=residstoadd_tau;								// add residstoadd to sum_of_tree_resids
                sum_of_trees_mat_tau[count]=mattoadd_tau;									// add mattoadd to sum_of_trees_mat
                count++;															// inremet the count variable.
                
                if(count==(size_mat-1)){											// If list size is too small
                  size_mat=size_mat*2;											// double the size
                  sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);				// double the length of sum_of_trees
                  sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,size_mat);	// double the length of sum_of_tree_resids
                  sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);		// double the length of sum_of_trees_mat
                }
              }
              
              if(other_tree_mu.size()>size_mat){											// If list size is too small
                sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,other_tree_mu.size());	// increase the length of sum_of_tree_resids
              }
              for(int f=0;f<other_tree_mu.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
                NumericVector resids_prevroundtemp_mu=other_tree_resids_mu[f];
                NumericVector residstoadd_mu=resids_prevroundtemp_mu-z*curr_round_preds_tau(_,k);						// create vector residstoadd equal to f+1^th element of other_tree_resids
                
                sum_of_tree_resids_mu[f]=residstoadd_mu;								// add residstoadd to sum_of_tree_resids
              }
              
              
              
              
              List temp1_prev_sum_tree = prev_sum_trees_mu[curr_round_parent[k]];
              List temp1_prev_sum_tree_mat = prev_sum_trees_mat_mu[curr_round_parent[k]];
              
              sum_of_trees_mu = resize_bcf(sum_of_trees_mu,temp1_prev_sum_tree.size());			
              sum_of_trees_mat_mu = resize_bcf(sum_of_trees_mat_mu,temp1_prev_sum_tree_mat.size());		
              
              sum_of_trees_mu = prev_sum_trees_mu[curr_round_parent[k]];
              sum_of_trees_mat_mu = prev_sum_trees_mat_mu[curr_round_parent[k]];
            }
          
          // Rcout << "GET TO LINE 6258 .\n";
          
          sum_of_trees_tau[count]=temp_sum_trees_tau[k];										// add k+1^th element of temp_sum_trees to sum_of_trees
          sum_of_tree_resids_tau[count]=temp_sum_tree_resids_tau[k];							// add k+1^th element of temp_sum_tree_resids to sum_of_tree_resids
          sum_of_trees_mat_tau[count]=temp_sum_trees_mat_tau[k];								// add k+1^th element of temp_sum_trees_mat to sum_of_trees_mat
          count++;																	// increment the count variable
          
          // Rcout << "GET TO LINE 6265 .\n";
          // Rcout << "LENGTH OF LIST SUM_OF_TREES_MU = " << sum_of_trees_mu.size() << ".\n";
          // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_MU = " << sum_of_trees_mat_mu.size() << ".\n";
          
          if(count==(size_mat-1)){													// If list size is too small
            size_mat=size_mat*2;													// double the size
            if(j==0) sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);						// double the length of sum_of_trees
            sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);						// double the length of sum_of_trees
            
            if(j==0) sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);				// double the length of sum_of_tree_resids
            sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);				// double the length of sum_of_tree_resids
            sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,size_mat);			// double the length of sum_of_trees_mat
          }
        }
        // Rcout << "count=" << count << ".\n";
        
        if(j==0) sum_of_trees_mu=resize_bcf(sum_of_trees_mu,count);										// remove spaces that are not filled in.
        sum_of_trees_tau=resize_bcf(sum_of_trees_tau,count);										// remove spaces that are not filled in.
        if(j==0) sum_of_trees_mat_mu=resize_bcf(sum_of_trees_mat_mu,count);								// remove spaces that are not filled in.
        sum_of_trees_mat_tau=resize_bcf(sum_of_trees_mat_tau,count);								// remove spaces that are not filled in.
        sum_of_tree_resids_tau=resize_bcf(sum_of_tree_resids_tau,count);							// remove spaces that are not filled in.
        
        List other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
        sum_of_tree_resids_mu=resize_bcf(sum_of_tree_resids_mu,other_tree_mu.size());							// remove spaces that are not filled in.
        
        
        // Rcout << "length of sum of trees = " << sum_of_trees_mu.size() << ".\n";
        // Rcout << "count = " << count << ".\n";
        // Rcout << "length of sum of treestau = " << sum_of_trees_tau.size() << ".\n";
        
        if(curr_round_parent[k]!=-1){															// If the k+1^th element of curr_round_parent is -1, (-1 is a terminal node?)
          overall_sum_trees_mu[overall_count]=sum_of_trees_mu;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_tau[overall_count]=sum_of_trees_tau;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
          overall_sum_tree_resids_mu[overall_count]=sum_of_tree_resids_mu;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
          overall_sum_tree_resids_tau[overall_count]=sum_of_tree_resids_tau;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_mat_mu[overall_count]=sum_of_trees_mat_mu;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_mat_tau[overall_count]=sum_of_trees_mat_tau;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
          overall_sum_BIC=temp_BIC;															// Let overall_sum_BIC equal temp_BIC, the vector of BICs.
          overall_sum_preds_outcome=Rcpp::as<arma::mat>(temp_preds_outcome);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          overall_sum_preds_mu=Rcpp::as<arma::mat>(temp_preds_mu);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          overall_sum_preds_tau=Rcpp::as<arma::mat>(temp_preds_tau);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          if(is_test_data==1){	// If there is test data, overall_sum_test_preds equal temp_test_preds, the matrix of out-of-sample predictions (columns correpond to different models?)
            //Rcout << "Get to Line 6327 in loop j = " << j  << ".\n";
            overall_sum_test_preds_outcome=Rcpp::as<arma::mat>(temp_test_preds_outcome);
            overall_sum_test_preds_mu=Rcpp::as<arma::mat>(temp_test_preds_mu);
            overall_sum_test_preds_tau=Rcpp::as<arma::mat>(temp_test_preds_tau);
          }					
          overall_count++;																	// increment overall_count
          if(overall_count==(overall_size-1)){												// If overall_size is too small
            overall_size=overall_size*2;													// double the size
            overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
            overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
            overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
            overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
            overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
            overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
          }
        }  
      }
    
      
      // Rcout << "overall - length of list of sum of tree models = " << overall_sum_trees_mu.size() << ".\n";
      List example_tree_tab = overall_sum_trees_mu[0];
      // Rcout << "overall - Length of input tree table list used to get best split sum = " << example_tree_tab.size() << ".\n";
      List example_tree_mat = overall_sum_trees_mat_mu[0];
      // Rcout << "overall - Length of input tree mat list used to get best split sum= " << example_tree_mat.size() << ".\n";
      
      //NumericMatrix example_tree_tab_mat = example_tree_tab[0];
      ////Rcout << "number of cols example mat = " << example_tree_tab_mat.ncol() << ".\n";
      
      //NumericMatrix example_tree_tab_mat2 = example_tree_tab[1];
      ////Rcout << "number of cols example mat2 = " << example_tree_tab_mat2.ncol() << ".\n";
      
      if(only_max_num_trees==0){
        
      //check if there were any trees from the previous round that didn't have daughter trees grown.
      //create vector to count number of possible parents for previous round
      if(j>0){																					// If not the first round of the outter loop
        IntegerVector prev_par_no_child=match(prev_par,curr_round_parent);						// create a vector equal to indices of the positions of the (the first matches of the) elements of prev_par in curr_round_parent.
        if(any(is_na(prev_par_no_child))){														// any of the vector of matches are NA (no match?)
          IntegerVector t4=ifelse(is_na(prev_par_no_child),1,0);								// create a vector t4 equal to 1 for the NA values, 0 otherwise.
          for(int h=0;h<prev_par_no_child.size();h++){										// for-loop of length equal to that of prev_par_no_child
            if(t4[h]==1){																	// If h+1^th element of vector of matches is NA
              if(prev_round_BIC2[h]-lowest_BIC<=log(c)){									// If the h+1^th model (from the previous round?) is in Occam's window
                SEXP s_mu = prev_sum_trees_mu[h];												// create a pointer to S expression type equal to the h+1^th element of prev_sum_trees (a tree table or list of tree tables from the previous round?)
                SEXP s_tau = prev_sum_trees_tau[h];												// create a pointer to S expression type equal to the h+1^th element of prev_sum_trees (a tree table or list of tree tables from the previous round?)
                //SEXP s_resid_mu = prev_sum_tree_resids_mu[h]; 
                SEXP s_resid_tau = prev_sum_tree_resids_tau[h];
                
                // if(is<List>(s_resid_mu)){	
                //   if(is<List>(s_resid_tau)){
                //     Rcout << "mu resid list and tau resid list. tau round j= " << j << ".\n"; 
                //   }else{
                //     Rcout << "mu resid list and tau resid vector?. tau round j= " << j << ".\n"; 
                //     
                //   }
                // }else{
                //   if(is<List>(s_resid_tau)){
                //     Rcout << "mu resid vector and tau resid list. tau round j= " << j << ".\n"; 
                //     
                //   }else{
                //     Rcout << "mu resid vector and tau resid vector. tau round j= " << j << ".\n"; 
                //     
                //   }
                // }
                
                
                if(is<List>(s_mu)){														// If prev_sum_trees[h] is a list

                  if(is<List>(s_tau)){
                    List tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    List tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                      NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }                    
                    List treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    List treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 6404 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }else{
                    List tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    NumericMatrix tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                      NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }                    
                    List treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    NumericMatrix treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 6439 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }


                }else{																// If prev_sum_trees[h] is NOT a list
                  if(is<List>(s_tau)){


                    NumericMatrix tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    List tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                      NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }                    
                    NumericMatrix treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    List treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 6480 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }else{
                    NumericMatrix tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                    NumericMatrix tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                    List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                    if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }else{
                      NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    }                    
                    NumericMatrix treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                    NumericMatrix treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                    overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                    overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                    overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                    overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                    overall_count++;													// increment overall_count

                    if(overall_count==(overall_size-1)){												// If overall_size is too small
                      overall_size=overall_size*2;													// double the size
                      overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                    }
                    double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                    overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                    overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                    if(is_test_data==1){
                      //Rcout << "Get to Line 6515 in loop j = " << j  << ".\n";
                      overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                    }
                  }

                }
              }
            }
          }
        }
      }
      }
      
      prev_round_preds_outcome=temp_preds_outcome;															// let prev_round_preds equal to the matrix of predictions.
      prev_round_preds_mu=temp_preds_mu;															// let prev_round_preds equal to the matrix of predictions.
      prev_round_preds_tau=temp_preds_tau;															// let prev_round_preds equal to the matrix of predictions.
      if(is_test_data==1){
        //Rcout << "Get to Line 6534 in loop j = " << j  << ".\n";
        prev_round_test_preds_outcome=temp_test_preds_outcome;								// if there is test data, let prev_round_test_preds equal the test data predictions
        prev_round_test_preds_mu=temp_test_preds_mu;								// if there is test data, let prev_round_test_preds equal the test data predictions
        prev_round_test_preds_tau=temp_test_preds_tau;								// if there is test data, let prev_round_test_preds equal the test data predictions
      }
      prev_round_BIC=temp_BIC;																// let prev_round_BIC equal the vector of BICs
      prev_round_BIC2=temp_BIC;																// let prev_round_BIC2 equal the vector of BICs
      prev_round_preds2_outcome=Rcpp::as<arma::mat>(temp_preds_outcome);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      prev_round_preds2_mu=Rcpp::as<arma::mat>(temp_preds_mu);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      prev_round_preds2_tau=Rcpp::as<arma::mat>(temp_preds_tau);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      if(is_test_data==1){
        //Rcout << "Get to Line 6545 in loop j = " << j  << ".\n";
        prev_round_test_preds2_outcome=Rcpp::as<arma::mat>(temp_test_preds_outcome);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
        prev_round_test_preds2_mu=Rcpp::as<arma::mat>(temp_test_preds_mu);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
        prev_round_test_preds2_tau=Rcpp::as<arma::mat>(temp_test_preds_tau);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
      }
      resids=temp_resids;																		// let resids equal the matrix of residuals
      parent=temp_parent;																		// let parent equal the parent vector
      overall_sum_trees_mu=resize_bcf(overall_sum_trees_mu,overall_count);								// remove spaces that are not filled in.
      overall_sum_trees_tau=resize_bcf(overall_sum_trees_tau,overall_count);								// remove spaces that are not filled in.
      overall_sum_tree_resids_mu=resize_bcf(overall_sum_tree_resids_mu,overall_count);					// remove spaces that are not filled in.
      overall_sum_tree_resids_tau=resize_bcf(overall_sum_tree_resids_tau,overall_count);					// remove spaces that are not filled in.
      overall_sum_trees_mat_mu=resize_bcf(overall_sum_trees_mat_mu,overall_count);						// remove spaces that are not filled in.
      overall_sum_trees_mat_tau=resize_bcf(overall_sum_trees_mat_tau,overall_count);						// remove spaces that are not filled in.
      
      // Rcout << "Check how far get in first round for tau. j = " << j <<" and first round = " << first_round << ".\n";
      
      if(j==0){																		// if in the first round of the outer for-loop (j==0)
        
         //Rcout << "INSIDE IF STATEMENT. j = " << j <<" and first round = " << first_round << ".\n";
        
        if(ntree_control==1){

          List temp_for_prev_tree(temp_sum_trees_tau.size());
          List temp_for_prev_mat(temp_sum_trees_mat_tau.size());
          List temp_for_prev_resids(temp_sum_trees_mat_tau.size());

          for(int p=0;p<temp_sum_trees_tau.size();p++){															// for-loop of length equal to number of models outputted in final round.
            List temp_tree_list(1);
            List temp_mat_list(1);
            List temp_resid_list(1);

            temp_tree_list[0]= temp_sum_trees_tau[p];
            temp_mat_list[0]=temp_sum_trees_mat_tau[p];
            temp_resid_list[0]=temp_sum_tree_resids_tau[p];

            temp_for_prev_tree[p] = temp_tree_list;
            temp_for_prev_mat[p] = temp_mat_list;
            temp_for_prev_resids[p] =temp_resid_list;

          }

            prev_sum_trees_tau = temp_for_prev_tree;
            prev_sum_trees_mat_tau = temp_for_prev_mat;
            prev_sum_tree_resids_tau = temp_for_prev_resids;

            prev_sum_trees_mu=overall_sum_trees_mu;														// let prev_sum_trees equal the list of tree tables (from the current round)
            prev_sum_trees_mat_mu= overall_sum_trees_mat_mu;												// let prev_sum_trees_mat equal the list of tree matrice
            prev_sum_tree_resids_mu=overall_sum_tree_resids_mu;											// let prev_sum_tree_resids equal the list of residual vectors (from the current round)
            
            overall_sum_trees_mu=prev_sum_trees_mu;		// not sure about this					// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
            overall_sum_trees_tau=prev_sum_trees_tau; // not sure about this													// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
            overall_sum_tree_resids_mu=prev_sum_tree_resids_mu;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)
            overall_sum_tree_resids_tau=prev_sum_tree_resids_tau;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)        overall_sum_trees_mat_mu=prev_sum_trees_mat_mu;										// let overall_sum_trees_mat equal the RESIZED list of tree matrice
            overall_sum_trees_mat_mu=prev_sum_trees_mat_mu;
            overall_sum_trees_mat_tau=prev_sum_trees_mat_tau;											// let overall_sum_trees_mat equal the RESIZED list of tree matrice
            if(is_test_data==1){ // if there is test data, let overall_sum_test_preds be an arma mat copy of the test data predictions
              //Rcout << "Get to Line 6601 in loop j = " << j  << ".\n";
              overall_sum_test_preds_outcome= Rcpp::as<arma::mat>(temp_test_preds_outcome);
              overall_sum_test_preds_mu= Rcpp::as<arma::mat>(temp_test_preds_mu);
              overall_sum_test_preds_tau= Rcpp::as<arma::mat>(temp_test_preds_tau);
            }
        }else{
        
        prev_sum_trees_mu=overall_sum_trees_mu;														// let prev_sum_trees equal the list of tree tables (from the current round)
        //ADDING TO MU TREES, therefore nothing yet added to tau trees
        prev_sum_trees_tau = temp_sum_trees_tau;
        prev_sum_trees_mat_mu= overall_sum_trees_mat_mu;												// let prev_sum_trees_mat equal the list of tree matrice
        prev_sum_trees_mat_tau = temp_sum_trees_mat_tau;
        
        
        // Rcout << "length of list of sum of tree models = " << prev_sum_trees_mu.size() << ".\n";
        //List example_tree_tab = prev_sum_trees_mu[0];
        //Rcout << "Length of input tree table list used to get best split sum = " << example_tree_tab.size() << ".\n";
        //List example_tree_mat = prev_sum_trees_mat_mu[0];
        //Rcout << "Length of input tree mat list used to get best split sum= " << example_tree_mat.size() << ".\n";
        
        //NumericMatrix example_tree_tab_mat = example_tree_tab[0];
        //Rcout << "number of cols example mat = " << example_tree_tab_mat.ncol() << ".\n";
        
        //NumericMatrix example_tree_tab_mat2 = example_tree_tab[1]; //CAUSES FATAL ERROR
        //Rcout << "number of cols example mat2 = " << example_tree_tab_mat2.ncol() << ".\n";

        
        prev_sum_tree_resids_tau=temp_sum_tree_resids_tau;											// let prev_sum_tree_resids equal the list of residual vectors (from the current round)
        prev_sum_tree_resids_mu=overall_sum_tree_resids_mu;											// let prev_sum_tree_resids equal the list of residual vectors (from the current round)
        //NumericMatrix test=prev_sum_trees[0];												// create a matrix test equal to the first element of temp_sum_trees (first obtained tree table)
        overall_sum_trees_mu=resize_bcf(prev_sum_trees_mu,temp_sum_trees_tau.size());						// remove spaces that are not filled in.
        overall_sum_trees_tau=resize_bcf(temp_sum_trees_tau,temp_sum_trees_tau.size());// INTENTIONALLY USING SIZE OF MU						// remove spaces that are not filled in.
        overall_sum_trees_mu=prev_sum_trees_mu;		// not sure about this					// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
        overall_sum_trees_tau=temp_sum_trees_tau; // not sure about this													// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
        overall_sum_tree_resids_mu=resize_bcf(prev_sum_tree_resids_tau,temp_sum_tree_resids_tau.size());	// remove spaces that are not filled in.
        overall_sum_tree_resids_mu=prev_sum_tree_resids_mu;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)
        overall_sum_tree_resids_tau=resize_bcf(temp_sum_tree_resids_tau,temp_sum_tree_resids_tau.size());	// remove spaces that are not filled in.
        overall_sum_tree_resids_tau=temp_sum_tree_resids_tau;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)        overall_sum_trees_mat_mu=prev_sum_trees_mat_mu;										// let overall_sum_trees_mat equal the RESIZED list of tree matrice
        overall_sum_trees_mat_tau=temp_sum_trees_mat_tau;											// let overall_sum_trees_mat equal the RESIZED list of tree matrice
        overall_sum_trees_mat_mu=resize_bcf(prev_sum_trees_mat_mu,temp_sum_trees_mat_tau.size());				// remove spaces that are not filled in.
        overall_sum_trees_mat_tau=resize_bcf(temp_sum_trees_mat_tau,temp_sum_trees_mat_tau.size());				// remove spaces that are not filled in.
        overall_sum_BIC=temp_BIC;															// let overall_sum_BIC equal the vector of BICs
        overall_sum_preds_outcome= Rcpp::as<arma::mat>(temp_preds_outcome);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
        overall_sum_preds_mu= Rcpp::as<arma::mat>(temp_preds_mu);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
        overall_sum_preds_tau= Rcpp::as<arma::mat>(temp_preds_tau);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
        if(is_test_data==1){ // if there is test data, let overall_sum_test_preds be an arma mat copy of the test data predictions
          //Rcout << "Get to Line 6647 in loop j = " << j  << ".\n";
          overall_sum_test_preds_outcome= Rcpp::as<arma::mat>(temp_test_preds_outcome);
          overall_sum_test_preds_mu= Rcpp::as<arma::mat>(temp_test_preds_mu);
          overall_sum_test_preds_tau= Rcpp::as<arma::mat>(temp_test_preds_tau);
        }
        }
      }else{																					// if not in the first round of the outer for-loop. (i.e. j>0)
        prev_sum_trees_mu=overall_sum_trees_mu;													// let prev_sum_trees equal the list of tree tables (up to the current round??)
        prev_sum_trees_tau=overall_sum_trees_tau;													// let prev_sum_trees equal the list of tree tables (up to the current round??)
        prev_sum_tree_resids_mu=overall_sum_tree_resids_mu;										// let prev_sum_tree_resids equal the list of residual vectors (up tp the current round)
        prev_sum_tree_resids_tau=overall_sum_tree_resids_tau;										// let prev_sum_tree_resids equal the list of residual vectors (up tp the current round)
        prev_sum_trees_mat_mu=overall_sum_trees_mat_mu;											// let prev_sum_trees_mat equal the list of tree matrice (up to the current round??)
        prev_sum_trees_mat_tau=overall_sum_trees_mat_tau;											// let prev_sum_trees_mat equal the list of tree matrice (up to the current round??)
        prev_round_BIC2=overall_sum_BIC;													// let prev_round_BIC2 equal the vector of BICs (up to the current round??)
        prev_round_preds2_outcome=overall_sum_preds_outcome;												// let prev_round_preds2 equal the predictions matrix
        prev_round_preds2_mu=overall_sum_preds_mu;												// let prev_round_preds2 equal the predictions matrix
        prev_round_preds2_tau=overall_sum_preds_tau;												// let prev_round_preds2 equal the predictions matrix
        if(is_test_data==1){
          //Rcout << "Get to Line 6665 in loop j = " << j  << ".\n";
          prev_round_test_preds2_outcome=overall_sum_test_preds_outcome;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
          prev_round_test_preds2_mu=overall_sum_test_preds_mu;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
          prev_round_test_preds2_tau=overall_sum_test_preds_tau;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
        }
      }
      
      overall_overall_sum_trees_mu[oo_count]=overall_sum_trees_mu;									// Add the current outer loop's table list to overall_overall_sum_trees (list of lists?? OR list of lists of lists??)
      overall_overall_sum_trees_tau[oo_count]=overall_sum_trees_tau;									// Add the current outer loop's table list to overall_overall_sum_trees (list of lists?? OR list of lists of lists??)
      overall_overall_sum_tree_resids_mu[oo_count]=overall_sum_tree_resids_mu;						// Add the current outer loop's list of residual vectors to the list overall_overall_sum_tree_resids (list of lists)
      overall_overall_sum_tree_resids_tau[oo_count]=overall_sum_tree_resids_tau;						// Add the current outer loop's list of residual vectors to the list overall_overall_sum_tree_resids (list of lists)
      overall_overall_sum_trees_mat_mu[oo_count]=overall_sum_trees_mat_mu;							// Add the current outer loop's list of matrices to the list overall_overall_sum_trees_mat (list of lists)
      overall_overall_sum_trees_mat_tau[oo_count]=overall_sum_trees_mat_tau;							// Add the current outer loop's list of matrices to the list overall_overall_sum_trees_mat (list of lists)
      overall_overall_sum_BIC[oo_count]=overall_sum_BIC;										// Add the vector of BICs to the list overall_overall_sum_BIC. (list of vectors?)
      oo_count ++;																					// increment the count 
      if(oo_count==(oo_size-1)){																		// If lists are not large enough
        oo_size=oo_size*2;																			// double the size.
        overall_overall_sum_trees_mu=resize_bigger_bcf(overall_overall_sum_trees_mu,oo_size);					// double the length of overall_overall_sum_trees
        overall_overall_sum_trees_tau=resize_bigger_bcf(overall_overall_sum_trees_tau,oo_size);					// double the length of overall_overall_sum_trees
        overall_overall_sum_tree_resids_mu=resize_bigger_bcf(overall_overall_sum_tree_resids_mu,oo_size);		// double the length of overall_overall_sum_tree_resids
        overall_overall_sum_tree_resids_tau=resize_bigger_bcf(overall_overall_sum_tree_resids_tau,oo_size);		// double the length of overall_overall_sum_tree_resids
        overall_overall_sum_trees_mat_mu=resize_bigger_bcf(overall_overall_sum_trees_mat_mu,oo_size);			// double the length of overall_overall_sum_trees_mat
        overall_overall_sum_trees_mat_tau=resize_bigger_bcf(overall_overall_sum_trees_mat_tau,oo_size);			// double the length of overall_overall_sum_trees_mat
        overall_overall_sum_BIC=resize_bigger_bcf(overall_overall_sum_BIC,oo_size);						// double the length of overall_overall_sum_BIC
      }    
      overall_overall_sum_preds_outcome=overall_sum_preds_outcome;											// let overall_overall_sum_preds equal the prediction matrix.
      overall_overall_sum_preds_mu=overall_sum_preds_mu;											// let overall_overall_sum_preds equal the prediction matrix.
      overall_overall_sum_preds_tau=overall_sum_preds_tau;											// let overall_overall_sum_preds equal the prediction matrix.
      
      if(is_test_data==1){
        //Rcout << "Get to Line 6695 in loop j = " << j  << ".\n";
        overall_overall_sum_test_preds_outcome=overall_sum_test_preds_outcome;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
        overall_overall_sum_test_preds_mu=overall_sum_test_preds_mu;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
        overall_overall_sum_test_preds_tau=overall_sum_test_preds_tau;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
      }
      //overall_trees_mu[j]=curr_round_trees_mu;														// let the j+1^th element of the list overall_trees be the list of trees obtained in the current round of the outer loop (list of lists of tree table matrices)
      //overall_trees_tau[j]=curr_round_trees_tau;														// let the j+1^th element of the list overall_trees be the list of trees obtained in the current round of the outer loop (list of lists of tree table matrices)
      //overall_mat_mu.push_back(curr_round_mat_mu);													// append the list of current round of tree matrices to the list overall_mat (list of lists of matrices)
      //overall_mat_tau.push_back(curr_round_mat_tau);													// append the list of current round of tree matrices to the list overall_mat (list of lists of matrices)
      overall_lik.push_back(curr_round_lik);													// append the current round's vector of BICs to the list overall_lik (list of vectors of BICs).
      prev_par=seq_len(overall_sum_trees_tau.size())-1;											// let prev_par equal a sequence 0,1,2,3,..., up to length of overall_sum_trees minus one
      
      
      
    }	// END OF tau(x) TREE CODE	(if-statement)
    
    
    // Rcout << "Get to end of loop j = " << j << ".\n";
    
  } //	END OF OUTER LOOP
  
  if(oo_count==0){
    throw std::range_error("BCF-BMA did not find any suitable model for the data. Maybe limit for Occam's window is too small. Maybe use more observations or change parameter values.");
  }
  
  overall_overall_sum_trees_mu=resize_bcf(overall_overall_sum_trees_mu,oo_count);						// remove spaces that are not filled in.
  overall_overall_sum_trees_tau=resize_bcf(overall_overall_sum_trees_tau,oo_count);						// remove spaces that are not filled in.
  overall_overall_sum_tree_resids_mu=resize_bcf(overall_overall_sum_tree_resids_mu,oo_count);			// remove spaces that are not filled in.
  overall_overall_sum_tree_resids_tau=resize_bcf(overall_overall_sum_tree_resids_tau,oo_count);			// remove spaces that are not filled in.
  overall_overall_sum_trees_mat_mu=resize_bcf(overall_overall_sum_trees_mat_mu,oo_count);				// remove spaces that are not filled in.
  overall_overall_sum_trees_mat_tau=resize_bcf(overall_overall_sum_trees_mat_tau,oo_count);				// remove spaces that are not filled in.
  overall_overall_sum_BIC=resize_bcf(overall_overall_sum_BIC,oo_count);							// remove spaces that are not filled in.
  NumericVector end_BIC=overall_overall_sum_BIC[overall_overall_sum_BIC.size()-1] ;			// final element of overall_overall_sum_BIC (vector of BICs from final round)
  NumericMatrix overallpreds_outcome(n,end_BIC.size());												// create a vector of dimensions: number of training obs by number of models outptted by the final round.
  NumericMatrix overallpreds_mu(n,end_BIC.size());												// create a vector of dimensions: number of training obs by number of models outptted by the final round.
  NumericMatrix overallpreds_tau(n,end_BIC.size());												// create a vector of dimensions: number of training obs by number of models outptted by the final round.
  NumericMatrix overall_test_preds_outcome(test_data.nrow(),end_BIC.size());							// create a vector of dimensions: number of test obs by number of models outptted by the final round.
  NumericMatrix overall_test_preds_mu(test_data.nrow(),end_BIC.size());							// create a vector of dimensions: number of test obs by number of models outptted by the final round.
  NumericMatrix overall_test_preds_tau(test_data.nrow(),end_BIC.size());							// create a vector of dimensions: number of test obs by number of models outptted by the final round.
  NumericVector post_weights(end_BIC.size());													// create a vector of length equal to number of models outputted by final round.
  for(int k=0;k<end_BIC.size();k++){															// for-loop of length equal to number of models outputted in final round.
    NumericMatrix oosp_outcome=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_preds_outcome));				// create a matrix equal to overall_overall_sum_preds, the training prediction matrix
    NumericMatrix oosp_mu=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_preds_mu));				// create a matrix equal to overall_overall_sum_preds, the training prediction matrix
    NumericMatrix oosp_tau=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_preds_tau));				// create a matrix equal to overall_overall_sum_preds, the training prediction matrix
    NumericVector temp_preds_outcome=oosp_outcome(_,k);														// let temp_preds equal the predictions from the k+1^th model
    NumericVector temp_preds_mu=oosp_mu(_,k);														// let temp_preds equal the predictions from the k+1^th model
    NumericVector temp_preds_tau=oosp_tau(_,k);														// let temp_preds equal the predictions from the k+1^th model

        
    
    NumericVector temp_test_preds_outcome;															// create a vector called temp_test_preds
    NumericVector temp_test_preds_mu;															// create a vector called temp_test_preds
    NumericVector temp_test_preds_tau;															// create a vector called temp_test_preds
    if(is_test_data==1){																	// if there is test data
      //Rcout << "Get to Line 6749  "  << ".\n";
      NumericMatrix oostp_outcome=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_test_preds_outcome));		// let oostp equal the test data prediction matrix
      NumericMatrix oostp_mu=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_test_preds_mu));		// let oostp equal the test data prediction matrix
      NumericMatrix oostp_tau=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_test_preds_tau));		// let oostp equal the test data prediction matrix
      temp_test_preds_outcome=oostp_outcome(_,k);																// let temp_test_preds equal the out of sample predictions from the k+1^th model.
      temp_test_preds_mu=oostp_mu(_,k);																// let temp_test_preds equal the out of sample predictions from the k+1^th model.
      temp_test_preds_tau=oostp_tau(_,k);																// let temp_test_preds equal the out of sample predictions from the k+1^th model.
    }
    NumericVector orig_temp_preds_outcome=get_original_bcf(min(y),max(y),-0.5,0.5,temp_preds_outcome) ;			// Rescale the in-sample predictions back to the original scale of the outcome. Defined on line 2216
    NumericVector orig_temp_preds_mu=get_original_bcf(min(y),max(y),-0.5,0.5,temp_preds_mu) ;			// Rescale the in-sample predictions back to the original scale of the outcome. Defined on line 2216
    NumericVector orig_temp_preds_tau=get_original_bcf(min(y),max(y),-0.5,0.5,temp_preds_tau) ;			// Rescale the in-sample predictions back to the original scale of the outcome. Defined on line 2216
    NumericVector BICi=-0.5*end_BIC;														// create a vector of the BICs multiplied by -0.5
    //Rcout << "end_BIC = "<< end_BIC << " .\n";
    // Rcout << "end_BIC[0] = "<< end_BIC[0] << " .\n";
    
    // Rcout << "TEST LINE 6689 weight =  "<< weight << " .\n";
    
    double max_BIC=max(BICi);																// set the variable max_BIC equal to the maximum of the (negative 0.5 times the) BICs
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));						// create the weight for the k+1^th model
    
    post_weights[k]=weight;																	// Let the k+1^th element of post_weights be the weight of the k+1^th model
    overallpreds_outcome(_,k) = temp_preds_outcome*weight;													// Let the k+1^th element of overallpreds be the predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    overallpreds_mu(_,k) = temp_preds_mu*weight;													// Let the k+1^th element of overallpreds be the predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    overallpreds_tau(_,k) = temp_preds_tau*weight;													// Let the k+1^th element of overallpreds be the predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    if(is_test_data==1){																	// if there is test data
      //Rcout << "Get to Line 6774  "  << ".\n";
      overall_test_preds_outcome(_,k) = temp_test_preds_outcome*weight;											// Let the k+1^th element of overall_test_preds be the out-of-sample predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
      overall_test_preds_mu(_,k) = temp_test_preds_mu*weight;											// Let the k+1^th element of overall_test_preds be the out-of-sample predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
      overall_test_preds_tau(_,k) = temp_test_preds_tau*weight;											// Let the k+1^th element of overall_test_preds be the out-of-sample predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    }
  }
  arma::mat M1_outcome(overallpreds_outcome.begin(), overallpreds_outcome.nrow(), overallpreds_outcome.ncol(), false);						// M1 is an arma mat copy of overallpreds (entry i,j gives contribution of j^th model to prediction of i^th observation)
  arma::mat M1_mu(overallpreds_mu.begin(), overallpreds_mu.nrow(), overallpreds_mu.ncol(), false);						// M1 is an arma mat copy of overallpreds (entry i,j gives contribution of j^th model to prediction of i^th observation)
  arma::mat M1_tau(overallpreds_tau.begin(), overallpreds_tau.nrow(), overallpreds_tau.ncol(), false);						// M1 is an arma mat copy of overallpreds (entry i,j gives contribution of j^th model to prediction of i^th observation)
  predicted_values_outcome=sum(M1_outcome,1);																					// predicted_values is a vector of predictions for each observation in the training data (before inverse scaling). (row sums of overallpreds)
  predicted_values_mu=sum(M1_mu,1);																					// predicted_values is a vector of predictions for each observation in the training data (before inverse scaling). (row sums of overallpreds)
  predicted_values_tau=sum(M1_tau,1);																					// predicted_values is a vector of predictions for each observation in the training data (before inverse scaling). (row sums of overallpreds)
 
  arma::mat M2_outcome(overall_test_preds_outcome.begin(), overall_test_preds_outcome.nrow(), overall_test_preds_outcome.ncol(), false);		// M2 is arma copy of overall_test_preds (entry i,j gives contribution of j^th model to prediction of i^th test observation)
  arma::mat M2_mu(overall_test_preds_mu.begin(), overall_test_preds_mu.nrow(), overall_test_preds_mu.ncol(), false);		// M2 is arma copy of overall_test_preds (entry i,j gives contribution of j^th model to prediction of i^th test observation)
  arma::mat M2_tau(overall_test_preds_tau.begin(), overall_test_preds_tau.nrow(), overall_test_preds_tau.ncol(), false);		// M2 is arma copy of overall_test_preds (entry i,j gives contribution of j^th model to prediction of i^th test observation)
  if(is_test_data==1){
    //Rcout << "Get to Line 6791  "  << ".\n";
    predicted_test_values_outcome=sum(M2_outcome,1);														// if there is test data, predicted_test_values is the vector of final test data predictions (before inverse scaling).
    predicted_test_values_mu=sum(M2_mu,1);														// if there is test data, predicted_test_values is the vector of final test data predictions (before inverse scaling).
    predicted_test_values_tau=sum(M2_tau,1);														// if there is test data, predicted_test_values is the vector of final test data predictions (before inverse scaling).
  }
  if(overall_lik.size()==0){																					// if the length of overall_lik is zero 
    throw std::range_error("BART-BMA didnt find any suitable model for the data. Maybe limit for Occam's window is too small.");
  }else{																										// if the length of overall_lik is greater than zero
    NumericVector orig_preds_outcome=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_values_outcome)) ;					// inverse scale in-sample predictions to original scale
    NumericVector orig_preds_mu=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_values_mu)) ;					// inverse scale in-sample predictions to original scale
    NumericVector orig_preds_tau=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_values_tau)) ;					// inverse scale in-sample predictions to original scale
    NumericVector orig_test_preds_outcome;																			// create vector
    NumericVector orig_test_preds_mu;																			// create vector
    NumericVector orig_test_preds_tau;																			// create vector
    if(is_test_data==1){																					// if have test data
      //Rcout << "Get to Line 6806  "  << ".\n";
      orig_test_preds_outcome=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_test_values_outcome)) ;					// inverse scale out-of-sample predictions to original scale
      orig_test_preds_mu=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_test_values_mu)) ;					// inverse scale out-of-sample predictions to original scale
      orig_test_preds_tau=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_test_values_tau)) ;					// inverse scale out-of-sample predictions to original scale
    }
    NumericVector minmax(2);					// create vector minmax of length 2
    minmax[0]=min(y);							// set first element equal to min(y)
    minmax[1]=max(y);							// set first element equal to max(y)
    if(is_test_data==1){						// if there is test data
      List ret(13);									// list of length 6.
      ret[0] = orig_preds_outcome;							// The first element is a vector of in-sample predictions
      ret[1] = orig_preds_mu;							// The first element is a vector of in-sample predictions
      ret[2] = orig_preds_tau;							// The first element is a vector of in-sample predictions
      ret[3] = overall_overall_sum_trees_mu;			// the second element is a list of lists (of lists?) of tree tables
      ret[4] = overall_overall_sum_trees_tau;			// the second element is a list of lists (of lists?) of tree tables
      ret[5] =overall_overall_sum_trees_mat_mu;		// the third element is a list of lists of tree matrices.
      ret[6] =overall_overall_sum_trees_mat_tau;		// the third element is a list of lists of tree matrices.
      ret[7] = end_BIC;								// the fourth element is the vector of BICs of sum-of-tree-models
      ret[8] = orig_test_preds_outcome;						// the fifth element is the vector of out-of-sample predictions
      ret[9] = orig_test_preds_mu;						// the fifth element is the vector of out-of-sample predictions
      ret[10] = orig_test_preds_tau;						// the fifth element is the vector of out-of-sample predictions
      ret[11] =overall_overall_sum_tree_resids_mu;		// the sixth element if a list of residual vectors
      ret[12] =overall_overall_sum_tree_resids_tau;		// the sixth element if a list of residual vectors
      return(ret);							// return the list
    }else{										// if there is no test data
      List ret(10);									// list of length 5.
      ret[0] = orig_preds_outcome;							// The first element is a vector of in-sample predictions
      ret[1] = orig_preds_mu;							// The first element is a vector of in-sample predictions
      ret[2] = orig_preds_tau;							// The first element is a vector of in-sample predictions
      ret[3] = overall_overall_sum_trees_mu;				// the second element is a list of lists (of lists?) of tree tables
      ret[4] = overall_overall_sum_trees_tau;				// the second element is a list of lists (of lists?) of tree tables
      ret[5] = overall_overall_sum_trees_mat_mu;			// the third element is a list of lists of tree matrices.
      ret[6] = overall_overall_sum_trees_mat_tau;			// the third element is a list of lists of tree matrices.
      ret[7] = end_BIC;								// the fourth element is the vector of BICs of sum-of-tree-models
      ret[8] =overall_overall_sum_tree_resids_mu;		// the sixth element if a list of residual vectors
      ret[9] =overall_overall_sum_tree_resids_tau;		// the sixth element if a list of residual vectors
      return(ret);									// return the list
    }
  }
}
//###########################################################################################################################//
// [[Rcpp::export]]
IntegerVector fuse_intvecs(IntegerVector a, IntegerVector b) {
  std::vector<int> x;                     
  x.insert( x.end(), a.begin(), a.end() );
  x.insert( x.end(), b.begin(), b.end() );
  return(wrap(x));
}
//###########################################################################################################################//
// [[Rcpp::export]]
NumericVector fuse_numvecs(NumericVector a, NumericVector b) {
  std::vector<double> x;                     
  x.insert( x.end(), a.begin(), a.end() );
  x.insert( x.end(), b.begin(), b.end() );
  return(wrap(x));
}
//###########################################################################################################################//
// [[Rcpp::export]]
List fuse_lists(List a, List b) {
  int xlength = a.size()+b.size();
  List x(xlength);
  for(int k=0;k< a.size();k++){
    x[k]=a[k];
  }
  for(int k=0;k< b.size();k++){
    x[a.size()+k]=b[k];
  }
  return(wrap(x));
}
//###########################################################################################################################//
// [[Rcpp::export]]
List fuse_listsof_intvecs(List a, List b) {
  std::list<IntegerVector> x;                     
  x.insert( x.end(), a.begin(), a.end() );
  x.insert( x.end(), b.begin(), b.end() );
  return(wrap(x));
}
//###########################################################################################################################//
// [[Rcpp::export]]
NumericMatrix mmult1(NumericMatrix a, NumericMatrix b) {
  int acoln = a.ncol();
  int bcoln = b.ncol();
  NumericMatrix out = no_init_matrix(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      out(_, j) = b(_, j - acoln);
    }
  }
  return out;
}
//###########################################################################################################################//

//' @title Obtain BCFBMA predictions, trees, BICs etc. to be called by R functions. Add mu or tau tree at each step
//' @export
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List BCF_BMA_sumLikelihood_add_mu_or_tau(NumericMatrix data,NumericVector y, NumericVector z, NumericMatrix pihat,
                           double a_mu,double a_tau,double mu_mu,double mu_tau,double nu,double lambda,double c,
                           double sigma_mu_mu,double sigma_mu_tau,double pen_mu,double pen_tau,int num_cp_mu,int num_cp_tau,
                           NumericMatrix test_data,NumericVector test_z,NumericMatrix test_pihat,
                           int ntree_control,int ntree_moderate,
                           double alpha_mu,double alpha_tau,double beta_mu,double beta_tau,bool split_rule_node,bool gridpoint,int maxOWsize,
                           int num_splits_mu,int num_splits_tau,int gridsize_mu, int gridsize_tau, int include_pi2,
                           bool zero_split, bool only_max_num_trees,int separate_tree_numbers
                                           ){
  // Rcout << "LINE 9172.\n";
  
  bool is_test_data=0;					// create bool is_test_data. Initialize equal to 0.
  if(test_data.nrow()>0){					// If test data has non-zero number of rows.
    is_test_data=1;						// set is_test_data equal to 1.
  }
  if(y.size() !=data.nrow()){				// If the length of input vector y is not equal to the nunber of rows in the input data (covariates)
    if(y.size()<data.nrow()){			// If the length of y is less than the number of rows in data
      throw std::range_error("Response length is smaller than the number of observations in the data"); 
    }else{								// If the length of y is greater than the number of rows in data
      throw std::range_error("Response length is greater than the number of observations in the data"); 
    }
  }
  //Rcout << "LINE 7678.\n";
  
  if(z.size() !=data.nrow()){				// If the length of input vector z is not equal to the nunber of rows in the input data (covariates)
    if(z.size()<data.nrow()){			// If the length of z is less than the number of rows in data
      throw std::range_error("Treatment indicator vector length is smaller than the number of observations in the data"); 
    }else{								// If the length of z is greater than the number of rows in data
      throw std::range_error("Treatment indicator vector length is greater than the number of observations in the data"); 
    }
  }
  if(pihat.nrow() !=data.nrow()){				// If the nunber of rows in the input matrix pihat is not equal to the nunber of rows in the input data (covariates)
    if(pihat.nrow()<data.nrow()){			// If the nunber of rows in the input matrix pihat is less than the number of rows in data
      throw std::range_error("The nunber of rows in the input matrix pihat is smaller than the number of observations in the data"); 
    }else{								// If the nunber of rows in the input matrix pihat is greater than the number of rows in data
      throw std::range_error("The nunber of rows in the input matrix pihat is greater than the number of observations in the data"); 
    }
  }
  //Rcout << "LINE 7694.\n";
  
  //check test data has the same number of variables as training data
  if(test_data.nrow()>0 && (data.ncol() != test_data.ncol())){	// If the number of rows in the test data is >0 AND the number of columns (variables) is not equal to that of data (the training data)
    throw std::range_error("Test data and training data must have the same number of variables. BART BMA assumes variables are in the same order."); 
  }
  if(test_z.size() != test_data.nrow()){	// If the number of rows in the test data covariate matrix is not equal to that of the test data treatment indicator variable
    throw std::range_error("Test data covariates and test data treatment indicator variable must have the same number of observations."); 
  }
  if(test_data.nrow() != test_pihat.nrow()){	// If the number of rows in the test data covariate matrix is not equal to that of the test data propensity score estimates matrix
    throw std::range_error("Test data covariates and test data propensity score estimates must have the same number of observations."); 
  }
  if(test_pihat.nrow()>0 && (pihat.ncol() != test_pihat.ncol())){	// If the number of rows in the test data propensity score estimates is >0 AND the number of columns (variables) is not equal to that of the training data propensity score estimates
    throw std::range_error("Test data propensity score estimates and training data propensity score estimates must have the same number of columns. BART BMA assumes variables are in the same order."); 
  }
  //check value of c is greater than 1!
  //	if(c<1){
  //		throw std::range_error("Value of Occam's Window has to be greater than 0."); 
  //	}
  if(num_cp_mu<0 || num_cp_mu>100){		// If input num_cp_mu is <0 or >100
    throw std::range_error("Value of num_cp_mu should be a value between 1 and 100."); 
  }
  if(num_cp_tau<0 || num_cp_tau>100){		// If input num_cp_tau is <0 or >100
    throw std::range_error("Value of num_cp_tau should be a value between 1 and 100."); 
  }
  // Now add propensity score estimates matrix as new leftmost column of data matrix. Call the resulting matrix x_control (to be consistent with terminology used by bcf package).
  arma::mat D1(data.begin(), data.nrow(), data.ncol(), false);				// copy the covariate data matrix into an arma mat
  arma::mat pihat_a(pihat.begin(), pihat.nrow(), pihat.ncol(), false);				// copy the pihat matrix into an arma mat
  arma::mat x_control_a=D1;				// create a copy of data arma mat called x_control_a
  if((include_pi2==0) | (include_pi2==2) ){
    if(pihat.nrow()>0 ){
      x_control_a.insert_cols(0,pihat_a);		// add propensity scores as new leftmost columns of x_control_a
    }
  }
   //Rcout << "Number of columns of matrix" << x_control_a.n_cols << ".\n";
  NumericMatrix x_control=wrap(x_control_a);	// convert x_control_a to a NumericMatrix called x_control
  // Name the matrix without the estimated propensity scores x_moderate.[CAN REMOVE THE DUPLICATION AND ADD x_control, x_moderate, and include_pi as input parameters later]
  //NumericMatrix x_moderate = data;	// x_moderate matrix is the covariate data without the propensity scores
  arma::mat x_moderate_a=D1;			// create arma mat copy of x_moderate.
  if((include_pi2==1)| (include_pi2==2) ){
    if(pihat.nrow()>0 ){
      x_moderate_a.insert_cols(0,pihat_a);		// add propensity scores as new leftmost columns of x_control_a
    }
  }
  NumericMatrix x_moderate=wrap(x_moderate_a);	// convert x_control_a to a NumericMatrix called x_control
   //Rcout << "Get to Line 7739  "  << ".\n";
  // Add test propensity scores to test data matrix
  arma::mat T1(test_data.begin(), test_data.nrow(), test_data.ncol(), false);				// copy the covariate test_data matrix into an arma mat
  arma::mat pihat_a_test(test_pihat.begin(), test_pihat.nrow(), test_pihat.ncol(), false);				// copy the test_pihat matrix into an arma mat
  arma::mat x_control_test_a=T1;				// create a copy of test_data arma mat called x_control_test_a
  if((include_pi2==0)| (include_pi2==2) ){
    if(test_pihat.nrow()>0 ){
      x_control_test_a.insert_cols(0,pihat_a_test);		// add propensity scores as new leftmost columns of x_control_test_a
    }
  }
  NumericMatrix x_control_test=wrap(x_control_test_a);	// convert x_control_test_a to a NumericMatrix called x_control_test
  // Name the matrix without the estimated propensity scores x_moderate_test.[CAN REMOVE THE DUPLICATION AND ADD x_control_test, x_moderate_test, and include_pi as input parameters later]
  //NumericMatrix x_moderate_test = test_data;	// x_moderate_test matrix is the covariate test_data without the propensity scores
  arma::mat x_moderate_test_a=T1;			// create arma mat copy of x_moderate_test.
  if((include_pi2==1)| (include_pi2==2) ){
    if(test_pihat.nrow()>0 ){
      x_moderate_test_a.insert_cols(0,pihat_a_test);		// add propensity scores as new leftmost columns of x_control_a
    }
  }
  NumericMatrix x_moderate_test=wrap(x_moderate_test_a);	// convert x_control_test_a to a NumericMatrix called x_control_test
  
   //Rcout << "Get to Line 7754  "  << ".\n";
  
  // NOT SURE IF SEPARATE INITIAL TREE MATRIX REQUIRED FOR mu(x) and tau(x)
  // BUT STILL DESIRABLE TO END UP WITH SEPARATE LISTS AND MATRICES FOR mu(x) and tau(x) trees
  //
  
  
  //	NumericMatrix treetable=start_tree_bcf(mu,sigma_mu);								// create matrix treetable (defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
  //	NumericMatrix treemat=start_matrix_bcf(data.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  //initialize the tree table and matrix
  
  NumericMatrix treetable_mu=start_tree_bcf(mu_mu,sigma_mu_mu);								// create matrix treetable (defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
  NumericMatrix treemat_mu=start_matrix_bcf(x_control.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  
  //	MIGHT NEED JUST NUMBER OF TREATED OBSERVATIONS for treemat_tau
  NumericMatrix treetable_tau=start_tree_bcf(mu_tau,sigma_mu_tau);								// create matrix treetable (defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
  
  // 	
  
  // create matrix of data for treated households only (not sure if this will be used)
  //	arma::mat x_mod_treat_a = x_moderate_a.rows(arma::find(z_ar==1));	
  //	NumericMatrix x_mod_treat = wrap(x_mod_treat_a);	
  //	NumericMatrix treemat_tau=start_matrix_bcf(x_mod_treat.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  
  // FOR NOW TRY TO PUT ALL OBSERVATIONS IN THE TAU TREE MATRIX OUTPUT
  NumericMatrix treemat_tau=start_matrix_bcf(x_moderate.nrow());								// line 623. Returns a data.nrow() by 1 matrix with all elements equal to 1
  
  
  //	MAYBE INCLUDE treetable_tau and treemat_tau HERE
  
  NumericVector y_scaled=scale_response_bcf(min(y),max(y),-0.5,0.5,y);		// re-scale the outcome variable
  double n=D1.n_rows;					// number of rows of the data matrix (number of training observations)
  //	CHANGE ABOVE LINE IF REMOVE data and D1
  
  // BELOW FUNCTIONS AND LINES ONLY USED DIRECTLY HERE FOR OBTAINING lowest_BIC
  //	//	double lik=likelihood_function_bcf(y_scaled,treetable,treemat,a,mu,nu,lambda);
  //	double lik=likelihood_function_bcf(y_scaled,treetable,treemat,a,mu,nu,lambda);	// defined on line 201. Returns log likelihood of initial model with no trees?
  //	double tree_prior=get_tree_prior_bcf(treetable,treemat,alpha,beta);				// defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees? but the tree is empty  at this stage?)
  //	double lowest_BIC=-2*(lik+log(tree_prior))+1*log(n);						// BIC actually not the BIC because add 2*log(prior). Closer to the Bayes factor? (Not noted in paper).
  
  //	THIS IS PROBABLY WRONG: SHOULD START WITH MODEL Y=a+bZ, WHERE THERE IS A mu(x) stump AND tau(x) STUMP.]
  // 	currently: will add tau stump after first mu(x) tree grown
  double lik=likelihood_function_mu_bcf(y_scaled,treetable_mu,treemat_mu,a_mu,mu_mu,nu,lambda);	// defined on line 201. Returns log likelihood of initial model with no trees?
  if(treetable_mu.ncol()<5) throw std::range_error("Line 4081");
  double tree_prior=get_tree_prior_bcf(treetable_mu,treemat_mu,alpha_mu,beta_mu);				// defined on line 566. Presumably returns a prior probability. (prior for single tree or sum of trees? but the tree is empty  at this stage?)
  double lowest_BIC=-2*(lik+log(tree_prior))+1*log(n);						// BIC actually not the BIC because add 2*log(prior). Closer to the Bayes factor? (Not noted in paper).
  
  // create initial tree table lists for mu(x) and tau(x)
  // Need to check if doing this correctly for tau(x)
  List tree_table_mu;								// create tree table.
  List tree_mat_mu;									// create tree matrix.
  tree_table_mu.push_back(treetable_mu);				// append treetable to list tree_table. (adding the model with no trees as the first element of the list).
  tree_mat_mu.push_back(treemat_mu);					// append treemat to list tree_mat. (adding the model with no trees as the first element of the list).
  
  List tree_table_tau;								// create tree table.
  List tree_mat_tau;									// create tree matrix.
  // these two lines might be unnecessary
  tree_table_tau.push_back(treetable_tau);				// append treetable to list tree_table. (adding the model with no trees as the first element of the list).
  tree_mat_tau.push_back(treemat_tau);					// append treemat to list tree_mat. (adding the model with no trees as the first element of the list).
  
  
  //	Not sure if should create separate CART_BMA_mu and CART_BMA_tau lists or just write over when moving from mu trees to tau trees
  // Probably unnecessary to have separate BART_BMA, but can edit this later
  List CART_BMA_mu;									// create list.
  List CART_BMA_tau;									// create list.
  
  
  arma::mat r;									// create matrix r.
  arma::colvec yarma=clone(y_scaled);				// create arma colvec copy, yarma, of scaled outcome vector.
  r.insert_cols(0,yarma);							// let the first column of r be yarma (the scaled outcome vector).
  NumericMatrix resids=wrap(r);					// create NumericMatrix copy of r called resids.
  
  
  // FOR NOW, TRYING TO LEAVE ALL SUBSETTING TO TREATED RESIDUALS TO FUNCTIONS (e.g. get_best_trees_sum_tau_bcf)	
  // create vector of treated outcome observations only?
  //	arma::mat yarma_treated = yarma.elem(arma::find(z_ar==1));	
  //	create vector of treated resids only?
  //	perhaps unnecessary if first should remove prediction from first mu tree
  //	arma::mat r_treated;									// create matrix r.
  //	r_treated.insert_cols(0,yarma_treated);					// let the first column of r be yarma (the scaled outcome vector).
  //	NumericMatrix resids=wrap(r_treated);					// create NumericMatrix copy of r called resids.
  
  
  
  //int first_round;								// create a variable first_round. Not initialized.
  
  //	NOT SURE IF SHOULD KEEP ALL OF THESE
  //	List overall_trees(max(ntree_control,ntree_moderate));					// create a list of length equal to the input value num_rounds.
  //	List overall_mat;								// create a list.
  List overall_lik;								// create a list.
  
  List overall_trees_mu(ntree_control);					// create a list of length equal to the input value num_rounds.
  List overall_mat_mu;								// create a list.
  //	List overall_lik_mu;								// create a list.
  
  List overall_trees_tau(ntree_moderate);					// create a list of length equal to the input value num_rounds.
  List overall_mat_tau;								// create a list.
  //	List overall_lik_tau;								// create a list.
  
  
  NumericMatrix prev_round_preds_outcome;  				// create a matrix.
  NumericMatrix prev_round_preds_mu;  				// create a matrix.
  NumericMatrix prev_round_preds_tau;  				// create a matrix.
  
  NumericVector prev_round_BIC;					// create a vector.
  NumericVector prev_round_BIC2;					// create a vector.
  arma::mat prev_round_preds2_outcome;					// create an arma mat.
  arma::mat prev_round_preds2_mu;					// create an arma mat.
  arma::mat prev_round_preds2_tau;					// create an arma mat.
  NumericMatrix prev_round_test_preds_outcome;			// create a matrix. (default values are zeroes?)
  NumericMatrix prev_round_test_preds_mu;			// create a matrix. (default values are zeroes?)
  NumericMatrix prev_round_test_preds_tau;			// create a matrix. (default values are zeroes?)
  
  arma::mat prev_round_test_preds2_outcome;				// create an arma mat
  arma::mat prev_round_test_preds2_mu;				// create an arma mat
  arma::mat prev_round_test_preds2_tau;				// create an arma mat
  
  arma::mat overall_overall_sum_test_preds_outcome;		// create an arma mat
  arma::mat overall_overall_sum_test_preds_mu;		// create an arma mat
  arma::mat overall_overall_sum_test_preds_tau;		// create an arma mat
  
  arma::colvec predicted_test_values_outcome;				// create an arma colvec
  arma::colvec predicted_test_values_mu;				// create an arma colvec
  arma::colvec predicted_test_values_tau;				// create an arma colvec
  
  
  //	Not sure if should keep these, or separately have lists for mu and tau trees, resids, mat	
  //	List prev_sum_trees;							// create a list
  // List prev_sum_tree_resids;						// create a list
  //	List prev_sum_trees_mat;						// create a list
  //	List cp_mat_list;								// create a list
  
  List prev_sum_trees_mu;								// create a list
  List prev_sum_tree_resids_mu;						// create a list
  List prev_sum_trees_mat_mu;							// create a list
  List cp_mat_list_mu;								// create a list
  
  // Create a vector recording the number of mu trees in all models
  // This vector should be of the same length as prev_sum_trees_mu
  // The purpose is to determine when the maximum number of mu trees is reached 
  // then will obtain a subset of models with less than the maximum before getting new mu trees to add 
  IntegerVector prev_mu_trees_count;
  
  List prev_sum_trees_tau;							// create a list
  List prev_sum_tree_resids_tau;						// create a list
  List prev_sum_trees_mat_tau;						// create a list
  List cp_mat_list_tau;								// create a list	

  // Create a vector recording the number of tau trees in all models
  // This vector should be of the same length as prev_sum_trees_mu
  // The purpose is to determine when the maximum number of tau trees is reached 
  // then will obtain a subset of models with less than the maximum before getting new tau trees to add   
  IntegerVector prev_tau_trees_count;
  
  
  
  
  int oo_size=300;								// create a variable. Initialized equal to 300.
  List overall_overall_sum_tree_resids_mu(oo_size);	// create a list of length 300.
  List overall_overall_sum_tree_resids_tau(oo_size);	// create a list of length 300.
  
  List overall_overall_sum_BIC(oo_size);			// create a list of length 300.
  int oo_count=0;									// create a variable, Initialized equal to 0.
  
  List overall_overall_sum_trees_mu(oo_size);			// create a list of length 300.
  List overall_overall_sum_trees_mat_mu(oo_size);		// create a list of length 300.
  
  List overall_overall_sum_trees_tau(oo_size);		// create a list of length 300.
  List overall_overall_sum_trees_mat_tau(oo_size);	// create a list of length 300.
  
  arma::mat overall_overall_sum_preds_outcome;			// create an arma mat.
  arma::mat overall_overall_sum_preds_mu;			// create an arma mat.
  arma::mat overall_overall_sum_preds_tau;			// create an arma mat.
  
  //	parent indexes a whole sum-of-tree model including mu(x) and tau(x) trees
  IntegerVector prev_par;							// create a vector
  
  arma::colvec predicted_values_outcome;					// create an arma colvec
  arma::colvec predicted_values_mu;					// create an arma colvec
  arma::colvec predicted_values_tau;					// create an arma colvec
  
  
  for(int j=0;j<ntree_control+ntree_moderate;j++){					// create a for-loop of length equal to the input value num_rounds.
      // Rcout << "Beginning of loop number = " << j << ".\n";
      //Rcout << "ntree_control = " << ntree_control << ".\n";
      //Rcout << "ntree_moderate = " << ntree_moderate << ".\n";
    
    //if(j<ntree_control){
      if(j>0){
        if(only_max_num_trees==1){
          lowest_BIC=100000;
        }
      }
      int overall_size=300;						// create a variable overall_size. Initialized equal to 0.
      List overall_sum_tree_resids_mu(overall_size);	// create a list of length 300
      List overall_sum_tree_resids_tau(overall_size);	// create a list of length 300
      
      List overall_sum_trees_mu(overall_size);		// create a list of length 300
      List overall_sum_trees_mat_mu(overall_size);	// create a list of length 300
      
      List overall_sum_trees_tau(overall_size);		// create a list of length 300
      List overall_sum_trees_mat_tau(overall_size);	// create a list of length 300
      
        //Rcout << "Get to Line 7957 in loop j = " << j << ".\n";
      int overall_count=0;						// set overall_count equal to 0.
      
      
      
      // CREATE VECTORS TO RECORD IF ADDED TREE IS A MU TREE OR TAU TREE FOR EACH MODEL 
      // EQUALS 1 IF MU TREE, 2 IF TAU TREE (maybe 0 should throw an error)
      
      // VECTOR FOR MODELS WITH A MU TREE ADDED
      IntegerVector curr_round_added_mu_addmu;
      // VECTOR FOR MODELS WITH A TAU TREE ADDED
      IntegerVector curr_round_added_mu_addtau;
      // VECTOR FOR ALL PROPOSED MODELS
      //IntegerVector curr_round_added_mu_addall;
      
      // CREATE NEW VECTORS TO RECORD COUNTS OF MU AND TAU TREES FOR PROPOSED SUM OF TREE MODELS
      // VECTORS FOR MODELS WITH A MU TREE ADDED
      IntegerVector new_mu_tree_counts_addmu;
      IntegerVector new_tau_tree_counts_addmu;
      // VECTORS FOR MODELS WITH A TAU TREE ADDED
      IntegerVector new_mu_tree_counts_addtau;
      IntegerVector new_tau_tree_counts_addtau;
      // VECTORS FOR ALL PROPOSED MODELS
      IntegerVector new_mu_tree_counts_addall;
      IntegerVector new_tau_tree_counts_addall;
      
      //		parent indexes the whole models to which the a tree can be appended
      IntegerVector parent;						// create vector.
      
      
      

      
      
      // CREATE VECTORS FOR OUTPUR OF get best trees
      
      // FOR THE MODELS WITH MU TREES ADDED
      NumericVector curr_round_lik_addmu;				// create vector	// To be filled with BICs for whole models suggested after a tree appended
      List curr_round_trees_mu_addmu;						// create list.
      List curr_round_trees_tau_addmu;						// create list.
      
      List curr_round_mat_mu_addmu;						// create list.
      List curr_round_mat_tau_addmu;						// create list.
      
      NumericVector curr_BIC_addmu;						// create vector.
      //		next line probably shouldn't need to be duplicated for curr_round_parent_mu and curr_round_parent_tau
      IntegerVector curr_round_parent_addmu;			// create vector.
      
      
      // FOR THE MODELS WITH TAU TREES ADDED
      NumericVector curr_round_lik_addtau;				// create vector	// To be filled with BICs for whole models suggested after a tree appended
      List curr_round_trees_mu_addtau;						// create list.
      List curr_round_trees_tau_addtau;						// create list.
      
      List curr_round_mat_mu_addtau;						// create list.
      List curr_round_mat_tau_addtau;						// create list.
      
      NumericVector curr_BIC_addtau;						// create vector.
      //		next line probably shouldn't need to be duplicated for curr_round_parent_mu and curr_round_parent_tau
      IntegerVector curr_round_parent_addtau;			// create vector.
      
      
      
      
      NumericVector overall_sum_BIC;				// create vector.
      
      arma::mat overall_sum_preds_outcome;				// create matrix.
      arma::mat overall_sum_preds_mu;				// create matrix.
      arma::mat overall_sum_preds_tau;				// create matrix.
      
      arma::mat overall_sum_test_preds_outcome;			// create matrix.
      arma::mat overall_sum_test_preds_mu;			// create matrix.
      arma::mat overall_sum_test_preds_tau;			// create matrix.
      
      
      
      
      
      
       //Rcout << "Get to Line 8062 in loop j = " << j << ".\n";
      if(j==0){									// If in the first round of the for-loop.
        parent.push_back(0);					// append a 0 to the end of the parent vector. (first and only element of parent vector so far).
        //first_round=1;							// set the variable first_round equal to 1.
      }else{										// If not in the first round of the for-loop.
        //first_round=0;							// set the variable first_round equal to 0.
      }
      //		The _mu	in resids_cp_mat_mu is probably unnecessary, but including it to remove ambiguity.
      //		Replace with just resids_cp_mat for memory efficiency after code is all working.
      //		similarly, cp_mat_list_mu is used (perhaps unnecessarily) here instead of cp_mat_list
      List resids_cp_mat_mu(resids.ncol());												// create a list of length equal to the number of colmns of resids
      int resids_count=0;																// create a variable resids_count. Initialize equal to zero.
      std::vector<int> err_list(resids.ncol());										// create a vector err_list of length equal to the number of colmns of resids
      std::vector<int> no_more_mu_trees(resids.ncol());										// create a vector err_list of length equal to the number of colmns of resids
      
      //get best splits
      for(int f=0;f<resids.ncol();f++){												// for-loop of length equal to the unmber of columns of resids
         
        if(j>0 && separate_tree_numbers==1){
          List temp_mutrees = prev_sum_trees_mu[f];
          if(temp_mutrees.size()==ntree_control){
            NumericMatrix temp_mat_err(0,3);
            resids_cp_mat_mu[resids_count]= temp_mat_err;
            err_list[resids_count]= 1;
            no_more_mu_trees[resids_count]= 1;
            
          }else{
            if(gridpoint==0){															// If input gridpoint equals 0. i.e. the PELT method will be used.
              cp_mat_list_mu=make_pelt_cpmat_mu_bcf(x_control,resids(_,f),pen_mu,num_cp_mu);				// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
            }else{																		// If input gridpoint equals 1. i.e. the PELT method will be used.
              cp_mat_list_mu=make_gridpoint_cpmat_mu_bcf(x_control,resids(_,f),gridsize_mu,num_cp_mu);			// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
            }
            resids_cp_mat_mu[resids_count]=cp_mat_list_mu[0];									// let the resids_count+1^th element be a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points.
            err_list[resids_count]=cp_mat_list_mu[1];
            no_more_mu_trees[resids_count]= 0;
          }
        }else{
        if(gridpoint==0){															// If input gridpoint equals 0. i.e. the PELT method will be used.
          cp_mat_list_mu=make_pelt_cpmat_mu_bcf(x_control,resids(_,f),pen_mu,num_cp_mu);				// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }else{																		// If input gridpoint equals 1. i.e. the PELT method will be used.
          cp_mat_list_mu=make_gridpoint_cpmat_mu_bcf(x_control,resids(_,f),gridsize_mu,num_cp_mu);			// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
        }
        resids_cp_mat_mu[resids_count]=cp_mat_list_mu[0];									// let the resids_count+1^th element be a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. 
        err_list[resids_count]=cp_mat_list_mu[1];										// let the resids_count+1^th element be a list that records error(s?)
        no_more_mu_trees[resids_count]= 0;
        }
        resids_count++;																// increment resids_count by 1.
      }
      resids_cp_mat_mu=resize_bcf(resids_cp_mat_mu,resids_count);								// remove elements of resids_cp_mat_mu that are not fileld in. How is  it possible that elements are not filled in??
      err_list.resize(resids_count);													// remove elements of err_list that are not fileld in. How is  it possible that elements are not filled in??
      no_more_mu_trees.resize(resids_count);													// remove elements of err_list that are not fileld in. How is  it possible that elements are not filled in??
      parent=seq_len(tree_table_mu.size())-1;											// set parent equal to a vector 0,1,2,3,..., up to the length of tree_table minus one.
      // Rcout << "ORIGINAL parent = " << parent << ".\n"; 
      // Rcout << "ORIGINAL tree_table_mu.size() = " << tree_table_mu.size() << ".\n"; 
      // Rcout << "ORIGINAL prev_sum_trees_mu.size() = " << prev_sum_trees_mu.size() << ".\n"; 
      // Rcout << "ORIGINAL prev_sum_trees_tau.size() = " << prev_sum_trees_tau.size() << ".\n"; 
      
      
      
      // if(is_true(all(as<IntegerVector>(wrap(err_list))==1))){							// if all elements of err_list equal 1.
      //   if(j==0){																	// If in the first round of the for-loop.
      //     throw std::range_error("No split points could be found to grow trees");
      //   }else{																		// If not in the first round of the for-loop.
      //     Rcout << "j= " << j << ".\n";
      //     throw std::range_error("No Mu trees can be grown for the number of iterations desired, as no splits were found.Please try fewer iterations.");
      //   }
      // } 
      
      
      
      
       // Rcout << "Get to Line 9620 in loop j = " << j << ".\n";
      //get current set of trees.
      if(j==0){						// If in the first round of the for-loop.
        CART_BMA_mu=get_best_trees_mu_bcf_2(x_control_a, x_moderate_a,z,resids,
                                          a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                          tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                          lowest_BIC,//first_round,
                                          parent,resids_cp_mat_mu,as<IntegerVector>(wrap(err_list)),
                                          x_control_test,x_moderate_test,test_z,
                                          alpha_mu,alpha_tau,beta_mu,beta_tau,
                                          is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                          split_rule_node,gridpoint,maxOWsize,num_splits_mu,num_splits_tau,gridsize_mu,zero_split,
                                          as<IntegerVector>(wrap(no_more_mu_trees)));
        

      
      }else{							// If not in the first round of the for-loop.
        //if j >0 then sum of trees become a list so need to read in list and get likelihood for each split point and terminal node
        
        CART_BMA_mu=get_best_trees_sum_mu_bcf_2(x_control_a, x_moderate_a,z,resids,
                                              a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                              tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                              lowest_BIC,//first_round,
                                              parent,resids_cp_mat_mu,as<IntegerVector>(wrap(err_list)),
                                              x_control_test,x_moderate_test,test_z,
                                              alpha_mu,alpha_tau,beta_mu,beta_tau,
                                              is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                              split_rule_node,gridpoint,maxOWsize,
                                              prev_sum_trees_mu,prev_sum_trees_tau,prev_sum_trees_mat_mu,
                                              prev_sum_trees_mat_tau,y_scaled,num_splits_mu,num_splits_tau,gridsize_mu,zero_split,
                                              as<IntegerVector>(wrap(no_more_mu_trees)));	// function defined on line 1953.
      }
       // Rcout << "Get to after get best trees in line 9650 in loop j = " << j << ".\n";
      
      curr_round_lik_addmu=CART_BMA_mu[0];							// vector of BICs (for whole sum-of tree-models after suggested trees added). Should be ordered ascending
      curr_round_trees_mu_addmu=CART_BMA_mu[1];						// list of tree tables
      //curr_round_trees_tau=CART_BMA_mu[2];
      curr_round_mat_mu_addmu=CART_BMA_mu[2];							// list of tree matrices
      //curr_round_mat_tau=CART_BMA_mu[4];							// list of tree matrices
      curr_round_parent_addmu=CART_BMA_mu[3];						// vector of tree parent numbers
      NumericMatrix curr_round_preds_mu=CART_BMA_mu[4];			// (in-sample single tree predictions for trees to add) matrix rows correspond to different units/individuals, columns corresponds to predictions from different (single or sums-of?) trees.
      curr_BIC_addmu=CART_BMA_mu[5];								// lowest BIC among trees?
      NumericMatrix curr_round_test_preds_mu=CART_BMA_mu[6];	// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspond to predictions from different (single or sums-of?) trees.
      
      
     
      
      
      curr_round_added_mu_addmu = rep(1,curr_round_lik_addmu.size());

      
      
    
      // INSERT CODE FOR ADDING TO MU TREE COUNTS HERE
      
      // VECTORS FOR MODELS WITH A MU TREE ADDED
      // new_mu_tree_counts_addmu;
      // new_tau_tree_counts_addmu;
      
      
      
      // INSERT CODE TO 
      // REMOVE MODELS WITH TOO MANY TAU TREES
      
      
      
      
      if(j==0){									// If in the first round of the for-loop.
        parent.push_back(0);					// append a 0 to the end of the parent vector. (first and only element of parent vector so far).
        //first_round=1;							// set the variable first_round equal to 1.
      }else{										// If not in the first round of the for-loop.
        //first_round=0;							// set the variable first_round equal to 0.
      }
      //		The _tau in resids_cp_mat_tau is probably unnecessary, but including it to remove ambiguity.
      //		Replace with just resids_cp_mat for memory efficiency after code is all working.
      //		similarly, cp_mat_list_tau is used (perhaps unnecessarily) here instead of cp_mat_list
      List resids_cp_mat_tau(resids.ncol());												// create a list of length equal to the number of colmns of resids
      int resids_count2=0;																// create a variable resids_count2. Initialize equal to zero.
      std::vector<int> err_list2(resids.ncol());										// create a vector err_list2 of length equal to the number of colmns of resids
      std::vector<int> no_more_tau_trees(resids.ncol());										// create a vector err_list of length equal to the number of colmns of resids
      //get best splits
      for(int f=0;f<resids.ncol();f++){												// for-loop of length equal to the unmber of columns of resids
        
        if(j>0 && separate_tree_numbers==1){
            List temp_tautrees = prev_sum_trees_tau[f];
            if(temp_tautrees.size()==ntree_moderate){
              NumericMatrix temp_mat_err(0,3);
              resids_cp_mat_tau[resids_count2]= temp_mat_err;
              err_list2[resids_count2]= 1;
              no_more_tau_trees[resids_count2]= 1;
            }else{
              if(gridpoint==0){															// If input gridpoint equals 0. i.e. the PELT method will be used.
                cp_mat_list_tau=make_pelt_cpmat_tau_bcf(x_moderate,resids(_,f),gridsize_tau,num_cp_tau,z);				// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
              }else{																		// If input gridpoint equals 1. i.e. the PELT method will be used.
                cp_mat_list_tau=make_gridpoint_cpmat_tau_bcf(x_moderate,resids(_,f),gridsize_tau,num_cp_tau,z);			// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
              }
              resids_cp_mat_tau[resids_count2]=cp_mat_list_tau[0];									// let the resids_count2+1^th element be a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. 
              err_list2[resids_count2]=cp_mat_list_tau[1];										// let the resids_count2+1^th element be a list that records error(s?)
              no_more_tau_trees[resids_count2]= 0;
              
            }
          }else{
        
            if(gridpoint==0){															// If input gridpoint equals 0. i.e. the PELT method will be used.
              cp_mat_list_tau=make_pelt_cpmat_tau_bcf(x_moderate,resids(_,f),gridsize_tau,num_cp_tau,z);				// make_pelt_cpmat defined on line 1612. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
            }else{																		// If input gridpoint equals 1. i.e. the PELT method will be used.
              cp_mat_list_tau=make_gridpoint_cpmat_tau_bcf(x_moderate,resids(_,f),gridsize_tau,num_cp_tau,z);			// make_gridpoint_cpmat defined on line 1550. First element of list is a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. Second element of list records error.
            }
            resids_cp_mat_tau[resids_count2]=cp_mat_list_tau[0];									// let the resids_count2+1^th element be a matrix, where the first column is a list of variable numbers for potential splits (column of variable in data), and the second column gives values of the covariates for the split points. 
            err_list2[resids_count2]=cp_mat_list_tau[1];										// let the resids_count2+1^th element be a list that records error(s?)
            no_more_tau_trees[resids_count2]= 0;
            
          }
        resids_count2++;																// increment resids_count2 by 1.
      }
      resids_cp_mat_tau=resize_bcf(resids_cp_mat_tau,resids_count2);								// remove elements of resids_cp_mat_tau that are not fileld in. How is  it possible that elements are not filled in??
      err_list2.resize(resids_count2);													// remove elements of err_list2 that are not fileld in. How is  it possible that elements are not filled in??
      no_more_tau_trees.resize(resids_count2);													// remove elements of err_list2 that are not fileld in. How is  it possible that elements are not filled in??
      parent=seq_len(tree_table_tau.size())-1;											// set parent equal to a vector 0,1,2,3,..., up to the length of tree_table minus one.
      
      
      // if(is_true(all(as<IntegerVector>(wrap(err_list2))==1))){							// if all elements of err_list2 equal 1.
      //   if(j==0){																	// If in the first round of the for-loop.
      //     throw std::range_error("No split points could be found to grow trees");
      //   }else{																		// If not in the first round of the for-loop.
      //     throw std::range_error("No Tau trees can be grown for the number of iterations desired, as no splits were found.Please try fewer iterations.");
      //   }
      // }
      
      
      
        // Rcout << "Get to 9749 in loop j = " << j << ".\n";
      
      //get current set of trees.
      if(j==0){						// If in the first round of the for-loop.
        //CART_BMA_tau=get_best_trees_tau(x_control_a, x_moderate_a,z,resids,
        //a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
        //tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
        //lowest_BIC,first_round,parent,resids_cp_mat_tau,as<IntegerVector>(wrap(err_list2)),
        //x_control_test,x_moderate_test,test_z,
        //alpha_mu,alpha_tau,beta_mu,beta_tau,
        //is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
        //split_rule_node,gridpoint,maxOWsize);
        CART_BMA_tau=get_best_trees_tau_bcf(x_control_a, x_moderate_a,z,resids,
                                            a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,
                                            sigma_mu_mu,sigma_mu_tau,
                                            tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                            lowest_BIC,//first_round,
                                            parent,resids_cp_mat_tau,as<IntegerVector>(wrap(err_list2)),
                                            x_control_test,x_moderate_test,test_z,
                                            alpha_mu,alpha_tau,beta_mu,beta_tau,
                                            is_test_data,
                                            pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                            split_rule_node,gridpoint,maxOWsize,
                                            num_splits_mu,num_splits_tau,gridsize_tau,zero_split,
                                            as<IntegerVector>(wrap(no_more_tau_trees)));	// function defined on line 1953.
      }else{							// If not in the first round of the for-loop.
        //if j >0 then sum of trees become a list so need to read in list and get likelihood for each split point and terminal node
        CART_BMA_tau=get_best_trees_sum_tau_bcf_2(x_control_a, x_moderate_a,z,resids,
                                                a_mu,a_tau,mu_mu,mu_tau,nu,lambda,c,sigma_mu_mu,sigma_mu_tau,
                                                tree_table_mu,tree_mat_mu,tree_table_tau,tree_mat_tau,
                                                lowest_BIC,//first_round,
                                                parent,resids_cp_mat_tau,as<IntegerVector>(wrap(err_list2)),
                                                x_control_test,x_moderate_test,test_z,
                                                alpha_mu,alpha_tau,beta_mu,beta_tau,
                                                is_test_data,pen_mu,num_cp_mu,pen_tau,num_cp_tau,	// some of these arguments are probably unnecessary
                                                split_rule_node,gridpoint,maxOWsize,
                                                prev_sum_trees_mu,prev_sum_trees_tau,
                                                prev_sum_trees_mat_mu,prev_sum_trees_mat_tau,y_scaled,
                                                num_splits_mu,num_splits_tau,gridsize_tau,zero_split,
                                                as<IntegerVector>(wrap(no_more_tau_trees)));	// function defined on line 1953.
      }
        // Rcout << "Get to 9789 in loop j = " << j << ".\n";
      
      curr_round_lik_addtau=CART_BMA_tau[0];							// vector of BICs (for whole sum-of tree-models after suggested trees added). Should be ordered ascending
      //curr_round_trees_mu_addtau=CART_BMA_tau[1];						// list of tree tables
      curr_round_trees_tau_addtau=CART_BMA_tau[1];					
      //curr_round_mat_mu_addtau=CART_BMA_tau[3];							// list of tree matrices
      curr_round_mat_tau_addtau=CART_BMA_tau[2];							// list of tree matrices
      curr_round_parent_addtau=CART_BMA_tau[3];						// vector of tree parent numbers
      NumericMatrix curr_round_preds_tau=CART_BMA_tau[4];			// (in-sample single tree predictions for trees to add) matrix rows correspond to different units/individuals, columns corresponds to predictions from different (single or sums-of?) trees.
      curr_BIC_addtau=CART_BMA_tau[5];								// lowest BIC among trees?
      NumericMatrix curr_round_test_preds_tau=CART_BMA_tau[6];	// (out-of-sample tree predictions?) matrix rows correspond to different units/individuals, columns correspond to predictions from different (single or sums-of?) trees.
      
      
      
      curr_round_added_mu_addtau = rep(2,curr_round_lik_addtau.size());
      
    
      // INSERT CODE FOR TAU TREE COUNTS HERE
      
      

      // VECTORS FOR MODELS WITH A TAU TREE ADDED
      // new_mu_tree_counts_addtau;
      // new_tau_tree_counts_addtau;
      
      
      // INSERT CODE TO 
      // REMOVE MODELS WITH TOO MANY TAU TREES
      
      
      
      
      
      
      
      if(curr_round_lik_addmu.size()+curr_round_lik_addtau.size()==0) {						// If number of sum of tree models is zero?
         // Rcout << "curr_round_lik.size()==0 BREAK in mu round in loop j = " << j << ".\n";
        //REMOVE THIS ERROR IF WANT TO ALLOW LESS THAN MAX NUMBER OF TREES
        //throw std::range_error("No mu trees chosen in round");
        
        break;											// break out of for-loop
      } 
      
      
      
       //Rcout << "Get to Line 8272 in loop j = " << j << ".\n";
      if(curr_BIC_addmu[0]<lowest_BIC){							// If the lowest BIC obtained by get_best_trees_sum is less than the currently saved lowest value
        lowest_BIC=curr_BIC_addmu[0];							// reset lowest_BIC to the new lowest value
      }
      
      if(curr_BIC_addtau[0]<lowest_BIC){							// If the lowest BIC obtained by get_best_trees_sum is less than the currently saved lowest value
        lowest_BIC=curr_BIC_addtau[0];							// reset lowest_BIC to the new lowest value
      }      
      
      
      
      // MERGE THE OUTPUT FROM get_best mu and tau trees
      
      
      //int num_total_newtrees = curr_round_lik_addmu.size() + curr_round_lik_addtau.size();
      
      //NumericVector curr_round_lik(num_total_newtrees);				// create vector	// To be filled with BICs for whole models suggested after a tree appended
      //List curr_round_trees_mu(num_total_newtrees);						// create list.
      //List curr_round_trees_tau(num_total_newtrees);						// create list.

      //List curr_round_mat_mu(num_total_newtrees);						// create list.
      //List curr_round_mat_tau(num_total_newtrees);						// create list.
      
      //NumericVector curr_BIC(num_total_newtrees);						// create vector.
      //		next line probably shouldn't need to be duplicated for curr_round_parent_mu and curr_round_parent_tau
      //IntegerVector curr_round_parent(num_total_newtrees);			// create vector.
      
      
      IntegerVector curr_round_added_mu_addall = fuse_intvecs(curr_round_added_mu_addmu,curr_round_added_mu_addtau);
      
      NumericVector curr_round_lik = fuse_numvecs(curr_round_lik_addmu,curr_round_lik_addtau);				// create vector	// To be filled with BICs for whole models suggested after a tree appended
      

      
      
      NumericVector curr_BIC = fuse_numvecs(curr_BIC_addmu,curr_BIC_addtau);						// create vector.
      IntegerVector curr_round_parent = fuse_intvecs(curr_round_parent_addmu,curr_round_parent_addtau) ;			// create vector.
      
      // There does not appear to be LISTS to fuse
      // unless put all mu and tau trees into the same lists
      //fuse_lists
      List curr_round_trees = fuse_lists(curr_round_trees_mu_addmu,curr_round_trees_tau_addtau);						// create list.
      List curr_round_mat = fuse_lists(curr_round_mat_mu_addmu,curr_round_mat_tau_addtau);
      
      // Check if the following lines are faster than the two above
      //List curr_round_trees = fuse_listsof_numericmats(curr_round_trees_mu_addmu,curr_round_trees_tau_addtau);						// create list.
      //List curr_round_mat = fuse_listsof_numericmats(curr_round_mat_mu_addmu,curr_round_mat_tau_addtau);
      
      //THERE ARE PERHAPS FASTER METHODS FOR COMBINING LISTS
      //BUT NEED TO USE RCPPARMADILLO OR STANDARD LIBRARY
      
      
      NumericMatrix curr_round_preds = cbind(curr_round_preds_mu,curr_round_preds_tau);
      NumericMatrix curr_round_test_preds = cbind(curr_round_test_preds_mu,curr_round_test_preds_tau);
      //CHECK IF 2 LINES BELOW ARE FASTER than those above
      //NumericMatrix curr_round_preds = mmult1(curr_round_preds_mu,curr_round_preds_tau);
      //NumericMatrix curr_round_test_preds = mmult1(curr_round_test_preds_mu,curr_round_test_preds_tau);
      
      
      
      
      
      // Create matriceas and lists to be used in the next two for-loops of length curr_round_lik.size()
      
      tree_table_mu=List();									// reset tree_table to an empty list
      tree_mat_mu=List();									// reset tree_mat to an empty list.
      tree_table_tau=List();									// reset tree_table to an empty list
      tree_mat_tau=List();									// reset tree_mat to an empty list.
      int lsize=curr_round_lik.size();					// create a variable equal to the number of sum of tree models returned by get_best_trees_sum
      tree_table_mu=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_mat_mu=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_table_tau=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      tree_mat_tau=List(lsize);								// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_outcome(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_mu(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_preds_tau(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      
      NumericMatrix temp_test_preds_outcome(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_test_preds_mu(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      NumericMatrix temp_test_preds_tau(test_data.nrow(),curr_round_lik.size());	// create a matrix of dimensions: number of test observations by the number of sum of tree models returned by get_best_trees_sum
      
      NumericMatrix temp_resids(n,curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericVector temp_parent(curr_round_lik.size());	// create a matrix of dimensions: n (number of training observations) by the number of sum of tree models returned by get_best_trees_sum
      NumericVector temp_BIC(curr_round_lik.size());		// create a vector of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      //List temp_sum_trees_mu(lsize);							// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      //List temp_sum_trees_tau(lsize);							// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      List temp_sum_tree_resids(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      //List temp_sum_tree_resids_mu(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      //List temp_sum_tree_resids_tau(lsize);					// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      //List temp_sum_trees_mat_mu(lsize);						// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      //List temp_sum_trees_mat_tau(lsize);						// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum

      List temp_sum_trees(lsize);							// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      List temp_sum_trees_mat(lsize);						// create a list of length equal to the number of sum of tree models returned by get_best_trees_sum
      
      
       //Rcout << "Get to Line 8371 in loop j = " << j << ".\n";
      
      
      
      
      
      
      
      
      // THE FOLLOWING for-loop CREATES temporary 
      // matrices of predictions and lists of trees, matrices, and partial residuals
      // to be used in the next for-loop of length curr_round_lik.size()
      
      
      int count=0; 										// create a count variable. Initialized equal to zero.
      for(int k=0;k<curr_round_lik.size();k++){
        if(curr_round_added_mu_addall[k]==1){
          
          
        // create a for-loop of length equal to the number of sum of tree models returned by get_best_trees_sum
        tree_table_mu[count]=start_tree_bcf(mu_mu,sigma_mu_mu);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
        tree_mat_mu[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
        tree_table_tau[count]=start_tree_bcf(mu_tau,sigma_mu_tau);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
        tree_mat_tau[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
        if(j==0){										// if in first round of for-loop
          temp_preds_outcome(_,k)=curr_round_preds(_,k);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
          temp_preds_mu(_,k)=curr_round_preds(_,k);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
          NumericVector zerovec1(n,0);
          temp_preds_tau(_,k)=zerovec1;
          if(is_test_data==1){       
            //Rcout << "Get to Line 8401 in loop j = " << j << "k= "<< k << ".\n";
            temp_test_preds_outcome(_,k)=curr_round_test_preds(_,k);
            temp_test_preds_mu(_,k)=curr_round_test_preds(_,k);
            NumericVector zerovectest(test_data.nrow(), 0.0);
            temp_test_preds_tau(_,k)=zerovectest;
          }
          // If there is test data, let the k+1^th column of temp_test_preds be the k+1^th column of curr_round_test_preds. These are the out-of-sample predictions of from the k+1^th model.
          temp_resids(_,k)=y_scaled-temp_preds_outcome(_,k);	// Let the k+1^th column of temp_resids be the outcome minus the predictons from the k+1^th model 
          temp_parent[k]=-1;							// Let the K=1^th element of temp_parent be -1.
          temp_BIC[k]=curr_round_lik[k];				// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
          
          
          //REMOVE ONE OF THE TWO LINES BELOW
   //       temp_sum_trees_mu[count]=curr_round_trees[k];	// Add the tree table to the list temp_sum_trees
          temp_sum_trees[count]=curr_round_trees[k];	// Add the tree table to the list temp_sum_trees
          
          //temp_sum_trees_tau[count]=curr_round_trees_tau[k];	// Add the tree table to the list temp_sum_trees
          
          //REMOVE ONE OF THE TWO LINES BELOW
    //      temp_sum_trees_mat_mu[count]=curr_round_mat[k];	// Add the tree matrix to temp_sum_trees_mat
          temp_sum_trees_mat[count]=curr_round_mat[k];	// Add the tree matrix to temp_sum_trees_mat
          
          
          //temp_sum_trees_mat_tau[count]=curr_round_mat_tau[k];	// Add the tree matrix to temp_sum_trees_mat
          //NOT SURE ABOUT THE FOLLOWING LINE
          
          temp_sum_tree_resids[count]=resids(_,0);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
          
          //temp_sum_tree_resids_mu[count]=resids(_,0);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
          // Rcout << "LENGTH OF resids(_,0)= " << resids(_,0).size() << ".\n";
        }else{											// If not in the first round of the for-loop.
          NumericVector curr_temp_pred_outcome=curr_round_preds(_,k) + 
            prev_round_preds_mu(_,curr_round_parent[k])+ 
            z*prev_round_preds_tau(_,curr_round_parent[k]);	// curr_temp_pred is the sum of the current round predictions and the previous round predictions?? Each round is for one tree? and explain more of the residuals in each round?
          NumericVector curr_temp_pred_mu=curr_round_preds(_,k) + 
            prev_round_preds_mu(_,curr_round_parent[k]);					
          NumericVector curr_temp_pred_tau=prev_round_preds_tau(_,curr_round_parent[k]);	
          
          NumericVector curr_temp_test_pred_outcome;			// create a vector
          NumericVector curr_temp_test_pred_mu;			// create a vector
          NumericVector curr_temp_test_pred_tau;			// create a vector
          
          if(is_test_data==1) {						// If there is test data.
             //Rcout << "Get to Line 8444 in loop j = " << j << "k= "<< k << ".\n";
            curr_temp_test_pred_outcome=curr_round_test_preds(_,k) + prev_round_test_preds_mu(_,curr_round_parent[k])+test_z*prev_round_test_preds_tau(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            curr_temp_test_pred_mu=curr_round_test_preds(_,k) + prev_round_test_preds_mu(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            curr_temp_test_pred_tau=prev_round_test_preds_tau(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
            
            temp_test_preds_outcome(_,k) = curr_temp_test_pred_outcome;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            temp_test_preds_mu(_,k) = curr_temp_test_pred_mu;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            temp_test_preds_tau(_,k) = curr_temp_test_pred_tau;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
            //Rcout << "Get to Line 8452 in loop j = " << j << "k= "<< k << ".\n";
          }
          temp_BIC[k]=curr_round_lik[k];									// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
          temp_preds_outcome(_,k)=curr_temp_pred_outcome;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          temp_preds_mu(_,k)=curr_temp_pred_mu;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          temp_preds_tau(_,k)=curr_temp_pred_tau;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
          
          temp_resids(_,k)=y_scaled-curr_temp_pred_outcome;						// Let the k+1^th column of temp_resids be the new residuals after appending the k+1^th
          temp_parent[k] = k;												// Let k+1^th element of temp_parent equal k.
          
          temp_sum_trees[count]=curr_round_trees[k];						// Add the tree table to the list temp_sum_trees.
          
          //SHOULD PROBABLY REMOVE LINE BELOW
     //     temp_sum_trees_mu[count]=curr_round_trees[k];						// Add the tree table to the list temp_sum_trees.
          //temp_sum_trees_tau[count]=curr_round_trees_tau[k];						// Add the tree table to the list temp_sum_trees.
          
          
          temp_sum_trees_mat[count]=curr_round_mat[k];					// Add the tree matrix to temp_sum_trees_mat.
          
          //SHOULD PROBABLY REMOVE LINE BELOW
      //    temp_sum_trees_mat_mu[count]=curr_round_mat[k];					// Add the tree matrix to temp_sum_trees_mat.
          //temp_sum_trees_mat_tau[count]=curr_round_mat_tau[k];					// Add the tree matrix to temp_sum_trees_mat.
          temp_sum_tree_resids[count]=resids(_,curr_round_parent[k]);		// Add the curr_round_parent[k]+1^th column of resids to temp_sum_tree_resids. resids contains predictions from previous rounds?
          
          //REMOVE LINE BELOW?? RELY ON vector indicating whether a mu tree or a tau tree added?
          //temp_sum_tree_resids_mu[count]=resids(_,curr_round_parent[k]);		// Add the curr_round_parent[k]+1^th column of resids to temp_sum_tree_resids. resids contains predictions from previous rounds?
          
        }
        
        }else{
        if(curr_round_added_mu_addall[k]==2){
          
          
          
          tree_table_mu[count]=start_tree_bcf(mu_mu,sigma_mu_mu);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
          tree_mat_mu[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
          tree_table_tau[count]=start_tree_bcf(mu_tau,sigma_mu_tau);		// add 1 by 7 matrix to the list tree_table. (start_tree_bcf defined on line 602). Returns 1 by 7 matrix with columns ("left daughter","right daughter","split var","split point","status","mean","std dev")) and 1st row values (0,0,0,0,-1,rand,0)
          tree_mat_tau[count]=start_matrix_bcf(n);				// Add matrix to list. start_matrix_bcf defined on line 623. Returns a n (number of obs) by 1 matrix with all elements equal to 1
          if(j==0){										// if in first round of for-loop
            //temp_preds_outcome(_,k)=prev_round_preds_mu(_,curr_round_parent[k])+z*curr_round_preds(_,k);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
            temp_preds_outcome(_,k)=z*curr_round_preds(_,k);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
            NumericVector zerovec2(n,0);
            temp_preds_mu(_,k)=zerovec2;
            //temp_preds_mu(_,k)=prev_round_preds_mu(_,curr_round_parent[k]);		// let k+1^th column of temp_preds equal k+1^th column of curr_round_preds. This is the in-sample predictions from the k+1^th (sum-of-trees) model?
            temp_preds_tau(_,k)=curr_round_preds(_,k);
            
            if(is_test_data==1){
              //Rcout << "Get to Line 8499 in loop j = " << j  << ".\n";
              //temp_test_preds_outcome(_,k)=prev_round_test_preds_mu(_,curr_round_parent[k])+z*curr_round_test_preds(_,k);
              temp_test_preds_outcome(_,k)=z*curr_round_test_preds(_,k);
              NumericVector zerovectest(test_data.nrow(), 0.0);
              temp_test_preds_mu(_,k)=zerovectest;
              //temp_test_preds_mu(_,k)=prev_round_test_preds_mu(_,curr_round_parent[k]);
              temp_test_preds_tau(_,k)=curr_round_test_preds(_,k);
            }
            // If there is test data, let the k+1^th column of temp_test_preds be the k+1^th column of curr_round_test_preds. These are the out-of-sample predictions of from the k+1^th model.
            
            temp_resids(_,k)=y_scaled-temp_preds_outcome(_,k);	// Let the k+1^th column of temp_resids be the outcome minus the predictons from the k+1^th model 
            temp_parent[k]=-1;							// Let the K=1^th element of temp_parent be -1.
            temp_BIC[k]=curr_round_lik[k];				// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
            //temp_sum_trees_mu[count]=curr_round_trees_mu[k];	// Add the tree table to the list temp_sum_trees
            
            temp_sum_trees[count]=curr_round_trees[k];						// Add the tree table to the list temp_sum_trees.
            
            //SHOULD PROBABLY REMOVE THE LINE BELOW
        //    temp_sum_trees_tau[count]=curr_round_trees[k];	// Add the tree table to the list temp_sum_trees
            
            //temp_sum_trees_mat_mu[count]=curr_round_mat_mu[k];	// Add the tree matrix to temp_sum_trees_mat
            
            temp_sum_trees_mat[count]=curr_round_mat[k];					// Add the tree matrix to temp_sum_trees_mat.
            
            //SHOULD PROBABLY REMOVE LINE BELOW
        //    temp_sum_trees_mat_tau[count]=curr_round_mat[k];	// Add the tree matrix to temp_sum_trees_mat
            
            //NOT SURE ABOUT THE FOLLOWING LINE
            //temp_sum_tree_resids[count]=resids(_,0);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
            //temp_sum_tree_resids_mu[count]=y_scaled-z*curr_round_preds(_,k);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
            
            
            temp_sum_tree_resids[count]=resids(_,0);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
            
            //REMOVE LINE BELOW?? RELY ON vector indicating whether a mu tree or a tau tree added?
            //temp_sum_tree_resids_tau[count]=resids(_,curr_round_parent[k]);	// Add to the list temp_sum_tree_resids the first column of resids from the start of the loop. (which, for j=0, is empty? No dimensions?)
            
            
          }else{											// If not in the first round of the for-loop.
            NumericVector curr_temp_pred_outcome=prev_round_preds_mu(_,curr_round_parent[k]) + 
              z*prev_round_preds_tau(_,curr_round_parent[k]) + 
              z*curr_round_preds(_,k);	// curr_temp_pred is the sum of the current round predictions and the previous round predictions?? Each round is for one tree? and explain more of the residuals in each round?
            NumericVector curr_temp_pred_mu=prev_round_preds_mu(_,curr_round_parent[k]);					
            NumericVector curr_temp_pred_tau=prev_round_preds_tau(_,curr_round_parent[k])+curr_round_preds(_,k);	
            
            NumericVector curr_temp_test_pred_outcome;			// create a vector
            NumericVector curr_temp_test_pred_mu;			// create a vector
            NumericVector curr_temp_test_pred_tau;			// create a vector
            
            if(is_test_data==1) {						// If there is test data.
              //Rcout << "Get to Line 8549 in loop j = " << j  << ".\n";
              curr_temp_test_pred_outcome=prev_round_test_preds_mu(_,curr_round_parent[k])+
                test_z*prev_round_test_preds_tau(_,curr_round_parent[k])+
                test_z*curr_round_test_preds(_,k);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
              curr_temp_test_pred_mu=prev_round_test_preds_mu(_,curr_round_parent[k]);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
              curr_temp_test_pred_tau=prev_round_test_preds_tau(_,curr_round_parent[k])+
                curr_round_test_preds(_,k);	// curr_temp_test_pred is the sum of the current round out-of-sample predictions and the previous round out-of-sample predictions?? Each round is for one tree? and explain more of the residuals in each round?
              
              temp_test_preds_outcome(_,k) = curr_temp_test_pred_outcome;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
              temp_test_preds_mu(_,k) = curr_temp_test_pred_mu;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
              temp_test_preds_tau(_,k) = curr_temp_test_pred_tau;	// Set the k+1^th column of temp_test_preds equal to curr_temp_test_pred (new out of sample predictinos from appending k+1^th tree to the sum of tree model?).
              
            }
            temp_BIC[k]=curr_round_lik[k];									// Let the k+1^th element of temp_BIC be the BIC of the k+1^th model.
             //Rcout << "TEST LINE 8563 temp_BIC[k]=  "<< temp_BIC[k] << " .\n";
            
            
            temp_preds_outcome(_,k)=curr_temp_pred_outcome;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
            temp_preds_mu(_,k)=curr_temp_pred_mu;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
            temp_preds_tau(_,k)=curr_temp_pred_tau;									// Let the k+1^th column of temp_preds be the in-samle predictions from adding the k+1^th tree
            
            temp_resids(_,k)=y_scaled-curr_temp_pred_outcome;						// Let the k+1^th column of temp_resids be the new residuals after appending the k+1^th
            temp_parent[k] = k;												// Let k+1^th element of temp_parent equal k.
            //temp_sum_trees_mu[count]=curr_round_trees_mu[k];						// Add the tree table to the list temp_sum_trees.
            
            temp_sum_trees[count]=curr_round_trees[k];						// Add the tree table to the list temp_sum_trees.
            
            //SHOULD PROBABLY REMOVE THE LINE BELOW
        //    temp_sum_trees_tau[count]=curr_round_trees[k];	// Add the tree table to the list temp_sum_trees

            //temp_sum_trees_mat[count]=curr_round_mat_mu[k];					// Add the tree matrix to temp_sum_trees_mat.
            
            temp_sum_trees_mat[count]=curr_round_mat[k];					// Add the tree matrix to temp_sum_trees_mat.
            
            //SHOULD PROBABLY REMOVE LINE BELOW
        //   temp_sum_trees_mat_tau[count]=curr_round_mat[k];	// Add the tree matrix to temp_sum_trees_mat
            
            temp_sum_tree_resids[count]=resids(_,curr_round_parent[k]);		// Add the curr_round_parent[k]+1^th column of resids to temp_sum_tree_resids. resids contains predictions from previous rounds?
            
            //REMOVE LINE BELOW?? RELY ON vector indicating whether a mu tree or a tau tree added?
            //temp_sum_tree_resids_tau[count]=resids(_,curr_round_parent[k]);		// Add the curr_round_parent[k]+1^th column of resids to temp_sum_tree_resids. resids contains predictions from previous rounds?
          }
          
          
          
          
        }else{
          throw std::range_error("Added tree not recorded as mu or tau tree");
        }
        }

        count++;													// Increment the count by 1.(Note this is within the innermost for-loop).
      }
      
      
      if(curr_round_lik.size()==0){									// if no new trees outputted in current round (by get_best_trees_sum? Why not throw this error earlier, at line 2350)
        throw std::range_error("No trees chosen in last round");
      }
        //Rcout << "Get to line 8607 in loop j = " << j << ".\n";
      
      
      
      
      
      
      
      

      // THE FOLLOWING LOOP CREATES THE SUM_OF_TREES LISTS (And matrices and partial resids)
      // FOR EACH OF THE MODELS PROPOSED IN THE ROUND OF THE OUTER LOOP
      // THEN ADDS THEN TO THE overall_sum lists
      
      for(int k=0;k<curr_round_lik.size();k++){	// create a for-loop of length equal to the number of sum of tree models returned by get_best_trees_sum
        
        if(curr_round_added_mu_addall[k]==1){
          
        int size_mat=300;						// create a variable, Initializd equal to 300.
        List sum_of_trees_mu(size_mat);			// create a list of length 300.
        List sum_of_tree_resids_mu(size_mat);		// create a list of length 300.
        List sum_of_tree_resids_tau(size_mat);		// create a list of length 300.
        List sum_of_trees_mat_mu(size_mat);		// create a list of length 300.
        int count=0;							// create a variable count equal to 0. (Count was already define, so could remove "int" at start of this line and just reset cound to 0).
        
         //Rcout << "Get to line 8632 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        List sum_of_trees_tau(size_mat);			// create a list of length 300.
        List sum_of_trees_mat_tau(size_mat);		// create a list of length 300.
        //if (j==1) List sum_of_trees_tau(1);			// create a list of length 300.
        //if (j==1) List sum_of_trees_mat_tau(1);		// create a list of length 300.
        //if (j>1) List sum_of_trees_tau = prev_sum_trees_tau[curr_round_parent[k]];
        //if (j>1) List sum_of_trees_mat_tau = prev_sum_trees_mat_tau[curr_round_parent[k]];
          //Rcout << "Get to line 8639 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        
        if(curr_round_parent[k]==-1){			// If the k+1^th element of curr_round_parent is -1, do nothing. (-1 is a terminal node?)
           //Rcout << "line 8642 curr_round_parent[k]==-1 in loop j = " << j << " and inner loop k = "<< k <<".\n";
           //Rcout << "overall_count =  = " << overall_count << ".\n";
          
          sum_of_trees_mu[count]= temp_sum_trees[k];
          sum_of_tree_resids_mu[count]=temp_sum_tree_resids[k];							// add k+1^th element of temp_sum_tree_resids to sum_of_tree_resids
          sum_of_trees_mat_mu[count]=temp_sum_trees_mat[k];
          count=1;
          
        }else{									// If the k+1^th element of curr_round_parent is not equal to -1.
          
          
          //REMOVED j==1 block of code.
          //Need all prev_sum_trees (and similar things) to be lists of lists of matrices
          //i.e. not just lists of matrices (even if only lists of lists of one matrix each)
          //Also, need to include if-statements for when prev_sum_trees lists are of length zero
          
          // NEED TO THINK MORE ABOUT j==1 CASE. Also need the prev_sum_trees to be correct
          // if(j==1){							// If in the SECOND round of the outter for-loop??
          //   
          //   sum_of_trees_tau = resize_bcf(sum_of_trees_tau,1);			// create a list of length 300.
          //   sum_of_trees_mat_tau = resize_bcf(sum_of_trees_mat_tau,1);		// create a list of length 300.
          //   
          //   // Rcout << "LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
          //   // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
          //   
          //   // Rcout << "Get to 5449 in mu round in loop j = " << j << ".\n";
          //   List other_tree_mulist=prev_sum_trees_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          //   NumericMatrix other_tree_mu=other_tree_mulist[0];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          //   
          //   //NumericMatrix other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          //   // Rcout << "Get to 5454 in mu round in loop j = " << j << ".\n";
          //   
          //   NumericMatrix other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          //   
          //   //NumericVector other_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];	// create vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
          //   
          //   //MIGHT NEED TO CHANGE THIS LINE
          //   List other_resids_mulist=prev_sum_tree_resids_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          //   NumericVector other_resids_mu=other_resids_mulist[0];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          //   // Rcout << "Line 5465 LENGTH OF other_resids_mu = " << other_resids_mu.size() << ".\n";
          //   NumericVector other_resids_tau=prev_sum_tree_resids_tau[curr_round_parent[k]];	// create vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
          //   
          //   sum_of_trees_mu[count]= other_tree_mu;										// add other_tree to the list sum_of_trees
          //   sum_of_trees_tau[count]= other_tree_tau;										// add other_tree to the list sum_of_trees
          //   sum_of_tree_resids_mu[count]=other_resids_mu-curr_round_preds(_,k);									// add other_resids to the list sum_of_tree_resids
          //   sum_of_tree_resids_tau[count]=other_resids_tau-curr_round_preds(_,k);									// add other_resids to the list sum_of_tree_resids
          //   
          //   
          //   List other_tree_mu_mat_list=prev_sum_trees_mat_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          //   NumericMatrix other_mat_mu=other_tree_mu_mat_list[0];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
          //   NumericMatrix other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
          //   sum_of_trees_mat_mu[count]=other_mat_mu;										// create a add other_mat to the list sum_of_trees_mat
          //   sum_of_trees_mat_tau[count]=other_mat_tau;										// create a add other_mat to the list sum_of_trees_mat
          //   count++;																// increment the count variable.
          //   // Rcout << "LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
          //   // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
          //   
          //   if(count==(size_mat-1)){												// If list size is too small
          //     size_mat=size_mat*2;												// double the size
          //     sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);					// double the length of sum_of_trees
          //     //sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);					// double the length of sum_of_trees
          //     sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,size_mat);		// double the length of sum_of_tree_resids
          //     sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,size_mat);		// double the length of sum_of_tree_resids
          //     sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);			// double the length of sum_of_trees_mat
          //     //sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);			// double the length of sum_of_trees_mat
          //   }
          //   // Rcout << "on  line 5490 LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
          //   // Rcout << "on  line 5491 LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
          //   // Rcout << "loop number" << j << "\n,";
          // }else{
          
          
          
          
            List other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            List other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
            List other_tree_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
            List other_tree_resids_tau=prev_sum_tree_resids_tau[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
            List other_mat_mu=prev_sum_trees_mat_mu[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            //List other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
            for(int f=0;f<other_tree_mu.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
              if(is<NumericMatrix>(other_tree_mu[f])){								// if f+1^th element of other_tree is a matrix, do nothing
              }else{																// if f+1^th element of other_tree is not a matrix
                throw std::range_error("Line 5505 tree is not a numeric matrix!");		// throw an error
              }
              NumericMatrix treetoadd_mu=other_tree_mu[f];								// create matrix treetoadd equal to f+1^th element of other_tree
              if(is<NumericVector>(other_tree_resids_mu[f])){						// if f+1^th element of other_tree_resids is a NumericVector, do nothing
              }else{																// if f+1^th element of other_tree_resids is not a NumericVector
                throw std::range_error("other resids not a numeric matrix!");	// throw an error
              }
              
              // UPDATE PARTIAL RESIDUALS FOR EACH PREVIOUS MU TREE IN THE SUM-OF-TREE MODEL
              NumericVector resids_prevroundtemp_mu=other_tree_resids_mu[f];
              NumericVector residstoadd_mu=resids_prevroundtemp_mu-curr_round_preds(_,k);						// create vector residstoadd equal to f+1^th element of other_tree_resids
              if(is<NumericMatrix>(other_mat_mu[f])){								// if f+1^th element of other_mat is a NumericMatrix, do nothing
                
              }else{																// if f+1^th element of other_mat is not a NumericMatrix
                throw std::range_error(" other mat not a numeric matrix!");		// throw an error
              }
              NumericMatrix mattoadd_mu=other_mat_mu[f];								// create matrix mattoadd equal to f+1^th element of other_mat
              
              sum_of_trees_mu[count]=treetoadd_mu;										// add treetoadd to sum_of_trees
              sum_of_tree_resids_mu[count]=residstoadd_mu;								// add residstoadd to sum_of_tree_resids
              sum_of_trees_mat_mu[count]=mattoadd_mu;									// add mattoadd to sum_of_trees_mat
              count++;															// inremet the count variable.
              
              if(count==(size_mat-1)){											// If list size is too small
                size_mat=size_mat*2;											// double the size
                sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);				// double the length of sum_of_trees
                sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,size_mat);	// double the length of sum_of_tree_resids
                sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);		// double the length of sum_of_trees_mat
              }
            }
            
            if(other_tree_tau.size()>size_mat){											// If list size is too small
              sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,other_tree_tau.size());	// increase the length of sum_of_tree_resids
            }
            for(int f=0;f<other_tree_tau.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
              // UPDATE PARTIAL RESIDUALS FOR EACH TAU TREE IN THE SUM-OF-TREE MODEL
              NumericVector resids_prevroundtemp_tau=other_tree_resids_tau[f];
              NumericVector residstoadd_tau=resids_prevroundtemp_tau-curr_round_preds(_,k);						// create vector residstoadd equal to f+1^th element of other_tree_resids
              
              sum_of_tree_resids_tau[f]=residstoadd_tau;								// add residstoadd to sum_of_tree_resids
            }
            
            
            
            // when a mu tree is added, the list of tau trees is unchanged,
            // therefore the following several lines just keep the same
            // lists of tau trees from the parent sum-of-tree model
            List temp1_prev_sum_tree = prev_sum_trees_tau[curr_round_parent[k]];
            List temp1_prev_sum_tree_mat = prev_sum_trees_tau[curr_round_parent[k]];
            
            sum_of_trees_tau = resize_bcf(sum_of_trees_tau,temp1_prev_sum_tree.size());			// create a list of length 300.
            sum_of_trees_mat_tau = resize_bcf(sum_of_trees_mat_tau,temp1_prev_sum_tree_mat.size());		// create a list of length 300.
            sum_of_trees_tau = prev_sum_trees_tau[curr_round_parent[k]];
            sum_of_trees_mat_tau = prev_sum_trees_mat_tau[curr_round_parent[k]];
              //Rcout << "on  line 8800 LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
              //Rcout << "on  line 8800 LENGTH OF LIST temp1_prev_sum_tree = " << temp1_prev_sum_tree.size() << ".\n";
              //Rcout << "on  line 8800 LENGTH OF LIST temp1_prev_sum_tree_mat = " << temp1_prev_sum_tree_mat.size() << ".\n";
              
          //} REMOVED j==1 block of code. commented out end of block on this line.
          
          
           //Rcout << "on  line 8779 LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
           //Rcout << "on  line 8806 LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
           //Rcout << "loop number" << j << "\n,";
          
          
          
          // ADD TO the MU sum_of lists because mu tree added to for current value of k
          
          sum_of_trees_mu[count]=temp_sum_trees[k];										// add k+1^th element of temp_sum_trees to sum_of_trees
          sum_of_tree_resids_mu[count]=temp_sum_tree_resids[k];							// add k+1^th element of temp_sum_tree_resids to sum_of_tree_resids
          sum_of_trees_mat_mu[count]=temp_sum_trees_mat[k];								// add k+1^th element of temp_sum_trees_mat to sum_of_trees_mat
          count++;																	// increment the count variable
          
          // Rcout << "LENGTH OF LIST SUM_OF_TREES_TAU = " << sum_of_trees_tau.size() << ".\n";
          // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_TAU = " << sum_of_trees_mat_tau.size() << ".\n";
          
          
          if(count==(size_mat-1)){													// If list size is too small
            size_mat=size_mat*2;													// double the size
            sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);						// double the length of sum_of_trees
            if(j==0) sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);						// double the length of sum_of_trees
            
            sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);				// double the length of sum_of_trees_mat_mu
            if(j==0) sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);				// double the length of sum_of_trees_mat_tau
            sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,size_mat);			// double the length of sum_of_trees_mat
          }
        }
         //Rcout << "Get to line 8806 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        sum_of_trees_mu=resize_bcf(sum_of_trees_mu,count);										// remove spaces that are not filled in.
        if(j==0) sum_of_trees_tau=resize_bcf(sum_of_trees_tau,0);
        // remove spaces that are not filled in.
        sum_of_trees_mat_mu=resize_bcf(sum_of_trees_mat_mu,count);								// remove spaces that are not filled in.
        if(j==0) sum_of_trees_mat_tau=resize_bcf(sum_of_trees_mat_tau,0);								// remove spaces that are not filled in.
        sum_of_tree_resids_mu=resize_bcf(sum_of_tree_resids_mu,count);							// remove spaces that are not filled in.
        if(j==0) sum_of_tree_resids_tau=resize_bcf(sum_of_tree_resids_tau,0);
         //Rcout << "Get to line 8813 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        if(j>0){
          List other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
          sum_of_tree_resids_tau=resize_bcf(sum_of_tree_resids_tau,other_tree_tau.size());							// remove spaces that are not filled in.
        }
         //Rcout << "Get to line 8846 in loop j = " << j << " and inner loop k = "<< k <<".\n";
        
        
        // NOW ADD TO overall_ LISTS
        // removed if(curr_round_parent[k]!=-1){
        // so that everything is a list
        
        //Rcout << "Line 8851.\n";
        //Rcout << "length of sum of trees = " << sum_of_trees_mu.size() << ".\n";
        //Rcout << "count = " << count << ".\n";
        //Rcout << "overall_count = " << overall_count << ".\n";
        
        //Rcout << "length of sum of treestau = " << sum_of_trees_tau.size() << ".\n";
        //Rcout << "curr_round_lik.size() = " << curr_round_lik.size() << ".\n";
        //Rcout << "k = " << k << ".\n";
        //Rcout << "j = " << j << ".\n";
        //Rcout << "sum_of_tree_resids_mu.size()  = " << sum_of_tree_resids_mu.size()  << ".\n";
        //Rcout << "sum_of_tree_resids_tau.size()  = " << sum_of_tree_resids_tau.size()  << ".\n";
        
        //Rcout << "Line 8865.\n";
        
        
        
        //if(curr_round_parent[k]!=-1){															// If the k+1^th element of curr_round_parent is -1, (-1 is a terminal node?)
          
          overall_sum_trees_mu[overall_count]=sum_of_trees_mu;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_tau[overall_count]=sum_of_trees_tau;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
          // Rcout <<"ADD ELEMENT TO LIST OF SUM TREE LISTS OF LENGTH"<< sum_of_trees_tau.size() << ".\n";
          overall_sum_tree_resids_mu[overall_count]=sum_of_tree_resids_mu;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
          overall_sum_tree_resids_tau[overall_count]=sum_of_tree_resids_tau;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_mat_mu[overall_count]=sum_of_trees_mat_mu;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
          overall_sum_trees_mat_tau[overall_count]=sum_of_trees_mat_tau;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
          
          //Rcout << "line 8879 overall_sum_trees_mu.size() = " << overall_sum_trees_mu.size() << ".\n";
          //Rcout << "overall_sum_trees_tau.size() = " << overall_sum_trees_tau.size() << ".\n";
          //Rcout << "overall_sum_tree_resids_mu.size() = " << overall_sum_tree_resids_mu.size() << ".\n";
          //Rcout << "overall_sum_tree_resids_tau.size() = " << overall_sum_tree_resids_tau.size() << ".\n";
          //Rcout << "line 8883 overall_count = " << overall_count << ".\n";
          
          
          // THE FOLLOWING LINES (before overall_count++) DON'T NEED TO BE IN THE LOOP OVER k.
          // CAN MOVE THESE LINES OUTSIDE THE LOOP OVER k
          overall_sum_BIC=temp_BIC;															// Let overall_sum_BIC equal temp_BIC, the vector of BICs.
          overall_sum_preds_outcome=Rcpp::as<arma::mat>(temp_preds_outcome);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          overall_sum_preds_mu=Rcpp::as<arma::mat>(temp_preds_mu);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          overall_sum_preds_tau=Rcpp::as<arma::mat>(temp_preds_tau);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
          if(is_test_data==1){	// If 8174 is test data, overall_sum_test_preds equal temp_test_preds, the matrix of out-of-sample predictions (columns correpond to different models?)
             //Rcout << "Get to Line 8843 in loop j = " << j << ". k= "<< k << " .\n";
            overall_sum_test_preds_outcome=Rcpp::as<arma::mat>(temp_test_preds_outcome);
            overall_sum_test_preds_mu=Rcpp::as<arma::mat>(temp_test_preds_mu);
            overall_sum_test_preds_tau=Rcpp::as<arma::mat>(temp_test_preds_tau);
          }
          
          
           //Rcout << " before increment 8850 overall_count=" << overall_count << ". k = "<< k << ".\n";
          overall_count++;																	// increment overall_count
          // Rcout << " after increment 8824 overall_count=" << overall_count << ". k = "<< k << ".\n";
          
          if(overall_count==(overall_size-1)){												// If overall_size is too small
            overall_size=overall_size*2;													// double the size
            overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
            overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
            overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
            overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
            overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
            overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
          }
        //}
        
        
        }else{
          
          //NOW TAU VERSION OF SECOND LOOP OF LENGTH curr_round_lik.size()
          if(curr_round_added_mu_addall[k]==2){
            
            
            
            int size_mat=300;						// create a variable, Initializd equal to 300.
            List sum_of_trees_mu(size_mat);			// create a list of length 300.
            List sum_of_trees_tau(size_mat);			// create a list of length 300.
            List sum_of_tree_resids_mu(size_mat);		// create a list of length 300.
            List sum_of_tree_resids_tau(size_mat);		// create a list of length 300.
            List sum_of_trees_mat_mu(size_mat);		// create a list of length 300.
            List sum_of_trees_mat_tau(size_mat);		// create a list of length 300.
            int count=0;							// create a variable count equal to 0. (Count was already define, so could remove "int" at start of this line and just reset cound to 0).
            
            if(curr_round_parent[k]==-1){			// If the k+1^th element of curr_round_parent is -1, do nothing. (-1 is a terminal node?)
              // Rcout << "line 8855 curr_round_parent[k]==-1 in loop j = " << j << " and inner loop k = "<< k <<".\n";
              sum_of_trees_tau[count]= temp_sum_trees[k];
              sum_of_tree_resids_tau[count]=temp_sum_tree_resids[k];							// add k+1^th element of temp_sum_tree_resids to sum_of_tree_resids
              sum_of_trees_mat_tau[count]=temp_sum_trees_mat[k];
              count=1;
              
            }else{									// If the k+1^th element of curr_round_parent is not equal to -1.
              // NEED TO THINK MORE ABOUT j==1 CASE. Also need the prev_sum_trees to be correct
              
              // if(j==0){
              //   sum_of_trees_mu = resize_bcf(sum_of_trees_mu,1);			// create a list of length 300.
              //   sum_of_trees_mat_mu = resize_bcf(sum_of_trees_mat_mu,1);		// create a list of length 300.
              //   
              //   // Rcout << "LENGTH OF LIST SUM_OF_TREES_MU = " << sum_of_trees_mu.size() << ".\n";
              //   // Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_MU = " << sum_of_trees_mat_mu.size() << ".\n";
              //   
              //   NumericMatrix other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
              //   //NumericMatrix other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];			// create matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
              //   NumericVector other_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];	// create vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
              //   // Rcout << "LENGTH OF LIST other_resids_mu = " << other_resids_mu.size() << ".\n";
              //   sum_of_trees_mu[count]= other_tree_mu;										// add other_tree to the list sum_of_trees
              //   //sum_of_trees_tau[count]= other_tree_tau;										// add other_tree to the list sum_of_trees
              //   sum_of_tree_resids_mu[count]=other_resids_mu-z*curr_round_preds(_,k);									// add other_resids to the list sum_of_tree_resids
              //   NumericMatrix other_mat_mu=prev_sum_trees_mat_mu[curr_round_parent[k]];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
              //   //NumericMatrix other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];		// create a matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
              //   sum_of_trees_mat_mu[count]=other_mat_mu;										// create a add other_mat to the list sum_of_trees_mat
              //   //sum_of_trees_mat_tau[count]=other_mat_tau;										// create a add other_mat to the list sum_of_trees_mat
              //   //count++;																// increment the count variable.
              //   
              //   //if(count==(size_mat-1)){												// If list size is too small
              //   //  size_mat=size_mat*2;												// double the size
              //   //  sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);					// double the length of sum_of_trees
              //   //  sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);					// double the length of sum_of_trees
              //   //  sum_of_tree_resids=resize_bigger_bcf(sum_of_tree_resids,size_mat);		// double the length of sum_of_tree_resids
              //   //  sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);			// double the length of sum_of_trees_mat
              //   //  sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);			// double the length of sum_of_trees_mat
              //   //}
              //   
              // }else{
                List other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
                List other_tree_tau=prev_sum_trees_tau[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
                List other_tree_resids_mu=prev_sum_tree_resids_mu[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
                List other_tree_resids_tau=prev_sum_tree_resids_tau[curr_round_parent[k]];		// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? vector equal to the curr_round_parent[k]+1^th element of the list prev_sum_tree_resids.
                //List other_mat_mu=prev_sum_trees_mat_mu[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
                List other_mat_tau=prev_sum_trees_mat_tau[curr_round_parent[k]];				// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the matrix prev_sum_trees_mat.
                for(int f=0;f<other_tree_tau.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
                  if(is<NumericMatrix>(other_tree_tau[f])){								// if f+1^th element of other_tree is a matrix, do nothing
                  }else{																// if f+1^th element of other_tree is not a matrix
                    throw std::range_error("Line 6211 tree is not a numeric matrix!");		// throw an error
                  }
                  NumericMatrix treetoadd_tau=other_tree_tau[f];								// create matrix treetoadd equal to f+1^th element of other_tree
                  if(is<NumericVector>(other_tree_resids_tau[f])){						// if f+1^th element of other_tree_resids is a NumericVector, do nothing
                  }else{																// if f+1^th element of other_tree_resids is not a NumericVector
                    throw std::range_error("other resids not a numeric matrix!");	// throw an error
                  }
                  NumericVector resids_prevroundtemp_tau=other_tree_resids_tau[f];						// create vector residstoadd equal to f+1^th element of other_tree_resids
                  NumericVector residstoadd_tau=resids_prevroundtemp_tau-z*curr_round_preds(_,k);
                  
                  if(is<NumericMatrix>(other_mat_tau[f])){								// if f+1^th element of other_mat is a NumericMatrix, do nothing
                    
                  }else{																// if f+1^th element of other_mat is not a NumericMatrix
                    throw std::range_error(" other mat not a numeric matrix!");		// throw an error
                  }
                  NumericMatrix mattoadd_tau=other_mat_tau[f];								// create matrix mattoadd equal to f+1^th element of other_mat
                  
                  sum_of_trees_tau[count]=treetoadd_tau;										// add treetoadd to sum_of_trees
                  sum_of_tree_resids_tau[count]=residstoadd_tau;								// add residstoadd to sum_of_tree_resids
                  sum_of_trees_mat_tau[count]=mattoadd_tau;									// add mattoadd to sum_of_trees_mat
                  count++;															// inremet the count variable.
                  
                  if(count==(size_mat-1)){											// If list size is too small
                    size_mat=size_mat*2;											// double the size
                    sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);				// double the length of sum_of_trees
                    sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,size_mat);	// double the length of sum_of_tree_resids
                    sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);		// double the length of sum_of_trees_mat
                  }
                }
                
                if(other_tree_mu.size()>size_mat){											// If list size is too small
                  sum_of_tree_resids_mu=resize_bigger_bcf(sum_of_tree_resids_mu,other_tree_mu.size());	// increase the length of sum_of_tree_resids
                }
                for(int f=0;f<other_tree_mu.size();f++){									// for-loop of length equal to that of other_tree?? length is one if list??
                  NumericVector resids_prevroundtemp_mu=other_tree_resids_mu[f];
                  NumericVector residstoadd_mu=resids_prevroundtemp_mu-z*curr_round_preds(_,k);						// create vector residstoadd equal to f+1^th element of other_tree_resids
                  
                  sum_of_tree_resids_mu[f]=residstoadd_mu;								// add residstoadd to sum_of_tree_resids
                }
                
                
                
                // when a tau tree is added, the list of mu trees is unchanged,
                // therefore the following several lines just keep the same
                // lists of mu trees from the parent sum-of-tree model
                List temp1_prev_sum_tree = prev_sum_trees_mu[curr_round_parent[k]];
                List temp1_prev_sum_tree_mat = prev_sum_trees_mat_mu[curr_round_parent[k]];
                
                sum_of_trees_mu = resize_bcf(sum_of_trees_mu,temp1_prev_sum_tree.size());			
                sum_of_trees_mat_mu = resize_bcf(sum_of_trees_mat_mu,temp1_prev_sum_tree_mat.size());		
                
                sum_of_trees_mu = prev_sum_trees_mu[curr_round_parent[k]];
                sum_of_trees_mat_mu = prev_sum_trees_mat_mu[curr_round_parent[k]];
                //Rcout << "on  line 9024 LENGTH OF LIST SUM_OF_TREES_MU = " << sum_of_trees_mu.size() << ".\n";
                //Rcout << "on  line 9025 LENGTH OF LIST (a list of mu trees) temp1_prev_sum_tree = " << temp1_prev_sum_tree.size() << ".\n";
                //Rcout << "on  line 9026 LENGTH OF LIST (a list of mu trees) temp1_prev_sum_tree_mat = " << temp1_prev_sum_tree_mat.size() << ".\n";
                
              
              
              //} Removed j==0 block of code, this line commented out because this was end of the else-block
              
              
                //Rcout << "GET TO LINE 8984 .\n";
              
              
              // ADD TO the TAU sum_of lists because tau tree added to for current value of k
              
              sum_of_trees_tau[count]=temp_sum_trees[k];										// add k+1^th element of temp_sum_trees to sum_of_trees
              sum_of_tree_resids_tau[count]=temp_sum_tree_resids[k];							// add k+1^th element of temp_sum_tree_resids to sum_of_tree_resids
              sum_of_trees_mat_tau[count]=temp_sum_trees_mat[k];								// add k+1^th element of temp_sum_trees_mat to sum_of_trees_mat
              count++;																	// increment the count variable
              
                //Rcout << "GET TO LINE 9038 .\n";
               //Rcout << "LENGTH OF LIST SUM_OF_TREES_MU = " << sum_of_trees_mu.size() << ".\n";
               //Rcout << "LENGTH OF LIST SUM_OF_TREES_MAT_MU = " << sum_of_trees_mat_mu.size() << ".\n";
              
              if(count==(size_mat-1)){													// If list size is too small
                size_mat=size_mat*2;													// double the size
                if(j==0) sum_of_trees_mu=resize_bigger_bcf(sum_of_trees_mu,size_mat);						// double the length of sum_of_trees
                sum_of_trees_tau=resize_bigger_bcf(sum_of_trees_tau,size_mat);						// double the length of sum_of_trees
                
                if(j==0) sum_of_trees_mat_mu=resize_bigger_bcf(sum_of_trees_mat_mu,size_mat);				// double the length of sum_of_tree_resids
                sum_of_trees_mat_tau=resize_bigger_bcf(sum_of_trees_mat_tau,size_mat);				// double the length of sum_of_tree_resids
                sum_of_tree_resids_tau=resize_bigger_bcf(sum_of_tree_resids_tau,size_mat);			// double the length of sum_of_trees_mat
              }
            }
             //Rcout << "count=" << count << ".\n";
            
            if(j==0) sum_of_trees_mu=resize_bcf(sum_of_trees_mu,0);										// remove spaces that are not filled in.
            sum_of_trees_tau=resize_bcf(sum_of_trees_tau,count);										// remove spaces that are not filled in.
            if(j==0) sum_of_trees_mat_mu=resize_bcf(sum_of_trees_mat_mu,0);								// remove spaces that are not filled in.
            sum_of_trees_mat_tau=resize_bcf(sum_of_trees_mat_tau,count);								// remove spaces that are not filled in.
            sum_of_tree_resids_tau=resize_bcf(sum_of_tree_resids_tau,count);							// remove spaces that are not filled in.
            if(j==0) sum_of_tree_resids_mu = resize_bcf(sum_of_tree_resids_mu,0);
            
            if(j>0){
              List other_tree_mu=prev_sum_trees_mu[curr_round_parent[k]];					// create List?? (Maybe curr_round_parent[k] is a vector/list of indices)?? matrix equal to the curr_round_parent[k]+1^th element of the list prev_sum_trees.
              sum_of_tree_resids_mu=resize_bcf(sum_of_tree_resids_mu,other_tree_mu.size());							// remove spaces that are not filled in.
            }
            
            //Rcout << "Line 9061.\n";
            //Rcout << "length of sum of trees = " << sum_of_trees_mu.size() << ".\n";
            //Rcout << "count = " << count << ".\n";
            //Rcout << "overall_count = " << overall_count << ".\n";
            //Rcout << "length of sum of treestau = " << sum_of_trees_tau.size() << ".\n";
            //Rcout << "curr_round_lik.size() = " << curr_round_lik.size() << ".\n";
            //Rcout << "k = " << k << ".\n";
            //Rcout << "j = " << j << ".\n";
            //Rcout << "sum_of_tree_resids_mu.size()  = " << sum_of_tree_resids_mu.size()  << ".\n";
            //Rcout << "sum_of_tree_resids_tau.size()  = " << sum_of_tree_resids_tau.size()  << ".\n";
            
            //Rcout << "Line 9071.\n";
            
      
      // NOW ADD TO overall_ LISTS
      // removed if(curr_round_parent[k]!=-1){
      // so that everything is a list
            
      //      if(curr_round_parent[k]!=-1){															// If the k+1^th element of curr_round_parent is -1, (-1 is a terminal node?)
              overall_sum_trees_mu[overall_count]=sum_of_trees_mu;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
              overall_sum_trees_tau[overall_count]=sum_of_trees_tau;								// overall_count+1^th element of overall_sum_trees is sum_of_trees (which is itself a list... therefore have a list of lists?)
              overall_sum_tree_resids_mu[overall_count]=sum_of_tree_resids_mu;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
              overall_sum_tree_resids_tau[overall_count]=sum_of_tree_resids_tau;						// overall_count+1^th element of overall_sum_tree_resids is sum_of_tree_resids (which is itself a list... therefore have a list of lists?)
              overall_sum_trees_mat_mu[overall_count]=sum_of_trees_mat_mu;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
              overall_sum_trees_mat_tau[overall_count]=sum_of_trees_mat_tau;						// overall_count+1^th element of overall_sum_trees_mat is sum_of_trees_mat (which is itself a list... therefore have a list of lists?)
              
              
              //Rcout << "line 9097 overall_sum_trees_mu.size() = " << overall_sum_trees_mu.size() << ".\n";
              //Rcout << "overall_sum_trees_tau.size() = " << overall_sum_trees_tau.size() << ".\n";
              //Rcout << "overall_sum_tree_resids_mu.size() = " << overall_sum_tree_resids_mu.size() << ".\n";
              //Rcout << "overall_sum_tree_resids_tau.size() = " << overall_sum_tree_resids_tau.size() << ".\n";
              //Rcout << "line 9101 overall_count = " << overall_count << ".\n";
              
              
              // THE FOLLOWING LINES (before overall_count++) DON'T NEED TO BE IN THE LOOP OVER k.
              // CAN MOVE THESE LINES OUTSIDE THE LOOP OVER k
              overall_sum_BIC=temp_BIC;															// Let overall_sum_BIC equal temp_BIC, the vector of BICs.
              overall_sum_preds_outcome=Rcpp::as<arma::mat>(temp_preds_outcome);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
              overall_sum_preds_mu=Rcpp::as<arma::mat>(temp_preds_mu);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
              overall_sum_preds_tau=Rcpp::as<arma::mat>(temp_preds_tau);									// Let overall_sum_preds equal temp_let preds, the matrix of predictions (columns correspond to different modes?).
              if(is_test_data==1){	// If there is test data, overall_sum_test_preds equal temp_test_preds, the matrix of out-of-sample predictions (columns correpond to different models?)
                //Rcout << "Get to Line 6327 in loop j = " << j  << ".\n";
                overall_sum_test_preds_outcome=Rcpp::as<arma::mat>(temp_test_preds_outcome);
                overall_sum_test_preds_mu=Rcpp::as<arma::mat>(temp_test_preds_mu);
                overall_sum_test_preds_tau=Rcpp::as<arma::mat>(temp_test_preds_tau);
              }
              
               //Rcout << " before increment 9052 overall_count=" << overall_count << ". k = "<< k << ".\n";
              
              overall_count++;																	// increment overall_count
              // Rcout << " after increment 9027 overall_count=" << overall_count << ". k = " << k << ".\n";
              
              if(overall_count==(overall_size-1)){												// If overall_size is too small
                overall_size=overall_size*2;													// double the size
                overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
              }
        //    }  
            
            
            
            
            
          }else{
            throw std::range_error("Added tree not recorded as mu or tau tree");
            
          }
        }
      
      }// End of second loop of length curr_round_lik.size()
      // The above loop creates the sum_of lists      
      
       //Rcout << "overall_sum_trees_mu.size() = " << overall_sum_trees_mu.size() << ".\n";
       //Rcout << "overall_count = " << overall_count << ".\n";
      
      
      
      
      

      
      
      
      
      
      
      // THE FOLLOWING CODE (optional) KEEPS SUM-OF-TREE MODELS FROM PREVIOUS ROUNDS
      // THAT DID NOT HAVE DAUGHTER TREES GROWN
      // i.e. SUM OF TREE MODELS WITH LESS THAN j(+1) TREES.
      
       //Rcout << "Get to line 9099 in loop j = " << j << ".\n";
      if(only_max_num_trees==0){
        
        //check if there were any trees from the previous round that didn't have daughter trees grown.
        //create vector to count number of possible parents for previous round
        
        //Note that this code is only spplied after the first round
        if(j>0){																					// If not the first round of the outter loop
          IntegerVector prev_par_no_child=match(prev_par,curr_round_parent);						// create a vector equal to indices of the positions of the (the first matches of the) elements of prev_par in curr_round_parent.
          if(any(is_na(prev_par_no_child))){														// any of the vector of matches are NA (no match?)
            IntegerVector t4=ifelse(is_na(prev_par_no_child),1,0);								// create a vector t4 equal to 1 for the NA values, 0 otherwise.
            for(int h=0;h<prev_par_no_child.size();h++){										// for-loop of length equal to that of prev_par_no_child
              if(t4[h]==1){																	// If h+1^th element of vector of matches is NA
                if(prev_round_BIC2[h]-lowest_BIC<=log(c)){									// If the h+1^th model (from the previous round?) is in Occam's window
                  
                  // EVERYTHING SHOULD BE A LIST. POSSIBLY LIST OF SIZE ZERO.
                  // CODE BELOW WITH IF STATEMENTS IS FOR TESTING.
                  // REMOVE AFTER TESTING
                  
                  SEXP s_mu = prev_sum_trees_mu[h];												// create a pointer to S expression type equal to the h+1^th element of prev_sum_trees (a tree table or list of tree tables from the previous round?)
                  SEXP s_tau = prev_sum_trees_tau[h];												// create a pointer to S expression type equal to the h+1^th element of prev_sum_trees (a tree table or list of tree tables from the previous round?)
                  SEXP s_resid_mu = prev_sum_tree_resids_mu[h]; 
                  SEXP s_resid_tau = prev_sum_tree_resids_tau[h];
                  
                  if(is<List>(s_resid_mu)){
                    if(is<List>(s_resid_tau)){
                      
                      //Rcout << "mu resid list and tau resid list. j= " << j << ".\n";
                      
                    }else{
                      // Rcout << "mu resid list and tau resid vector?. j= " << j << ".\n";
                      throw std::range_error("mu resid list and tau resid vector?");
                      
                    }
                  }else{
                    if(is<List>(s_resid_tau)){
                      // Rcout << "mu resid vector and tau resid list. j= " << j << ".\n";
                      throw std::range_error("mu resid vector and tau resid list");
                      
                    }else{
                      // Rcout << "mu resid vector and tau resid vector. j= " << j << ".\n";
                      throw std::range_error("mu resid vector and tau resid vector");
                      
                    }
                  }
                  
                  if(is<List>(s_mu)){														// If prev_sum_trees[h] is a list
                    
                    if(is<List>(s_tau)){
                      List tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                      List tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                      List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                      //if(is<List>(s_resid_tau)){
                      List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      //}else{
                      //  NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      //  overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      //}
                      List treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                      List treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                      overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                      overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                      overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                      overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                      overall_count++;													// increment overall_count
                      
                      if(overall_count==(overall_size-1)){												// If overall_size is too small
                        overall_size=overall_size*2;													// double the size
                        overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                        overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                        overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                        overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                        overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                        overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                      }
                      double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                      overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                      overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      if(is_test_data==1){
                        //Rcout << "Get to Line 5673 in loop j = " << j  << ".\n";
                        overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                        overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                        overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      }
                    }else{
                      throw std::range_error("prev_sum_trees_tau NOT A LIST");
                      
                      //EVERYTHING SHOULD BE A LIST. 
                      //REMOVE THE COMMENTED OUT CODE BELOW AFTER TESTING
                      
                      // List tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                      // NumericMatrix tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                      // List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                      // if(is<List>(s_resid_tau)){
                      //   List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      //   overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // }else{
                      //   NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      //   overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // }                    
                      // List treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                      // NumericMatrix treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                      // overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                      // overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                      // overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                      // overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                      // overall_count++;													// increment overall_count
                      // 
                      // if(overall_count==(overall_size-1)){												// If overall_size is too small
                      //   overall_size=overall_size*2;													// double the size
                      //   overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      //   overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      //   overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      //   overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      //   overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      //   overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                      // }
                      // double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                      // overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                      // overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // if(is_test_data==1){
                      //   //Rcout << "Get to Line 5708 in loop j = " << j  << ".\n";
                      //   overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      //   overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      //   overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      // }
                      
                      
                      
                    }
                    
                    
                  }else{																// If prev_sum_trees[h] is NOT a list
                    if(is<List>(s_tau)){
                      
                      throw std::range_error("prev_sum_trees_mu NOT A LIST");
                      
                      //EVERYTHING SHOULD BE A LIST. 
                      //REMOVE THE COMMENTED OUT CODE BELOW AFTER TESTING
                      
                      // NumericMatrix tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                      // List tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                      // List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                      // if(is<List>(s_resid_tau)){
                      //   List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      //   overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // }else{
                      //   NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      //   overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // }                    
                      // NumericMatrix treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                      // List treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                      // overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                      // overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                      // overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                      // overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                      // overall_count++;													// increment overall_count
                      // 
                      // if(overall_count==(overall_size-1)){												// If overall_size is too small
                      //   overall_size=overall_size*2;													// double the size
                      //   overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      //   overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      //   overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      //   overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      //   overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      //   overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                      // }
                      // double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                      // overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                      // overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // if(is_test_data==1){
                      //   //Rcout << "Get to Line 5749 in loop j = " << j  << ".\n";
                      //   overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      //   overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      //   overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      // }
                      
                      
                      
                    }else{
                      
                      throw std::range_error("prev_sum_trees_mu and prev_sum_trees_tau NOT LISTS");
                      
                      //EVERYTHING SHOULD BE A LIST. 
                      //REMOVE THE COMMENTED OUT CODE BELOW AFTER TESTING
                      
                      
                      // NumericMatrix tree_no_child_mu=prev_sum_trees_mu[h];								// create a list equal to prev_sum_trees[h]
                      // NumericMatrix tree_no_child_tau=prev_sum_trees_tau[h];								// create a list equal to prev_sum_trees[h]
                      // List resids_no_child_mu=prev_sum_tree_resids_mu[h];						// create a list equal to prev_sum_tree_resids[h]
                      // if(is<List>(s_resid_tau)){
                      //   List resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      //   overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // }else{
                      //   NumericVector resids_no_child_tau=prev_sum_tree_resids_tau[h];						// create a list equal to prev_sum_tree_resids[h]
                      //   overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // }                    
                      // NumericMatrix treemat_no_child_mu=prev_sum_trees_mat_mu[h];						// create a list equal to prev_sum_trees_mat[h]
                      // NumericMatrix treemat_no_child_tau=prev_sum_trees_mat_tau[h];						// create a list equal to prev_sum_trees_mat[h]
                      // overall_sum_trees_mu[overall_count]=tree_no_child_mu;						// add prev_sum_trees[h] to overall_sum_trees
                      // overall_sum_trees_tau[overall_count]=tree_no_child_tau;						// add prev_sum_trees[h] to overall_sum_trees
                      // overall_sum_tree_resids_mu[overall_count]=resids_no_child_mu;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // //overall_sum_tree_resids_tau[overall_count]=resids_no_child_tau;				// add prev_sum_tree_resids[h] to overall_sum_tree_resids
                      // overall_sum_trees_mat_mu[overall_count]=treemat_no_child_mu;				// add  to overall_sum_trees_mat
                      // overall_sum_trees_mat_tau[overall_count]=treemat_no_child_tau;				// add  to overall_sum_trees_mat
                      // overall_count++;													// increment overall_count
                      // 
                      // if(overall_count==(overall_size-1)){												// If overall_size is too small
                      //   overall_size=overall_size*2;													// double the size
                      //   overall_sum_trees_mu=resize_bigger_bcf(overall_sum_trees_mu,overall_size);				// double the length of overall_sum_trees
                      //   overall_sum_trees_tau=resize_bigger_bcf(overall_sum_trees_tau,overall_size);				// double the length of overall_sum_trees
                      //   overall_sum_tree_resids_mu=resize_bigger_bcf(overall_sum_tree_resids_mu,overall_size);	// double the length of overall_sum_tree_resids
                      //   overall_sum_tree_resids_tau=resize_bigger_bcf(overall_sum_tree_resids_tau,overall_size);	// double the length of overall_sum_tree_resids
                      //   overall_sum_trees_mat_mu=resize_bigger_bcf(overall_sum_trees_mat_mu,overall_size);		// double the length of overall_sum_trees_mat
                      //   overall_sum_trees_mat_tau=resize_bigger_bcf(overall_sum_trees_mat_tau,overall_size);		// double the length of overall_sum_trees_mat
                      // }
                      // double BIC_to_add=prev_round_BIC2[h];												// let BIC_to_add equal the h+1^th BIC
                      // overall_sum_BIC.push_back(BIC_to_add);												// append the the h+1^th BIC to overall_sum_BIC
                      // overall_sum_preds_outcome.insert_cols(overall_sum_preds_outcome.n_cols,prev_round_preds2_outcome.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // overall_sum_preds_mu.insert_cols(overall_sum_preds_mu.n_cols,prev_round_preds2_mu.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // overall_sum_preds_tau.insert_cols(overall_sum_preds_tau.n_cols,prev_round_preds2_tau.col(h));	// add the predictions of the h+1^th model to overall_sum_preds as the last (rightmost) column
                      // if(is_test_data==1){
                      //   //Rcout << "Get to Line 5784 in loop j = " << j  << ".\n";
                      //   overall_sum_test_preds_outcome.insert_cols(overall_sum_test_preds_outcome.n_cols,prev_round_test_preds2_outcome.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      //   overall_sum_test_preds_mu.insert_cols(overall_sum_test_preds_mu.n_cols,prev_round_test_preds2_mu.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      //   overall_sum_test_preds_tau.insert_cols(overall_sum_test_preds_tau.n_cols,prev_round_test_preds2_tau.col(h));	// If there is test data, add the out-of-sample predictions of the h+1^th model to overall_sum_test_preds as the last (rightmost) column.
                      // }
                      
                      
                      
                    }//closes else (tau list or not)
                    
                  }//closes else (mu list or not)
                } // closes if(prev_round_BIC2[h]-lowest_BIC<=log(c))
              }// closes if(t4[h]==1){
            } // closes for(int h=0;h<prev_par_no_child.size();h++)
          }// closes if(any(is_na(prev_par_no_child)))
        } // closes if(j>0){
      } // closes if(only_max_num_trees==0)
      
      
      
       //Rcout << "Get to Line 9354 in loop j = " << j  << ".\n";
      
      
      prev_round_preds_outcome=temp_preds_outcome;															// let prev_round_preds equal to the matrix of predictions.
      prev_round_preds_mu=temp_preds_mu;															// let prev_round_preds equal to the matrix of predictions.
      prev_round_preds_tau=temp_preds_tau;															// let prev_round_preds equal to the matrix of predictions.
      if(is_test_data==1){
        //Rcout << "Get to Line 5802 in loop j = " << j  << ".\n";
        prev_round_test_preds_outcome=temp_test_preds_outcome;								// if there is test data, let prev_round_test_preds equal the test data predictions
        prev_round_test_preds_mu=temp_test_preds_mu;								// if there is test data, let prev_round_test_preds equal the test data predictions
        prev_round_test_preds_tau=temp_test_preds_tau;								// if there is test data, let prev_round_test_preds equal the test data predictions
      }
      prev_round_BIC=temp_BIC;																// let prev_round_BIC equal the vector of BICs
      prev_round_BIC2=temp_BIC;																// let prev_round_BIC2 equal the vector of BICs
      prev_round_preds2_outcome=Rcpp::as<arma::mat>(temp_preds_outcome);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      prev_round_preds2_mu=Rcpp::as<arma::mat>(temp_preds_mu);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      prev_round_preds2_tau=Rcpp::as<arma::mat>(temp_preds_tau);										// let prev_round_preds equal arma mat copy of the matrix of predictions.
      if(is_test_data==1){
        //Rcout << "Get to Line 5813 in loop j = " << j  << ".\n";
        prev_round_test_preds2_outcome=Rcpp::as<arma::mat>(temp_test_preds_outcome);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
        prev_round_test_preds2_mu=Rcpp::as<arma::mat>(temp_test_preds_mu);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
        prev_round_test_preds2_tau=Rcpp::as<arma::mat>(temp_test_preds_tau);		// if there is test data, let prev_round_test_preds2 be an arma mat copy of the test data predictions
      }
      resids=temp_resids;																		// let resids equal the matrix of residuals
      parent=temp_parent;																		// let parent equal the parent vector
      overall_sum_trees_mu=resize_bcf(overall_sum_trees_mu,overall_count);								// remove spaces that are not filled in.
      overall_sum_trees_tau=resize_bcf(overall_sum_trees_tau,overall_count);								// remove spaces that are not filled in.
      overall_sum_tree_resids_mu=resize_bcf(overall_sum_tree_resids_mu,overall_count);					// remove spaces that are not filled in.
      overall_sum_tree_resids_tau=resize_bcf(overall_sum_tree_resids_tau,overall_count);					// remove spaces that are not filled in.
      overall_sum_trees_mat_mu=resize_bcf(overall_sum_trees_mat_mu,overall_count);						// remove spaces that are not filled in.
      overall_sum_trees_mat_tau=resize_bcf(overall_sum_trees_mat_tau,overall_count);						// remove spaces that are not filled in.
      
       //Rcout << "overall_count = " << overall_count << ".\n";
      // Rcout << "AFTER RESIZE overall_sum_trees_mu.size() = " << overall_sum_trees_mu.size() << ".\n";
      // Rcout << "AFTER RESIZE overall_sum_trees_tau.size() = " << overall_sum_trees_tau.size() << ".\n";
       //Rcout << "curr_round_lik.size() =" << curr_round_lik.size() << ".\n";
      
      
      
      // commented out if(first_round==1) because everything should be a list regardless of the round
      
      // if(first_round==1){																		// if in the first round of the outer for-loop (j==0)
      //   prev_sum_trees_mu=temp_sum_trees_mu;														// let prev_sum_trees equal the list of tree tables (from the current round)
      //   //ADDING TO MU TREES, therefore nothing yet added to tau trees
      //   List prev_sum_trees_tau(prev_sum_trees_mu.size());
      //   prev_sum_trees_mat_mu=temp_sum_trees_mat_mu;												// let prev_sum_trees_mat equal the list of tree matrice
      //   List prev_sum_trees_mat_tau(prev_sum_trees_mu.size());
      //   
      //   for(int p=0;p<curr_round_lik.size();p++){
      //     prev_sum_trees_tau[p] = start_tree_bcf(0,0); // maybe should be empty list rather than list of empty trees?
      //     prev_sum_trees_mat_tau[p] = start_matrix_bcf(n); // maybe should be empty list? .. model shouldn't even have stub tree
      //   }
      //   prev_sum_tree_resids_mu=temp_sum_tree_resids_mu;											// let prev_sum_tree_resids equal the list of residual vectors (from the current round)
      //   //prev_sum_tree_resids_tau=;											// let prev_sum_tree_resids equal the list of residual vectors (from the current round)
      //   
      //   //NumericMatrix test=prev_sum_trees[0];												// create a matrix test equal to the first element of temp_sum_trees (first obtained tree table)
      //   overall_sum_trees_mu=resize_bcf(temp_sum_trees_mu,temp_sum_trees_mu.size());						// remove spaces that are not filled in.
      //   overall_sum_trees_tau=resize_bcf(prev_sum_trees_tau,temp_sum_trees_mu.size());// INTENTIONALLY USING SIZE OF MU						// remove spaces that are not filled in.
      //   overall_sum_trees_mu=temp_sum_trees_mu;													// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
      //   overall_sum_trees_tau=prev_sum_trees_tau; // not sure about this													// let overall_sum_trees equal the RESIZED list of tree tables (from the current round)
      //   overall_sum_tree_resids_mu=resize_bcf(temp_sum_tree_resids_mu,temp_sum_tree_resids_mu.size());	// remove spaces that are not filled in.
      //   overall_sum_tree_resids_mu=temp_sum_tree_resids_mu;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)
      //   overall_sum_tree_resids_tau=resize_bcf(prev_sum_tree_resids_tau,prev_sum_tree_resids_tau.size());	// remove spaces that are not filled in.
      //   overall_sum_tree_resids_tau=prev_sum_tree_resids_tau;										// let overall_sum_tree_resids equal the RESIZED list of tree residual vectors (from the current round)
      //   overall_sum_trees_mat_mu=temp_sum_trees_mat_mu;											// let overall_sum_trees_mat equal the RESIZED list of tree matrice
      //   overall_sum_trees_mat_tau=prev_sum_trees_mat_tau;											// let overall_sum_trees_mat equal the RESIZED list of tree matrice
      //   overall_sum_trees_mat_mu=resize_bcf(temp_sum_trees_mat_mu,temp_sum_trees_mat_mu.size());				// remove spaces that are not filled in.
      //   overall_sum_trees_mat_tau=resize_bcf(prev_sum_trees_mat_tau,temp_sum_trees_mat_mu.size());				// remove spaces that are not filled in.
      //   overall_sum_BIC=temp_BIC;															// let overall_sum_BIC equal the vector of BICs
      //   overall_sum_preds_outcome= Rcpp::as<arma::mat>(temp_preds_outcome);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
      //   overall_sum_preds_mu= Rcpp::as<arma::mat>(temp_preds_mu);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
      //   overall_sum_preds_tau= Rcpp::as<arma::mat>(temp_preds_tau);									// let overall_sum_preds equal arma mat copy of the matrix of predictions.
      //   if(is_test_data==1){ // if there is test data, let overall_sum_test_preds be an arma mat copy of the test data predictions
      //     //Rcout << "Get to Line 5860 in loop j = " << j  << ".\n";
      //     overall_sum_test_preds_outcome= Rcpp::as<arma::mat>(temp_test_preds_outcome);
      //     overall_sum_test_preds_mu= Rcpp::as<arma::mat>(temp_test_preds_mu);
      //     overall_sum_test_preds_tau= Rcpp::as<arma::mat>(temp_test_preds_tau);
      //   }
      // }else{																					// if not in the first round of the outer for-loop. (i.e. j>0)
        
        prev_sum_trees_mu=overall_sum_trees_mu;													// let prev_sum_trees equal the list of tree tables (up to the current round??)
        prev_sum_trees_tau=overall_sum_trees_tau;													// let prev_sum_trees equal the list of tree tables (up to the current round??)
        prev_sum_tree_resids_mu=overall_sum_tree_resids_mu;										// let prev_sum_tree_resids equal the list of residual vectors (up tp the current round)
        prev_sum_tree_resids_tau=overall_sum_tree_resids_tau;										// let prev_sum_tree_resids equal the list of residual vectors (up tp the current round)
        prev_sum_trees_mat_mu=overall_sum_trees_mat_mu;											// let prev_sum_trees_mat equal the list of tree matrice (up to the current round??)
        prev_sum_trees_mat_tau=overall_sum_trees_mat_tau;											// let prev_sum_trees_mat equal the list of tree matrice (up to the current round??)
        
        
        prev_round_BIC2=overall_sum_BIC;													// let prev_round_BIC2 equal the vector of BICs (up to the current round??)
        prev_round_preds2_outcome=overall_sum_preds_outcome;												// let prev_round_preds2 equal the predictions matrix
        prev_round_preds2_mu=overall_sum_preds_mu;												// let prev_round_preds2 equal the predictions matrix
        prev_round_preds2_tau=overall_sum_preds_tau;												// let prev_round_preds2 equal the predictions matrix
        if(is_test_data==1){
          //Rcout << "Get to Line 5880 in loop j = " << j  << ".\n";
          prev_round_test_preds2_outcome=overall_sum_test_preds_outcome;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
          prev_round_test_preds2_mu=overall_sum_test_preds_mu;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
          prev_round_test_preds2_tau=overall_sum_test_preds_tau;					// if there is test data, let prev_round_test_preds2 equal the out-of-sample predictions matrix.
        }
      //}  commented out if(first_round==1) because everything should be a list regardless of the round
      
      
      
      
      overall_overall_sum_trees_mu[oo_count]=overall_sum_trees_mu;									// Add the current outer loop's table list to overall_overall_sum_trees (list of lists?? OR list of lists of lists??)
      overall_overall_sum_trees_tau[oo_count]=overall_sum_trees_tau;									// Add the current outer loop's table list to overall_overall_sum_trees (list of lists?? OR list of lists of lists??)
      overall_overall_sum_tree_resids_mu[oo_count]=overall_sum_tree_resids_mu;						// Add the current outer loop's list of residual vectors to the list overall_overall_sum_tree_resids (list of lists)
      overall_overall_sum_tree_resids_tau[oo_count]=overall_sum_tree_resids_tau;						// Add the current outer loop's list of residual vectors to the list overall_overall_sum_tree_resids (list of lists)
      overall_overall_sum_trees_mat_mu[oo_count]=overall_sum_trees_mat_mu;							// Add the current outer loop's list of matrices to the list overall_overall_sum_trees_mat (list of lists)
      overall_overall_sum_trees_mat_tau[oo_count]=overall_sum_trees_mat_tau;							// Add the current outer loop's list of matrices to the list overall_overall_sum_trees_mat (list of lists)
      overall_overall_sum_BIC[oo_count]=overall_sum_BIC;										// Add the vector of BICs to the list overall_overall_sum_BIC. (list of vectors?)
      oo_count ++;																					// increment the count 
      if(oo_count==(oo_size-1)){																		// If lists are not large enough
        oo_size=oo_size*2;																			// double the size.
        overall_overall_sum_trees_mu=resize_bigger_bcf(overall_overall_sum_trees_mu,oo_size);					// double the length of overall_overall_sum_trees
        overall_overall_sum_trees_tau=resize_bigger_bcf(overall_overall_sum_trees_tau,oo_size);					// double the length of overall_overall_sum_trees
        overall_overall_sum_tree_resids_mu=resize_bigger_bcf(overall_overall_sum_tree_resids_mu,oo_size);		// double the length of overall_overall_sum_tree_resids
        overall_overall_sum_tree_resids_tau=resize_bigger_bcf(overall_overall_sum_tree_resids_tau,oo_size);		// double the length of overall_overall_sum_tree_resids
        overall_overall_sum_trees_mat_mu=resize_bigger_bcf(overall_overall_sum_trees_mat_mu,oo_size);			// double the length of overall_overall_sum_trees_mat
        overall_overall_sum_trees_mat_tau=resize_bigger_bcf(overall_overall_sum_trees_mat_tau,oo_size);			// double the length of overall_overall_sum_trees_mat
        overall_overall_sum_BIC=resize_bigger_bcf(overall_overall_sum_BIC,oo_size);						// double the length of overall_overall_sum_BIC
      }    
      overall_overall_sum_preds_outcome=overall_sum_preds_outcome;											// let overall_overall_sum_preds equal the prediction matrix.
      overall_overall_sum_preds_mu=overall_sum_preds_mu;											// let overall_overall_sum_preds equal the prediction matrix.
      overall_overall_sum_preds_tau=overall_sum_preds_tau;											// let overall_overall_sum_preds equal the prediction matrix.
      
      if(is_test_data==1){
        //Rcout << "Get to Line 5909 in loop j = " << j  << ".\n";
        overall_overall_sum_test_preds_outcome=overall_sum_test_preds_outcome;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
        overall_overall_sum_test_preds_mu=overall_sum_test_preds_mu;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
        overall_overall_sum_test_preds_tau=overall_sum_test_preds_tau;				// If there is test data, let overall_overall_sum_test_preds equal the out-of-sample prediction matrix
      }
      //overall_trees_mu[j]=curr_round_trees_mu;														// let the j+1^th element of the list overall_trees be the list of trees obtained in the current round of the outer loop (list of lists of tree table matrices)
      //overall_trees_tau[j]=curr_round_trees_tau;														// let the j+1^th element of the list overall_trees be the list of trees obtained in the current round of the outer loop (list of lists of tree table matrices)
      //overall_mat_mu.push_back(curr_round_mat_mu);													// append the list of current round of tree matrices to the list overall_mat (list of lists of matrices)
      //overall_mat_tau.push_back(curr_round_mat_tau);													// append the list of current round of tree matrices to the list overall_mat (list of lists of matrices)
      overall_lik.push_back(curr_round_lik);													// append the current round's vector of BICs to the list overall_lik (list of vectors of BICs).
      prev_par=seq_len(overall_sum_trees_mu.size())-1;											// let prev_par equal a sequence 0,1,2,3,..., up to length of overall_sum_trees minus one
      
      
      //Rcout << "Line 9654 overall_sum_trees_mu.size() = " << overall_sum_trees_mu.size() << ".\n";
      //Rcout << "overall_sum_trees_tau.size() = " << overall_sum_trees_tau.size() << ".\n";
      //Rcout << "overall_sum_tree_resids_mu.size() = " << overall_sum_tree_resids_mu.size() << ".\n";
      //Rcout << "overall_sum_tree_resids_tau.size() = " << overall_sum_tree_resids_tau.size() << ".\n";
      //Rcout << "curr_round_lik.size() = " << curr_round_lik.size() << ".\n";
      //Rcout << "overall_lik.size() = " << overall_lik.size() << ".\n";
      //Rcout << "overall_overall_sum_trees_mu.size() = " << overall_overall_sum_trees_mu.size() << ".\n";
      //Rcout << "overall_overall_sum_trees_tau.size() = " << overall_overall_sum_trees_tau.size() << ".\n";
      
      //Rcout << "overall_count = " << overall_count << ".\n";
      //Rcout << "Line 9574 oo_count = " << oo_count << ".\n";
      
      
    //}
    //	END OF mu(x) TREE CODE
    

    
    
      //Rcout << "Get to end of loop j = " << j << ".\n";
    
  } //	END OF OUTER LOOP
  
   //Rcout << "Get to outside outer loop.\n";
  
  
  if(oo_count==0){
    throw std::range_error("BCF-BMA did not find any suitable model for the data. Maybe limit for Occam's window is too small. Maybe use more observations or change parameter values.");
  }
  
  overall_overall_sum_trees_mu=resize_bcf(overall_overall_sum_trees_mu,oo_count);						// remove spaces that are not filled in.
  overall_overall_sum_trees_tau=resize_bcf(overall_overall_sum_trees_tau,oo_count);						// remove spaces that are not filled in.
  overall_overall_sum_tree_resids_mu=resize_bcf(overall_overall_sum_tree_resids_mu,oo_count);			// remove spaces that are not filled in.
  overall_overall_sum_tree_resids_tau=resize_bcf(overall_overall_sum_tree_resids_tau,oo_count);			// remove spaces that are not filled in.
  overall_overall_sum_trees_mat_mu=resize_bcf(overall_overall_sum_trees_mat_mu,oo_count);				// remove spaces that are not filled in.
  overall_overall_sum_trees_mat_tau=resize_bcf(overall_overall_sum_trees_mat_tau,oo_count);				// remove spaces that are not filled in.
  
  //Rcout << "Line 9601 overall_overall_sum_trees_mu.size() = " << overall_overall_sum_trees_mu.size() << ".\n";
  //Rcout << "overall_overall_sum_trees_tau.size() = " << overall_overall_sum_trees_tau.size() << ".\n";
  
  
  overall_overall_sum_BIC=resize_bcf(overall_overall_sum_BIC,oo_count);							// remove spaces that are not filled in.
  NumericVector end_BIC=overall_overall_sum_BIC[overall_overall_sum_BIC.size()-1] ;			// final element of overall_overall_sum_BIC (vector of BICs from final round)
  NumericMatrix overallpreds_outcome(n,end_BIC.size());												// create a vector of dimensions: number of training obs by number of models outptted by the final round.
  NumericMatrix overallpreds_mu(n,end_BIC.size());												// create a vector of dimensions: number of training obs by number of models outptted by the final round.
  NumericMatrix overallpreds_tau(n,end_BIC.size());												// create a vector of dimensions: number of training obs by number of models outptted by the final round.
  NumericMatrix overall_test_preds_outcome(test_data.nrow(),end_BIC.size());							// create a vector of dimensions: number of test obs by number of models outptted by the final round.
  NumericMatrix overall_test_preds_mu(test_data.nrow(),end_BIC.size());							// create a vector of dimensions: number of test obs by number of models outptted by the final round.
  NumericMatrix overall_test_preds_tau(test_data.nrow(),end_BIC.size());							// create a vector of dimensions: number of test obs by number of models outptted by the final round.
  NumericVector post_weights(end_BIC.size());													// create a vector of length equal to number of models outputted by final round.
  for(int k=0;k<end_BIC.size();k++){															// for-loop of length equal to number of models outputted in final round.
    NumericMatrix oosp_outcome=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_preds_outcome));				// create a matrix equal to overall_overall_sum_preds, the training prediction matrix
    NumericMatrix oosp_mu=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_preds_mu));				// create a matrix equal to overall_overall_sum_preds, the training prediction matrix
    NumericMatrix oosp_tau=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_preds_tau));				// create a matrix equal to overall_overall_sum_preds, the training prediction matrix
    NumericVector temp_preds_outcome=oosp_outcome(_,k);														// let temp_preds equal the predictions from the k+1^th model
    NumericVector temp_preds_mu=oosp_mu(_,k);														// let temp_preds equal the predictions from the k+1^th model
    NumericVector temp_preds_tau=oosp_tau(_,k);														// let temp_preds equal the predictions from the k+1^th model
    
    
    
    NumericVector temp_test_preds_outcome;															// create a vector called temp_test_preds
    NumericVector temp_test_preds_mu;															// create a vector called temp_test_preds
    NumericVector temp_test_preds_tau;															// create a vector called temp_test_preds
    if(is_test_data==1){																	// if there is test data
      //Rcout << "Get to Line 6749  "  << ".\n";
      NumericMatrix oostp_outcome=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_test_preds_outcome));		// let oostp equal the test data prediction matrix
      NumericMatrix oostp_mu=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_test_preds_mu));		// let oostp equal the test data prediction matrix
      NumericMatrix oostp_tau=Rcpp::as<NumericMatrix>(wrap(overall_overall_sum_test_preds_tau));		// let oostp equal the test data prediction matrix
      temp_test_preds_outcome=oostp_outcome(_,k);																// let temp_test_preds equal the out of sample predictions from the k+1^th model.
      temp_test_preds_mu=oostp_mu(_,k);																// let temp_test_preds equal the out of sample predictions from the k+1^th model.
      temp_test_preds_tau=oostp_tau(_,k);																// let temp_test_preds equal the out of sample predictions from the k+1^th model.
    }
    NumericVector orig_temp_preds_outcome=get_original_bcf(min(y),max(y),-0.5,0.5,temp_preds_outcome) ;			// Rescale the in-sample predictions back to the original scale of the outcome. Defined on line 2216
    NumericVector orig_temp_preds_mu=get_original_bcf(min(y),max(y),-0.5,0.5,temp_preds_mu) ;			// Rescale the in-sample predictions back to the original scale of the outcome. Defined on line 2216
    NumericVector orig_temp_preds_tau=get_original_bcf(min(y),max(y),-0.5,0.5,temp_preds_tau) ;			// Rescale the in-sample predictions back to the original scale of the outcome. Defined on line 2216
    NumericVector BICi=-0.5*end_BIC;														// create a vector of the BICs multiplied by -0.5
    //Rcout << "Line 9640 end_BIC = "<< end_BIC << " .\n";
    // Rcout << "end_BIC[0] = "<< end_BIC[0] << " .\n";
    
    // Rcout << "TEST LINE 6689 weight =  "<< weight << " .\n";
    
    double max_BIC=max(BICi);																// set the variable max_BIC equal to the maximum of the (negative 0.5 times the) BICs
    double weight=exp(BICi[k]-(max_BIC+log(sum(exp(BICi-max_BIC)))));						// create the weight for the k+1^th model
    
    post_weights[k]=weight;																	// Let the k+1^th element of post_weights be the weight of the k+1^th model
    overallpreds_outcome(_,k) = temp_preds_outcome*weight;													// Let the k+1^th element of overallpreds be the predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    overallpreds_mu(_,k) = temp_preds_mu*weight;													// Let the k+1^th element of overallpreds be the predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    overallpreds_tau(_,k) = temp_preds_tau*weight;													// Let the k+1^th element of overallpreds be the predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    if(is_test_data==1){																	// if there is test data
      //Rcout << "Get to Line 6774  "  << ".\n";
      overall_test_preds_outcome(_,k) = temp_test_preds_outcome*weight;											// Let the k+1^th element of overall_test_preds be the out-of-sample predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
      overall_test_preds_mu(_,k) = temp_test_preds_mu*weight;											// Let the k+1^th element of overall_test_preds be the out-of-sample predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
      overall_test_preds_tau(_,k) = temp_test_preds_tau*weight;											// Let the k+1^th element of overall_test_preds be the out-of-sample predictions of the k+1^th model multiplied by the model's weight (i.e, the contributions of the k+1^th model to the predictions).
    }
  }
  arma::mat M1_outcome(overallpreds_outcome.begin(), overallpreds_outcome.nrow(), overallpreds_outcome.ncol(), false);						// M1 is an arma mat copy of overallpreds (entry i,j gives contribution of j^th model to prediction of i^th observation)
  arma::mat M1_mu(overallpreds_mu.begin(), overallpreds_mu.nrow(), overallpreds_mu.ncol(), false);						// M1 is an arma mat copy of overallpreds (entry i,j gives contribution of j^th model to prediction of i^th observation)
  arma::mat M1_tau(overallpreds_tau.begin(), overallpreds_tau.nrow(), overallpreds_tau.ncol(), false);						// M1 is an arma mat copy of overallpreds (entry i,j gives contribution of j^th model to prediction of i^th observation)
  predicted_values_outcome=sum(M1_outcome,1);																					// predicted_values is a vector of predictions for each observation in the training data (before inverse scaling). (row sums of overallpreds)
  predicted_values_mu=sum(M1_mu,1);																					// predicted_values is a vector of predictions for each observation in the training data (before inverse scaling). (row sums of overallpreds)
  predicted_values_tau=sum(M1_tau,1);																					// predicted_values is a vector of predictions for each observation in the training data (before inverse scaling). (row sums of overallpreds)
  
  arma::mat M2_outcome(overall_test_preds_outcome.begin(), overall_test_preds_outcome.nrow(), overall_test_preds_outcome.ncol(), false);		// M2 is arma copy of overall_test_preds (entry i,j gives contribution of j^th model to prediction of i^th test observation)
  arma::mat M2_mu(overall_test_preds_mu.begin(), overall_test_preds_mu.nrow(), overall_test_preds_mu.ncol(), false);		// M2 is arma copy of overall_test_preds (entry i,j gives contribution of j^th model to prediction of i^th test observation)
  arma::mat M2_tau(overall_test_preds_tau.begin(), overall_test_preds_tau.nrow(), overall_test_preds_tau.ncol(), false);		// M2 is arma copy of overall_test_preds (entry i,j gives contribution of j^th model to prediction of i^th test observation)
  if(is_test_data==1){
    //Rcout << "Get to Line 6791  "  << ".\n";
    predicted_test_values_outcome=sum(M2_outcome,1);														// if there is test data, predicted_test_values is the vector of final test data predictions (before inverse scaling).
    predicted_test_values_mu=sum(M2_mu,1);														// if there is test data, predicted_test_values is the vector of final test data predictions (before inverse scaling).
    predicted_test_values_tau=sum(M2_tau,1);														// if there is test data, predicted_test_values is the vector of final test data predictions (before inverse scaling).
  }
  if(overall_lik.size()==0){																					// if the length of overall_lik is zero 
    throw std::range_error("BART-BMA didnt find any suitable model for the data. Maybe limit for Occam's window is too small.");
  }else{																										// if the length of overall_lik is greater than zero
    NumericVector orig_preds_outcome=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_values_outcome)) ;					// inverse scale in-sample predictions to original scale
    NumericVector orig_preds_mu=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_values_mu)) ;					// inverse scale in-sample predictions to original scale
    NumericVector orig_preds_tau=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_values_tau)) ;					// inverse scale in-sample predictions to original scale
    NumericVector orig_test_preds_outcome;																			// create vector
    NumericVector orig_test_preds_mu;																			// create vector
    NumericVector orig_test_preds_tau;																			// create vector
    if(is_test_data==1){																					// if have test data
      //Rcout << "Get to Line 6806  "  << ".\n";
      orig_test_preds_outcome=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_test_values_outcome)) ;					// inverse scale out-of-sample predictions to original scale
      orig_test_preds_mu=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_test_values_mu)) ;					// inverse scale out-of-sample predictions to original scale
      orig_test_preds_tau=get_original_bcf(min(y),max(y),-0.5,0.5,wrap(predicted_test_values_tau)) ;					// inverse scale out-of-sample predictions to original scale
    }
    NumericVector minmax(2);					// create vector minmax of length 2
    minmax[0]=min(y);							// set first element equal to min(y)
    minmax[1]=max(y);							// set first element equal to max(y)
    if(is_test_data==1){						// if there is test data
      List ret(13);									// list of length 6.
      ret[0] = orig_preds_outcome;							// The first element is a vector of in-sample predictions
      ret[1] = orig_preds_mu;							// The first element is a vector of in-sample predictions
      ret[2] = orig_preds_tau;							// The first element is a vector of in-sample predictions
      ret[3] = overall_overall_sum_trees_mu;			// the second element is a list of lists (of lists?) of tree tables
      ret[4] = overall_overall_sum_trees_tau;			// the second element is a list of lists (of lists?) of tree tables
      ret[5] =overall_overall_sum_trees_mat_mu;		// the third element is a list of lists of tree matrices.
      ret[6] =overall_overall_sum_trees_mat_tau;		// the third element is a list of lists of tree matrices.
      ret[7] = end_BIC;								// the fourth element is the vector of BICs of sum-of-tree-models
      ret[8] = orig_test_preds_outcome;						// the fifth element is the vector of out-of-sample predictions
      ret[9] = orig_test_preds_mu;						// the fifth element is the vector of out-of-sample predictions
      ret[10] = orig_test_preds_tau;						// the fifth element is the vector of out-of-sample predictions
      ret[11] =overall_overall_sum_tree_resids_mu;		// the sixth element if a list of residual vectors
      ret[12] =overall_overall_sum_tree_resids_tau;		// the sixth element if a list of residual vectors
      return(ret);							// return the list
    }else{										// if there is no test data
      List ret(10);									// list of length 5.
      ret[0] = orig_preds_outcome;							// The first element is a vector of in-sample predictions
      ret[1] = orig_preds_mu;							// The first element is a vector of in-sample predictions
      ret[2] = orig_preds_tau;							// The first element is a vector of in-sample predictions
      ret[3] = overall_overall_sum_trees_mu;				// the second element is a list of lists (of lists?) of tree tables
      ret[4] = overall_overall_sum_trees_tau;				// the second element is a list of lists (of lists?) of tree tables
      ret[5] = overall_overall_sum_trees_mat_mu;			// the third element is a list of lists of tree matrices.
      ret[6] = overall_overall_sum_trees_mat_tau;			// the third element is a list of lists of tree matrices.
      ret[7] = end_BIC;								// the fourth element is the vector of BICs of sum-of-tree-models
      ret[8] =overall_overall_sum_tree_resids_mu;		// the sixth element if a list of residual vectors
      ret[9] =overall_overall_sum_tree_resids_tau;		// the sixth element if a list of residual vectors
      return(ret);									// return the list
    }
  }
}