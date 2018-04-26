#pragma once

#include <RcppEigen.h>

using namespace std; 
using namespace boost; 
using namespace Eigen; 

// Will only actually work for double... couldn't get Boost approach to work
template<class scalar_type>
Matrix<scalar_type,Dynamic,Dynamic> randn_matrix(int rows, int cols,scalar_type mean=scalar_type(0.0),scalar_type sd=scalar_type(1.0))
{
  MatrixXd result(rows,cols);
  result << Rcpp::as<VectorXd>(Rcpp::rnorm(rows * cols,mean,sd));
  return(result);
}

template<class scalar_type>
Matrix<scalar_type,Dynamic,Dynamic> gamrnd_matrix(int rows, int cols,scalar_type shape,scalar_type rate)
{
  MatrixXd result(rows,cols);
  result << Rcpp::as<VectorXd>(Rcpp::rgamma(rows * cols,shape,rate));
  return(result);
}
