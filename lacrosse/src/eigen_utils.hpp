#pragma once

#include <RcppEigen.h>

#include <math.h>
#include <assert.h>

#include "settings.hpp"
#include "random_matrix.hpp"

using namespace Eigen;
using namespace std;

class EigenUtils {

  static int bool_to_int(bool x)
  {
    return x ? 1 : 0;
  }

  static bool _double_to_bool(double x){
    return x != 0.0; 
  }

public:  
  static int binary_sum(const Matrix<bool,Dynamic,1> &x)
  {
    return x.unaryExpr(ptr_fun(bool_to_int)).sum();
  }

  static Matrix<bool,Dynamic,1> double_to_bool(VectorXd &x){
    return x.unaryExpr(pointer_to_unary_function<double,bool>(_double_to_bool)); 
  }

  template<class scalar_type>
  static Matrix<scalar_type,Dynamic,Dynamic> logical_index_rows(const Matrix<scalar_type,Dynamic,Dynamic> &input, const Matrix<bool,Dynamic,1> &index)
  {
    int len = binary_sum(index); 
    Matrix<scalar_type,Dynamic,Dynamic> result(len,input.cols());
    int counter=0; 
    for (int i=0;i<input.rows();i++){
      if (index(i)){
	result.row(counter)=input.row(i); 
	counter++; 
      }
    }
    return result;
  }

  template<class scalar_type>
  static Matrix<scalar_type,Dynamic,Dynamic> logical_index_cols(const Matrix<scalar_type,Dynamic,Dynamic> &input, const Matrix<bool,Dynamic,1> &index)
  {
    int len = binary_sum(index); 
    Matrix<scalar_type,Dynamic,Dynamic> result(input.rows(),len);
    int counter=0; 
    for (int i=0;i<input.cols();i++){
      if (index(i)){
	result.col(counter)=input.col(i); 
	counter++; 
      }
    }
    return result;
  }

  template<class scalar_type>
  static Matrix<scalar_type,Dynamic,1> logical_index(const Matrix<scalar_type,Dynamic,1> &input, const Matrix<bool,Dynamic,1> &index)
  {
    int len = binary_sum(index); 

    Matrix<scalar_type,Dynamic,1> result(len,1);
    int counter=0; 
    for (int i=0;i<input.rows();i++){
      if (index(i)){
	result(counter)=input(i); 
	counter++; 
      }
    }
    return result;
  }

  template<class scalar_type>
  static void logical_index_set(Matrix<scalar_type,Dynamic,Dynamic> &to_set, int row_index, const Matrix<bool,Dynamic,1> &index, const Matrix<scalar_type,Dynamic,1> &x)
  {
    int counter=0; 
    for (int i=0;i<to_set.cols();i++){
      if (index(i)){
	to_set(row_index,i)=x(counter); 
	counter++; 
      }
    }
  }

  template<class scalar_type>
  static scalar_type log_determinant(const LDLT<Matrix<scalar_type,Dynamic,Dynamic> > &chol)
  {
    return chol.vectorD().unaryExpr(pointer_to_unary_function<double,double>(log)).sum();
  }

  template<class scalar_type>
  static scalar_type log_determinant(const LLT<Matrix<scalar_type,Dynamic,Dynamic> > &chol)
  {
    return 2.0*chol.matrixLLT().diagonal().unaryExpr(pointer_to_unary_function<double,double>(log)).sum();
  }

};
