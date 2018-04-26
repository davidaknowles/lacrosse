#pragma once
#define _USE_MATH_DEFINES
// #include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include "settings.hpp"
#include "random_matrix.hpp"
#include "eigen_utils.hpp"
 
using namespace Eigen;
using namespace std; 

template<class real_type>
class params {
public:
  real_type alpha, beta, Bprior, factor_precision_rate, noise_precision_rate, loglike, test_loglike, test_rmse; 
  Matrix<real_type,Dynamic,Dynamic> Z, G, X, B; 
  Matrix<real_type,Dynamic,1> lambda_x, lambda_e;
  VectorXd FF;

  params(settings<real_type> &s)
  {
    alpha=s.alpha; 
    beta=s.beta; 
    if (s.lasso)
      Bprior=1.0; 
    else
      Bprior=min(0.8,5.0/(double)s.P); 
    factor_precision_rate=s.d; 
    noise_precision_rate=s.b; 
    //int initial_size=s.truncate ? s.D : 0; 
    Z=MatrixXd::Zero(s.D,s.initial_size); 
    G=MatrixXd::Zero(s.D,s.initial_size); 
    B=MatrixXd::Zero(s.initial_size,s.P); 
    lambda_x.setConstant(s.initial_size,s.lambda_x); 
    lambda_e.setConstant(s.D,s.lambda_e); 
    X=MatrixXd::Zero(s.initial_size,s.N); 
  }
  
  static real_type sqr(real_type x)
  {
    return x*x; 
  }

  static real_type greater_than_zero(real_type x)
  {
    return x > 0.0 ? 1.0 : 0.0; 
  }
 
  void delete_features(VectorXd &XX, MatrixXd &BF, VectorXd &m, Matrix<bool,Dynamic,1> &to_keep)
  {
    Z=EigenUtils::logical_index_cols(Z,to_keep); 
    X=EigenUtils::logical_index_rows(X,to_keep); 
    XX=EigenUtils::logical_index(XX,to_keep); 
    m=EigenUtils::logical_index(m,to_keep); 
    G=EigenUtils::logical_index_cols(G,to_keep);
    B=EigenUtils::logical_index_rows(B,to_keep); 
    BF=EigenUtils::logical_index_rows(BF,to_keep); 
    lambda_x=EigenUtils::logical_index(lambda_x,to_keep); 
  }

  int num_features(settings<real_type> &s){
  //return s.truncate ? Z.colwise().sum().unaryExpr(ptr_fun(greater_than_zero)).sum() : Z.cols(); 
  return Z.cols(); 
  }

  void one_MCMC_iteration(settings<real_type> &s, int iteration)
  {
    FF= s.F.unaryExpr(ptr_fun(sqr)).rowwise().sum(); 
    MatrixXd BF=B*s.F; 
    sample_factor_loadings(s,BF);
    sample_parametric_part(s);
    if (Z.cols()>0){
      sample_B(s,BF,B);
      if (!s.lasso)
	sample_parametric_part_B(s);
      sample_factors(s,BF);
    }
    if (s.sample_noise_variance)
      sample_noise_hypers(s); 
    if (s.sample_alpha)
      sample_alpha(s); 
    if (s.sample_Bprior)
      sample_Bprior(s); 
    if (s.sample_mixing_matrix_hypers)
      sample_mixing_matrix_hypers(s,BF); 
    if (s.sample_beta)
      Rcpp::stop( "ERROR: sampling beta not implemented yet" ); 
    if (s.missing_values)
      impute_Y(s); 
    loglike = log_likelihood(s); 
    //    double loglikecheck = slow_ll(s); 
    // Rcpp::Rcout << loglikecheck  << endl; 
    stringstream ss;
    double nonZeros=B.unaryExpr(ptr_fun(binarise)).sum(); 
    ss << iteration << "/" << s.iterations << " K=" << num_features(s) << " mZ=" << Z.rowwise().sum().mean() << " nB="<< nonZeros << " ll=" << loglike;
    ss << " alpha=" << alpha << " Bprior " << Bprior << " lambda_e[1]=" << lambda_e[0]; 
    if (s.predictive_performance)
      ss << " testl=" << test_loglike << " rmse=" << test_rmse; 
    ss << endl;
    string st = ss.str(); 
    if (s.verbose)
      Rcpp::Rcout << st; 
    if (!s.log_file.empty()){
      ofstream logfile(s.log_file.c_str(),ios::app);
      logfile << st;
      logfile.close();
    }
  }

  real_type log_likelihood(settings<real_type> &s)
  {
    const static real_type log2pi= log((real_type)(2.0f*M_PI)); 
    VectorXi row_counts(s.D) ;
    VectorXd log_lambda_2pi(s.D);
    for (int i=0;i<s.D;i++){
      row_counts(i)=s.N- (s.missing_values ? s.missing.row(i).count() : 0); 
      log_lambda_2pi(i)=log(lambda_e(i))-log2pi; 
    }

    MatrixXd GX = G*X; 
    MatrixXd E = (*(s.Y_original))-GX; 
    MatrixXd E_obs = E; 
    if (s.missing_values)
      E_obs=s.missing.select(0.0,E); 
    E_obs=E_obs.unaryExpr( ptr_fun(sqr) ); 
    if (s.predictive_performance){
      MatrixXd E_heldout = s.missing.select(*(s.Y_original)-GX, 0.0);
      E_heldout=E_heldout.unaryExpr( ptr_fun(sqr)); 
      VectorXi row_counts_missing=VectorXi::Constant(s.D,s.N) - row_counts; 
      test_loglike=.5f*row_counts_missing.cast<real_type>().dot(log_lambda_2pi) - .5f*E_heldout.rowwise().sum().dot(lambda_e); 

      test_rmse = sqrt(E_heldout.sum() / (real_type)s.missing.count()); 
      test_loglike /= (real_type)s.missing.count(); 
    }

    return .5f*row_counts.cast<real_type>().dot(log_lambda_2pi) - .5f * E_obs.rowwise().sum().dot(lambda_e);
  }

  real_type slow_ll(settings<real_type> &s){
    const static real_type log2pi= log((real_type)(2.0*M_PI)); 
    MatrixXd GX= G*X;
    real_type ll = 0.0;
    test_loglike=0.0;
    test_rmse=0.0;

    for (int d=0;d<s.D;d++)
      for (int n=0;n<s.N;n++){
	real_type e2 = (*s.Y_original)(d,n) - GX(d,n); 
	e2 *= e2; 
	real_type l = .5*( log(lambda_e(d)) - log2pi - lambda_e(d)*e2 ); 
	assert( isnormal(l));

	if (s.missing_values && s.missing(d,n)){
	  test_loglike+=l; 
	  test_rmse += e2; 
	}
	else
	  ll += l; 
      }
    test_loglike /= (real_type)s.missing.count(); 
    test_rmse = sqrt( test_rmse / (real_type)s.missing.count() ); 
    return ll; 
  }

  void impute_Y(settings<real_type> &s) {
    MatrixXd r=randn_matrix<real_type>(s.D,s.N);
    VectorXd l=lambda_e.unaryExpr(pointer_to_unary_function<real_type,real_type>(sqrt)); 
    MatrixXd noise = r.array().colwise() / l.array(); 
    MatrixXd Ypred(G*X); 
    Ypred += noise; 
    *(s.Y)=s.missing.select(Ypred,*(s.Y)); 
  }

  void sample_factors(settings<real_type> &s,MatrixXd &BF) {
    int K=G.cols(); 
    MatrixXd lambdaG = G.array().colwise() * lambda_e.array(); // TODO: check this
    assert( lambdaG.cols()==G.cols()); 
    MatrixXd prec=G.transpose()*lambdaG;
    prec+=lambda_x.asDiagonal(); 
    LLT<MatrixXd> llt(prec); 
    MatrixXd mp = lambdaG.transpose()*(*(s.Y))+lambda_x.asDiagonal()*BF; 
    X=llt.solve(mp);
    X+=llt.matrixL().solve(randn_matrix<real_type>(K,s.N));
    assert( ! isnan(X(0,0)));
  }

  void sample_one_element(real_type m,int di,int ki,VectorXd &XX, MatrixXd &E,settings<real_type> &s)
  {
    RowVectorXd ek(s.N); 
    if (Z(di,ki)==1.0)
      ek=E.row(di)+G(di,ki)*X.row(ki); 
    else
      ek=E.row(di); 
    real_type lambda=lambda_e(di)*XX(ki)+1.0; 
    // TODO most time seems to be spent on the following line, so cache lambda_e XE? 
    // or could experiment with row-major vs col-major? 
    real_type mu=lambda_e(di)*X.row(ki).dot(ek)/lambda; 
    real_type K=(real_type)Z.cols(); 
    real_type logrp= log(m + (s.truncate ? alpha/K : 0.0)) - log(beta+s.D-1.0-m) ; 
    logrp+=.5*(lambda*mu*mu-log(lambda)+log(1.0)); 

    // contribution from mrf
    if (s.use_mrf)
      for (int i=0;i<s.mrf_connections[di].size();i++)
	      if (Z(s.mrf_connections[di][i]-1,ki)==1.0)
	        logrp+=s.mrf_strength[di][i];

    real_type prob_zis1=1.0/(1.0+exp(-logrp)); 
    assert( !isnan(prob_zis1));
    bool z=s.rand()<prob_zis1;
	     
    Z(di,ki)=z ? 1.0 : 0.0; 
    G(di,ki)=z ? mu+s.randn()/sqrt(lambda) : 0.0; 
    E.row(di)=ek-G(di,ki)*X.row(ki);     
  }

  static bool not_op(bool x) { return !x; } 

  void sample_new_features(int di,VectorXd &XX,MatrixXd &E,MatrixXd &BF, VectorXd &mall,settings<real_type> &s) 
  {
    int Knew=s.rpoiss(alpha*beta/(beta+(real_type)s.D-1)); 
    VectorXd lambda_x_new; 
    if (s.sample_mixing_matrix_hypers) {
      if (s.perfactor_variance)
	lambda_x_new=gamrnd_matrix<double>(Knew,1,s.c,factor_precision_rate);
      else {
	if (Z.cols()==0)
	  lambda_x_new=VectorXd::Constant(Knew,s.lambda_x);
	else
	  lambda_x_new=VectorXd::Constant(Knew,lambda_x(0));
      }
    }
    else
      lambda_x_new=VectorXd::Constant(Knew,s.lambda_x);
    VectorXd m = mall-Z.row(di).transpose(); 
    
    Matrix<bool,Dynamic,1> singletons = m.array()==0.0; 
    int current_kappa=EigenUtils::binary_sum(singletons); 
    // TODO: if current_kappa==Knew just return? 
    VectorXd Gtemp=G.row(di); 
    VectorXd current_g=EigenUtils::logical_index(Gtemp,singletons); 
    VectorXd current_lambda_x = EigenUtils::logical_index(lambda_x,singletons); 
    MatrixXd current_bf = EigenUtils::logical_index_rows(BF,singletons); 
    Gtemp = singletons.select(0.0,Gtemp); 
    real_type logrprop = 0.0; 
    RowVectorXd Ed = s.Y->row(di); 
    MatrixXd prec,M;
    LLT<MatrixXd> llt;
    if (current_kappa > 0 | Knew > 0)
      Ed -= Gtemp.transpose()*X; 
    if (current_kappa > 0) {
      prec = lambda_e(di)*current_g*current_g.transpose();
      prec += current_lambda_x.asDiagonal();
      llt.compute(prec); 
      M=lambda_e(di)*(llt.solve(current_g)*Ed)+llt.solve(current_lambda_x.asDiagonal()*current_bf); 
      logrprop = .5*s.N*EigenUtils::log_determinant(llt)-.5*(M.transpose()*prec*M).diagonal().sum(); 
    }
    VectorXd g;
    MatrixXd new_b,new_bf; 
    if (Knew > 0){
      g = randn_matrix<real_type>(Knew,1).array();
      // TODO do this properly 
      new_b=MatrixXd::Zero(Knew,s.P); 
      new_bf=new_b*s.F; 
      
      prec=lambda_e(di)*(g*g.transpose());
      prec+=lambda_x_new.asDiagonal(); 
      llt.compute(prec); 
      M=lambda_e(di)*llt.solve(g)*Ed+llt.solve(lambda_x_new.asDiagonal()*new_bf);
      logrprop-=.5*s.N*EigenUtils::log_determinant(llt)-.5*(M.transpose()*prec*M).diagonal().sum(); 
    }
    assert( !isnan(logrprop) );
    // decide whether to accept the change
    if (s.rand()<exp(logrprop)){
      if (current_kappa>0){
	E.row(di)+=current_g.transpose()*EigenUtils::logical_index_rows(X,singletons); 
	Matrix<bool,Dynamic,1> to_keep=singletons.unaryExpr(ptr_fun(not_op)); 
	delete_features(XX,BF,mall,to_keep);
      }
      if (Knew>0){
	int K=Z.cols(); 
	Z.conservativeResize(s.D,K+Knew);
	Z.block(0,K,s.D,Knew)=MatrixXd::Zero(s.D,Knew); 
	Z.block(di,K,1,Knew)=MatrixXd::Ones(1,Knew); 
	X.conservativeResize(K+Knew,s.N);
	X.block(K,0,Knew,s.N)=M+llt.matrixL().solve(randn_matrix<real_type>(Knew,s.N));
	XX.conservativeResize(K+Knew); 
	mall.conservativeResize(K+Knew); 
	mall.segment(K,Knew)=VectorXd::Ones(Knew); 
	XX.segment(K,Knew)=X.block(K,0,Knew,s.N).unaryExpr(ptr_fun(sqr)).rowwise().sum(); 
	G.conservativeResize(s.D,K+Knew); 
	G.block(0,K,s.D,Knew)=MatrixXd::Zero(s.D,Knew); 
	G.block(di,K,1,Knew)=g.transpose(); 
	B.conservativeResize(K+Knew,s.P); 
	B.block(K,0,Knew,s.P)=new_b; 
	BF.conservativeResize(K+Knew,s.N); 
	BF.block(K,0,Knew,s.N)=new_bf; 
	lambda_x.conservativeResize(K+Knew); 
	lambda_x.segment(K,Knew)=lambda_x_new;
	E.row(di)-=g.transpose()*X.block(K,0,Knew,s.N); 
	K=K+Knew; 
	// TODO: sample new g again here?
      }
    }
  }

  void sample_factor_loadings(settings<real_type> &s,MatrixXd &BF)
  {
    VectorXd XX = X.unaryExpr(ptr_fun(sqr)).rowwise().sum(); 
    MatrixXd E = (*(s.Y)) - G*X; 
    VectorXd m=Z.colwise().sum();
    for (int di=0;di<s.D;di++){
      int K=Z.cols(); 
      for (int ki=0;ki<K;ki++){
	m(ki)-=Z(di,ki); 
	if (s.truncate || (m(ki)>0))
	  sample_one_element(m(ki),di,ki,XX,E,s);
	m(ki)+=Z(di,ki); 
	if (s.debug)
	  assert( (E.row(di) - (s.Y->row(di)-G.row(di)*X) ).norm() < 1e-6 ); 
      }
      if (!s.truncate)
	sample_new_features(di,XX,E,BF,m,s); 
      if (s.debug){
	assert( (E.row(di) - (s.Y->row(di)-G.row(di)*X) ).norm() < 1e-6 ); 
	assert( (Z.colwise().sum()-m.transpose()).norm() < 1e-6 );
	assert( (X.unaryExpr(ptr_fun(sqr)).rowwise().sum() - XX).norm() < 1e-6 );
      }
    }
  }

  void sample_parametric_part(settings<real_type> &s){
    MatrixXd XX=X * X.transpose(); 
    MatrixXd XY=X * (*(s.Y)).transpose(); 
    for (int di=0;di<s.D;di++){
      VectorXd t=Z.row(di).transpose(); 
      if (t.sum()==0.0)
	continue; 
      Matrix<bool,Dynamic,1> z= EigenUtils::double_to_bool(t);
      MatrixXd prec=EigenUtils::logical_index_cols(EigenUtils::logical_index_rows(XX,z),z)*lambda_e[di];
      int numEl=prec.cols(); 
      if (numEl==0)
	continue; 
      prec+=MatrixXd::Identity(numEl,numEl);   // EigenUtils::logical_index(lambda_g,z).asDiagonal(); 
      VectorXd xy=XY.col(di); 
      VectorXd mean_prec=EigenUtils::logical_index(xy,z)*lambda_e[di];
      LLT<MatrixXd> llt(prec); 

      VectorXd temp=llt.solve(mean_prec); 
      temp+=llt.matrixL().solve(randn_matrix<real_type>(numEl,1));
      EigenUtils::logical_index_set(G,di,z,temp); 
    }
  }
  
  real_type soft_threshold(real_type x, real_type l){
    if (x > l)
      return x-l; 
    else if (x < -l)
      return x+l;
    else return 0.0; 
  }

  void sample_one_element_B(int ki,int pi,MatrixXd &BF,MatrixXd &B,settings<real_type> &s)
  {
    if (B(ki,pi)!=0.0)
      BF.row(ki) -= B(ki,pi)*s.F.row(pi); 
    RowVectorXd ek = X.row(ki) - BF.row(ki); 

    if (s.lasso){
      real_type betahat=s.F.row(pi).dot(ek);
      B(ki,pi)=soft_threshold(betahat, Bprior/lambda_x(ki))/FF(pi);
    }
    else {
      real_type lambda=lambda_x(ki)*FF(pi)+1.0; 
      real_type mu=lambda_x(ki)*s.F.row(pi).dot(ek)/lambda; 
      real_type logrp=log(Bprior)-log(1.0-Bprior); // -log((real_type)s.P); // log(m+(1.0/(double)s.P)) - log(s.N-1.0+m); 
      assert( !isnan(logrp) ); 
      logrp+=.5*(lambda*mu*mu-log(lambda)+log(1.0)); 
      
      if (s.use_feature_mrf)
	for (int i=0;i<s.feature_mrf_connections[pi].size();i++)
	  if (B(ki,s.feature_mrf_connections[pi][i]-1)!=0.0)
	    logrp+=s.feature_mrf_strength[pi][i];
      
      real_type prob_zis1=1.0/(1.0+exp(-logrp)); 
      assert( !isnan(prob_zis1));
      bool z=s.rand()<prob_zis1;
	     
      B(ki,pi)=z ? mu+s.randn()/sqrt(lambda) : 0.0; 
    }

    if (B(ki,pi)!=0.0)
      BF.row(ki) += B(ki,pi)*s.F.row(pi); 
  }

  void sample_B(settings<real_type> &s,MatrixXd &BF,MatrixXd &B)
  {
    for (int ki=0;ki<X.rows();ki++)
      for (int pi=0;pi<s.P;pi++)
	sample_one_element_B(ki,pi,BF,B,s);
    //if (s.debug) 
    // assert( ( B.row(ki)*s.F - BF ).norm() < 1e-6); 
  }

  void sample_parametric_part_B(settings<real_type> &s){
    for (int ki=0;ki<Z.cols();ki++){
      VectorXd t=B.row(ki).transpose(); 
      if (t.sum()==0.0)
	continue; 
      Matrix<bool,Dynamic,1> z=EigenUtils::double_to_bool(t);
      /* if (EigenUtils::binary_sum(z)>1000) {
	      Rcpp::Rcout << "Big! : " << EigenUtils::binary_sum(z) << " ki " << ki << endl;
	      Rcpp::Rcout << "B[ki,]=" << B.row(ki) << endl; 
      } */
      MatrixXd f=EigenUtils::logical_index_rows(s.F,z); 
      MatrixXd prec=f*f.transpose()*lambda_x[ki];
      int numEl=prec.cols(); 
      if (numEl==0)
	continue; 
      prec+=MatrixXd::Identity(numEl,numEl);   // EigenUtils::logical_index(lambda_g,z).asDiagonal(); 
      VectorXd mean_prec=f * X.row(ki).transpose() * lambda_x[ki]; 
      assert(mean_prec.size()==numEl); 
      LLT<MatrixXd> llt(prec); 
      VectorXd temp=llt.solve(mean_prec); 
      temp+=llt.matrixL().solve(randn_matrix<real_type>(numEl,1));
      EigenUtils::logical_index_set(B,ki,z,temp); 
    }
  }

  void sample_noise_hypers(settings<real_type> &s)
  {
    MatrixXd E = (*(s.Y)) - G*X;
    // TODO: E = E.array() * mvmask; total_observations=mvmask.sum()
    int total_observations = s.N*s.D; 
    if (s.isotropic_noise){
      lambda_e=VectorXd::Ones(s.D).array()*s.gamrnd(s.a+.5*(real_type)total_observations,noise_precision_rate+.5*E.squaredNorm() );
    }
    else {
      for (int d=0;d<s.D;d++) {
	lambda_e(d)=s.gamrnd(s.a+.5*(real_type)s.N,noise_precision_rate+.5*E.row(d).squaredNorm()); 
      }
      if (s.hierarchical_noise_prior){
	noise_precision_rate=s.gamrnd(s.a0+s.D*s.a,s.b0+lambda_e.sum());
      }
    }
  }

  void sample_mixing_matrix_hypers(settings<real_type> &s,MatrixXd &BF)
  {
    int K=Z.cols();
    MatrixXd err=X-BF; 
    if (s.perfactor_variance){
      for (int k=0;k<K;k++){
	lambda_x(k)=s.gamrnd(s.c+.5*(real_type)s.N,factor_precision_rate+.5*err.row(k).squaredNorm()); 
      }
      if (s.learn_scale){
	factor_precision_rate=s.gamrnd(s.c0+s.c*(real_type)K,s.d0+lambda_x.sum()); 
      }
    }
    else {
      lambda_x=VectorXd::Ones(K)*s.gamrnd(s.c+.5*(real_type)K*s.N,s.d+.5*err.squaredNorm()); 
    }
  }

  static real_type HD(real_type x, int D)
  {
    real_type res=0.0;
    for (int i=0;i<D;i++)
      res += x/(x+(real_type)i); 
    return res; 
  }

  static double binarise(double x){
    return x==0.0 ? 0.0 : 1.0; 
  }

  void sample_alpha(settings<real_type> &s)
  {
    int K=s.truncate ? Z.colwise().sum().unaryExpr(ptr_fun(greater_than_zero)).sum() : Z.cols(); 
    alpha=s.gamrnd(s.e+(double)K,s.f+HD(beta,s.D)); 

  }

  static double reciprocal(double x){ return 1.0/x; } 

  double marginaliseX(settings<real_type> &s,MatrixXd &BF){
    /* MatrixXd temp(BF-X);
      VectorXd err=temp.unaryExpr(ptr_fun(sqr)).rowwise().sum();
      if (err.size()!=X.rows()) throw 1; 
      return (err.array() * lambda_x.array()).sum(); */ 
    MatrixXd lambdaG = G.array().rowwise() * lambda_x.unaryExpr(ptr_fun(reciprocal)).transpose().array();      
    MatrixXd Lambda = lambdaG * G.transpose(); 
    Lambda += lambda_e.unaryExpr(ptr_fun(reciprocal)).asDiagonal(); 
    LLT<MatrixXd> llt(Lambda);
    MatrixXd GBF(G*BF);
    GBF-=(*(s.Y)); 
    return ( GBF.array() * llt.solve(GBF).array() ).sum(); 
  }

  double marginaliseX1(settings<real_type> &s,MatrixXd &BF){
    MatrixXd lambdaG = G.array().rowwise() * lambda_x.unaryExpr(ptr_fun(reciprocal)).transpose().array();      
    MatrixXd Lambda = lambdaG * G.transpose(); 
    Lambda += lambda_e.unaryExpr(ptr_fun(reciprocal)).asDiagonal(); 
    LLT<MatrixXd> llt(Lambda);
    MatrixXd GBF(G*BF);
    GBF-=(*(s.Y)); 
    return -2.0*( (*(s.Y)).array()*llt.solve(GBF).array()).sum(); 
  }

  double aic(settings<real_type> &s,MatrixXd &BF,MatrixXd &B){
    double nonZeros=B.unaryExpr(ptr_fun(binarise)).sum(); 
    if (s.debug) Rcpp::Rcout << "nz: " << nonZeros << endl; 
    return marginaliseX(s,BF) + log((double)X.cols()) * nonZeros / (double)X.cols();
  }

  double fit_quadratic(double x[], double y[], double &quad_term){
    // double g=(y[1]-y[0])/(x[1]-x[0]); 
    //double a2=(y[2]-y[0]-g*x[2])/((x[2]-x[0])*(x[2]-x[1]));
    double a2= ((y[1]-y[0])*(x[0]-x[2]) + (y[2]-y[0])*(x[1]-x[0]))/((x[0]-x[2])*(x[1]*x[1]-x[0]*x[0]) + (x[1]-x[0])*(x[2]*x[2]-x[0]*x[0])); 
    //    double a1=-x[0]*a2-x[1]*a2+g; 
    double a1= ((y[1] - y[0]) - a2*(x[1]*x[1] - x[0]*x[0])) / (x[1]-x[0]);
    quad_term=a2; 
    return a1; 
  }

  void test_quad_fit(){
    double x[3] = { -1.0, 0.0, 1.0 } ; 
    double a2=3.4; 
    double a1=-2.3; 
    double a0=0.9, q; 
    double y[3]; 
    for (int i=0;i<3;i++) y[i]=a2*x[i]*x[i]+a1*x[i]+a0; 
    double a1i=fit_quadratic(x,y,q);
    if ((abs(a1-a1i)>1e-4)||(abs(a2-q)>1e-4)) Rcpp::stop("Quadratic fitting failed.");
  }

  void check_calcs(settings<real_type> &s,MatrixXd &BF,MatrixXd &B,int ki,int pi){
    double y[3], x[3] = { -1.0, 0.0, 1.0 } ; 
    for (int i=0;i<3;i++){
      BF.row(ki) += x[i]*s.F.row(pi); 
      y[i]=marginaliseX(s,BF); 
      BF.row(ki) -= x[i]*s.F.row(pi); 
    }
    double quad; 
    double lin=fit_quadratic(x,y,quad); 
    Rcpp::Rcout << "Quad: "<< quad << " lin: "<< -.5*lin << endl; 
  }

  double opt_B(settings<real_type> &s,MatrixXd &B)
  {
    if (s.debug) test_quad_fit();
    MatrixXd BF(B*s.F); 
    MatrixXd lambdaG = G.array().rowwise() * lambda_x.unaryExpr(ptr_fun(reciprocal)).transpose().array();      
    MatrixXd Lambda = lambdaG * G.transpose(); 
    Lambda += lambda_e.unaryExpr(ptr_fun(reciprocal)).asDiagonal(); 
    LLT<MatrixXd> llt(Lambda);
    MatrixXd LG=llt.solve(G); 
    MatrixXd YLG=(*(s.Y)).transpose()*LG; 
    MatrixXd A=G.transpose()*LG;
    double old_aic=aic(s,BF,B); 

    for (int i=0;i<30;i++){
      for (int ki=0;ki<X.rows();ki++)
	for (int pi=0;pi<s.P;pi++){
	  double old_lk; 
	  double old_b=B(ki,pi); 
	  if (s.debug) old_lk = .5 * marginaliseX(s,BF) + Bprior * abs(B(ki,pi)); 
	  if (B(ki,pi)!=0.0)
	    BF.row(ki) -= B(ki,pi)*s.F.row(pi); 
	  if (A(ki,ki)!=0.0){
	    double quad_term=FF(pi)*A(ki,ki); 
	    double lin_term=YLG.col(ki).dot(s.F.row(pi)); 
	    lin_term-=s.F.row(pi).dot(A.row(ki)*BF); // note BF has B(ki,pi) set to 0 already, so don't need to remove quadratic term	    
	    /*if (s.debug){
	      Rcpp::Rcout << "Analytic quad: " << quad_term << " lin: " << lin_term << endl; 
	      check_calcs(s,BF,B,ki,pi); 
	      }*/
	    
	    B(ki,pi)=soft_threshold( lin_term, Bprior) / quad_term; 
	    /*RowVectorXd ek = X.row(ki) - BF.row(ki); 
	      real_type betahat=s.F.row(pi).dot(ek)/FF(pi);
	      B(ki,pi)=soft_threshold(betahat, Bprior/lambda_x(ki)); */
	  } else {
	    B(ki,pi)=0.0; 
	  }
	  if (B(ki,pi)!=0.0)
	    BF.row(ki) += B(ki,pi)*s.F.row(pi); 
	  if (s.debug){
	    double new_lk=.5 * marginaliseX(s,BF) + Bprior * abs(B(ki,pi));  
	    if (new_lk-0.0001 > old_lk)
	      Rcpp::Rcout << "new: " << new_lk << " old: " << old_lk << " oldB " << old_b << " newB " << B(ki,pi) << endl; 
	  }
	}
      double new_aic=aic(s,BF,B);
      if (s.debug) Rcpp::Rcout << "AIC: " << new_aic << endl; 
      if ((i>1) && (abs(new_aic-old_aic)<1.0)){
	old_aic=new_aic; 
	if (s.debug) Rcpp::Rcout << "Converged after " << i << " its" << endl; 
	break; 
      }
      old_aic=new_aic; 
    }
    return old_aic; 
    
  }
  

  double convergeB(settings<real_type> &s,MatrixXd &B){
      MatrixXd BF(B*s.F); 
      double old_aic=aic(s,BF,B); 
      for (int i=0;i<30;i++){
	sample_B(s,BF,B); 
	double new_aic=aic(s,BF,B);
	if (s.debug) Rcpp::Rcout << "AIC: " << new_aic << endl; 
	if ((i>1) && (abs(new_aic-old_aic)<1.0)){
	  old_aic=new_aic; 
	  if (s.debug) Rcpp::Rcout << "Converged after " << i << " its" << endl; 
	  break; 
	}
	old_aic=new_aic; 
      }
      return old_aic; 
  }

  void sample_Bprior(settings<real_type> &s){
    Matrix<bool,Dynamic,1> temp= Z.colwise().sum().array() != 0; 
    if (s.lasso){
      double m=0.0,p=.01; 
      double current_aic=opt_B(s,B) + sqr(log(Bprior)-m)*p*.5; 
      double old_Bprior=Bprior; 
      Bprior=exp(log(Bprior)+3.0*s.randn());
      MatrixXd Bprime(B);
      double new_aic=opt_B(s,Bprime) + sqr(log(Bprior)-m)*p*.5; 
      if (new_aic < current_aic){
	B=Bprime; 
      } else {
	Bprior=old_Bprior; 
      }
      //      VectorXd l1s=B.unaryExpr(pointer_to_unary_function<real_type,real_type>(abs)).rowwise().sum(); 
      //VectorXd relevant_l1s=EigenUtils::logical_index(l1s,temp); 
      //double onSumL1=relevant_l1s.sum();
      //Bprior = s.gamrnd((double)s.P*EigenUtils::binary_sum(temp)+1.0,onSumL1+1.0);
    } else {
      VectorXd counts=B.unaryExpr(ptr_fun(binarise)).rowwise().sum(); 
      VectorXd relevant_counts=EigenUtils::logical_index(counts,temp); 
      double onSum=relevant_counts.sum();
      Bprior=s.betarnd(onSum+1.0/(double)s.P,(double)s.P*(double)relevant_counts.size()-onSum+1.0); 
      assert( Bprior > 0.0 );
      Bprior=max(Bprior, .3/(double)s.P); 
    }
  }
};
