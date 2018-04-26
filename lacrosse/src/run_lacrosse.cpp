#include "run_lacrosse.h"
#include <assert.h>
using namespace Rcpp ;
using namespace Eigen; 

void setup_mrf(SEXP mrf_connections_s, SEXP mrf_strength_s, settings<double> &s){
  List mrf_connections_list(mrf_connections_s); 
  List mrf_strength_list(mrf_strength_s); 

  if (mrf_connections_list.size()==0){
    s.use_mrf=false; 
  } else {
    s.use_mrf=true; 
    int D=mrf_connections_list.size();
    s.mrf_connections.resize(D); 
    s.mrf_strength.resize(D); 
    for (int d=0;d<D;d++){
      NumericVector nv((SEXP)mrf_strength_list[d]);
      s.mrf_strength[d]=as<vector<double> >(nv); 
      IntegerVector iv((SEXP)mrf_connections_list[d]); 
      s.mrf_connections[d]=as<vector<int> >(iv); 
    }
  }
}

void setup_feature_mrf(SEXP mrf_connections_s, SEXP mrf_strength_s, settings<double> &s){
  List mrf_connections_list(mrf_connections_s); 
  List mrf_strength_list(mrf_strength_s); 
  if (mrf_connections_list.size()==0){
    s.use_feature_mrf=false; 
  } else {
    s.use_feature_mrf=true; 
    int K=mrf_connections_list.size(); 
    s.feature_mrf_connections.resize(K); 
    s.feature_mrf_strength.resize(K); 
    for (int d=0;d<K;d++){
      NumericVector nv((SEXP)mrf_strength_list[d]);
      s.feature_mrf_strength[d]=as<vector<double> >(nv); 
      IntegerVector iv((SEXP)mrf_connections_list[d]); 
      s.feature_mrf_connections[d]=as<vector<int> >(iv); 
    }
  }
}
 
settings<double> parse_settings(SEXP settings_sexp){
  Rcpp::List rparam(settings_sexp);
  settings<double> s; 
  s.iterations=as<int>(rparam["iterations"]); 
  s.truncate=as<bool>(rparam["truncate"]); 
  s.a=as<double>(rparam["a"]);
  s.b=as<double>(rparam["b"]);
  s.a0=as<double>(rparam["a0"]);
  s.b0=as<double>(rparam["b0"]);
  s.c=as<double>(rparam["c"]);
  s.c0=as<double>(rparam["c0"]);
  s.d=as<double>(rparam["d"]);
  s.d0=as<double>(rparam["d0"]);
  s.e=as<double>(rparam["e"]);
  s.f=as<double>(rparam["f"]);
  s.alpha=as<double>(rparam["alpha"]);
  s.beta=as<double>(rparam["beta"]);
  s.lambda_e=as<double>(rparam["lambda.e"]);
  s.lambda_x=as<double>(rparam["lambda.x"]);
  s.verbose=as<bool>(rparam["verbose"]);
  s.lasso=as<bool>(rparam["lasso"]);
  s.sample_noise_variance=as<bool>(rparam["sample.noise.variance"]);
  s.sample_mixing_matrix_hypers=as<bool>(rparam["sample.mixing.matrix.hypers"]);
  s.sample_beta=as<bool>(rparam["sample.beta"]);
  s.sample_alpha=as<bool>(rparam["sample.alpha"]);
  s.sample_Bprior=as<bool>(rparam["sample.Bprior"]); 
  s.isotropic_noise=as<bool>(rparam["isotropic.noise"]);
  s.perfactor_variance=as<bool>(rparam["perfactor.variance"]);
  s.learn_scale=as<bool>(rparam["learn.scale"]);
  s.hierarchical_noise_prior=as<bool>(rparam["hierarchical.noise.prior"]);
  s.predictive_performance=as<bool>(rparam["predictive.performance"]); 
  s.log_file=as<string>(rparam["log.file"]);
  s.samples_dir=as<string>(rparam["samples.dir"]); 
  s.debug=as<bool>(rparam["debug"]);
  s.store_samples=as<bool>(rparam["store.samples"]); 
  s.thinout=as<int>(rparam["thin.out"]); 
  s.burnin=as<int>(rparam["burn.in"]);
  s.initial_size=as<int>(rparam["initial_size"]);
  setup_mrf(rparam["mrf.connections"],rparam["mrf.strengths"],s); 
  setup_feature_mrf(rparam["feature.mrf.connections"],rparam["feature.mrf.strengths"],s); 
								  
  return s; 
}

bool double_to_bool(double x){ return x!=0.0; }

List params_to_list(params<double> &p){
  SparseMatrix<double> sparseB=p.B.sparseView();
  SparseMatrix<double> sparseG=p.G.sparseView();
  return List::create(Named("factor.loadings") = sparseG,
		      Named("factors") = p.X,
		      Named("B")=sparseB, //p.B,
		      Named("alpha") = p.alpha,
		      Named("beta") = p.beta,
		      Named("Bprior") = p.Bprior,
		      Named("lambda.x") = p.lambda_x,
		      Named("lambda.e") = p.lambda_e);
}


SEXP params_to_sexp(params<double> &p){
  return params_to_list(p); 
}

params<double> sexp_to_params(SEXP sexp,settings<double> &s){
  params<double> p(s);
  Rcpp::List rparam(sexp);
  Map<MatrixXd> G = as<Map<MatrixXd> >(rparam["factor.loadings"]);
  p.G = G; 
  p.Z = G.unaryExpr(pointer_to_unary_function<double,bool>(double_to_bool)).cast<double>(); 
  Map<MatrixXd> B = as<Map<MatrixXd> >(rparam["B"]);
  p.B = B; 
  p.X = as<Map<MatrixXd> >(rparam["factors"]); 
  p.alpha = as<double>(rparam["alpha"]); 
  p.beta = as<double>(rparam["beta"]); 
  p.lambda_x = as<Map<VectorXd> >(rparam["lambda.x"]); 
  p.lambda_e = as<Map<VectorXd> >(rparam["lambda.e"]);
  p.Bprior = as<double>(rparam["Bprior"]); 
  return p; 
}

SEXP run_helper(settings<double> &s, SEXP init_param)
{
  params<double> p = sexp_to_params(init_param,s);
  params<double> best_p(s); 
  best_p.loglike=-1e10; 
  VectorXf loglike(s.iterations); 
  VectorXf test_loglike(s.iterations); 
  VectorXf test_rmse(s.iterations); 
  VectorXi Ks(s.iterations); 
  if (s.verbose)
    cout << "Running lacrosse: K " << p.Z.cols() << " debug " << s.debug << endl; 

  Function save=Environment::base_env()["save"]; 
  for (int i=0;i<s.iterations;i++){
    p.one_MCMC_iteration(s,i+1); 
    if (p.loglike>best_p.loglike)
      best_p=p;
    loglike(i)=p.loglike;
    test_loglike(i)=p.test_loglike;
    test_rmse(i)=p.test_rmse;
    Ks(i)=p.num_features(s); 
    if (s.store_samples && (i>s.burnin) && (i % s.thinout == 0)){
      char buff[200]; 
      sprintf(buff, "%s/sample%d.RData", s.samples_dir.c_str(), i); 
      string mystring=buff; 
      cout << "Saving sample to " << mystring << endl; 
      List l=params_to_list(p); 
      Environment::global_env().assign("param",l); 
      save("param",Named("file",mystring));
    }
    R_CheckUserInterrupt();
  }
  
  if (s.missing_values) {
    best_p.impute_Y(s); 
    return List::create(Named("best.params") = params_to_sexp(best_p),
			Named("final.params") = params_to_sexp(p),
			Named("num.features") = Ks,
			Named("loglikelihood") = loglike,
			Named("test.loglike") = test_loglike,
			Named("test.rmse") = test_rmse,
			Named("imputed.Y") = *s.Y);
  }
  else
    return List::create(Named("best.params") = params_to_sexp(best_p),
			Named("final.params") = params_to_sexp(p),
			Named("num.features") = Ks,
			Named("loglikelihood") = loglike);
}

RcppExport SEXP initial_param(SEXP N, SEXP D, SEXP P, SEXP settings_sexp){
  settings<double> s = parse_settings(settings_sexp); 
  s.N=as<int>(N); 
  s.D=as<int>(D); 
  s.P=as<int>(P);
  params<double> p(s); 
  return params_to_sexp(p); 
}

RcppExport SEXP run_lacrosse(SEXP Ys,SEXP Fs,SEXP settings_sexp,SEXP init_param){
  settings<double> s = parse_settings(settings_sexp); 
  MatrixXd Y(as<Map<MatrixXd> >(Ys)); 
  s.N=Y.cols();
  s.D=Y.rows();
  MatrixXd F(as<Map<MatrixXd> >(Fs)); 
  s.F=F; 
  s.P=F.rows(); 
  assert( Y.cols() == F.cols()); 
  s.setup_Y(Y,false);
  return run_helper(s,init_param); 
}

RcppExport SEXP run_lacrosse_missing(SEXP Ys,SEXP Fs,SEXP missing,SEXP settings_sexp,SEXP init_param){
  settings<double> s = parse_settings(settings_sexp); 
  MatrixXd Y(as<Map<MatrixXd> >(Ys)); 
  Map<MatrixXd> mv(as<Map<MatrixXd> >(missing));
  s.missing =  mv.unaryExpr(ptr_fun(double_to_bool));  
  MatrixXd F(as<Map<MatrixXd> >(Fs)); 
  s.F=F; 
  s.N=Y.cols();
  s.D=Y.rows();
  s.P=F.rows(); 
  assert( Y.cols() == F.cols()); 
  s.setup_Y(Y,true); 
  return run_helper(s,init_param); 
}

RcppExport SEXP test_Bprior(SEXP n,SEXP its){
  int N=as<int>(n);
  int iterations=as<int>(its); 
  VectorXd x(N); 
  NumericVector Bpriors(its); 
  double Bprior=1.0/(double)N; 
  settings<double> s; 
  for (int i=0; i<iterations; i++){
    for (int j=0; j<N; j++){
      double logrp=log(Bprior)-log(1.0-Bprior); 
      double prob_zis1=1.0/(1.0+exp(-logrp)); 
      assert( !isnan(prob_zis1));
      bool z=s.rand()<prob_zis1;
      x[j]=z ? 1.0 : 0.0; 
    }
    double sum=x.sum(); 
    Bprior=s.betarnd(sum+1.0/(double)N,(double)N-sum+1.0); 
    Bpriors[i]=Bprior; 
    cout << "Bprior " << Bprior << endl; 
  }
  return wrap(Bpriors); 
}
