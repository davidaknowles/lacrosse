#pragma once
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
// #include <Eigen/Dense>
#include <RcppEigen.h>
using namespace boost; 
using namespace std;
using namespace Eigen;  

template<class real_type>
class settings {

public:
  // need this to be public to pass to random_matrix
  boost::random::mt19937 gen;

  // see constructor for hyperparameter descriptions
  real_type alpha,beta,a,b,a0,b0,c,d,c0,d0,e,f,lambda_e,lambda_x;
  bool verbose,store_samples,sample_noise_variance,sample_mixing_matrix_hypers,sample_alpha,sample_beta,isotropic_noise,perfactor_variance,learn_scale,hierarchical_noise_prior,missing_values,predictive_performance, truncate, sample_Bprior; 
  int iterations,thinout,burnin,N,D,P,initial_size;
  Matrix<bool,Dynamic,Dynamic> missing; 
  bool debug;
  bool lasso; 

  // this is a pointer so it can either be the original matrix or, in the case of
  // missing values, an imputed version
  Matrix<real_type,Dynamic,Dynamic> *Y, *Y_original; 
  Matrix<real_type,Dynamic,Dynamic> F;

  bool use_mrf, use_feature_mrf; 
  vector<vector<int> > mrf_connections, feature_mrf_connections; 
  vector<vector<real_type> > mrf_strength, feature_mrf_strength;

  string log_file,samples_dir; 

  void setup_Y(MatrixXd &Y_o,bool mv){
    missing_values=mv; 
    Y_original = &Y_o; 
    Y= missing_values ? new MatrixXd(Y_o) : &Y_o; 
  }

  // sample uniform between fMin and fMax
  real_type rand(real_type fMin = 0.0, real_type fMax = 1.0)
  {
    boost::uniform_real<> dist(fMin, fMax);
    return dist(gen);
  }

  real_type randn(real_type mu=0.0, real_type sd=1.0)
  {
    boost::normal_distribution<> dist(mu,sd); 
    return dist(gen); 
  }

  int rpoiss(real_type rate)
  {
    boost::poisson_distribution<> dist(rate); 
    return dist(gen); 
  }

  real_type gamrnd(real_type shape, real_type rate)
  {
    boost::gamma_distribution<> dist(shape, 1.0/rate) ;
    return dist(gen); 
  }

  real_type betarnd(real_type a, real_type b)
  {
    double x=gamrnd(a,1.0); 
    double y=gamrnd(b,1.0); 
    return x/(x+y); 
  }

  settings(){
    log_file.clear();
    use_mrf=false;
    debug=false;
    missing_values=false;
    predictive_performance=false;
    initial_size=0;
    a=1; // noise variance hyperparameters
    b=1;
    a0=1; // hyperparameters on b
    b0=1;
    d=1; // factor variance hyperparameters
    c0=1;
    d0=1;
    e=1; // IBP alpha hyperparamters
    f=1;
    alpha=1; // IBP hyperparameters
    beta=1;
    burnin=500;
    thinout=5;
    iterations=1000;
    verbose=false;
    lasso=false; 
    store_samples=false; 
    sample_noise_variance=true; 
    lambda_e=10; // noise precision
    sample_mixing_matrix_hypers=true; 
    lambda_x=1; // factor precision
    sample_beta=false;
    sample_Bprior=false; 
    // Interesting settings
    sample_alpha=true; // whether to sample the IBP concentration parameter
    isotropic_noise=true; // enforce isotropic noise?
    perfactor_variance=true; // Per factor covariances? With c=1 corresponds to ARD prior
    learn_scale=true; // Hierarchical prior on the factor covariances? Only valid if perfactor is true
    c=1; // 1 corresponds to an ARD prior. 2 corresponds to not. 
    hierarchical_noise_prior=false; // share power across noise dimensions? 
  }

};
