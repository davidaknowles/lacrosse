#pragma once
#ifndef _lacrosse_H
#define _lacrosse_H

#include <RcppEigen.h>
#include "params.hpp"
#include "settings.hpp"
// [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
// [[Rcpp::export]]

RcppExport SEXP run_lacrosse(SEXP Ys,SEXP Fs,SEXP settings_sexp, SEXP ip);
RcppExport SEXP run_lacrosse_missing(SEXP Ys,SEXP Fs,SEXP missing,SEXP settings_sexp,SEXP ip);
RcppExport SEXP initial_param(SEXP N, SEXP D, SEXP P,SEXP settings_sexp);

#endif
