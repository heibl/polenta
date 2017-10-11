#include <Rcpp.h>
#include <cstdlib>
using namespace std;
using namespace Rcpp;

// // [[Rcpp::export]]
// NumericVector arrayC(NumericVector input, IntegerVector dim) {
//   input.attr("dim") = dim;
//   return input;
// }

// [[Rcpp::export]]
NumericVector rep_C(NumericVector x, NumericVector y) {
  int n = y.size();
  NumericVector myvector(sum(y));
  int ind = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < y(i); ++j) {
      myvector(ind) = x[i];
      ind = ind + 1;
    }
  }
  return myvector;
}

// [[Rcpp::export]]
NumericVector rep_example( NumericVector x, int n){
  NumericVector res = rep_each(x, n) ;
  return res ;
}

// [[Rcpp::export]]
NumericVector res_pair_hit_3d(NumericMatrix mat){

  int n_row = mat.nrow();
  int n_col = mat.ncol();

  NumericVector dim = NumericVector::create(n_row, n_row, n_col);
  NumericVector x = 1;
  int n = n_row*n_row*n_col;
  NumericVector respairhit = rep_each(x, n) ;
  respairhit.attr("dim") = dim;
  Rcout << respairhit(1,2) << endl;
  return(respairhit);
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cube array2cube( SEXP myArray ) {

  Rcpp::NumericVector vecArray(myArray);
  Rcpp::IntegerVector arrayDims = vecArray.attr("dim");

  arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

  return(cubeArray);

}

//   // int count = 0;
//   for(int col=0; col <= n_col-1; col++){
//
//     for(int row1 = 0; row1 <= n_row-2; row1++){
//       int resup = row1+1;
//       for(int row2 = resup; row2 <= n_row-1; row2++){
//
//         int res1 = col2res(row1, col);
//         int res2 = col2res(row2, col);
//
//         int oddornot1 = res1 % 2;
//         int oddornot2 = res2 % 2;
//
//
//         // if both residues are bases => 0
//         // if one of the basis is a gab => -1
//         int hit;
//         if ( oddornot1 == 1 && oddornot2 == 1 ) {
//           respairhit(1,2) = 0;
//         }else{
//           respairhit(row1, row2, col) = -1;
//         }
//       }
//     }
//   }
//   // colnames(respairhit) = CharacterVector::create("col", "row1", "row2","score");
//   return(respairhit);
// }
