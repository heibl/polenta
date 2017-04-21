// cmatrix recoding
// Franz-Sebastian Krah
// [2017-04.17]
#include <Rcpp.h>
#include <cstdlib>
// #include <iostream>
using namespace std;
using namespace Rcpp;

// For more info about the Cmatrix recoding see page 3 on the supplemantary material of:
//   Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot:
//   Bayesian alignment and phylogenetic footprinting with MCMC,
//   BMC Evolutionary Biology, 9, 217.
// http://www.biomedcentral.com/content/supplementary/1471-2148-9-217-s1.pdf
//   */


// this is code if a CharacterMatrix is provided
// this code, however, does not work because Rcpp
// seems to have problems with if(X == "-")
// // [[Rcpp::export]]
// NumericMatrix cmatrixC(CharacterMatrix msa){
//   int n_row = msa.nrow();
//   int n_col = msa.ncol();
//   NumericMatrix msa_recode(n_row, n_col);
//
//   for(int row=0; row<=n_row; row++){
//     int LastNonGap = -1;
//     int res = 0;
//
//     for(int col=0; col<=n_col; col++){
//       if (msa(row,col) == '-'){
//         msa_recode(row,col) = LastNonGap+1;
//       } else {
//         msa_recode(row,col) = (2*res)+1;
//         LastNonGap = (2*res)+1;
//
//         res++;
//       }
//     }
//   }
//   return (msa_recode);
// }

// this is code if a NumericMatrix is provided (0 = gab, 1 = base)
// works fine, much faster than R version
// [[Rcpp::export]]
NumericMatrix Cmatrix(NumericMatrix msa){
  int n_row = msa.nrow();
  int n_col = msa.ncol();
  NumericMatrix msa_recode(n_row, n_col);

  for(int row=0; row<=n_row; row++){
    int LastNonGap = -1;
    int res = 0;

    for(int col=0; col<=n_col; col++){
      if (msa(row,col) == 0){
        msa_recode(row,col) = LastNonGap+1;
      } else {
        msa_recode(row,col) = (2*res)+1;
        LastNonGap = (2*res)+1;

        res++;
      }
    }
  }
  return (msa_recode);
}
