// Calculate residue pair score
// Franz-Sebastian Krah
// [2017-04-22]
#include <Rcpp.h>
#include <cstdlib>
using namespace std;
using namespace Rcpp;

//   For more info about the Cmatrix recoding see page 3 on the supplemantary material of:
//   Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot:
//   Bayesian alignment and phylogenetic footprinting with MCMC,
//   BMC Evolutionary Biology, 9, 217.
//   http://www.biomedcentral.com/content/supplementary/1471-2148-9-217-s1.pdf
// [[Rcpp::export]]
List msa_recode(NumericMatrix msa){
  // msa = transpose(msa);
  // Rcout << "transpose done"<< endl;
  int n_row = msa.nrow();
  int n_col = msa.ncol();
  // Rcout << "numbers  done" << n_row << " " << n_col << endl;

  NumericMatrix col2res(n_row, n_col);
  NumericMatrix res2col(n_row, n_col);

  // Rcout << "initializing done"<< endl;

  for(int row=0; row<n_row; row++){
    int lastnongap = -1;
    int res = 0;

    for(int col=0; col < n_col; col++){
      if (msa(row,col) == 0){
        col2res(row,col) = lastnongap+1;
      } else {
        res2col(row,res) = col;
        col2res(row,col) = (2*res)+1;
        lastnongap = (2*res)+1;
        res++;
        // Rcout << count << endl;
      }
    }
  }
  // msa_recode = transpose(msa_recode);
  // Rcout << "transpose done"<< endl;
  // return (col2res);
  return List::create(Rcpp::Named("col2res") = col2res,
    Rcpp::Named("res2col") = res2col);
}
