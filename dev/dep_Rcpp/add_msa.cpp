// Calculate residue pair score
// Franz-Sebastian Krah
// [2017-08-16]
#include <Rcpp.h>
#include <cstdlib>
// #include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
// using namespace arma;
// [[Rcpp::export]]
NumericMatrix add_msa(NumericMatrix ref_col2res, NumericMatrix alt_col2res,
  NumericMatrix alt_res2col, NumericMatrix hit_mat){

  // Rcout << "initializing done"<< endl;
  int ncol = ref_col2res.ncol();
  int nrow = ref_col2res.nrow();

  int count = 0;
  for(int col=0; col <= ncol-1; col++){

    for(int row1 = 0; row1 <= nrow-2; row1++){

      int res_Cmatrix = ref_col2res(row1,col);

      int res;
      if(res_Cmatrix%2  == 1){
        res=0.5*(res_Cmatrix+1)-1;
      }
      else{ res = -1;}

      int alt_col;
      if(res > -1){
        alt_col = alt_res2col(row1, res);
        // Rcout << "alt_col" << alt_col << endl;

        int resup = row1+1;
        for(int row2 = resup; row2 <= nrow-1; row2++){

          if(hit_mat(count,3)>-1 && (ref_col2res(row2,col) == alt_col2res(row2, alt_col))){
            hit_mat(count,3)++;
          }
          count++;
        }
      }
      else{
        int resup = row1+1;
        for(int row2 = resup; row2 <= nrow-1; row2++){
          count++;
        }
      }
    }
  }
  // colnames(respairhit) = CharacterVector::create("col", "row1", "row2","hit");
  return(hit_mat);
}

