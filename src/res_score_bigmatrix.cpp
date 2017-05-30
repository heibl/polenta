#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec add_msa_score(const arma::mat& ref, const arma::mat& com){

  int nr = ref.nrow();
  int nc = ref.ncol();
  int rpsc_n_row = nChoosek(nr,2);
  arma::vec rpsc(rpsc_n_row*nc);

  // Rcout << "intitializing done"<< endl;

  int count = 0;
  for(int col=0; col <= nc-1; col++){
    for(int row1 = 0; row1 <= nr-2; row1++){
      int resup = row1+1;
      for(int row2=resup; row2 <= nr-1; row2++){

        // reference residue pair
        int ref_rp1 = ref(row1, col);
        int ref_rp2 = ref(row2, col);
        // Rcout << "1: residue nr1 " << ref_rp1 << " and residue number2 " << ref_rp2 << std::endl;
        // if gap => NA

        if(ref_rp1%2 == 0 || ref_rp2%2 == 0){
          rpsc(count) = NA_REAL;
        }else{ //no gap
          // if the two aligned residues from the REF
          // are aligned in one column in the ALT => 1, else 0
          int upbase = which_true2(com(row1,_) == ref_rp1);
          int downbase = which_true2(com(row2,_) == ref_rp2);
          // Rcout << "2: base 1 " << upbase << " and base 2 " << upbase<< std::endl;

          if(upbase == downbase){
            rpsc(count) = 1;
          }else{
            rpsc(count) = 0;
          }
        } // row2
        count++;
        // Rcout << count << endl;
      } // row1
    } // col
  }
  return(rpsc);
}
