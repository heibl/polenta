// Calculate residue pair score
// Franz-Sebastian Krah
// [2017-04-22]
#include <Rcpp.h>
#include <cstdlib>
using namespace std;
using namespace Rcpp;

// For more info about the Cmatrix recoding see page 3 on the supplemantary material of:
//   Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot:
//   Bayesian alignment and phylogenetic footprinting with MCMC,
//   BMC Evolutionary Biology, 9, 217.
// http://www.biomedcentral.com/content/supplementary/1471-2148-9-217-s1.pdf
//   */
// [[Rcpp::export]]
NumericMatrix Cmatrix(NumericMatrix msa){
  msa = transpose(msa);
  int n_row = msa.nrow();
  int n_col = msa.ncol();
  NumericMatrix msa_recode(n_row, n_col);
  Rcout << "initializing done"<< endl;

  for(int col=0; col<=n_col; col++){
    int lastnongap = -1;
    int count = 0;

    for(int row=0; row<=n_row; row++){
      if (msa(row,col) == 0){
        msa_recode(row,col) = lastnongap+1;
      } else {
        msa_recode(row,col) = (2*count)+1;
        lastnongap = (2*count)+1;
        count++;
        Rcout << count << endl;
      }
    }
  }
  msa_recode = transpose(msa_recode);
  // this is experimental: R quits unexpectedly, could have to do with memory issues
  // this might be a way to fix it
  // Rcpp::Environment G = Rcpp::Environment::global_env();
  // Rcpp::Function gc = G["gc"];
  return (msa_recode);
}
//[[Rcpp::export]]
int nChoosek( int n, int k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}
// [[Rcpp::export]]
int which_true(LogicalVector x) {
  int counter = 0;
  int out;
  for(int i = 0; i < x.size(); i++) {
    counter++;
    if(x[i] == TRUE) {
      out = counter;
    }
  }
  return out;
}
// [[Rcpp::export]]
NumericMatrix add_msa(NumericMatrix ref, NumericMatrix com){

  int ref_n_row = ref.nrow();
  int ref_n_col = ref.ncol();
  int rpsc_n_row = nChoosek(ref_n_row,2);
  NumericMatrix rpsc(rpsc_n_row*ref_n_col,5);
  NumericVector ref_rp;

  Rcout << "intitializing done"<< endl;

  int count = 0;
  for(int col=0; col <= ref_n_col-1; col++){
    for(int row1 = 0; row1 <= ref_n_row-2; row1++){
      int resup = row1+1;
      for(int row2=resup; row2 <= ref_n_row-1; row2++){


        // reference residue pair
        int ref_rp1 = ref(row1, col);
        int ref_rp2 = ref(row2, col);
        // Rcout << "1: residue nr1 " << ref_rp1 << " and residue number2 " << ref_rp2 << std::endl;
        // if gap => NA
        if(ref_rp1%2 == 0 || ref_rp2%2 == 0){

          rpsc(count,0) = col;
          rpsc(count,1) = row1;
          rpsc(count,2) = row2;
          rpsc(count,3) = NA_REAL;
          rpsc(count,4) = count;

        }else{ //no gap
          // if the two aligned residues from the REF
          // are aligned in one column in the ALT => 1, else 0
          int upbase = which_true(com(row1,_) == ref_rp1);
          int downbase = which_true(com(row2,_) == ref_rp2);
          // Rcout << "2: base 1 " << upbase << " and base 2 " << upbase<< std::endl;

          if(upbase == downbase){
            rpsc(count,0) = col;
            rpsc(count,1) = row1;
            rpsc(count,2) = row2;
            rpsc(count,3) = 1;
            rpsc(count,4) = count;
          }else{
            rpsc(count,0) = col;
            rpsc(count,1) = row1;
            rpsc(count,2) = row2;
            rpsc(count,3) = NA_REAL;
            // rpsc(count,4) = 0;}

          }
        } // row2
        count++;
        Rcout << count << endl;
      } // row1
    } // col

    //   colnames(rps) =  c("col", "row1", "row2", "score")
  }
  return(rpsc);
}
