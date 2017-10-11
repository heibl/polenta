// Franz-Sebastian Krah [2017-08-16]
#include <Rcpp.h>
#include <cstdlib>
using namespace std;
using namespace Rcpp;

//[[Rcpp::export]]
NumericMatrix msa_recode(NumericMatrix score){

  NumericVector cols = unique(score(,0));
  int n = cols.size();

  for(int j=0; j<n; j++){

    df = gr_scr(gr_scr(,1)==j, );
    rows <- unique(c(unique(df$row1), unique(df$row2)))

    for(int col=0; col < n_col; col++){
      if (msa(row,col) == 0){
        col2res(row,col) = lastnongap+1;
      } else {
        // 'res2col': indicies of the characters in the alignment
        res2col(row,res) = col;
        // 'col2res': characters are represented by odd numbers
        //            and gaps by even numbers
        col2res(row,col) = (2*res)+1;
        lastnongap = (2*res)+1;
        res++;
        // Rcout << count << endl;
      }
    }
  }
  return List::create(Rcpp::Named("col2res") = col2res,
    Rcpp::Named("res2col") = res2col);
}
