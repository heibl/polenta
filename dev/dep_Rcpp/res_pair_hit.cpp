// Calculate residue pair score
// Franz-Sebastian Krah
// [2017-04-22]
#include <Rcpp.h>
#include <cstdlib>
// #include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
// using namespace arma;

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
NumericMatrix res_pair_hit(NumericMatrix col2res){
  // msa = transpose(msa);
  // Rcout << "transpose done"<< endl;
  int n_row = col2res.nrow();
  int n_col = col2res.ncol();

  int nrows = nChoosek(n_row, 2)*n_col;
  NumericMatrix respairhit(nrows, 4);

  // Rcout << "initializing done"<< endl;

  int count = 0;
  for(int col=0; col <= n_col-1; col++){

    for(int row1 = 0; row1 <= n_row-2; row1++){
      int resup = row1+1;
      for(int row2 = resup; row2 <= n_row-1; row2++){

        int res1 = col2res(row1, col);
        int res2 = col2res(row2, col);

        int oddornot1 = res1 % 2;
        int oddornot2 = res2 % 2;

        respairhit(count,0) = col+1;
        respairhit(count,1) = row1+1;
        respairhit(count,2) = row2+1;

        int hit;
        if ( oddornot1 == 1 && oddornot2 == 1 ) {
          hit = 0;
        }else{
          hit = -1;
        }
        respairhit(count,3) = hit;
        // Rcout << "col:" << col << "row1:" << row << "row2:" << row2 << "hit:"<< hit << endl;
        count++;
      }
    }
  }
  colnames(respairhit) = CharacterVector::create("col", "row1", "row2","hit");
  return(respairhit);
}
