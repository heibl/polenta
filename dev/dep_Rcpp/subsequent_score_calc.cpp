#include <Rcpp.h>
#include <cstdlib>
using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix RPS(NumericMatrix scores){

  int ncol = scores.ncol();
  int nrow = scores.nrow();

  NumericVector col = scores(_,1);
  NumericVector ucol = unique(col);
  NumericVector row1 = scores(_,2);
  NumericVector urow1 = unique(row1);
  NumericMatrix rps(ucol.size()*urow1.size(),3);

  NumericVector col_v = scores( _, 0);
  NumericVector row1_v = scores( _, 1);
  NumericVector row2_v = scores( _, 2);
  NumericVector scr = scores( _, 3);




  double hit = 0;
  for(int col; col< ucol.size(); col++){
    int count = 1;
    for(int row; row < urow1.size(); row++){

      hit = scr[col_v > 2 ];
      hit = cumsum(hit);
      hit = hit/count;
      count++;
    }
  }
}


