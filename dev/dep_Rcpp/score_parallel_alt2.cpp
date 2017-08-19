// Residue pair score function
// Franz Krah
// 06-24-2017
// Status: running, bug after first iteration
// examples: http://gallery.rcpp.org/articles/parallel-distance-matrix/
// and http://gallery.rcpp.org/articles/parallel-matrix-transform/


#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>


// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct add_msa_score : public Worker
{
  // source matrix
  const RMatrix<double> msa;
  const RMatrix<double> ref;

  // destination matrix
  RMatrix<double> rpsc;

  // initialize with source and destination
  // msa and ref are outputs from Cmatrix (odds are bases, evens are gaps)
  add_msa_score(const NumericMatrix msa, const NumericMatrix ref, NumericMatrix rpsc)
    : msa(msa), ref(ref), rpsc(rpsc) {}

  // calc Cmatrix
  void operator()(std::size_t begin, std::size_t end) {

    // double count = 0;
    double nr = ref.nrow();
    // column loop


    for(std::size_t col=begin; col < end; col++){

      // count of row number for res matrix
      RVector<int> count = 0;

      // first residue
      for(std::size_t row1 = 0; row1 <= nr-2; row1++){

        // second residue
        double raise = row1+1;
        for(std::size_t row2=raise; row2 <= nr-1; row2++){

          // reference residue pair
          int ref_rp1 = ref(row1, col);
          int ref_rp2 = ref(row2, col);

          int oddornot1 = ref_rp1 % 2;
          int oddornot2 = ref_rp2 % 2;

          // if one residue is gap (even number) => NA
          // if ( (ref_rp1 % 2 == 0) || (ref_rp2 % 2 == 0) ){
          double hit;

          if ( oddornot1 != 1 || oddornot2 != 1 ){

            hit = NA_REAL;

          } else { // no gap => hit or no hit

            // if the two aligned residues from the REF
            // are aligned in one column in the ALT (MSA) => 1, else 0

            // rows we will operate on
            RMatrix<double>::Row r1 = msa.row(row1);
            RMatrix<double>::Row r2 = msa.row(row2);

            // find positions of residues in the ALT MSA
            // *which* function was not recognized as function
            double hit1;
            for(std::size_t counter = 0; counter < msa.ncol(); counter++){
              if(r1[counter] == ref_rp1){
                hit1 = counter;
              }
            }
            // while (r1[counter] != ref_rp1) {
            //   counter++ ;
            // }
            // hit1 = counter;
            double hit2;
            for(std::size_t counter = 0; counter < msa.ncol(); counter++){
              if(r2[counter] == ref_rp2){
                hit2 = counter;
              }
            }

            // score
            if(hit1 == hit2){
              hit = 1;
            }
            if(hit1 != hit2){
              hit = 0;
            }
          }
          rpsc(count, col) = hit;
          count++;
        } // row2
      } // row1
    } // col
  }
};
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
// // [[Rcpp::export]]
////[[Franz; 06-24-2017]] would like to use this in the main function,
// however, this is not working; need template form?
// int which_true(NumericVector x, int y) {
//   int counter = 0;
//   while (x[counter] != y) {
//     counter++ ;
//     if(x[counter] == y){
//       break;
//     }
//   }
//   return (counter+1);
// }

// [[Rcpp::export]]
NumericVector add_msa_sc(NumericMatrix msa, NumericMatrix ref) {

  // allocate the output matrix
  NumericMatrix rpsc(nChoosek(ref.nrow(),2), ref.ncol());

  // SquareRoot functor (pass input and output matrixes)
  add_msa_score addmsa(msa, ref, rpsc);

  // call parallelFor to do the work
  parallelFor(0, ref.ncol(), addmsa);

  // return the output matrix
  return rpsc;
}
