#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>


// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;



struct add_msa_par : public Worker
{
  // source matrix
  const RMatrix<double> ref_col2res;
  const RMatrix<double> ref_col2res_odd;
  const RMatrix<double> alt_col2res;
  const RMatrix<double> alt_res2col;
  RMatrix<double> hit_mat;

  // // // return matrix
  // RMatrix<double> res_mat;

  // initialize with source and destination
  add_msa_par(const NumericMatrix ref_col2res,
    const NumericMatrix ref_col2res_odd,
    const NumericMatrix alt_col2res,
    const NumericMatrix alt_res2col,
    NumericMatrix hit_mat)
    : ref_col2res(ref_col2res),
      ref_col2res_odd(ref_col2res_odd),
      alt_col2res(alt_col2res),
      alt_res2col(alt_res2col),
      hit_mat(hit_mat){}

  // calc
  void operator()(std::size_t begin, std::size_t end) {

    // double ncol = ref_col2res.ncol();
    // double nrow = ref_col2res.nrow();
    // RMatrix<double>::Column hit_v = hit_mat.column(4);

    for(std::size_t col=begin; col < end; col++){

      double count = 0;
      for(std::size_t row1 = 0; row1 <= ref_col2res.nrow()-2; row1++){

        // double res_Cmatrix = ref_col2res(row1,col);

        double res;
        // in the Cmatrix (col2res) the value for each non gap residue is 2*res+1
        // => this trasformes it back to the residue number
        // We need this number to look up the position of the residue in the
        // alternative MSA (by using res2col)
        if(ref_col2res_odd(row1,col) == 1){
          res=0.5*(ref_col2res(row1,col)+1)-1;
        }
        else{ res = -1;}
        // double alt_col = alt_res2col(row1, res);

        double raise = row1+1;
        for(std::size_t row2=raise; row2 <= ref_col2res.nrow()-1; row2++){

          // if the residue pair in the REF MSA is not a gab
          // and the second residue associated with the first residue is equal
          // => hit and rise the hit matrix (res_pair_hit) by one for this residue pair
          if(res > -1){
          if(hit_mat(count,col)>-1 && (ref_col2res(row2,col) == alt_col2res(row2, alt_res2col(row1, res)))){
            hit_mat(count, col)++;
            }
          }
          count++;
        }
      }
    }
  }
};


// [[Rcpp::export]]
NumericMatrix add_msa_parallel(NumericMatrix ref_col2res, NumericMatrix alt_col2res,
  NumericMatrix alt_res2col, NumericMatrix hit_mat, NumericMatrix ref_col2res_odd) {

  // // allocate the output matrix
  // NumericMatrix res_mat(hit_mat.nrow(), hit_mat.ncol());

  // functor (pass input and output matrixes)
  add_msa_par addmsa(ref_col2res, alt_col2res,
    alt_res2col, ref_col2res_odd, hit_mat);

  // call parallelFor to do the work
  parallelFor(0, ref_col2res.ncol(), addmsa);

  // return the output matrix
  return hit_mat;
}
