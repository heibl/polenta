// examples: http://gallery.rcpp.org/articles/parallel-distance-matrix/
// and http://gallery.rcpp.org/articles/parallel-matrix-transform/

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::depends(RcppParallel)]]
struct Cmatrix_par : public Worker
{
  // source matrix
  const RMatrix<double> mat;

  // destination matrix
  RMatrix<double> cmat;

  // initialize with source and destination
  Cmatrix_par(const NumericMatrix mat, NumericMatrix cmat)
    : mat(mat), cmat(cmat) {}

  // calc Cmatrix
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t col = begin; col < end; col++){
      double lastnongap = -1;
      double count = 0;

      for(std::size_t row=0; row <= mat.nrow() - 1; row++){
        if (mat(row,col) == 0){
          cmat(row,col) = lastnongap+1;
        } else {
          cmat(row,col) = (2*count)+1;
          lastnongap = (2*count)+1;
          count++;
          // Rcout << count << endl;
        }
      }
    }
  }
  // cmat = transpose(cmat);
};


// [[Rcpp::export]]
NumericMatrix Cmatrix_p(NumericMatrix mat) {

  // mat = transpose(mat);

  // allocate the output matrix
  NumericMatrix cmat(mat.nrow(), mat.ncol());

  // SquareRoot functor (pass input and output matrixes)
  Cmatrix_par cmatrix_par(mat, cmat);

  // call parallelFor to do the work
  parallelFor(0, mat.ncol(), cmatrix_par);

  // return the output matrix
  // cmat = transpose(cmat);
  return cmat;
}
