<<<<<<< HEAD:src/res_pair_score.cpp
// // Calculate residue pair score
// // Franz-Sebastian Krah
// // [2017-04-22]
// #include <Rcpp.h>
// #include <cstdlib>
// using namespace std;
// using namespace Rcpp;
//
// //   For more info about the Cmatrix recoding see page 3 on the supplemantary material of:
// //   Satija R, Novak A., Mikls I., Lyngs R., and Hein J. (2009) BigFoot:
// //   Bayesian alignment and phylogenetic footprinting with MCMC,
// //   BMC Evolutionary Biology, 9, 217.
// //   http://www.biomedcentral.com/content/supplementary/1471-2148-9-217-s1.pdf
// // [[Rcpp::export]]
// NumericMatrix Cmatrix(NumericMatrix msa){
//   msa = transpose(msa);
//   // Rcout << "transpose done"<< endl;
//   int n_row = msa.nrow();
//   int n_col = msa.ncol();
//   // Rcout << "numbers  done" << n_row << " " << n_col << endl;
//
//   NumericMatrix msa_recode(n_row, n_col);
//   // Rcout << "initializing done"<< endl;
//
//   for(int col=0; col<=n_col-1; col++){
//     int lastnongap = -1;
//     int count = 0;
//
//     for(int row=0; row <= n_row-1; row++){
//       if (msa(row,col) == 0){
//         msa_recode(row,col) = lastnongap+1;
//       } else {
//         msa_recode(row,col) = (2*count)+1;
//         lastnongap = (2*count)+1;
//         count++;
//         // Rcout << count << endl;
//       }
//     }
//   }
//   msa_recode = transpose(msa_recode);
//   // Rcout << "transpose done"<< endl;
//   return (msa_recode);
// }
// //[[Rcpp::export]]
// int nChoosek( int n, int k )
// {
//   if (k > n) return 0;
//   if (k * 2 > n) k = n-k;
//   if (k == 0) return 1;
//
//   int result = n;
//   for( int i = 2; i <= k; ++i ) {
//     result *= (n-i+1);
//     result /= i;
//   }
//   return result;
// }
=======
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
NumericMatrix Cmatrix(NumericMatrix msa){
  msa = transpose(msa);
  // Rcout << "transpose done"<< endl;
  int n_row = msa.nrow();
  int n_col = msa.ncol();
  // Rcout << "numbers  done" << n_row << " " << n_col << endl;

  NumericMatrix msa_recode(n_row, n_col);
  // Rcout << "initializing done"<< endl;

  for(int col=0; col<=n_col-1; col++){
    int lastnongap = -1;
    int count = 0;

    for(int row=0; row <= n_row-1; row++){
      if (msa(row,col) == 0){
        msa_recode(row,col) = lastnongap+1;
      } else {
        msa_recode(row,col) = (2*count)+1;
        lastnongap = (2*count)+1;
        count++;
        // Rcout << count << endl;
      }
    }
  }
  msa_recode = transpose(msa_recode);
  // Rcout << "transpose done"<< endl;
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


>>>>>>> 7254b7c08b135f842d2ba3b215ff9a5fa9f1db6a:dev/dep_Rcpp/res_pair_score.cpp
// // [[Rcpp::export]]
// int which_true2(LogicalVector x) {
//   int counter = 0;
//   int out;
//   for(int i = 0; i < x.size(); i++) {
//     counter++;
//     if(x[i] == TRUE) {
//       out = counter;
//     }
//   }
//   return out;
// }
//
// // // [[Rcpp::export]]
// // int which_true2(LogicalVector x) {
// //   int counter = 0;
// //   while (x[counter] != TRUE) {
// //     counter++ ;
// //     if(x[counter] == TRUE){
// //       break;
// //     }
// //   }
// //   return (counter+1);
// // }
// // This matrix version is very slow // [[Rcpp::export]]
// // NumericMatrix add_msa(NumericMatrix ref, NumericMatrix com){
// //
// //   int nr = ref.nrow();
// //   int nc = ref.ncol();
// //   int rpsc_n_row = nChoosek(nr,2);
// //   NumericMatrix rpsc(rpsc_n_row*nc,4);
// //
// //   // Rcout << "intitializing done"<< endl;
// //
// //   int count = 0;
// //   for(int col=0; col <= nc-1; col++){
// //     for(int row1 = 0; row1 <= nr-2; row1++){
// //       int resup = row1+1;
// //       for(int row2=resup; row2 <= nr-1; row2++){
// //
// //
// //         // reference residue pair
// //         int ref_rp1 = ref(row1, col);
// //         int ref_rp2 = ref(row2, col);
// //         // Rcout << "1: residue nr1 " << ref_rp1 << " and residue number2 " << ref_rp2 << std::endl;
// //         // if gap => NA
// //
// //         rpsc(count,0) = col+1;
// //         rpsc(count,1) = row1+1;
// //         rpsc(count,2) = row2+1;
// //
// //         if(ref_rp1%2 == 0 || ref_rp2%2 == 0){
// //           rpsc(count,3) = NA_REAL;
// //           // rpsc(count,4) = count+1;
// //
// //         }else{ //no gap
// //           // if the two aligned residues from the REF
// //           // are aligned in one column in the ALT => 1, else 0
// //           int upbase = which_true(com(row1,_) == ref_rp1);
// //           int downbase = which_true(com(row2,_) == ref_rp2);
// //           // Rcout << "2: base 1 " << upbase << " and base 2 " << upbase<< std::endl;
// //
// //           if(upbase == downbase){
// //             rpsc(count,3) = 1;
// //             // rpsc(count,4) = count+1;
// //           }else{
// //             rpsc(count,3) = 0;
// //             // rpsc(count,4) = count+1;
// //           }
// //         } // row2
// //         count++;
// //         // Rcout << count << endl;
// //       } // row1
// //     } // col
// //   }
// //   return(rpsc);
// // }
// // [[Rcpp::export]]
// NumericVector add_msa_score(NumericMatrix ref, NumericMatrix com){
//
//   int nr = ref.nrow();
//   int nc = ref.ncol();
//   int rpsc_n_row = nChoosek(nr,2);
//   NumericVector rpsc(rpsc_n_row*nc);
//
//   // Rcout << "intitializing done"<< endl;
//
//   int count = 0;
//   for(int col=0; col <= nc-1; col++){
//     for(int row1 = 0; row1 <= nr-2; row1++){
//       int resup = row1+1;
//       for(int row2=resup; row2 <= nr-1; row2++){
//
//         // reference residue pair
//         int ref_rp1 = ref(row1, col);
//         int ref_rp2 = ref(row2, col);
//         // Rcout << "1: residue nr1 " << ref_rp1 << " and residue number2 " << ref_rp2 << std::endl;
//         // if gap => NA
//
//         if(ref_rp1%2 == 0 || ref_rp2%2 == 0){
//           rpsc(count) = NA_REAL;
//         }else{ //no gap
//           // if the two aligned residues from the REF
//           // are aligned in one column in the ALT => 1, else 0
//           int upbase = which_true2(com(row1,_) == ref_rp1);
//           int downbase = which_true2(com(row2,_) == ref_rp2);
//           // Rcout << "2: base 1 " << upbase << " and base 2 " << upbase<< std::endl;
//
//           if(upbase == downbase){
//             rpsc(count) = 1;
//           }else{
//             rpsc(count) = 0;
//           }
//         } // row2
//         count++;
//         // Rcout << count << endl;
//       } // row1
//     } // col
//   }
//   return(rpsc);
// }
<<<<<<< HEAD:src/res_pair_score.cpp
//
// // [[Rcpp::export]]
// NumericMatrix rps_mat_maker(int nr, int nc){
//
//   int rpsc_n_row = nChoosek(nr, 2);
//   NumericMatrix rpsc(rpsc_n_row*nc, 3);
//
//   int count = 0;
//   for(int col=0; col <= nc-1; col++){
//     for(int row1 = 0; row1 <= nr-2; row1++){
//       int resup = row1+1;
//       for(int row2=resup; row2 <= nr-1; row2++){
//
//         rpsc(count,0) = col+1;
//         rpsc(count,1) = row1+1;
//         rpsc(count,2) = row2+1;
//
//         count++;
//       } // row1
//     } // col
//   }
//   return(rpsc);
// }
=======
// [[Rcpp::export]]
NumericVector add_msa_score(NumericMatrix ref, NumericMatrix com){

  int nr = ref.nrow();
  int nc = ref.ncol();
  int rpsc_n_row = nChoosek(nr,2);
  NumericVector rpsc(rpsc_n_row*nc);

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

// [[Rcpp::export]]
NumericMatrix rps_mat_maker(int nr, int nc){

  int rpsc_n_row = nChoosek(nr, 2);
  NumericMatrix rpsc(rpsc_n_row*nc, 3);

  int count = 0;
  for(int col=0; col <= nc-1; col++){
    for(int row1 = 0; row1 <= nr-2; row1++){
      int resup = row1+1;
      for(int row2=resup; row2 <= nr-1; row2++){

        rpsc(count,0) = col+1;
        rpsc(count,1) = row1+1;
        rpsc(count,2) = row2+1;

        count++;
      } // row1
    } // col
}
  return(rpsc);
}
>>>>>>> 7254b7c08b135f842d2ba3b215ff9a5fa9f1db6a:dev/dep_Rcpp/res_pair_score.cpp
// // [[Rcpp::export]]
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
