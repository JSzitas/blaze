#ifndef SIMPLE_OLS
#define SIMPLE_OLS

namespace small_lm {
// assume both have equal length
template <typename T> T inner_product( const std::vector<T> &x,
                                       const std::vector<T> &y) {
  T result = 0;
  for( int i = 0; i < x.size(); i++ ) {
    result += (x[i] * y[i]);
  }
  return result;
}

template <typename T> struct Matrix {
  Matrix( std::vector<T> data, int i, int j ) {
    this->X = data;
    this->n_row = i;
    this->n_col = j;
    this->size = i*j;
  };
  std::vector<T> col( int j ){
    // iterate over data and find the corresponding indices which match with j
    std::vector<T> result( this->n_row );
    for(int i = 0; i < this->n_row; i++ ) {
      result[i] = this->X[i + (j*this->n_row) ];
    }
    return result;
  }
  std::vector<T> row( int i) {
    // iterate over data and find the corresponding indices which match with j
    // this is actually not a contingent block - this skips nrow rows each time
    std::vector<T> result(this->n_col);
    for( int j=0; j < this->n_col; j++ ) {
      // this does index i + however many offsets we need to get to the next row
      result[j] = this->X[ i + (j*this->n_row) ];
    }
    return result;
  }
  T elem( int ind ) {
    return this->X[ind];
  }
  void t() {
    // figure out how to do transposition nicely
    // this->X
  }
  int nrow(){ return this->n_row; }
  int ncol(){ return this->n_col; }
  int size(){ return this->size(); }
private:
  std::vector<T> X;
  int n_col;
  int n_row;
};

template <typename T> Matrix<T> operator *( Matrix<T> X, Matrix<T> Y ) {

  std::vector<T> temp(X.nrow() * Y.ncol());
  int p = 0;
  for( int i=0; i < X.nrow(); i++ ) {
    for( int j=0; j < Y.ncol(); j++ ) {
      temp[p] = inner_product( X.row(i), Y.col(j));
      p++;
    }
  }
  Matrix<T> result( temp,
                    X.nrow(),
                    Y.ncol());
  return result;
}

template <typename T> Matrix<T> cholesky( Matrix<T> A ) {
  // int p = 0; // index of current element that we are updating

  // int rows = A.nrow();

  // result[0] = sqrt( A.elem(0) );

  // for (int j = 0; j < A.ncol(); j++) {
  //   // iterate over columns - this has a funky offsetting rule
  //   // in the first column, there is no offset, and then the offset deepens by 1
  //
  //   // we compute the first row of a particular column like this
  //   result[p] = sqrt( A.elem(j * rows) );
  //
  //   for( int k=j+1; k < rows; k++ ) {
  //     // this iterates over rows - but only the rows left in the current column -
  //     // since this is lower triangular, we know that k is at most j
  //
  //
  //
  //   }
  //
  //
  //
  // }
    // sum = 0;
    // for (int k = 0; k < j; k++) {
      // sum += L[j][k] * L[j][k];
    // }
    // L[j][j] = sqrt(A[j][j] - sum);

  //   for (i = j + 1; i < A.ncol(); i++) {
  //     sum = 0;
  //     for (k = 0; k < j; k++) {
  //       sum += L[i][k] * L[j][k];
  //     }
  //     L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
  //   }
  // }
}



// quick algorithm to compute ordinary least squares solution Ax = b
template <typename T> std::vector<T> lm( std::vector<T> & y,
                                         Matrix<T> &X) {




}

//
// template <typename T> std::vector<T> projection( const std::vector<T> & u,
//                                                  const std::vector<T> & a) {
//   T upper = inner_product(u,a);
//   T lower = inner_product(u,u);
//
//   T ratio = upper/lower;
//   std::vector<T> result(u.size());
//   int i = 0;
//   for( auto &item:u ) {
//     result[i] = item * ratio;
//     i++;
//   }
//
//   return result;
// }
//
//
// template <typename T> std::vector<T> QR_gram_schmidt( Matrix<T> &X ) {
//
//
// }

}

#endif
