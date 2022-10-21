#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
void tseries_arma ( std::vector<double> &x,
                    std::vector<double> &u,
                    std::vector<double> &a,
                    std::vector<int> ar,
                    std::vector<int> ma,
                    int arl,
                    int mal,
                    int max,
                    int n,
                    int intercept)
  /* compute conditional sum of squares */
{
  int i, j;
  double sum;

  for (i=max; i<n; i++) {
    if( intercept) sum = a[mal + arl];
    else sum = 0.0;
    for (j=0; j<arl; j++)
      sum += a[j]*x[i-ar[j]];
    for (j=0; j<mal; j++)
      sum += a[j+arl] * u[i-ma[j]];
    u[i]=x[i]-sum;
  }
}
