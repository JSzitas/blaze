#ifndef FORECAST_HEADER
#define FORECAST_HEADER

#include "vector"

struct forecast_result {
  forecast_result(int h) {
    this->forecast = std::vector<double>(h);
    this->se = std::vector<double>(h);
  }
  void add(int i, double forecast, double se) {
    this->forecast[i] = forecast;
    this->se[i] = se;
  }
  std::vector<double> forecast;
  std::vector<double> se;
};

/* Forecasts based on state space representation of ARIMA via
 * the kalman filter.
 */
forecast_result kalman_forecast( int n_ahead,
                                 std::vector<double> & Z,
                                 std::vector<double> & a,
                                 std::vector<double> & P,
                                 std::vector<double> & T,
                                 std::vector<double> & V,
                                 double h,
                                 bool update = false ){
  int p = a.size();
  std::vector<double> anew(p);
  std::vector<double> Pnew(p*p);
  std::vector<double> mm(p*p);

  forecast_result res(n_ahead);

  double fc, tmp;
  for (int l = 0; l < n_ahead; l++) {
    fc = 0.0;
    for (int i = 0; i < p; i++) {
      tmp = 0.0;
      for (int k = 0; k < p; k++) {
        tmp += T[i + p * k] * a[k];
      }
      anew[i] = tmp;
      fc += tmp * Z[i];
    }
    for (int i = 0; i < p; i++) {
      a[i] = anew[i];
    }
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        tmp = 0.0;
        for (int k = 0; k < p; k++) {
          tmp += T[i + p * k] * P[k + p * j];
        }
        mm[i + p * j] = tmp;
      }
    }
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        tmp = V[i + p * j];
        for (int k = 0; k < p; k++) {
          tmp += mm[i + p * k] * T[j + p * k];
        }
        Pnew[i + p * j] = tmp;
      }
      tmp = h;
    }
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        P[i + j * p] = Pnew[i + j * p];
        tmp += Z[i] * Z[j] * P[i + j * p];
      }
      res.add(l, fc, tmp );
    }
  }
  return res;
}

#endif
