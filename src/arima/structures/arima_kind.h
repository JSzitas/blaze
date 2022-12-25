#ifndef ARIMA_KIND
#define ARIMA_KIND

// defines the arima structure
struct arima_kind {
  arima_kind(){};
  arima_kind(int p, int d, int q, int P, int D, int Q, int s_period) {
    this->arma_p = p;
    this->diff_d = d;
    this->arma_q = q;
    this->sarma_P = P;
    this->seas_diff_D = D;
    this->sarma_Q = Q;
    this->s_period = s_period;
  }
  arima_kind(std::vector<int> arima_def) {
    this->arma_p = arima_def[0];
    this->diff_d = arima_def[1];
    this->arma_q = arima_def[2];
    this->sarma_P = arima_def[3];
    this->seas_diff_D = arima_def[4];
    this->sarma_Q = arima_def[5];
    this->s_period = arima_def[6];
  }
  const int p() const { return this->arma_p; }
  const int d() const { return this->diff_d; }
  const int q() const { return this->arma_q; }
  const int P() const { return this->sarma_P; }
  const int D() const { return this->seas_diff_D; }
  const int Q() const { return this->sarma_Q; }
  const int period() const { return this->s_period; }

private:
  int arma_p, diff_d, arma_q, sarma_P, seas_diff_D, sarma_Q, s_period;
};

#endif
