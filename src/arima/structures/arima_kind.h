#ifndef ARIMA_KIND
#define ARIMA_KIND

// defines the arima structure
struct arima_kind {
  arima_kind(){
    this->arma_p = 0;
    this->diff_d = 0;
    this->arma_q = 0;
    this->sarma_P = 0;
    this->seas_diff_D = 0;
    this->sarma_Q = 0;
    this->s_period = 0;
  };
  arima_kind(const size_t p,
             const size_t d,
             const size_t q,
             const size_t P,
             const size_t D,
             const size_t Q,
             const size_t s_period) : arma_p(p),
             diff_d(d), arma_q(q), sarma_P(P),
             seas_diff_D(D), sarma_Q(Q), s_period(s_period) {}
  arima_kind(const std::vector<size_t> &arima_def) {
    this->arma_p = arima_def[0];
    this->diff_d = arima_def[1];
    this->arma_q = arima_def[2];
    this->sarma_P = arima_def[3];
    this->seas_diff_D = arima_def[4];
    this->sarma_Q = arima_def[5];
    this->s_period = arima_def[6];
  }
  const size_t p() const { return this->arma_p; }
  const size_t d() const { return this->diff_d; }
  const size_t q() const { return this->arma_q; }
  const size_t P() const { return this->sarma_P; }
  const size_t D() const { return this->seas_diff_D; }
  const size_t Q() const { return this->sarma_Q; }
  const size_t period() const { return this->s_period; }

private:
  size_t arma_p, diff_d, arma_q, sarma_P, seas_diff_D, sarma_Q, s_period;
};

#endif
