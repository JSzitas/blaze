#ifndef ARIMA_CSS_LIKELIHOOD
#define ARIMA_CSS_LIKELIHOOD

#include "arima/structures/structural_model.h"
#include "arima/structures/arima_kind.h"

#include "third_party/eigen.h"

#include "utils/utils.h"

/* arma is p, q, sp, sq, ns, d, sd
 * Note that this function is very similar to the one that follows it
 * the main difference is in what they return -
 */
double arima_css_ssq(const Eigen::VectorXd &y, const Eigen::VectorXd &pars,
                     const arima_kind &kind, const int n_cond,
                     std::vector<double> &resid) {

  const size_t n = y.size(), p = kind.p() + kind.period() * kind.P(),
            q = kind.q() + kind.period() * kind.Q();
  // prepare the residuals - possibly move this out and never allocate here?
  int ma_offset, nu = n-n_cond;
  double ssq = 0.0, tmp = 0.0;
  for (size_t l = n_cond; l < n; l++) {
    ma_offset = min(l - n_cond, q);
    tmp = y[l];
    for (size_t j = 0; j < p; j++) {
      tmp -= pars[j] * y[l - j - 1];
    }
    // to offset that this is all in one vector, we need to
    // start at p and go to p + q
    for (size_t j = 0; j < ma_offset; j++) {
      tmp -= pars[p + j] * resid[l - j - 1];
    }
    resid[l] = tmp;
    if (!isnan(tmp)) {
      ssq += tmp * tmp;
    }
    else {
      nu--;
    }
  }
  return ssq / nu;
}

template <const size_t update_point = 0,
          class T> void arima_steady_state(const T &y,
                                           structural_model<double> &model) {
  const size_t n = y.size(), rd = model.a.size(), p = model.phi.size(),
    q = model.theta.size(), d = model.delta.size(), r = rd - d;

    std::vector<double> anew(rd);
    std::vector<double> M(rd);
    std::vector<double> mm( (d > 0) * rd * rd);

    for (size_t l = 0; l < n; l++) {
      for (size_t i = 0; i < r; i++) {
        double tmp = (i < r - 1) ? model.a[i + 1] : 0.0;
        if (i < p) tmp += model.phi[i] * model.a[0];
        anew[i] = tmp;
      }
      if (d > 0) {
        for (size_t i = r + 1; i < rd; i++) anew[i] = model.a[i - 1];
        double tmp = model.a[0];
        for (size_t i = 0; i < d; i++) tmp += model.delta[i] * model.a[r + i];
        anew[r] = tmp;
      }
      if(l > update_point) {
        if (d == 0) {
          for (size_t i = 0; i < r; i++) {
            double vi = 0.0;
            if (i == 0) vi = 1.0; else if (i - 1 < q) vi = model.theta[i - 1];
            for (size_t j = 0; j < r; j++) {
              double tmp = 0.0;
              if (j == 0) tmp = vi; else if (j - 1 < q) tmp = vi * model.theta[j - 1];
              if (i < p && j < p) tmp += model.phi[i] * model.phi[j] * model.P[0];
              if (i < r - 1 && j < r - 1) tmp += model.P[i + 1 + r * (j + 1)];
              if (i < p && j < r - 1) tmp += model.phi[i] * model.P[j + 1];
              if (j < p && i < r - 1) tmp += model.phi[j] * model.P[i + 1];
              model.Pn[i + r * j] = tmp;
            }
          }
        } else {
          /* mm = TP */
          for (size_t i = 0; i < r; i++)
            for (size_t j = 0; j < rd; j++) {
              double tmp = 0.0;
              if (i < p) tmp += model.phi[i] * model.P[rd * j];
              if (i < r - 1) tmp += model.P[i + 1 + rd * j];
              mm[i + rd * j] = tmp;
            }
            for (size_t j = 0; j < rd; j++) {
              double tmp = model.P[rd * j];
              for (size_t k = 0; k < d; k++)
                tmp += model.delta[k] * model.P[r + k + rd * j];
              mm[r + rd * j] = tmp;
            }
            for (size_t i = 1; i < d; i++)
              for (size_t j = 0; j < rd; j++)
                mm[r + i + rd * j] = model.P[r + i - 1 + rd * j];

          /* Pnew = mmT' */
          for (size_t i = 0; i < r; i++)
            for (size_t j = 0; j < rd; j++) {
              double tmp = 0.0;
              if (i < p) tmp += model.phi[i] * mm[j];
              if (i < r - 1) tmp += mm[rd * (i + 1) + j];
              model.Pn[j + rd * i] = tmp;
            }
            for (size_t j = 0; j < rd; j++) {
              double tmp = mm[j];
              for (size_t k = 0; k < d; k++)
                tmp += model.delta[k] * mm[rd * (r + k) + j];
              model.Pn[rd * r + j] = tmp;
            }
            for (size_t i = 1; i < d; i++)
              for (size_t j = 0; j < rd; j++)
                model.Pn[rd * (r + i) + j] = mm[rd * (r + i - 1) + j];
          /* Pnew <- Pnew + (1 theta) %o% (1 theta) */
          for (size_t i = 0; i <= q; i++) {
            double vi = (i == 0) ? 1. : model.theta[i - 1];
            for (size_t j = 0; j <= q; j++)
              model.Pn[i + rd * j] += vi * ((j == 0) ? 1. : model.theta[j - 1]);
          }
        }
      }
      if (!isnan(y[l])) {
        double resid = y[l] - anew[0];
        for (size_t i = 0; i < d; i++)
          resid -= model.delta[i] * anew[r + i];

        for (size_t i = 0; i < rd; i++) {
          double tmp = model.Pn[i];
          for (size_t j = 0; j < d; j++)
            tmp += model.Pn[i + (r + j) * rd] * model.delta[j];
          M[i] = tmp;
        }

        double gain = M[0];
        for (size_t j = 0; j < d; j++) gain += model.delta[j] * M[r + j];
        for (size_t i = 0; i < rd; i++)
          model.a[i] = anew[i] + M[i] * resid / gain;
        for (size_t i = 0; i < rd; i++)
          for (size_t j = 0; j < rd; j++)
            model.P[i + j * rd] = model.Pn[i + j * rd] - M[i] * M[j] / gain;
      } else {
        for (size_t i = 0; i < rd; i++) model.a[i] = anew[i];
        for (size_t i = 0; i < rd * rd; i++) model.P[i] = model.Pn[i];
      }
    }
}

template <const size_t update_point = 0,
          typename T> std::array<double,2> arima_likelihood(
              const T &y,
              structural_model<double> &model) {
  const size_t rd = model.a.size(), d = model.delta.size();

  std::vector<double> anew(rd);
  std::vector<double> M(rd);
  std::vector<double> mm( (d > 0) * rd * rd);
  return arima_likelihood_impl(y, model, anew, M, mm);
}

template <const size_t update_point = 0,
          typename T> std::array<double,2> arima_likelihood(
              const T &y, structural_model<double> &model,
              std::vector<double> &anew,
              std::vector<double> &M,
              std::vector<double> &mm) {
  return arima_likelihood_impl(y, model, anew, M, mm);
}

template <const size_t update_point = 0,
          typename T> std::array<double,2> arima_likelihood_impl(
              const T &y,
              structural_model<double> &model,
              std::vector<double> &anew,
              std::vector<double> &M,
              std::vector<double> &mm) {
  const size_t n = y.size(), rd = model.a.size(), p = model.phi.size(),
    q = model.theta.size(), d = model.delta.size(), r = rd - d;

  double ssq = 0, sumlog = 0;
  size_t nu = n;
  for (size_t l = 0; l < n; l++) {
    for (size_t i = 0; i < r; i++) {
      double tmp = (i < r - 1) ? model.a[i + 1] : 0.0;
      if (i < p) tmp += model.phi[i] * model.a[0];
      anew[i] = tmp;
    }
    if (d > 0) {
      for (size_t i = r + 1; i < rd; i++) anew[i] = model.a[i - 1];
      double tmp = model.a[0];
      for (size_t i = 0; i < d; i++) tmp += model.delta[i] * model.a[r + i];
      anew[r] = tmp;
    }
    if(l > update_point) {
      if (d == 0) {
        for (size_t i = 0; i < r; i++) {
          double vi = 0.0;
          if (i == 0) vi = 1.0; else if (i - 1 < q) vi = model.theta[i - 1];
          for (size_t j = 0; j < r; j++) {
            double tmp = 0.0;
            if (j == 0) tmp = vi; else if (j - 1 < q) tmp = vi * model.theta[j - 1];
            if (i < p && j < p) tmp += model.phi[i] * model.phi[j] * model.P[0];
            if (i < r - 1 && j < r - 1) tmp += model.P[i + 1 + r * (j + 1)];
            if (i < p && j < r - 1) tmp += model.phi[i] * model.P[j + 1];
            if (j < p && i < r - 1) tmp += model.phi[j] * model.P[i + 1];
            model.Pn[i + r * j] = tmp;
          }
        }
      } else {
        /* mm = TP */
        for (size_t i = 0; i < r; i++)
          for (size_t j = 0; j < rd; j++) {
            double tmp = 0.0;
            if (i < p) tmp += model.phi[i] * model.P[rd * j];
            if (i < r - 1) tmp += model.P[i + 1 + rd * j];
            mm[i + rd * j] = tmp;
          }
          for (size_t j = 0; j < rd; j++) {
            double tmp = model.P[rd * j];
            for (size_t k = 0; k < d; k++)
              tmp += model.delta[k] * model.P[r + k + rd * j];
            mm[r + rd * j] = tmp;
          }
          for (size_t i = 1; i < d; i++)
            for (size_t j = 0; j < rd; j++)
              mm[r + i + rd * j] = model.P[r + i - 1 + rd * j];

        /* Pnew = mmT' */
        for (size_t i = 0; i < r; i++)
          for (size_t j = 0; j < rd; j++) {
            double tmp = 0.0;
            if (i < p) tmp += model.phi[i] * mm[j];
            if (i < r - 1) tmp += mm[rd * (i + 1) + j];
            model.Pn[j + rd * i] = tmp;
          }
          for (size_t j = 0; j < rd; j++) {
            double tmp = mm[j];
            for (size_t k = 0; k < d; k++)
              tmp += model.delta[k] * mm[rd * (r + k) + j];
            model.Pn[rd * r + j] = tmp;
          }
          for (size_t i = 1; i < d; i++)
            for (size_t j = 0; j < rd; j++)
              model.Pn[rd * (r + i) + j] = mm[rd * (r + i - 1) + j];
        /* Pnew <- Pnew + (1 theta) %o% (1 theta) */
        for (size_t i = 0; i <= q; i++) {
          double vi = (i == 0) ? 1. : model.theta[i - 1];
          for (size_t j = 0; j <= q; j++)
            model.Pn[i + rd * j] += vi * ((j == 0) ? 1. : model.theta[j - 1]);
        }
      }
    }
    if (!isnan(y[l])) {
      double resid = y[l] - anew[0];
      for (size_t i = 0; i < d; i++)
        resid -= model.delta[i] * anew[r + i];

      for (size_t i = 0; i < rd; i++) {
        double tmp = model.Pn[i];
        for (size_t j = 0; j < d; j++)
          tmp += model.Pn[i + (r + j) * rd] * model.delta[j];
        M[i] = tmp;
      }

      double gain = M[0];
      for (size_t j = 0; j < d; j++) gain += model.delta[j] * M[r + j];
      if (gain < 1e4) {
        ssq += resid * resid / gain;
        sumlog += log(gain);
      } else {
        nu--;
      }
      for (size_t i = 0; i < rd; i++)
        model.a[i] = anew[i] + M[i] * resid / gain;
      for (size_t i = 0; i < rd; i++)
        for (size_t j = 0; j < rd; j++)
          model.P[i + j * rd] = model.Pn[i + j * rd] - M[i] * M[j] / gain;
    } else {
      for (size_t i = 0; i < rd; i++) model.a[i] = anew[i];
      for (size_t i = 0; i < rd * rd; i++) model.P[i] = model.Pn[i];
    }
  }
  // finally, compute likelihood and return
  const double s2 = ssq/(double)nu;
  const double loglik = 0.5 * (log(s2) + (sumlog/(double)nu));
  return std::array<double,2>{s2, loglik};
}

#endif
