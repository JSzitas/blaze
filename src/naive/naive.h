#ifndef BLAZE_NAIVE
#define BLAZE_NAIVE

#include "utils/utils.h"
#include "utils/xreg.h"
#include "utils/seasonality_utils.h"

#include "common/forecast_result.h"
#include "structures/circulant.h"

template <typename scalar_t> scalar_t sum_positive_squares(const size_t n) {
  return n*(n+1)*((2*n) + 1)/6;
}

template <typename scalar_t> scalar_t drift_coef(
    const std::vector<scalar_t> & x) {
  const scalar_t scale = inv_sum_positive_squares(x.size());
  scalar_t res = 0.0;
  for(size_t i = 0; i < x.size(); i++) res += (i + 1) * x[i];
  return res * scale;
}

template <typename scalar_t> class RW {
private:
  // data
  std::vector<scalar_t> y;
  size_t lags;
  bool drift;
  scalar_t sigma2;
  std::vector<scalar_t> fitted_vals, residuals;
  scalar_t drift_val, drift_se;
  bool fitted;
public:
  RW<scalar_t>(){};
  RW<scalar_t>(const std::vector<scalar_t> &y,
               const size_t lags = 1,
               const bool drift = true) : 
    y(y), lags(lags), drift(drift), sigma2(0.0),
    fitted_vals(std::vector<scalar_t>(y.size(), 0.0)), 
    residuals(std::vector<scalar_t>(y.size(), 0.0)),
    drift_val(0.0), drift_se(0.0), fitted(false){};
  void fit() {
    this->fitted_vals = std::vector<scalar_t>(this->y.size() - this->lags);
    this->drift_val = 0.0;
    this->drift_se = 0.0;
    scalar_t drift_temp = 0.0;
    // compute fit and residuals 
    for(size_t j = this->lags; j < this->y.size(); j++) {
      auto offset = j - this->lags;
      auto curr_fit = y[j] - y[offset];
      this->fitted_vals[offset] = curr_fit;
      auto resid = y[j] - curr_fit;
      this->residuals[offset] = resid;
      // also compute optional drift 
      if(this->drift) this->drift_val += (j + 1) * resid;
    }
    // potentially adjust residuals for drift? 
    // compute sigma2
    this->sigma2 = sum_of_squares(this->residuals)/
      (this->y.size() - this->lags() - this->drift);
    if(this->drift) {
      // drift computation 
      // this is a fast way to compute what is normally (X'X)^-1 - since this 
      // reduces to an inverse of a single number (as X is a vector), and the 
      // product X'X is just sum( x_i^2 for i in 1:n ) - which has a fast formula
      // finally, since we are dropping the first lags observations, the appropriate
      // way to account for this is to actually subtract their sum to have a drift 
      // coefficient which is index consistent
      const scalar_t scale = 1/(sum_positive_squares(this->y.size()) - 
                                sum_positive_squares(this->lags));
      this->drift_val *= scale;
      // also estimate drift standard error
      this->drift_se = std::sqrt(this->sigma2 * scale); 
    }
    this->fitted = true;
  };
  forecast_result<scalar_t> forecast(
      const size_t h = 10) {
    if (!this->fitted)
      return forecast_result<scalar_t>(0);
    std::vector<scalar_t> fcst(h), std_err(h);
    for(size_t i = 0; i < h; i++) {
      fcst[i] = this->y[this->y.size() - h + i] + (this->drift_val * i);
    }
    for(size_t i = 0; i < h; i++) {
      std_err[i] = std::sqrt(this->sigma2 * i + std::pow(i * this->drift_se));
    }
    return forecast_result<scalar_t>(fcst, std_err);
  };
  const std::vector<scalar_t> get_residuals() const { return this->residuals; }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_vals; }
  const scalar_t get_sigma2() const { return this->sigma2; }
  const bool is_fitted() const { return this->fitted; }
};

template <typename scalar_t> class AverageTile {
private:
  // data
  std::vector<scalar_t> y;
  std::vector<size_t> seasons;
  size_t last_season;
  bool fitted;
public:
  AverageTile<scalar_t>(){};
  AverageTile<scalar_t>(const std::vector<scalar_t> &y,
                        const std::vector<size_t> seasonalities = {},
                        const size_t max_seas_iter = 5,
                        const size_t max_seas_size = 1500) : 
    y(y){
    if(seasonalities.size() == 0) {
      this->seasons = find_seasonalities(y, max_seas_iter, max_seas_size);
    }
    else {
      this->seasons = seasonalities;
    }
    this->last_season = 0;
    this->fitted = false;
  };
  void fit(){
    // prepare tiles - we keep a unique tile for each averaged seasonal mean
    // across all seasonalities
    std::vector<scalar_t> tiles(this->seasons[this->seasons.size()-1]);
    for(auto &season:seasons) {
      // for the current season, keep a number of running means, and 
      // update them using 
    }
    
    
    
    
    
  }
  forecast_result<scalar_t> forecast(
      const size_t h = 10) {
    if (!this->fitted)
      return forecast_result<scalar_t>(0);
    // carry out respective tiling using last known periods 
    
    
  };
  const std::vector<scalar_t> get_residuals() const { return this->residuals; }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_vals; }
  const scalar_t get_sigma2() const { return this->sigma2; }
  const bool is_fitted() const { return this->fitted; }
};


// average_tile <- function(x, max_seasonalities = 5, tiler = mean, trend_prefilter = FALSE) {
//   seasonalities <- sort(find_seasonalities(x), TRUE)#[seq_len(max_seasonalities)]
// 
//   if(trend_prefilter) {
//     total_avg = mean(x)
//     
//     
//   }
//   else {
//     res <- purrr::map(seasonalities, function(seas) seq_len(seas) )
//     df <- purrr::map( res, ~ rep(.x, ceiling(length(x)/length(.x)))[seq_len(length(x))] )
//     df <- do.call(cbind, df)
//     colnames(df) <- paste0( "seas_", seasonalities )
//     df <- as.data.frame(cbind( y = x, df ))
//     df <- dplyr::group_by(df, !!!rlang::syms(paste0( "seas_", seasonalities)))
//     res <- dplyr::summarise( df, dplyr::across( .fns = tiler), .groups = "keep")
//     res <- dplyr::ungroup(res)
//     res <- dplyr::arrange( res, !!!rlang::syms(paste0( "seas_", seasonalities )))
//     df <- dplyr::ungroup(df)
//     last_row <- dplyr::select( df[ length(x),], tidyselect::all_of(
//         paste0( "seas_", seasonalities )))
//   }
//   return(structure(
//       list( seasonalities = seasonalities, last_seas = last_row, tiles = res ),
//       class = "average_tile")
//            )
// }
// 
// rep_df <- function(x, k = 2) {
//   result <- data.frame()
//   for(i in seq_len(k)) {
//     result <- rbind(result, x)
//   }
//   return(result)
// }
// 
// 
// forecast.average_tile <- function(object, h = 10,... ) {
// 
//   last_obs <- object$last_seas
//   tiles <- object$tiles
//   seasonalities <- object$seasonalities
//   first_pred <- c( unlist(last_obs + 1)) %% seasonalities
// 
// 
//   new_df <- rep_df( tiles, ceiling( h/nrow(tiles) )+1 )
//   for( row in seq_len(nrow(new_df)) ) {
// 
//     res <- new_df[row, intersect( colnames( last_obs), colnames(new_df) ) ] - first_pred
//     if( sum(unlist(res)) == 0 ) {
//       res <- new_df[ row:(row + h-1),]
//       break;
//     }
//   }
//   return(res)
// }

#endif
