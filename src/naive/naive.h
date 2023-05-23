#ifndef BLAZE_NAIVE
#define BLAZE_NAIVE

#include "utils/utils.h"
#include "utils/xreg.h"

#include "common/forecast_result.h"

template <typename scalar_t> class RW {
private:
  // data
  std::vector<scalar_t> y;
  size_t lags;
  bool drift;
  std::vector<scalar_t> model, fitted_vals, residuals;
  scalar_t drift_val;
  bool fitted;
public:
  RW<scalar_t>(){};
  RW<scalar_t>(const std::vector<scalar_t> &y,
               const size_t lags = 1,
               const bool drift = true) : 
    y(y), lags(lags), drift(drift){
  };
  void fit() {
    this->model = std::vector<scalar_t>(this->lags);
    this->fitted_vals = std::vector<scalar_t>(this->y.size() - this->lags);
    this->drift_val = 0;
    // take lagged values  
    for(size_t j = 0; j < this->lags; j++) {
      this->model[j] = this->y[this->y.size() - 1 - this->lags + j];
    }
    // compute fit and residuals 
    for(size_t j = this->lags; j < this->y.size(); j++) {
      auto offset = j - this->lags;
      auto curr_fit = y[j] - y[offset];
      this->fitted_vals[offset] = curr_fit;
      auto resid = y[j] - curr_fit;
      this->residuals[offset] = resid;
      // also compute optional drift 
      this->drift_val += this->drift * resid;
    }
    // drift is average residual
    this->drift_val /= this->y.size() - this->lags();
    this->sigma2 = sum_of_squares(this->residuals)/this->y.size();
    this->fitted = true;
  };
  forecast_result<scalar_t> forecast(
      const size_t h = 10) {
    if (!this->fitted)
      return forecast_result<scalar_t>(0);
    
  };
  const std::vector<scalar_t> get_residuals() const { return this->residuals; }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_vals; }
  const scalar_t get_sigma2() const { return this->sigma2; }
  const bool is_fitted() const { return this->fitted; }
};



// average_tile <- function(x, tiler = mean) {
//   seasonalities <- sort(find_seasonalities(x), TRUE)
//   
//   res <- purrr::map(seasonalities, function(seas) seq_len(seas) )
//   df <- purrr::map( res, ~ rep(.x, ceiling(length(x)/length(.x)))[seq_len(length(x))] )
//   df <- do.call(cbind, df)
//     colnames(df) <- paste0( "seas_", seasonalities )
//     df <- as.data.frame(cbind( y = x, df ))
//     
//     df <- dplyr::group_by(df, !!!rlang::syms(paste0( "seas_", seasonalities )))
//     res <- dplyr::summarise( df, dplyr::across( .fns = tiler), .groups = "keep" )
//     res <- dplyr::ungroup(res)
//     res <- dplyr::arrange( res,
//                            !!!rlang::syms(paste0( "seas_", seasonalities ))
//     )
//     df <- dplyr::ungroup(df)
//     
//     last_row <- dplyr::select( df[ length(x),],
//                                tidyselect::all_of(
//                                  paste0( "seas_", seasonalities )
//                                )
//     )
//     
//     return( structure(
//         list( seasonalities = seasonalities,
//               last_seas = last_row,
//               tiles = res ),
//               class = "average_tile")
//     )
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
