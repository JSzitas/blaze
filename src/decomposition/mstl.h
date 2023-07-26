#ifndef BLAZE_MSTL
#define BLAZE_MSTL

#include "decomposition/stl.h"
#include "utils/seasonality_utils.h"

#include "decomposition/decompose_result.h"

template <typename scalar_t> class MSTL {
  const int trend_size, low_pass_size, seasonal_degree, trend_degree, 
  low_pass_degree, seasonal_jump, trend_jump, low_pass_jump;
  const bool robust;
  std::vector<scalar_t> deseasonalized, trend, remainder;
  std::vector<int> s_periods, s_windows;
  std::vector<std::vector<scalar_t>> seasonalities;
public:
  MSTL<scalar_t>(){}
  MSTL<scalar_t>(const std::vector<scalar_t> &y,
                 const std::vector<int> &periods = {}, 
                 const std::vector<int> &seasonal = {},
                 const int trend_size = 0, 
                 const int low_pass_size = 0,
                 const int seasonal_degree = 1,
                 const int trend_degree = 1,
                 const int low_pass_degree = 1,
                 const int seasonal_jump = 1,
                 const int trend_jump = 1,
                 const int low_pass_jump = 1,
                 const bool robust = false) : trend_size(trend_size),
                 low_pass_size(low_pass_size), seasonal_degree(seasonal_degree),
                 trend_degree(trend_degree), low_pass_degree(low_pass_degree),
                 seasonal_jump(seasonal_jump), trend_jump(trend_jump),
                 low_pass_jump(low_pass_jump), robust(robust),
                 deseasonalized(y), trend(std::vector<scalar_t>(y.size(), 0.0)),
                 remainder(std::vector<scalar_t>(y.size(), 0.0)) {
    if(periods.size() == 0) {
      this->s_periods = find_seasonalities<scalar_t, int>(y);
    } else{
      this->s_periods = periods;
    }
    if(seasonal.size() == 0) {
      this->s_windows = std::vector<int>(this->s_periods.size(), 7);
      int increment = 4;
      for(auto &window:this->s_windows) {
        window += increment;
        increment += 4;
      }
    } else {
      this->s_windows = seasonal;
    }
    this->seasonalities = std::vector<std::vector<scalar_t>>(
      this->s_periods.size(), std::vector<scalar_t>(y.size(), 0.0));
  }
  void fit() {
    std::cout << "Fitting MSTL model" << std::endl;
    size_t i = 0; 
    for(; i < this->s_periods.size()-1; i++) fit_seasonality(i);
    // this is only the last seasonality - we save the trend as well here
    if(i < this->s_periods.size()) fit_seasonality<true>(i);
    for(size_t i = 0; i < deseasonalized.size(); i++) {
      this->remainder[i] = this->deseasonalized[i] - this->trend[i];
    }
    return;
  }
  const std::vector<scalar_t> get_trend() const { return this->trend; }
  const std::vector<scalar_t> get_season(const size_t which_season) const {
    size_t intended_season = which_season;
    if(which_season < this->seasonalities.size()) {
      std::cout << "You asked for more seasons than were available in get_season " <<
        " you specified you wanted season (which_season): " << which_season << 
          " but only " << (this->seasonalities.size()-1) << " seasons are available."
          << " Returning last available season instead" << std::endl;
      intended_season = this->seasonalities.size()-1;    
    }
    return this->seasonalities[intended_season]; 
  }
  const std::vector<scalar_t> get_remainder() const { return this->residual; }
  const MultipleDecompositionResult<scalar_t> get_decomposition_result() const {
    return MultipleDecompositionResult<scalar_t>(
      this->trend, this->seasonalities, this->remainder);
  }
  void print_summary() {
    std::cout << "Seasonal periods: ";
    print_vector(this->s_periods);
    std::cout << "Seasonal windows: ";
    print_vector(this->s_windows);
  }
private:
  template <const bool save_trend = false>
  void fit_seasonality(const size_t seasonal_period) {
    std::cout << "Fitting to seasonality: " << seasonal_period << std::endl;
    auto model = STL(this->deseasonalized, this->s_periods[seasonal_period],
                     this->s_windows[seasonal_period], this->trend_size,
                     this->low_pass_size, this->seasonal_degree, 
                     this->trend_degree, this->low_pass_degree, 
                     this->seasonal_jump, this->trend_jump, this->low_pass_jump,
                     this->robust);
    model.fit();
    // recover seasonal component
    const std::vector<scalar_t> seasonal_component = model.get_season();
    // save seasonality 
    this->seasonalities[seasonal_period] = seasonal_component;
    // deseasonalize 
    for(size_t i = 0; i < this->deseasonalized.size(); i++)
      deseasonalized[i] -= seasonal_component[i];
    // save trend - technically only necessary for the final period, i.e the last trend
    if constexpr(save_trend) {
      std::cout << "Saving trend" << std::endl;
      this->trend = model.get_trend(); 
    }
    return;
  }
};

#endif
