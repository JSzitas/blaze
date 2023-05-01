#ifndef FORECAST_RESULT
#define FORECAST_RESULT

template <typename scalar_t = float> struct forecast_result {
  forecast_result<scalar_t>( std::vector<scalar_t> &forecasts,
                             std::vector<scalar_t> &std_errs ) :
  forecast(std::move(forecasts)), std_err(std::move(std_errs)){}
  forecast_result<scalar_t>( std::vector<scalar_t> &&forecasts,
                             std::vector<scalar_t> &&std_errs ) :
  forecast(std::move(forecasts)), std_err(std::move(std_errs)){}
  forecast_result<scalar_t>(size_t h) {
    this->forecast = std::vector<scalar_t>(h);
    this->std_err = std::vector<scalar_t>(h);
  }
  void add(size_t i, scalar_t fcst, scalar_t se) {
    this->forecast[i] = fcst;
    this->std_err[i] = se;
  }
  std::vector<scalar_t> forecast;
  std::vector<scalar_t> std_err;
};


#endif
