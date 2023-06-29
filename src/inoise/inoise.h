#ifndef INTEGRATED_NOISE
#define INTEGRATED_NOISE

#include <random>

#include "utils/boxcox.h"
#include "utils/utils.h"

#include "common/forecast_result.h"

template <typename scalar_t> class sqrtStabilizer{
public:
  sqrtStabilizer<scalar_t>(){};
  sqrtStabilizer<scalar_t>(const std::vector<scalar_t> &y){};
  std::vector<scalar_t> transform(const std::vector<scalar_t> &x) {
    std::vector<scalar_t> result(x.size());
    for(size_t i = 0; i < x.size(); i++) {
      result[i] = std::sqrt(x[i]);
    }
    return result;
  }
  std::vector<scalar_t> inverse_transform(const std::vector<scalar_t> &x) {
    std::vector<scalar_t> result(x.size());
    for(size_t i = 0; i < x.size(); i++) {
      result[i] = std::pow(x[i],2);
    }
    return result;
  }
  const scalar_t inverse_transform(const scalar_t y) {
    return std::pow(y,2);
  }
  template <typename ForwardIt> void inverse_transform(
      ForwardIt first, ForwardIt last) {
    for(; first!= last; first++) {
      *first = std::move(std::pow(*first, 2));
    }
  }
};

template <typename scalar_t> class logStabilizer{
public:
  logStabilizer<scalar_t>(){};
  logStabilizer<scalar_t>(const std::vector<scalar_t> &y){};
  std::vector<scalar_t> transform(const std::vector<scalar_t> &x) {
    std::vector<scalar_t> result(x.size());
    for(size_t i = 0; i < x.size(); i++) {
      result[i] = std::log(x[i]);
    }
  }
  std::vector<scalar_t> inverse_transform(const std::vector<scalar_t> &x) {
    std::vector<scalar_t> result(x.size());
    for(size_t i = 0; i < x.size(); i++) {
      result[i] = std::exp(x[i]);
    }
    return result;
  }
  const scalar_t inverse_transform(const scalar_t y) {
    return std::exp(y);
  }
  template <typename ForwardIt> void inverse_transform(
      ForwardIt first, ForwardIt last) {
    for(; first!= last; first++) {
      *first = std::move(std::exp(*first));
    }
  }
};

template <typename scalar_t> class boxCoxStabilizer{
private:
  BoxCoxTransformer<scalar_t> transformer;
public:
  boxCoxStabilizer<scalar_t>(const std::vector<scalar_t> &y){
    this->transformer = BoxCoxTransformer<scalar_t>(y);
  };
  std::vector<scalar_t> transform(const std::vector<scalar_t> &x) {
    return this->transformer.transform(x);
  }
  std::vector<scalar_t> inverse_transform(const std::vector<scalar_t> &x) {
    return this->transformer.inverse_transform(x);
  }
  template <typename ForwardIt> void inverse_transform(
      ForwardIt first, ForwardIt last) {
    for(; first!= last; first++) {
      *first = this->transformer.inverse_transform(*first);
    }
  }
};

template <typename scalar_t> class identityStabilizer{
public:
  identityStabilizer<scalar_t>(){};
  identityStabilizer<scalar_t>(const std::vector<scalar_t> &y){};
  std::vector<scalar_t> transform(
      const std::vector<scalar_t> &x) {return x;}
  std::vector<scalar_t> inverse_transform(
      const std::vector<scalar_t> &x) {return x;}
  template <typename ForwardIt> void inverse_transform(
      ForwardIt first, ForwardIt last) {}
};

using mt19937_size_t = std::mersenne_twister_engine<size_t, 64, 312, 156, 31,
                                                    0xb5026f5aa96619e9, 29,
                                                    0x5555555555555555, 17,
                                                    0x71d67fffeda60000, 37,
                                                    0xfff7eee000000000, 43, 6364136223846793005>;


template <typename scalar_t,
          typename Stabilizer=sqrtStabilizer<scalar_t>> class INoise {
private:
  std::vector<scalar_t> y, diff_y, diff_y_s;
  size_t d, s_d, season;
  Stabilizer stabilizer;
  mt19937_size_t twister;
  scalar_t last_y;
  bool fitted;
public:
  INoise<scalar_t, Stabilizer>() {};
  INoise<scalar_t, Stabilizer>(const std::vector<scalar_t> &y, 
                               const size_t differences = 1,
                               const size_t seasonal_difference = 10,
                               const bool stabilize = true) :
  y(y), d(differences) {
    this->stabilizer = Stabilizer(y);
    this->s_d = seasonal_difference;
    this->twister = mt19937_size_t{std::random_device{}()};
    // if(seasonal_difference < 0) {
      // find seasonality instead
    // }
    this->last_y = 0;
    this->fitted = false;
  }
  void fit() {
    // apply variance stabilizing transformation
    this->diff_y = this->stabilizer.transform(this->y);
    this->last_y = this->diff_y[this->diff_y.size()-1];
    // take difference 
    if(this->d > 0) {
      this->diff_y = diff(this->diff_y, 1, this->d);
    }
    this->diff_y_s = std::vector<scalar_t>(this->diff_y.size() - this->s_d, 0);
    // take seasonal difference (if any)
    if(this->s_d > 1) {
      this->diff_y_s = diff(this->diff_y, this->s_d, 1);
    }
    this->fitted = true;
  }
  forecast_result<scalar_t> forecast(
      const size_t h = 10,
      const bool produce_se = true,
      const size_t samples = 10000) {
    // keep track of average and standard errors 
    std::vector<scalar_t> avg(h), std_err(h);
    if(produce_se) {
      // intermediate result temporary
      std::vector<scalar_t> temp(h);
      // h * 2 as this will hold intermediaries for both the averages and 
      // standard errors when running welfords algorithm
      std::vector<scalar_t> welfords_temp(h*2);//, 0.0);
      // initialize welfords algorithm from first sample
      for(size_t j = 0; j < h; j++) {
        // draw and create temp
        temp[j] = this->diff_y[this->diff_y.size() - this->s_d + (j % this->s_d)] + 
          draw_from_vec(this->diff_y_s, this->twister);
      }
      // run cummulative sum over current range 
      *(temp.begin()) += this->last_y;
      cummulative_sum(temp.begin(), temp.end());
      // apply inverse transform
      this->stabilizer.inverse_transform(temp.begin(), temp.end());
      // initialize welfords algorithm
      for(size_t j = 0; j < h; j++) {
        welfords_temp[2*j] = temp[j];
        welfords_temp[(2*j)+1] = 0.0;
      }
      // finish initialization, run actual sampling
      for(size_t i = 1; i < samples; i++) {
        for(size_t j = 0; j < h; j++) {
          temp[j] = this->diff_y[this->diff_y.size() - this->s_d + (j % this->s_d)] + 
            draw_from_vec(this->diff_y_s, this->twister);
        }
        *(temp.begin()) += this->last_y;
        cummulative_sum(temp.begin(), temp.end());
        this->stabilizer.inverse_transform(temp.begin(), temp.end());
        // run welfords algorithm
        for(size_t j = 0; j < h; j++) {
          const scalar_t current_val = temp[j], prev_mean = welfords_temp[2*j];
          welfords_temp[2*j] += (current_val - prev_mean)/static_cast<scalar_t>(i);
          welfords_temp[(2*j)+1] += (current_val - welfords_temp[2*j]) * (current_val - prev_mean);
        }
      }
      // move results out of welfords temp 
      for(size_t j = 0; j < h; j++) {
        avg[j] = welfords_temp[2*j];
        std_err[j] = std::sqrt(welfords_temp[(2*j)+1]/(samples-1));
      }
    }
    else{
      // the seasonal mean is correct
      scalar_t mean_seas_val = mean(this->diff_y_s);
      // repeat pattern from last s_d observations until we have h/s_d repeats
      for(size_t i = 0; i < h; i++) {
        /* i % s_d allows you to correctly fold over seasonal periods 
         * because we might need e.g. 20 repeats of the last 12 values, so repeat 
         * all 12 values once, and then get the first 8 values in that series again
         * after looping around */
        avg[i] = this->diff_y[this->diff_y.size() - this->s_d + (i % this->s_d)] + mean_seas_val;
      }
      *(avg.begin()) += this->last_y;
      cummulative_sum(avg.begin(), avg.end());
      this->stabilizer.inverse_transform(avg.begin(), avg.end());
      // for std_err just get observed values and compute standard error based on them
      auto sigma2 = welfords_algorithm(this->y)[1];
      std::fill(std_err.begin(), std_err.end(), sigma2);
    }
    return forecast_result(avg, std_err);      
  }
  const bool is_fitted() const { return this->fitted; }
};


#endif
