#ifndef TIMER_STOPWATCH
#define TIMER_STOPWATCH

#include <chrono>
#include <type_traits>
#include <iostream>

// graciously taken from: https://stackoverflow.com/a/61881422
template <typename Resolution = std::chrono::duration<double,std::micro>>
class Stopwatch {
  typedef std::chrono::steady_clock Clock;
private:
  std::chrono::time_point<Clock> last;
public:
  void reset() noexcept {
    last = Clock::now();
  }   
  Stopwatch() noexcept {
    reset();
  }   
  auto operator()() const noexcept {// returns time in Resolution
    return Resolution(Clock::now() - last).count();
  }
  ~Stopwatch() {
    std::cout << "This code took: " << (*this)() * 1e-6 << " seconds.\n";
  }
};



#endif
