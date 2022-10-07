#ifndef SOLVER_HEADER
#define SOLVER_HEADER

#include <cmath>
#include "rng.h"
#include "utils.h"

// get secant value (aka linearly interpolate)
template <class ObjectiveFun, typename T> T secant( const ObjectiveFun & f,
                                                    T& new_val,
                                                    T& old_val,
                                                    T& alternative,
                                                    T tol = 0.0001) {
  // secant
  T new_obj = f(new_val);
  T old_obj = f(old_val);
  if( abs(new_obj - old_obj) < tol ) {
    return alternative;
  }
  return new_val - ((new_val - old_val)/( new_obj - old_obj )*new_obj);
}

template <typename T> bool between( T &x, T &a, T &b ) {
  if( a < b ) {
     return x < b && x > a;
  }
  return x > b && x < a;
}

template <typename T> void nudge( std::vector<T> &particles,
                                  const std::vector<T>& scores ) {

  T best = particles[0];
  T best_score = scores[0];
  // find best particle
  for( int i=1; i < particles.size(); i++ ) {
    // check if this particle has a better score than the best score
    // if yes, update best particle
    if( scores[i] < best_score ) {
      best = particles[i];
      best_score = scores[i];
    }
  }
  // nudge all particles towards best particle - get them "halfway there"
  for( auto &particle:particles ) {
    particle = (particle + best)/2;
  }
}

template <class ObjectiveFun,
          typename NumType=double> NumType particle_drainer( ObjectiveFun &f,
                                                             std::vector<NumType> &x,
                                                             NumType lower,
                                                             NumType upper,
                                                             int n_particles = 20,
                                                             NumType tol = 0.00001,
                                                             int max_iter = 200) {
    xoshiro gen;
    std::vector<NumType> particles(n_particles);
    for( int i=0; i < n_particles; i++ ) {
      particles[i] = runif<NumType, xoshiro>(lower, upper, gen);
    }
    int iter = 0;
    std::vector<NumType> scores(n_particles);
    while(iter <= max_iter) {
      // get particle scores
      scores = map(particles, f);
      // for( auto &particle:scores ) {
      //   std::cout << particle << " ";
      // }
      // std::cout << " " << std::endl;

      for( auto &particle:particles ) {
        std::cout << particle << " ";
      }
      std::cout << " - pre nudge. " << std::endl;

      // update particles by "nudging" them towards the best particle
      // this works by reference - hence particles get updated
      nudge( particles, scores);
      for( auto &particle:particles ) {
        std::cout << particle << " ";
      }
      std::cout << " - post nudge." << std::endl;
      if( all_const( particles, tol ) ) {
        std::cout << "Are all same." << std::endl;
        break;
      }
      iter++;
    }
  return mean( particles );
}



// template <class ObjectiveFun,
//           typename NumType> NumType dekker( const ObjectiveFun & f,
//                                             NumType lower_bound,
//                                             NumType upper_bound) {
//     NumType contrapoint, guess, bisection, old_guess, new_guess;
//     contrapoint = lower_bound, guess = upper_bound, old_guess = lower_bound;
//     while( true ) {
//       // bisection result:
//       bisection = (guess + contrapoint)/2;
//       // secant with bisection as alternative result
//       new_guess = secant( f, guess, old_guess, bisection );
//       if( between(new_guess, guess, bisection) ) {
//         guess = new_guess;
//       }
//       else {
//         guess = bisection;
//       }
//
//
//     }
//
//
//
//
//
//
//
// }














#endif
