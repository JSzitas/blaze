// Copyright https://github.com/PatWie/CppNumericalSolvers, MIT license
#ifndef INCLUDE_CPPOPTLIB_FUNCTION_H_
#define INCLUDE_CPPOPTLIB_FUNCTION_H_

#include <optional>
#include <type_traits>
#include "utils/derivatives.h"

namespace cppoptlib::function {

// Specifies a current function state.
template <class scalar_t, class vector_t, class matrix_t>
struct State {
  int dim;
  int order;

  scalar_t value = 0;               // The objective value.
  vector_t x;                       // The current input value in x.
  vector_t gradient;                // The gradient in x.
  std::optional<matrix_t> hessian;  // The Hessian in x;

  // TODO(patwie): There is probably a better way.
  State() : dim(-1), order(-1) {}

  State(const int dim, const int order)
    : dim(dim),
      order(order),
      x(vector_t::Zero(dim)),
      gradient(vector_t::Zero(dim)) {
    if (order > 1) {
      hessian = std::optional<matrix_t>(matrix_t::Zero(dim, dim));
    }
  }

  State(const State<scalar_t, vector_t, matrix_t> &rhs) { CopyState(rhs); }

  State operator=(const State<scalar_t, vector_t, matrix_t> &rhs) {
    CopyState(rhs);
    return *this;
  }

  void CopyState(const State<scalar_t, vector_t, matrix_t> &rhs) {
    assert(rhs.order > -1);
    dim = rhs.dim;
    order = rhs.order;
    value = rhs.value;
    x = rhs.x.eval();
    if (order >= 1) {
      gradient = rhs.gradient.eval();
    }
    if ((order >= 2) && rhs.hessian) {
      hessian = std::optional<matrix_t>(rhs.hessian->eval());
    }
  }
};


// replace previous implementation of Function with CRTP
template <typename TScalar, typename DerivedFunction, int TDim = Eigen::Dynamic>
class Function {
public:
  static constexpr int Dim = TDim;
  using scalar_t = TScalar;
  using vector_t = Eigen::Matrix<TScalar, Dim, 1>;
  using hessian_t = Eigen::Matrix<TScalar, Dim, Dim>;
  using matrix_t = Eigen::Matrix<TScalar, Eigen::Dynamic, Eigen::Dynamic>;
  using index_t = typename vector_t::Index;

  using state_t = function::State<scalar_t, vector_t, hessian_t>;

  // template<class T>  using has_grad =
  // decltype(std::declval<T&>().Gradient(std::declval<const vector_t&, vector_t*>()));
  // template<class T>  using has_eval =
  // decltype(std::declval<T&>().Eval(std::declval<const vector_t&, const int>()));

public:
  Function() = default;
  ~Function() = default;
  // Computes the value of a function.
  scalar_t operator()(const vector_t &x) {
    return static_cast<DerivedFunction*>(this)->operator()(x);
  }
  // gradient
  void Gradient(const vector_t &x, vector_t *grad) {
    // static_assert(std::experimental::is_detected<void, has_grad, this>);
    return static_cast<DerivedFunction*>(this)->Gradient(x, grad);
  }
  // hessian
  void Hessian(const vector_t &x,  hessian_t *hessian) {
    return static_cast<DerivedFunction*>(this)->Hessian(x, hessian);
  }
  // evaluate entire state
  State<scalar_t, vector_t, hessian_t> Eval(const vector_t &x,
                                            const int order = 2) {
    // static_assert(std::experimental::is_detected<
    //   State<scalar_t, vector_t, hessian_t>, has_eval, this>);
    return static_cast<DerivedFunction*>(this)->Eval(x, order);
  }
};

template <typename TScalar, int TDim = Eigen::Dynamic>
class DefaultFunction : public Function<TScalar, DefaultFunction<TScalar, TDim>, TDim> {

  using scalar_t = TScalar;
  using vector_t = Eigen::Matrix<TScalar, TDim, 1>;
  using hessian_t = Eigen::Matrix<TScalar, TDim, TDim>;
  using matrix_t = Eigen::Matrix<TScalar, Eigen::Dynamic, Eigen::Dynamic>;
  using index_t = typename vector_t::Index;

  scalar_t operator()(const vector_t &x) = 0;
  // Computes the gradient of a function.
  void Gradient(const vector_t &x, vector_t *grad) {
    utils::ComputeFiniteGradient(*this, x, grad);
  }
  // Computes the Hessian of a function.
  void Hessian(const vector_t &x, hessian_t *hessian) {
    utils::ComputeFiniteHessian(*this, x, hessian);
  }
  // For improved performance, this function will return the state directly.
  // Override this method if you can compute the objective value, gradient and
  // Hessian simultaneously.
  State<scalar_t, vector_t, hessian_t> Eval(const vector_t &x,
                                            const int order = 2) { //const
    State<scalar_t, vector_t, hessian_t> state(x.rows(), order);
    state.value = this->operator()(x);
    state.x = x;
    if (order >= 1) {
      this->Gradient(x, &state.gradient);
    }
    if ((order >= 2) && (state.hessian)) {
      this->Hessian(x, &*(state.hessian));
    }
    return state;
  }
};

}  // namespace cppoptlib::function

#endif  // INCLUDE_CPPOPTLIB_FUNCTION_H_
