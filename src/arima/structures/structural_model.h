#ifndef STRUCT_MODEL
#define STRUCT_MODEL

// Kalman filtering model structure
template <typename U = float> struct structural_model {
  structural_model<U>() {}
  structural_model<U>(std::vector<U> phi, std::vector<U> theta,
                      std::vector<U> delta, std::vector<U> Z, std::vector<U> a,
                      std::vector<U> P, std::vector<U> T, std::vector<U> V,
                      std::vector<U> Pn, U h)
      : phi(std::move(phi)), theta(std::move(theta)), delta(std::move(delta)),
        Z(std::move(Z)), a(std::move(a)), P(std::move(P)), T(std::move(T)),
        V(std::move(V)), Pn(std::move(Pn)), h(h){};
  structural_model<U>(std::vector<U> &&phi, std::vector<U> &&theta,
                      std::vector<U> &&delta, std::vector<U> &&Z,
                      std::vector<U> &&a, std::vector<U> &&P,
                      std::vector<U> &&T, std::vector<U> &&V,
                      std::vector<U> &&Pn, U h)
      : phi(std::move(phi)), theta(std::move(theta)), delta(std::move(delta)),
        Z(std::move(Z)), a(std::move(a)), P(std::move(P)), T(std::move(T)),
        V(std::move(V)), Pn(std::move(Pn)), h(h) {};
  void set(structural_model<U> &model) {
    this->phi = std::move(model.phi);
    this->theta = std::move(model.theta);
    this->delta = std::move(model.delta);
    this->Z = std::move(model.Z);
    this->a = std::move(model.a);
    this->P = std::move(model.P);
    this->T = std::move(model.T);
    this->V = std::move(model.V);
    this->Pn = std::move(model.Pn);
    this->h = model.h;
  }
  // private:
  std::vector<U> phi, theta, delta, Z, a, P, T, V, Pn;
  U h;
};

#endif
