#ifndef RNG
#define RNG

#include <stdint.h>
#define MAX_SIZE_64_BIT_UINT (18446744073709551615U)

uint64_t bitwise_rotate(uint64_t x, int bits, int rotate_bits) {
  return (x << rotate_bits) | (x >> (bits - rotate_bits));
}

struct halton {
  halton(){
    b=2;
    y = 1;
    n=0;
    d=1;
    x = 1;
  };
  float yield(){
    x = d-n;
    if(x == 1){
      n = 1;
      d *= b;
    }
    else {
      y = d;
      while(x <= y) {
        y /= b;
        n = (b + 1) * y - x;
      }
    }
    return (float)(n/d);
  };
  void set( float base ) {
    b = base;
    y = 1;
    n = 0;
    d = 1;
    x = 1;
  };
  void reset(){
    b=2;
    y = 1;
    n = 0;
    d = 1;
    x = 1;
  };
  private:
    float b, y, n, d, x;
};


struct recurrent {
  recurrent(){
    seed = 0.5;
    alpha = 0.618034;
    z = alpha + seed;
    z -= (double)(int)(z);
  };
  recurrent( double init_seed ) {
    seed = init_seed;
    alpha = 0.618034;
    z = alpha + seed;
    z -= (double)(int)(z);
  };
  double yield() {
    z = (z+alpha);
    // a slightly evil way to do z % 1 with floats
    z -= (double)(int)(z);
    return z;
  };
  void reset(){
    alpha = 0.618034;
    seed = 0.5;
    z = 0;
  };
  private:
    double alpha = 0.618034, seed = 0.5, z = 0;
  //   float phi(int d = 1) {
  //     float x = 2.0;
  //     for( int i=0; i< 15; i++) {
  //       x = pow((x+1),1/(d+1));
  //     }
  //     return x;
  //   }
};

struct splitmix {
  splitmix(){
    s = 12374563468;
  };
  float yield() {
    uint64_t result = (s += 0x9E3779B97f4A7C15);
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    return (float)(result ^ (result >> 31))/(float)MAX_SIZE_64_BIT_UINT;
  };
  int yield_init() {
    uint64_t result = (s += 0x9E3779B97f4A7C15);
    result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
    result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
    return result ^ (result >> 31);
  };
  void set( uint64_t x ) {
    s = x;
  };
  private:
    uint64_t s;
};

struct xoshiro {
  xoshiro(){
    splitmix gn;
    s[0] = gn.yield_init();
    s[1] = s[0] >> 32;

    s[2] = gn.yield();
    s[3] = s[2] >> 32;
  };
  double yield(){
    uint64_t const result = s[0] + s[3];
    uint64_t const t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = bitwise_rotate(s[3], 64, 45);

    return (double)result/(double)MAX_SIZE_64_BIT_UINT;
  }
  void reset(){
    splitmix gn;
    s[0] = gn.yield_init();
    s[1] = s[0] >> 32;

    s[2] = gn.yield();
    s[3] = s[2] >> 32;
  };
  void set( uint64_t x,
            uint64_t y,
            uint64_t z,
            uint64_t t) {
    s[0] = x;
    s[1] = y;
    s[2] = z;
    s[3] = t;
  };
  private:
    uint64_t rol64(uint64_t x, int k)
    {
      return (x << k) | (x >> (64 - k));
    }
    uint64_t s[4];
};

struct xorshift {
  xorshift() {
    splitmix gn;
    x[0] = gn.yield_init();
    x[1] = x[0] >> 32;
  };
  double yield() {
    uint64_t t = x[0];
    uint64_t const s = x[1];
    x[0] = s;
    t ^= t << 23;		// a
    t ^= t >> 18;		// b -- Again, the shifts and the multipliers are tunable
    t ^= s ^ (s >> 5);	// c
    x[1] = t;
    return (double)(t + s)/(double)MAX_SIZE_64_BIT_UINT;
  };
  void reset(){
    splitmix gn;
    x[0] = gn.yield_init();
    x[1] = x[0] >> 32;
  };
  void set( uint64_t y, uint64_t z ) {
    x[0] = y;
    x[1] = z;
  };
  private:
    uint64_t x[2];
};

#endif
