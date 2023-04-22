#ifndef MAP_NAN
#define MAP_NAN

#include <cstdint>
#include <cstring>

double map_inf_and_nan_to_zero(const double orig_val) {
  auto d = orig_val;
  auto u64 = std::uint64_t{0};

  // Compile time sanity check:
  static_assert(sizeof(uint64_t) == sizeof(double));

  // Expose the bitwise representation of the value.
  std::memcpy((void *)&u64, (void *)&d, sizeof(uint64_t));

  // Extract exponent.
  const auto exp = u64 & 0x7FF0'0000'0000'0000;

  // Add 1 to the exponent and extract the carry out.
  const auto carry = ((exp + (uint64_t{1} << 52)) >> 63) & 1;

  // Produce a mask of 0xFFFF'FFFF'FFFF'FFFF if the carry was 0,
  // or 0x0000'0000'0000'0000 if the carry was 1.
  const auto mask = ~(-carry);

  // Apply the bitmask to the uint64_t representation of the
  // original value. If it was Inf or NaN, it will get forced to
  // +0. Otherwise, the original value survives.
  u64 &= mask;

  // Move the bits back into a double.
  std::memcpy((void *)&d, (void *)&u64, sizeof(double));
  // Return the resulting double.
  return d;
}

float map_inf_and_nan_to_zero(const float orig_val) {
  auto d = orig_val;
  auto u32 = std::uint32_t{0};

  // Compile time sanity check:
  static_assert(sizeof(uint32_t) == sizeof(float));

  // Expose the bitwise representation of the value.
  std::memcpy((void *)&u32, (void *)&d, sizeof(uint32_t));

  // Extract exponent.
  const auto exp = u32 & 0x7FF0'0000'0000'0000;

  // Add 1 to the exponent and extract the carry out.
  const auto carry = ((exp + (uint32_t{1} << 20)) >> 31) & 1;

  // Produce a mask of 0xFFFF'FFFF'FFFF'FFFF if the carry was 0,
  // or 0x0000'0000'0000'0000 if the carry was 1.
  const auto mask = ~(-carry);

  // Apply the bitmask to the uint64_t representation of the
  // original value. If it was Inf or NaN, it will get forced to
  // +0. Otherwise, the original value survives.
  u32 &= mask;

  // Move the bits back into a double.
  std::memcpy((void *)&d, (void *)&u32, sizeof(float));
  // Return the resulting double.
  return d;
}


#endif

