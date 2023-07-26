#ifndef BLAZE_DECOMPOSE_RESULT
#define BLAZE_DECOMPOSE_RESULT

template <typename scalar_t> struct DecompositionResult{
  std::vector<scalar_t> trend, season, remainder;
  DecompositionResult<scalar_t>(
    const std::vector<scalar_t> &trend, 
    const std::vector<scalar_t> &season,
    const std::vector<scalar_t> &remainder) : trend(trend), season(season),
    remainder(remainder) {}
};

template <typename scalar_t> struct MultipleDecompositionResult {
  std::vector<scalar_t> trend, remainder;
  std::vector<std::vector<scalar_t>> seasonalities;
  MultipleDecompositionResult<scalar_t>(
    const std::vector<scalar_t> &trend, 
    const std::vector<std::vector<scalar_t>> &seasonalities,
    const std::vector<scalar_t> &remainder) : trend(trend), 
    seasonalities(seasonalities), remainder(remainder) {}
};

#endif
