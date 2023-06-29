#ifndef MEDIAN_IMPL
#define MEDIAN_IMPL

#include <queue>
#include <algorithm>

template <typename T> T median(const std::vector<T> &x) {
  // take copy of x as sort_heap would be destructive otherwise
  std::vector<T> temp = x;
  const size_t n = temp.size();
  // construct heap in O(3N) time 
  std::make_heap(temp.begin(), temp.end());
  // sort heap in O(2Nlog(N)) time 
  std::sort_heap(temp.begin(), temp.end());
  // figure out if we are taking the mean of the two values in the middle 
  if(n % 2) {
    return (temp[n/2] + temp[(n-1)/2])/2;
  }
  // otherwise return the value in the middle 
  return temp[n/2];
}

// streaming median requires two medians 
template <typename T> struct StreamingMedian{
  std::priority_queue<T, std::vector<T>, std::greater<T>> right;
  std::priority_queue<T, std::vector<T>, std::less<T>> left;
public:
  StreamingMedian<T>() {
    this->right = std::priority_queue<T, std::vector<T>, std::greater<T>>();
    this->left = std::priority_queue<T, std::vector<T>, std::less<T>>();
  }
  void push_back(const T x) {
    // push onto left heap 
    this->left.push(x);
    // periodically call reorder 
    if(this->left.size() > (this->right.size() + 15)) {
      reorder();
    }
  }
  const T value() const {
    reorder();
    if(left.size() == right.size())
      return (this->left.top() + this->right.top())/2;
    // otherwise I know the left heap holds the median
    return this->left.top();
  }
private:
  // reorder elements between heaps 
  void reorder() {
    // this moves elements from the left heap to the right heap 
    // if left + right is even, we will take an average of two tops 
    // so we need to do the correct number of reorderings 
    // since we only ever push to left heap, this is really simple
    if((this->left.size() + this->right.size()) % 2) {
      while(this->left.size() != this->right.size()) {
        this->right.push(this->left.pop());
      }
      return;
    }
    // otherwise we need to do fewer pop-push pairs :) 
    while(this->left.size() != (this->right.size()+1)) {
      this->right.push(this->left.pop());
    }
  }
};

template <typename T> class FudgyStreamingMedian{
  T average, median, eta;
public:
  FudgyStreamingMedian<T>(const T eta = 0.1) {
    this->average = 0;
    this->median = 0;
    this->eta = eta;
  }
  void update(const T x) {
      this->average += (x - average) * eta; // rough running average.
      this->median += std::copysign(this->average * eta, x - median);
  }
  void reset() {
    this->average = 0;
    this->median = 0;
  }
  const T value() const {
    return this->median; 
  }
};

#endif
