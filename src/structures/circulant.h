#ifndef CIRCULANT
#define CIRCULANT

#include <iostream>
#include <vector>
#include <cstddef>
#include <array>

template <typename T> struct circulant{
private:
  std::vector<T> data;
  size_t circle_index, size_;
public:  
  circulant<T>(const size_t size) {
    this->data = std::vector<T>(size);
    this->size_ = size;
    this->circle_index = 0;
  };
  circulant<T>(const std::vector<T> &x) {
    this->data = x;
    this->size_ = x.size();
    this->circle_index = 0;
  }
  T& operator [] (const size_t i){
    return this->data[(this->circle_index + i) % this->size_];
  }
  void push_back(const T item) {
    this->data[this->circle_index] = item;
    this->circle_index = (this->circle_index + 1) % this->size_;
  }
  void print() {
    for(size_t j = 0; j < this->size_; j++) {
      std::cout << this->data[j] << ", ";
    }
    std::cout << std::endl;
  }
  void print(const size_t i) {
    std::cout << "Index " << i << " value: " << this->data[(this->circle_index + i) % this->size_] << std::endl;
  }
  void print_in_order() {
    for(size_t i=0; i < this->size_; i++) {
      std::cout << this->data[(this->circle_index + i) % this->size_] << ", ";
    }
    std::cout << std::endl;
  }
  const size_t size() const {
    return this->size_;
  }
  const T& back() const {
    return this->data[(this->circle_index - 1) % this->size_];
  }
};

template <typename T, const size_t size> struct fixed_circulant{
private:
  std::array<T, size> data = std::array<T,size>();
  size_t circle_index = 0, size_ = size;
public:
  fixed_circulant<T, size>(circulant<T> &x) {
    for(size_t i= 0; i < x.size(); i++) {
      this->data[i] = x[i];
    }
    this->circle_index = 0;
    this->size_ = x.size();
  }
  fixed_circulant<T, size>(const size_t current_size) {
    this->size_ = current_size;
  }
  T& operator [] (const size_t i){
    return this->data[(this->circle_index + i) % size_];
  }
  void push_back(const T item) {
    this->data[this->circle_index] = item;
    this->circle_index = (this->circle_index + 1) % size_;
  }
  void print() {
    for(size_t j = 0; j < size_; j++) {
      std::cout << this->data[j] << ", ";
    }
    std::cout << std::endl;
  }
  void print(const size_t i) {
    std::cout << "Index " << i << " value: " << this->data[(this->circle_index + i) % size_] << std::endl;
  }
};


#endif
