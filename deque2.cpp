#include <algorithm>
#include <cstddef>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

template <typename T, typename Allocator = std::allocator<T>>
class Deque {
 public:
  using value_type = T;
  using allocator_type = Allocator;
  using alloc_traits = std::allocator_traits<allocator_type>;
  Deque();
  Deque(const Allocator& alloc);
  Deque(const Deque& deq);
  Deque(Deque&& deq);
  Deque(size_t count, const Allocator& alloc = Allocator());
  Deque(size_t count, const T& value, const Allocator& alloc = Allocator());
  Deque(std::initializer_list<T> init, const Allocator& alloc = Allocator());
  ~Deque();
  Deque<T, Allocator>& operator=(const Deque& deq);
  Deque<T, Allocator>& operator=(Deque&& other);
  template <typename... Args>
  void create_one_with_args(Args&&... args);
  void create_one_with_move(T&& value);
  allocator_type get_allocator() const;
  void reserve(Deque& newdeq);
  void swap(Deque& other);
  void swap_deque(Deque& newdeq);
  void move_right(size_t index_val);
  size_t size() const;
  bool empty() const;
  T& operator[](size_t num);
  const T& operator[](size_t num) const;
  T& at(size_t num);
  const T& at(size_t num) const;
  void push_back(const T& val);
  void push_back(T&& val);
  void pop_back();
  template <typename... Args>
  void emplace_back(Args&&... args);
  void push_front(const T& val);
  void push_front(T&& val);
  template <typename... Args>
  void emplace_front(Args&&... args);
  void pop_front();

  template <bool IsConst>
  class CommonIterator
      : public std::iterator<std::random_access_iterator_tag,
                             std::conditional_t<IsConst, const T, T>> {
   public:
    using ref = std::conditional_t<IsConst, const T&, T&>;
    using point = std::conditional_t<IsConst, const T*, T*>;
    CommonIterator();
    CommonIterator(
        point ptr, size_t shift_ptr,
        std::conditional_t<IsConst, const std::vector<T*>&, std::vector<T*>&>
            deq);
    CommonIterator(const CommonIterator& itt);
    CommonIterator& operator=(const CommonIterator& itt);
    CommonIterator operator++(int val);
    CommonIterator& operator++();
    CommonIterator operator--(int val);
    CommonIterator& operator--();
    CommonIterator operator+(int val) const;
    CommonIterator& operator+=(int val);
    CommonIterator operator-(int val) const;
    CommonIterator& operator-=(int val);
    bool operator<(const CommonIterator& itt) const;
    bool operator==(const CommonIterator& itt) const;
    bool operator<=(const CommonIterator& itt) const;
    bool operator>=(const CommonIterator& itt) const;
    bool operator!=(const CommonIterator& itt) const;
    bool operator>(const CommonIterator& itt) const;
    std::ptrdiff_t operator-(const CommonIterator& itt) const;
    ref operator*() const;
    point operator->() const;
    size_t get_shift() const;

   private:
    std::conditional_t<IsConst, const T*, T*> ptr_;
    size_t shift_ptr_;
    std::vector<T*> deq_;
  };

  using iterator = CommonIterator<false>;
  using const_iterator = CommonIterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin();
  const_iterator begin() const;
  const_iterator cbegin() const;
  iterator end();
  const_iterator end() const;
  const_iterator cend() const;
  reverse_iterator rbegin();
  reverse_iterator rend();
  const_reverse_iterator crbegin() const;
  const_reverse_iterator crend() const;
  void insert(iterator itt, const T& value);
  template <typename... Args>
  void emplace(iterator itt, Args&&... args);
  void erase(iterator itt);

 private:
  std::vector<T*> deq_;
  size_t size_;
  size_t shift_;
  allocator_type alloc_;
};

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque() : size_(0), shift_(0) {}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(const Allocator& alloc) : alloc_(alloc) {}

template <typename T, typename Allocator>
Deque<T, Allocator>::~Deque() {
  const size_t kCol = 6;
  for (size_t i = shift_; i < shift_ + size_; ++i) {
    alloc_traits::destroy(alloc_, &deq_[i / kCol][i % kCol]);
  }
  for (size_t i = 0; i < deq_.size(); ++i) {
    alloc_traits::deallocate(alloc_, deq_[i], kCol);
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(size_t count, const Allocator& alloc)
    : size_(count), alloc_(alloc) {
  const size_t kCol = 6;
  size_ = count;
  if (count % kCol != 0) {
    count += kCol - (count % kCol);
  }
  size_t capacity = count * 3;
  deq_.resize(capacity / kCol, nullptr);
  for (size_t i = 0; i < capacity / kCol; ++i) {
    deq_[i] = alloc_traits::allocate(alloc_, kCol);
  }
  shift_ = (capacity / 3);
  for (size_t i = shift_; i < shift_ + size_; ++i) {
    try {
      alloc_traits::construct(alloc_, &deq_[i / kCol][(i % kCol)]);
    } catch (...) {
      for (size_t j = shift_; j < i; ++j) {
        alloc_traits::destroy(alloc_, &deq_[j / kCol][j % kCol]);
      }
      for (size_t i = 0; i < deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kCol);
      }
      size_ = 0;
      shift_ = 0;
      throw;
    }
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(size_t count, const T& value, const Allocator& alloc)
    : size_(count), alloc_(alloc) {
  const size_t kCol = 6;
  size_ = count;
  if (count % kCol != 0) {
    count += kCol - (count % kCol);
  }
  size_t capacity = count * 3;
  deq_.resize(capacity / kCol, nullptr);
  for (size_t i = 0; i < capacity / kCol; ++i) {
    deq_[i] = alloc_traits::allocate(alloc_, kCol);
  }
  shift_ = (capacity / 3);
  for (size_t i = shift_; i < shift_ + size_; ++i) {
    try {
      alloc_traits::construct(alloc_, &deq_[i / kCol][(i % kCol)], value);
    } catch (...) {
      for (size_t j = shift_; j < i; ++j) {
        alloc_traits::destroy(alloc_, &deq_[j / kCol][j % kCol]);
      }
      for (size_t i = 0; i < deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kCol);
      }
      size_ = 0;
      shift_ = 0;
      throw;
    }
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(const Deque<T, Allocator>& deq)
    : size_(deq.size_), shift_(deq.shift_) {
  size_t count = deq.size_;
  alloc_ = alloc_traits::select_on_container_copy_construction(deq.alloc_);
  const size_t kCol = 6;
  if (count % kCol != 0) {
    count += kCol - (count % kCol);
  }
  size_t capacity = count * 3;
  deq_.resize(capacity / kCol, nullptr);
  for (size_t i = 0; i < capacity / kCol; ++i) {
    deq_[i] = alloc_traits::allocate(alloc_, kCol);
  }
  shift_ = (capacity / 3);
  for (size_t i = shift_; i < shift_ + size_; ++i) {
    try {
      alloc_traits::construct(alloc_, &deq_[i / kCol][(i % kCol)],
                              deq.deq_[i / kCol][i % kCol]);
    } catch (...) {
      for (size_t j = shift_; j < i; ++j) {
        alloc_traits::destroy(alloc_, &deq_[j / kCol][j % kCol]);
      }
      for (size_t i = 0; i < deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kCol);
      }
      size_ = 0;
      shift_ = 0;
      throw;
    }
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(Deque<T, Allocator>&& deq)
    : size_(deq.size_), shift_(deq.shift_), deq_(deq.deq_), alloc_(deq.alloc_) {
  deq.size_ = 0;
  deq.shift_ = 0;
  for (size_t i = 0; i < deq.deq_.size(); ++i) {
    deq.deq_[i] = nullptr;
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>::Deque(std::initializer_list<T> init,
                           const Allocator& alloc)
    : size_(init.size()), alloc_(alloc) {
  const size_t kCol = 6;
  size_t count = init.size();
  if (count % kCol != 0) {
    count += kCol - (count % kCol);
  }
  size_t capacity = count * 3;
  deq_.resize(capacity / kCol, nullptr);
  for (size_t i = 0; i < capacity / kCol; ++i) {
    deq_[i] = alloc_traits::allocate(alloc_, kCol);
  }
  shift_ = (capacity / 3);
  size_t index = shift_;
  for (auto itt = init.begin(); itt != init.end(); ++itt, ++index) {
    try {
      alloc_traits::construct(alloc_, &deq_[index / kCol][(index % kCol)],
                              *itt);
    } catch (...) {
      for (size_t i = shift_; i < index; ++i) {
        alloc_traits::destroy(alloc_, &deq_[i / kCol][i % kCol]);
      }
      for (size_t i = 0; i < deq_.size(); ++i) {
        alloc_traits::deallocate(alloc_, deq_[i], kCol);
      }
      size_ = 0;
      shift_ = 0;
      throw;
    }
  }
}

template <typename T, typename Allocator>
Deque<T, Allocator>& Deque<T, Allocator>::operator=(
    const Deque<T, Allocator>& deq) {
  if (this == &deq) {
    return *this;
  }
  Deque<T, Allocator> newdeq;
  size_t count = deq.size_;
  newdeq.size_ = deq.size_;
  newdeq.shift_ = deq.shift_;
  newdeq.alloc_ = alloc_;
  if (alloc_traits::propagate_on_container_copy_assignment::value) {
    newdeq.alloc_ = deq.alloc_;
  }
  const size_t kCol = 6;
  if (count % kCol != 0) {
    count += kCol - (count % kCol);
  }
  size_t capacity = count * 3;
  newdeq.deq_.resize(capacity / kCol, nullptr);
  for (size_t i = 0; i < capacity / kCol; ++i) {
    newdeq.deq_[i] = alloc_traits::allocate(newdeq.alloc_, kCol);
  }
  newdeq.shift_ = (capacity / 3);
  for (size_t i = newdeq.shift_; i < newdeq.shift_ + newdeq.size_; ++i) {
    try {
      alloc_traits::construct(newdeq.alloc_, newdeq.deq_[i / kCol] + (i % kCol),
                              deq.deq_[i / kCol][i % kCol]);
    } catch (...) {
      for (size_t j = newdeq.shift_; j < i; ++j) {
        alloc_traits::destroy(newdeq.alloc_, &newdeq.deq_[j / kCol][j % kCol]);
      }
      for (size_t i = 0; i < newdeq.deq_.size(); ++i) {
        alloc_traits::deallocate(newdeq.alloc_, newdeq.deq_[i], kCol);
      }
      newdeq.deq_.clear();
      newdeq.size_ = 0;
      newdeq.shift_ = 0;
      throw;
    }
  }
  if (!empty()) {
    for (size_t i = shift_; i < shift_ + size_; ++i) {
      alloc_traits::destroy(alloc_, &deq_[i / kCol][(i % kCol)]);
    }
    for (size_t i = 0; i < deq_.size(); ++i) {
      alloc_traits::deallocate(alloc_, deq_[i], kCol);
    }
  }
  size_ = 0;
  shift_ = 0;
  deq_.clear();
  swap(newdeq);
  return *this;
}

template <typename T, typename Allocator>
Deque<T, Allocator>& Deque<T, Allocator>::operator=(
    Deque<T, Allocator>&& other) {
  Deque<T, Allocator> newdeq = std::move(other);
  swap(newdeq);
  return *this;
}

template <typename T, typename Allocator>
template <typename... Args>
void Deque<T, Allocator>::create_one_with_args(Args&&... args) {
  const size_t kCol = 6;
  size_t count = size_ = 1;
  if (count % kCol != 0) {
    count += kCol - (count % kCol);
  }
  size_t capacity = count * 3;
  deq_.resize(capacity / kCol, nullptr);
  for (size_t i = 0; i < capacity / kCol; ++i) {
    deq_[i] = alloc_traits::allocate(alloc_, kCol);
  }
  shift_ = (capacity / 3);
  try {
    alloc_traits::construct(alloc_, &deq_[shift_ / kCol][(shift_ % kCol)],
                            std::forward<Args>(args)...);
  } catch (...) {
    for (size_t i = 0; i < deq_.size(); ++i) {
      alloc_traits::deallocate(alloc_, deq_[i], kCol);
    }
    throw;
  }
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::allocator_type
Deque<T, Allocator>::get_allocator() const {
  return alloc_;
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::reserve(Deque<T, Allocator>& newdeq) {
  const size_t kCol = 6;
  newdeq.alloc_ = alloc_;
  shift_ += (deq_.size() * kCol);
  for (size_t i = 0; i < newdeq.deq_.size(); ++i) {
    if (i >= deq_.size() && i < deq_.size() * 2) {
      newdeq.deq_[i] = deq_[i - deq_.size()];
    } else {
      newdeq.deq_[i] = alloc_traits::allocate(newdeq.alloc_, kCol);
    }
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::create_one_with_move(T&& value) {
  const size_t kCol = 6;
  size_t count = size_ = 1;
  if (count % kCol != 0) {
    count += kCol - (count % kCol);
  }
  size_t capacity = count * 3;
  deq_.resize(capacity / kCol, nullptr);
  for (size_t i = 0; i < capacity / kCol; ++i) {
    deq_[i] = alloc_traits::allocate(alloc_, kCol);
  }
  shift_ = (capacity / 3);
  try {
    alloc_traits::construct(alloc_, &deq_[shift_ / kCol][(shift_ % kCol)],
                            std::move(value));
  } catch (...) {
    for (size_t i = 0; i < deq_.size(); ++i) {
      alloc_traits::deallocate(alloc_, deq_[i], kCol);
    }
    throw;
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::swap(Deque<T, Allocator>& other) {
  std::swap(size_, other.size_);
  std::swap(shift_, other.shift_);
  std::swap(deq_, other.deq_);
  std::swap(alloc_, other.alloc_);
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::swap_deque(Deque<T, Allocator>& newdeq) {
  for (size_t i = 0; i < deq_.size(); ++i) {
    deq_[i] = nullptr;
  }
  deq_.resize(newdeq.deq_.size(), nullptr);
  for (size_t i = 0; i < deq_.size(); ++i) {
    deq_[i] = newdeq.deq_[i];
    newdeq.deq_[i] = nullptr;
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::move_right(size_t index_val) {
  const size_t kCol = 6;
  size_t index = shift_ + size_;
  index_val += shift_;
  std::swap(deq_[index_val / kCol][index_val % kCol],
            deq_[(index - 1) / kCol][(index - 1) % kCol]);
  for (size_t i = index - 1; i > index_val + 1; --i) {
    deq_[i / kCol][i % kCol] = deq_[(i - 1) / kCol][(i - 1) % kCol];
  }
}

template <typename T, typename Allocator>
size_t Deque<T, Allocator>::size() const {
  return size_;
}

template <typename T, typename Allocator>
bool Deque<T, Allocator>::empty() const {
  return (size_ == 0);
}

template <typename T, typename Allocator>
T& Deque<T, Allocator>::operator[](size_t num) {
  const size_t kCol = 6;
  num += shift_;
  return deq_[num / kCol][num % kCol];
}

template <typename T, typename Allocator>
const T& Deque<T, Allocator>::operator[](size_t num) const {
  const size_t kCol = 6;
  num += shift_;
  return deq_[num / kCol][num % kCol];
}

template <typename T, typename Allocator>
T& Deque<T, Allocator>::at(size_t num) {
  if (num >= size_) {
    const size_t kCol = 6;
    throw std::out_of_range("at_out_of_range");
  }
  const size_t kCol = 6;
  num += shift_;
  return deq_[num / kCol][num % kCol];
}

template <typename T, typename Allocator>
const T& Deque<T, Allocator>::at(size_t num) const {
  if (num >= size_) {
    throw std::out_of_range("at_out_of_range");
  }
  const size_t kCol = 6;
  num += shift_;
  return deq_[num / kCol][num % kCol];
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_back(const T& val) {
  emplace_back(val);
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_back(T&& val) {
  emplace_back(std::move(val));
}

template <typename T, typename Allocator>
template <typename... Args>
void Deque<T, Allocator>::emplace_back(Args&&... args) {
  const size_t kCol = 6;
  if (empty()) {
    try {
      create_one_with_args(std::forward<Args>(args)...);
    } catch (...) {
      throw;
    }
  } else {
    if (shift_ + size_ + 1 == (deq_.size() * kCol)) {
      Deque<T, Allocator> newdeq;
      newdeq.deq_.resize(deq_.size() * 3, nullptr);
      reserve(newdeq);
      size_t index = shift_ + size_;
      try {
        alloc_traits::construct(newdeq.alloc_,
                                &newdeq.deq_[index / kCol][(index % kCol)],
                                std::forward<Args>(args)...);
      } catch (...) {
        for (size_t i = 0; i < newdeq.deq_.size(); ++i) {
          if (!(i >= deq_.size() && i < deq_.size() * 2)) {
            alloc_traits::deallocate(newdeq.alloc_, newdeq.deq_[i], kCol);
          }
        }
        throw;
      }
      ++size_;
      swap_deque(newdeq);
    } else {
      size_t index = shift_ + size_;
      try {
        alloc_traits::construct(alloc_, &deq_[index / kCol][(index % kCol)],
                                std::forward<Args>(args)...);
      } catch (...) {
        throw;
      }
      ++size_;
    }
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::pop_back() {
  const size_t kCol = 6;
  size_t index = shift_ + size_ - 1;
  alloc_traits::destroy(alloc_, deq_[index / kCol] + (index % kCol));
  --size_;
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_front(const T& val) {
  emplace_front(val);
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::push_front(T&& val) {
  emplace_front(std::move(val));
}

template <typename T, typename Allocator>
template <typename... Args>
void Deque<T, Allocator>::emplace_front(Args&&... args) {
  const size_t kCol = 6;
  if (empty()) {
    try {
      create_one_with_args(std::forward<Args>(args)...);
    } catch (...) {
      throw;
    }
  } else {
    if (shift_ == 0) {
      Deque<T, Allocator> newdeq;
      newdeq.deq_.resize(deq_.size() * 3, nullptr);
      reserve(newdeq);
      try {
        --shift_;
        alloc_traits::construct(newdeq.alloc_,
                                &newdeq.deq_[shift_ / kCol][(shift_ % kCol)],
                                std::forward<Args>(args)...);
      } catch (...) {
        for (size_t i = 0; i < newdeq.deq_.size(); ++i) {
          if (!(i >= deq_.size() && i < deq_.size() * 2)) {
            alloc_traits::deallocate(newdeq.alloc_, newdeq.deq_[i], kCol);
          }
        }
        throw;
      }
      ++size_;
      swap_deque(newdeq);
    } else {
      --shift_;
      try {
        alloc_traits::construct(alloc_, &deq_[shift_ / kCol][(shift_ % kCol)],
                                std::forward<Args>(args)...);
      } catch (...) {
        ++shift_;
        throw;
      }
      ++size_;
    }
  }
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::pop_front() {
  const size_t kCol = 6;
  alloc_traits::destroy(alloc_, &deq_[shift_ / kCol][(shift_ % kCol)]);
  ++shift_;
  --size_;
}

template <typename T, typename Allocator>
template <bool IsConst>
Deque<T, Allocator>::CommonIterator<IsConst>::CommonIterator()
    : ptr_(nullptr), shift_ptr_(0) {}

template <typename T, typename Allocator>
template <bool IsConst>
Deque<T, Allocator>::CommonIterator<IsConst>::CommonIterator(
    point ptr, size_t shift_ptr,
    std::conditional_t<IsConst, const std::vector<T*>&, std::vector<T*>&> deq) {
  ptr_ = ptr;
  shift_ptr_ = shift_ptr;
  deq_ = deq;
}

template <typename T, typename Allocator>
template <bool IsConst>
Deque<T, Allocator>::CommonIterator<IsConst>::CommonIterator(
    const CommonIterator& itt) {
  ptr_ = itt.ptr_;
  shift_ptr_ = itt.shift_ptr_;
  deq_ = itt.deq_;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>&
Deque<T, Allocator>::CommonIterator<IsConst>::operator=(
    const CommonIterator& itt) {
  ptr_ = itt.ptr_;
  shift_ptr_ = itt.shift_ptr_;
  deq_ = itt.deq_;
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>
Deque<T, Allocator>::CommonIterator<IsConst>::operator++(int val) {
  if (ptr_ == nullptr) {
    return *this;
  }
  val = 0;
  iterator newitt = *this;
  (*this) += val + 1;
  return newitt;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>&
Deque<T, Allocator>::CommonIterator<IsConst>::operator++() {
  if (ptr_ == nullptr) {
    return *this;
  }
  const size_t kCol = 6;
  ++shift_ptr_;
  ptr_ = &deq_[shift_ptr_ / kCol][shift_ptr_ % kCol];
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>
Deque<T, Allocator>::CommonIterator<IsConst>::operator--(int val) {
  if (ptr_ == nullptr) {
    return *this;
  }
  val = 0;
  CommonIterator<IsConst> newitt = *this;
  (*this) -= val + 1;
  return newitt;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>&
Deque<T, Allocator>::CommonIterator<IsConst>::operator--() {
  if (ptr_ == nullptr) {
    return *this;
  }
  const size_t kCol = 6;
  --shift_ptr_;
  ptr_ = &deq_[shift_ptr_ / kCol][shift_ptr_ % kCol];
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>
Deque<T, Allocator>::CommonIterator<IsConst>::operator+(int val) const {
  if (ptr_ == nullptr) {
    return *this;
  }
  CommonIterator<IsConst> newitt = *this;
  newitt += val;
  return newitt;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>&
Deque<T, Allocator>::CommonIterator<IsConst>::operator+=(int val) {
  if (ptr_ == nullptr) {
    return *this;
  }
  const size_t kCol = 6;
  shift_ptr_ += val;
  ptr_ = &deq_[shift_ptr_ / kCol][shift_ptr_ % kCol];
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>
Deque<T, Allocator>::CommonIterator<IsConst>::operator-(int val) const {
  if (ptr_ == nullptr) {
    return *this;
  }
  CommonIterator<IsConst> newitt = *this;
  newitt -= val;
  return newitt;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>&
Deque<T, Allocator>::CommonIterator<IsConst>::operator-=(int val) {
  if (ptr_ == nullptr) {
    return *this;
  }
  const size_t kCol = 6;
  shift_ptr_ -= val;
  ptr_ = &deq_[shift_ptr_ / kCol][shift_ptr_ % kCol];
  return *this;
}

template <typename T, typename Allocator>
template <bool IsConst>
bool Deque<T, Allocator>::CommonIterator<IsConst>::operator<(
    const CommonIterator& itt) const {
  return (ptr_ < itt.ptr_);
}

template <typename T, typename Allocator>
template <bool IsConst>
bool Deque<T, Allocator>::CommonIterator<IsConst>::operator==(
    const CommonIterator& itt) const {
  return (ptr_ == itt.ptr_);
}

template <typename T, typename Allocator>
template <bool IsConst>
bool Deque<T, Allocator>::CommonIterator<IsConst>::operator<=(
    const CommonIterator& itt) const {
  return (*this < itt) && (*this == itt);
}

template <typename T, typename Allocator>
template <bool IsConst>
bool Deque<T, Allocator>::CommonIterator<IsConst>::operator>=(
    const CommonIterator& itt) const {
  return !(*this < itt);
}

template <typename T, typename Allocator>
template <bool IsConst>
bool Deque<T, Allocator>::CommonIterator<IsConst>::operator!=(
    const CommonIterator& itt) const {
  return !(*this == itt);
}

template <typename T, typename Allocator>
template <bool IsConst>
bool Deque<T, Allocator>::CommonIterator<IsConst>::operator>(
    const CommonIterator& itt) const {
  return (itt < *this);
}

template <typename T, typename Allocator>
template <bool IsConst>
std::ptrdiff_t Deque<T, Allocator>::CommonIterator<IsConst>::operator-(
    const CommonIterator& itt) const {
  return (shift_ptr_ - itt.shift_ptr_);
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>::ref
Deque<T, Allocator>::CommonIterator<IsConst>::operator*() const {
  return *ptr_;
}

template <typename T, typename Allocator>
template <bool IsConst>
typename Deque<T, Allocator>::template CommonIterator<IsConst>::point
Deque<T, Allocator>::CommonIterator<IsConst>::operator->() const {
  return ptr_;
}

template <typename T, typename Allocator>
template <bool IsConst>
size_t Deque<T, Allocator>::CommonIterator<IsConst>::get_shift() const {
  return shift_ptr_;
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::iterator Deque<T, Allocator>::begin() {
  const size_t kCol = 6;
  if (empty()) {
    return iterator();
  }
  return iterator(&deq_[shift_ / kCol][shift_ % kCol], shift_, deq_);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_iterator Deque<T, Allocator>::begin()
    const {
  const size_t kCol = 6;
  if (empty()) {
    return const_iterator();
  }
  return const_iterator(&deq_[shift_ / kCol][shift_ % kCol], shift_, deq_);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_iterator Deque<T, Allocator>::cbegin()
    const {
  const size_t kCol = 6;
  if (empty()) {
    return const_iterator();
  }
  return const_iterator(&deq_[shift_ / kCol][shift_ % kCol], shift_, deq_);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::iterator Deque<T, Allocator>::end() {
  const size_t kCol = 6;
  if (empty()) {
    return iterator();
  }
  return iterator(&deq_[(shift_ + size_) / kCol][(shift_ + size_) % kCol],
                  shift_ + size_, deq_);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_iterator Deque<T, Allocator>::end() const {
  const size_t kCol = 6;
  if (empty()) {
    return const_iterator();
  }
  return const_iterator(&deq_[(shift_ + size_) / kCol][(shift_ + size_) % kCol],
                        shift_ + size_, deq_);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_iterator Deque<T, Allocator>::cend() const {
  const size_t kCol = 6;
  if (empty()) {
    return const_iterator();
  }
  return const_iterator(&deq_[(shift_ + size_) / kCol][(shift_ + size_) % kCol],
                        shift_ + size_, deq_);
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::reverse_iterator Deque<T, Allocator>::rbegin() {
  const size_t kCol = 6;
  if (empty()) {
    return reverse_iterator();
  }
  return reverse_iterator(
      iterator(&deq_[(shift_ + size_) / kCol][(shift_ + size_) % kCol],
               shift_ + size_, deq_));
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::reverse_iterator Deque<T, Allocator>::rend() {
  const size_t kCol = 6;
  if (empty()) {
    return reverse_iterator();
  }
  return reverse_iterator(
      iterator(&deq_[(shift_) / kCol][(shift_) % kCol], shift_, deq_));
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_reverse_iterator
Deque<T, Allocator>::crbegin() const {
  const size_t kCol = 6;
  if (empty()) {
    return const_reverse_iterator();
  }
  return const_reverse_iterator(
      const_iterator(&deq_[(shift_ + size_) / kCol][(shift_ + size_) % kCol],
                     shift_ + size_, deq_));
}

template <typename T, typename Allocator>
typename Deque<T, Allocator>::const_reverse_iterator
Deque<T, Allocator>::crend() const {
  const size_t kCol = 6;
  if (empty()) {
    return const_reverse_iterator();
  }
  return const_reverse_iterator(
      const_iterator(&deq_[(shift_) / kCol][(shift_) % kCol], shift_, deq_));
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::insert(iterator itt, const T& value) {
  emplace(itt, value);
}

template <typename T, typename Allocator>
template <typename... Args>
void Deque<T, Allocator>::emplace(iterator itt, Args&&... args) {
  size_t index_val = itt.get_shift() - shift_;
  emplace_back(std::forward<Args>(args)...);
  move_right(index_val);
}

template <typename T, typename Allocator>
void Deque<T, Allocator>::erase(iterator itt) {
  const size_t kCol = 6;
  size_t index_val = itt.get_shift();
  size_t index = shift_ + size_;
  for (size_t i = index_val; i < index - 1; ++i) {
    deq_[i / kCol][i % kCol] = deq_[(i + 1) / kCol][(i + 1) % kCol];
  }
  alloc_traits::destroy(alloc_, &deq_[(index - 1) / kCol][(index - 1) % kCol]);
  --size_;
}

struct Accountant {
  // Some field of strange size
  char arr[40];

  static size_t ctor_calls;
  static size_t dtor_calls;

  // NO LINT
  static void reset() {
    ctor_calls = 0;
    dtor_calls = 0;
  }

  Accountant() {
    ++ctor_calls;
  }

  // NO LINT
  Accountant(const Accountant&) {
    ++ctor_calls;
  }

  // NO LINT
  Accountant& operator=(const Accountant&) {
    // Actually, when it comes to assign one list to another,
    // list can use element-wise assignment instead of destroying nodes and creating new ones
    ++ctor_calls;
    ++dtor_calls;
    return *this;
  }

  Accountant(Accountant&&) = default;
  Accountant& operator=(Accountant&&) = default;

  ~Accountant() {
    ++dtor_calls;
  }
};

struct ThrowingAccountant: public Accountant {
  static bool need_throw;

  int value = 0;

  // NO LINT
  ThrowingAccountant(int value = 0): Accountant(), value(value) {
    if (need_throw && ctor_calls % 5 == 4)
      throw std::string("Ahahahaha you have been cocknut");
  }

  // NO LINT
  ThrowingAccountant(const ThrowingAccountant& other): Accountant(), value(other.value) {
    if (need_throw && ctor_calls % 5 == 4)
      throw std::string("Ahahahaha you have been cocknut");
  }

  // NO LINT
  ThrowingAccountant& operator=(const ThrowingAccountant& other) {
    value = other.value;
    ++ctor_calls;
    ++dtor_calls;
    if (need_throw && ctor_calls % 5 == 4)
      throw std::string("Ahahahaha you have been cocknut");
    return *this;
  }

};


size_t Accountant::ctor_calls = 0;
size_t Accountant::dtor_calls = 0;

bool ThrowingAccountant::need_throw = true;

int main() {

  ThrowingAccountant::need_throw = true;
  try {
    Deque<ThrowingAccountant> d(8);
  } catch(...) {
    std::cout << (Accountant::ctor_calls == 4);
    std::cout << (Accountant::dtor_calls == 4);
  }

  return 0;
}