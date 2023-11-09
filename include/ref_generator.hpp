#pragma once

#include <coroutine>
#include <exception>
#include <utility>

template<typename T>
struct generator {
  struct promise_type;
  using handle_type = std::coroutine_handle<promise_type>;

  struct promise_type {
    T value_;
    std::exception_ptr exception_;

    generator get_return_object() {
      return generator(handle_type::from_promise(*this));
    }
    std::suspend_always initial_suspend() { return {}; }
    std::suspend_always final_suspend() noexcept { return {}; }
    void unhandled_exception() { exception_ = std::current_exception(); }
    template<std::convertible_to<T> From> // C++20 concept
    std::suspend_always yield_value(From &&from) {
      value_ = std::forward<From>(from);
      return {};
    }
    void return_void() {}
  };

  handle_type h_;

  generator(handle_type h) : h_(h) {}
  generator(const generator &) = delete;
  ~generator() { h_.destroy(); }
  explicit operator bool() {
    fill();
    return !h_.done();
  }
  T operator()() {
    fill();
    full_ = false;
    return std::move(h_.promise().value_);
  }

  // begin
    struct Iter {
        handle_type h_;
        bool done_;
    
        void operator++() {
        h_.resume();
        done_ = h_.done();
        }
        T operator*() const { return h_.promise().value_; }
        bool operator==(Iter const &rhs) const {
        return done_ == rhs.done_;
        }
        bool operator!=(Iter const &rhs) const { return !(*this == rhs); }
    };
    Iter begin() {
        fill();
        return Iter{h_, h_.done()};
    }
    Iter end() { return Iter{h_, true}; }

private:
  bool full_ = false;

  void fill() {
    if (!full_) {
      h_();
      if (h_.promise().exception_)
        std::rethrow_exception(h_.promise().exception_);
      full_ = true;
    }
  }
};

