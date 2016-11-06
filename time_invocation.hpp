#pragma once

#include <chrono>
#include <cstddef>
#include <utility>

template<class Duration, class Clock, class Function, class... Args>
Duration time_invocation_in(const Clock& clock, std::size_t num_trials, Function&& f, Args&&... args)
{
  auto start = clock.now();
  for(std::size_t i = 0;
      i < num_trials;
      ++i)
  {
    std::forward<Function>(f)(std::forward<Args&&>(args)...);
  }
  auto end = clock.now();

  // return mean duration
  return std::chrono::duration_cast<Duration>(end - start) / num_trials;
}

template<class Function, class... Args>
std::size_t time_invocation_in_nanoseconds(std::size_t num_trials, Function&& f, Args&&... args)
{
  return ::time_invocation_in<std::chrono::nanoseconds>(std::chrono::high_resolution_clock(), num_trials, std::forward<Function>(f), std::forward<Args>(args)...).count();
}

template<class Function, class... Args>
std::size_t time_invocation_in_microseconds(std::size_t num_trials, Function&& f, Args&&... args)
{
  return ::time_invocation_in<std::chrono::microseconds>(std::chrono::high_resolution_clock(), num_trials, std::forward<Function>(f), std::forward<Args>(args)...).count();
}

template<class Function, class... Args>
std::size_t time_invocation_in_milliseconds(std::size_t num_trials, Function&& f, Args&&... args)
{
  return ::time_invocation_in<std::chrono::milliseconds>(std::chrono::high_resolution_clock(), num_trials, std::forward<Function>(f), std::forward<Args>(args)...).count();
}

template<class Function, class... Args>
std::size_t time_invocation_in_seconds(std::size_t num_trials, Function&& f, Args&&... args)
{
  return ::time_invocation_in<std::chrono::seconds>(std::chrono::system_clock(), num_trials, std::forward<Function>(f), std::forward<Args>(args)...).count();
}

