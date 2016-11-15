#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include <type_traits>

template<class T, class Bounder>
struct memoized_bounder
{
  // we're memoizing the result of Bounder, which returns a bounding_box_type
  using bounding_box_type = std::result_of_t<Bounder(const T&)>;

  // avoid copying this thing unintentionally
  memoized_bounder(memoized_bounder&&) = delete;

  template<class ContiguousRange>
  memoized_bounder(const ContiguousRange& elements, Bounder bounder)
    : elements_(&*elements.begin())
  {
    bounding_boxes_.reserve(elements.size());
    std::transform(elements.begin(), elements.end(), std::back_inserter(bounding_boxes_), bounder);
  }

  bounding_box_type operator()(const T& element) const
  {
    size_t element_idx = &element - elements_;
    return bounding_boxes_[element_idx];
  }

  const T* elements_;
  std::vector<bounding_box_type> bounding_boxes_;
};

