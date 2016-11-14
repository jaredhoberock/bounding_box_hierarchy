#pragma once

#include <utility>
#include <iterator>
#include <type_traits>
#include <algorithm>


struct partition_largest_axis_at_median_element
{
  template<class BoundingBox>
  static std::array<float,3> centroid(const BoundingBox& box)
  {
    std::array<float,3> result{(box[1][0] + box[0][0])/2,
                               (box[1][1] + box[0][1])/2,
                               (box[1][2] + box[0][2])/2};
    return result;
  }


  template<typename Bounder>
  struct sort_bounding_boxes_by_axis
  {
    sort_bounding_boxes_by_axis(const size_t axis_, Bounder bounder_)
      :axis(axis_),bounder(bounder_)
    {}

    template<class Element>
    bool operator()(const Element& lhs, const Element& rhs) const
    {
      auto lhs_val = centroid(bounder(lhs))[axis];
      auto rhs_val = centroid(bounder(rhs))[axis];

      return lhs_val < rhs_val;
    }

    size_t axis;
    Bounder bounder;
  };


  template<class BoundingBox>
  static size_t largest_axis(const BoundingBox& box)
  {
    // find the largest dimension of the box
    size_t axis = 0;
    float largest_length = -std::numeric_limits<float>::infinity();
    for(size_t i = 0; i < 3; ++i)
    {
      float length = box[1][i] - box[0][i];
      if(length > largest_length)
      {
        largest_length = length;
        axis = i;
      }
    }

    return axis;
  }


  template<class Iterator, class BoundingBox, class Bounder>
  Iterator operator()(Iterator first, Iterator last, const BoundingBox& box, Bounder bounder) const
  {
    // create an ordering
    sort_bounding_boxes_by_axis<Bounder> compare(largest_axis(box), bounder);

    // split at median element
    Iterator split = first + (last - first) / 2;

    std::nth_element(first, split, last, compare);

    return split;
  }
};

