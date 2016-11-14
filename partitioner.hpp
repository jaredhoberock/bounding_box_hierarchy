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
  Iterator operator()(Iterator first, Iterator middle, Iterator last, const BoundingBox& box, Bounder bounder) const
  {
    // create an ordering
    sort_bounding_boxes_by_axis<Bounder> compare(largest_axis(box), bounder);

    std::nth_element(first, middle, last, compare);

    return middle;
  }


  template<class Iterator, class BoundingBox, class Bounder>
  Iterator operator()(Iterator first, Iterator last, const BoundingBox& box, Bounder bounder) const
  {
    // split at median element
    Iterator split = first + (last - first) / 2;

    return operator()(first, split, last, box, bounder);
  }
};


struct minimize_surface_area_heuristic
{
  template<class Iterator, class Bounder>
  static auto bounding_box(Iterator first, Iterator last, Bounder bounder)
  {
    using bounding_box_type = std::result_of_t<Bounder(typename std::iterator_traits<Iterator>::reference)>;

    float inf = std::numeric_limits<float>::infinity();
    bounding_box_type result{{{inf, inf, inf}, {-inf, -inf, -inf}}};
        
    for(; first != last; ++first)
    {
      auto bounding_box = bounder(*first);

      for(int i = 0; i < 3; ++i)
      {
        result[0][i] = std::min(result[0][i], bounding_box[0][i]);
        result[1][i] = std::max(result[1][i], bounding_box[1][i]);
      }
    }

    return result;
  }


  template<class BoundingBox>
  static float surface_area(const BoundingBox& box)
  {
    float width  = box[1][0] - box[0][0];
    float height = box[1][1] - box[0][1];
    float depth  = box[1][2] - box[0][2];

    return 2.f * (width * height + width * depth + height * depth);
  }


  template<class Iterator, class BoundingBox, class Bounder>
  Iterator operator()(Iterator first, Iterator last, const BoundingBox& box, Bounder bounder) const
  {
    size_t num_elements = last - first;

    size_t best_candidate = 0;
    float smallest_cost = 0;

    partition_largest_axis_at_median_element partitioner;

    std::vector<typename std::iterator_traits<Iterator>::value_type> elements;

    size_t stride = std::max(1ul, num_elements / 1000);
    for(size_t i = 0; i < (last - first); i += stride)
    {
      elements.assign(first, last);

      auto candidate = partitioner(elements.begin(), elements.begin() + i, elements.end(), box, bounder);

      auto left = bounding_box(elements.begin(), candidate, bounder);
      auto right = bounding_box(candidate, elements.end(), bounder);

      float left_surface_area = surface_area(left);
      float right_surface_area = surface_area(right);

      float left_cost = float(candidate - elements.begin()) * (left_surface_area / right_surface_area);
      float right_cost = float(elements.end() - candidate) * (right_surface_area / left_surface_area);

      if(left_cost + right_cost < smallest_cost)
      {
        smallest_cost = left_cost + right_cost;
        best_candidate = candidate - elements.begin();
      }
    }

    return partitioner(first, first + best_candidate, last, box, bounder);
  }
};

