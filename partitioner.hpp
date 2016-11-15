#pragma once

#include <utility>
#include <iterator>
#include <type_traits>
#include <algorithm>


template<class BoundingBox>
BoundingBox empty_box()
{
  float inf = std::numeric_limits<float>::infinity();
  BoundingBox result{{{inf, inf, inf}, {-inf, -inf, -inf}}};
  return result;
}


template<class BoundingBox>
BoundingBox combine_bounding_boxes(const BoundingBox& a, const BoundingBox& b)
{
  BoundingBox result = a;

  for(int i = 0; i < 3; ++i)
  {
    result[0][i] = std::min(result[0][i], b[0][i]);
    result[1][i] = std::max(result[1][i], b[1][i]);
  }

  return result;
}


template<class BoundingBox, class Point>
BoundingBox add_point_to_bounding_box(const BoundingBox& b, const Point& p)
{
  BoundingBox result = b;

  for(int i = 0; i < 3; ++i)
  {
    result[0][i] = std::min(result[0][i], p[i]);
    result[1][i] = std::max(result[1][i], p[i]);
  }

  return result;
}


template<class BoundingBox>
static std::array<float,3> centroid(const BoundingBox& box)
{
  std::array<float,3> result{(box[1][0] + box[0][0])/2,
                             (box[1][1] + box[0][1])/2,
                             (box[1][2] + box[0][2])/2};
  return result;
}


template<class BoundingBox>
size_t largest_axis(const BoundingBox& box)
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


struct partition_largest_axis_at_median_element
{
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

    bounding_box_type result = empty_box<bounding_box_type>();
        
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

    // compute the bounding box of elements' centroids
    BoundingBox centroid_bounding_box = empty_box<BoundingBox>();
    for(Iterator i = first; i != last; ++i)
    {
      centroid_bounding_box = add_point_to_bounding_box(centroid_bounding_box, centroid(bounder(*i)));
    }

    // select an axis to split
    int axis = largest_axis(box);

    float buckets_min = centroid_bounding_box[0][axis];
    float buckets_max = centroid_bounding_box[1][axis];

    constexpr size_t num_buckets = 10;
    float bucket_width = (buckets_max - buckets_min) / num_buckets;

    struct bucket
    {
      float cost;
      float centroid;
      size_t num_elements_in_left_partition;
      BoundingBox left_box;
      BoundingBox right_box;

      bucket()
        : cost{}, centroid{}, num_elements_in_left_partition(0),
          left_box(empty_box<BoundingBox>()),
          right_box(empty_box<BoundingBox>())
      {}

      bool operator<(const bucket& other)
      {
        return cost < other.cost;
      }
    };

    // initialize buckets' centroids
    std::array<bucket, num_buckets> buckets;
    buckets[0].centroid = buckets_min + bucket_width/2;
    for(int i = 1; i < buckets.size(); ++i)
    {
      // each bucket's centroid is at an offset bucket_width from the previous
      buckets[i].centroid = buckets[i-1].centroid + bucket_width;
    }

    if(buckets.front().centroid == buckets.back().centroid)
    {
      // in this degenerate case, we weren't able to subdivide the space
      // as a consequence, we will not be able to partition elements into exactly two sets
      // so use a different partitioning strategy
      return partition_largest_axis_at_median_element()(first, last, box, bounder);
    }

    for(Iterator i = first; i != last; ++i)
    {
      auto this_box = bounder(*i);
      auto this_centroid = centroid(this_box);

      // for each bucket, find which side of centroid element i falls into
      for(bucket& b : buckets)
      {
        if(this_centroid[axis] < b.centroid)
        {
          b.left_box = combine_bounding_boxes(b.left_box, this_box);
          ++b.num_elements_in_left_partition;
        }
        else
        {
          b.right_box = combine_bounding_boxes(b.right_box, this_box);
        }
      }
    }

    // compute the cost of partitioning the elements at each bucket's centroid 
    for(bucket& b : buckets)
    {
      size_t num_elements_in_right_partition = num_elements - b.num_elements_in_left_partition; 

      float left_area = surface_area(b.left_box);
      float right_area = surface_area(b.right_box);

      // compute the surface area heuristic cost of the proposed split
      b.cost = left_area * float(b.num_elements_in_left_partition) + right_area * float(num_elements_in_right_partition);
    }

    // find the bucket which produces non-empty partitions with the smallest traversal cost
    std::sort(buckets.begin(), buckets.end());
    int bucket = -1;
    for(int i = 0; i < buckets.size(); ++i)
    {
      if(buckets[i].num_elements_in_left_partition > 0 && buckets[i].num_elements_in_left_partition != num_elements)
      {
        bucket = i;
        break;
      }
    }

    if(bucket == -1)
    {
      // we weren't able to partition the elements into exactly two subsets, so use a different partitioning strategy
      return partition_largest_axis_at_median_element()(first, last, box, bounder);
    }

    // partition the elements based on whether their centroids are on the left or the right of the selected bucket's centroid
    return std::partition(first, last, [&](const auto& element)
    {
      return centroid(bounder(element))[axis] < buckets[bucket].centroid;
    });
  }
};

