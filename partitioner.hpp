#pragma once

#include <array>
#include <algorithm>
#include <limits>


struct partition_largest_axis_at_middle_element
{
  template<class BoundingBox>
  static std::array<float,3> centroid(const BoundingBox& box)
  {
    std::array<float,3> result{(box[1][0] + box[0][0])/2,
                               (box[1][1] + box[0][1])/2,
                               (box[1][2] + box[0][2])/2};
    return result;
  }


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
  Iterator operator()(Iterator first, Iterator last, const BoundingBox& box, Bounder bounder) const
  {
    // split at middle element
    Iterator middle = first + (last - first) / 2;

    // create an ordering
    sort_bounding_boxes_by_axis<Bounder> compare(largest_axis(box), bounder);

    std::nth_element(first, middle, last, compare);

    return middle;
  }
};


struct minimize_surface_area_heuristic
{
  template<class BoundingBox>
  static std::array<float,3> centroid(const BoundingBox& box)
  {
    std::array<float,3> result{(box[1][0] + box[0][0])/2,
                               (box[1][1] + box[0][1])/2,
                               (box[1][2] + box[0][2])/2};
    return result;
  }


  template<class BoundingBox>
  static BoundingBox empty_box()
  {
    float inf = std::numeric_limits<float>::infinity();
    BoundingBox result{{{inf, inf, inf}, {-inf, -inf, -inf}}};
    return result;
  }


  template<class BoundingBox>
  static BoundingBox combine_bounding_boxes(const BoundingBox& a, const BoundingBox& b)
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
  static BoundingBox add_point_to_bounding_box(const BoundingBox& b, const Point& p)
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
    // compute the bounding box of elements' centroids
    BoundingBox centroid_bounding_box = empty_box<BoundingBox>();
    for(Iterator i = first; i != last; ++i)
    {
      centroid_bounding_box = add_point_to_bounding_box(centroid_bounding_box, centroid(bounder(*i)));
    }

    struct bucket
    {
      float cost;
      int axis;
      float centroid;
      size_t num_elements_in_left_partition;
      BoundingBox left_box;
      BoundingBox right_box;

      bucket()
        : num_elements_in_left_partition(0),
          left_box(empty_box<BoundingBox>()),
          right_box(empty_box<BoundingBox>())
      {}

      bool operator<(const bucket& other) const
      {
        return cost < other.cost;
      }
    };

    // initialize buckets' axes and centroids
    constexpr size_t num_buckets_per_axis = 10;
    std::array<bucket, 3 * num_buckets_per_axis> buckets;

    for(int axis = 0; axis < 3; ++axis)
    {
      int axis_begin = axis * num_buckets_per_axis;
      int axis_end = axis_begin + num_buckets_per_axis;

      float buckets_min = centroid_bounding_box[0][axis];
      float buckets_max = centroid_bounding_box[1][axis];
      float bucket_width = (buckets_max - buckets_min) / num_buckets_per_axis;

      buckets[axis_begin].axis = axis;
      buckets[axis_begin].centroid = buckets_min + bucket_width/2;
      for(int i = axis_begin + 1; i < axis_end; ++i)
      {
        // each bucket's centroid is at an offset bucket_width from the previous
        buckets[i].axis = axis;
        buckets[i].centroid = buckets[i-1].centroid + bucket_width;
      }
    }

    for(Iterator i = first; i != last; ++i)
    {
      auto this_box = bounder(*i);
      auto this_centroid = centroid(this_box);

      // for each bucket, find which side of its centroid element i falls into
      for(bucket& b : buckets)
      {
        if(this_centroid[b.axis] < b.centroid)
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

    size_t num_elements = last - first;

    // compute the cost of partitioning the elements at each bucket's centroid 
    for(bucket& b : buckets)
    {
      size_t num_elements_in_right_partition = num_elements - b.num_elements_in_left_partition; 

      float left_area = surface_area(b.left_box);
      float right_area = surface_area(b.right_box);

      // compute the surface area heuristic cost of the proposed split
      b.cost = left_area * float(b.num_elements_in_left_partition) + right_area * float(num_elements_in_right_partition);
    }

    // remove buckets which have a NaN cost or which produce partitions with empty sets
    auto valid_buckets_end = std::remove_if(buckets.begin(), buckets.end(), [=](const bucket& b)
    {
      return (b.cost != b.cost) || (b.num_elements_in_left_partition == 0) || (b.num_elements_in_left_partition == num_elements);
    });

    if(buckets.begin() == valid_buckets_end)
    {
      // we weren't able to partition the elements into exactly two subsets, so use a different partitioning strategy
      return partition_largest_axis_at_middle_element()(first, last, box, bounder);
    }

    // select the bucket with minimal splitting cost
    auto selected_bucket = std::min_element(buckets.begin(), valid_buckets_end);

    // partition the elements based on whether their centroids are on the left or the right of the selected bucket's centroid
    return std::partition(first, last, [&](const auto& element)
    {
      int axis = selected_bucket->axis;
      return centroid(bounder(element))[axis] < selected_bucket->centroid;
    });
  }
};

