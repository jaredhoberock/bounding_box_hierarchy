#pragma once

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <utility>
#include <functional>
#include <type_traits>
#include <cassert>

#include "optional.hpp"


template<class T, class Bounder>
struct memoized_bounder
{
  // we're memoizing the result of Bounder, which returns a bounding_box_type
  // that may be different than the type of bounding box used by bounding_volume_hierarchy
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


template<typename T>
class bounding_volume_hierarchy
{
  private:
    struct call_member_intersect
    {
      template<class... Args>
      auto operator()(const T& element, Args&&... args) const
      {
        return element.intersect(std::forward<Args>(args)...);
      }
    };


    struct call_member_bounding_box
    {
      auto operator()(const T& element) const
      {
        return element.bounding_box();
      }
    };


    struct convert_to_float
    {
      template<class U>
      float operator()(const U& x) const
      {
        return static_cast<float>(x);
      }
    };


  public:
    using element_type = T;


    // a bounding box is an array of two arrays of three floats
    // XXX if T::bounding_box() exists, we should use its result for our bounding_box_type
    using bounding_box_type = std::array<std::array<float,3>,2>;


    template<class ContiguousRange,
             class Bounder = call_member_bounding_box>
    bounding_volume_hierarchy(const ContiguousRange& elements,
                              Bounder bounder = call_member_bounding_box(),
                              float epsilon = std::numeric_limits<float>::epsilon())
      : elements_(&*elements.begin()), nodes_(make_tree(elements, bounder, epsilon))
    {}


    bounding_box_type bounding_box() const
    {
      return bounding_box(root_node());
    }


    /// Finds the nearest intersection, if any, between a ray and the elements of this bounding_volume_hierarchy.
    /// \param origin The ray's origin.
    /// \param direction The ray's direction.
    /// \param interval The parametric ray interval to consider. Defaults to [0,1).
    /// \param intersector A function to test for intersection between a ray and element.
    ///                    Signature is expected to be result_type intersector(T,Point,Direction,Interval). Defaults to `T::intersect()`.
    /// \param hit_time A function to return the value of the ray parameter at the time of an intersection.
    ///                 Signature is expected to be `float hit_time(result_type)`, where `result_type` is the type returned by `intersector`.
    ///                 By default, this function converts `result_type` to `float`.
    /// \return The nearest result of `intersector` along the ray, if an intersection exists. Empty, otherwise.
    //
    // XXX maybe the default hit_time should be to call result_type::hit_time() if result_type is not convertible to float?
    // XXX another default could be to call std::get<float>() in case the intersection type is tuple-like
    template<class Point, class Vector,
             class Interval = std::array<float,2>,
             class Function1 = call_member_intersect,
             class Function2 = convert_to_float>
    std::result_of_t<Function1(T,Point,Vector,Interval)>
      intersect(Point origin, Vector direction,
                Interval interval = Interval{0.f, 1.f},
                Function1 intersector = call_member_intersect(),
                Function2 hit_time = convert_to_float()) const
    {
      Vector one_over_direction = {1.f/direction[0], 1.f/direction[1], 1.f/direction[2]};

      std::result_of_t<Function1(T,Point,Vector,Interval)> result;

      const node* current_node = root_node();

      while(current_node != nullptr)
      {
        bool hit_current_node = false;

        if(!is_leaf(current_node))
        {
          hit_current_node = intersect_box(bounding_box(current_node), origin, one_over_direction, interval);
        }
        else
        {
          auto current_result = intersector(element(current_node), origin, direction, interval);
          
          if(current_result)
          {
            // shorten interval
            interval[1] = hit_time(*current_result);

            // update result
            result = *current_result;
          }
        }

        current_node = hit_current_node ? current_node->hit_node_ : current_node->miss_node_;
      }

      return result;
    }


  private:
    template<class Point, class Vector, class Interval>
    static bool intersect_box(const bounding_box_type& box,
                              Point origin,
                              Vector one_over_direction,
                              Interval interval)
    {
      Point t_min3, t_max3;
      for(int i = 0; i < 3; ++i)
      {
        t_min3[i] = (box[0][i] - origin[i]) * one_over_direction[i];
        t_max3[i] = (box[1][i] - origin[i]) * one_over_direction[i];
      }

      Point t_near3{std::min(t_min3[0], t_max3[0]), std::min(t_min3[1], t_max3[1]), std::min(t_min3[2], t_max3[2])};
      Point t_far3{ std::max(t_min3[0], t_max3[0]), std::max(t_min3[1], t_max3[1]), std::max(t_min3[2], t_max3[2])};

      auto t_near = std::max(std::max(t_near3[0], t_near3[1]), t_near3[2]);
      auto t_far  = std::min(std::min(t_far3[0],  t_far3[1]),  t_far3[2]);

      bool hit = t_near <= t_far;
      return hit && interval[0] <= t_far && t_near <= interval[1];
    }


    template<class ContiguousRange, class Bounder>
    static bounding_box_type bounding_box(const std::vector<size_t>::iterator begin,
                                          const std::vector<size_t>::iterator end,
                                          const ContiguousRange& elements,
                                          Bounder bounder,
                                          float epsilon)
    {
      float inf = std::numeric_limits<float>::infinity();
      bounding_box_type result{{{inf, inf, inf}, {-inf, -inf, -inf}}};
          
      for(std::vector<size_t>::iterator t = begin; t != end; ++t)
      {
        auto bounding_box = bounder(elements[*t]);

        for(int i = 0; i < 3; ++i)
        {
          result[0][i] = std::min(result[0][i], bounding_box[0][i]);
          result[1][i] = std::max(result[1][i], bounding_box[1][i]);
        }
      }

      // widen the bounding box by 2*epsilon
      // XXX might want to instead ensure that each side is at least epsilon in width
      // this ensures that axis-aligned elements always
      // lie strictly within the bounding box
      for(size_t i = 0; i != 3; ++i)
      {
        result[0][i] -= epsilon;
        result[1][i] += epsilon;
      }

      return result;
    }


    template<typename Bounder>
    struct sort_bounding_boxes_by_axis
    {
      template<class ContiguousRange>
      sort_bounding_boxes_by_axis(const size_t axis_,
                                  const ContiguousRange& elements_,
                                  Bounder bounder_)
        :axis(axis_),elements(&*elements_.begin()),bounder(bounder_)
      {}

      bool operator()(const size_t lhs, const size_t rhs) const
      {
        auto lhs_val = bounder(elements[lhs])[0][axis];
        auto rhs_val = bounder(elements[rhs])[0][axis];

        return lhs_val < rhs_val;
      }

      size_t axis;
      const T* elements;
      Bounder bounder;
    };


    static size_t largest_axis(const bounding_box_type& box)
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


    struct node
    {
      const node* hit_node_;
      const node* miss_node_;
      bounding_box_type bounding_box_;

      node() = default;

      node(const node* hit_node,
           const node* miss_node,
           const bounding_box_type& bounding_box)
        : hit_node_(hit_node),
          miss_node_(miss_node),
          bounding_box_(bounding_box)
      {}

      node(const node* hit_node, const node* miss_node, size_t primitive_index)
        : hit_node_(hit_node),
          miss_node_(miss_node)
      {}
    };


    template<class ContiguousRange, class Bounder>
    static const node* make_tree_recursive(std::vector<node>& tree,
                                           const node* miss_node,
                                           const node* right_brother,
                                           std::vector<size_t>::iterator begin,
                                           std::vector<size_t>::iterator end,
                                           const ContiguousRange &elements,
                                           Bounder bounder,
                                           float epsilon)
    {
      if(begin + 1 == end)
      {
        // a right leaf's hit index is always the miss index
        // a left leaf's hit index is always the right brother
        const node* hit_node = right_brother == nullptr ? miss_node : right_brother;

        // leaves come at the beginning of the tree
        node* result = &tree[*begin];
        *result = node(hit_node, miss_node, *begin);
        return result;
      }
      else
      {
        // find the bounding box of the elements
        bounding_box_type box = bounding_box(begin, end, elements, bounder, epsilon);

        // create an ordering
        sort_bounding_boxes_by_axis<Bounder> compare(largest_axis(box),elements,bounder);
        
        // sort the median
        std::vector<size_t>::iterator split = begin + (end - begin) / 2;

        std::nth_element(begin, split, end, compare);

        // build the right subtree first
        const node* right_child = make_tree_recursive(tree, miss_node, nullptr, split, end, elements, bounder, epsilon);

        // we have to pass the right child to the left subtree build
        const node* left_child = make_tree_recursive(tree, right_child, right_child, begin, split, elements, bounder, epsilon);

        // create a new node
        tree.emplace_back(left_child, miss_node, box);
        return &tree.back();
      }
    }


    template<class ContiguousRange, class Bounder>
    static std::vector<node> make_tree(const ContiguousRange& elements, Bounder bounder, float epsilon)
    {
      // we will sort an array of indices
      std::vector<size_t> indices(elements.size());
      std::iota(indices.begin(), indices.end(), 0);
    
      // reserve 2*n - 1 nodes to ensure that no iterators are invalidated during construction
      std::vector<node> tree;
      tree.reserve(2 * elements.size() - 1);

      // create the leaves first so that they come first in the array
      tree.resize(elements.size());
    
      // memoize the bound function
      memoized_bounder<T,Bounder> memoized_bounder(elements,bounder);
    
      // recurse
      make_tree_recursive(tree,
                          nullptr,
                          nullptr,
                          indices.begin(),
                          indices.end(),
                          elements,
                          std::ref(memoized_bounder),
                          epsilon);
    
      assert(tree.size() == 2 * elements.size() - 1);

      return tree;
    }

    
    const T& element(const node* leaf) const
    {
      // the index of the element corresponding to the leaf node is the same as the leaf node's index
      size_t element_idx = leaf - nodes_.data();
      return elements_[element_idx];
    }

    const bounding_box_type& bounding_box(const node* n) const
    {
      return n->bounding_box_;
    }


    const node* leaves_end() const
    {
      return nodes_.data() + (nodes_.size() + 1) / 2;
    }


    bool is_leaf(const node* n) const
    {
      return n < leaves_end();
    }


    const node* root_node() const
    {
      return &nodes_.back();
    }


    const T* elements_;
    std::vector<node> nodes_;
};

