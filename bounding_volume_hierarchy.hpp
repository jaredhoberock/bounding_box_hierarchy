#pragma once

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <utility>
#include <tuple>
#include <type_traits>
#include <cassert>

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

    using point = std::array<float,3>;

  public:
    using element_type = T;

    // XXX Bounder should take const T& and return pair<Point,Point>
    // XXX Bounder should default to call_member_bounding_box
    //
    // XXX bounder should come after epsilon? to match the parameter order of similar parameters of intersect()?
    template<class ContiguousRange, class Bounder>
    bounding_volume_hierarchy(const ContiguousRange& elements, Bounder bounder, float epsilon = std::numeric_limits<float>::epsilon())
      : elements_(&*elements.begin()), nodes_(make_tree(elements, bounder, epsilon))
    {}

    auto bounding_box() const
    {
      return std::make_pair(root_node().min_corner_, root_node().max_corner_);
    }

    template<class Point, class Vector, typename Interval = std::array<float,2>, class Function = call_member_intersect>
    bool intersect(Point origin, Vector direction, Interval interval = Interval{0.f, 1.f}, Function intersector = call_member_intersect()) const
    {
      point one_over_direction = {1.f/direction[0], 1.f/direction[1], 1.f/direction[2]};

      const node* current_node = root_node();
      bool hit = false;
      bool result = false;
      auto t = interval[1];
      while(current_node != nullptr)
      {
        if(!is_leaf(current_node))
        {
          hit = intersect_box(origin, one_over_direction, current_node->min_corner_, current_node->max_corner_, interval);
        }
        else
        {
          // the index of the element corresponding to the leaf node is the same as the leaf node's index
          size_t element_idx = current_node - nodes_.data();

          hit = intersector(elements_[element_idx], origin, direction, t) && interval[0] < t && t < interval[1];
          result |= hit;
          if(hit)
            interval[1] = std::min(t, interval[1]);
        }

        current_node = hit ? current_node->hit_node_ : current_node->miss_node_;
      }

      return result;
    }

  private:
    template<class Point, class Vector, class Interval>
    static bool intersect_box(Point origin, Vector one_over_direction, const point &min_corner, const point &max_corner, Interval interval)
    {
      point t_min3, t_max3;
      for(int i = 0; i < 3; ++i)
      {
        t_min3[i] = (min_corner[i] - origin[i]) * one_over_direction[i];
        t_max3[i] = (max_corner[i] - origin[i]) * one_over_direction[i];
      }

      point t_near3{std::min(t_min3[0], t_max3[0]), std::min(t_min3[1], t_max3[1]), std::min(t_min3[2], t_max3[2])};
      point  t_far3{std::max(t_min3[0], t_max3[0]), std::max(t_min3[1], t_max3[1]), std::max(t_min3[2], t_max3[2])};

      auto t_near = std::max(std::max(t_near3[0], t_near3[1]), t_near3[2]);
      auto t_far  = std::min(std::min(t_far3[0],  t_far3[1]),  t_far3[2]);

      bool hit = t_near <= t_far;
      return hit && interval[0] <= t_far && t_near <= interval[1];
    }


    // XXX this should probably be external to this class
    template<class Bounder>
    struct memoized_bounder
    {
      using pair_of_points = std::result_of_t<Bounder(const T&)>;
      using min_corner_type = std::tuple_element_t<0,pair_of_points>;
      using max_corner_type = std::tuple_element_t<1,pair_of_points>;

      // avoid copying this thing unintentionally
      memoized_bounder(memoized_bounder&&) = delete;

      template<class ContiguousRange>
      memoized_bounder(const ContiguousRange& elements, Bounder bounder)
        : elements_(&*elements.begin())
      {
        min_corners.resize(elements.size());
        max_corners.resize(elements.size());

        size_t i = 0;
        for(auto element = elements.begin(); element != elements.end(); ++element, ++i)
        {
          // XXX need to just call bounder once and destructure the result with tie()
          min_corners[i][0] = bounder(*element, 0, true);
          min_corners[i][1] = bounder(*element, 1, true);
          min_corners[i][2] = bounder(*element, 2, true);

          max_corners[i][0] = bounder(*element, 0, false);
          max_corners[i][1] = bounder(*element, 1, false);
          max_corners[i][2] = bounder(*element, 2, false);
        }
      }

      // XXX operator()() needs to be:
      // result_of_t<Bounder(T)> operator()(const T& element) const

      float operator()(const size_t axis, const bool min, size_t element_idx)
      {
        if(min)
        {
          return min_corners[element_idx][axis];
        }

        return max_corners[element_idx][axis];
      }

      pair_of_points operator()(const T& element) const
      {
        size_t element_idx = &element - elements_;
        return pair_of_points(min_corners[element_idx], max_corners[element_idx]);
      }

      const T* elements_;
      std::vector<min_corner_type> min_corners;
      std::vector<max_corner_type> max_corners;
    };


    template<typename Bounder>
    static std::pair<point,point> bounding_box(const std::vector<size_t>::iterator begin,
                                               const std::vector<size_t>::iterator end,
                                               const std::vector<T> &elements,
                                               memoized_bounder<Bounder>& bounder,
                                               float epsilon)
    {
      float inf = std::numeric_limits<float>::infinity();
      point min_corner{ inf,  inf,  inf};
      point max_corner{-inf, -inf, -inf};

      float x;
          
      for(std::vector<size_t>::iterator t = begin;
          t != end;
          ++t)
      {
        auto bounding_box = bounder(elements[*t]);

        for(int i = 0; i < 3; ++i)
        {
          min_corner[i] = std::min(min_corner[i], std::get<0>(bounding_box)[i]);
          max_corner[i] = std::max(max_corner[i], std::get<1>(bounding_box)[i]);
        }
      }

      // widen the bounding box by 2*epsilon
      // XXX might want to instead ensure that each side is at least epsilon in width
      // this ensures that axis-aligned elements always
      // lie strictly within the bounding box
      for(size_t i = 0; i != 3; ++i)
      {
        min_corner[i] -= epsilon;
        max_corner[i] += epsilon;
      }

      return std::make_pair(min_corner, max_corner);
    }

    template<typename Bounder>
      struct element_sorter
    {
      template<class ContiguousRange>
      element_sorter(const size_t axis_,
                     const ContiguousRange& elements_,
                     Bounder &bounder_)
        :axis(axis_),elements(&*elements_.begin()),bounder(bounder_)
      {
        ;
      }

      bool operator()(const size_t lhs, const size_t rhs) const
      {
        return bounder(axis, true, lhs) < bounder(axis, true, rhs);
      }

      size_t axis;
      const T* elements;
      Bounder& bounder; // XXX make this a value
    };


    static size_t findPrincipalAxis(const point &min, const point &max);


    struct node
    {
      const node* hit_node_;
      const node* miss_node_;
      point min_corner_;
      point max_corner_;

      node() = default;

      node(const node* hit_node,
           const node* miss_node,
           const point& min_corner,
           const point& max_corner)
        : hit_node_(hit_node),
          miss_node_(miss_node),
          min_corner_(min_corner),
          max_corner_(max_corner)
      {}

      template<class To, class From>
      static To coerce(const From& from)
      {
        union {From from; To to;} u;
        u.from = from;
        return u.to;
      }

      node(const node* hit_node, const node* miss_node, size_t primitive_index)
        : hit_node_(hit_node),
          miss_node_(miss_node),
          min_corner_{coerce<float>(primitive_index), 0, 0}
      {}

      size_t element_index() const
      {
        return coerce<size_t>(min_corner_[0]);
      }
    };


    template<class ContiguousRange, class Bounder>
    static const node* make_tree_recursive(std::vector<node>& tree,
                                           const node* miss_node,
                                           const node* right_brother,
                                           std::vector<size_t>::iterator begin,
                                           std::vector<size_t>::iterator end,
                                           const ContiguousRange &elements,
                                           Bounder& bounder,
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
        // find the bounds of the elements
        std::pair<point,point> bounds = bounding_box(begin, end, elements, bounder, epsilon);

        size_t axis = findPrincipalAxis(bounds.first, bounds.second);

        // create an ordering
        element_sorter<Bounder> sorter(axis,elements,bounder);
        
        // sort the median
        std::vector<size_t>::iterator split = begin + (end - begin) / 2;

        std::nth_element(begin, split, end, sorter);

        // build the right subtree first
        const node* right_child = make_tree_recursive(tree, miss_node, nullptr, split, end, elements, bounder, epsilon);

        // we have to pass the right child to the left subtree build
        const node* left_child = make_tree_recursive(tree, right_child, right_child, begin, split, elements, bounder, epsilon);

        // create a new node
        tree.emplace_back(left_child, miss_node, bounds.first, bounds.second);
        return &tree.back();
      }
    }


    template<class ContiguousRange, class Bounder>
    static std::vector<node> make_tree(const ContiguousRange &elements, Bounder bounder, float epsilon)
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
      memoized_bounder<Bounder> memoized_bounder(elements,bounder);
    
      // recurse
      make_tree_recursive(tree,
                          nullptr,
                          nullptr,
                          indices.begin(),
                          indices.end(),
                          elements,
                          memoized_bounder,
                          epsilon);
    
      assert(tree.size() == 2 * elements.size() - 1);

      return tree;
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

#include "bounding_volume_hierarchy.inl"

