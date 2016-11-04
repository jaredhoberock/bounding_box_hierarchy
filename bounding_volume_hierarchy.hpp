#pragma once

#include <vector>
#include <gpcpu/Vector.h>
#include <algorithm>
#include <numeric>
#include <cassert>

template<typename PrimitiveType, typename PointType, typename RealType = float>
  class bounding_volume_hierarchy
{
  public:
    typedef PrimitiveType Primitive;
    typedef PointType Point;
    typedef RealType Real;

    static const float EPS;

    template<class Bounder>
    bounding_volume_hierarchy(const std::vector<Primitive>& primitives, Bounder& bound)
      : nodes_(make_tree(primitives, bound))
    {}

    template<typename Intersector>
    bool intersect(const Point &o, const Point &d, Real tMin, Real tMax, Intersector &intersect) const
    {
      Point invDir;
      invDir[0] = Real(1.0) / d[0];
      invDir[1] = Real(1.0) / d[1];
      invDir[2] = Real(1.0) / d[2];

      const node* current_node = root_node();
      bool hit = false;
      bool result = false;
      float t = tMax;
      while(current_node != nullptr)
      {
        if(!is_leaf(current_node))
        {
          hit = intersectBox(o, invDir,
                             current_node->min_corner_,
                             current_node->max_corner_,
                             tMin, tMax);
        }
        else
        {
          hit = intersect(o,d,current_node->primitive_index(),t) && t < tMax && t > tMin;
          result |= hit;
          if(hit)
            tMax = std::min(t, tMax);
        }

        current_node = hit ? current_node->hit_node_ : current_node->miss_node_;
      }

      return result;
    }

    static bool intersectBox(const Point &o, const Point &invDir,
                             const Point &minBounds, const Point &maxBounds, 
                             const Real &tMin, const Real &tMax);

  protected:
    /*! The idea of this class is to wrap Bounder
     *  and accelerate build() by caching the results
     *  of Bounder.
     *
     *  This gives about a 10x build speedup on the bunny
     *  in Cornell box scene.
     */
    template<typename Bounder>
      class CachedBounder
    {
      public:
        inline CachedBounder(Bounder &bound,
                             const std::vector<Primitive> &primitives);

        inline float operator()(const size_t axis,
                                const bool min,
                                const size_t primIndex)
        {
          if(min)
          {
            return mPrimMinBounds[axis][primIndex];
          }

          return mPrimMaxBounds[axis][primIndex];
        } // end operator()()

      protected:
        std::vector<RealType> mPrimMinBounds[3];
        std::vector<RealType> mPrimMaxBounds[3];
    }; // end CachedBounder

    template<typename Bounder>
      static void findBounds(const std::vector<size_t>::iterator begin,
                             const std::vector<size_t>::iterator end,
                             const std::vector<Primitive> &primitives,
                             CachedBounder<Bounder> &bound,
                             Point &m, Point &M);

    template<typename Bounder>
      struct PrimitiveSorter
    {
      inline PrimitiveSorter(const size_t axis,
                             const std::vector<PrimitiveType> &primitives,
                             Bounder &bound)
        :mAxis(axis),mPrimitives(primitives),mBound(bound)
      {
        ;
      } // end PrimitiveSorter()

      inline bool operator()(const size_t lhs, const size_t rhs) const
      {
        return mBound(mAxis, true, lhs) < mBound(mAxis, true, rhs);
      } // end operator<()

      size_t mAxis;
      const std::vector<PrimitiveType> &mPrimitives;
      Bounder &mBound;
    }; // end PrimitiveSorter

    static size_t findPrincipalAxis(const Point &min, const Point &max);

    struct node
    {
      const node* hit_node_;
      const node* miss_node_;
      Point min_corner_;
      Point max_corner_;

      node() = default;

      node(const node* hit_node,
           const node* miss_node,
           const Point& min_corner,
           const Point& max_corner)
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
          min_corner_(coerce<float>(primitive_index), 0, 0)
      {}

      size_t primitive_index() const
      {
        return coerce<size_t>(min_corner_[0]);
      }
    };


    template<typename Bounder>
    static const node* make_tree_recursive(std::vector<node>& tree,
                                           const node* miss_node,
                                           const node* right_brother,
                                           std::vector<size_t>::iterator begin,
                                           std::vector<size_t>::iterator end,
                                           const std::vector<PrimitiveType> &primitives,
                                           Bounder &bound)
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
        // find the bounds of the primitives
        Point m, M;
        findBounds(begin, end, primitives, bound, m, M);

        size_t axis = findPrincipalAxis(m, M);

        // create an ordering
        PrimitiveSorter<Bounder> sorter(axis,primitives,bound);
        
        // sort the median
        std::vector<size_t>::iterator split = begin + (end - begin) / 2;

        std::nth_element(begin, split, end, sorter);

        // build the right subtree first
        const node* right_child = make_tree_recursive(tree, miss_node, nullptr, split, end, primitives, bound);

        // we have to pass the right child to the left subtree build
        const node* left_child = make_tree_recursive(tree, right_child, right_child, begin, split, primitives, bound);

        // create a new node
        tree.emplace_back(left_child, miss_node, m, M);
        return &tree.back();
      }
    }


    template<typename Bounder>
    static std::vector<node> make_tree(const std::vector<Primitive> &primitives, Bounder &bound)
    {
      // we will sort an array of indices
      std::vector<size_t> indices(primitives.size());
      std::iota(indices.begin(), indices.end(), 0);
    
      // reserve 2*n - 1 nodes to ensure that no iterators are invalidated during construction
      std::vector<node> tree;
      tree.reserve(2 * primitives.size() - 1);

      // create the leaves first so that they come first in the array
      tree.resize(primitives.size());
    
      // create a CachedBounder
      CachedBounder<Bounder> cachedBound(bound,primitives);
    
      // recurse
      make_tree_recursive(tree,
                          nullptr,
                          nullptr,
                          indices.begin(),
                          indices.end(),
                          primitives,
                          cachedBound);
    
      assert(tree.size() == 2 * primitives.size() - 1);

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

    std::vector<node> nodes_;
};

#include "bounding_volume_hierarchy.inl"

