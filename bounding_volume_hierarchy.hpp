#pragma once

#include <vector>
#include <gpcpu/Vector.h>
#include <cassert>

template<typename PrimitiveType, typename PointType, typename RealType = float>
  class bounding_volume_hierarchy
{
  public:
    typedef PrimitiveType Primitive;
    typedef PointType Point;
    typedef RealType Real;

    static const float EPS;

    template<typename Bounder>
      void build(const std::vector<Primitive> &primitives,
                 Bounder &bound);

    template<typename Intersector>
      bool intersect(const Point &o, const Point &d,
                     Real tMin, Real tMax,
                     Intersector &intersect) const;

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

    template<typename Bounder>
    void build(const size_t miss_index,
               const size_t right_brother_index,
               std::vector<size_t>::iterator begin,
               std::vector<size_t>::iterator end,
               const std::vector<PrimitiveType> &primitives,
               Bounder &bound);

    static size_t findPrincipalAxis(const Point &min,
                                          const Point &max);

    struct node
    {
      size_t left_child_index_;
      size_t right_child_index_;
      size_t hit_index_;
      size_t miss_index_;
      Point min_corner_;
      Point max_corner_;

      node(size_t left_child_index,
           size_t right_child_index,
           size_t hit_index,
           size_t miss_index,
           const Point& min_corner,
           const Point& max_corner)
        : left_child_index_(left_child_index),
          right_child_index_(right_child_index),
          hit_index_(hit_index),
          miss_index_(miss_index),
          min_corner_(min_corner),
          max_corner_(max_corner)
      {}

      node(size_t hit_index, size_t miss_index, size_t primitive_index)
        : left_child_index_(null_node),
          right_child_index_(primitive_index),
          hit_index_(hit_index),
          miss_index_(miss_index)
      {}

      bool is_leaf() const
      {
        return left_child_index_ == null_node;
      }

      size_t primitive_index() const
      {
        assert(is_leaf());
        return right_child_index_;
      }
    };

    std::vector<node> nodes_;
    static const size_t null_node;
};

#include "bounding_volume_hierarchy.inl"

