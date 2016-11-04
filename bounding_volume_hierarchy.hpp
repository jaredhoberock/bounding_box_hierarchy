#pragma once

#include <vector>
#include <gpcpu/Vector.h>

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
    /*! This method adds a new Node to this hierarchy.
     *  \param parent The parent of the Node to add.
     *  \return The index of the node.
     */
    size_t addNode(const size_t parent);

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
    size_t build(const size_t parent,
                 std::vector<size_t>::iterator begin,
                 std::vector<size_t>::iterator end,
                 const std::vector<PrimitiveType> &primitives,
                 Bounder &bound);

    static size_t findPrincipalAxis(const Point &min,
                                          const Point &max);

    /*! This method computes the index of the next node in a
     *  depth first traversal of this tree, from node i.
     *  \param i The Node of interest.
     *  \return The index of the next Node from i, if it exists;
     *          NULL_NODE, otherwise.
     */
    size_t computeHitIndex(const size_t i) const;

    /*! This method computes the index of the next node in a
     *  depth first traversal of this tree, from node i.
     *  \param i The Node of interest.
     *  \return The index of the next Node from i, if it exists;
     *          NULL_NODE, otherwise.
     */
    size_t computeMissIndex(const size_t i) const;

    /*! This method computes the index of a Node's brother to the right,
     *  if it exists.
     *  \param i The index of the Node of interest.
     *  \return The index of Node i's brother to the right, if it exists;
     *          NULL_NODE, otherwise.
     */
    size_t computeRightBrotherIndex(const size_t i) const;

    struct node
    {
      size_t parent_index_;
      size_t left_child_index_;
      size_t right_child_index_;
      Point min_corner_;
      Point max_corner_;
      size_t hit_index_;
      size_t miss_index_;

      node() {}

      node(size_t parent)
        : parent_index_(parent),
          left_child_index_(null_node),
          right_child_index_(null_node)
      {}

      bool is_leaf() const
      {
        return left_child_index_ == null_node;
      }
    };

    std::vector<node> nodes_;
    static const size_t null_node;
    size_t root_index_;
};

#include "bounding_volume_hierarchy.inl"

