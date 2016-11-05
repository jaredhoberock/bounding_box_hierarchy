#include "bounding_volume_hierarchy.hpp"
#include <limits>

template<typename PrimitiveType, typename RealType>
const float bounding_volume_hierarchy<PrimitiveType,RealType>::EPS = 0.00005f;


template<typename PrimitiveType,
         typename RealType>
  template<typename Bounder>
    void bounding_volume_hierarchy<PrimitiveType, RealType>
      ::findBounds(const std::vector<size_t>::iterator begin,
                   const std::vector<size_t>::iterator end,
                   const std::vector<Primitive> &primitives,
                   CachedBounder<Bounder> &bound,
                   point &min_corner, point &max_corner)
{
  Real inf = std::numeric_limits<Real>::infinity();
  min_corner = point{ inf,  inf,  inf};
  max_corner = point{-inf, -inf, -inf};

  Real x;
      
  for(std::vector<size_t>::iterator t = begin;
      t != end;
      ++t)
  {
    for(size_t i =0;
        i < 3;
        ++i)
    {
      x = bound(i, true, *t);

      if(x < min_corner[i])
      {
        min_corner[i] = x;
      } // end if

      x = bound(i, false, *t);

      if(x > max_corner[i])
      {
        max_corner[i] = x;
      } // end if
    } // end for j
  } // end for t

  // always widen the bounding box
  // this ensures that axis-aligned primitives always
  // lie strictly within the bounding box
  for(size_t i = 0; i != 3; ++i)
  {
    min_corner[i] -= EPS;
    max_corner[i] += EPS;
  } // end for i
} // end bounding_volume_hierarchy::findBounds()


template<typename PrimitiveType,
         typename RealType>
  size_t bounding_volume_hierarchy<PrimitiveType, RealType>
    ::findPrincipalAxis(const point &min_corner,
                        const point &max_corner)
{
  // find the principal axis of the points
  size_t axis = 4;
  float maxLength = -std::numeric_limits<Real>::infinity();
  float temp;
  for(size_t i = 0; i < 3; ++i)
  {
    temp = max_corner[i] - min_corner[i];
    if(temp > maxLength)
    {
      maxLength = temp;
      axis = i;
    } // end if
  } // end for

  return axis;
} // end bounding_volume_hierarchy::findPrincipalAxis()


template<typename PrimitiveType,
         typename RealType>
  template<typename Bounder>
    bounding_volume_hierarchy<PrimitiveType,RealType>::CachedBounder<Bounder>
      ::CachedBounder(Bounder &bound,
                      const std::vector<Primitive> &primitives)
{
  mPrimMinBounds[0].resize(primitives.size());
  mPrimMinBounds[1].resize(primitives.size());
  mPrimMinBounds[2].resize(primitives.size());

  mPrimMaxBounds[0].resize(primitives.size());
  mPrimMaxBounds[1].resize(primitives.size());
  mPrimMaxBounds[2].resize(primitives.size());

  size_t i = 0;
  for(typename std::vector<PrimitiveType>::const_iterator prim = primitives.begin();
      prim != primitives.end();
      ++prim, ++i)
  {
    mPrimMinBounds[0][i] = bound(0, true, *prim);
    mPrimMinBounds[1][i] = bound(1, true, *prim);
    mPrimMinBounds[2][i] = bound(2, true, *prim);

    mPrimMaxBounds[0][i] = bound(0, false, *prim);
    mPrimMaxBounds[1][i] = bound(1, false, *prim);
    mPrimMaxBounds[2][i] = bound(2, false, *prim);
  } // end for prim
} // end CachedBounder::CachedBounder()

