#include "bounding_volume_hierarchy.hpp"
#include <limits>

template<typename T>
  size_t bounding_volume_hierarchy<T>
    ::findPrincipalAxis(const point &min_corner,
                        const point &max_corner)
{
  // find the principal axis of the points
  size_t axis = 4;
  float maxLength = -std::numeric_limits<float>::infinity();
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


template<typename T>
  template<typename Bounder>
    bounding_volume_hierarchy<T>::CachedBounder<Bounder>
      ::CachedBounder(Bounder &bound,
                      const std::vector<T> &primitives)
{
  mPrimMinBounds[0].resize(primitives.size());
  mPrimMinBounds[1].resize(primitives.size());
  mPrimMinBounds[2].resize(primitives.size());

  mPrimMaxBounds[0].resize(primitives.size());
  mPrimMaxBounds[1].resize(primitives.size());
  mPrimMaxBounds[2].resize(primitives.size());

  size_t i = 0;
  for(typename std::vector<T>::const_iterator prim = primitives.begin();
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

