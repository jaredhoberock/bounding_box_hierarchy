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

