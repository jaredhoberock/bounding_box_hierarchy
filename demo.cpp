#include <array>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

#include "bounding_volume_hierarchy.hpp"

using point = std::array<float,3>;

point operator-(const point& lhs, const point& rhs)
{
  return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

using vector = point;

float dot(const vector& lhs, const vector& rhs)
{
  return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), 0.f);
}

vector cross(const vector& lhs, const vector& rhs)
{
  vector result;
  
  vector subtract_me{lhs[2]*rhs[1], lhs[0]*rhs[2], lhs[1]*rhs[0]};
  
  result[0]  = (lhs[1] * rhs[2]);
  result[0] -= subtract_me[0];
  result[1]  = (lhs[2] * rhs[0]);
  result[1] -= subtract_me[1];
  result[2]  = (lhs[0] * rhs[1]);
  result[2] -= subtract_me[2];
  
  return result;
}


struct triangle : std::array<point,3>
{
  std::experimental::optional<float> intersect(const point& origin, const point& direction, const std::array<float,2>& interval) const
  {
    const point& p0 = (*this)[0];
    const point& p1 = (*this)[1];
    const point& p2 = (*this)[2];

    vector e1 = p1 - p0;
    vector e2 = p2 - p0;
    vector s1 = cross(direction,e2);
    float divisor = dot(s1,e1);
    if(divisor == 0.f)
    {
      return std::experimental::nullopt;
    }

    float inv_divisor = 1.f / divisor;

    // compute barycentric coordinates 
    vector d = origin - p0;
    float b0 = dot(d,s1) * inv_divisor;
    if(b0 < 0.f || b0 > 1.f)
    {
      return std::experimental::nullopt;
    }

    vector s2 = cross(d,e1);
    float b1 = dot(direction, s2) * inv_divisor;
    if(b1 < 0.f || b0 + b1 > 1.f)
    {
      return std::experimental::nullopt;
    }

    // compute t
    float t = inv_divisor * dot(e2,s2);

    return (interval[0] <= t && t < interval[1]) ? std::experimental::make_optional(t) : std::experimental::nullopt;
  }

  std::array<point,2> bounding_box() const
  {
    point min_corner;

    min_corner[0] = std::min((*this)[0][0], std::min((*this)[1][0], (*this)[2][0]));
    min_corner[1] = std::min((*this)[0][1], std::min((*this)[1][1], (*this)[2][1]));
    min_corner[2] = std::min((*this)[0][2], std::min((*this)[1][2], (*this)[2][2]));

    point max_corner;

    max_corner[0] = std::max((*this)[0][0], std::max((*this)[1][0], (*this)[2][0]));
    max_corner[1] = std::max((*this)[0][1], std::max((*this)[1][1], (*this)[2][1]));
    max_corner[2] = std::max((*this)[0][2], std::max((*this)[1][2], (*this)[2][2]));

    return {min_corner, max_corner};
  }
};


std::vector<triangle> random_triangles_in_unit_cube(size_t n, int seed = 0)
{
  std::mt19937 rng(seed);
  std::uniform_real_distribution<float> unit_interval(0,1);

  std::vector<triangle> result(n);
  for(triangle& tri : result)
  {
    for(int i = 0; i < 3; ++i)
    {
      tri[i] = {unit_interval(rng), unit_interval(rng), unit_interval(rng)};
    }
  }

  return result;
}


std::vector<std::pair<point,vector>> random_rays_in_unit_cube(size_t n, int seed = 13)
{
  std::default_random_engine rng(seed);
  std::uniform_real_distribution<float> unit_interval(0,1);

  std::vector<std::pair<point,vector>> result(n);
  for(auto& ray : result)
  {
    point from;
    point to;

    from = {unit_interval(rng), unit_interval(rng), unit_interval(rng)};
    to = {unit_interval(rng), unit_interval(rng), unit_interval(rng)};

    vector direction = to - from;

    ray.first = from;
    ray.second = direction;
  }

  return result;
}


void test(size_t num_triangles, size_t num_rays, size_t seed = 0)
{
  // generate some random triangles
  auto triangles = random_triangles_in_unit_cube(num_triangles, seed);

  // build bvh
  bounding_volume_hierarchy<triangle> bvh(triangles);

  // generate some random rays
  auto rays = random_rays_in_unit_cube(num_rays, seed + 1);

  using intersection_type = std::pair<const triangle*,float>;

  std::vector<intersection_type> bvh_intersections;
  for(int i = 0; i < rays.size(); ++i)
  {
    auto& ray = rays[i];

    // use a custom intersection functor to return a pointer to the triangle and the hit time
    auto intersection = bvh.intersect(ray.first, ray.second, {0,1}, [](const auto& tri, const auto& o, const auto& d, const auto& i)
    {
      std::experimental::optional<float> intermediate_result = tri.intersect(o,d,i);

      if(intermediate_result)
      {
        auto result = intersection_type(&tri, *intermediate_result);
        return std::experimental::make_optional(result);
      }

      return std::experimental::optional<intersection_type>();
    },
    [](const intersection_type& result)
    {
      // the hit time is the intersection's second half
      return result.second;
    });

    if(intersection)
    {
      bvh_intersections.push_back(*intersection);
    }
  }

  std::vector<intersection_type> reference;
  for(int i = 0; i < rays.size(); ++i)
  {
    auto& ray = rays[i];

    std::experimental::optional<intersection_type> nearest_intersection;

    for(const auto& tri : triangles)
    {
      std::experimental::optional<float> current_intersection = tri.intersect(ray.first, ray.second, {0,1});

      if(current_intersection)
      {
        if(!nearest_intersection || *current_intersection < nearest_intersection->second)
        {
          nearest_intersection = intersection_type(&tri,*current_intersection);
        }
      }
    }

    if(nearest_intersection)
    {
      reference.push_back(*nearest_intersection);
    }
  }

  if(bvh_intersections != reference)
  {
    std::cerr << "bvh_intersections.size(): " << bvh_intersections.size() << std::endl;
    std::cerr << "reference.size(): " << reference.size() << std::endl;
  }

  assert(bvh_intersections == reference);
}


int main()
{
  for(size_t i = 0; i < 20; ++i)
  {
    std::mt19937 rng(i);
    std::uniform_int_distribution<> dist(500, 1000);

    int m = dist(rng);
    int n = dist(rng);

    std::cout << "testing " << m << " " << n << std::endl;

    test(m, n, rng());
  }

  std::cout << "OK" << std::endl;

  return 0;
}

