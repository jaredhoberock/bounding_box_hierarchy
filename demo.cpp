#include <array>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>

#include "bounding_volume_hierarchy.hpp"
#include "BoundingVolumeHierarchy.h"
#include "optional.hpp"

struct point : std::array<float,3>
{
  point() = default; 
  point(const point&) = default; 
  point(float x, float y, float z) : std::array<float,3>{x,y,z} {}
};

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
  
  vector subtract_me(lhs[2]*rhs[1], lhs[0]*rhs[2], lhs[1]*rhs[0]);
  
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


struct triangle_bounding_box
{
  auto operator()(const triangle& tri) const
  {
    return tri.bounding_box();
  }

  float operator()(const triangle& tri, int axis, bool min) const
  {
    if(min)
    {
      return std::min(tri[0][axis], std::min(tri[1][axis], tri[2][axis]));
    }

    return std::max(tri[0][axis], std::max(tri[1][axis], tri[2][axis]));
  }

  float operator()(int axis, bool min, const triangle& tri) const
  {
    return operator()(tri, axis, min);
  }
};


struct triangle_intersect
{
  const std::vector<triangle>& triangles;

  bool operator()(const triangle& tri, const point& origin, const point& direction, float& t) const
  {
    std::array<float,2> interval{0, 1};
    auto result = tri.intersect(origin, direction, interval);
    if(result) t = *result;
    return bool(result);
  }

  bool operator()(int idx, const point& origin, const point& direction, float& t) const
  {
    return operator()(triangles[idx], origin, direction, t);
  }

  bool operator()(const point& origin, const point& direction, int idx, float& t) const
  {
    return operator()(idx, origin, direction, t);
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

  // build bvhs
  BoundingVolumeHierarchy<triangle, point> old_bvh;
  triangle_bounding_box bound;
  old_bvh.build(triangles, bound);

  bounding_volume_hierarchy<triangle> new_bvh(triangles);

  // generate some random rays
  auto rays = random_rays_in_unit_cube(num_rays, seed + 1);

  std::vector<int> old_intersections;
  for(int i = 0; i < rays.size(); ++i)
  {
    auto& ray = rays[i];

    triangle_intersect intersector{triangles};
    if(old_bvh.intersect(ray.first, ray.second, 0.f, 1.f, intersector))
    {
      old_intersections.push_back(i);
    }
  }

  std::vector<int> new_intersections;
  for(int i = 0; i < rays.size(); ++i)
  {
    auto& ray = rays[i];

    if(new_bvh.intersect(ray.first, ray.second))
    {
      new_intersections.push_back(i);
    }
  }

  if(new_intersections != old_intersections)
  {
    std::cerr << "old_intersections.size(): " << old_intersections.size() << std::endl;
    std::cerr << "new_intersections.size(): " << new_intersections.size() << std::endl;
  }

  assert(new_intersections == old_intersections);
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

