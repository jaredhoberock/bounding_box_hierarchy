#include "bounding_volume_hierarchy.hpp"
#include "BoundingVolumeHierarchy.h"
#include <array>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>

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

using triangle = std::array<point,3>;

// this should just return a pair of points (min, max)
float bound_triangle(int axis, bool min, const triangle& tri)
{
  if(min)
  {
    return std::min(tri[0][axis], std::min(tri[1][axis], tri[2][axis]));
  }

  return std::max(tri[0][axis], std::min(tri[1][axis], tri[2][axis]));
}


struct intersect_triangle
{
  const std::vector<triangle>& triangles;

  bool operator()(const point& origin, const point& direction, int idx, float& t) const
  {
    const triangle& tri = triangles[idx];
    const point& p0 = tri[0];
    const point& p1 = tri[1];
    const point& p2 = tri[2];

    vector e1 = p1 - p0;
    vector e2 = p2 - p0;
    vector s1 = cross(direction,e2);
    float divisor = dot(s1,e1);
    if(divisor == 0.f)
    {
      return false;
    }

    float inv_divisor = 1.f / divisor;

    // compute barycentric coordinates 
    vector d = origin - p0;
    float b0 = dot(d,s1) * inv_divisor;
    if(b0 < 0.f || b0 > 1.f)
    {
      return false;
    }

    vector s2 = cross(d,e1);
    float b1 = dot(direction, s2) * inv_divisor;
    if(b1 < 0.f || b0 + b1 > 1.f)
    {
      return false;
    }

    // compute t
    t = inv_divisor * dot(e2,s2);

    return true;
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
  old_bvh.build(triangles, bound_triangle);

  bounding_volume_hierarchy<triangle, point> new_bvh;
  new_bvh.build(triangles, bound_triangle);

  // generate some random rays
  auto rays = random_rays_in_unit_cube(num_rays, seed + 1);

  std::vector<int> old_intersections;
  for(int i = 0; i < rays.size(); ++i)
  {
    auto& ray = rays[i];

    intersect_triangle intersector{triangles};
    if(old_bvh.intersect(ray.first, ray.second, 0.f, 1.f, intersector))
    {
      old_intersections.push_back(i);
    }
  }

  std::vector<int> new_intersections;
  for(int i = 0; i < rays.size(); ++i)
  {
    auto& ray = rays[i];

    intersect_triangle intersector{triangles};
    if(new_bvh.intersect(ray.first, ray.second, 0.f, 1.f, intersector))
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

