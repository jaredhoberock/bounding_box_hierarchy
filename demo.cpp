#include <array>
#include <random>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <chrono>

#include "bounding_volume_hierarchy.hpp"
#include "bounding_box_hierarchy.hpp"
#include "exhaustive_searcher.hpp"
#include "time_invocation.hpp"

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
  float intersect(const point& origin, const point& direction, float nearest) const
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
      return nearest;
    }

    float inv_divisor = 1.f / divisor;

    // compute barycentric coordinates 
    vector d = origin - p0;
    float b0 = dot(d,s1) * inv_divisor;
    if(b0 < 0.f || b0 > 1.f)
    {
      return nearest;
    }

    vector s2 = cross(d,e1);
    float b1 = dot(direction, s2) * inv_divisor;
    if(b1 < 0.f || b0 + b1 > 1.f)
    {
      return nearest;
    }

    // compute t
    float t = inv_divisor * dot(e2,s2);
    if(t < 0)
    {
      return nearest;
    }

    return std::min(nearest, t);
  }

  std::experimental::optional<float> intersect(const point& origin, const point& direction, const std::array<float,2>& interval) const
  {
    float t = intersect(origin, direction, interval[1]);
    return t < interval[1] ? std::experimental::make_optional(t) : std::experimental::nullopt;
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


using ray = std::pair<point,vector>;


std::vector<ray> random_rays_in_unit_cube(size_t n, int seed = 13)
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


using intersection_type = std::pair<float, const triangle*>;

template<class Searcher>
std::vector<intersection_type> find_intersections(const std::vector<triangle>& triangles, const std::vector<ray>& rays)
{
  // build searcher
  Searcher searcher(triangles);

  std::vector<intersection_type> intersections;
  for(int i = 0; i < rays.size(); ++i)
  {
    auto& ray = rays[i];

    intersection_type init(1.f, nullptr);

    // use a custom intersection functor to return the hit time and a pointer to the triangle
    auto intersection = searcher.intersect(ray.first, ray.second, init, [](const auto& tri, const auto& o, const auto& d, intersection_type nearest)
    {
      return intersection_type(tri.intersect(o,d,nearest.first), &tri);
    });

    if(intersection.first < 1.f)
    {
      intersections.push_back(intersection);
    }
  }

  return std::move(intersections);
}


template<class Searcher>
bool test(const std::vector<triangle>& triangles, const std::vector<ray>& rays)
{
  return find_intersections<Searcher>(triangles, rays) == find_intersections<exhaustive_searcher<triangle>>(triangles, rays);
}


template<class Hierarchy>
double measure_performance(const Hierarchy& hierarchy, const std::vector<ray>& rays)
{
  // warm up
  for(const auto& r : rays)
  {
    hierarchy.intersect(r.first, r.second, 1.f);
  }

  size_t milliseconds = time_invocation_in_milliseconds(20, [&]
  {
    for(const auto& r : rays)
    {
      hierarchy.intersect(r.first, r.second, 1.f);
    }
  });

  return double(1000 * rays.size()) / milliseconds;
}


int main()
{
  for(size_t i = 0; i < 20; ++i)
  {
    std::mt19937 rng(i);
    std::uniform_int_distribution<> dist(500, 1000);

    int m = dist(rng);
    int n = dist(rng);

    auto triangles = random_triangles_in_unit_cube(dist(rng), rng());
    auto rays = random_rays_in_unit_cube(dist(rng), rng());

    std::cout << "testing " << m << " " << n << std::endl;

    assert(test<bounding_box_hierarchy<triangle>>(triangles, rays));
  }

  size_t num_triangles = 100000;
  size_t num_rays = 1 << 10;

  std::mt19937 rng;
  auto triangles = random_triangles_in_unit_cube(num_triangles, rng());
  auto rays = random_rays_in_unit_cube(num_rays, rng());

  std::cout << "timing exhaustive_searcher: " << std::endl;
  exhaustive_searcher<triangle> es(triangles);
  auto es_rays_per_second = measure_performance(es, rays);
  std::cout << "exhaustive_searcher: " << es_rays_per_second << " rays/s" << std::endl;

  std::cout << "timing bounding_box_hierarchy: " << std::endl;
  bounding_box_hierarchy<triangle> bbh(triangles);
  auto bbh_rays_per_second = measure_performance(bbh, rays);
  std::cout << "bounding_box_hierarchy: " << bbh_rays_per_second << " rays/s" << std::endl;

  std::cout << "OK" << std::endl;

  return 0;
}

