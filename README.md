# bounding_box_hierarchy
Spatial data structures such as *bounding box hierarchies* accelerate ray
intersection queries by disqualifying entire collections of objects from
consideration in bulk. If a ray does not intersect the *bounding box* of a
collection of objects, then those objects need not be considered individually.
Organizing bounding boxes into a *hierarchical* tree structure in turn allows
quick disqualification of entire subtrees.

The code in this repository implements a `bounding_box_hierarchy` suitable for accelerating
ray intersection queries against generic objects.

## Construction

For example, a `bounding_box_hierarchy` can be constructed from a collection of `triangle`s:

```
struct triangle
{
  ...
};

std::vector<triangle> triangles;

bounding_box_hierarchy<triangle> bbh(triangle, bounder);
```

A `bounder` is a function that returns a `triangle`'s bounding box. It
takes a reference to a `triangle` as a parameter and returns its bounding box:

```
using point = std::array<float,3>;
using bounding_box_type = std::array<point,2>;

bounding_box_type bounder(const triangle& tri)
{
  ...
}
```

In this example, a `bounding_box` is an array of two `point`s. The first
(second) element of the `bounding_box` is understood to be the "minimal"
("maximal") corner of the box.

In general, the `bounder` parameter of `bounding_box_hierarchy`'s constructor
can be any callable type that receives a `T` as a parameter and returns a
bounding box type that can be indexed as described above.

The `bounder` parameter may be omitted if the type `T` has a member function
named `bounding_box()` that returns a bounding box type:

```
using point = std::array<float,3>;
using bounding_box_type = std::array<point,2>;

struct fancy_triangle
{
  bounding_box_type bounding_box() const
  {
    ...
  }

  // other data here
  ...
};

std::vector<fancy_triangle> triangles;

bounding_box_hierarchy<fancy_triangle> bbh(triangles);

```

## Intersection

After construction, a `bounding_box_hierarchy` can be queried for intersections with rays with the `.intersect()` member function:

```
template<class Point, class Vector, class Result, class Intersector>
Result intersect(Point origin, Vector direction, Result init, Intersector intersector);
```

It returns the nearest intersection with the ray beginning at `origin` pointing in `direction`, if one exists. Otherwise, `.intersect()` returns `init`.

The `intersector` parameter is a callable type with the following signature:

```
Result intersector(const T& element, Point origin, Vector direction, Result init);
```

The role of `intersector` is to test for an intersection between the given
`element` and ray. If an intersection exists, it should return the details in
the result. Otherwise, it should return `init`. These details can be arbitrary, but they should somehow encode the parametric hit "time" along the ray with the element.

For example, a simple `intersector` can be defined as a lambda function and can simply return the hit time as a floating point number:

```
bounding_box_hierarchy<triangle> bbh = ...

using point  = std::array<float,3>;
using vector = std::array<float,3>;

// create a ray beginning at the origin pointing in the positive z direction
point  ray_origin{0,0,0};
vector ray_direction{0,0,1};

// only consider intersections in the parametric ray interval [0,1)
float init = 1.0f;

float hit_time = bbh.intersect(ray_origin, ray_direction, init, [](const triangle& tri, point o, vector d, float nearest_hit)
{
  // compute the "time" at which the ray defined by o & d hits tri
  float t = ...

  // return the nearer of the two hit times
  return std::min(t,nearest_hit)
});

if(hit_time < init)
{
  std::cout << "intersection at time " << hit_time << std::endl;
}
else
{
  std::cout << "no intersection" << std::endl;
}
```

### Custom Intersection Results

Usually, we need to know which triangle was intersected, not just the hit time. In this case, we can return both values from the `intersector`:

```
// define a custom intersection result type
using intersection = std::pair<float, const triangle*>;

// only consider intersections in the parametric ray interval [0,1), and initialize the nearest triangle to null
intersection init(1.f, nullptr);

intersection result = bbh.intersect(ray_origin, ray_direction, init, [](const triangle& tri, point o, vector d, const intersection_type& nearest_result)
{
  // compute the time at which the ray defined by o & d hits tri
  float t = ...

  // get the hit time of the nearest result encountered so far
  float nearest_t = nearest_result.first;

  // return the nearer of the two results
  return t < nearest_t ? nearest_result : intersection(t, &tri);
});

if(result.second)
{
  std::cout << "intersection with triangle " << result.second << " at time " << result.first << std::endl;
}
else
{
  std::cout << "no intersection" << std::endl;
}
```

When the `result` of `intersector` is not a floating point number, it tries to
call `std::get<float>(result)` to determine the hit time at the intersection.
`std::get<float>(result)` works for `std::pair<T1,T2>` as long as either `T1`
or `T2` is a `float` (but not both).

For other types of intersection results, we can use a custom `hit_time` function, which is passed as a parameter following the `intersector`:

```
// define a fancy intersection result type
struct fancy_intersection
{
  float hit_time;
  const triangle* hit_triangle;

  // other info might go here
  ...
};

// only consider intersections in the parametric ray interval [0,1), and initialize the nearest triangle to null
fancy_intersection fancy_init(1.f, nullptr, other data...)

fancy_intersection result = bbh.intersect(ray_origin, ray_direction, fancy_init, [](const triangle& tri, point o, vector d, const fancy_intersection& nearest_result)
{
  // compute the time at which the ray defined by o & d hits tri
  float t = ...

  // get the hit time of the nearest result encountered so far
  float nearest_t = nearest_result.hit_time;

  // return the nearer of the two results
  return t < nearest_t ? nearest_result : fancy_intersection(t, &tri, other data...);
},
[](const fancy_intersection& i)
{
  // the hit_time() function simply returns the .hit_time member of i
  return t.hit_time;
});

if(result.hit_triangle)
{
  std::cout << "intersection with triangle " << result.hit_triangle << " at time " << result.hit_time << std::endl;
}
else
{
  std::cout << "no intersection" << std::endl;
}
```

### Even Fancier

If our geometric object types are fancy enough, we can avoid passing the `intersector` and `hit_time` functions to `.intersect()`: 

```
struct fancier_triangle
{
  // define the type of .intersect()'s result
  struct intersection
  {
    // define a member function .hit_time() which returns the hit time of this intersection
    float hit_time() const { return time; }

    float time;
    const triangle* tri;

    // info might go here
    ...
  };

  // define a member function .intersect() which returns the intersection between this
  // fancier_triangle and the ray, if it exists and is nearer than nearest
  intersection intersect(const point& o, const vector& d, const intersection& nearest) const
  {
    // compute the time at which the ray defined by o & d hits tri
    float t = ...

    // get the hit time of the nearest result encountered so far
    float nearest_t = nearest.hit_time();

    // return the nearer of the two results
    return t < nearest_t ? nearest : intersection(t, this, other data...);
  }

  // define a member function .bounding_box() which returns this fancier_triangle's bounding box
  bounding_box_type bounding_box() const
  {
    ...
  }

  // other data here
  ...
};

std::vector<fancier_triangle> triangles;

bounding_box_hierarchy<fancier_triangle> bbh(triangles);

// only consider intersections in the parametric ray interval [0,1), and initialize the nearest triangle to null
fancier_triangle::intersection fancier_init(1.f, nullptr, other data...)

// because fancier_triangle has the member function .intersect() and the result type of .intersect() has the
// member function .hit_time(), we don't need to pass the intersector or hit_time functions to bbh.intersect()
auto result = bbh.intersect(ray_origin, ray_direction, fancier_init);

if(result.tri)
{
  std::cout << "intersection with triangle " << result.tri << " at time " << result.hit_time() << std::endl;
}
else
{
  std::cout << "no intersection" << std::endl;
}
```

The [demo](./demo.cpp) program demonstrates these techniques.

