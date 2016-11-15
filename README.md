# bounding_box_hierarchy
Spatial data structures such as *bounding box hierarchies* accelerate ray
intersection queries by disqualifying entire collections of objects from
consideration in bulk. If a ray does not intersect the *bounding box* of a
collection of objects, then those objects need not be considered individually.
Organizing bounding boxes into a *hierarchical* tree structure in turn allows
quick disqualification of entire subtrees.

The code in this repository implements a `bounding_box_hierarchy` suitable for accelerating
ray intersection queries against generic objects.

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
};

std::vector<fancy_triangle> triangles;

bounding_box_hierarchy<fancy_triangle> bbh(triangles);

```

TODO: explain how `.intersect()` works here

