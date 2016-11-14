#pragma once

#include <vector>
#include <array>
#include <stack>
#include <numeric>
#include <functional>
#include <algorithm>
#include <tuple>

#include "memoized_bounder.hpp"
#include "partitioner.hpp"


// XXX the template parameter list might go something like this:
//
// template<class T, class Volume = default_bounding_box_type, class Function = default_bounding_volume_intersector>
template<class T>
class bounding_box_hierarchy
{
  private:
    struct call_member_intersect
    {
      template<class... Args>
      auto operator()(const T& element, Args&&... args) const
      {
        return element.intersect(std::forward<Args>(args)...);
      }
    };


    struct call_member_bounding_box
    {
      auto operator()(const T& element) const
      {
        return element.bounding_box();
      }
    };


    struct default_projection
    {
      float operator()(float x) const
      {
        return x;
      }

      template<class... Types>
      auto operator()(const std::tuple<Types...>& t) const
      {
        return std::get<float>(t);
      }

      template<class T1, class T2>
      auto operator()(const std::pair<T1,T2>& p) const
      {
        return std::get<float>(p);
      }
    };


    struct select_bounding_box_type
    {
      // if U::bounding_box() exists, use its type as the bounding box type
      template<class U>
      static auto test(int) -> decltype(std::declval<U>().bounding_box());

      // otherwise, a bounding box is an array of two arrays of three floats
      template<class>
      static std::array<std::array<float,3>,2> test(...);

      using type = decltype(test<T>(0));
    };


  public:
    using element_type = T;

    using bounding_box_type = typename select_bounding_box_type::type;


    template<class ContiguousRange,
             class Bounder = call_member_bounding_box,
             class Partitioner = partition_largest_axis_at_median_element>
    bounding_box_hierarchy(const ContiguousRange& elements,
                           Bounder bounder = call_member_bounding_box(),
                           Partitioner partitioner = partition_largest_axis_at_median_element())
      : nodes_(make_tree(elements, bounder, partitioner))
    {}


    bounding_box_type bounding_box() const
    {
      return bounding_box(root_node());
    }


    template<class Point, class Vector, class U,
             class Function1 = call_member_intersect,
             class Function2 = default_projection>
    U intersect(Point origin, Vector direction, U init,
                Function1 intersector = call_member_intersect(),
                Function2 hit_time = default_projection()) const
    {
      U result = init;
      auto result_t = hit_time(result);

      Vector one_over_direction = {1.f/direction[0], 1.f/direction[1], 1.f/direction[2]};

      using stack_type = short_stack<const node*,64>;

      stack_type stack;
      stack.push(root_node());

      while(!stack.empty())
      {
        const node* current_node = stack.top();
        stack.pop();

        if(is_leaf(current_node))
        {
          // we pass result to intersector() to implement things like
          // * mailboxing
          // * ray intervals
          auto current_result = intersector(element(current_node), origin, direction, result);
          auto current_t = hit_time(current_result);
          if(current_t < result_t)
          {
            result_t = current_t;
            result = current_result;
          }
        }
        else
        {
          if(intersect_box(bounding_box(current_node), origin, one_over_direction, result_t))
          {
            // push children to stack 
            stack.push(current_node->left_child_);
            stack.push(current_node->right_child_);
          }
        }
      }

      return result;
    }


  private:
    template<class U, size_t N>
    class short_stack : private std::array<U,N>
    {
      public:
        short_stack()
          : top_(this->data())
        {}

        void push(const U& value)
        {
          ++top_;
          *top_ = value;
        }

        U&& pop()
        {
          return std::move(*top_--);
        }

        bool empty() const
        {
          return top_ == this->data();
        }

        U& top()
        {
          return *top_;
        }

        const U& top() const
        {
          return *top_;
        }

      private:
        U* top_;
    };


    template<class Point, class Vector>
    static bool intersect_box(const bounding_box_type& box,
                              Point origin,
                              Vector one_over_direction,
                              float t_bound)
    {
      for(int i = 0; i < 3; ++i)
      {
        float t_near = (box[0][i] - origin[i]) * one_over_direction[i];
        float t_far  = (box[1][i] - origin[i]) * one_over_direction[i];

        if(t_near > t_far) std::swap(t_near, t_far);
        float t0 = t_near > 0.f ? t_near : 0.f;
        float t1 = t_far < t_bound ? t_far : t_bound;
        if(t0 > t1) return false;
      }

      return true;
    }


    template<class IndirectBounder>
    static bounding_box_type bounding_box(const std::vector<size_t>::iterator begin,
                                          const std::vector<size_t>::iterator end,
                                          IndirectBounder bounder)
    {
      float inf = std::numeric_limits<float>::infinity();
      bounding_box_type result{{{inf, inf, inf}, {-inf, -inf, -inf}}};
          
      for(std::vector<size_t>::iterator e = begin; e != end; ++e)
      {
        auto bounding_box = bounder(*e);

        for(int i = 0; i < 3; ++i)
        {
          result[0][i] = std::min(result[0][i], bounding_box[0][i]);
          result[1][i] = std::max(result[1][i], bounding_box[1][i]);
        }
      }

      return result;
    }


    struct node
    {
      const node* left_child_;
      const node* right_child_;
      bounding_box_type bounding_box_;

      node(const node* left_child,
           const node* right_child,
           const bounding_box_type& bounding_box)
        : left_child_(left_child),
          right_child_(right_child),
          bounding_box_(bounding_box)
      {}
    };

    template<class Bounder>
    struct indirect_bounder
    {
      template<class ContiguousRange>
      indirect_bounder(Bounder bounder_, const ContiguousRange& elements)
        : bounder(bounder_),
          data(&*elements.begin())
      {}

      auto operator()(size_t element_idx) const
      {
        return bounder(data[element_idx]);
      }

      Bounder bounder;
      const T* data;
    };

    template<class Bounder, class ContiguousRange>
    static indirect_bounder<Bounder> make_indirect_bounder(Bounder bounder, const ContiguousRange& elements)
    {
      return indirect_bounder<Bounder>(bounder, elements);
    }


    template<class ContiguousRange, class IndirectBounder, class Partitioner>
    static const node* make_tree_recursive(std::vector<node>& tree,
                                           std::vector<size_t>::iterator begin,
                                           std::vector<size_t>::iterator end,
                                           const ContiguousRange& elements,
                                           IndirectBounder bounder,
                                           Partitioner partitioner)
    {
      if(begin + 1 == end)
      {
        // we've hit a leaf, so return a pointer to the element
        return reinterpret_cast<const node*>(&elements[*begin]);
      }

      // find the bounding box of the elements
      bounding_box_type box = bounding_box(begin, end, bounder);

      // partition the elements into two sets
      std::vector<size_t>::iterator split = partitioner(begin, end, box, bounder);

      // build subtrees
      const node* left_child  = make_tree_recursive(tree, begin, split, elements, bounder, partitioner);
      const node* right_child = make_tree_recursive(tree, split, end,   elements, bounder, partitioner);

      // create a new node
      tree.emplace_back(left_child, right_child, box);
      return &tree.back();
    }


    template<class ContiguousRange, class Bounder, class Partitioner>
    static std::vector<node> make_tree(const ContiguousRange& elements, Bounder bounder, Partitioner partitioner)
    {
      // we will sort an array of indices
      std::vector<size_t> indices(elements.size());
      std::iota(indices.begin(), indices.end(), 0);

      // reserve n - 1 nodes to ensure that no iterators are invalidated during construction
      std::vector<node> tree;
      tree.reserve(elements.size() - 1);

      // memoize the bound function
      memoized_bounder<T,Bounder> memoized_bounder(elements, bounder);

      // add indirection to the bound function
      auto indirect_bounder = make_indirect_bounder(std::ref(memoized_bounder), elements);

      // recurse
      make_tree_recursive(tree, indices.begin(), indices.end(), elements, indirect_bounder, partitioner);

      assert(tree.size() == elements.size() - 1);

      return tree;
    }


    const T& element(const node* leaf) const
    {
      return *reinterpret_cast<const T*>(leaf);
    }

    const bounding_box_type& bounding_box(const node* n) const
    {
      return n->bounding_box_;
    }

    bool is_leaf(const node* n) const
    {
      return n < &*nodes_.begin() || &*nodes_.end() <= n;
    }

    const node* root_node() const
    {
      return &nodes_.back();
    }

    std::vector<node> nodes_;
};

