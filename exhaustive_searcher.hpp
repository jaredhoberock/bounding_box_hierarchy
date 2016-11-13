#pragma once

#include <utility>
#include <tuple>
#include <array>
#include <limits>


template<class T>
class exhaustive_searcher
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

    template<class ContiguousRange, class Bounder = call_member_bounding_box>
    exhaustive_searcher(const ContiguousRange& elements,
                        Bounder bounder = call_member_bounding_box(),
                        float epsilon = std::numeric_limits<float>::epsilon())
      : begin_(&*elements.begin()),
        end_(&*elements.end()),
        bounding_box_(bounding_box(elements, bounder, epsilon))
    {}

    bounding_box_type bounding_box() const
    {
      return bounding_box_;
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

      for(const T* element = begin_; element != end_; ++element)
      {
        auto current_result = intersector(*element, origin, direction, result);
        auto current_t = hit_time(current_result);

        if(current_t < result_t)
        {
          result = current_result;
          result_t = current_t;
        }
      }

      return result;
    }

  private:
    template<class ContiguousRange, class Bounder>
    static bounding_box_type bounding_box(const ContiguousRange& elements, Bounder bounder, float epsilon)
    {
      float inf = std::numeric_limits<float>::infinity();
      bounding_box_type result{{{inf, inf, inf}, {-inf, -inf, -inf}}};
          
      for(const auto& element : elements)
      {
        auto bounding_box = bounder(element);

        for(int i = 0; i < 3; ++i)
        {
          result[0][i] = std::min(result[0][i], bounding_box[0][i]);
          result[1][i] = std::max(result[1][i], bounding_box[1][i]);
        }
      }

      // widen the bounding box by 2*epsilon
      // XXX might want to instead ensure that each side is at least epsilon in width
      // this ensures that axis-aligned elements always
      // lie strictly within the bounding box
      for(size_t i = 0; i != 3; ++i)
      {
        result[0][i] -= epsilon;
        result[1][i] += epsilon;
      }

      return result;
    }

    const T* begin_;
    const T* end_;
    bounding_box_type bounding_box_;
};

