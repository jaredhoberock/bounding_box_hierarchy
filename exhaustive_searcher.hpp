#pragma once

#include <utility>
#include <tuple>

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


  public:
    using element_type = T;

    template<class ContiguousRange>
    exhaustive_searcher(const ContiguousRange& elements)
      : begin_(&*elements.begin()),
        end_(&*elements.end())
    {}

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
    const T* begin_;
    const T* end_;
};

