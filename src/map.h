#ifndef __map_h__
#define __map_h__

template <size_t N> 
class map_c
{
  public:
    map_c() {}
    ~map_c() {}

    virtual int sample_free_space(double[N])
    {
      return 1;
    }
    virtual bool is_in_collision(const double[N])
    {
      return false;
    }
    virtual double get_state_cost(const double s[N])
    { 
      return 0;
    }
};

#endif
