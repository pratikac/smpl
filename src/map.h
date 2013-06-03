#ifndef __map_h__
#define __map_h__

#define RANDF   (rand()/(RAND_MAX+1.0))

template <size_t N> 
class map_c
{
  public:
    map_c() {}
    ~map_c() {}

    virtual int sample_free_space(float[N])
    {
      return 1;
    };
    virtual bool is_in_collision(const float[N])
    {
      return false;
    };
    virtual float get_state_cost(const float s[N])
    { 
      return 0;
    }
};

#endif
