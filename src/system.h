#ifndef __system_h__
#define __system_h__

#include "dynamical_system.h"
#include "map.h"
#include <cmath>
#include <cfloat>
#include <vector>
#include <ostream>
#include <string.h>
using namespace std;


template<size_t N>
class region_c
{
  public:
    float s[N];
    float c[N];

    region_c()
    {
      for(size_t i=0; i<N; i++)
        s[i] = c[i] = 0;
    }
    
    region_c(const float* cin, const float* sin)
    {
      for(size_t i=0; i< N; i++)
      {
        s[i] = sin[i];
        c[i] = cin[i];
      }
    }
    
    bool is_inside(const state_c<N>& sin) const
    {
      for(size_t i=0; i<N; i++)
      {
        if(fabs(sin.x[i] - c[i]) > s[i]/2.0)
          return false;
      }
      return true;
    }
};

class cost_c
{
  public:
    float val;
    cost_c() : val(FLT_MAX/2){}
    cost_c(float xin) : val(xin) {}
    cost_c(const cost_c& c2) : val(c2.val){}
    
    virtual ~cost_c(){}

    virtual cost_c& operator+=(const cost_c& rhs)
    {
      val += rhs.val;
      return *this;
    }
    virtual cost_c operator+(const cost_c& c2) const
    {
      cost_c toret = *this;
      toret.val += c2.val;
      return toret;
    }
    virtual bool operator<(const cost_c& rhs) const
    {
      return (val < rhs.val);
    }
    virtual bool operator>(const cost_c& rhs) const
    {
      return (val > rhs.val);
    }
    virtual float difference(const cost_c& c2) const
    {
      return fabs(val-c2.val);
    }
    virtual ostream& print(ostream& os=cout, const char* prefix=NULL, const char* suffix=NULL) const
    {
      if(prefix)
        os<<prefix;
      os<<val;
      if(suffix)
        os<<suffix;
      return os;
    }
};


template<class dynamical_system_tt, class map_tt, class cost_tt>
class system_c
{
  public:
    typedef dynamical_system_tt dynamical_system_t;
    typedef map_tt map_t;
    typedef cost_tt cost_t;

    typedef typename dynamical_system_t::state_t state;
    typedef typename dynamical_system_t::control_t control;
    typedef typename dynamical_system_t::opt_data_t opt_data_t;
    typedef typename dynamical_system_t::trajectory_t trajectory;
    const static size_t N = dynamical_system_t::state_t::N;
    typedef region_c<N> region_t;
  
    map_t obstacle_map;
    dynamical_system_t dynamical_system;

    region_c<N> operating_region;
    region_c<N> goal_region;

    system_c(){};
    ~system_c(){}

    virtual int get_key(const state& s, float* key)
    {
      for(size_t i=0; i<N; i++)
        key[i] = (s.x[i] - operating_region.c[i])/operating_region.s[i] + 0.5;
      //cout<<"key: "<< key[0]<<" "<<key[1]<<" "<<key[2]<<endl;
      return 0;
    }

    virtual bool is_in_collision(const state& s)
    {
      return obstacle_map.is_in_collision(s.x);
    }
    float get_goal_cost(const state& s)
    {
      state goal_state(goal_region.c);
      return s.dist(goal_state);
    }
    virtual bool is_in_goal(const state& s)
    {
      return goal_region.is_inside(s);
    }
    
    virtual int sample_state(state& s, bool sample_from_free=false)
    {
      if(sample_from_free)
      {
        if(obstacle_map.sample_free_space(s.x))
          return 1;
      }
      else
      {
        bool found_free_state = false;
        while(!found_free_state)
        {
          for(size_t i=0; i<N; i++)
            s.x[i] = operating_region.c[i] + (RANDF-0.5)*operating_region.s[i];
          found_free_state = !is_in_collision(s);
        }
      }
      return 0;
    }
    virtual int sample_in_goal(state& s)
    {
      bool found_free_state = false;
      while(!found_free_state)
      {
        for(size_t i=0; i<N; i++)
          s.x[i] = goal_region.c[i] + (RANDF-0.5)*goal_region.s[i];
        found_free_state = !is_in_collision(s);
      }
      return 0;
    }
    
    int copy_array(const float* xin, float* xout, int dim)
    {
      memcpy(xout, xin, sizeof(float)*dim);
      return 0;
    }
    virtual bool is_safe_trajectory(const trajectory& traj)
    {
      if(traj.states.empty())
        return true;

      state sm;
      int drop_counter=0;
      for(auto& x : traj.states)
      {
        if(!drop_counter)
        {
          sm = state(x);
          if(is_in_collision(sm))
            return false;
        }
        drop_counter++;
        if(drop_counter == 25)
          drop_counter = 0;
      }
      return true;
    }

    virtual int extend_to(const state& si, const state& sf, bool check_obstacles,
        trajectory& traj, opt_data_t& opt_data)
    {
      int res = dynamical_system.extend_to(si, sf, traj, opt_data);
      if(res)
      {
        traj.clear();
        return res;
      }
      if(check_obstacles)
      {
        if(is_safe_trajectory(traj))
          res = 0;
        else
        {
          traj.clear();
          res = 1;
        }
      }
      return res;
    }

    virtual int evaluate_extend_cost(const state& si, const state& sf,
        opt_data_t& opt_data, cost_t& extend_cost)
    {
      float total_variation = dynamical_system.evaluate_extend_cost(si, sf, opt_data);
      extend_cost = cost_t(total_variation);
      if(total_variation > 0)
        return 0;
      return 1;
    }

    virtual float get_state_cost(const state& s)
    {
      return obstacle_map.get_state_cost(s.x);
    };

    virtual cost_t evaluate_trajectory_cost(trajectory& traj)
    {
      float c=0;
      for(auto& ps : traj.states)
        c += get_state_cost(state(ps));
      return cost_t(c);
    }
   
    virtual cost_t get_zero_cost()
    {
      return cost_t(0);
    }
    virtual cost_t get_inf_cost()
    {
      return cost_t(FLT_MAX/2);
    }

    void test_extend_to()
    {
      trajectory traj;
      float zero[3] ={0};
      state origin(zero);
      for(int i=0; i<10; i++)
      {
        state sr;
        sample_in_goal(sr);
        sr.print(cout, "sampled:","\n");
        
        opt_data_t opt_data;
        extend_to(origin, sr, true, traj, opt_data);
        cout<<"cost: "<< traj.total_variation<<endl;
        traj.print();
        getchar();
      }
      traj.clear();
    }
};

#endif
