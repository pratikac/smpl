#ifndef __system_h__
#define __system_h__

#include "dynamical_system.h"
#include "map.h"
#include <cmath>
#include <cfloat>
#include <vector>
#include <ostream>
#include <string.h>
#include <cstdlib>
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
    
    virtual bool is_inside(const state_c<N>& sin) const
    {
      for(size_t i=0; i<N; i++)
      {
        if(fabs(sin.x[i] - c[i]) > s[i]/2.0)
          return false;
      }
      return true;
    }
};

template<size_t dim_t>
class cost_c
{
  public:
    const static size_t dim = dim_t;
    vector<float> val;
    
    cost_c(float x=FLT_MAX/2){
      val = vector<float>(dim, x);
    }
    cost_c(vector<float>& xin) : val(xin) {}
    cost_c(const cost_c& c2) : val(c2.val){}
    cost_c(float c, int d){
      val[d] = c;
    }
    virtual ~cost_c(){}

    float& operator[](int d){
      return val[d];
    }
    virtual cost_c& operator+=(const cost_c& rhs)
    {
      for(size_t i=0; i<dim; i++)
        val[i] += rhs.val[i];
      return *this;
    }
    virtual cost_c operator+(const cost_c& c2) const
    {
      cost_c toret = *this;
      for(size_t i=0; i<dim; i++)
        toret.val[i] += c2.val[i];
      return toret;
    }
    virtual bool operator<(const cost_c& rhs) const
    {
      for(size_t i=0; i<dim; i++){
        if(val[i] < rhs.val[i])
          return true;
        else if(val[i] > rhs.val[i])
          return false;
      }
      return true;
    }
    virtual bool operator>(const cost_c& rhs) const
    {
      return !(*this < rhs);
    }
    virtual float difference(const cost_c& c2) const
    {
      float t1 = 0;
      for(size_t i=0; i<dim; i++)
        t1 += (val[i]-c2.val[i])*(val[i]-c2.val[i]);
      return sqrt(t1);
    }
    virtual ostream& print(ostream& os=cout, const char* prefix=NULL, const char* suffix=NULL) const
    {
      if(prefix)
        os<<prefix;
      for(auto& vi : val)
        os<<vi<<",";
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
      extend_cost = cost_t(total_variation, extend_cost.dim-1);
      if(total_variation > 0)
        return 0;
      return 1;
    }

    virtual cost_t get_state_cost(const state& s)
    {
      cost_t t1;
      float& t2 = t1.val.back();
      t2 = obstacle_map.get_state_cost(s.x);
      return t1;
    };

    virtual cost_t evaluate_trajectory_cost(trajectory& traj)
    {
      cost_t c(0);
      for(auto& ps : traj.states)
        c += get_state_cost(state(ps));
      return c;
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
