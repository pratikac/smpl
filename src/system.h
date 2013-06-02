#ifndef __system_h__
#define __system_h__

#include "dynamical_system.h"
#include "map.h"
#include <cmath>
#include <cfloat>
#include <vector>
#include <ostream>
using namespace std;


template<size_t N>
class region_t
{
  public:
    float s[N];
    float c[N];

    region_t()
    {
      for(size_t i=0; i<N; i++)
        s[i] = c[i] = 0;
    }
    region_t(const float* cin, const float* sin)
    {
      for(size_t i=0; i< N; i++)
      {
        s[i] = sin[i];
        c[i] = cin[i];
      }
    }
    bool is_inside(const state_t<N>& sin) const
    {
      for(size_t i=0; i<N; i++)
      {
        if(fabs(sin.x[i] - c[i]) > s[i]/2.0)
          return false;
      }
      return true;
    }
};

// generic cost structure
class cost_t
{
  public:
    float val;
    cost_t() : val(FLT_MAX){}
    cost_t(float xin) : val(xin) {}
    cost_t(const cost_t& c2) : val(c2.val){}
    
    virtual ~cost_t(){
    }
    virtual cost_t& operator+=(const cost_t& rhs)
    {
      val += rhs.val;
      return *this;
    }
    virtual cost_t operator+(const cost_t& c2) const
    {
      cost_t toret = *this;
      toret.val += c2.val;
      return toret;
    }
    virtual bool operator<(const cost_t& rhs) const
    {
      return (val < rhs.val);
    }
    virtual bool operator>(const cost_t& rhs) const
    {
      return (val > rhs.val);
    }
    virtual float difference(const cost_t& c2) const
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

template<class dynamical_system_t, class map_t>
class system_t
{
  public:
    typedef dynamical_system_t::s_t s_t;
    typedef dynamical_system_t::c_t c_t;
    typedef dynamical_system_t::opt_data_t opt_data_t;
    typedef dynamical_system_t::traj_t traj_t;
    typedef state_t::N N;

    map_t* obstacle_map;
    dynamical_system_t<s_t, c_t, opt_data_t>* dynamical_system;

    region_t<N> operating_region;
    region_t<N> goal_region;

    system_t(){};
    ~system_t(){}

    virtual int get_key(const state& s, float* key)
    {
      for(size_t i=0; i<N; i++)
        key[i] = (s.x[i] - operating_region.c[i])/operating_region.s[i] + 0.5;
      //cout<<"key: "<< key[0]<<" "<<key[1]<<" "<<key[2]<<endl;
      return 0;
    }

    virtual bool is_in_collision(const state& s)
    {
      return obstacle_map->is_in_collision(s.x);
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
    
    virtual state* sample_state(bool sample_from_free=false)
    {
      //cout << "system_t::sample_state() called" << endl << flush;
      if(sample_from_free)
      {
        float ps[N] ={0};
        if(obstacle_map->sample_free_space(ps))
          return NULL;
        else
          return new state(ps);
      }
      else
      {
        //cout << "sample_from_free is FALSE" << endl << flush;
        state* spos = new state();
        bool found_free_state = false;
        while(!found_free_state)
        {
          //cout << "not in free space...retrying" << endl << flush;
          for(size_t i=0; i<N; i++)
            spos->x[i] = operating_region.c[i] + (RANDF-0.5)*operating_region.s[i];
          found_free_state = !is_in_collision(*spos);
        }
        //cout << "Now in free space...exiting.." << endl << flush;
        return spos;
      }
      return NULL;
    }
    virtual state* sample_in_goal()
    {
      state* ps = new state();
      bool found_free_state = false;
      while(!found_free_state)
      {
        for(size_t i=0; i<N; i++)
          ps->x[i] = goal_region.c[i] + (RANDF-0.5)*goal_region.s[i];
        found_free_state = !is_in_collision(*ps);
      }
      return ps;
    }
    
    int copy_array(const float* xin, float* xout, int dim)
    {
      for(int i=0; i<dim; i++)
        xout[i] = xin[i];
      return 0;
    }
    virtual bool is_safe_trajectory(const trajectory_t& traj)
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

    virtual int extend_to(const state* si, const state* sf, bool check_obstacles, trajectory_t& traj, optimization_data_t* opt_data)
    {
      int res = dynamical_system->extend_to(si, sf, traj, opt_data);
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

    virtual cost_t* evaluate_extend_cost(const state* si, const state* sf, optimization_data_t*& opt_data)
    {
      float total_variation = dynamical_system->evaluate_extend_cost(si, sf, opt_data);
      return new cost_t(total_variation);
    }

    virtual float get_state_cost(const state& s)
    {
      return obstacle_map->get_state_cost(s.x);
    };

    virtual cost_t evaluate_trajectory_cost(trajectory_t& traj)
    {
      float c=0;
      for(auto& ps : traj.states)
        c += get_state_cost(state(ps));
      return cost_t(c);
    }
    
    virtual cost_t* get_zero_cost(){
      return new cost_t(0);
    }

    virtual cost_t* get_infinite_cost(){
      return new cost_t(FLT_MAX);
    }

    void test_extend_to()
    {
      trajectory_t traj;
      float zero[3] ={0};
      state origin(zero);
      for(int i=0; i<10; i++)
      {
        state sr;
        sample_in_goal(sr);
        cout<<"sampled: "; sr.print();
        extend_to(&origin, &sr, true, traj);
        cout<<"cost: "<< traj.total_variation<<endl;
        traj.print();
        getchar();
      }
      traj.clear();
    }
};

#endif
