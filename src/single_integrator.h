#ifndef __single_integrator_h__
#define __single_integrator_h__

#include <float.h>
#include <iostream>
#include "dynamical_system.h"
using namespace std;

template<size_t N>
class single_integrator_c : public dynamical_system_c<state_c<N>, control_c<N>, optimization_data_c>
{
  public:
    typedef dynamical_system_c<state_c<N>, control_c<N>, optimization_data_c> dynamical_system_t;
    typedef typename dynamical_system_t::state_t state_t;
    typedef typename dynamical_system_t::control_t control_t;
    typedef typename dynamical_system_t::trajectory_t trajectory_t;
    typedef optimization_data_c single_integrator_opt_data_t;

    int normalize_diff(const state_t& s1, const state_t& s2, float factor, float* vout)
    {
      for(size_t i=0; i<N; i++)
        vout[i] = (s2.x[i]-s1.x[i])*factor;
      return 0;
    }


    float delta_distance;

    single_integrator_c() : delta_distance(0.05){};

    int get_plotter_state(const state_t& s, float* ps)
    {
      int m = min(3,(int)N);
      for(int i=0; i<m; i++)
        ps[i] = s[i];
      return 0;
    }

    int sample_state(float* center, float* size, float* s)
    {
      for(int i : range(0,N))
        s[i] = center[i] + (RANDF-0.5)*size[i];
      return 0;
    }

    int extend_to(const state_t& si, const state_t& sf,
        trajectory_t& traj, single_integrator_opt_data_t& opt_data)
    {
      traj.clear();
      float dist = evaluate_extend_cost(si, sf, opt_data);
      traj.total_variation = dist;
      traj.dt = delta_distance;     // assume velocity of 1 m/s
      int num_points = ceil(dist/delta_distance);

      float step[N];
      if(normalize_diff(si, sf, float(1.0/num_points), step)!=0)
        return 1;
  
      traj.states.reserve(num_points+1);
      traj.states.push_back(si);
      for(int i=1; i<num_points; i++)
      {
        state_t cur_state = traj.states.back() + state_t(step);
        traj.states.push_back(cur_state);
      }
      return 0;
    }

    float evaluate_extend_cost(const state_t& si, const state_t& sf,
        single_integrator_opt_data_t& opt_data)
    {
      return si.dist(sf);
    }
    
    void test_extend_to()
    {
      trajectory_t traj;
      float zero[N] ={0};
      state_t origin(zero);
      float goal[N] = {10};
      state_t sr(goal);
      sr.print(cout, "sampled:","\n");

      optimization_data_c opt_data;
      extend_to(origin, sr, traj, opt_data);
      cout<<"cost: "<< traj.total_variation<<endl;
      traj.print();
      traj.clear();
    }
};

#endif
