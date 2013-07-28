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
      for(int i=0; i<3; i++)
        ps[i] = 0;
      if(N < 4){
        for(int i=0; i<N-1; i++)
          ps[i] = s[i];
      }
      else{
        cout<<"plotting single integrator for N > 3"<<endl;
      }
      return 0;
    }

    int extend_to(const state_t& si, const state_t& sf,
        trajectory_t& traj, single_integrator_opt_data_t& opt_data)
    {
      traj.clear();
      float dist = evaluate_extend_cost(si, sf, opt_data);
      traj.total_variation = dist;
      int num_points = ceil(dist/delta_distance);

      float step[N];
      if(normalize_diff(si, sf, float(1.0/num_points), step)!=0)
        return 1;

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

};

#endif
