#ifndef __SINGLE_INTEGRATOR_H__
#define __SINGLE_INTEGRATOR_H__

#include <float.h>
#include <iostream>
#include "dynamical_system.h"

using namespace std;

template<size_t N>
class single_integrator_t : public dynamical_system_t<N>
{
  typedef state_t<N> state;

  int normalize_diff(const state* s1, const state* s2, float factor, float* vout)
  {
    for(size_t i=0; i<N; i++)
      vout[i] = (s2->x[i]-s1->x[i])*factor;
    return 0;
  }

  int add(const float* s1, const float* s2, float* res)
  {
    for(size_t i=0; i<N; i++)
      res[i] = s2[i]+s1[i];
    return 0;
  }
  public:

  float delta_distance;

  single_integrator_t() : delta_distance(0.5)
  {
  };

  int extend_to(const state* si, const state* sf, trajectory_t& traj, optimization_data_t* opt_data=NULL)
  {
    traj.clear();
    float dist = evaluate_extend_cost(si, sf);
    traj.total_variation = dist;
    int num_points = ceil(dist/delta_distance);

    float step[N];
    if(normalize_diff(si, sf, float(1.0/num_points), step)!=0)
      return 1;

    float* init_state = new float[N];
    for(size_t i=0; i<N; i++) 
      init_state[i] = si->x[i];

    traj.states.push_back(init_state);

    for(int i=1; i<num_points; i++){
      float* cur_state = new float[N];
      if(add(traj.states.back(), step, cur_state)) 
        return 1;
      traj.states.push_back(cur_state);
    }

    return 0;
  }

  float evaluate_extend_cost(const state* si, const state* sf, optimization_data_t* opt_data=NULL)
  {
    return si->dist(*sf);
  }

};

#endif
