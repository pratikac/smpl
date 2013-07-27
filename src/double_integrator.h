#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"

class double_integrator_opt_data_c : public optimization_data_c
{
  public:
    double_integrator_opt_data_c() {}
    ~double_integrator_opt_data_c(){}
};

class double_integrator_c : public dynamical_system_c<state_c<4>, control_c<2>, optimization_data_c>
{
  public:
    typedef dynamical_system_c<state_c<4>, control_c<2>, optimization_data_c> dynamical_system_t;
    typedef typename dynamical_system_t::state_t state_t;
    typedef typename dynamical_system_t::control_t control_t;
    typedef typename dynamical_system_t::trajectory_t trajectory_t;
    typedef double_integrator_opt_data_c double_integrator_opt_data_t;

    int extend_to(const state_t& si, const state_t& sf,
        trajectory_t& traj, double_integrator_opt_data_t& opt_data)
    {
      return 0;
    }
    
    float evaluate_extend_cost(const state_t& si, const state_t& sf,
        double_integrator_opt_data_t& opt_data)
    {
      float u1m,u1p,u2m, u2p;
      float vx1f = sf[0] - si[0];
      float vdx1f = sf[2]-si[2];
      float vx2f = sf[1] - si[1];
      float vdx2f = sf[3] - si[3];
      
      if(vdx1f > 0){
        u1m = 1;
        u1p = -1;
      }
      else{
        u1m = -1;
        u1p = 1;
      }
      if(vdx2f > 0){
        u2m = 1;
        u2p = -1;
      }
      else{
        u2m = -1;
        u2p = 1;
      }

      float T1 = (vdx1f - (u1m*sqrt(power(vdx1f,2) - 2*u1p*vx1f))/sqrt(u1m*(u1m - u1p)) + 
          (u1p*sqrt(power(vdx1f,2) - 2*u1p*vx1f))/sqrt(u1m*(u1m - u1p)))/u1p;
    
      float T2 = (vdx2f - (u2m*sqrt(power(vdx2f,2) - 2*u2p*vx2f))/sqrt(u2m*(u2m - u2p)) + 
            (u2p*sqrt(power(vdx2f,2) - 2*u2p*vx2f))/sqrt(u2m*(u2m - u2p)))/u2p;
      
      float T = max(T1, T2);
      if(T == T1)
      {

      }
      else
      {

      }
    }
}; 

#endif
