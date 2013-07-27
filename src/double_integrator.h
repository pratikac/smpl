#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"
#include "double_integrator_mathematica.h"

class double_integrator_opt_data_c : public optimization_data_c
{
  public:
    float u111, u121, u211, u221;
    float t11, t21, t221, t222, t121, t122; 
    float T1, T2, T;
    bool is_initialized;

    double_integrator_opt_data_c() : is_initialized(false) {}
    ~double_integrator_opt_data_c(){}
};

class double_integrator_c : public dynamical_system_c<state_c<4>, control_c<2>, double_integrator_opt_data_c>
{
  public:
    typedef dynamical_system_c<state_c<4>, control_c<2>, optimization_data_c> dynamical_system_t;
    typedef typename dynamical_system_t::state_t state_t;
    typedef typename dynamical_system_t::control_t control_t;
    typedef typename dynamical_system_t::trajectory_t trajectory_t;
    typedef double_integrator_opt_data_c double_integrator_opt_data_t;
    
    double_integrator_c() {}

    int extend_to(const state_t& si, const state_t& sf,
        trajectory_t& traj, double_integrator_opt_data_t& opt_data)
    {
      if(!opt_data.is_initialized)
        if(evaluate_extend_cost(si, sf, opt_data) < 0)
          return 1;
      
      return 0;
    }
    
    float evaluate_extend_cost(const state_t& si, const state_t& sf,
        double_integrator_opt_data_t& opt_data)
    {
      if(opt_data.is_initialized)
        return opt_data.T;

      float um=1;
      float up=-um;
      float u111, u121;
      float u211, u221;
      float x1f = sf[0] - si[0];
      float dx1f = sf[2]-si[2];
      float x2f = sf[1] - si[1];
      float dx2f = sf[3] - si[3];
      
      if(dx1f > 0){
        u111 = um;
        u121 = um;
      }
      else{
        u111 = up;
        u121 = up;
      }
      if(dx2f > 0){
        u211 = um;
        u221 = um;
      }
      else{
        u211 = up;
        u221 = up;
      }

      float T = max(mT1, mT2);
      //cout<<"T: "<< T << endl;
      if(T != T)
      {
        //cout<<"returned"<<endl;
        return -1;
      }
      opt_data.u111 = u111;
      opt_data.u121 = u121;
      opt_data.u211 = u211;
      opt_data.u221 = u221;

      opt_data.t11 = mt11;
      opt_data.t21 = mt21;
      
      opt_data.t121 = mt121;
      opt_data.t122 = mt122;
      opt_data.t221 = mt221;
      opt_data.t222 = mt222;

      opt_data.T1 = mT1;
      opt_data.T2 = mT2;
      opt_data.T = max(opt_data.T1, opt_data.T2);
      opt_data.is_initialized = true;
      //cout<<"T: "<< opt_data.T << endl;
      return  opt_data.T;
    }
}; 

#endif
