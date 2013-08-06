#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"
#include "double_integrator_mathematica.h"

class double_integrator_optimization_data_c : public optimization_data_c
{
  public:
    float u111, u121, u211, u221;
    float t11, t21, t221, t222, t121, t122; 
    float T1, T2, T;
    bool is_initialized;

    double_integrator_optimization_data_c() : is_initialized(false) {}
    ~double_integrator_optimization_data_c(){}
};

class double_integrator_c : public dynamical_system_c<state_c<4>, control_c<2>, double_integrator_optimization_data_c>
{
  public:
    typedef dynamical_system_c<state_c<4>, control_c<2>, optimization_data_c> dynamical_system_t;
    typedef typename dynamical_system_t::state_t state_t;
    typedef typename dynamical_system_t::control_t control_t;
    typedef typename dynamical_system_t::trajectory_t trajectory_t;
    typedef double_integrator_optimization_data_c double_integrator_opt_data_t;
    
    double_integrator_c() {}
    
    int get_plotter_state(const state_t& s, float* ps)
    {
      ps[0] = s[0];
      ps[1] = s[1];
      ps[2] = 0;
      return 0;
    }

    int extend_to(const state_t& si, const state_t& sf,
        trajectory_t& traj, double_integrator_opt_data_t& opt_data)
    {
      if(!opt_data.is_initialized)
        if(evaluate_extend_cost(si, sf, opt_data) < 0)
          return 1;
      
      traj.clear();
      traj.total_variation = opt_data.T;

      float u111 = opt_data.u111, u121 = opt_data.u121, u211 = opt_data.u211, u221 = opt_data.u221;
      float t11 = opt_data.t11, t21 = opt_data.t21, t221 = opt_data.t221, t222 = opt_data.t222, t121 = opt_data.t121, t122 = opt_data.t122; 
      
      float T = opt_data.T;
      float T1 = opt_data.T1;
      float T2 = opt_data.T2;
      float dt = 0.01;
      if(T == T1)
      {
        float t = 0;
        state_t sc = si;
        control_t cc;
        while(t < T)
        {
          if(t< t11){
            cc.x[0] = u111;
          }
          else{
            cc.x[0] = -u111;
          }
          sc.x[2] += cc[0]*dt;
          sc.x[0] += sc[2]*dt;

          if(t < t221){
            cc.x[1] = u221;
          }
          else if(t < t222){
            cc.x[1] = -u221;
          }
          else{
            cc.x[1] = u221;
          }
          sc.x[3] += cc[1]*dt;
          sc.x[1] += sc[3]*dt;
          
          traj.states.push_back(sc);
          traj.controls.push_back(cc);
          t += dt;
        }
      }
      else if(T == T2)
      {
        float t = 0;
        state_t sc = si;
        control_t cc;
        while(t < T)
        {
          if(t< t21){
            cc.x[1] = u211;
          }
          else{
            cc.x[1] = -u211;
          }
          sc.x[3] += cc[1]*dt;
          sc.x[1] += sc[3]*dt;

          if(t < t121){
            cc.x[0] = u121;
          }
          else if(t < t122){
            cc.x[0] = -u121;
          }
          else{
            cc.x[0] = u121;
          }
          sc.x[2] += cc[0]*dt;
          sc.x[0] += sc[2]*dt;
          
          traj.states.push_back(sc);
          traj.controls.push_back(cc);
          t += dt;
        }
      }
      else{
        cout<<"T is none of T1 and T2, return"<< endl;
        return 1;
      }
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
      
      float T1 = mT1;
      float T2 = mT2;
      if((T1 != T1) || (T2 != T2))
        return -1;
      float T = max(T1, T2);
      //cout<<"T: "<< T << endl;
      
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

      opt_data.T1 = T1;
      opt_data.T2 = T2;
      opt_data.T = T; 
      opt_data.is_initialized = true;
      //cout<<"T: "<< opt_data.T << endl;
      return  opt_data.T;
    }

    void test_extend_to()
    {
      trajectory_t traj;
      float zero[4] ={0};
      state_t origin(zero);
      float goal[4] = {-10, 5, 2, 1};
      state_t sr(goal);
      sr.print(cout, "sampled:","\n");

      double_integrator_optimization_data_c opt_data;
      if(extend_to(origin, sr, traj, opt_data))
        cout<<"could not connect"<<endl;
      cout<<"cost: "<< traj.total_variation<<endl;
      traj.print();
      traj.clear();
    }
}; 

#endif
