#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"
#include "double_integrator_mathematica.h"
#include <cassert>

class double_integrator_optimization_data_c : public optimization_data_c
{
  public:
    float t11, t12, t21, t22;
    float T1, T2, T;
    float u11, u21;
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
      
      float x10 = si[0] - sf[0];
      float dx10 = si[2]-sf[2];
      float x20 = si[1] - sf[1];
      float dx20 = si[3] - sf[3];
      float t1[2] = {x10, dx10};
      float t2[2] = {x20, dx20};
      
      float T = opt_data.T, T1 = opt_data.T1, T2 = opt_data.T2;
      float t11 = opt_data.t11, t12 = opt_data.t12, t21 = opt_data.t21, t22 = opt_data.t22;
      float u11 = opt_data.u11, u21 = opt_data.u21;
      
      float um = 1, g = 1;
      // dim 2 needs to slow down
      if( fabs(T-T1) < fabs(T-T2))
      {
        g = get_gain(t2, T, um, u21);
      }
      // dim 1 needs to slow down
      else
      {
        g = get_gain(t1, T, um, u11);
      }
      return 0;
    }
    
    /*
     * use control g*um to reach origin from x0[2] in T time units
     * calculates g using newton-raphson
     * u1 is used to check the direction left /right
     */
    void get_f_df(float x10, float dx10, float g, float um, float T, bool left, bool t1one, float& f, float& df)
    {
      if(left)
      {
        if(t1one){
          f = mt1l1(g) + mt2l(g) -T;
          df = mdt1l1(g) + mdt2l(g);
        }
        else{
          f = mt1l2(g) + mt2l(g) - T;
          df = mdt1l2(g) + mdt2l(g);
        }
      }
      else
      {
        float g1 = -g;
        if(t1one){
          f = mt1l1(g1) + mt2l(g1) -T;
          df = mdt1l1(g1) + mdt2l(g1);
        }
        else{
          f = mt1l2(g1) + mt2l(g1) -T;
          df = mdt1l2(g1) + mdt2l(g1);
        }
      }
    }

    float get_gain(float x0[2], float T, float um, float u1)
    {
      float x10 = x0[0];
      float dx10 = x0[1];
      bool left, t1one;
      if(u1 > 0)
      {
        left = true;
        if(mt1l1(1) > 0)
          t1one = true;
        else if(mt1l2(1) > 0)
          t1one = false;
      }
      else
      {
        left = false;
        if(mt1l1(1) > 0)
          t1one = true;
        else if(mt1l2(1) > 0)
          t1one = false;
      }
     
#if 0
      float g = 0.8;
      float gp = g;
      bool is_converged = false;
      while(!is_converged)
      {
        float f, df;
        get_f_df(x10, dx10, g, um, T, left, t1one, f, df);
        g = g - f/df;
        is_converged = fabs(g - gp) < 0.01;
        gp = g;

        cout<<"g: "<< g << endl;
        getchar();
      }
#else
      float g=0.5, gm=1e-10, gp=1;
      float f, fm, fp, df;
      get_f_df(x10, dx10, gm, um, T, left, t1one, fm, df);
      get_f_df(x10, dx10, gp, um, T, left, t1one, fp, df);
      assert(fm*fp < 0);

      bool is_converged = false;
      while(!is_converged)
      {
        g = (gm+gp)/2;
        get_f_df(x10, dx10, g, um, T, left, t1one, f, df);
        
        if(f*fm < 0)
          gp = g;
        else if(f*fp < 0)
          gm = g;

        is_converged = (gp-gm) < 0.05;
        cout<<"g: "<< g << endl;
        getchar();
      }
#endif
      return g;
    }

    /*
     * t1 = time till first switch
     * t2 = T - t1, i.e., time until reaches origin
     */
    float get_time(float x0[2], float um, float& t1, float& t2, float& u1)
    {
      float x10 = x0[0];
      float dx10 = x0[1];
      if(is_right(x10, dx10, 1))
      {
        u1 = -um;
        if(mt1l1(-1) > 0){
          t1 = mt1l1(-1);
        }
        else if(mt1l2(-1) > 0){
          t1 = mt1l2(-1);
        }
        else
        {

          cout<<"both t1 < 0"<<endl;
          getchar();
        }
        t2 = mt2l(-1);
        return t1+t2;
      }
      else if(is_left(x10, dx10, 1))
      {
        u1 = um;
        if(mt1l1(1) > 0){
          t1 = mt1l1(1);
        }
        else if(mt1l2(1) > 0){
          t1 = mt1l2(1);
        }
        else
        {

          cout<<"both t1 < 0"<<endl;
          getchar();
        }
        t2 = mt2l(1);
        return t1+t2;
      }
      else
      {
        cout<<"could not find region in state-space"<<endl;
        getchar();
      }
      return -1;
    }
    
    float evaluate_extend_cost(const state_t& si, const state_t& sf,
        double_integrator_opt_data_t& opt_data)
    {
      if(opt_data.is_initialized)
        return opt_data.T;

      float um=1;

      float x10 = si[0] - sf[0];
      float dx10 = si[2]-sf[2];
      float x20 = si[1] - sf[1];
      float dx20 = si[3] - sf[3];
      
      float t1[2] = {x10, dx10};
      float t2[2] = {x20, dx20};
      float t11, t12, t21, t22;
      float u11, u21;

      float T1 = get_time(t1, um, t11, t12, u11);
      float T2 = get_time(t2, um, t21, t22, u21);
      float T = max(T1, T2);

      opt_data.T = T;
      opt_data.T1 = T1;
      opt_data.T2 = T2;

      opt_data.t11 = t11;
      opt_data.t12 = t12;
      opt_data.t21 = t21;
      opt_data.t22 = t22;

      opt_data.u11 = u11;
      opt_data.u21 = u21;

      opt_data.is_initialized = true;
      
      cout<<"T: "<< T << endl;
      return T;
    }

    void test_extend_to()
    {
      trajectory_t traj;
      float zero[4] ={0};
      state_t origin(zero);
      float goal[4] = {10, 5, 2, 1};
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
