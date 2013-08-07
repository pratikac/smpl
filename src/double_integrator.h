#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"
#include <cassert>

class double_integrator_optimization_data_c : public optimization_data_c
{
  public:
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
      
      float x10 = si[0] - sf[0];
      float x20 = si[1] - sf[1];
      float dx10 = si[2]- sf[2];
      float dx20 = si[3] - sf[3];
      float t1[2] = {x10, dx10};
      float t2[2] = {x20, dx20};
      
      float T = opt_data.T, T1 = opt_data.T1, T2 = opt_data.T2;
      float t11, t12, t21, t22;

      float um = 1, g = 1;
      float u1, u2;
      
      get_time(t1, um, t11, t12);
      get_time(t2, um, t21, t22);
      if(is_right(t1[0], t1[1], um))
        u1 = -um;
      else
        u1 = um;
      if(is_right(t2[0], t2[1], um))
        u2 = -um;
      else
        u2 = um;

      // dim 2 needs to slow down
      if( fabs(T-T1) < 1e-6)
      {
        g = get_gain(t2, T, um);
        get_time(t2, g*um, t21, t22);
          
        if(is_right(t2[0], t2[1], g*um))
          u2 = -g*um;
        else
          u2 = g*um;

        if(g < 0)
          return -1;
      }
      // dim 1 needs to slow down
      else
      {
        g = get_gain(t1, T, um);
        get_time(t1, g*um, t11, t12); 
        
        if(is_right(t1[0], t1[1], g*um))
          u1 = -g*um;
        else
          u1 = g*um;
        
        if(g < 0)
          return -1;
      }
#if 1 
      float t = 0, dt = 0.005;
      state_t sc;
      sc.x[0] = x10; sc.x[1] = x20;
      sc.x[2] = dx10; sc.x[3] = dx20;
      control_t cc;
      while(t < T)
      {
        if(t < t11)
          cc.x[0] = u1;
        else
          cc.x[0] = -u1;
        if(t < t21)
          cc.x[1] = u2;
        else
          cc.x[1] = -u2;

        sc.x[2] += cc[0]*dt;
        sc.x[3] += cc[1]*dt;
        sc.x[0] += sc[2]*dt;
        sc.x[1] += sc[3]*dt;

        traj.states.push_back(sc+sf);
        traj.controls.push_back(cc);
        t += dt;
      }
#else
      traj.states.push_back(si);
      traj.states.push_back(sf);
#endif
      return 0;
    }
    
    /*
     * use control g*um to reach origin from x0[2] in T time units
     */
    float get_f(float x10, float dx10, float g, float um, float T)
    {
      return get_t1(x10, dx10, g*um) + get_t2(x10, dx10, g*um) - T;        
    }

    float get_gain(float x0[2], float T, float um)
    {
      float x10 = x0[0];
      float dx10 = x0[1];
      
      float g=0.5, gm=1e-10, gp=1;
      float f, fm, fp;
      fm = get_f(x10, dx10, gm, um, T);
      fp = get_f(x10, dx10, gp, um, T);
      if( (fm < -FLT_MAX/2) || (fp < -FLT_MAX/2) || (fm*fp > 0))
        return -1;

      bool is_converged = false;
      while(!is_converged)
      {
        g = (gm+gp)/2.0;
        f = get_f(x10, dx10, g, um, T);
        if(f < -FLT_MAX/2)
          return -1;
        
        if(f*fm < 0)
          gp = g;
        else if(f*fp < 0)
          gm = g;

        is_converged = (gp-gm) < 0.1;
      }
      return g;
    }

    bool is_right(float x10, float dx10, float um)
    {
      return (((x10 > dx10*dx10/2/um) || (dx10 > 0)) && (x10 > -dx10*dx10/2/um));
    }
    
    float get_t1(float x10, float dx10, float um)
    {
      if(is_right(x10, dx10, um))
      {
        float t1 = 2*dx10*um;
        float d1 = sqrt(2)*sqrt(dx10*dx10*um*um + 2*pow(um,3)*x10);
        float t1r1 = (t1 - d1)/2/um/um;
        float t1r2 = (t1 + d1)/2/um/um;
        assert(max(t1r1, t1r2) > 0);
        return max(t1r1, t1r2);
      }
      else
      {
        float t1 = -2*dx10*um;
        float d1 = sqrt(2)*sqrt(dx10*dx10*um*um - 2*pow(um,3)*x10);
        float t1l1 = (t1 - d1)/2/um/um;
        float t1l2 = (t1 + d1)/2/um/um;
        assert(max(t1l1, t1l2) > 0);
        return max(t1l1, t1l2);
      }
    }

    float get_t2(float x10, float dx10, float um)
    {
      if(is_right(x10, dx10, um))
        return sqrt(um*um*(dx10*dx10 + 2*um*x10)/2)/um/um;
      else
        return sqrt(-um*um*(-dx10*dx10 + 2*um*x10)/2)/um/um;
    }
    
    float get_time(float x0[2], float um, float& t1, float& t2)
    {
      float x10 = x0[0];
      float dx10 = x0[1];
      
      t1 = get_t1(x10, dx10, um);
      t2 = get_t2(x10, dx10, um);

      return t1+t2;
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

      float T1 = get_time(t1, um, t11, t12);
      float T2 = get_time(t2, um, t21, t22);
      float T = max(T1, T2);
      //cout<<"T1: "<< T1 << " T2: "<< T2 << " T: "<< T << endl;
      
      opt_data.T1 = T1;
      opt_data.T2 = T2;
      opt_data.T = T;

      opt_data.is_initialized = true;
      
      //cout<<"T: "<< T << endl;
      return T;
    }

    void test_extend_to()
    {
      trajectory_t traj;
      float zero[4] = {0};
      state_t origin(zero);
      float goal[4] = {10, 5, 2, 1};
      state_t sr(goal);
      sr.print(cout, "sampled:","\n");

      double_integrator_optimization_data_c opt_data;
      if(extend_to(origin, sr, traj, opt_data))
        cout<<"could not connect"<<endl;
      cout<<"cost: "<< traj.total_variation<<endl;
      traj.print();
    }
}; 

#endif
