#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"
#include <cassert>

class double_integrator_optimization_data_c : public optimization_data_c
{
  public:
    float T1, T2, T;
    float t11, t21;
    float u1, u2;
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
    
    float umm;
    double_integrator_c() : umm(1) {}
    
    int get_plotter_state(const state_t& s, float* ps)
    {
      ps[0] = s[0];
      ps[1] = s[1];
      ps[2] = 0; //sqrt(s[2]*s[2]+s[3]*s[3]);
      return 0;
    }

    int sample_state(float* center, float* size, float* s)
    {
      for(int i : range(0,4))
        s[i] = center[i] + (RANDF-0.5)*size[i];
      return 0;
    }
    
    int extend_to(const state_t& si, const state_t& sf,
        trajectory_t& traj, double_integrator_opt_data_t& opt_data)
    {
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
      float u1, u2;
      float um = umm;

      if(!opt_data.is_initialized)
      {
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

        float g = 1.0;
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
            return 1;
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
            return 1;
        }
        
        opt_data.t11 = t11;
        opt_data.t21 = t21;
        opt_data.u1 = u1;
        opt_data.u2 = u2;
        opt_data.is_initialized = true;
      } 
      else
      {
        t11 = opt_data.t11;
        t21 = opt_data.t21;
        u1 = opt_data.u1;
        u2 = opt_data.u2;
      }
#if 1 
      float t = 0, dt = 0.01;
      traj.dt = dt;
      state_t sc;
      sc.x[0] = x10; sc.x[1] = x20;
      sc.x[2] = dx10; sc.x[3] = dx20;
      
      control_t cc;
      int num_steps = T/dt;
      traj.states.reserve(num_steps);
      traj.controls.reserve(num_steps);
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
      float t1, t2;
      float t[2] = {x10, dx10};
      return get_time(t, fabs(g*um), t1, t2) - T;        
    }

    float get_gain(float x0[2], float T, float um)
    {
      float x10 = x0[0];
      float dx10 = x0[1];
      
      float g=0.5, gm=1e-12, gp=1-1e-12;
      float f, fm, fp;
      fm = get_f(x10, dx10, gm, um, T);
      fp = get_f(x10, dx10, gp, um, T);
      if(fabs(fp)<1e-6)
      {
        return 1;
      }
      if((fm < -FLT_MAX/2) || (fp < -FLT_MAX/2) || (fm*fp > 0))
      {
        cout<<"get_gain: 184, ret -1"<<endl;
        return -1;
      }

      bool is_converged = false;
      int c = 0;
      while(!is_converged)
      {
        g = (gm+gp)/2.0;
        f = get_f(x10, dx10, g, um, T);
        if((f < -FLT_MAX/2) || (f != f))
        {
          cout<<"get_gain: 196, ret -1"<<endl;
          return -1;
        }

        if(fabs(f) < 1e-6)
          return g;

        if(f*fm < 0){
          gp = g;
          fp = get_f(x10, dx10, gp, um, T);
        }
        else if(f*fp < 0){
          gm = g;
          fm = get_f(x10, dx10, gm, um, T);
        }
        
        c++;
        if(c > 100)
          cout<<"g: "<< g << " fm: "<< fm <<" fp: "<< fp <<" f: "<< f << endl;
        is_converged = (gp-gm) < 1e-6;
      }
      return g;
    }

    bool is_right(float x10, float dx10, float um)
    {
      return (((x10 > dx10*dx10/2/um) || (dx10 > 0)) && (x10 > -dx10*dx10/2/um));
    }
    
    float get_t1(float x10, float dx10, float um)
    {
      if(is_right(x10, dx10, um)){
        float t1 = 2*dx10*um;
        float d1 = sqrt(2)*sqrt(dx10*dx10*um*um + 2*pow(um,3)*x10);
        return (t1+d1)/2/um/um;
      }
      else{
        float t1 = -2*dx10*um;
        float d1 = sqrt(2)*sqrt(dx10*dx10*um*um - 2*pow(um,3)*x10);
        return (t1+d1)/2/um/um;
      }
    }

    float get_t2(float x10, float dx10, float um)
    {
      if(is_right(x10, dx10, um))
        return sqrt(um*um*(dx10*dx10 + 2*um*x10)/2)/um/um;
      else
        return sqrt(um*um*(dx10*dx10 - 2*um*x10)/2)/um/um;
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

      float um=umm;

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

      //cout<<"T: "<< T << endl;
      return T;
    }

    void test_extend_to()
    {
      trajectory_t traj1, traj2;
      float zero[4] = {1, 2, -1, 0.1};
      state_t origin(zero);
      float goal[4] = {10, 5, -2, 1};
      state_t sr(goal);
      sr.print(cout, "sampled:","\n");

      double_integrator_optimization_data_c opt_data;
      if(extend_to(origin, sr, traj1, opt_data))
        cout<<"could not connect"<<endl;
      //traj1.print();

      float g2[4] = {12, 4, 1, -1};
      state_t sg2(g2);
      sg2.print(cout, "second goal:","\n");
      
      double_integrator_optimization_data_c opt_data2;
      if(extend_to(sr, sg2, traj2, opt_data2))
        cout<<"could not connect"<<endl;
      traj1.append(traj2);
      traj1.print();
    }
}; 

#endif
