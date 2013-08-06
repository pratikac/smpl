#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"
#include "double_integrator_mathematica.h"

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

      return 0;
    }
    
    /*
     * t1 = time till first switch
     * t2 = T - t1, i.e., time until reaches origin
     */
    float get_time(float x0[2], float um, float& t1, float& t2)
    {
      float x10 = x0[0];
      float dx10 = x0[1];
      if( ((x10 > dx10*dx10/2/um) || (dx10 > 0)) && ((x10 > -dx10*dx10/2/um) || (x10 > 0)) )
      {
        if(t1r1 > 0){
          t1 = t1r1;
        }
        else if(t1r2 > 0){
          t1 = t1r2;
        }
        else
        {

          cout<<"both t1 < 0"<<endl;
          getchar();
        }
        t2 = t2r;
        return t1+t2;
      }
      else if( ((x10 < -dx10*dx10/2/um) || (dx10 < 0)) && ((x10 < dx10*dx10/2/um) || (x10 < 0)) )
      {
        if(t1l1 > 0){
          t1 = t1l1;
        }
        else if(t1l2 > 0){
          t1 = t1l2;
        }
        else
        {

          cout<<"both t1 < 0"<<endl;
          getchar();
        }
        t2 = t2l;
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

      float T1 = get_time(t1, um, t11, t12);
      float T2 = get_time(t2, um, t21, t22);
      float T = max(T1, T2);

      opt_data.T = T;
      opt_data.T1 = T1;
      opt_data.T2 = T2;
      
      //cout<<"T: "<< T << endl;
      return T;
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
