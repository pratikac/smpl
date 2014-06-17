#ifndef __reeds_shepp_h__
#define __reeds_shepp_h__

#include <float.h>
#include <iostream>
#include "dynamical_system.h"
#include "rs.h"
using namespace std;

class reeds_shepp_optimization_data_c : public optimization_data_c
{
  public:
    int _turning_radius;
    int _numero;
    double _t, _u, _v;
    
    reeds_shepp_optimization_data_c(){
      _turning_radius = -1;
    }
    reeds_shepp_optimization_data_c(int tr, int numero, double t, double u, double v)
    {
      _turning_radius = tr;
      _numero = numero;
      _t = t;
      _u = u;
      _v = v;
    }
    ~reeds_shepp_optimization_data_c(){}
    void print(ostream& os)
    {
      os<<"tr: "<< _turning_radius<< endl;
      os <<"numero: "<< _numero << endl;
      cout<< "t: "<< _t << " u: "<< _u << " v: " << _v << endl;
    }

};

/*
 * state:   (x,y,theta)
 * control: (numero of pattern, t, u, v)
 */
class reeds_shepp_c : public dynamical_system_c<state_c<3>, control_c<4>, reeds_shepp_optimization_data_c >
{
  public:
    typedef typename dynamical_system_c::state_t state_t;
    typedef typename dynamical_system_c::control_t control_t;
    typedef dynamical_system_c::trajectory_t trajectory_t;
    typedef reeds_shepp_optimization_data_c reeds_shepp_optimization_data_t;

    double delta_distance;

#define num_turning_radii   (1)
    float turning_radii[num_turning_radii];
    
    reeds_shepp_c()
    {
      delta_distance = 0.05;
      turning_radii[0] = 1;
      //turning_radii[1] = 6;
      //turning_radii[2] = 8;
    };

    int get_plotter_state(const state_t& s, float* ps)
    {
      ps[0] = s[0];
      ps[1] = s[1];
      ps[2] = 0;
      return 0;
    }

    int sample_state(float* center, float* size, float* s)
    {
      for(int i : range(0,3))
        s[i] = center[i] + (RANDF-0.5)*size[i];
      return 0;
    }

    void get_trajectory_from_numero(const state_t& si, const state_t& sf,
        int numero, double t, double u, double v, trajectory_t& traj)
    {
      int N = 10000;
      double px[N], py[N], pth[N];
      int traj_len = constRS(numero, t, u, v, si[0], si[1], si[2],
        delta_distance,
        px, py, pth);
      
      for(int i=0; i< traj_len; i++)
      {
        float s[3] = {(float)px[i], (float)py[i], (float)pth[i]};
        traj.states.push_back(state_t(s));
      }
      traj.total_variation = traj_len;
      return;
    }

    int extend_to(const state_t& si, const state_t& sf, trajectory_t& traj, reeds_shepp_optimization_data_t& opt_data)
    {
      if(opt_data._turning_radius < 0)
      {
        evaluate_extend_cost(si, sf, opt_data);
      }
      get_trajectory_from_numero(si, sf,
          opt_data._numero, opt_data._t, opt_data._u, opt_data._v,
          traj);
      return 0;
    }

    float evaluate_extend_cost(const state_t& si, const state_t& sf,
        reeds_shepp_optimization_data_t& opt_data)
    {
      double rho = 1;     // cost = len + rho/turning_radius
      double min_cost = FLT_MAX;

      for(int i=num_turning_radii-1; i >=0; i--)
      {
        float tr = turning_radii[i];
        int numero;
        double t, u, v;
        float len = min_length_rs(si[0], si[1], si[2],
            sf[0], sf[1], sf[2],
            &numero, &t, &u, &v);

        float cost = len + rho*1./tr; 
        if(cost < min_cost)
        {
          min_cost = cost;
          opt_data = reeds_shepp_optimization_data_t(tr, numero, t, u, v); 
        }

      }
      if((min_cost < 0) || (min_cost > FLT_MAX/2.))
        return -1;
      return min_cost;
    }


    void test_extend_to()
    {
      trajectory_t traj;
      float zero[3] ={0};
      state_t origin(zero);
      float goal[3] = {0,0, 175/180.0*M_PI};
      state_t sr(goal);
      sr.print(cout, "sampled:","\n");

      reeds_shepp_optimization_data_c opt_data;
      extend_to(origin, sr, traj, opt_data);
      cout<<"cost: "<< traj.total_variation<<endl;
      //traj.print();
      traj.clear();
      
      /*
      float zero[3] ={0};
      state_t origin(zero);
      float goal[3] = {0, -1, 1/180.*M_PI};
      state_t sr(goal);
      
      reeds_shepp_optimization_data_c opt_data;
      cout<<"cost: "<< evaluate_extend_cost(origin, sr, opt_data) << endl;
      opt_data.print(cout);
      */
    }
};
#endif

