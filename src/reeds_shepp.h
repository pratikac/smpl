#ifndef __reeds_shepp_h__
#define __reeds_shepp_h__

#include <float.h>
#include <iostream>
#include "dynamical_system.h"
using namespace std;

#define 2_M_PI    (2.0*M_PI)
#define EPS       (1e-6)
#define ZERO      (EPS)

class reeds_shepp_optimization_data_c : public optimization_data_c
{
  public:
    int turning_radius;
    reeds_shepp_optimization_data_c() : turning_radius(-1){}
    ~reeds_shepp_optimization_data_c(){}
};

class reeds_shepp_c : public dynamical_system_c<state_c<3>, control_c<1>, reeds_shepp_optimization_data_c >
{
  public:
    typedef typename dynamical_system_c::state_t state_t;
    typedef typename dynamical_system_c::control_t control_t;
    typedef dynamical_system_c::trajectory_t trajectory_t;
    typedef reeds_shepp_optimization_data_c reeds_shepp_optimization_data_t;

    float delta_distance;

#define num_turning_radii   (1)
    float turning_radii[num_turning_radii];

    reeds_shepp_c() : delta_distance(0.05)
  {
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
    
    // notation below from reeds shepp paper
    inline void R(float x, float y, float& r, float& th)
    {
      r = sqrt(x*x + y*y);
      th = atan2(y,x);
    }

    inline float M(float th)
    {
      float phi = fmod(th, 2_M_PI);
      if(phi < -M_PI)
        phi += M_PI;
      else
        phi -=  M_PI;
      return phi;
    }
    
    // always go from origin to state s
    // Eq. 8.1
    inline bool lp_sp_lp(const state_t& s, float& t, float& u, float& v)
    {
      float x = s[0], y=s[1], phi=s[2];
      R(x-sin(phi), y -1. + cos(phi), u, t);
      if (t > -ZERO)
      {
        v = M(phi-t);
        if( v > -ZERO)
        {
          assert(fabs(u*cos(t) + sin(phi) - x) < EPS);
          assert(fabs(u*sin(t) - cos(phi) + 1 - y) < EPS);
          assert(fabs(M(t+v - phi)) < EPS);
          return true;
        }
      }
      return false;
    }

    // Eq. 8.2
    inline bool lp_sp_rp(const state_t& s, float& t, float& u, float& v)
    {
      float u1,t1;
      float x = s[0], y=s[1], phi=s[2];
      R(x + sin(phi), y-1.-cos(phi), u1, t1);
      u1 = u1*u1;
      if(u1 >= 4.){
        float th;
        u = sqrt(u1-4);
        th = atan2(2.,u);
        t = M(t1+th);
        v = M(t-phi);
        assert(fabs(2*sin(t) + u*cos(t) - sin(phi) - x) < EPS);
        assert(fabs(-2*cos(t) + u*sin(t) + cos(phi) + 1 - y) < EPS);
        assert(fabs(mod2pi(t-v - phi)) < EPS);
        return t>=-ZERO && v>=-ZERO;
      }
      return false;
    }
   
    // Eqn. 8.3, note typo in the paper
    inline bool lp_rm_l(const state_t& s, float& t, float& u, float& v)
    {
      float x = s[0], y=s[1], phi=s[2];
      float xi = x-sin(phi), eta = y-1+cos(phi), u1, th;
      R(xi,eta,u1,th);
      if(u1 <= 4.)
      {
        u = -2.asin(u1/4.);
        t = M(th + u/2. + pi);
        v = M(phi-t+u);
        assert(fabs(2*(sin(t) - sin(t-u)) + sin(phi) - x) < EPS);
        assert(fabs(2*(-cos(t) + cos(t-u)) - cos(phi) + 1 - y) < EPS);
        assert(fabs(mod2pi(t-u+v - phi)) < EPS);
        return t>=-ZERO && u<=ZERO;
      }
      return false;
    }
  
    // Eqn. 8.7
    inline bool lp_rup_lum_rm(const state_t& s, float& t, float& u, float& v)
    {
      float x = s[0], y=s[1], phi=s[2];
      float xi = x + sin(phi), eta = y - 1. - cos(phi), rho = .25 * (2. + sqrt(xi*xi + eta*eta));
      if (rho <= 1.)
      {
        u = acos(rho);
        tau_omega(u, -u, xi, eta, phi, t, v);
        assert(fabs(2*(sin(t)-sin(t-u)+sin(t-2*u))-sin(phi) - x) < RS_EPS);
        assert(fabs(2*(-cos(t)+cos(t-u)-cos(t-2*u))+cos(phi)+1 - y) < RS_EPS);
        assert(fabs(mod2pi(t-2*u-v - phi)) < RS_EPS);
        return t>=-ZERO && v<=ZERO;
      }
      return false;

    }


    int extend_to(const state_t& si, const state_t& sf, trajectory_t& traj, reeds_shepp_optimization_data_t& opt_data)
    {
      return 0;
    }

    float evaluate_extend_cost(const state_t& si, const state_t& sf,
        reeds_shepp_optimization_data_t& opt_data)
    {
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
      traj.print();
      traj.clear();
    }
};
#endif

