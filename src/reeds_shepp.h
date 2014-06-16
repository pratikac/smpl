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
    
    enum reeds_shepp_segment_type{
      RS_NOP = 0,
      RS_LEFT = 1,
      RS_STRAIGHT = 2,
      RS_RIGHT = 3
    };
    static const reeds_shepp_segment_type reeds_shepp_path_type[18][5] = 
    {
      { RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP, RS_NOP },             // 0
      { RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP, RS_NOP },            // 1
      { RS_LEFT, RS_RIGHT, RS_LEFT, RS_RIGHT, RS_NOP },           // 2
      { RS_RIGHT, RS_LEFT, RS_RIGHT, RS_LEFT, RS_NOP },           // 3
      { RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP },        // 4
      { RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP },       // 5
      { RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP },        // 6
      { RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP },       // 7
      { RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP },       // 8
      { RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP },        // 9
      { RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_LEFT, RS_NOP },       // 10
      { RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_RIGHT, RS_NOP },        // 11
      { RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP },         // 12
      { RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP },         // 13
      { RS_LEFT, RS_STRAIGHT, RS_LEFT, RS_NOP, RS_NOP },          // 14
      { RS_RIGHT, RS_STRAIGHT, RS_RIGHT, RS_NOP, RS_NOP },        // 15
      { RS_LEFT, RS_RIGHT, RS_STRAIGHT, RS_LEFT, RS_RIGHT },      // 16
      { RS_RIGHT, RS_LEFT, RS_STRAIGHT, RS_RIGHT, RS_LEFT }       // 17
    };

    class reeds_shepp_path_c
    {
      public:
        float total_length;
        float segment_lengths[5];
        const reeds_shepp_segment_type* _type;

        reeds_shepp_path_c(const reeds_shepp_segment_type* type=reeds_shepp_path_type[0],
            float t=FLT_MAX/2., float u =0., float v=0.,
            float w=0., float x=0.);
    };
    
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

    inline void tau_omega(float u, float v, float xi, float eta, float phi,
        float& tau, float& omega)
    {
      float delta = M(u-v), A = sin(u) - sin(delta), B = cos(u) - cos(delta) - 1.;
      float t1 = atan2(eta*A - xi*B, xi*A + eta*B), t2 = 2. * (cos(delta) - cos(v) - cos(u)) + 3;
      tau = (t2<0) ? M(t1+pi) : M(t1);
      omega = M(tau - u + v - phi);
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
        assert(fabs(M(t-v - phi)) < EPS);
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
        assert(fabs(M(t-u+v - phi)) < EPS);
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
        assert(fabs(M(t-2*u-v - phi)) < RS_EPS);
        return t>=-ZERO && v<=ZERO;
      }
      return false;

    }

    // Eqn. 8.8
    inline bool lp_rum_lum_rp(const state_t& s, float& t, float& u, float& v)
    {
      float x = s[0], y=s[1], phi=s[2];
      float xi = x + sin(phi), eta = y - 1. - cos(phi), rho = (20. - xi*xi - eta*eta) / 16.;
      if (rho>=0 && rho<=1)
      {
        u = -acos(rho);
        if (u >= -.5 * pi)
        {
          tau_omega(u, u, xi, eta, phi, t, v);
          assert(fabs(4*sin(t)-2*sin(t-u)-sin(phi) - x) < RS_EPS);
          assert(fabs(-4*cos(t)+2*cos(t-u)+cos(phi)+1 - y) < RS_EPS);
          assert(fabs(M(t-v - phi)) < RS_EPS);
          return t>=-ZERO && v>=-ZERO;
        }
      }
      return false;
    }

    // Eqn. 8.9
    inline bool lp_rm_sm_lm(const state_t& s, float& t, float& u, float& v)
    {
      float x = s[0], y=s[1], phi=s[2];
      float xi = x - sin(phi), eta = y - 1. + cos(phi), rho, theta;
      R(xi, eta, rho, theta);
      if (rho >= 2.)
      {
        float r = sqrt(rho*rho - 4.);
        u = 2. - r;
        t = M(theta + atan2(r, -2.));
        v = M(phi - .5*pi - t);
        assert(fabs(2*(sin(t)-cos(t))-u*sin(t)+sin(phi) - x) < EPS);
        assert(fabs(-2*(sin(t)+cos(t))+u*cos(t)-cos(phi)+1 - y) < EPS);
        assert(fabs(M(t+pi/2+v-phi)) < EPS);
        return t>=-ZERO && u<=ZERO && v<=ZERO;
      }
      return false;
    }

    // Eqn. 8.10
    inline bool lp_rm_sm_rm(const state_t& s, float& t, float& u, float& v)
    {
      float x = s[0], y=s[1], phi=s[2];
      float xi = x + sin(phi), eta = y - 1. - cos(phi), rho, theta;
      R(-eta, xi, rho, theta);
      if (rho >= 2.)
      {
        t = theta;
        u = 2. - rho;
        v = M(t + .5*pi - phi);
        assert(fabs(2*sin(t)-cos(t-v)-u*sin(t) - x) < EPS);
        assert(fabs(-2*cos(t)-sin(t-v)+u*cos(t)+1 - y) < EPS);
        assert(fabs(M(t+pi/2-v-phi)) < EPS);
        return t>=-ZERO && u<=ZERO && v<=ZERO;
      }
      return false;
    }

    // Eqn. 8.11 *** TYPO IN PAPER ***
    inline bool lp_rm_slm_rp(const state_t& s, float& t, float& u, float& v)
    {
      float x = s[0], y=s[1], phi=s[2];
      float xi = x + sin(phi), eta = y - 1. - cos(phi), rho, theta;
      R(xi, eta, rho, theta);
      if (rho >= 2.)
      {
        u = 4. - sqrt(rho*rho - 4.);
        if (u <= ZERO)
        {
          t = M(atan2((4-u)*xi -2*eta, -2*xi + (u-4)*eta));
          v = M(t - phi);
          assert(fabs(4*sin(t)-2*cos(t)-u*sin(t)-sin(phi) - x) < EPS);
          assert(fabs(-4*cos(t)-2*sin(t)+u*cos(t)+cos(phi)+1 - y) < EPS);
          assert(fabs(M(t-v-phi)) < EPS);
          return t>=-ZERO && v>=-ZERO;
        }
      }
      return false;
    }   

    void CSC(float x, float y, float phi, reeds_shepp_path_c& path)
    {
        float t, u, v, Lmin = path.total_length, L;
        if (lp_sp_lp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[14], t, u, v);
            Lmin = L;
        }
        if (lp_sp_lp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[14], -t, -u, -v);
            Lmin = L;
        }
        if (lp_sp_lp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[15], t, u, v);
            Lmin = L;
        }
        if (lp_sp_lp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[15], -t, -u, -v);
            Lmin = L;
        }
        if (lp_sp_rp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[12], t, u, v);
            Lmin = L;
        }
        if (lp_sp_rp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[12], -t, -u, -v);
            Lmin = L;
        }
        if (lp_sp_rp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[13], t, u, v);
            Lmin = L;
        }
        if (lp_sp_rp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[13], -t, -u, -v);
    }
    
    void CCC(float x, float y, float phi, reeds_shepp_path_c& path)
    {
        float t, u, v, Lmin = path.length(), L;
        if (lp_rm_l(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[0], t, u, v);
            Lmin = L;
        }
        if (lp_rm_l(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[0], -t, -u, -v);
            Lmin = L;
        }
        if (lp_rm_l(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[1], t, u, v);
            Lmin = L;
        }
        if (lp_rm_l(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[1], -t, -u, -v);
            Lmin = L;
        }

        // backwards
        float xb = x*cos(phi) + y*sin(phi), yb = x*sin(phi) - y*cos(phi);
        if (lp_rm_l(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[0], v, u, t);
            Lmin = L;
        }
        if (lp_rm_l(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[0], -v, -u, -t);
            Lmin = L;
        }
        if (lp_rm_l(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[1], v, u, t);
            Lmin = L;
        }
        if (lp_rm_l(-xb, -yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[1], -v, -u, -t);
    }
    
    void CCCC(float x, float y, float phi, reeds_shepp_path_c& path)
    {
        float t, u, v, Lmin = path.length(), L;
        if (lp_rup_lum_rm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[2], t, u, -u, v);
            Lmin = L;
        }
        if (lp_rup_lum_rm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[2], -t, -u, u, -v);
            Lmin = L;
        }
        if (lp_rup_lum_rm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[3], t, u, -u, v);
            Lmin = L;
        }
        if (lp_rup_lum_rm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v))) // timeflip + reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[3], -t, -u, u, -v);
            Lmin = L;
        }

        if (lp_rum_lum_rp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[2], t, u, u, v);
            Lmin = L;
        }
        if (lp_rum_lum_rp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[2], -t, -u, -u, -v);
            Lmin = L;
        }
        if (lp_rum_lum_rp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[3], t, u, u, v);
            Lmin = L;
        }
        if (lp_rum_lum_rp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + 2.*fabs(u) + fabs(v))) // timeflip + reflect
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[3], -t, -u, -u, -v);
    }
    
    void CCSC(float x, float y, float phi, reeds_shepp_path_c& path)
    {
        float t, u, v, Lmin = path.length() - .5*pi, L;
        if (lp_rm_sm_lm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[4], t, -.5*pi, u, v);
            Lmin = L;
        }
        if (lp_rm_sm_lm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[4], -t, .5*pi, -u, -v);
            Lmin = L;
        }
        if (lp_rm_sm_lm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[5], t, -.5*pi, u, v);
            Lmin = L;
        }
        if (lp_rm_sm_lm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[5], -t, .5*pi, -u, -v);
            Lmin = L;
        }

        if (lp_rm_sm_rm(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[8], t, -.5*pi, u, v);
            Lmin = L;
        }
        if (lp_rm_sm_rm(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[8], -t, .5*pi, -u, -v);
            Lmin = L;
        }
        if (lp_rm_sm_rm(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[9], t, -.5*pi, u, v);
            Lmin = L;
        }
        if (lp_rm_sm_rm(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[9], -t, .5*pi, -u, -v);
            Lmin = L;
        }

        // backwards
        float xb = x*cos(phi) + y*sin(phi), yb = x*sin(phi) - y*cos(phi);
        if (lp_rm_sm_lm(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[6], v, u, -.5*pi, t);
            Lmin = L;
        }
        if (lp_rm_sm_lm(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[6], -v, -u, .5*pi, -t);
            Lmin = L;
        }
        if (lp_rm_sm_lm(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[7], v, u, -.5*pi, t);
            Lmin = L;
        }
        if (lp_rm_sm_lm(-xb, -yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[7], -v, -u, .5*pi, -t);
            Lmin = L;
        }

        if (lp_rm_sm_rm(xb, yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[10], v, u, -.5*pi, t);
            Lmin = L;
        }
        if (lp_rm_sm_rm(-xb, yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[10], -v, -u, .5*pi, -t);
            Lmin = L;
        }
        if (lp_rm_sm_rm(xb, -yb, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[11], v, u, -.5*pi, t);
            Lmin = L;
        }
        if (lp_rm_sm_rm(-xb, -yb, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[11], -v, -u, .5*pi, -t);
    }
    
    void CCSCC(float x, float y, float phi, reeds_shepp_path_c& path)
    {
        float t, u, v, Lmin = path.length() - pi, L;
        if (lp_rm_slm_rp(x, y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v)))
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[16], t, -.5*pi, u, -.5*pi, v);
            Lmin = L;
        }
        if (lp_rm_slm_rp(-x, y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[16], -t, .5*pi, -u, .5*pi, -v);
            Lmin = L;
        }
        if (lp_rm_slm_rp(x, -y, -phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // reflect
        {
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[17], t, -.5*pi, u, -.5*pi, v);
            Lmin = L;
        }
        if (lp_rm_slm_rp(-x, -y, phi, t, u, v) && Lmin > (L = fabs(t) + fabs(u) + fabs(v))) // timeflip + reflect
            path = reeds_shepp_path_c(
                reeds_shepp_path_type[17], -t, .5*pi, -u, .5*pi, -v);
    }
    
    reeds_shepp_path_c get_path(const state_t& s)
    {
      reeds_shepp_path_c path;
      float x = s[0], y=s[1], phi=s[2];
      CSC(x, y, phi, path);
      CCC(x, y, phi, path);
      CCCC(x, y, phi, path);
      CCSC(x, y, phi, path);
      CCSCC(x, y, phi, path);
      return path;
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

