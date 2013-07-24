#ifndef __dynamical_system_h__
#define __dynamical_system_h__

#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;

#define SQ(x)   ((x)*(x))

template<size_t N_t>
class state_c
{
  public:
    const static size_t N = N_t;
    float x[N];

    state_c()
    {
      for(size_t i=0; i<N; i++)
        x[i] = 0;
    }
    state_c(const state_c& s)
    {
      for(size_t i=0; i<N; i++)
        x[i] = s.x[i];
    }
    state_c(const float* sin)
    {
      for(size_t i=0; i<N; i++)
        x[i] = sin[i];
    }
    state_c& operator=(const state_c& sin)
    {
      if(this == &sin)
        return *this;
      for(size_t i=0; i<N; i++)
        x[i] = sin.x[i];
      return *this;
    }
    state_c operator+(const state_c& s2) const
    {
      state_c toret = *this;
      for(size_t i=0; i<N; i++)
        toret.x[i] += s2.x[i];
      return toret;
    }
    float operator[](const size_t i) const 
    {
      return x[i];
    }
    virtual ostream& print(ostream& os=cout, const char* prefix=NULL, const char* suffix=NULL) const
    {
      if(prefix)
        os<<prefix<<" [";
      else
        os<<" [";
      for(size_t i=0; i<N-1; i++)
        os<<x[i]<<",";
      os<<x[N-1]<<"]";
      if(suffix)
        os<<suffix;
      return os;
    }
    float dist(const state_c& s, bool only_xy=false) const
    {
      size_t len = N;
      if(only_xy)
        len = 2;
      float t=0;
      for(size_t i=0; i<len; i++)
        t = t + SQ(x[i]-s.x[i]);
      return sqrt(t);
    }
};

template<size_t M>
class control_c : public state_c<M>
{
  public:
    control_c() : state_c<M>() {}
    control_c(const control_c& c) : state_c<M>(c) {};
    control_c(const float* cin) : state_c<M>(cin) {};
};

template<class state_t, class control_t>
class trajectory_c
{
  public:
    vector<state_t> states;
    vector<control_t> controls;
    float total_variation;
     
    trajectory_c() : total_variation(0){}
    
    // does not clear memory, explicitly
    // call trajectory.clear() to deallocate memory
    ~trajectory_c(){}
    int clear()
    {
      total_variation = 0;
      states.clear();
      controls.clear();
      return 0;
    }
    int pop_front(int how_many)
    {
      states = vector<state_t>(states.begin()+how_many, states.end());
      controls = vector<control_t>(controls.begin()+how_many, controls.end());
      return 0;
    }
    int append(trajectory_c& t2)
    {
      total_variation += t2.total_variation;
      states.insert(states.end(), t2.states.begin(), t2.states.end());
      controls.insert(controls.end(), t2.controls.begin(), t2.controls.end());
      return 0;
    }
    void reverse()
    {
      std::reverse(states.begin(), states.end());
      std::reverse(controls.begin(), controls.end());
      return;
    }
    int print(const char* prefix="traj")
    {
      print_states("states");
      print_controls("controls");
      return 0;
    }
    int print_states(const char* prefix="")
    {
      cout<<prefix<<endl;
      for(auto& s : states)
        s.print(cout, "", "\n");
      return 0;
    }
    int print_controls(const char* prefix="")
    {
      cout<<prefix<<endl;
      for(auto& c : controls)
        c.print(cout, "", "\n");
      return 0;
    }
};

class optimization_data_c
{
  public:
    virtual ~optimization_data_c(){};
};

template<class state_tt, class control_tt, class opt_data_tt>
class dynamical_system_c
{
  public:
    typedef state_tt state_t;
    typedef control_tt control_t;
    typedef opt_data_tt opt_data_t;
    typedef trajectory_c<state_t, control_t> trajectory_t;

    int modulo_mpi_pi(float& th)
    {
      while(th < -M_PI)
        th = th + 2*M_PI;
      while(th > M_PI)
        th = th - 2*M_PI;
      return 0;
    }
    int modulo_zero_2pi(float& th)
    {
      while(th < 0)
        th = th + 2*M_PI;
      while(th > 2*M_PI)
        th = th - 2*M_PI;
      return 0;
    }
    
    virtual int extend_to(const state_t& si, const state_t& sf, trajectory_t& traj, opt_data_t& opt_data)=0;
    virtual float evaluate_extend_cost(const state_t& si, const state_t& sf, opt_data_t& opt_data)=0;
};

#endif
