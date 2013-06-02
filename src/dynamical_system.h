#ifndef __dynamical_system_h__
#define __dynamical_system_h__

#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;

#define SQ(x)   ((x)*(x))

template<size_t N>
class state_t
{
  public:
    float x[N];
    const static size_t dim = N;

    state_t()
    {
      for(size_t i=0; i<N; i++)
        x[i] = 0;
    }
    state_t(const state_t& s)
    {
      for(size_t i=0; i<N; i++)
        x[i] = s.x[i];
    }
    state_t(const float* sin)
    {
      for(size_t i=0; i<N; i++)
        x[i] = sin[i];
    }
    state_t& operator=(const state_t& sin)
    {
      if(this == &sin)
        return *this;
      for(size_t i=0; i<N; i++)
        x[i] = sin.x[i];
      return *this;
    }
    
    float operator[](const size_t i) const {return x[i];}
    
    virtual ostream& print(ostream& os=cout, const char* prefix=NULL, const char* suffix=NULL) const
    {
      if(prefix)
        os<<prefix<<" (";
      else
        os<<" (";
      for(size_t i=0; i<N-1; i++)
        os<<x[i]<<",";
      os<<x[N-1]<<")";
      if(suffix)
        os<<suffix;
      return os;
    }
    
    float dist(const state_t& s, bool only_xy=false) const
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
class control_t : public state_t<M>
{
  public:
    control_t() : state_t<M>() {}
    control_t(const control_t& c) : state_t<M>(c) {};
    control_t(const float* cin) : state_t<M>(cin) {};
};

template<class s_t, class c_t>
class trajectory_t
{
  public:
    vector<s_t> states;
    vector<c_t> controls;
    float total_variation;
     
    trajectory_t() : total_variation(0){}
    
    // does not clear memory, explicitly
    // call trajectory.clear() to deallocate memory
    ~trajectory_t(){}
    int clear()
    {
      total_variation = 0;
      states.clear();
      controls.clear();
      return 0;
    }
    int pop_front(int how_many)
    {
      states = vector<s_t>(states.begin()+how_many, states.end());
      controls = vector<s_t>(controls.begin()+how_many, controls.end());
      return 0;
    }
    int append(trajectory_t& t2)
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
        s.print();
      return 0;
    }
    int print_controls(const char* prefix="")
    {
      cout<<prefix<<endl;
      for(auto& c : controls)
        c.print();
      return 0;
    }
};

class optimization_data_t
{
  public:
    virtual ~optimization_data_t(){};
};

template<class s_t, class c_t, class opt_data_t>
class dynamical_system_t
{
  public:
    typedef trajectory_t<s_t, c_t> traj_t;
    typedef s_t state_t;
    typedef c_t control_t;

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
    
    virtual int extend_to(const s_t* si, const s_t* sf, traj_t& traj, opt_data_t* opt_data)=0;
    virtual float evaluate_extend_cost(const s_t* si, const s_t* sf, opt_data_t*& opt_data)=0;
};

#endif
