#include <iostream>
#include <ctime>

#include "dubins.h"
#include "single_integrator.h"
#include "double_integrator.h"
#include "map.h"
#include "system.h"
#include "rrts.h"
using namespace std;

#if 0
int main()
{
  //typedef system_c<single_integrator_c<3>, map_c<3>, cost_c<1> > system_t;
  //typedef system_c<dubins_c, map_c<3>, cost_c<1> > system_t;
  
  typedef system_t::state state;
  typedef typename system_t::control control;
  typedef typename system_t::trajectory trajectory;
  typedef typename system_t::region_t region;
  
  rrts_c<system_t> rrts;

  float zero[3] = {0};
  float size[3] = {100,100,2*M_PI};
  rrts.system.operating_region = region(zero, size);

  float gc[3] = {10,10,0};
  float gs[3] = {1,1,0.1*M_PI};
  state goal_state(gc);
  rrts.system.goal_region = region(gc,gs);
  
  state origin(zero);
  rrts.initialize(origin);

  time_t ts=time(0), te;
  int max_iterations = 1e4, diter=max_iterations/10;
  trajectory traj;
  for(int i=0; i<max_iterations; i++)
  {
    rrts.iteration();
    if(i%diter == 0)
      cout<<i<<" "<<rrts.get_best_cost().val[0]<<endl;
    //cout<<"check_tree: "<< rrts.check_tree() << endl;
#if 0
    if(rrts.system.is_in_goal(rrts.root->state))
      break;
    if(i % 100 == 0)
    { 
      cout<<i<<" "<<rrts.get_best_cost().val<<endl;
      rrts.switch_root(25, traj);
      cout<<"is safe: "<< rrts.system.is_safe_trajectory(traj)<<endl;
      traj.clear();
      cout<<"switched root: ";
      rrts.root->state.print();
      cout<<endl;
    }
#endif
  }
  cout<<rrts.get_best_cost().val[0]<<endl;
  trajectory best_traj;
  rrts.get_best_trajectory(best_traj);
  best_traj.print();
  
  cout<<"time: "<< difftime(time(0), ts)<<endl;

  return 0;
}

#else

int main()
{
  typedef system_c<double_integrator, map_c<4>, cost_c<1> > system_t;
  
  typedef system_t::state state;
  typedef typename system_t::control control;
  typedef typename system_t::trajectory trajectory;
  typedef typename system_t::region_t region;
  
  rrts_c<system_t> rrts;

  float zero[4] = {0};
  float size[4] = {100,100,5,5};
  rrts.system.operating_region = region(zero, size);

  float gc[4] = {10,10,1,0};
  float gs[4] = {1,1,0.1,0.1};
  state goal_state(gc);
  rrts.system.goal_region = region(gc,gs);
  
  state origin(zero);
  rrts.initialize(origin);

  time_t ts=time(0), te;
  int max_iterations = 1e4, diter=max_iterations/10;
  trajectory traj;
  for(int i=0; i<max_iterations; i++)
  {
    rrts.iteration();
    if(i%diter == 0)
      cout<<i<<" "<<rrts.get_best_cost().val[0]<<endl;
    //cout<<"check_tree: "<< rrts.check_tree() << endl;
  }
  cout<<rrts.get_best_cost().val[0]<<endl;
  trajectory best_traj;
  rrts.get_best_trajectory(best_traj);
  best_traj.print();
  
  cout<<"time: "<< difftime(time(0), ts)<<endl;

  return 0;
};
#endif
