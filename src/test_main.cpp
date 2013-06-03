#include <iostream>
#include <ctime>

#include "dubins.h"
#include "single_integrator.h"
#include "map.h"
#include "system.h"
#include "rrts.h"
using namespace std;

int main()
{
  //typedef system_c<single_integrator_c<3>, map_c<3>, cost_c> system_t;
  typedef system_c<dubins_c, map_c<3>, cost_c> system_t;
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
  int max_iterations = 1000, diter=max_iterations/10;
  trajectory traj;
  for(int i=0; i<max_iterations; i++)
  {
    cout<<i<<" "<<rrts.get_best_cost().val<<endl;
    //cout<<"check_tree: "<< rrts.check_tree() << endl;
    for(int j=0; j< 10; j++)
      rrts.iteration();
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
  cout<<rrts.get_best_cost().val<<endl;
  //traj.print();
  
  cout<<"time: "<< difftime(time(0), ts)<<endl;

  return 0;
}
