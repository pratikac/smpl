#include <iostream>
#include <ctime>

#include "dubins.h"
#include "dubins_velocity.h"
#include "single_integrator.h"
#include "double_integrator.h"
#include "rrts.h"
#include "brrts.h"
using namespace std;

int test_single_integrator()
{
    typedef system_c<single_integrator_c<3>, map_c<3>, region_c<3>, cost_c<1> > system_t;
    //typedef system_c<dubins_c, map_c<3>, region_c<3>, cost_c<1> > system_t;

    typedef system_t::state state;
    typedef typename system_t::control control;
    typedef typename system_t::trajectory trajectory;
    typedef typename system_t::region_t region;

    lcm_t *lcm          = bot_lcm_get_global(NULL);
    bot_lcmgl_t *lcmgl  = bot_lcmgl_init(lcm, "plotter");
    bot_lcmgl_line_width(lcmgl, 2.0);
    bot_lcmgl_switch_buffer(lcmgl);

    rrts_c<vertex_c<system_t>, edge_c<system_t> > rrts(lcmgl);

    float zero[3] = {0};
    float size[3] = {100,100,2*M_PI};
    rrts.system.operating_region = region(zero, size);

    float gc[3] = {10,10,0};
    float gs[3] = {1,1,0.1*M_PI};
    state goal_state(gc);
    rrts.system.goal_region = region(gc,gs);

    state origin(zero);
    rrts.initialize(origin);

    tt clock;
    clock.tic();
    int max_iterations = 1e3, diter=max_iterations/10;
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
    cout<<"time: "<< clock.toc() <<" [ms]"<<endl;

    rrts.plot_environment();
    rrts.plot_tree();
    rrts.plot_best_trajectory();
    bot_lcmgl_switch_buffer(lcmgl);
    cout<<rrts.get_best_cost().val[0]<<endl;


    return 0;
}

int test_double_integrator()
{
    typedef system_c<double_integrator_c, map_c<4>, region_c<4>, cost_c<1> > system_t;

    typedef system_t::state state;
    typedef typename system_t::control control;
    typedef typename system_t::trajectory trajectory;
    typedef typename system_t::region_t region;

    lcm_t *lcm          = bot_lcm_get_global(NULL);
    bot_lcmgl_t *lcmgl  = bot_lcmgl_init(lcm, "plotter");
    bot_lcmgl_switch_buffer(lcmgl);

    rrts_c<vertex_c<system_t>, edge_c<system_t> > rrts(lcmgl);

    float zero[4] = {0};
    float size[4] = {25, 25, 5, 5};
    rrts.system.operating_region = region(zero, size);
    rrts.goal_sample_freq = 0.01;

    float gc[4] = {0, 0, 0, 0};
    float gs[4] = {0.1,0.1,0.1,0.1};
    state goal_state(gc);
    rrts.system.goal_region = region(gc,gs);

    rrts.system.test_extend_to();
    return 0;

    float ss[4] = {5,5,5,5};
    state start_state(ss);
    rrts.initialize(start_state, lcmgl);

    tt clock;
    clock.tic();
    int max_iterations = 1e5, diter=1000;
    for(int i=0; i<max_iterations; i++)
    {
        rrts.iteration();
        if(i%diter == 0)
        {
            cout<<i<<" "<<rrts.get_best_cost().val[0]<<endl;
            rrts.plot_tree();
            rrts.plot_best_trajectory();
            bot_lcmgl_switch_buffer(lcmgl);
            getchar();
        }
    }
    rrts.plot_tree();
    rrts.plot_best_trajectory();
    bot_lcmgl_switch_buffer(lcmgl);
    cout<<"time: "<< clock.toc() <<" [ms]"<<endl;
    cout<<rrts.get_best_cost().val[0]<<endl;

    return 0;
}

int test_brrts()
{
    srand(time(NULL));

    //typedef system_c<single_integrator_c<3>, map_c<3>, region_c<3>, cost_c<1> > system_t;
    typedef system_c<double_integrator_c, map_c<4>, region_c<4>, cost_c<1> > system_t;

    typedef system_t::state state;
    typedef typename system_t::control control;
    typedef typename system_t::trajectory trajectory;
    typedef typename system_t::region_t region;

    lcm_t *lcm          = bot_lcm_get_global(NULL);
    bot_lcmgl_t *lcmgl  = bot_lcmgl_init(lcm, "plotter");
    bot_lcmgl_switch_buffer(lcmgl);

    brrts_c<bvertex_c<system_t>, bedge_c<system_t> > brrts(lcmgl);

    float zero[4] = {0};
    float size[4] = {25, 25, 10, 10};
    brrts.system.operating_region = region(zero, size);

    float gc[4] = {10,5,0,0};
    float gs[4] = {0.1,0.1,0.1,0.1};
    state goal_state(gc);
    brrts.system.goal_region = region(gc,gs);

    state origin(zero);
    brrts.initialize(origin, lcmgl);

    tt clock;
    clock.tic();
    int max_iterations = 1e2, diter=max_iterations/10;
    for(int i=0; i<max_iterations; i++)
    {
        brrts.iteration();
        if(i%diter == 0)
        {
            cout<<i<<" "<< brrts.get_best_cost().val[0]<<endl;
            brrts.plot_tree();
            brrts.plot_best_trajectory();
            bot_lcmgl_switch_buffer(lcmgl);
        }
    }
    cout<<"time: "<< clock.toc() <<" [ms]"<<endl;
    cout<< brrts.get_best_cost().val[0]<<endl;

    return 0;
}

int test_dubins_velocity()
{
    typedef system_c<dubins_velocity_c, map_c<4>, region_c<4>, cost_c<1> > system_t;

    typedef system_t::state state;
    typedef typename system_t::control control;
    typedef typename system_t::trajectory trajectory;
    typedef typename system_t::region_t region;

    lcm_t *lcm          = bot_lcm_get_global(NULL);
    bot_lcmgl_t *lcmgl  = bot_lcmgl_init(lcm, "plotter");
    bot_lcmgl_line_width(lcmgl, 2.0);
    bot_lcmgl_switch_buffer(lcmgl);

    rrts_c<vertex_c<system_t>, edge_c<system_t> > rrts(lcmgl);

    float zero[4] = {0};
    float size[4] = {25,25,2*M_PI, 5};
    rrts.system.operating_region = region(zero, size);

    float gc[4] = {10,10,0,1};
    float gs[4] = {1,1,0.1*M_PI,0.1};
    state goal_state(gc);
    rrts.system.goal_region = region(gc,gs);

    //rrts.system.test_extend_to();
    //return 0;

    state origin(zero);
    rrts.initialize(origin);

    tt clock;
    clock.tic();
    int max_iterations = 1e3, diter=1;
    trajectory traj;
    for(int i=0; i<max_iterations; i++)
    {
        rrts.iteration();
        if(i%diter == 0)
            cout<<i<<" "<<rrts.get_best_cost().val[0]<<endl;

        rrts.plot_tree();
        rrts.plot_best_trajectory();
        bot_lcmgl_switch_buffer(lcmgl);
        getchar();
    }
    cout<<"time: "<< clock.toc() <<" [ms]"<<endl;

    rrts.plot_tree();
    rrts.plot_best_trajectory();
    bot_lcmgl_switch_buffer(lcmgl);
    cout<<rrts.get_best_cost().val[0]<<endl;

    return 0;
}

int main()
{

    //test_single_integrator();
    test_double_integrator();
    //test_brrts();
    //test_dubins_velocity();
    return 0;
};
