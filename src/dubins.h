#ifndef __dubins_h__
#define __dubins_h__

#include <float.h>
#include <iostream>
#include "dynamical_system.h"
using namespace std;

#include "dbs.h"

class dubins_optimization_data_c : public optimization_data_c
{
    public:
        int turning_radius;
        dubins_optimization_data_c() : turning_radius(-1){}
        ~dubins_optimization_data_c(){}
};

class dubins_c : public dynamical_system_c<state_c<3>, control_c<1>, dubins_optimization_data_c >
{
    public:
        typedef typename dynamical_system_c::state_t state_t;
        typedef typename dynamical_system_c::control_t control_t;
        typedef dynamical_system_c::trajectory_t trajectory_t;
        typedef dubins_optimization_data_c dubins_optimization_data_t;

        double delta_distance;

        int num_turning_radii;
        vector<double> turning_radii;

        dubins_c()
        {
            num_turning_radii = 1;
            turning_radii.reserve(num_turning_radii);

            delta_distance = 0.05;
            turning_radii[0] = 10;
            //turning_radii[1] = 6;
            //turning_radii[2] = 8;
        };
        
        dubins_c(double* radii, int num_tr_)
        {
            num_turning_radii = num_tr_;
            turning_radii.reserve(num_turning_radii);
            for(int i=0; i< num_tr_; i++)
                turning_radii[i] = radii[i];

            delta_distance = 0.05;
        };

        int get_plotter_state(const state_t& s, double* ps)
        {
            ps[0] = s[0];
            ps[1] = s[1];
            ps[2] = 0;
            return 0;
        }

        int sample_state(double* center, double* size, double* s)
        {
            for(int i : range(0,3))
                s[i] = center[i] + (RANDF-0.5)*size[i];
            return 0;
        }

        static int dubins_path_callback(double q[3], double x, void* user_data)
        {
            state_t s(q);
            trajectory_t* traj = (trajectory_t*) user_data;
            
            traj->states.push_back(s);
            return 0;
        }

        // not populating controls for now
        trajectory_t dubins_path_get_discretized(DubinsPath* path)
        {
            trajectory_t traj;
            traj.clear();
            traj.dt = delta_distance;
            traj.t0 = 0;
            traj.total_variation = dubins_path_length(path);
            
            dubins_path_sample_many(path, dubins_path_callback, delta_distance, &traj);
            
            return traj;
        }

        int extend_to(const state_t& si, const state_t& sf, trajectory_t& traj,
                dubins_optimization_data_t& opt_data)
        {
            bool return_trajectory = true;
            if(opt_data.turning_radius < 0)
            {
                if(evaluate_extend_cost(si, sf, opt_data) < 0)
                    return 1;
            }
            double tr = turning_radii[opt_data.turning_radius];
            
            DubinsPath path;
            dubins_init(si.x, sf.x, tr, &path);
            traj = dubins_path_get_discretized(&path);
            
            return 0;
        }

        double evaluate_extend_cost(const state_t& si, const state_t& sf,
                dubins_optimization_data_t& opt_data)
        {
            trajectory_t traj;
            if(opt_data.turning_radius >= 0)
            {
                double tr = turning_radii[opt_data.turning_radius];
                DubinsPath path;
                dubins_init(si.x, sf.x, tr, &path);
                return dubins_path_length(&path);
            }
            else
            {
                double min_cost = DBL_MAX;

                int& best_turning_radius = opt_data.turning_radius;
                best_turning_radius = -1;

                bool return_trajectory = false;
                for(int i=num_turning_radii-1; i >=0; i--)
                {
                    double tr = turning_radii[i];
                    
                    DubinsPath path;
                    dubins_init(si.x, sf.x, tr, &path);
                    double cost = dubins_path_length(&path);

                    if(cost > 0)
                    {
                        if(cost < min_cost)
                        {
                            min_cost = cost;
                            best_turning_radius = i;
                        }
                    }
                }
                if((min_cost < 0) || (min_cost > DBL_MAX/2.0))
                    return -1;
                return min_cost;
            }
        }

        void test_extend_to()
        {
            trajectory_t traj;
            double zero[3] ={0};
            state_t origin(zero);
            double goal[3] = {10, 10, M_PI/2. + 0.01};
            state_t sr(goal);
            sr.print(cout, "sampled:","\n");

            dubins_optimization_data_c opt_data;
            extend_to(origin, sr, traj, opt_data);
            cout<<"cost: "<< traj.total_variation<<endl;
            traj.print();
            traj.clear();
        }
};

#endif
