#ifndef __system_h__
#define __system_h__

#include "dynamical_system.h"
#include "map.h"
#include <cmath>
#include <cfloat>
#include <vector>
#include <ostream>
#include <string.h>
#include <cstdlib>
using namespace std;


template<size_t N>
class region_c
{
    public:
        double s[N];
        double c[N];
        double color[4];

        region_c()
        {
            for(size_t i=0; i<N; i++)
                s[i] = c[i] = 0;

            color[0] = 1;
            color[1] = 0;
            color[2] = 0;
            color[3] = 0.01;
        }

        region_c(const double* cin, const double* sin, const double* color_in=NULL)
        {
            for(size_t i=0; i< N; i++)
            {
                s[i] = sin[i];
                c[i] = cin[i];
            }
            if(color_in)
            {
                for(int i : range(0,4))
                    color[i] = color_in[i];
            }
            else
            {
                color[0] = 1;
                color[1] = 0;
                color[2] = 0;
                color[3] = 0.01;
            }
        }

        virtual bool is_inside(const state_c<N>& sin, bool only_xy=false) const
        {
            int Nm = N;
            if(only_xy)
                Nm = 2;
            for(size_t i=0; i<N; i++)
            {
                if(fabs(sin.x[i] - c[i]) > s[i]/2.0)
                    return false;
            }
            return true;
        }

        virtual bool is_inside(const double sin[N], bool only_xy=false) const
        {
            int Nm = N;
            if(only_xy)
                Nm = 2;
            for(size_t i=0; i<Nm; i++)
            {
                if(fabs(sin[i] - c[i]) > s[i]/2.0)
                    return false;
            }
            return true;
        }


        virtual int get_plotter_state(double* cp, double* sp, int dim=3)
        {
            int d = min(dim,(int)N);
            memcpy(cp, c, d*sizeof(double));
            memcpy(sp, s, d*sizeof(double));
            return 0;
        }
};

template<size_t dim_t>
class cost_c
{
    public:
        const static size_t dim = dim_t;
        vector<double> val;

        cost_c(double x=FLT_MAX/2){
            val = vector<double>(dim, x);
        }
        cost_c(vector<double>& xin) : val(xin) {}
        cost_c(const cost_c& c2) : val(c2.val){}
        cost_c(double c, int d){
            val = vector<double>(dim, FLT_MAX/2);
            val[d] = c;
        }
        virtual ~cost_c(){}

        double& operator[](int d){
            return val[d];
        }
        virtual cost_c& operator+=(const cost_c& rhs)
        {
            for(size_t i=0; i<dim; i++)
                val[i] += rhs.val[i];
            return *this;
        }
        virtual cost_c operator+(const cost_c& c2) const
        {
            cost_c toret = *this;
            for(size_t i=0; i<dim; i++)
                toret.val[i] += c2.val[i];
            return toret;
        }
        virtual bool operator<(const cost_c& rhs) const
        {
            for(size_t i=0; i<dim; i++){
                if(val[i] < rhs.val[i])
                    return true;
                else if(val[i] > rhs.val[i])
                    return false;
            }
            return true;
        }
        virtual bool operator>(const cost_c& rhs) const
        {
            return !(*this < rhs);
        }
        virtual double difference(const cost_c& c2) const
        {
            double t1 = 0;
            for(size_t i=0; i<dim; i++)
                t1 += (val[i]-c2.val[i])*(val[i]-c2.val[i]);
            return sqrt(t1);
        }
        virtual ostream& print(ostream& os=cout, const char* prefix=NULL, const char* suffix=NULL) const
        {
            if(prefix)
                os<<prefix;
            for(auto& vi : val)
                os<<vi<<",";
            if(suffix)
                os<<suffix;
            return os;
        }
};


template<class dynamical_system_tt, class map_tt, class region_tt, class cost_tt>
class system_c
{
    public:
        typedef dynamical_system_tt dynamical_system_t;
        typedef map_tt map_t;
        typedef cost_tt cost_t;

        typedef typename dynamical_system_t::state_t state;
        typedef typename dynamical_system_t::control_t control;
        typedef typename dynamical_system_t::opt_data_t opt_data_t;
        typedef typename dynamical_system_t::trajectory_t trajectory;
        const static size_t N = dynamical_system_t::state_t::N;
        typedef region_tt region_t;

        map_t obstacle_map;
        dynamical_system_t dynamical_system;

        region_t operating_region;
        region_t goal_region;
        vector<region_t> heuristic_sampling_regions;
        double heuristic_sampling_probability;

        system_c(){
            heuristic_sampling_probability = 0.5;
        };
        ~system_c(){}

        virtual int get_key(const state& s, double* key)
        {
            for(size_t i=0; i<N; i++)
                key[i] = (s.x[i] - operating_region.c[i])/operating_region.s[i] + 0.5;
            //cout<<"key: "<< key[0]<<" "<<key[1]<<" "<<key[2]<<endl;
            return 0;
        }

        virtual int get_plotter_state(const state& s, double* ps)
        {
            return dynamical_system.get_plotter_state(s, ps);
        }

        virtual bool is_in_collision(const state& s)
        {
#if 1
            if(operating_region.is_inside(s, true))
                return obstacle_map.is_in_collision(s.x);
            else
                return true;
#else
            return obstacle_map.is_in_collision(s.x);
#endif
        }

        double get_goal_cost(const state& s)
        {
            state goal_state(goal_region.c);
            return s.dist(goal_state);
        }
        virtual bool is_in_goal(const state& s)
        {
            return goal_region.is_inside(s);
        }

        virtual int sample_state(state& s, bool sample_from_free=false)
        {
            if(sample_from_free)
            {
                if(obstacle_map.sample_free_space(s.x))
                    return 1;
            }
            else
            {
                bool found_free_state = false;
                while(!found_free_state)
                {
                    double p = RANDF;
                    region_t* r = &operating_region;
                    if(heuristic_sampling_regions.size())
                    {
                        if(p > heuristic_sampling_probability)
                            r = &operating_region;
                        else
                        {
                            int which = RANDF*heuristic_sampling_regions.size();
                            r = &(heuristic_sampling_regions[which]);
                        }
                    }

                    dynamical_system.sample_state(r->c, r->s, s.x);
                    found_free_state = !is_in_collision(s);
                }
            }
            return 0;
        }
        virtual int sample_in_goal(state& s)
        {
            bool found_free_state = false;
            while(!found_free_state)
            {
                for(size_t i=0; i<N; i++)
                    s.x[i] = goal_region.c[i] + (RANDF-0.5)*goal_region.s[i];
                found_free_state = !is_in_collision(s);
            }
            return 0;
        }

        int copy_array(const double* xin, double* xout, int dim)
        {
            memcpy(xout, xin, sizeof(double)*dim);
            return 0;
        }
        virtual bool is_safe_trajectory(const trajectory& traj)
        {
            if(traj.states.empty())
                return true;

            int drop_counter=0;
            for(auto& s : traj.states)
            {
                if(drop_counter %10 == 0)
                {
                    if(is_in_collision(s))
                        return false;
                }
                drop_counter++;
            }
            return true;
        }

        virtual int extend_to(const state& si, const state& sf, bool check_obstacles,
                trajectory& traj, opt_data_t& opt_data)
        {
            int res = dynamical_system.extend_to(si, sf, traj, opt_data);
            if(res)
            {
                traj.clear();
                return res;
            }
            if(check_obstacles)
            {
                if(is_safe_trajectory(traj))
                    res = 0;
                else
                {
                    traj.clear();
                    res = 1;
                }
            }
            return res;
        }

        virtual int evaluate_extend_cost(const state& si, const state& sf,
                opt_data_t& opt_data, cost_t& extend_cost)
        {
            double total_variation = dynamical_system.evaluate_extend_cost(si, sf, opt_data);
            if(total_variation < 0)
                return 1;
            extend_cost = cost_t();
            extend_cost[extend_cost.dim-1] = total_variation;
            return 0;
        }

        virtual cost_t get_state_cost(const state& s)
        {
            cost_t t1;
            double& t2 = t1.val.back();
            t2 = obstacle_map.get_state_cost(s.x);
            return t1;
        };

        virtual cost_t evaluate_trajectory_cost(trajectory& traj)
        {
            cost_t c(0);
            for(auto& ps : traj.states)
                c += get_state_cost(state(ps));
            return c;
        }

        virtual cost_t get_zero_cost()
        {
            return cost_t(0);
        }
        virtual cost_t get_inf_cost()
        {
            return cost_t(FLT_MAX/2);
        }

        void test_extend_to()
        {
            dynamical_system.test_extend_to();       
        }
};

#endif
