#ifndef __double_integrator_h__
#define __double_integrator_h__

#include "dynamical_system.h"
#include <cassert>

class double_integrator_optimization_data_c : public optimization_data_c
{
    public:
        double T1, T2, T;
        double ts1, ts2;
        double u1, u2;
        bool is_initialized;

        double_integrator_optimization_data_c() : is_initialized(false) {}
        ~double_integrator_optimization_data_c(){}
};

class double_integrator_c : public dynamical_system_c<state_c<4>, control_c<2>, double_integrator_optimization_data_c>
{
    public:
        typedef dynamical_system_c<state_c<4>, control_c<2>, optimization_data_c> dynamical_system_t;
        typedef typename dynamical_system_t::state_t state_t;
        typedef typename dynamical_system_t::control_t control_t;
        typedef typename dynamical_system_t::trajectory_t trajectory_t;
        typedef double_integrator_optimization_data_c double_integrator_opt_data_t;

        double umm;
        double_integrator_c() : umm(1) {}

        int get_plotter_state(const state_t& s, double* ps)
        {
            ps[0] = s[0];
            ps[1] = s[1];
            ps[2] = 0; //s[2];
            return 0;
        }

        int sample_state(double* center, double* size, double* s)
        {
            for(int i : range(0,4))
                s[i] = center[i] + (RANDF-0.5)*size[i];
            return 0;
        }
        
        int extend_to(const state_t& si, const state_t& sf,
                trajectory_t& traj, double_integrator_opt_data_t& opt_data)
        {
            if(evaluate_extend_cost(si, sf, opt_data) < 0)
                return 1;

            traj.clear();
            traj.total_variation = opt_data.T;

            double si1[2] = {si[0], si[2]};
            double sf1[2] = {sf[0], sf[2]};
            double si2[2] = {si[1], si[3]};
            double sf2[2] = {sf[1], sf[3]};

            double T = opt_data.T, T1 = opt_data.T1, T2 = opt_data.T2;
            double ts1, ts2;
            double um = umm;
            double u1, u2;

            if(!opt_data.is_initialized)
            {
                double g = 1.0;
                // dim 2 needs to slow down
                if( fabs(T-T1) < 1e-6)
                {
                    g = get_gain(si2, sf2, um, T);
                    if(g < 0)
                        return 1;
                    get_time(si2, sf2, fabs(g*um), u2, ts2);
                    get_time(si1, sf1, fabs(um), u1, ts1);
                }
                // dim 1 needs to slow down
                else
                {
                    g = get_gain(si1, sf1, um, T);
                    if(g < 0)
                        return 1;
                    get_time(si1, sf1, fabs(g*um), u1, ts1); 
                    get_time(si2, sf2, fabs(um), u2, ts2);
                }

                opt_data.ts1 = ts1;
                opt_data.ts2 = ts2;
                opt_data.u1 = u1;
                opt_data.u2 = u2;
                opt_data.is_initialized = true;
            } 
            else
            {
                ts1 = opt_data.ts1;
                ts2 = opt_data.ts2;
                u1 = opt_data.u1;
                u2 = opt_data.u2;
            }
            double t = 0, dt = 0.1;
            traj.dt = dt;
            state_t sc = si;

            control_t cc;
            int num_steps = T/dt;
            traj.states.reserve(num_steps);
            traj.controls.reserve(num_steps);

            traj.states.push_back(sc);
            while(t < T)
            {
                if(t < ts1)
                    cc.x[0] = u1;
                else
                    cc.x[0] = -u1;
                if(t < ts2)
                    cc.x[1] = u2;
                else
                    cc.x[1] = -u2;

                sc.x[2] += cc.x[0]*dt;
                sc.x[0] += sc.x[2]*dt;

                sc.x[3] += cc.x[1]*dt;
                sc.x[1] += sc.x[3]*dt;
                
                traj.states.push_back(sc);
                traj.controls.push_back(cc);
                t += dt;
            }
            return 0;
        }

        //use control g*um to reach origin from x0[2] in T time units
        double get_f(const double x0[2], const double xf[2],
                const double g, const double um, const double T)
        {
            double ts, u;
            return get_time(x0, xf, fabs(g*um), u, ts) - T; 
        }
        
        double get_gain(const double x0[2], const double xf[2], const double um,
                const double T, const double eps=1e-6)
        {
            double g=0.5, gm=eps, gp=1-gm;
            double f, fm, fp;
            fm = get_f(x0, xf, gm, um, T);
            fp = get_f(x0, xf, gp, um, T);
            if(fabs(fp)<eps)
            {
                return 1;
            }
            if(fm*fp > 0)
            {
                //cout<<"new eps: "<< eps/2<<endl;
                //printf("x0: %.2f, %.2f, xf: %.2f, %.2f, T: %.2f, fm: %.2f, fp: %.2f\n",
                //        x0[0], x0[1], xf[0], xf[1], T, fm, fp);
                return get_gain(x0,xf,um,T, eps/2);
            }
            
            bool is_converged = false;
            int c = 0;
            while(!is_converged)
            {
                g = (gm+gp)/2.0;
                f = get_f(x0, xf, g, um, T);
                assert(f == f);

                if(fabs(f) < eps)
                    return g;

                if(f*fm < 0){
                    gp = g;
                    fp = get_f(x0, xf, gp, um, T);
                }
                else if(f*fp < 0){
                    gm = g;
                    fm = get_f(x0, xf, gm, um, T);
                }

                c++;
                if(c > 100)
                    cout<<"g: "<< g << " fm: "<< fm <<" fp: "<< fp <<" f: "<< f << endl;
                is_converged = (gp-gm) < eps;
            }
            return g;
        }

        // check to see if to right of switching curve
        bool is_right(const double x0[2], const double xf[2], double um)
        {
            assert(um > 0);
            double u = um;
            if(x0[1] > xf[1])
                u = -um;
            if ( x0[0] - xf[0] - 1/2/u*(xf[1]*xf[1] - x0[1]*x0[1]) > 0)
                return true;
            return false;
        }

        // control is: 0 --u-- ts --(-u1)-- T
        double get_time(const double x0[2], const double xf[2], const double um,
                double& u, double& ts)
        {
            u = um;
            if(is_right(x0,xf,um))
                u = -um;

            double D = um*um*(x0[1]*x0[1] + xf[1]*xf[1] + \
                                2*um*x0[0] - 2*um*xf[0]);
            if (D < 0)
                return -1;
            double sqD = sqrt(D);
            ts = max((x0[1]*um - 0.7071*sqD)/um/um, (x0[1]*um + 0.7071*sqD)/um/um);
            double T;
            T = max((x0[1]*um+ xf[1]*um - 1.4142*sqD)/um/um, (x0[1]*um+ xf[1]*u + 1.4142*sqD)/um/um);
            return T;
        }

        double evaluate_extend_cost(const state_t& si, const state_t& sf,
                double_integrator_opt_data_t& opt_data)
        {
            if(opt_data.is_initialized)
                return opt_data.T;

            double um=umm;

            const double si1[2] = {si[0], si[2]};
            const double sf1[2] = {sf[0], sf[2]};
            const double si2[2] = {si[1], si[3]};
            const double sf2[2] = {sf[1], sf[3]};
            double ts1, ts2;
            double u1;

            double T1 = get_time(si1, sf1, um, u1, ts1);
            double T2 = get_time(si2, sf2, um, u1, ts2);
            double T = max(T1, T2);
            if(T < 0)
                return -1;
            
            //cout<<endl;
            //si.print(cout, "si: ", "\n");
            //sf.print(cout, "sf: ", "\n");
            //cout<<"T1: "<< T1 << " T2: "<< T2 << " T: "<< T << endl;

            opt_data.T1 = T1;
            opt_data.T2 = T2;
            opt_data.T = T;

            //cout<<"T: "<< T << endl;
            return T;
        }

        void test_extend_to()
        {
            trajectory_t traj1, traj2;

            double s1t[4] = {0};
            state_t s1(s1t);
            double s2t[4] = {10, 10, 10, 10};
            state_t s2(s2t);
            s2.print(cout, "sampled:","\n");

            double_integrator_optimization_data_c opt_data;
            evaluate_extend_cost(s2, s1, opt_data);
            
            if(extend_to(s2, s1, traj1, opt_data))
                cout<<"could not connect"<<endl;
            traj1.print();
            cout<<"first trajectory finished"<<endl;
            
            /*
               double g2[4] = {12, 4, 1, -1};
               state_t sg2(g2);
               sg2.print(cout, "second goal:","\n");

               double_integrator_optimization_data_c opt_data2;
               if(extend_to(sr, sg2, traj2, opt_data2))
               cout<<"could not connect"<<endl;
               traj1.append(traj2);
               traj1.print();
               */
        }
}; 

#endif
