#include <bot_core/bot_core.h>

#include "../dubins.h"
#include <bot_vis/gl_util.h>
#include <bot_lcmgl_client/lcmgl.h>
#include <unistd.h>

int main()
{
    lcm_t *lcm          = bot_lcm_get_global(NULL);
    bot_lcmgl_t* lcmgl = bot_lcmgl_init(lcm, "plotter");
    bot_lcmgl_line_width(lcmgl, 2.0);
    bot_lcmgl_switch_buffer(lcmgl);

    double r = 10;
    double th = M_PI/4.;
    double eps = 0.1;

    while(eps < th - eps)
    {
        double sth = sin(th);
        double cth = cos(th);
        double L = 1.0;

        double turning_radii[] = {r};
        dubins_c dubins(turning_radii, 1);

        state_c<3> s0;
        for(int i=0; i< 3; i++)
            s0.x[i] = 0;

        state_c<3> a, ap, b, bn, sf, sf2;
        a.x[0] = s0.x[0] + r*sth;
        a.x[1] = s0.x[1] + r*(1-cth);
        a.x[2] = s0.x[2] + th;

        ap.x[0] = s0.x[0] + r*sin(th-eps);
        ap.x[1] = s0.x[1] + r*(1-cos(th-eps));
        ap.x[2] = s0.x[2] + (th-eps);

        b.x[0] = a.x[0] + L*sth;
        b.x[1] = a.x[1] + L*cth;
        b.x[2] = a.x[2];

        bn.x[0] = b.x[0] + r*(cos(th) - cos(th+eps));
        bn.x[1] = b.x[1] + r*(sin(th+eps) - sin(th));
        bn.x[2] = b.x[2] - eps;

        sf.x[0] = b.x[0] + r*sth;
        sf.x[1] = b.x[1] + r*(1-cth);
        sf.x[2] = s0.x[0];

        sf2.x[0] = r;
        sf2.x[1] = r;
        sf2.x[2] = 2*th;

        auto z0 = ap;
        auto zf = sf2;

        dubins_optimization_data_c opt_data;
        dubins_c::trajectory_t traj;
        dubins.extend_to(z0, zf, traj, opt_data);
        
        bot_lcmgl_color4f(lcmgl, 1, 1, 0, 1);
        bot_lcmgl_point_size(lcmgl, 2);
        bot_lcmgl_begin(lcmgl, GL_POINTS);
        for(int i=0; i < traj.states.size(); i++)
        {
            auto c = traj.states[i];
            bot_lcmgl_vertex3d(lcmgl, c.x[0], c.x[1], 0);
            //printf("c: (%f,%f)\n", c[0], c[1]);
        }
        bot_lcmgl_end(lcmgl);
        bot_lcmgl_switch_buffer(lcmgl);
   
        traj.clear();

        bot_lcmgl_color4f(lcmgl, 0, 1, 0, 1);
        bot_lcmgl_text(lcmgl, s0.x, "s0");
        bot_lcmgl_text(lcmgl, a.x, "a");
        bot_lcmgl_text(lcmgl, ap.x, "ap");
        bot_lcmgl_text(lcmgl, b.x, "b");
        bot_lcmgl_text(lcmgl, bn.x, "bn");
        bot_lcmgl_text(lcmgl, sf.x, "sf");
        bot_lcmgl_text(lcmgl, sf2.x, "sf2");

        usleep(0.2e6);
        eps += 0.02;
    }
};
