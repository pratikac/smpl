#ifndef __dbs_h__
#define __dbs_h__

// Copyright (c) 2008-2014, Andrew Walker
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include <math.h>
#include <assert.h>
#include <glib.h>

#define DEBUG       (1)

#define LSL (0)
#define LSR (1)
#define RSL (2)
#define RSR (3)
#define RLR (4)
#define LRL (5)

// The three segment types a path can be made up of
#define L_SEG (0)
#define S_SEG (1)
#define R_SEG (2)

// Error codes
#define EDUBOK        (0)   // No error
#define EDUBCOCONFIGS (1)   // Colocated configurations
#define EDUBPARAM     (2)   // Path parameterisitation error
#define EDUBBADRHO    (3)   // the rho value is invalid
#define EDUBNOPATH    (4)   // no connection between configurations with this word

// The segment types for each of the Path types
const int DIRDATA[][3] = {
    { L_SEG, S_SEG, L_SEG },
    { L_SEG, S_SEG, R_SEG },
    { R_SEG, S_SEG, L_SEG },
    { R_SEG, S_SEG, R_SEG },
    { R_SEG, L_SEG, R_SEG },
    { L_SEG, R_SEG, L_SEG }
};

// The various types of solvers for each of the path types
int dubins_LSL(double alpha, double beta, double d, double* outputs);
int dubins_LSR(double alpha, double beta, double d, double* outputs);
int dubins_RSL(double alpha, double beta, double d, double* outputs);
int dubins_RSR(double alpha, double beta, double d, double* outputs);
int dubins_RLR(double alpha, double beta, double d, double* outputs);
int dubins_LRL(double alpha, double beta, double d, double* outputs);

typedef int (*DubinsWord)(double, double, double, double* );
DubinsWord dubins_words[] = {
    dubins_LSL,
    dubins_LSR,
    dubins_RSL,
    dubins_RSR,
    dubins_RLR,
    dubins_LRL,
};

typedef struct
{
    double qi[3];       // the initial configuration
    double param[3];    // the lengths of the three segments
    double rho;         // model forward velocity / model angular velocity
    int type;           // path type. one of LSL, LSR, ... 
} DubinsPath;


#define UNPACK_INPUTS(alpha, beta)     \
    double sa = sin(alpha);            \
double sb = sin(beta);             \
double ca = cos(alpha);            \
double cb = cos(beta);             \
double c_ab = cos(alpha - beta);   \

#define PACK_OUTPUTS(outputs)       \
    outputs[0]  = t;                \
outputs[1]  = p;                \
outputs[2]  = q;

/**
 * Callback function for path sampling
 *
 * @note the q parameter is a configuration
 * @note the t parameter is the distance along the path
 * @note the user_data parameter is forwarded from the caller
 * @note return non-zero to denote sampling should be stopped
 */
typedef int (*DubinsPathSamplingCallback)(double q[3], double t, void* user_data);

#define DUBINS_ZERO     (-1e-9)
#define DUBINS_EPS      (1e-6)
static double dbsmod2pi(double theta)
{
    if(theta < 0 && theta > DUBINS_ZERO)
        return 0;
    return theta - 2*M_PI * floor(theta/2./M_PI);
}

int dubins_init_normalised( double alpha, double beta, double d, DubinsPath* path)
{
    if(d < DUBINS_EPS && fabs(alpha-beta) < DUBINS_EPS)
    {
        double params[3];
        int err = dubins_words[0](0, 0, d, params);
        if(err == EDUBOK)
        {
            path->param[0] = params[0];
            path->param[1] = params[1];
            path->param[2] = params[2];
            path->type = 0;
        }
        return err;
    }

    double best_cost = INFINITY;
    int    best_word;
    int    i;

    best_word = -1;
    for( i = 0; i < 6; i++ )
    {
        double params[3];
        int err = dubins_words[i](alpha, beta, d, params);
        if(err == EDUBOK)
        {
            double cost = params[0] + params[1] + params[2];
            //printf("i: %d, cost: %f\n", i, cost);
            if(cost < best_cost)
            {
                best_word = i;
                best_cost = cost;
                path->param[0] = params[0];
                path->param[1] = params[1];
                path->param[2] = params[2];
                path->type = i;
            }
        }
    }

    if(best_word == -1)
    {
        return EDUBNOPATH;
    }
    path->type = best_word;
    return EDUBOK;
}

int dubins_init(const double q0[3], const double q1[3], double rho, DubinsPath* path)
{
    int i;
    double dx = q1[0] - q0[0];
    double dy = q1[1] - q0[1];
    double D = sqrt(dx * dx + dy * dy);
    double d = D / rho;
    if(rho <= 0.)
    {
        return EDUBBADRHO;
    }
    double theta = dbsmod2pi(atan2(dy, dx));
    double alpha = dbsmod2pi(q0[2] - theta);
    double beta  = dbsmod2pi(q1[2] - theta);
    for(i = 0; i < 3; i ++)
    {
        path->qi[i] = q0[i];
    }
    path->rho = rho;

    return dubins_init_normalised(alpha, beta, d, path);
}

int dubins_LSL(double alpha, double beta, double d, double* outputs)
{
    UNPACK_INPUTS(alpha, beta);
    double tmp0 = d+sa-sb;
    double p_squared = 2 + (d*d) -(2*c_ab) + (2*d*(sa - sb));
    if(p_squared < DUBINS_ZERO)
    {
        return EDUBNOPATH;
    }
    double tmp1 = atan2( (cb-ca), tmp0);
    if(fabs(cb-ca) < DUBINS_EPS)
        tmp1 = 0;
    double t = dbsmod2pi(-alpha + tmp1);
    double p = sqrt(MAX(p_squared, 0));
    double q = dbsmod2pi(beta - tmp1);
    PACK_OUTPUTS(outputs);

#if DEBUG
    assert(fabs(p*cos(alpha + t) - sa + sb - d) < DUBINS_EPS);
    assert(fabs(p*sin(alpha + t) + ca - cb) < DUBINS_EPS);
    assert(dbsmod2pi(alpha + t + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
#endif

    return EDUBOK;
}

int dubins_RSR( double alpha, double beta, double d, double* outputs )
{
    UNPACK_INPUTS(alpha, beta);
    double tmp0 = d-sa+sb;
    double p_squared = 2 + (d*d) -(2*c_ab) + (2*d*(sb-sa));
    if(p_squared < DUBINS_ZERO)
    {
        return EDUBNOPATH;
    }
    double tmp1 = atan2((ca-cb), tmp0);
    if(fabs(cb-ca) < DUBINS_EPS)
        tmp1 = 0;
    double t = dbsmod2pi( alpha - tmp1);
    double p = sqrt(MAX(p_squared, 0));
    double q = dbsmod2pi( -beta + tmp1);
    PACK_OUTPUTS(outputs);

#if DEBUG
    assert(fabs(p*cos(alpha - t) + sa - sb - d) < DUBINS_EPS);
    assert(fabs(p*sin(alpha - t) - ca + cb) < DUBINS_EPS);
    assert(dbsmod2pi(alpha - t - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
#endif

    return EDUBOK;
}

int dubins_LSR( double alpha, double beta, double d, double* outputs )
{
    UNPACK_INPUTS(alpha, beta);
    double p_squared = -2 + (d*d) + (2*c_ab) + (2*d*(sa+sb));
    if(p_squared < DUBINS_ZERO)
    {
        return EDUBNOPATH;
    }
    double p    = sqrt(MAX(p_squared, 0));
    double tmp2 = atan2((-ca-cb), (d+sa+sb)) - atan2(-2.0, p);
    if(fabs(-ca-cb) < DUBINS_EPS)
        tmp2 = -atan2(-2., p);
    double t    = dbsmod2pi(-alpha + tmp2);
    double q    = dbsmod2pi(-dbsmod2pi(beta) + tmp2);
    PACK_OUTPUTS(outputs);

#if DEBUG
    assert(fabs(p*cos(alpha + t) + 2. * sin(alpha + t) - sa - sb - d) < DUBINS_EPS);
    assert(fabs(p*sin(alpha + t) - 2. * cos(alpha + t) + ca + cb) < DUBINS_EPS);
    assert(dbsmod2pi(alpha + t - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
#endif
    return EDUBOK;
}

int dubins_RSL( double alpha, double beta, double d, double* outputs )
{
    UNPACK_INPUTS(alpha, beta);
    double p_squared = (d*d) -2 + (2*c_ab) - (2*d*(sa+sb));
    if(p_squared < DUBINS_ZERO)
    {
        return EDUBNOPATH;
    }
    double p    = sqrt(MAX(p_squared, 0));
    double tmp2 = atan2((ca+cb), (d-sa-sb)) - atan2(2.0, p);
    if(fabs(ca + cb) < DUBINS_EPS)
        tmp2 = 0 - atan2(2.0, p);
    double t    = dbsmod2pi(alpha - tmp2);
    double q    = dbsmod2pi(beta - tmp2);
    PACK_OUTPUTS(outputs);

#if DEBUG
    assert(fabs(p*cos(alpha - t) - 2. * sin(alpha - t) + sa + sb - d) < DUBINS_EPS);
    assert(fabs(p*sin(alpha - t) + 2. * cos(alpha - t) - ca - cb) < DUBINS_EPS);
    assert(dbsmod2pi(alpha - t + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
#endif
    return EDUBOK;
}

int dubins_RLR( double alpha, double beta, double d, double* outputs )
{
    UNPACK_INPUTS(alpha, beta);
    double tmp_rlr = (6. - d*d + 2*c_ab + 2*d*(sa-sb)) / 8.;
    if( fabs(tmp_rlr) > 1)
    {
        return EDUBNOPATH;
    }
    double tmp1 = atan2(ca-cb, d-sa+sb);
    if(fabs(ca - cb) < DUBINS_EPS)
        tmp1 = 0;
    double p = dbsmod2pi(2*M_PI - acos(tmp_rlr ));
    double t = dbsmod2pi(alpha - tmp1 + dbsmod2pi(p/2.));
    double q = dbsmod2pi(alpha - beta - t + dbsmod2pi(p));
    PACK_OUTPUTS(outputs);

#if DEBUG
    assert(fabs( 2.*sin(alpha - t + p) - 2. * sin(alpha - t) - d + sa - sb) < DUBINS_EPS);
    assert(fabs(-2.*cos(alpha - t + p) + 2. * cos(alpha - t) - ca + cb) < DUBINS_EPS);
    assert(dbsmod2pi(alpha - t + p - q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
#endif
    return EDUBOK;
}

int dubins_LRL( double alpha, double beta, double d, double* outputs )
{
    UNPACK_INPUTS(alpha, beta);
    double tmp_lrl = (6. - d*d + 2*c_ab + 2*d*(- sa + sb)) / 8.;
    if( fabs(tmp_lrl) > 1)
    {
        return EDUBNOPATH;
    }
    double tmp1 = atan2(ca-cb, d+sa-sb);
    if(fabs(ca-cb) < DUBINS_EPS)
        tmp1 = 0;
    double p = dbsmod2pi(2*M_PI - acos(tmp_lrl));
    double t = dbsmod2pi(-alpha - tmp1 + p/2.);
    double q = dbsmod2pi(dbsmod2pi(beta) - alpha -t + dbsmod2pi(p));
    PACK_OUTPUTS(outputs);

#if DEBUG
    assert(fabs(-2.*sin(alpha + t - p) + 2. * sin(alpha + t) - d - sa + sb) < DUBINS_EPS);
    assert(fabs( 2.*cos(alpha + t - p) - 2. * cos(alpha + t) + ca - cb) < DUBINS_EPS);
    assert(dbsmod2pi(alpha + t - p + q - beta + .5 * DUBINS_EPS) < DUBINS_EPS);
#endif
    return EDUBOK;
}

double dubins_path_length( DubinsPath* path )
{
    double length = 0.;
    length += path->param[0];
    length += path->param[1];
    length += path->param[2];
    length = length * path->rho;
    return length;
}

int dubins_path_type( DubinsPath* path )
{
    return path->type;
}

void dubins_segment( double t, double qi[3], double qt[3], int type)
{
    assert( type == L_SEG || type == S_SEG || type == R_SEG );

    if( type == L_SEG ) {
        qt[0] = qi[0] + sin(qi[2]+t) - sin(qi[2]);
        qt[1] = qi[1] - cos(qi[2]+t) + cos(qi[2]);
        qt[2] = qi[2] + t;
    }
    else if( type == R_SEG ) {
        qt[0] = qi[0] - sin(qi[2]-t) + sin(qi[2]);
        qt[1] = qi[1] + cos(qi[2]-t) - cos(qi[2]);
        qt[2] = qi[2] - t;
    }
    else if( type == S_SEG ) {
        qt[0] = qi[0] + cos(qi[2]) * t;
        qt[1] = qi[1] + sin(qi[2]) * t;
        qt[2] = qi[2];
    }
}

int dubins_path_sample( DubinsPath* path, double t, double q[3] )
{
    if( t < 0 || t >= dubins_path_length(path) ) {
        // error, parameter out of bounds
        return EDUBPARAM;
    }

    // tprime is the normalised variant of the parameter t
    double tprime = t / path->rho;

    // In order to take rho != 1 into account this function needs to be more complex
    // than it would be otherwise. The transformation is done in five stages.
    //
    // 1. translate the components of the initial configuration to the origin
    // 2. generate the target configuration
    // 3. transform the target configuration
    //      scale the target configuration
    //      translate the target configration back to the original starting point
    //      normalise the target configurations angular component

    // The translated initial configuration
    double qi[3] = { 0, 0, path->qi[2] };

    // Generate the target configuration
    const int* types = DIRDATA[path->type];
    double p1 = path->param[0];
    double p2 = path->param[1];
    double q1[3]; // end-of segment 1
    double q2[3]; // end-of segment 2
    dubins_segment(p1,      qi,    q1, types[0]);
    dubins_segment(p2,      q1,    q2, types[1]);
    if(tprime < p1)
    {
        dubins_segment( tprime, qi, q, types[0]);
    }
    else if(tprime < (p1+p2))
    {
        dubins_segment(tprime-p1, q1, q, types[1]);
    }
    else
    {
        dubins_segment(tprime-p1-p2, q2, q, types[2]);
    }

    // scale the target configuration, translate back to the original starting point
    q[0] = q[0] * path->rho + path->qi[0];
    q[1] = q[1] * path->rho + path->qi[1];
    q[2] = dbsmod2pi(q[2]);

    return 0;
}

int dubins_path_sample_many(DubinsPath* path, DubinsPathSamplingCallback cb, double stepSize, void* user_data)
{
    double x = 0.0;
    double length = dubins_path_length(path);
    while(x <  length)
    {
        double q[3];
        dubins_path_sample(path, x, q);
        int retcode = cb(q, x, user_data);
        if(retcode != 0)
        {
            return retcode;
        }
        x += stepSize;
    }
    return 0;
}

int dubins_path_endpoint( DubinsPath* path, double q[3] )
{
    // TODO - introduce a new constant rather than just using DUBINS_EPS
    return dubins_path_sample( path, dubins_path_length(path) - DUBINS_EPS, q );
}

int dubins_extract_subpath( DubinsPath* path, double t, DubinsPath* newpath )
{
    // calculate the true parameter
    double tprime = t / path->rho;

    // copy most of the data
    newpath->qi[0] = path->qi[0];
    newpath->qi[1] = path->qi[1];
    newpath->qi[2] = path->qi[2];
    newpath->rho   = path->rho;
    newpath->type  = path->type;

    // fix the parameters
    newpath->param[0] = fmin( path->param[0], tprime );
    newpath->param[1] = fmin( path->param[1], tprime - newpath->param[0]);
    newpath->param[2] = fmin( path->param[2], tprime - newpath->param[0] - newpath->param[1]);
    return 0;
}


#endif
