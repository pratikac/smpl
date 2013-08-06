#ifndef __double_integrator_mathematica_h__
#define __double_integrator_mathematica_h__

#define t1l1  (-2*dx10*um - sqrt(2)*sqrt(pow(dx10,2)*pow(um,2) - \
                2*pow(um,3)*x10))/(2.*pow(um,2))

#define t1l2  (-2*dx10*um + sqrt(2)*sqrt(pow(dx10,2)*pow(um,2) - \
                2*pow(um,3)*x10))/(2.*pow(um,2))

#define t2l   (pow(dx10,2) - 2*um*x10)/(4.*pow(um,2))

#define t1r1  (2*dx10*um - sqrt(2)*sqrt(pow(dx10,2)*pow(um,2) + \
                2*pow(um,3)*x10))/(2.*pow(um,2))

#define t1r2  (2*dx10*um + sqrt(2)*sqrt(pow(dx10,2)*pow(um,2) + \
                2*pow(um,3)*x10))/(2.*pow(um,2))

#define t2r   (pow(dx10,2) + 2*um*x10)/(4.*pow(um,2))

#endif
