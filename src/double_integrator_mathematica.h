#ifndef __double_integrator_mathematica_h__
#define __double_integrator_mathematica_h__

#define mt1l1(g)    (-2*dx10*(g)*um - sqrt(2)*sqrt(pow(dx10,2)*pow((g)*um,2) - \
                    2*pow((g)*um,3)*x10))/(2.*pow((g)*um,2))

#define mt1l2(g)    (-2*dx10*(g)*um + sqrt(2)*sqrt(pow(dx10,2)*pow((g)*um,2) - \
                    2*pow((g)*um,3)*x10))/(2.*pow((g)*um,2))

#define mt2l(g)     (pow(dx10,2) - 2*(g)*um*x10)/(4.*pow((g)*um,2))

#define mt1r1(g)    (2*dx10*(g)*um - sqrt(2)*sqrt(pow(dx10,2)*pow((g)*um,2) + \
                    2*pow((g)*um,3)*x10))/(2.*pow((g)*um,2))

#define mt1r2(g)    (2*dx10*(g)*um + sqrt(2)*sqrt(pow(dx10,2)*pow((g)*um,2) + \
                    2*pow((g)*um,3)*x10))/(2.*pow((g)*um,2))

#define mt2r(g)     (pow(dx10,2) + 2*(g)*um*x10)/(4.*pow((g)*um,2))

#define is_right(x10, dx10, g)  ( ((x10 > dx10*dx10/2/(g)/um) || (dx10 > 0)) && ((x10 > -dx10*dx10/2/(g)/um) || (x10 > 0)) )
#define is_left(x10, dx10, g)   ( ((x10 < -dx10*dx10/2/(g)/um) || (dx10 < 0)) && ((x10 < dx10*dx10/2/(g)/um) || (x10 < 0)) )

#endif
