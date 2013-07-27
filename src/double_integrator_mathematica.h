#ifndef __double_integrator_mathematica_h__
#define __double_integrator_mathematica_h__

#define mT1   (-(dx1f*u111) + sqrt(2)*sqrt(pow(dx1f,2)*pow(u111,2) + 2*pow(u111,3)*x1f))/pow(u111,2)

#define mt11  sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f))/(sqrt(2)*pow(u111,2))

#define mT2   (-(dx2f*u211) + sqrt(2)*sqrt(pow(dx2f,2)*pow(u211,2) + 2*pow(u211,3)*x2f))/pow(u211,2)

#define mt21  sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))/(sqrt(2)*pow(u211,2))

#define mt221 (-pow(dx2f,2) - 2*dx2f*mT1*u221 + pow(mT1,2)*pow(u221,2) + 4*u221*x2f)/(4.*u221*(-dx2f + mT1*u221))

#define mt222 (-pow(dx2f,2) + 6*dx2f*mT1*u221 - 3*pow(mT1,2)*pow(u221,2) - 4*u221*x2f)/(4.*u221*(dx2f - mT1*u221))

#define mt121 (-pow(dx1f,2) - 2*dx1f*mT2*u121 + pow(mT2,2)*pow(u121,2) + 4*u121*x1f)/(4.*u121*(-dx1f + mT2*u121))

#define mt122 (-pow(dx1f,2) + 6*dx1f*mT2*u121 - 3*pow(mT2,2)*pow(u121,2) - 4*u121*x1f)/(4.*u121*(dx1f - mT2*u121))

#endif
