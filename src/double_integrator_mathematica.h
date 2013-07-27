#ifndef __double_integrator_mathematica_h__
#define __double_integrator_mathematica_h__

#define mT1   (-(dx1f*u111) - sqrt(2)*sqrt(pow(dx1f,2)*pow(u111,2) + 2*pow(u111,3)*x1f))/pow(u111,2)

#define mt11  sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f))/(sqrt(2)*pow(u111,2))

#define mT2   (-(dx2f*u211) - sqrt(2)*sqrt(pow(dx2f,2)*pow(u211,2) + 2*pow(u211,3)*x2f))/pow(u211,2)

#define mt21  sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))/(sqrt(2)*pow(u211,2))

#define mt221 (pow(dx2f,3)*pow(u111,4) - pow(dx2f,2)*pow(u111,2)*u221*(dx1f*u111 + 3*sqrt(2)*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f))) - \
              dx2f*u111*u221*(u111*u221*(pow(dx1f,2) - 4*u111*x1f) + 2*sqrt(2)*dx1f*u221*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f)) + 4*pow(u111,3)*x2f) + \
              pow(u221,2)*(pow(dx1f,3)*u111*u221 + sqrt(2)*pow(dx1f,2)*u221*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f)) + 4*dx1f*pow(u111,2)*(u221*x1f - u111*x2f) + \
              4*sqrt(2)*u111*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f))*(u221*x1f + u111*x2f)))/ \
              (4.*pow(u111,2)*u221*(pow(dx2f,2)*pow(u111,2) + 2*dx1f*dx2f*u111*u221 - pow(u221,2)*(pow(dx1f,2) + 4*u111*x1f)))

#define mt222 (-(pow(dx2f,3)*pow(u111,4) + pow(dx2f,2)*pow(u111,2)*u221*(7*dx1f*u111 + 5*sqrt(2)*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f))) +  \
              dx2f*u111*u221*(3*u111*u221*(pow(dx1f,2) - 4*u111*x1f) + 6*sqrt(2)*dx1f*u221*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f)) + 4*pow(u111,3)*x2f) - \
              pow(u221,2)*(3*pow(dx1f,3)*u111*u221 + 3*sqrt(2)*pow(dx1f,2)*u221*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f)) - 4*dx1f*pow(u111,2)*(-3*u221*x1f + u111*x2f) + \
              4*sqrt(2)*u111*sqrt(pow(u111,2)*(pow(dx1f,2) + 2*u111*x1f))*(3*u221*x1f + u111*x2f)))/ \
              (4.*pow(u111,2)*u221*(pow(dx2f,2)*pow(u111,2) + 2*dx1f*dx2f*u111*u221 - pow(u221,2)*(pow(dx1f,2) + 4*u111*x1f))))

#define mt121 (-(2*u121*u211*(dx2f*u121*u211 + dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))) + \
              sqrt(2)*sqrt(u121*pow(u211,3)*(3*pow(dx2f,2)*pow(u121,2)*(u121 + u1221)*u211 + 2*dx2f*u121*(u121 + u1221)*(dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))) + \
              u211*(pow(dx1f,2)*(3*u121 - u1221)*pow(u211,2) + 4*u121*u211*((-u121 + u1221)*u211*x1f + u121*(u121 + u1221)*x2f) + \
              2*sqrt(2)*dx1f*u121*(u121 + u1221)*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))))))/(2.*u121*(u121 - u1221)*pow(u211,3)), \
              (-2*u121*u211*(dx2f*u121*u211 + dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))) + \
              sqrt(2)*sqrt(u121*pow(u211,3)*(3*pow(dx2f,2)*pow(u121,2)*(u121 + u1221)*u211 + 2*dx2f*u121*(u121 + u1221)*(dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))) +  \
              u211*(pow(dx1f,2)*(3*u121 - u1221)*pow(u211,2) + 4*u121*u211*((-u121 + u1221)*u211*x1f + u121*(u121 + u1221)*x2f) + \
              2*sqrt(2)*dx1f*u121*(u121 + u1221)*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))))))/(2.*u121*(u121 - u1221)*pow(u211,3)))

#define mt122 (-((3*u121 - u1221)*u211*(dx2f*u121*u211 + dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))) + \
              sqrt(2)*sqrt(u121*pow(u211,3)*(3*pow(dx2f,2)*pow(u121,2)*(u121 + u1221)*u211 + 2*dx2f*u121*(u121 + u1221)*(dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))) + \
              u211*(pow(dx1f,2)*(3*u121 - u1221)*pow(u211,2) + 4*u121*u211*((-u121 + u1221)*u211*x1f + u121*(u121 + u1221)*x2f) + \
              2*sqrt(2)*dx1f*u121*(u121 + u1221)*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))))))/(2.*u121*(u121 - u1221)*pow(u211,3)),\
              (-((3*u121 - u1221)*u211*(dx2f*u121*u211 + dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f)))) + \
              sqrt(2)*sqrt(u121*pow(u211,3)*(3*pow(dx2f,2)*pow(u121,2)*(u121 + u1221)*u211 + 2*dx2f*u121*(u121 + u1221)*(dx1f*pow(u211,2) + sqrt(2)*u121*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))) + \
              u211*(pow(dx1f,2)*(3*u121 - u1221)*pow(u211,2) + 4*u121*u211*((-u121 + u1221)*u211*x1f + u121*(u121 + u1221)*x2f) + \
              2*sqrt(2)*dx1f*u121*(u121 + u1221)*sqrt(pow(u211,2)*(pow(dx2f,2) + 2*u211*x2f))))))/(2.*u121*(u121 - u1221)*pow(u211,3)))

#endif
