from sympy import *
init_session(quiet=True)

u, um = symbols('u, um')
x0,dx0 = symbols('x0, dx0')
xf,dxf = symbols('xf, dxf')
t,ts,T = symbols('t, ts, T')
x, dx = symbols('x, dx')
xs, dxs = symbols('xs, dxs')

if 1:
    u = 1
    x0,dx0 = -10,10
    xf,dxf = 0,0

dx = -1*u*ts + dx0
x = -0.5*u*ts**2 + dx0*ts + x0

e1 = x - xf - 1/u*(0.5*dx**2 - 0.5*dxf**2)
ans = solve(e1, ts)
ts_ans = [simplify(a) for a in ans]

# get T
dxs_ans = [simplify(dx.subs(ts, a)) for a  in ts_ans]
#xs_ans = [simplify(x.subs(ts, a)) for a  in ts_ans]

T_ans = [simplify(-(dxs_ans[i] - dxf)/u + ts_ans[i]) for i in xrange(2)]
