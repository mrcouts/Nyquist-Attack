from Serial import *

R = Serial("Rx", '', Matrix([['x','x'],['y','y']]).T)

x_ = Matrix([R.qh_,R.ph_])
xn = lambda X:x_.subs([(x_[i],X[i]) for i in range(len(x_)-1,-1,-1)]).evalf()

m = Matrix([1.0,1.0])
l = Matrix([0.2,0.2])
lg = Matrix([0.1,0.1])
Jx = Matrix([0.00333333333333,0.00333333333333])
g = Matrix([9.8])

PI = Matrix([m,l,lg,Jx,g])
symPI = Matrix([Matrix(R.m),Matrix(R.l),Matrix(R.lg),Matrix(R.Jx), Matrix([symbols('g')]) ])
PIrep = [(symPI[i], PI[i]) for i in range(len(symPI))]

r = Matrix([sin(t),sin(t)])
Kp = 100.0*eye(len(R.qh_))
Kv = 20.0*eye(len(R.qh_))
u = ( R.vh_ + R.gh_ + R.Mh_*(r.diff(t,2) + Kv*(r.diff(t) - R.ph_ ) + Kp*(r - R.qh_ )  ) )

f = simplify(Matrix([R.ph_, R.Mh_**-1 * (u - R.vh_ - R.gh_ )]))
f2 = simplify(f.subs(PIrep).evalf())
f3 = Matrix([f2[i] for i in range(len(f2))])
fn = lambda t,X:f3.subs([(x_[i],X[i]) for i in range(len(x_)-1,-1,-1)]).subs(symbols('t'),t).evalf()

if(True):   
    from RK import *
    import time     
    t0 = 0
    tf = 10
    n = 100
    start_time = time.time()
    Y = TR('Euler','Euler').TRX(fn, t0, Matrix([0,0,0,0]), n, tf, tol=1e-5, nmax_gnr=50,nmax_gss=100)
    #Y = RK('RK5').RKX(fn, t0, Matrix([0,0,0,0]), n, tf)
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.linspace(t0, tf, n)
    y = x.copy()
    for i in np.arange(np.size(x)):
        y[i] = Y[0,i]
    plt.figure()
    plt.plot(x, y, 'r')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('title')
    plt.show()