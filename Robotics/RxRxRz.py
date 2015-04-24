from Serial import *

R = Serial("RxRxRz", '', Matrix([['x','x','z'],['y','y','y']]).T)
mot  = Motor("mot1",'4')
mot2 = Motor("mot2",'5')


ph_ = R.ph_
po_ = Matrix([mot.p_,mot2.p_])
p_ = Matrix([ph_,po_])
M_ =    diag(R.Mh_sb_, mot.M_, mot2.M_)
v_ = Matrix([R.vh_sb_, mot.v_, mot2.v_])
g_ = Matrix([R.gh_sb_, mot.g_, mot2.g_])
phi_ = po_ - Matrix([R.vw_[0][0]+ph_[1],0,0, R.vw_[0][0]+symbols('beta')*ph_[1],0,0]).subs(R.StaticBal)
Ah_ = phi_.jacobian(ph_)
Ao_ = phi_.jacobian(po_)
C_ = Matrix([eye(3),simplify(-Ao_**-1 * Ah_)])
Mh_ = simplify( (C_.T * M_ * C_).subs(R.Jy[2], R.Jx[2]) )
vh_ = simplify(  (C_.T * ( M_ * C_.diff(t)*ph_ + v_ ) ).subs(R.Jy[2], R.Jx[2]))
gh_ = simplify(C_.T *g_)
Sol = solve( Matrix([Mh_[1,0]]) , [symbols('beta')] )
Mh_db_ = simplify(Mh_.subs(Sol))
vh_db_ = simplify(vh_.subs(Sol))
gh_db_ = simplify(gh_.subs(Sol))

pprint(p_)
pprint(M_)
pprint(v_)
pprint(g_)
pprint(phi_)
pprint(Ah_)
pprint(Ao_)
pprint(C_)
pprint(Mh_)
pprint(vh_)
pprint(gh_)
pprint(Sol)
pprint(Mh_db_)
pprint(vh_db_)
pprint(gh_db_)
