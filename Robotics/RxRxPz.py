from Serial import *

R = Serial("RRP", '', Matrix([['x','x','0'],['y','y','x']]).T)
mot = Motor("mot",'4')

R.StaticBal = solve(R.gh_, R.lg)

R.Mh_sb_ = simplify(R.Mh_.subs(R.StaticBal))
R.vh_sb_ = simplify(R.vh_.subs(R.StaticBal))
R.gh_sb_ = simplify(R.gh_.subs(R.StaticBal))

ph_ = R.ph_
po_ = mot.p_
p_ = Matrix([ph_,po_])
M_ = diag(R.Mh_sb_, mot.M_)
v_ = Matrix([R.vh_sb_, mot.v_])
g_ = Matrix([R.gh_sb_, mot.g_])
phi_ = po_ - Matrix([R.vw_[0][0]+symbols('beta')*ph_[1],0,0,R.vv_[0][0],R.vv_[0][1],R.vv_[0][2]]).subs(R.StaticBal)
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
