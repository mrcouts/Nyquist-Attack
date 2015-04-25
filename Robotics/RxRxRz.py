from Serial import *

R = Serial("RxRxRz", '', Matrix([['x','x','z'],['y','y','y']]).T)
n_mot = 2
mot = [Motor("mot"+str(i+1),str(R.dof+i+1)) for i in range(n_mot)]

ph_ = R.ph_
po_ = Matrix([mot[i].p_ for i in range(n_mot)])
p_ = Matrix([ph_,po_])

M_ = R.Mh_sb_
for i in range(n_mot):
	M_ = diag(M_, mot[i].M_)
v_ = Matrix([R.vh_sb_] + [mot[i].v_ for i in range(n_mot)])
g_ = Matrix([R.gh_sb_] + [mot[i].g_ for i in range(n_mot)])

phi_ = po_ - Matrix([R.vw_[0][0]+ph_[1],0,0, R.vw_[0][0]+symbols('beta')*ph_[1],0,0]).subs(R.StaticBal)

Ah_ = phi_.jacobian(ph_)
Ao_ = phi_.jacobian(po_)
C_ = Matrix([eye(R.dof),simplify(-Ao_**-1 * Ah_)])
Mh_ = simplify( (C_.T * M_ * C_) )
vh_ = simplify(  (C_.T * ( M_ * C_.diff(t)*ph_ + v_ ) ))
gh_ = simplify(C_.T *g_)

Sol = solve( Matrix([Mh_[i,j] for i in range(1,R.dof) for j in range(i)]) , [symbols('beta')] )
Mh_db_ = simplify(Mh_.subs(Sol).subs(R.Jy[2], R.Jx[2]))
vh_db_ = simplify(vh_.subs(Sol).subs(R.Jy[2], R.Jx[2]))
gh_db_ = simplify(gh_.subs(Sol).subs(R.Jy[2], R.Jx[2]))

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
