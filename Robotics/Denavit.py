from sympy import *
#import math
#from operator import concat
init_printing(use_unicode=True)

###############################################################################

t = symbols('t')

Rotx = lambda T: Matrix([[1, 0, 0],[0, cos(T),-sin(T)],[0, sin(T), cos(T)]])
Roty = lambda T: Matrix([[cos(T), 0, sin(T)],[0, 1, 0],[-sin(T), 0, cos(T)]])
Rotz = lambda T: Matrix([[cos(T), -sin(T), 0],[sin(T), cos(T), 0], [0, 0, 1]])

Hx = lambda T,dx,dy,dz: Rotx(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
Hy = lambda T,dx,dy,dz: Roty(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
Hz = lambda T,dx,dy,dz: Rotz(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def H(axis, T, dx, dy, dz):
    if axis == 'x' :
        return Hx(T,dx,dy,dz)
    elif axis == 'y' :
        return Hy(T,dx,dy,dz)
    elif axis == 'z' :
        return Hz(T,dx,dy,dz)
    elif axis == '0':
        return Hz(0,dx,dy,dz)
    else:
        raise ValueError("Soh pode ser string x, y, z ou 0!")
        
###############################################################################

def H_d(mat):
    _H_d = lambda a, alpha, d, theta: H('z',theta,0,0,0)*H('0',0,0,0,d)*H('0',0,a,0,0)*H('x',alpha,0,0,0)
    return [_H_d(*mat[i,0:4]) for i in xrange(mat.rows)]

class Serial(object):
    """Serial robots dynamics."""
    def __init__(self, name, ID, dof, fDH_, g_dir_, l_, lg_, m_, I__, c_, gamma_, n_):
        
        z = list(zeros(1,dof))
        if fDH_(z,z,z).rows != dof:
            raise ValueError("DH_ tem que ter dof colunas!")
                
        #Identificadores:
        self.name = name
        self.ID = ID
        Id = '' if ID == '' else '_' + str(ID)
        
        #Numero de graus de liberdade do sistema:
        self.dof = dof
        dof_ = xrange(dof)
        
        #Parametros geometricos:
        self.l_ = l_
        self.lg_ = lg_
        
        #Parametros de inercia:
        self.m_  = m_
        self.I__ = I__
        
        #Parametros de atrito:
        self.c_  = c_
        self.gamma_ = gamma_
        self.n_ = n_
        
        #Direcao da gravidade
        self.g_dir_ = g_dir_/g_dir_.norm()
        
        #Coordenadas Generalizadas
        self.q_ = Matrix([Function('theta_'+str(i+1)+Id)(t) if str(fDH_(z,z,z)[i,7]) == 'R' else 
        	               Function('d_'    +str(i+1)+Id)(t) for i in dof_ ])   
                        
        #Quasi-velocidades independentes
        self.dq_ = self.q_.diff(t)
                        
        #Matriz dos Parametros de Denavit-Hartemberg
        self.DH_ = fDH_(self.q_, l_, lg_)
                        
        #Componentes de velocidades e velocidades angulares de cada corpo rigido do sistema
        self.vel__ = [Matrix([Function('vx'+str(i+1)+Id)(t),Function('vy'+str(i+1)+Id)(t),Function('vz'+str(i+1)+Id)(t)]) for i in dof_]
        self.w__   = [Matrix([Function('wx'+str(i+1)+Id)(t),Function('wy'+str(i+1)+Id)(t),Function('wz'+str(i+1)+Id)(t)]) for i in dof_]
        self.w_ = Matrix([self.w__[i] for i in dof_])
        
        #Quasi-velocidades        
        p_ = Matrix([self.dq_, Matrix([ Matrix([self.vel__[i],self.w__[i]]) for i in dof_]) ])
        
        #-----------------------------Cinematica------------------------------#
        
        #Matrizes de transformacao homogenea relativas (B_i | B_{i+1} )
        self.Hr__ = H_d(self.DH_[:,0:4])
        
        #Matrizes de transformacao homogenea absolutas (N | B_i)
        self.H__ = []
        self.H__.append(self.Hr__[0])
        for i in range(1,dof):
            self.H__.append(simplify(self.H__[i-1]*self.Hr__[i]))
        
        self.z__ = [Matrix([0,0,1])] + [self.H__[i][0:3,2] for i in dof_] #[k]_N
        self.o__ = [Matrix([0,0,0])] + [self.H__[i][0:3,3] for i in dof_] #[o]_N
        self.og__= [simplify(Matrix([(self.H__[i]*Matrix([self.DH_[i,4:7].T,[1]]))[0:3]]).T) for i in dof_] #[g]_N
        
        #Jacobianos dos centros de massa
        self.Jv__ = [simplify(Matrix([self.z__[i].cross(self.og__[j]- self.o__[i]).T if str(self.DH_[i,7]) == 'R' and i <= j else 
                                     (self.z__[i].T if i <= j else 
                                      zeros(1,3) ) for i in dof_]).T ) for j in dof_]
                                          
        self.Jw__ = [simplify( self.H__[j][0:3,0:3].T*Matrix([self.z__[i].T if str(self.DH_[i,7]) == 'R' and i <= j else
                                                              zeros(1,3) for i in dof_]).T ) for j in dof_]
        
        self.J_ = Matrix( [Matrix([self.Jv__[i],self.Jw__[i]]) for i in dof_] )
        self.Jw_ = Matrix([self.Jw__[i] for i in dof_])
        
        #Jacobianos do efetuador
        self.Jv_n_ = simplify(Matrix( [self.z__[i].cross(self.o__[dof] - self.o__[i]).T if str(self.DH_[i,7]) == 'R' else
                                       self.z__[i].T for i in dof_ ]).T )
                                           
        self.Jw_n_ = simplify( self.H__[dof-1][0:3,0:3].T*Matrix( [self.z__[i].T if str(self.DH_[i,7]) == 'R' else 
                                                                   zeros(1,3) for i in dof_ ]).T )
                                   
        self.J_n_ = Matrix([self.Jv_n_,self.Jw_n_])
        
        #Dinamica
        C_ = Matrix([eye(dof),self.J_])

        non_null_p_index = []
        null_p_index = []
        for i in xrange(C_.rows):
        	if C_[i,:] != zeros(1,dof):
        		non_null_p_index.append(i)
        	else:
        		null_p_index.append(i)
          
        null_p_ = p_.extract(null_p_index,[0])
        replace = [(null_p_[i],0) for i in xrange(null_p_.rows)]
          
        self.C_ = C_.extract(non_null_p_index, range(C_.cols))
        self.p_ = p_.extract(non_null_p_index, range(p_.cols))

        self.A_ = Matrix([self.C_[dof:,:].T,-eye(self.C_.rows-dof).T]).T
        self.b_ = simplify( -self.A_.diff(t)*self.p_)
        
        self.M__ = [diag(eye(3)*m_[i], I__[i]) for i in dof_]
        self.v__ = [simplify(Matrix([zeros(3,1),self.w__[i].cross(I__[i]*self.w__[i])])) for i in dof_]
        self.g__ = [Matrix([-m_[i]*symbols('g')*self.g_dir_,zeros(3,1)]) for i in dof_]
        
        M_ = diag(zeros(dof),*self.M__)
        v_ = Matrix([zeros(dof,1), Matrix([self.v__[i] for i in dof_]) ]).subs(replace)
        g_ = Matrix([zeros(dof,1), Matrix([self.g__[i] for i in dof_]) ])
        f_ = Matrix([Matrix([self.c_[i]*self.dq_[i] + self.gamma_[i]*tanh(self.n_[i]*self.dq_[i]) for i in dof_]), zeros(6*dof,1) ])
        
        self.M_ = M_.extract(non_null_p_index, non_null_p_index)
        self.v_ = v_.extract(non_null_p_index, range(v_.cols))
        self.g_ = g_.extract(non_null_p_index, range(g_.cols))
        self.f_ = f_.extract(non_null_p_index, range(f_.cols))
        
        v_aux_ = simplify(self.v_.subs([(self.w_[i],(self.Jw_[i,:]*self.dq_)[0] ) for i in xrange(3*dof)]))
        
        self.Mh_ = simplify(self.C_.T*self.M_*self.C_)
        self.vh_ = simplify(self.C_.T*( v_aux_ + self.M_*self.C_.diff(t)*self.dq_ ))
        self.gh_ = simplify(self.C_.T*self.g_)
        self.fh_ = simplify(self.C_.T*self.f_)