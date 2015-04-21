from sympy import *
import math
from operator import concat
init_printing(use_unicode=True)

t = symbols('t')

def Rotx(T):
    return Matrix([[1, 0, 0],[0, cos(T),-sin(T)],[0, sin(T), cos(T)]])

def Roty(T):
    return Matrix([[cos(T), 0, sin(T)],[0, 1, 0],[-sin(T), 0, cos(T)]])

def Rotz(T):
    return Matrix([[cos(T), -sin(T), 0],[sin(T), cos(T), 0], [0, 0, 1]])

def Hx(T,dx,dy,dz):
    return Rotx(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def Hy(T,dx,dy,dz):
    return Roty(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def Hz(T,dx,dy,dz):
    return Rotz(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
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
        
def H2(axis_rot, axis_trans,T, l):
    if axis_trans == 'x':
        return H(axis_rot,T,l,0,0)
    elif axis_trans == 'y' :
        return H(axis_rot,T,0,l,0)
    elif axis_trans == 'z' :
        return H(axis_rot,T,0,0,l)
    else:
        raise ValueError("Soh pode ser string x, y, z!")
        
def MassCenter(axis, l):
    if axis == 'x':
        return Matrix([l,0,0])
    elif axis == 'y':
        return Matrix([0,l,0])
    elif axis == 'z' :
        return Matrix([0,0,l])
    else:
        print("Tah bem loko?")        
    
def vec2h(v):
    return Matrix([v,Matrix([1])])
    
def h2vec(hv):
    return hv[0:3,:]
    
def h2rot(M):
    return M[0:3,0:3]
    
def nrows(M):
    return M.shape[0]
    
def ncols(M):
    return M.shape[1]

def Jacobian(M,v):
    n = nrows(v)
    J = [M.diff(v[i]) for i in range(n)]
    return J

def JacDotVec(J, v):
    n = nrows(J[0])
    m = ncols(J[0])
    Sum = zeros(n,m)
    for i in range(len(J)):
        Sum += J[i]*v[i]
    return Sum         
    
class Serial(object):
    """Serial robots dynamics."""
    def __init__(self, name, ID, code):
        self.name = name
        self.ID = ID
        self.dof = nrows(code)
        self.code = code

        Id = '' if ID == '' else '_' + str(ID)

        self.m  = [symbols('m_' + str(i+1) + Id) for i in range(self.dof)]
        self.Jx = [symbols('Jx' + str(i+1) + Id) for i in range(self.dof)]
        self.Jy = [symbols('Jy' + str(i+1) + Id) for i in range(self.dof)]
        self.Jz = [symbols('Jz' + str(i+1) + Id) for i in range(self.dof)]
        self.Jx = [symbols('Jx' + str(i+1) + Id) for i in range(self.dof)]
        self.l  = [symbols('l_' + str(i+1) + Id) for i in range(self.dof)]
        self.lg = [symbols('lg_'+ str(i+1) + Id) for i in range(self.dof)]

        self.x  = [Function('x_' + str(i+1) + Id )(t) for i in range(self.dof)]
        self.y  = [Function('y_' + str(i+1) + Id )(t) for i in range(self.dof)]
        self.z  = [Function('z_' + str(i+1) + Id )(t) for i in range(self.dof)]

        self.wx = [Function('wx_' + str(i+1) + Id )(t) for i in range(self.dof)]
        self.wy = [Function('wy_' + str(i+1) + Id )(t) for i in range(self.dof)]
        self.wz = [Function('wz_' + str(i+1) + Id )(t) for i in range(self.dof)]

        self.vx = [Function('vx_' + str(i+1) + Id )(t) for i in range(self.dof)]
        self.vy = [Function('vy_' + str(i+1) + Id )(t) for i in range(self.dof)]
        self.vz = [Function('vz_' + str(i+1) + Id )(t) for i in range(self.dof)]

        #Verificando se o mecanismo tem juntas prismaticas
        self.pflag = any(str(code[i,0]) == "0" for i in range(self.dof))
        
        #Coordenadas Generalizadas
        self.qh_ = Matrix([Function('theta_'+str(i+1)+Id)(t) if str(code[i,0]) != '0' else 
        	               Function('d_'    +str(i+1)+Id)(t) for i in range(self.dof) ])
                
        self.qo_ = Matrix(reduce(concat, [ [self.x[i], self.y[i], self.z[i]] for i in range(self.dof) ]) )
        self.q_ = Matrix([self.qh_,self.qo_])
                
        #VelocIDades Generalizadas
        self.ph_ = self.qh_.diff(t)
        
        pw_ = Matrix(reduce(concat, [ [self.wx[i], self.wy[i], self.wz[i]] for i in range(self.dof) ]) )
        pv_ = Matrix(reduce(concat, [ [self.vx[i], self.vy[i], self.vz[i]] for i in range(self.dof) ]) )
        po_ = Matrix([pw_,pv_])
        p_ = Matrix([self.ph_, po_])
        
        #Matrizes de transformacao homogenea relativas  
        vH_ = ([H2(str(code[0,0]),str(code[0,1]),self.qh_[0],0)] + 
               [H2(str(code[i,0]),str(code[i,1]),self.qh_[i],self.l[i-1]) if str(code[i-1,0]) != '0' else 
                H2(str(code[i,0]),str(code[i,1]),self.qh_[i],self.qh_[i-1]) for i in range(1,self.dof)])
                
        #Centros de massa nos S.C. das barras:
        vgb_ = [MassCenter(str(code[i,1]),self.lg[i]) if str(code[i,0]) != '0' else 
                MassCenter(str(code[i,1]),self.qh_[i] - self.l[i] +  self.lg[i]) 
                for i in range(self.dof)]

        #Posicao do efetuador no S.C. da ultima barra:
        self.xb_ = (MassCenter(str(code[self.dof-1,1]),self.l[self.dof-1]) if str(code[self.dof-1,0]) != '0' else 
        	        MassCenter(str(code[i,1]),self.qh_[i]))
            
        #Matrizes de transformacao homogenea absolutas
        vI_ = []
        vI_.append(vH_[0])
        for i in range(1,self.dof):
            vI_.append(simplify(vI_[i-1]*vH_[i]))

        #Orientacao do efetuador
        self.R_ = h2rot(vI_[self.dof-1])
            
        #Centros de massa no S.C. N:
        vgn_ = [h2vec(simplify(vI_[i]*vec2h(vgb_[i]))) for i in range(self.dof)]

        #Posicao do efetuador no S.C. N:
        self.xn_ = h2vec(simplify(vI_[self.dof-1]*vec2h(self.xb_)))
            
        #Velocidades dos centros de massa:
        self.vv_ = [simplify(vgn_[i].diff(t)) for i in range(self.dof)]
            
        #Matriz S:
        vS_ = [simplify(h2rot(vI_[i]).T * h2rot(vI_[i]).diff(t)) for i in range(self.dof)]
            
        #Velocidades angulares das barras:
        self.vw_ = [Matrix([vS_[i][2,1],vS_[i][0,2],vS_[i][1,0]]) for i in range(self.dof)]
            
        #Vetor po_ em funcao de ph_:
        w_ = self.vw_[0]
        for i in range(1,self.dof):
        	w_ = Matrix([w_,self.vw_[i]])

        v_ = self.vv_[0]
        for i in range(1,self.dof):
        	v_ = Matrix([v_,self.vv_[i]])
            
        _Po_ = Matrix([w_,v_])

        #Indices dos elementos de _Po_ nao nulos      
        nonzerolist = [i for i in range(6*self.dof) if _Po_[i] != 0] 
        
        #Indices dos elementos de _p_ nao nulos   
        nonzerolist2 = range(self.dof) + [i + self.dof for i in nonzerolist]

        _po_ = _Po_[nonzerolist, :]
        self.po_ = po_[nonzerolist, :]
        self.p_ = Matrix([self.ph_,self.po_])
                
        #Matriz C:
        self.C_ = Matrix([eye(self.dof),_po_.jacobian(self.ph_)])
        
        #Matriz A:
        self.A_ = Matrix([_po_.jacobian(self.ph_).T, -eye(len(self.po_))]).T
        
        #Matriz b:
        self.b_ = simplify(-self.A_.diff(t)*self.p_)
        
        #energia de aceleracoes:
        s = sum([
           	(self.m[i]*( self.vx[i].diff(t)**2 + self.vy[i].diff(t)**2 + self.vz[i].diff(t)**2 ) + 
           	 self.Jx[i]*self.wx[i].diff(t)**2 + 
           	 self.Jy[i]*self.wy[i].diff(t)**2 + 
           	 self.Jz[i]*self.wz[i].diff(t)**2 + 
           	 2*self.wx[i].diff(t)*(self.Jz[i]-self.Jy[i])*self.wz[i]*self.wy[i] + 
           	 2*self.wy[i].diff(t)*(self.Jx[i]-self.Jz[i])*self.wx[i]*self.wz[i] + 
           	 2*self.wz[i].diff(t)*(self.Jy[i]-self.Jx[i])*self.wy[i]*self.wx[i] )/2 
           	 for i in range(self.dof) ])
            
        #energia potencial
        ep = sum([self.m[i]*symbols('g')*self.z[i] for i in range(self.dof)])

        #energia cinetica
        self.ec = sum([
           	(self.m[i]*( self.vv_[i][0]**2 + self.vv_[i][1]**2 + self.vv_[i][2]**2 ) + 
           	 self.Jx[i]*self.vw_[i][0]**2 + 
           	 self.Jy[i]*self.vw_[i][1]**2 + 
           	 self.Jz[i]*self.vw_[i][2]**2 )/2 
           	 for i in range(self.dof) ])
    
        #Matrizes da dinamica
        self.M_ = (Matrix([s]).jacobian(self.p_.diff(t))).jacobian(self.p_.diff(t))
        self.v_ = ( (Matrix([s]).jacobian(self.p_.diff(t))).T - (self.M_)*self.p_.diff(t) ).jacobian(self.p_)*self.p_ / 2
        g_ = Matrix([self.ph_,pv_]).jacobian(p_).T * Matrix([ep]).jacobian(self.q_).T
        self.g_ = g_[nonzerolist2, :]

        v_aux_ = simplify(JacDotVec(Jacobian(self.v_.jacobian(self.p_),self.p_),self.C_*self.ph_)*(self.C_*self.ph_)/2)
        
        #Matrizes da dinamica hashtag
        self.Mh_ = simplify(self.C_.T * self.M_ * self.C_)
        self.vh_ = simplify(self.C_.T * (v_aux_ + self.M_ * self.C_.diff(t) * self.ph_ ) )
        self.gh_ = simplify(self.C_.T * self.g_)

        self.Mh2_ = simplify( (Matrix([self.ec]).jacobian(self.ph_).T).jacobian(self.ph_) )
        self.vh2_ = simplify( (Matrix([self.ec]).jacobian(self.ph_).T).diff(t) - (Matrix([self.ec]).jacobian(self.qh_).T) -  self.Mh2_*self.ph_.diff(t)  )

        
        #Balanceamento estatico
        if self.pflag == False:
            self.StaticBal = solve(self.gh_, self.lg)
            
            self.gh_sb_ = simplify(self.gh_.subs(self.StaticBal))
            self.Mh_sb_ = simplify(self.Mh_.subs(self.StaticBal))
            self.vh_sb_ = simplify(self.vh_.subs(self.StaticBal))
            self.ec_sb  = simplify(self.ec.subs(self.StaticBal))
        else:
            self.StaticBal = None
            self.gh_sb_ = None
            self.Mh_sb_ = None
            self.vh_sb_ = None
    
    def description(self):
        print "Sou um robo %s, de %d graus de liberdade, com ID = %s." % (self.name, self.dof, str(self.ID))
        print "qh_ = "
        pprint(self.qh_)
        print " "
        print "p_ = "
        pprint(self.p_)
        print " "
        print "C_ = "
        pprint(self.C_)
        print " "
        print "A_ = "
        pprint(self.A_)
        print " "
        print "b_ = "
        pprint(self.b_)
        print " "
        print "M_ = "
        pprint(self.M_)
        print " "
        print "v_ = "
        pprint(self.v_)
        print " "
        print "g_ = "
        pprint(self.g_)
        print " "
        print "Mh_ = "
        pprint(self.Mh_)
        print " "
        print "vh_ = "
        pprint(self.vh_)
        print " "
        print "gh_ = "
        pprint(self.gh_)
        print " "
        print "Solucao do Balanceamento Estatico:"
        pprint(self.StaticBal)
        print " "
        print "Mh_sb_ = "
        pprint(self.Mh_sb_)
        print " "
        print "vh_sb_ = "
        pprint(self.vh_sb_)
        print " "
        print "gh_sb_ = "
        pprint(self.gh_sb_)
        #for i in range(self.dof):
        #    pprint(self.vv_[i])
        #for i in range(self.dof):
        #    pprint(self.vw_[i])
        
class Motor(object):
    """Motor dynamics."""
    dof = 6 
    def __init__(self, name, ID):
        self.name = name
        self.ID = ID
        Id = '' if ID == '' else '_' + str(ID)

        self.Jx = symbols('Jx'+Id)
        self.Jy = symbols('Jy'+Id)
        self.Jz = symbols('Jz'+Id)

        self.wx = Function('wx'+Id)(t)
        self.wy = Function('wy'+Id)(t)
        self.wz = Function('wz'+Id)(t)
        self.vx = Function('vx'+Id)(t)
        self.vy = Function('vy'+Id)(t)
        self.vz = Function('vz'+Id)(t)

        self.theta = Function('theta'+Id)(t)

        #self.qh_ = Matrix([self.theta])
        #self.ph_ = self.qh_.diff(t)
        #self.po_ = Matrix([self.wx , self.wy, self.wz, self.vx, self.vy, self.vz])
        #self.p_  = Matrix([self.ph_,self.po_])
        self.p_ = Matrix([self.wx , self.wy, self.wz, self.vx, self.vy, self.vz])

        self.M_ = diag(self.Jx, self.Jy, self.Jz, 0, 0, 0)
        self.v_ = zeros(6,1)
        self.g_ = zeros(6,1)
    def description(self):
        print "Sou um motor chamado %s, com ID = %s." % (self.name, str(self.ID))
        #print "qh_ = "
        #pprint(self.qh_)
        #print " "
        print "p_ = "
        pprint(self.p_)
        print " "
        print "M_ = "
        pprint(self.M_)
        print " "
        print "v_ = "
        pprint(self.v_)
        print " "
        print "g_ = "
        pprint(self.g_)
        print " "