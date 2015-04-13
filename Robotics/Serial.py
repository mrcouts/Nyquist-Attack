from sympy import *
import math
init_printing(use_unicode=True)

t = symbols('t')

def Rotx(T):
    return Matrix([[1, 0, 0],[0, cos(T),-sin(T)],[0, sin(T), cos(T)]])

def Roty(T):
    return Matrix([[ cos(T), 0, sin(T)],[0, 1, 0],[-sin(T), 0, cos(T)]])

def Rotz(T):
    return Matrix([[cos(T),-sin(T),0],[sin(T), cos(T),0], [0, 0, 1]])

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
        print("Tah bem loko?")
        
def H2(axis_rot, axis_trans,T, l):
    if axis_trans == 'x':
        return H(axis_rot,T,l,0,0)
    elif axis_trans == 'y' :
        return H(axis_rot,T,0,l,0)
    elif axis_trans == 'z' :
        return H(axis_rot,T,0,0,l)
    else:
        print("Tah bem loko?")
        
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
    return v.row_insert(3,Matrix([[1]]))
    
def h2vec(hv):
    return hv[range(3),0]
    
def h2rot(M):
    return M[range(3),range(3)]
    
def nrows(M):
    return M.shape[0]
    
def ncols(M):
    return M.shape[1]
         
    
class Serial(object):
    """Serial robots dynamics."""
    def __init__(self, name, id, code):
        self.name = name
        self.id = id
        self.dof = nrows(code)
        self.code = code
        
        #Coordenadas Generalizadas
        self.qh_ = zeros(self.dof,1)
        for i in range(self.dof):
            if str(code[i,0]) != '0' :
                self.qh_[i] = Function('theta_' + str(i+1) )(t)
            else:
                self.qh_[i] = Function('d_' + str(i+1) )(t)
                
        self.qo_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            self.qo_[3*i+0] = Function('x_' + str(i+1) )(t)
            self.qo_[3*i+1] = Function('y_' + str(i+1) )(t)
            self.qo_[3*i+2] = Function('z_' + str(i+1) )(t)
            
        self.q_ = Matrix([self.qh_,self.qo_])
                
        #Velocidades Generalizadas
        self.ph_ = self.qh_.diff(t)
        
        pw_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            pw_[3*i+0] = Function('wx_' + str(i+1) )(t)
            pw_[3*i+1] = Function('wy_' + str(i+1) )(t)
            pw_[3*i+2] = Function('wz_' + str(i+1) )(t)
            
        pv_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            pv_[3*i+0] = Function('vx_' + str(i+1) )(t)
            pv_[3*i+1] = Function('vy_' + str(i+1) )(t)
            pv_[3*i+2] = Function('vz_' + str(i+1) )(t)
            
        po_ = Matrix([pw_,pv_])
        p_ = Matrix([self.ph_, po_])
        
        #Matrizes de transformacao homogenea relativas  
        vH_ = []
        vH_.append(H2(str(code[0,0]),str(code[0,1]),self.qh_[0],0))
        for i in range(1,self.dof):
            if code[i-1,0] != 0 :
                vH_.append(H2(str(code[i,0]),str(code[i,1]),self.qh_[i],symbols('l_'+str(i))))
            else:
                vH_.append(H2(str(code[i,0]),str(code[i,1]),self.qh_[i],self.qh_[i-1]))
                
        #Centros de massa nos S.C. das barras:
        vgb_ = []
        for i in range(self.dof):
            if code[i,0] != 0 :
                vgb_.append(MassCenter(str(code[i,1]),symbols('lg_'+str(i+1))))
            else:
                vgb_.append(MassCenter(str(code[i,1]),self.qh_[i] - symbols('l_'+str(i+1)) +  symbols('lg_'+str(i+1))))
            
        #Matrizes de transformacao homogenea absolutas
        vI_ = []
        vI_.append(vH_[0])
        for i in range(1,self.dof):
            vI_.append(0)
            vI_[i] = simplify(vI_[i-1]*vH_[i])
            
        #Centros de massa no S.C. N:
        vgn_ = []
        for i in range(self.dof):
            vgn_.append(0)
            vgn_[i] = h2vec(simplify(vI_[i]*vec2h(vgb_[i])))
            
        #Velocidades dos centros de massa:
        vv_ = []
        for i in range(self.dof):
            vv_.append(0)
            vv_[i] = simplify(vgn_[i].diff(t))
            
        #Matriz S:
        vS_ = []
        for i in range(self.dof):
            vS_.append(0)
            vS_[i] = simplify(h2rot(vI_[i]).T * h2rot(vI_[i]).diff(t))
            
        #Velocidades angulares das barras:
        vw_ = []
        for i in range(self.dof):
            vw_.append(0)
            vw_[i] = Matrix([vS_[i][2,1],vS_[i][0,2],vS_[i][1,0]])
            
        #Vetor po_ em funcao de ph_:
        w_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            w_[3*i+0] = vw_[i][0]
            w_[3*i+1] = vw_[i][1]
            w_[3*i+2] = vw_[i][2]
            
        v_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            v_[3*i+0] = vv_[i][0]
            v_[3*i+1] = vv_[i][1]
            v_[3*i+2] = vv_[i][2]
            
        _Po_ = Matrix([w_,v_])
        
        self.nonzerolist = []
        for i in range(6*self.dof):
            if _Po_[i] != 0:
                self.nonzerolist.append(i)
                
        _po_ = zeros(len(self.nonzerolist),1)
        for i in range(len(self.nonzerolist)):
            _po_[i] = _Po_[self.nonzerolist[i]]
                
        self.po_ = zeros(len(self.nonzerolist),1)
        for i in range(len(self.nonzerolist)):
            self.po_[i] = po_[self.nonzerolist[i]]
            
        self.p_ = Matrix([self.ph_,self.po_])
                
        
        #Matriz C:
        self.C_ = Matrix([eye(self.dof),_po_.jacobian(self.ph_)])
        
        #Matriz A:
        self.A_ = Matrix([_po_.jacobian(self.ph_).T, -eye(len(self.po_))]).T
        
        #Matriz b:
        self.b_ = simplify(-self.A_.diff(t)*self.p_)
        
        #energia de aceleracoes:
        self.s = 0
        for i in range(self.dof):
            self.s += ( symbols('m_' + str(i+1))*( pv_[3*i].diff(t)**2 + pv_[3*i+1].diff(t)**2 + pv_[3*i+2].diff(t)**2 ) + symbols('Jx'+str(i+1))*pw_[3*i].diff(t)**2 + symbols('Jy'+str(i+1))*pw_[3*i+1].diff(t)**2 + symbols('Jz'+str(i+1))*pw_[3*i+2].diff(t)**2 + 2*pw_[3*i].diff(t)*(symbols('Jz'+str(i+1))-symbols('Jy'+str(i+1)))*pw_[3*i+2]*pw_[3*i+1] + 2*pw_[3*i+1].diff(t)*(symbols('Jx'+str(i+1))-symbols('Jz'+str(i+1)))*pw_[3*i]*pw_[3*i+2] + + 2*pw_[3*i+2].diff(t)*(symbols('Jy'+str(i+1))-symbols('Jx'+str(i+1)))*pw_[3*i+1]*pw_[3*i]  )/2
            
        #energia potencial
        self.ep = 0
        for i in range(self.dof):
            self.ep += symbols('m_' + str(i+1))*symbols('g')*self.qo_[3*i+2]
    
        self.M_ = (Matrix([self.s]).jacobian(self.p_.diff(t))).jacobian(self.p_.diff(t))
        self.v_ = ( (Matrix([self.s]).jacobian(self.p_.diff(t))).T - (self.M_)*self.p_.diff(t) ).jacobian(self.p_)*self.p_ / 2
        g_ = Matrix([self.ph_,pv_]).jacobian(p_).T * Matrix([self.ep]).jacobian(self.q_).T
        
        self.nonzerolist2 = range(self.dof)
        for i in range(len(self.nonzerolist)):
            self.nonzerolist2.append(0)
            self.nonzerolist2[i+self.dof] = self.nonzerolist[i] + self.dof
        
        self.g_ = zeros(len(self.nonzerolist2),1)
        for i in range(len(self.nonzerolist2)):
            self.g_[i] = g_[self.nonzerolist2[i]]
    
    def description(self):
        print "Sou um robo %s, de %d graus de liberdade, com id = %d." % (self.name, self.dof, self.id)

RR = Serial("RR", 0, Matrix([['z','x','x'],['z','y','y']]).T)
pprint(RR.q_)
pprint(RR.p_)

RR.description()
 
pprint(RR.C_)
pprint(RR.A_)
pprint(RR.b_)
#pprint(RR.s)
pprint(RR.p_)
pprint(RR.M_)
pprint(RR.v_)
pprint(RR.g_)
print(RR.nonzerolist)
print(RR.nonzerolist2)
    
T = symbols('T', cls=Function)
t = symbols('t')
theta = symbols('theta', cls=Function)
x = symbols('x', cls=Function)
y = symbols('y', cls=Function)
z = symbols('z', cls=Function)