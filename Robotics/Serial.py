from sympy import *
import math
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
    return Matrix([v,Matrix([1])])
    
def h2vec(hv):
    return hv[0:3,:]
    
def h2rot(M):
    return M[0:3,0:3]
    
def nrows(M):
    return M.shape[0]
    
def ncols(M):
    return M.shape[1]
         
    
class Serial(object):
    """Serial robots dynamics."""
    def __init__(self, name, ID, code):
        self.name = name
        self.ID = ID
        self.dof = nrows(code)
        self.code = code

        #Verificando se o mecanismo tem juntas prismaticas
        self.pflag = False
        for i in range(self.dof):
            self.pflag |= (str(code[i,0]) == '0')

        if ID == '':
            Id = ''
        else:
            Id = '_' + str(ID)
        
        #Coordenadas Generalizadas
        self.qh_ = zeros(self.dof,1)
        for i in range(self.dof):
            if str(code[i,0]) != '0' :
                self.qh_[i] = Function('theta_' + str(i+1) + Id )(t)
            else:
                self.qh_[i] = Function('d_' + str(i+1) + Id )(t)
                
        self.qo_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            self.qo_[3*i+0] = Function('x_' + str(i+1) + Id )(t)
            self.qo_[3*i+1] = Function('y_' + str(i+1) + Id )(t)
            self.qo_[3*i+2] = Function('z_' + str(i+1) + Id )(t)
            
        self.q_ = Matrix([self.qh_,self.qo_])
                
        #VelocIDades Generalizadas
        self.ph_ = self.qh_.diff(t)
        
        pw_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            pw_[3*i+0] = Function('wx_' + str(i+1) + Id )(t)
            pw_[3*i+1] = Function('wy_' + str(i+1) + Id )(t)
            pw_[3*i+2] = Function('wz_' + str(i+1) + Id )(t)
            
        pv_ = zeros(3*self.dof,1)
        for i in range(self.dof):
            pv_[3*i+0] = Function('vx_' + str(i+1) + Id )(t)
            pv_[3*i+1] = Function('vy_' + str(i+1) + Id )(t)
            pv_[3*i+2] = Function('vz_' + str(i+1) + Id )(t)
            
        po_ = Matrix([pw_,pv_])
        p_ = Matrix([self.ph_, po_])
        
        #Matrizes de transformacao homogenea relativas  
        vH_ = []
        vH_.append(H2(str(code[0,0]),str(code[0,1]),self.qh_[0],0))
        for i in range(1,self.dof):
            if str(code[i-1,0]) != '0' :
                vH_.append(H2(str(code[i,0]),str(code[i,1]),self.qh_[i],symbols('l_'+str(i) + Id)))
            else:
                vH_.append(H2(str(code[i,0]),str(code[i,1]),self.qh_[i],self.qh_[i-1]))
                
        #Centros de massa nos S.C. das barras:
        vgb_ = []
        for i in range(self.dof):
            if str(code[i,0]) != '0':
                vgb_.append(MassCenter(str(code[i,1]),symbols('lg_'+str(i+1) + Id)))
            else:
                vgb_.append(MassCenter(str(code[i,1]),self.qh_[i] - symbols('l_'+str(i+1) + Id) +  symbols('lg_'+str(i+1) + Id)))

        #Posicao do efetuador no S.C. da ultima barra:
        if str(code[self.dof-1,0]) != '0' :
            self.xb_ = MassCenter(str(code[self.dof-1,1]),symbols('l_'+str(self.dof) + Id))
        else:
            self.xb_ = MassCenter(str(code[i,1]),self.qh_[i])
            
        #Matrizes de transformacao homogenea absolutas
        vI_ = []
        vI_.append(vH_[0])
        for i in range(1,self.dof):
            vI_.append(simplify(vI_[i-1]*vH_[i]))

        #Orientacao do efetuador
        self.R_ = h2rot(vI_[self.dof-1])
            
        #Centros de massa no S.C. N:
        vgn_ = []
        for i in range(self.dof):
            vgn_.append(h2vec(simplify(vI_[i]*vec2h(vgb_[i]))))

        #Posicao do efetuador no S.C. N:
        self.xn_ = h2vec(simplify(vI_[self.dof-1]*vec2h(self.xb_)))
            
        #Velocidades dos centros de massa:
        self.vv_ = []
        for i in range(self.dof):
            self.vv_.append(simplify(vgn_[i].diff(t)))
            
        #Matriz S:
        vS_ = []
        for i in range(self.dof):
            vS_.append(simplify(h2rot(vI_[i]).T * h2rot(vI_[i]).diff(t)))
            
        #Velocidades angulares das barras:
        self.vw_ = []
        for i in range(self.dof):
            self.vw_.append(Matrix([vS_[i][2,1],vS_[i][0,2],vS_[i][1,0]]))
            
        #Vetor po_ em funcao de ph_:
        w_ = self.vw_[0]
        for i in range(1,self.dof):
        	w_ = Matrix([w_,self.vw_[i]])

        v_ = self.vv_[0]
        for i in range(1,self.dof):
        	v_ = Matrix([v_,self.vv_[i]])
            
        _Po_ = Matrix([w_,v_])

        #Indices dos elementos de _Po_ nao nulos      
        nonzerolist = []
        for i in range(6*self.dof):
            if _Po_[i] != 0:
                nonzerolist.append(i)
        
        #Indices dos elementos de _p_ nao nulos   
        nonzerolist2 = range(self.dof)
        for i in nonzerolist:
            nonzerolist2.append(i + self.dof)
                
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
        s = 0
        for i in range(self.dof):
            s += ( 
            	symbols('m_' + str(i+1) + Id)*( pv_[3*i].diff(t)**2 + pv_[3*i+1].diff(t)**2 + pv_[3*i+2].diff(t)**2 ) + 
            	symbols('Jx'+str(i+1) + Id)*pw_[3*i+0].diff(t)**2 + 
            	symbols('Jy'+str(i+1) + Id)*pw_[3*i+1].diff(t)**2 + 
            	symbols('Jz'+str(i+1) + Id)*pw_[3*i+2].diff(t)**2 + 
            	2*pw_[3*i+0].diff(t)*(symbols('Jz'+str(i+1) + Id)-symbols('Jy'+str(i+1) + Id))*pw_[3*i+2]*pw_[3*i+1] + 
            	2*pw_[3*i+1].diff(t)*(symbols('Jx'+str(i+1) + Id)-symbols('Jz'+str(i+1) + Id))*pw_[3*i+0]*pw_[3*i+2] + 
            	2*pw_[3*i+2].diff(t)*(symbols('Jy'+str(i+1) + Id)-symbols('Jx'+str(i+1) + Id))*pw_[3*i+1]*pw_[3*i+0] )/2
            
        #energia potencial
        ep = 0
        for i in range(self.dof):
            ep += symbols('m_' + str(i+1) + Id)*symbols('g')*self.qo_[3*i+2]
    
        #Matrizes da dinamica
        self.M_ = (Matrix([s]).jacobian(self.p_.diff(t))).jacobian(self.p_.diff(t))
        self.v_ = ( (Matrix([s]).jacobian(self.p_.diff(t))).T - (self.M_)*self.p_.diff(t) ).jacobian(self.p_)*self.p_ / 2
        g_ = Matrix([self.ph_,pv_]).jacobian(p_).T * Matrix([ep]).jacobian(self.q_).T
        self.g_ = g_[nonzerolist2, :]
        
        #Matrizes da dinamica hashtag
        self.Mh_ = simplify(self.C_.T * self.M_ * self.C_)
        self.vh_ = simplify(self.C_.T * (self.v_ + self.M_ * self.C_.diff(t) * self.ph_ ) )
        self.gh_ = simplify(self.C_.T * self.g_)
        
        #Balanceamento estatico
        if self.pflag == False:
            lgx = []
            for i in range(self.dof):
                lgx.append(symbols('lg_'+str(i+1) + Id))
            self.StaticBal = solve(self.gh_, lgx)
            
            self.gh_sb_ = simplify(self.gh_.subs(self.StaticBal))
            self.Mh_sb_ = simplify(self.Mh_.subs(self.StaticBal))
            self.vh_sb_ = simplify(self.vh_.subs(self.StaticBal))
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
       

RR = Serial("RR", '', Matrix([['x','x'],['y','y']]).T)
RR.description()