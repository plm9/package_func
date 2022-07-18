import numpy as np
from sympy import *

theta_12= 33.45  #degrees
theta_23= 42.1   #degrees
theta_13= 8.62   #degrees

delta_cp= 270    #degrees

Dm_21   = 7.42*10**(-5) #eV^2

#Normal ordering
Dm_31_no= 2.51*10**(-3) #eV^2
Dm_32_no= 2.51*10**(-3) #eV^2

#Inverse ordering
Dm_31_io= -2.49*10**(-3) #eV^2
Dm_32_io= -2.49*10**(-3) #eV^2

#Defining symbols
s_12=Symbol("s_12",real=True)
s_23=Symbol("s_23",real=True)
s_13=Symbol("s_13",real=True)

c_12=Symbol("c_12",real=True)
c_23=Symbol("c_23",real=True)
c_13=Symbol("c_13",real=True)

d_cp=Symbol("d_cp",real=True)

D_m_21=Symbol("D_m_12",real=True)
D_m_31=Symbol("D_m_13",real=True)
D_m_32=Symbol("D_m_23",real=True)

L=Symbol("L",real=True,positive=True)
E_nu=Symbol("E_nu",real=True)

#Constants
G_F=1.1163787*10**(-5)#GeV^-2

""" def Y_e(L):
    if L<15:
        return 
    else:
        return """

rho=5.51  #g/cm^3
#N_e=Y_e(L)*rho

def PMNS_param_matrix():
    m=Matrix([
        [c_12*c_13,s_12*c_13,s_13*exp(-d_cp*1j)],
        [-s_12*c_23-c_12*s_23*s_13*exp(d_cp*1j),c_12*c_23-s_12*s_23*s_13*exp(d_cp*1j),s_23*c_13],
        [s_12*s_23-c_12*c_23*s_13*exp(d_cp*1j),-c_12*s_23-s_12*c_23*s_13*exp(d_cp*1j),c_23*c_13]])

    return m

def deg_to_rad(theta):
    return (theta/180)*np.pi

def rad_to_deg(theta):
    return (theta/np.pi)*180

def flavor_to_index(a,b):
    if a == "e":
        index_a=0
    elif a == "mu":
        index_a=1
    elif a == "tau":
        index_a=2
    else:
        raise ValueError("Non existing neutrino flavor for initial state.")
    
    if b == "e":
        index_b=0
    elif b == "mu":
        index_b=1
    elif b == "tau":
        index_b=2
    else:
        raise ValueError("Non existing neutrino flavor for final state.")
    
    return index_a,index_b


def D_mass(i,j,ordering="NO"): #for ordering: NO is normal ordering and IO inverse ordering
    if i==j:
        return 0
    elif i==1:
        if j==2:
            return -Dm_21
        elif j==3:
            if ordering=="NO":
                return -Dm_31_no
            else:
                return -Dm_31_io
    elif i==2 :
        if j==1:
            return Dm_21
        elif j==3:
            if ordering=="NO":
                return -Dm_32_no
            else:
                return -Dm_32_io
    elif i==3:
        if j==1:
            if ordering=="NO":
                return Dm_31_no
            else:
                return Dm_31_io
        elif j==2:
            if ordering=="NO":
                return Dm_32_no
            else:
                return Dm_32_io

def D_mass_param(i,j):
    if i==j:
        return 0
    elif i==1:
        if j==2:
            return -D_m_21
        elif j==3:
            return -D_m_31
    elif i==2 :
        if j==1:
            return D_m_21
        elif j==3:
            return -D_m_32
    elif i==3:
        if j==1:
            return D_m_31
        elif j==2:            
            return D_m_32

def norm(z):
    return sqrt(re(z)**2+im(z)**2)

def theta(i,j):
    if i==1 or j==1:
        if i==2 or j==2:
            return theta_12
        elif i==3 or j==3:
            return theta_13
    elif i==2 or j==2:
        if i==3 or j==3:
            return theta_23
    
def delta_Kro(a,b):
    if a==b:
        return 1
    else:
        return 0

def MSW_Dmass(j,i):
    A=0 #THIS IS 0 FOR VACUUM BUT NOT FOR MATTER. FOR MATTER A=2sqrt(2)G_F N_e E_\nu
    term=D_mass_param(j,i)*sqrt((cos(2*theta(i,j))-(A/D_mass_param(i,j)))**2+(sin(2*theta(i,j)))**2)
    return term 

def MSW_angle(i,j):
    A=0 #THIS IS 0 FOR VACUUM BUT NOT FOR MATTER. FOR MATTER A=2sqrt(2)G_F N_e E_\nu
    term=0.5*asin((sin(2*theta(i,j))**2)/((cos(2*theta(i,j))-(A/D_mass_param(i,j)))**2+sin(2*theta(i,j))**2))
    return term


def theta_MSW(i,j):
    if i==1 or j==1:
        if i==2 or j==2:
            return theta_12
        elif i==3 or j==3:
            return MSW_angle(1,3)
    elif i==2 or j==2:
        if i==3 or j==3:
            return theta_23