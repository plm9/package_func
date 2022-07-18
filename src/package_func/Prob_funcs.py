import numpy as np
import package_func.my_functions as mf

#Constants
R_earth=6371 #km

#angles
theta_12=33.45 #deg
theta_23=42.1
theta_13=8.62

d_cp=270 #deg -> this can change as much as we want

s_12=np.sin(theta_12)
s_23=np.sin(theta_23)
s_13=np.sin(theta_13)

c_12=np.cos(theta_12)
c_23=np.cos(theta_23)
c_13=np.cos(theta_13)


def PMNS_matrix():
    m=np.array([
        [c_12*c_13,s_12*c_13,s_13*np.exp(-d_cp*1j)],
        [-s_12*c_23-c_12*s_23*s_13*np.exp(d_cp*1j),c_12*c_23-s_12*s_23*s_13*np.exp(d_cp*1j),s_23*c_13],
        [s_12*s_23-c_12*c_23*s_13*np.exp(d_cp*1j),-c_12*s_23-s_12*c_23*s_13*np.exp(d_cp*1j),c_23*c_13]])

    return m

def Long(cosZenith):
    return -2*R_earth*cosZenith

def Prob_calc_General(a,b,cz,En):
    index_a,index_b=mf.flavor_to_index(a,b)

    U=PMNS_matrix()

    re_sum=0
    im_sum=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                term=U[index_b,i]*U[index_a,i].conjugate()*U[index_b,j].conjugate()*U[index_a,j]
                re_sum+=term.real*(np.sin((mf.D_mass(i+1,j+1))*Long(cz)/(4*En)))**2
                im_sum+=term.imag*(np.sin((mf.D_mass(i+1,j+1))*Long(cz)/(2*En)))
                del term
    
    return mf.delta_Kro(a,b)-4*re_sum-2*im_sum

def Prob_a_to_b_Gen_MSW(a,b,cz,En):
    index_a,index_b=mf.flavor_to_index(a,b)

    U=PMNS_matrix()

    re_sum=0
    im_sum=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                if i+1==1 and j+1==3:
                    mass_term=mf.MSW_Dmass(i+1,j+1)
                else:
                    mass_term=mf.D_mass(i+1,j+1)
                term=U[index_b,i]*U[index_a,i].conjugate()*U[index_b,j].conjugate()*U[index_a,j]
                re_sum+=term.real*(np.sin((mass_term)*Long(cz)/(4*En)))**2
                im_sum+=term.imag*(np.sin((mass_term)*Long(cz)/(2*En)))
                del term
    returner=mf.delta_Kro(a,b)-4*re_sum-2*im_sum
    return returner