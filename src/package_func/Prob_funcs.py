import numpy as np
import package_func.my_functions as mf

#Constants
R_earth=6371 #km

#angles
theta_12=33.45 #deg
theta_23=42.1
theta_13=8.62

#make this angles into rad for numpy
theta_12_r=mf.deg_to_rad(theta_12)
theta_23_r=mf.deg_to_rad(theta_23)
theta_13_r=mf.deg_to_rad(theta_13)

delta_cp=90

d_cp=mf.deg_to_rad(delta_cp) #deg -> this can change as much as we want

s_12=np.sin(theta_12_r)
s_23=np.sin(theta_23_r)
s_13=np.sin(theta_13_r)

c_12=np.cos(theta_12_r)
c_23=np.cos(theta_23_r)
c_13=np.cos(theta_13_r)


def PMNS_matrix():
    m=np.array([
        [c_12*c_13,s_12*c_13,s_13*np.exp(-d_cp*1j)],
        [-s_12*c_23-c_12*s_23*s_13*np.exp(d_cp*1j),c_12*c_23-s_12*s_23*s_13*np.exp(d_cp*1j),s_23*c_13],
        [s_12*s_23-c_12*c_23*s_13*np.exp(d_cp*1j),-c_12*s_23-s_12*c_23*s_13*np.exp(d_cp*1j),c_23*c_13]])

    return m

def Long(cosZenith):
    return abs(-2*R_earth*cosZenith)

def Prob_calc_General(a,b,cz,En,type="normal",ordering="NO"):
    index_a,index_b=mf.flavor_to_index(a,b)

    if type!="normal":
        U=PMNS_matrix().conjugate()
    else:
        U=PMNS_matrix()

    re_sum=0
    im_sum=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                term=U[index_b,i]*U[index_a,i].conjugate()*U[index_b,j].conjugate()*U[index_a,j]
                re_sum+=term.real*(np.sin(1.27*(mf.D_mass(i+1,j+1,ordering))*Long(cz)/(2*En)))**2#
                im_sum+=term.imag*(np.sin(1.27*(mf.D_mass(i+1,j+1,ordering))*Long(cz)/(En)))#*1.27*
                del term
    
    return mf.delta_Kro(a,b)-4*re_sum-2*im_sum

def Prob_a_to_b_Gen_MSW(a,b,cz,En,type="normal",ordering="NO"): ## DOES NOT WORK FOR En<10**0 for e to e 
    index_a,index_b=mf.flavor_to_index(a,b)

    if type!="normal":
        U=PMNS_matrix().conjugate()
    else:
        U=PMNS_matrix()

    re_sum=0
    im_sum=0 
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                if i+1==1 and j+1==3:
                    mass_term=mf.MSW_Dmass(i+1,j+1,cz,En,ordering)
                else:
                    mass_term=mf.D_mass(i+1,j+1,ordering)
                term=U[index_b,i]*U[index_a,i].conjugate()*U[index_b,j].conjugate()*U[index_a,j]
                re_sum+=term.real*(np.sin(1.27*(mass_term)*Long(cz)/(2*En)))**2#
                im_sum+=term.imag*(np.sin(1.27*(mass_term)*Long(cz)/(En)))#
                del term
    returner=mf.delta_Kro(a,b)-4*re_sum-2*im_sum
    return returner

def Prob_a_to_b_alt(a,b,cz,En,type="normal"):
    index_a,index_b=mf.flavor_to_index(a,b)
    
    fst_rhs=0
    if type!="normal":
        U=PMNS_matrix().conjugate()
    else:
        U=PMNS_matrix()

    for i in [0,1,2]:
        fst_rhs+=(mf.norm(U[index_b,i])**2)* mf.norm(U[index_a,i])**2

    scd_rhs=0
    for i in [0,1,2]:
        for j in [0,1,2]:
            if i>j:
                scd_rhs+=U[index_b,i]*(U[index_a,i].conjugate())*(U[index_b,j].conjugate())*U[index_a,j]*np.cos((mf.D_mass(i+1,j+1))*Long(cz)/(2*En))

    return fst_rhs+(2*scd_rhs.real)