"""Exact symbolic 50-field GH PRINCIPAL symbol; no full Einstein residual.

Lindblom et al. (2006) gamma1=-1, gamma2=0. The 50 fields are ten
symmetric spacetime metric components, each with (psi,Pi,Phi_x,Phi_y,Phi_z).
"""
import sympy as sp


def blocks(N,beta,spatial_inverse):
    gamma=sp.Matrix(spatial_inverse)
    if gamma.shape!=(3,3) or gamma!=gamma.T:
        raise ValueError("symmetric spatial inverse required")
    if len(beta)!=3 or N<=0:
        raise ValueError("positive lapse and three shifts required")
    S5=sp.diag(1,1,1,1,1)
    for i in range(3):
        for j in range(3):
            S5[i+2,j+2]=gamma[i,j]
    A5=[]
    for k in range(3):
        A=sp.zeros(5)
        A[1,1]=-beta[k]
        for j in range(3):
            A[1,j+2]=N*gamma[k,j]
            A[j+2,1]=N*(1 if j==k else 0)
            A[j+2,j+2]=-beta[k]
        A5.append(A)
    I10=sp.eye(10)
    return sp.kronecker_product(I10,S5),[
        sp.kronecker_product(I10,A) for A in A5]


def schwarzschild_ks_axis(M=sp.Integer(1),r=sp.Integer(2)):
    if M<=0 or r<=0:raise ValueError("M,r positive")
    h=2*M/r
    gamma=sp.diag(1/(1+h),1,1)
    N=1/sp.sqrt(1+h)
    beta=[h/(1+h),0,0]
    return N,beta,gamma


def exact_axis_report(M=sp.Integer(1),r=sp.Integer(2)):
    N,beta,ginv=schwarzschild_ks_axis(M,r)
    S,As=blocks(N,beta,ginv)
    if any(sp.simplify(S*A-A.T*S)!=sp.zeros(50) for A in As):
        raise ValueError("GH 50-field principal symmetrizer failed")
    if ginv.det()<=0 or N<=0:
        raise ValueError("hyperbolicity failed")
    eigmin=min(sp.Rational(1),*[ginv[i,i] for i in range(3)])
    assert eigmin>0
    return {
        "field_count":50,"symmetric_spacetime_components":10,
        "field_per_pair":5,"lapse":N,"shift_x":beta[0],
        "positive_symmetrizer_min_eigen":eigmin,
        "principal_symmetric_all_directions":True,
        "schwarzschild_ks_horizon_r_over_M":r/M,
        "full_einstein_source_and_constraints_verified":False
    }
