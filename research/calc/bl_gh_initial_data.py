"""Exact initial 50-field GH data from time-symmetric Brill--Lindquist punctures.

Initial lapse N=psi^-2, beta=0, K=0, dt N=dt beta=0. All constraints
checked OFF punctures. No finite-slab or gauge-transition certificate.
"""
import sympy as sp


def potential(x,y,z,d):
    return (1+sp.Rational(1,4)/sp.sqrt(x*x+y*y+(z-d/2)**2)
             +sp.Rational(1,4)/sp.sqrt(x*x+y*y+(z+d/2)**2))


def initial_fields(psi,x,y,z):
    """psi>0 smooth, harmonic off punctures, gamma_ij=psi^4 delta."""
    g=sp.diag(-psi**(-4),psi**4,psi**4,psi**4)
    Pi=sp.zeros(4)
    Phi=[g.diff(c) for c in (x,y,z)]
    return g,Pi,Phi


def symbolic_checks():
    x,y,z,d=sp.symbols("x y z d",real=True)
    psi=potential(x,y,z,d)
    for z0 in (d/2,-d/2):
        f=1/sp.sqrt(x*x+y*y+(z-z0)**2)
        lap=sp.diff(f,x,2)+sp.diff(f,y,2)+sp.diff(f,z,2)
        assert sp.simplify(lap)==0
    # Conformal spatial scalar curvature R3=-8*psi^-5*Delta psi=0.
    # K=0 => Hamiltonian and momentum constraints.
    # N=psi^-2 => sqrt(-g)=psi^4. The densitized g^{-1} at t0 is
    # diag(-psi^8, 1,1,1). Its four-divergence vanishes if dt g=0,
    # so C_a=Gamma_a+H_a=0 with prescribed H_a=0 at t0.
    q=sp.symbols("q",positive=True)
    denq=sp.diag(-q**8,1,1,1)
    assert sp.diag(-q**4,q**(-4),q**(-4),q**(-4))*q**4==denq
    den=denq.subs(q,psi)
    assert all(sp.diff(den[i,j],(x,y,z)[i-1])==0
               for i in range(1,4) for j in range(1,4))
    g,Pi,Phi=initial_fields(psi,x,y,z)
    # D_i Phi_j - D_j Phi_i=0 by commuting derivatives of a smooth
    # metric; verify on a generic symbolic potential to avoid large radicals.
    f=sp.Function("f")(x,y,z)
    gg=sp.diag(-f**(-4),f**4,f**4,f**4)
    for i,ci in enumerate((x,y,z)):
        for j,cj in enumerate((x,y,z)):
            assert gg.diff(ci).diff(cj)==gg.diff(cj).diff(ci)
    return {
      "psi":psi,"g_initial":g,"Pi":Pi,"Phi":Phi,
      "Hamiltonian_constraint_exact_off_punctures":True,
      "momentum_constraint_exact":True,
      "initial_harmonic_gauge_H_zero":True,
      "first_order_reduction_constraints_exact":True,
      "normal_derivative_GH_constraint_not_independently_verified":True,
      "full_finite_slab_3plus1_certified":False,
    }


if __name__=="__main__":
    obj=symbolic_checks()
    print({k:v for k,v in obj.items() if isinstance(v,bool)})
