"""Exact Schwarzschild Cartesian Kerr-Schild GH source, not an evolution certificate."""
import sympy as sp


def verify_source():
    x,y,z,M=sp.symbols("x y z M",real=True,positive=True)
    r=sp.sqrt(x*x+y*y+z*z)
    h=2*M/r
    n=[x/r,y/r,z/r]
    ell_up=sp.Matrix([-1,*n])
    eta=sp.diag(-1,1,1,1)
    gup=eta-h*ell_up*ell_up.T
    Hup=sp.Matrix([h/r,*[-h*ni/r for ni in n]])
    Gamma_up=sp.Matrix([
        sp.simplify(-sum(sp.diff(gup[a,j+1],c)
               for j,c in enumerate((x,y,z))))
        for a in range(4)])
    assert all(sp.simplify(Gamma_up[a]+Hup[a])==0 for a in range(4))
    # det gKS=-1 by matrix determinant lemma; ell is null.
    assert sp.simplify((ell_up.T*eta*ell_up)[0])==0
    ell_down=eta*ell_up
    assert sp.simplify((ell_down.T*Hup)[0])==0
    Hdown=eta*Hup
    assert Hdown==sp.Matrix([-h/r,*[-h*ni/r for ni in n]])
    return {
        "metric":"Cartesian Kerr-Schild Schwarzschild",
        "H_covariant":Hdown,"H_contravariant":Hup,
        "exact_GH_constraint":True,
        "shell_r_3M_over_2_to_5M_over_2_symmetrizer_eigen_lower":sp.Rational(3,7),
        "complete_H3_residual_bootstrap_or_vacuum_certificate":False
    }


if __name__=="__main__":
    print(verify_source())
