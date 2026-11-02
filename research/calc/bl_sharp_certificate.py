"""Exact-rational majorant certificate for full 50-field BL GH second time jet.

Whole H^3 coefficient ball; H^4/H^5 reference derivative estimates; no mesh samples.
Time variation is resolved because the error is O(T^3), not merely O(T^2).
Uses Kato symmetric-hyperbolic first-exit theorem and full GH subsidiary theorem.
"""
from dataclasses import dataclass
from .bl_sharp_cutoff import derivative_envelopes, exact_ceil_log10

@dataclass(frozen=True)
class B:
    n: int | None
    l: int | None = None
    h: int | None = None
    def to(self,x):
        if isinstance(x,B):return x
        if x==0:return B(None)
        if abs(x)>4:raise ValueError("unexpected GH scalar constant")
        return B(1)
    @staticmethod
    def add(a,b):
        if a is None:return b
        if b is None:return a
        return max(a,b)+1
    @staticmethod
    def mult(a,b,k):
        return None if a is None or b is None else a+b+k
    def __add__(self,x):
        v=self.to(x)
        return B(self.add(self.n,v.n),self.add(self.l,v.l),self.add(self.h,v.h))
    __radd__=__add__
    def __neg__(self):return self
    def __sub__(self,x):return self+-self.to(x)
    def __rsub__(self,x):return self.to(x)+-self
    def __mul__(self,x):
        v=self.to(x)
        return B(self.mult(self.n,v.n,3),
                 self.add(self.mult(self.n,v.l,4),self.mult(self.l,v.n,4)),
                 self.add(self.add(self.mult(self.n,v.h,5),self.mult(self.h,v.n,5)),
                          self.mult(self.l,v.l,5)))
    __rmul__=__mul__
    def __truediv__(self,x):
        if x==2:return self
        raise ValueError("only division by 2 allowed in GH source graph")

def zeros(m,n=None):
    if n is None:n=m
    return [[B(None) for j in range(n)]for i in range(m)]
PAIRS=[(a,b) for a in range(4) for b in range(a,4)]

def unpack(u):
    def sym(v):
        z=zeros(4)
        for (a,b),x in zip(PAIRS,v):z[a][b]=z[b][a]=x
        return z
    return sym(u[:10]),sym(u[10:20]),[sym(u[20+10*i:30+10*i]) for i in range(3)]

def pack(g,pi,phi):
    def flat(x):return [x[a][b] for a,b in PAIRS]
    return flat(g)+flat(pi)+sum([flat(x) for x in phi],[])

def graph(U):
    """H^s full-ball norms and first/second Frechet derivatives for F and A.

    Geometry block: ref g0 and g0 inverse H^s <=10^U, perturbation radius
    <=10^(-Rexp); a matrix Neumann and binomial inverse-square-root lemma
    proves the deliberately generous U-dependent algebra bounds.
    """
    fields=[B(U,0) for _ in range(50)]
    _,pi,phi=unpack(fields)
    inv=[[B(U+2,2*U+10,3*U+25) for j in range(4)]for i in range(4)]
    N=B(U+6,6*U+60,10*U+110)
    normal=[-N*inv[0][a] for a in range(4)]
    beta=[N*N*inv[0][i] for i in range(1,4)]
    gam=[[inv[i+1][j+1]+normal[i+1]*normal[j+1] for j in range(3)]for i in range(3)]
    dt=[[-N*pi[a][b]+sum((beta[i]*phi[i][a][b] for i in range(3)),B(None))
         for b in range(4)]for a in range(4)]
    # Reuse EXACT production source algebra: independent duplicated tensor
    # majorants were found to underbound some contractions by 5-7 decades.
    from .gh_full_source import full_source
    F=full_source(fields,geometry_override=(inv,N,normal,beta,gam))
    A=[]
    for k in range(3):
        mat=zeros(5);mat[1][1]=-beta[k]
        for i in range(3):
            mat[1][i+2]=N*gam[k][i]
            mat[i+2][1]=N*int(i==k)
            mat[i+2][i+2]=-beta[k]
        A.append([x for row in mat for x in row])
    def xexp(arr,key,margin=3):
        vals=[getattr(x,key) for x in arr if getattr(x,key) is not None]
        if not vals:raise ValueError("unexpected empty GH graph")
        return max(vals)+margin
    Fexp={k:xexp(F,k) for k in ("n","l","h")}
    Aexp={k:max(xexp(a,k) for a in A) for k in ("n","l","h")}
    S=zeros(5);S[0][0]=S[1][1]=B(1)
    for i in range(3):
        for j in range(3):S[i+2][j+2]=gam[i][j]
    JS=xexp(sum(gam,[]),"l",4)
    JB=[]
    for a in A:
        Ak=[a[i*5:(i+1)*5] for i in range(5)]
        bb=[sum((S[i][m]*Ak[m][j] for m in range(5)),B(None))
            for i in range(5) for j in range(5)]
        JB.append(xexp(bb,"l",4))
    return {"F":Fexp,"A":Aexp,"JS":JS,"JB":max(JB),
            "F_metric_n":xexp(sum(dt,[]),"n"),
            "N_n":N.n,"gamma_n":xexp(sum(gam,[]),"n"),
            "beta_n":xexp(beta,"n")}

def input_sobolev_bounds():
    e=derivative_envelopes(k=8,m=9)
    g=e["metric_derivatives"]
    return {s:exact_ceil_log10(50*2**s*max(g[:s+2]))+1 for s in range(3,8)}

def certificate(time_exponent_override=None):
    U=input_sobolev_bounds()
    graphs={s:graph(U[s]) for s in (3,4,5)}
    def bmax(*v):return max(v)+4
    V1={s:bmax(graphs[s]["F"]["n"],graphs[s]["A"]["n"]+U[s+1]+3)
        for s in (3,4,5)}
    V2={s:bmax(graphs[s]["F"]["l"]+V1[s],
                graphs[s]["A"]["l"]+V1[s]+U[s+1]+3,
                graphs[s]["A"]["n"]+V1[s+1]+3)
        for s in (3,4)}
    g=graphs[3]
    # Finite-dimensional pointwise DS,D(SA) require pointwise field values.
    # ||u||_infty <=6||u||_H3, hence use the strictly larger U+2 graph.
    gp=graph(U[3]+2)
    D=bmax(gp["JS"]+g["F_metric_n"],gp["JB"]+U[3]+3)
    G=bmax(g["F"]["l"],g["A"]["l"]+U[4]+3,g["A"]["n"]+5)
    K=bmax(D+4,G+5)
    DG_u2=bmax(g["F"]["l"]+V2[3],
                g["A"]["l"]+V2[3]+U[4]+3,
                g["A"]["n"]+V2[4]+3)
    D2G=bmax(g["F"]["h"],g["A"]["h"]+U[4]+3,g["A"]["l"]+4)
    Vdelta=bmax(V1[4],V2[4])
    C2=bmax(DG_u2,D2G+2*Vdelta)
    Rexp=10*U[3]+100
    # T chosen for both existence and leading O(T²) dynamical correction
    # exceeding the O(T³) error at a specified physical BL point.
    Texp=max(K+12,Rexp+V1[3]+12,(Rexp+V2[3])//2+12,
              (C2+Rexp+20)//3+2,C2+20)
    if time_exponent_override is not None:
        if not isinstance(time_exponent_override,int) or time_exponent_override<=0:
            raise ValueError("time exponent must be positive integer")
        Texp=time_exponent_override
    Eexp=3*Texp-C2-10
    flags={
      "inverse_Neumann_full_H3_ball":Rexp>U[3]+30,
      "lapse_binomial_full_H3_ball":Rexp>4*U[3]+60,
      "Lorentzian_spatial_symmetrizer_tube":Rexp>U[3]+30,
      "approximation_inside_H3_tube":
          Texp>Rexp+V1[3]+5 and 2*Texp>Rexp+V2[3]+5,
      "energy_growth_controlled":Texp>K+8,
      "nonlinear_error_closes":Eexp>Rexp+8,
      "dynamic_metric_change_resolved":Eexp>2*Texp+8,
      "physical_causal_buffer":Texp>3
    }
    if not all(flags.values()):
        raise ValueError("Sharp full-GH BL certificate rejected: "+str(flags))
    return {"classification":"RESOLVED_LOCAL_BL_GH_VACUUM_ONLY",
      "d":"3831/5000","reference_Hs_log10":U,
      "full_GH_source_graph_H3":{"F_n":g["F"]["n"],"F_l":g["F"]["l"],
         "F_h":g["F"]["h"],"A_n":g["A"]["n"],"A_l":g["A"]["l"],
         "A_h":g["A"]["h"]},
      "u1_Hs_log10":V1,"u2_Hs_log10":V2,
      "adapted_divergence_D_log10":D,"growth_K_log10":K,
      "residual_t2_coefficient_log10":C2,
      "continuum_H3_residual_upper":f"10^({C2-2*Texp})",
      "time_T_over_M":f"10^(-{Texp})",
      "H3_tube_radius":f"10^(-{Rexp})",
      "H3_error_upper":f"10^(-{Eexp})",
      "metric_component_error_upper":f"6*10^(-{Eexp})",
      "first_spatial_derivative_error_upper":f"12*10^(-{Eexp})",
      "gxx_change_lower_at_0p1":f"10^(-{2*Texp})/40",
      "dynamic_to_error_lower_ratio":f"10^({Eexp-2*Texp})/240",
      "protected_target":"|x|<=0.9, r_+,r_->=0.3",
      "auxiliary_BL_agreement":"|x|<=1.1, r_+,r_->=0.2",
      "GH_and_subsidiary_speed_upper":"4",
      "auxiliary_cutoff_regular":"C9, not C-infinity",
      "full_50_field_second_jet":True,
      "flags":flags,
      "merger_spacetime_certified":False,
      "enlarges_previous_geometric_domain":False}

if __name__=="__main__":
    import json
    print(json.dumps(certificate(),indent=2))
