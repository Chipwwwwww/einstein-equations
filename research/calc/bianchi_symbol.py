"""Principal symbols of two Einstein-Bianchi reductions in Gauss gauge (notes/19 §2).

(F) Frame form (Friedrich-Rendall 2000, §4.4; unknowns e_a^mu, Gamma_b^i_j, E, B): frame and connection
    obey ODEs along u, so the principal part is only the (E, H) block
        d_t E = curl H + ...,  d_t H = -curl E + ...,
    with curl taken in the orthonormal frame. Symbol: P(xi) = [[0, C_xi], [-C_xi, 0]] on STF + STF.
(C) Coordinate form (unknowns gamma_ij, K_ij, E_ij, H_ij; Christoffels computed from gamma, as in
    bl_bianchi.py): the lower-order terms of curl contain d gamma, so the symbol has an extra block
    H-row <- gamma-column (and E-row <- gamma-column when H != 0).

The script computes: C_xi is antisymmetric (so (F) is symmetric hyperbolic with the identity
symmetriser); the speeds of (F) are 0, +-1/2, +-1 (Friedrich-Rendall's characteristic polynomial
xi_0^K (xi_0^2 + |xi|^2/4)^L (g(xi, xi))^N); and whether the symbol of (C) is diagonalisable.
"""

import itertools

import numpy as np

EPS = np.zeros((3, 3, 3))
for (i, j, k), s in {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.items():
    EPS[i, j, k] = s

SYM = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]


def to_mat(v):
    A = np.zeros((3, 3))
    for c, (i, j) in enumerate(SYM):
        A[i, j] = A[j, i] = v[c]
    return A


def to_vec(A):
    return np.array([A[i, j] for (i, j) in SYM])


def stf_basis():
    """Orthonormal basis (Frobenius) of symmetric trace-free 3x3 matrices."""
    B = []
    for (i, j) in SYM:
        A = np.zeros((3, 3))
        A[i, j] = A[j, i] = 1
        B.append(A - np.trace(A) / 3 * np.eye(3))
    M = np.array([A.ravel() for A in B]).T
    U, s, _ = np.linalg.svd(M, full_matrices=False)
    return [U[:, c].reshape(3, 3) for c in range(5)]


def curl_symbol(xi, A):
    """sym(eps_a^{cd} xi_c A_db) (flat orthonormal frame)."""
    T = np.einsum("acd,c,db->ab", EPS, xi, A)
    return (T + T.T) / 2


def frame_symbol(xi):
    Bs = stf_basis()
    C = np.array([[np.sum(Bs[r] * curl_symbol(xi, Bs[c])) for c in range(5)] for r in range(5)])
    P = np.block([[np.zeros((5, 5)), C], [-C, np.zeros((5, 5))]])
    return C, P


def christoffel_symbol(xi, dg):
    """Gamma^m_kl from d_k gamma_ab = xi_k dg_ab (flat background gamma = delta)."""
    G = np.zeros((3, 3, 3))
    for m, k, l in itertools.product(range(3), repeat=3):
        G[m, k, l] = 0.5 * (xi[k] * dg[m, l] + xi[l] * dg[m, k] - xi[m] * dg[k, l])
    return G


def curl_gamma_term(xi, dg, A):
    """Principal (d gamma) part of curl A: -sym(eps_i^{kl} Gamma^m_kj A_lm)."""
    G = christoffel_symbol(xi, dg)
    T = -np.einsum("ikl,mkj,lm->ij", EPS, G, A)
    return (T + T.T) / 2


def coordinate_symbol(xi, E, H):
    """Symbol of the coordinate form on (gamma, K, E, H), 6 components each, background gamma = delta."""
    n = 6
    P = np.zeros((4 * n, 4 * n))
    for c in range(n):
        e = np.zeros(n)
        e[c] = 1
        A = to_mat(e)
        # columns: gamma (0), K (1), E (2), H (3); rows likewise
        P[2 * n:3 * n, 3 * n + c] = to_vec(curl_symbol(xi, A))           # E <- curl H
        P[3 * n:4 * n, 2 * n + c] = -to_vec(curl_symbol(xi, A))          # H <- -curl E
        P[2 * n:3 * n, 0 * n + c] = to_vec(curl_gamma_term(xi, A, H))     # E <- d gamma (via curl H)
        P[3 * n:4 * n, 0 * n + c] = -to_vec(curl_gamma_term(xi, A, E))    # H <- d gamma (via curl E)
    return P


def jordan_test(P, tol=1e-9):
    """(algebraic multiplicity of eigenvalue 0, geometric multiplicity, rank P, rank P^2)."""
    ev = np.linalg.eigvals(P)
    alg = int(np.sum(np.abs(ev) < 1e-6))
    r1 = np.linalg.matrix_rank(P, tol)
    r2 = np.linalg.matrix_rank(P @ P, tol)
    return alg, P.shape[0] - r1, r1, r2


def main():
    rng = np.random.default_rng(0)
    xi = rng.normal(size=3)
    xi /= np.linalg.norm(xi)
    C, P = frame_symbol(xi)
    out = {"C_antisymmetric": float(np.abs(C + C.T).max()),
           "frame_speeds": sorted(np.round(np.linalg.eigvals(P).real, 10).tolist())}
    E = to_mat(rng.normal(size=6))
    E -= np.trace(E) / 3 * np.eye(3)
    out["coordinate_H0"] = jordan_test(coordinate_symbol(xi, E, np.zeros((3, 3))))
    out["coordinate_flat"] = jordan_test(coordinate_symbol(xi, np.zeros((3, 3)), np.zeros((3, 3))))
    return out


if __name__ == "__main__":
    print(main())
