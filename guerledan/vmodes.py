import xarray as xr
import numpy as np

import scipy.sparse as sp
import scipy.sparse.linalg as linalgs
from scipy.sparse import bmat, eye
import scipy.linalg as linalg

sparse_format = "lil"

g = 9.81

def compute_vmodes(N2, z, nmodes=10, U=None, rho0=None, kdv=False, nlinear=False):
    """ Compute vertical normal modes.
    Solve for buoyancy modes.
    
    $
    \frac{d}{dz} \Big [ \rho_0 (c - U(z))^2  \frac{d\Phi}{dz} \Big ] + \rho_0 N^2 \Phi = 0, \\
    \Phi(z=-h) = 0, \\
    (c-U(z=0))^2 \frac{d\Phi}{dz}(z=0) - g \Phi(z=0) = 0
    $
    
    Grimshaw et al. 2002

    Returns
    -------
    vm_pos: xr.Dataset
        dataset of positive eigenvalues
    vm_neg: xr.Dataset
        dataset of negative eigenvalues
    """

    # assemble into a single dataset for easier preprocessing manipulations:
    ds = xr.merge([N2.rename("N2"), z.rename("z")])
    ds = ds.where(~np.isnan(ds.N2), drop=True)
    ds = ds.sortby("z", ascending=True)
    # cell are from bottom to top, sea floor is at zf[0], sea level at zf[-1]
    zf = ds["z"].values
    zc = (zf[:-1] + zf[1:]) * 0.5
    N2 = ds["N2"].values
    # N2 = N2*0 + 1e-4 # dev, test

    N = len(zf)
    dzc = np.diff(zc)
    dzf = np.diff(zf)

    # should eventually be passed as input parameters
    if U is None:
        U = np.zeros(zf.shape)
    if rho0 is None:
        rho0 = np.ones(zf.shape)
        
    phi, dphi, c = solve_vmodes_linear(N2, dzf, dzc, rho0, U, nmodes)
    
    vm = xr.Dataset(
        dict(
            phi=(("z", "mode"), phi),
            c=("mode", c),
        ),
        coords=dict(z=zf),
    )
    vm["dphidz"] = (("z", "mode"), dphi)
    vm = vm.assign_coords(
        U=(("z"), U), 
        rho0=(("z"), rho0),
        N2=(("z"), N2),
        H=np.abs(zf).max(),
        N_max=np.nanmax(np.sqrt(N2)),
    )

    if kdv:
        co = get_kdv_coeffs(vm)
        alpha, beta = co.alpha.values, co.beta.values
        vm = xr.merge([vm, co])
        
    if nlinear:
    
        Tn, dTn, Td, dTd = solve_vmodes_nlinear_corrections(N2, dzf, dzc, rho0, U, phi, dphi, c, alpha, beta)
        vm["Tn"] = (("z", "mode"), Tn)
        vm["dTndz"] = (("z", "mode"), dTn)
        vm["Td"] = (("z", "mode"), Td)
        vm["dTddz"] = (("z", "mode"), dTd)

        if kdv:
            co = get_kdv_coeffs(vm, nlinear=True)
            vm = xr.merge([vm, co])

    # drop bad modes
    #vm = vm.where((~np.isnan(vm.c)) & (np.abs(vm.c) <= 400), drop=True)

    # sort and split into positive/negative eigenvalues
    vm_p = (
        vm
        .where(vm.c > 0, drop=True)
        .sortby("c", ascending=False)
    )
    vm_p = vm_p.assign_coords(mode=np.arange(vm_p.mode.size))
    vm_n = (
        vm
        .where(vm.c < 0, drop=True)
        .sortby("c", ascending=True)
    )
    vm_n = vm_n.assign_coords(mode=np.arange(vm_n.mode.size))

    return vm_p, vm_n
    
    
def solve_vmodes_linear(N2, dzf, dzc, rho0, U, nmodes):

    N = len(dzf)+1

    # assemble quadratic eigenvalue problem
    A0 = sp.lil_matrix((N, N))
    A1 = sp.lil_matrix((N, N))
    A2 = sp.lil_matrix((N, N))
    # interior
    s = np.mean(dzf) ** 2 / np.max(N2)  # scaling coefficient
    for i in range(1, N - 1):
        # equation at zf i+1
        rho_up = (rho0[i] + rho0[i + 1]) * 0.5
        rho_down = (rho0[i] + rho0[i - 1]) * 0.5
        U_up = (U[i] + U[i + 1]) * 0.5
        U_down = (U[i] + U[i - 1]) * 0.5
        idz_up = 1 / dzf[i] / dzc[i - 1]
        idz_down = 1 / dzf[i - 1] / dzc[i - 1]
        #
        A2[i, i + 1] = rho_up * idz_up * s
        A2[i, i] = -(rho_up * idz_up + rho_down * idz_down) * s
        A2[i, i - 1] = rho_down * idz_down * s
        #
        A1[i, i + 1] = -2 * rho_up * U_up * idz_up * s
        A1[i, i] = 2 * (rho_up * U_up * idz_up + rho_down * U_down * idz_down) * s
        A1[i, i - 1] = -2 * rho_down * U_down * idz_down * s
        #
        A0[i, i + 1] = rho_up * U_up**2 * idz_up * s
        A0[i, i] = (
            -(rho_up * U_up**2 * idz_up + rho_down * U_down**2 * idz_down)
            + rho0[i] * N2[i]
        ) * s
        A0[i, i - 1] = rho_down * U_down**2 * idz_down * s
    # lower boundary condition is phi=0
    i = 0
    A0[i, i] = 1.0
    # upper boundary condition
    i = N - 1
    s = dzf[i - 1]
    A2[i, i] = 1 / dzf[i - 1] * s
    A2[i, i - 1] = -1 / dzf[i - 1] * s
    A1[i, i] = 2 * U[i] / dzf[i - 1] * s
    A1[i, i - 1] = -2 * U[i] / dzf[i - 1] * s
    A0[i, i] = (U[i] ** 2 / dzf[i - 1] - g) * s
    A0[i, i - 1] = -U[i] ** 2 / dzf[i - 1] * s

    # compute modes
    # ev,ef = la.eigs(L.tocsc(),nmodes+1,sigma=sigma)
    phi, c = polyeig((A0, A1, A2), nmodes, e_max=1000, e_split=True)
    c = np.real(c)
    phi = np.real(phi)

    # massage c, sort appropriately
    # ii = c.argsort()
    # c=c[ii]              # c value (c-U with a background flow)
    # phi=X[:,ii]         # phi

    # normalizes such that phi_max = 1
    for m in range(phi.shape[1]):
        imax = np.abs(phi[:,m]).argmax()
        phi[:, m] = phi[:, m] / phi[imax,m]

    # computes vertical derivative manually
    dphi = phi * 0
    for i in range(1, N - 1):
        # equation at zf i+1
        dphi[i,:] = (
            (phi[i + 1,:] - phi[i, :]) / dzf[i] + (phi[i,:] - phi[i - 1,:]) / dzf[i - 1]
        ) * 0.5
    dphi[-1, :] = g * phi[-1, :] / (c - U[-1]) ** 2
    dphi[0, :] = dphi[1, :]
        
    return phi, dphi, c


def solve_vmodes_nlinear_corrections(N2, dzf, dzc, rho0, U, phi, dphi, c, alpha, beta):

    N, nmodes = phi.shape

    # assemble L operator
    L0 = sp.lil_matrix((N, N))  # c**0 terms
    L1 = sp.lil_matrix((N, N))  # c**1 terms
    L2 = sp.lil_matrix((N, N))  # c**2 terms
    RHSn = np.zeros((N,))
    RHSd = np.zeros((N,))
    # interior
    for i in range(1, N - 1):
        # equation at zf i+1
        rho_up = (rho0[i] + rho0[i + 1]) * 0.5
        rho_down = (rho0[i] + rho0[i - 1]) * 0.5
        U_up = (U[i] + U[i + 1]) * 0.5
        U_down = (U[i] + U[i - 1]) * 0.5
        idz_up = 1 / dzf[i] / dzc[i - 1]
        idz_down = 1 / dzf[i - 1] / dzc[i - 1]
        #
        L2[i, i + 1] = rho_up * idz_up
        L2[i, i] = -(rho_up * idz_up + rho_down * idz_down)
        L2[i, i - 1] = rho_down * idz_down
        #
        L1[i, i + 1] = -2 * rho_up * U_up * idz_up
        L1[i, i] = 2 * (rho_up * U_up * idz_up + rho_down * U_down * idz_down)
        L1[i, i - 1] = -2 * rho_down * U_down * idz_down
        #
        L0[i, i + 1] = rho_up * U_up**2 * idz_up
        L0[i, i] = (
            -(rho_up * U_up**2 * idz_up + rho_down * U_down**2 * idz_down)
            + rho0[i] * N2[i]
        )
        L0[i, i - 1] = rho_down * U_down**2 * idz_down
    # lower boundary condition is phi=0
    i = 0
    L0[i, i] = 1.0
    # upper boundary condition is updated for each mode
    
    Tn = np.zeros_like(phi)
    Td = np.zeros_like(phi)
    
    # !! should loop around m here
    for m in range(nmodes):
    
        _alpha, _beta, _c = alpha[m], beta[m], c[m]
        # interior
        L = _c**2 * L2 + _c * L1 + L0
        RHSn[:] = 0.
        RHSd[:] = 0.
        for i in range(1, N - 1):
            rho_up = (rho0[i] + rho0[i + 1]) * 0.5
            rho_down = (rho0[i] + rho0[i - 1]) * 0.5
            cU_up = _c - (U[i] + U[i + 1]) * 0.5
            cU_down = _c - (U[i] + U[i - 1]) * 0.5
            dphi_up = (phi[i+1,m]-phi[i,m])/ dzf[i]
            dphi_down = (phi[i,m]-phi[i-1,m])/ dzf[i-1]
            RHSn[i] = (
                -_alpha * 
                    (
                        rho_up*cU_up*dphi_up 
                        - rho_down*cU_down*dphi_down
                    ) / dzc[i - 1]
                + 3/2 * (
                        rho_up*cU_up**2*dphi_up**2 
                        - rho_down*cU_down**2*dphi_down**2
                    ) / dzc[i - 1]
            )
            RHSd[i] = (
                -2*_beta * 
                    (
                        rho_up*cU_up*dphi_up
                        - rho_down*cU_down*dphi_down
                    ) / dzc[i - 1]
                - rho0[i] * (_c - U[i])**2 * phi[i, m]
            )
        # fix upper boundary operator and condition and solve
        i = N - 1
        cU = _c - U[i]
        L[i, i] = g - cU**2 / dzf[i - 1]
        L[i, i - 1] = cU**2 / dzf[i - 1]
        RHSn[i] = _alpha * cU * dphi[i, m] + 3/2 * cU**2 * dphi[i, m]**2
        RHSd[i] = 2* _beta * cU * dphi[i, m]

        # solve
        Tn[:,m] = linalgs.spsolve(L, RHSn)
        Td[:,m] = linalgs.spsolve(L, RHSd)
        # tested with non sparse arrays
        #_L = L.toarray()
        #Tn[:,m] = linalg.solve(_L, RHSn)
        #Td[:,m] = linalg.solve(_L, RHSd)
    
    # normalizes such that phi_max = 1
    for m in range(phi.shape[1]):
        imax = np.abs(phi[:,m]).argmax()
        Tn[:,m] = Tn[:,m] - Tn[imax,m] * phi[:,m]
        Td[:,m] = Td[:,m] - Td[imax,m] * phi[:,m]
    
    # computes vertical derivatives manually
    dTn = phi * 0
    dTd = phi * 0
    for i in range(1, N - 1):
        # equation at zf i+1
        dTn[i, :] = (
            (Tn[i + 1, :] - Tn[i, :]) / dzf[i] 
            + (Tn[i, :] - Tn[i - 1,:]) / dzf[i - 1]
        ) * 0.5
        dTd[i, :] = (
            (Td[i + 1, :] - Td[i, :]) / dzf[i] 
            + (Td[i, :] - Td[i - 1,:]) / dzf[i - 1]
        ) * 0.5
    dTn[-1, :] = dTn[-2, :]
    dTn[0, :] = dTn[1, :]
    dTd[0, :] = dTd[1, :]
    dTd[-1, :] = dTd[-2, :]
    
    return Tn, dTn, Td, dTd


def polyeig(A, nmodes, sigma=None, linearization=0, sparse=False, e_max=None, e_split=False):
    """
    Solve the polynomial eigenvalue problem:
        (A0 + e A1 +...+  e**p Ap)x=0
    We consider only quadratic form in fact: (K + C e + M e**2) x = 0
        X,e = polyeig((K,C),M)
        https://en.wikipedia.org/wiki/Quadratic_eigenvalue_problem

    Return the eigenvectors [x_i] and eigenvalues [e_i] that are solutions.
    """
    if len(A) <= 0:
        raise Exception("Provide at least one matrix")
    for Ai in A:
        if Ai.shape[0] != Ai.shape[1]:
            raise Exception("Matrices must be square")
        if Ai.shape != A[0].shape:
            raise Exception("All matrices must have the same shapes")

    n = A[0].shape[0]
    l = len(A) - 1
    assert l == 2, "Only quadratic eigenvalue problem implemented for now"

    # Assemble matrices for generalized problem
    if linearization == 0:
        # C = bmat([[-A[i] for i in range(l-1,-1,-1)], [eye(n*(l-1)), None] ])
        # D = bmat([[A[-1], None],[None, eye(n*(l-1))]]);
        C = bmat([[None, eye(n)], [-A[0], -A[1]]])
        D = bmat([[eye(n), None], [None, A[2]]])
        # D = D+eye(D.shape[0])*1e-6 # not working
    elif linearization == 1:
        C = bmat([[-A[0], None], [A[1], A[2]]])
        D = bmat([[None, eye(n)], [eye(n), None]])
    elif linearization == 2:
        C = bmat([[A[0], None], [A[1], A[0]]])
        D = bmat([[None, A[0]], [-A[2], None]])
    elif linearization == 3:
        C = bmat([[None, -A[0]], [A[2], None]])
        D = bmat([[A[2], A[1]], [None, A[2]]])
    # Solve generalized eigenvalue problem
    if sparse:
        # this is not working at the moment, most likely because M does not meet eigs requirements
        e, X = linalgs.eigs(
            C, M=D, k=2*nmodes, sigma=sigma, which="LM", maxiter=10000
        )
    else:
        # looses benefit of sparsity but sparse solver fails at the moment
        e, X = linalg.eig(C.toarray(), D.toarray())
        if not e_split:
            idx = np.argsort(np.abs(e))[::-1][:nmodes]
            X = X[:, idx]
            e = e[idx]
        else:
            # select nmodes largest positive eigenvalues and nmodes largest negative eigenvalues
            e = np.real(e) # e has to be real
            if e_max is not None:
                idx_p = np.where( (e>0) & (np.abs(e)<e_max) )[0]
                idx_n = np.where( (e<0) & (np.abs(e)<e_max) )[0]
            else:
                idx_p = np.where( (e>0) )[0]
                idx_n = np.where( (e<0) )[0]
            ep = e[idx_p]
            en = e[idx_n]
            idx_p1 = np.argsort(ep)[::-1][:nmodes]
            idx_n1 = np.argsort(en)[:nmodes]
            ep = e[idx_p][idx_p1]
            en = e[idx_n][idx_n1]
            e = np.hstack([ep, en])
            X = np.hstack(
                [X[:, idx_p][:,idx_p1], X[:, idx_n][:,idx_n1]]
            )
    if np.all(np.isreal(e)):
        e = np.real(e)
    X = X[:n, :]

    return X, e



def get_kdv_coeffs(vm, U=None, rho0=None, nlinear=False):
    """get KdV coefficients"""
    if U is None:
        if "U" in vm:
            U = vm["U"]
        else:
            U = 0.0
    if rho0 is None:
        if "rho0" in vm:
            rho0 = vm["rho0"]
        else:
            rho0 = 1.0
    phi, dphi, c = vm.phi, vm.dphidz, vm.c
    I = 2 * (rho0 * (c - U) * dphi**2).integrate("z")
    alpha = 3 * (rho0 * (c - U) ** 2 * dphi**3).integrate("z") / I
    beta = (rho0 * (c - U) ** 2 * phi**2).integrate("z") / I
    if not nlinear:
        co = xr.Dataset(dict(alpha=alpha, beta=beta))
        return co
    # nonlinear coefficients
    if "Td" in vm and "Tn" in vm:
        Td, dTd = vm.Td, vm.dTddz
        Tn, dTn = vm.Tn, vm.dTndz
        alpha1 = (
            rho0 * 
             (
                 3 * (c - U) ** 2 *  ( 3*dTn - 2*dphi**2 ) * dphi**2
                 - alpha**2 * dphi**2 
                 + alpha*(c-U) *( 5 * dphi**2 - 4* dTn ) * dphi
             )
        ).integrate("z") / I
        beta1 = (
            rho0 * 
              (
                  (c - U) ** 2 * phi * Td 
                  - beta**2 * dphi**2
                  + 2*beta*(c-U)*(phi**2 - dphi*dTd)
              )
        ).integrate("z") / I
        gamma1 = -(
            rho0 * (
                2*alpha*beta*dphi**2 - 2*alpha*(c-U)*phi**2
                +(c-U)**2*phi**2*dphi 
                - (c-U)**2 * (
                    2*Tn*phi+3*dTd*dphi**2
                )
                + 2*(c-U) * (
                    alpha*dTd+2*beta*dTn
                )*dphi
                -4*beta*(c-U)*dphi**3
            )
        ).integrate("z") / I
        gamma2 = (
            rho0 * (
                (c-U)*(2*beta*dphi**3 + 6*alpha*phi**2)
                - 3*alpha*beta*dphi**2 - 2*(c-U)**2 * (phi**2*dphi -3*Tn*phi)
                -6*alpha*(c-U)*dTd*dphi
                +3*(c-U)**2*dTd*dphi**2
            )
        ).integrate("z") / I
    co = xr.Dataset(dict(alpha1=alpha1, beta1=beta1, gamma1=gamma1, gamma2=gamma2))
        
    return co
