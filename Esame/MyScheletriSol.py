import math
from math import copysign
import numpy as np
import scipy.linalg as spl

from SolveTriang import LSolve
from SolveTriang import USolve

########


def sign(x):
    return copysign(1, x)


def Bisezione(fname, a, b, tolx):
    fa = fname(a)
    fb = fname(b)

    if sign(fa) * sign(fb) >= 0:  ########
        # Non è possibile applicare il metodo della bisezione
        return None, None, None

    it = 0
    vxk = []

    maxit = math.ceil(math.log((b - a) / tolx) / math.log(2)) - 1

    while abs(b - a) > tolx:  ########
        xk = a + (b - a) / 2  ########
        vxk.append(xk)
        it += 1
        fxk = fname(xk)
        if fxk == 0:
            return xk, it, vxk

        if sign(fa) * sign(fxk) > 0:
            a = xk  ########
            fa = fxk  ########
        elif sign(fb) * sign(fxk) > 0:
            b = xk  ########
            fb = fxk  ########

    return xk, it, vxk


def Falsi(fname, a, b, tolx, tolf, nmax):
    fa = fname(a)
    fb = fname(b)

    if sign(fa) * sign(fb) >= 0:  ########
        # Non è possibile applicare il metodo della falsa posizione
        return None, None, None

    it = 0
    vxk = []
    fxk = 10

    while it < nmax and abs(b - a) > tolx and abs(fxk > tolf):  ########
        xk = a - fa * (b - a) / (fb - fa)  ########
        vxk.append(xk)
        it += 1
        fxk = fname(xk)

        if fxk == 0:
            return xk, it, vxk

        if sign(fa) * sign(fxk) > 0:
            a = xk  ########
            fa = fxk  ########
        elif sign(fb) * sign(fxk) > 0:
            b = xk  ########
            fb = fxk  ########

    return xk, it, vxk


def Corde(fname, m, x0, tolx, tolf, nmax):
    vxk = []
    fx0 = fname(x0)  ########
    d = fx0 / m  ########
    x1 = x0 - d  ########
    fx1 = fname(x1)
    vxk.append(x1)
    it = 1

    while it < nmax and abs(fx1 / m) >= tolx * abs(x1) and abs(fx1) >= tolf:  ########
        x0 = x1  ########
        fx0 = fname(x0)  ########
        d = fx0 / m  ########
        """
        x1= ascissa del punto di intersezione tra  la retta che passa per il punto
        (xi,f(xi)) e ha pendenza uguale a m  e l'asse x
        """
        x1 = x0 - d  ########
        fx1 = fname(x1)
        it += 1
        vxk.append(x1)

        if it == nmax:
            # Max Num di iterazioni
            return x1, it, vxk
    return x1, it, vxk


def Secanti(fname, xm1, x0, tolx, tolf, nmax):
    vxk = []
    fxm1 = fname(xm1)  ########
    fx0 = fname(x0)  ########
    d = fx0 * (x0 - xm1) / (fx0 - fxm1)  ########
    x1 = x0 - d  ########
    fx1 = fname(x1)
    vxk.append(x1)
    it = 1

    while it <= nmax and abs(d) >= tolx * abs(x1) and abs(fx1) >= tolf:
        xm1 = x0  ########
        x0 = x1  ########
        fxm1 = fname(xm1)  ########
        fx0 = fname(x0)  ########
        d = fx0 * (x0 - xm1) / (fx0 - fxm1)  ########
        x1 = x0 - d  ########

        fx1 = fname(x1)
        vxk.append(x1)
        it += 1

    return x1, it, vxk


def Newton(fname, fpname, x0, tolx, tolf, nmax):
    vxk = []
    fx0 = fname(x0)
    fpx0 = fpname(x0)
    if abs(fpx0) <= np.spacing(1):  ########
        # derivata prima nulla
        return None, None, None

    d = fx0 / fpx0  ########
    x1 = x0 - d  ########

    fx1 = fname(x1)
    vxk.append(x1)
    it = 1

    while it <= nmax and abs(fx1) >= tolf and abs(d) >= tolx * abs(x1):  ########
        x0 = x1  ########
        fx0 = fname(x0)  ########
        fpx0 = fpname(x0)
        if abs(fpx0) <= np.spacing(1):  ########
            # derivata prima nulla
            return None, None, None

        d = fx0 / fpx0  ########
        x1 = x0 - d  ########
        fx1 = fname(x1)
        it += 1

        vxk.append(x1)

    return x1, it, vxk


def NewtonMod(fname, fpname, m, x0, tolx, tolf, nmax):
    vxk = []
    fx0 = fname(x0)  ########
    fpx0 = fpname(x0)
    if abs(fpx0) <= np.spacing(1):  ########
        # derivata prima nulla in x0
        return None, None, None

    d = fx0 / fpx0  ########
    x1 = x0 - m * d  ########

    fx1 = fname(x1)  ########
    vxk.append(x1)
    it = 1

    while it <= nmax and abs(fx1) >= tolf and abs(d) >= tolx * abs(x1):  ########
        x0 = x1  ########
        fx0 = fname(x0)  ########
        fpx0 = fpname(x0)
        if abs(fpx0) <= np.spacing(1):  ########
            # derivata prima nulla in x0
            return None, None, None
        d = fx0 / fpx0  ########
        """
           x1= ascissa del punto di intersezione tra  la retta che passa per il punto
           (xi,f(xi)) ed è tangente alla funzione f(x) nel punto (xi.f(xi))  e l'asse x
        """
        x1 = x0 - m * d  ########
        fx1 = fname(x1)
        it = it + 1

        vxk.append(x1)

    return x1, it, vxk


def StimaOrdine(xk, it):
    k = it - 4
    ordine = np.log(abs(xk[k + 2] - xk[k + 3]) / abs(xk[k + 1] - xk[k + 2])) / np.log(
        abs(xk[k + 1] - xk[k + 2]) / abs(xk[k] - xk[k + 1])
    )
    return ordine


def NewtonSys(fun, jac, x0, tolx, tolf, nmax):
    matjac = jac(x0)
    if np.linalg.det(matjac) == 0:  ########
        return None, None, None

    s = -np.linalg.solve(matjac, fun(x0))  ########

    it = 1
    x1 = x0 + s  ########
    fx1 = fun(x1)
    vxm = [np.linalg.norm(s, 1) / np.linalg.norm(x1, 1)]

    while (
        it <= nmax
        and np.linalg.norm(fx1, 1) >= tolf
        and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1)
    ):  ########
        x0 = x1  ########
        it += 1

        matjac = jac(x0)
        if np.linalg.det(matjac) == 0:  ########
            return None, None, None

        s = -np.linalg.solve(matjac, fun(x0))  ########

        x1 = x0 + s  ########
        fx1 = fun(x1)
        vxm.append(np.linalg.norm(s, 1) / np.linalg.norm(x1, 1))

    return x1, it, vxm


def NewtonSysCorde(fun, jac, x0, tolx, tolf, nmax):
    matjac = jac(x0)
    if np.linalg.det(matjac) == 0:  ########
        return None, None, None

    s = -np.linalg.solve(matjac, fun(x0))  ########
    it = 1
    x1 = x0 + s  ########
    fx1 = fun(x1)

    vxm = [np.linalg.norm(s, 1) / np.linalg.norm(x1, 1)]

    while (
        it <= nmax
        and np.linalg.norm(fx1, 1) >= tolf
        and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1)
    ):  ########
        x0 = x1  ########
        it += 1

        if np.linalg.det(matjac) == 0:  ########
            return None, None, None

        s = -np.linalg.solve(matjac, fun(x0))  ########

        x1 = x0 + s  ########
        fx1 = fun(x1)
        vxm.append(np.linalg.norm(s, 1) / np.linalg.norm(x1, 1))

    return x1, it, vxm


def NewtonSysSham(fun, jac, x0, tolx, tolf, nmax):
    matjac = jac(x0)
    if np.linalg.det(matjac) == 0:  ########
        return None, None, None

    s = -np.linalg.solve(matjac, fun(x0))  ########

    it = 1
    x1 = x0 + s  ########
    fx1 = fun(x1)

    vxm = [np.linalg.norm(s, 1) / np.linalg.norm(x1, 1)]
    update = 10

    while (
        it <= nmax
        and np.linalg.norm(fx1, 1) >= tolf
        and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1)
    ):  ########
        x0 = x1  ########
        it += 1
        if it % update == 0:
            matjac = jac(x0)  ########

            if np.linalg.det(matjac) == 0:  ########
                return None, None, None
        else:
            s = -np.linalg.solve(matjac, fun(x0))  ########

        x1 = x0 + s  ########
        fx1 = fun(x1)
        vxm.append(np.linalg.norm(s, 1) / np.linalg.norm(x1, 1))

    return x1, it, vxm


def NewtonMin(fgrad, Hess, x0, tolx, tolf, nmax):
    matHess = Hess(x0)  ########
    if np.linalg.det(matHess) == 0:  ########
        # la matrice non è a rango max
        return None, None, None

    grad_x0 = fgrad(x0)
    s = -np.linalg.solve(matHess, grad_x0)  ########

    it = 1
    x1 = x0 + s  ########
    grad_x1 = fgrad(x1)
    vxk = [np.linalg.norm(s, 1)]

    while (
        it <= nmax
        and np.linalg.norm(grad_x1, 1) >= tolf
        and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1)
    ):  ########
        x0 = x1  ########
        it += 1
        matHess = Hess(x0)  ########
        grad_x0 = grad_x1  ########

        if np.linalg.det(Hess) == 0:
            # la matrice non è a rango max
            return None, None, None

        s = -np.linalg.solve(matHess, grad_x0)  ########
        x1 = x0 + s  ########
        grad_x1 = fgrad(x1)
        vxk.append(np.linalg.norm(s, 1))
    return x1, it, vxk


def NewtonMinMod(fgrad, Hess, x0, tolx, tolf, nmax):
    matHess = np.array(
        [
            [Hess[0, 0](x0[0], x0[1]), Hess[0, 1](x0[0], x0[1])],
            [Hess[1, 0](x0[0], x0[1]), Hess[1, 1](x0[0], x0[1])],
        ]
    )
    grad_x0 = np.array([fgrad[0](x0[0], x0[1]), fgrad[1](x0[0], x0[1])])

    if np.linalg.det(matHess) == 0:  ########
        # la matrice non è a rango max
        return None, None, None

    s = -np.linalg.solve(matHess, grad_x0)  ########

    it = 1
    x1 = x0 + s  ########
    grad_x1 = np.array([fgrad[0](x1[0], x1[1]), fgrad[1](x1[0], x1[1])])
    vxk = [np.linalg.norm(s, 1)]

    while (
        it <= nmax
        and np.linalg.norm(grad_x1, 1) >= tolf
        and np.linalg.norm(s, 1) >= tolx * np.linalg.norm(x1, 1)
    ):  ########
        x0 = x1
        grad_x0 = grad_x1
        it += 1

        matHess = np.array(
            [
                [Hess[0, 0](x0[0], x0[1]), Hess[0, 1](x0[0], x0[1])],
                [Hess[1, 0](x0[0], x0[1]), Hess[1, 1](x0[0], x0[1])],
            ]
        )  ########

        if np.linalg.det(Hess) == 0:
            # la matrice non è a rango max
            return None, None, None

        s = -np.linalg.solve(matHess, grad_x0)  ########
        x1 = x0 + s  ########
        grad_x1 = np.array([fgrad[0](x1[0], x1[1]), fgrad[1](x1[0], x1[1])])

        vxk.append(np.linalg.norm(s, 1))
    return x1, it, vxk


def Jacobi(A, b, x0, toll, itmax):
    err = 1000
    d = np.diag(A)
    n = A.shape[0]
    invM = np.diag(1 / d)
    E = np.tril(A, -1)  ########
    F = np.triu(A, 1)  ########
    N = -(E + F)  ########
    T = np.dot(invM, N)  ########

    autoval = np.linalg.eigvals(T)
    raggiospettrale = np.max(np.abs(autoval))  ########
    it = 0
    verr = []

    while it < itmax and err >= toll:  ########
        x = (b + np.dot(N, x0)) / d.reshape(n, 1)  ########
        err = np.linalg.norm(x - x0) / np.linalg.norm(x)
        verr.append(err)
        x0 = x.copy()
        it += 1
    return x, it, verr


def GaussSeidel(A, b, x0, toll, itmax):
    err = 1000
    d = np.diag(A)

    D = np.diag(d)  ########
    E = np.tril(A, -1)  ########
    F = np.tiru(A, 1)  ########

    M = D + E  ########
    N = -F  ########
    T = np.dot(np.linalg.inv(M), N)  ########

    autoval = np.linalg.eigvals(T)
    raggiospettrale = np.max(np.abs(autoval))  ########

    it = 0
    verr = []

    while it < itmax and err >= toll:  ########
        tmp = b - F @ x0  ########
        x, _ = LSolve(M, tmp)  ########
        err = np.linalg.norm(x - x0) / np.linalg.norm(x)
        verr.append(err)
        x0 = x.copy()
        it += 1

    return x, it, verr


def GaussSeidelSor(A, b, x0, toll, itmax, omega):
    err = 1000
    d = np.diag(A)  ########
    D = np.diag(d)  ########
    Dinv = np.diag(1 / d)  ########
    E = np.tril(A, -1)  ########
    F = np.tril(A, 1)  ########

    Momega = D + omega * E
    Nomega = (1 - omega) * D - omega * F
    T = np.dot(np.linalg.inv(Momega), Nomega)  ########

    autoval = np.linalg.eigvals(T)
    raggiospettrale = np.max(np.abs(autoval))

    M = D + E
    N = -F
    it = 0
    xold = x0.copy()
    xnew = x0.copy()
    verr = []

    while it < itmax and err >= toll:  ########
        tmp = b - np.dot(F, xold)  ########
        xtilde, flag = LSolve(M, tmp)  ########
        xnew = (1 - omega) * xold + omega * xtilde  ########

        err = np.linalg.norm(xnew - xold) / np.linalg.norm(xnew)
        verr.append(err)
        xold = xnew.copy()
        it += 1

    return xnew, it, verr


def SteepestDescent(A, b, x0, itmax, toll):
    n, m = A.shape
    if n != m:
        # matrice non quadrata
        return None, None, None, None

    x = x0
    r = A @ x - b
    p = -r  ########
    it = 0

    nb = np.linalg.norm(b)
    err = np.linalg.norm(r) / nb
    vsol = [x]
    verr = [err]

    while it < itmax and err >= toll:  ########
        it += 1
        Ap = A @ p  ########

        alpha = -(r.T @ p) / (p.T @ Ap)  ########
        x = x + alpha * p  ########
        vsol.append(x)

        r = r + alpha * Ap  ########

        err = np.linalg.norm(r) / nb
        verr.append(err)
        # Direzione opposta alla direzione del gradiente
        p = -r  ########

    return x, verr, vsol, it


def ConjugateGradient(A, b, x0, itmax, toll):
    n, m = A.shape
    if n != m:
        # matrice non quadrata
        return None, None, None, None

    x = x0
    r = A @ x - b
    p = -r
    it = 0

    nb = np.linalg.norm(b)
    err = np.linalg.norm(r) / nb
    vsol = [x]
    verr = [err]

    while it < itmax and err > toll:
        it += 1
        Ap = A @ p  ########

        alpha = -(r.T @ p) / (p.T @ Ap)  ########
        x = x + alpha * p  ########
        vsol.append(x)

        rtr = r.T @ r
        r += alpha * Ap
        gamma = r.T @ r / rtr  ########
        err = np.linalg.norm(r) / nb
        verr.append(err)

        # La nuova direzione appartiene al piano individuato da -r e p. gamma è scelto in maniera tale che la nuova direzione
        # sia coniugata rispetto alla direzione precedente( che geometricamente significa che punti verso il centro)
        p = -r + gamma * p  ########

    return x, verr, vsol, it


def eqnorm(A, b):
    G = A.T @ A  ########
    f = A.T @ b  ########

    L = spl.cholesky(G, lower=True)  ########
    U = L.T  ########

    z, flag = LSolve(L, f)  ########
    if flag == 0:  # optional
        x, flag = USolve(U, z)  ########
    return x


def QRLS(A, b):
    n = A.shape[1]

    Q, R = spl.qr(A)
    h = Q.T @ b  ########

    x, flag = USolve(R[0:n, :], h[0:n])  ########
    residuo = np.linalg.norm(h[n:]) ** 2
    return x, residuo


def SVDLS(A, b):
    m, n = A.shape
    U, s, VT = spl.svd(A)
    V = VT.T
    thresh = np.spacing(1) * m * s[0]
    k = np.count_nonzero(s > thresh)

    d = U.T @ b ########
    d1 = d[:k].reshape(k, 1)  ########
    s1 = s[:k].reshape(k, 1)  ########

    c = d1 / s1  ########
    x = V[:, :k] @ c
    residuo = np.linalg.norm(d[k:]) ** 2

    return x, residuo


def Plagr(xnodi, j):
    xzeri = np.zeros_like(xnodi)
    n = xnodi.size

    if j == 0:
        xzeri = xnodi[1:n]
    else:
        xzeri = np.append(xnodi[0:j], xnodi[j + 1 : n])  ########

    num = np.poly(xzeri)  ########
    den = np.polyval(num, xnodi[j])  ########

    p = num / den
    return p


def InterpL(x, y, xx):
    n = x.size
    m = xx.size
    L = np.zeros((m, n))

    for j in range(n):
        p = Plagr(x, j)  ########
        L[:, j] = np.polyval(p, xx)  ########

    return L @ y
