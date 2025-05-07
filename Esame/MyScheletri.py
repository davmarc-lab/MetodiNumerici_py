import math
from math import copysign
import numpy as np
import scipy.linalg as spl

from SolveTriang import LSolve
from SolveTriang import USolve


def sign(x):
    return copysign(1, x)


def Bisezione(fname, a, b, tolx):
    fa = fname(a)
    fb = fname(b)

    if:
        # Non è possibile applicare il metodo della bisezione
        return None, None, None

    it = 0
    vxk = []

    maxit = math.ceil(math.log((b - a) / tolx) / math.log(2)) - 1

    while : 
        xk = 
        vxk.append(xk)
        it += 1
        fxk = fname(xk)
        if fxk == 0:
            return xk, it, vxk

        if sign(fa) * sign(fxk) > 0:
            a = 
            fa = 
        elif sign(fb) * sign(fxk) > 0:
            b = 
            fb = 

    return xk, it, vxk


def Falsi(fname, a, b, tolx, tolf, nmax):
    fa = fname(a)
    fb = fname(b)

    if :
        # Non è possibile applicare il metodo della falsa posizione
        return None, None, None

    it = 0
    vxk = []
    fxk = 10

    while :  
        xk = 
        vxk.append(xk)
        it += 1
        fxk = fname(xk)

        if fxk == 0:
            return xk, it, vxk

        if sign(fa) * sign(fxk) > 0:
            a = 
            fa = 
        elif sign(fb) * sign(fxk) > 0:
            b = 
            fb = 

    return xk, it, vxk


def Corde(fname, m, x0, tolx, tolf, nmax):
    vxk = []
    fx0 = 
    d = 
    x1 = 
    fx1 = fname(x1)
    vxk.append(x1)
    it = 1

    while :
        x0 = 
        fx0 = 
        d = 
        """
        x1= ascissa del punto di intersezione tra  la retta che passa per il punto
        (xi,f(xi)) e ha pendenza uguale a m  e l'asse x
        """
        x1 = 
        fx1 = fname(x1)
        it += 1
        vxk.append(x1)

        if it == nmax:
            # Max Num di iterazioni
            return x1, it, vxk
    return x1, it, vxk


def Secanti(fname, xm1, x0, tolx, tolf, nmax):
    vxk = []
    fxm1 = 
    fx0 = 
    d = 
    x1 = 
    fx1 = fname(x1)
    vxk.append(x1)
    it = 1

    while :
        xm1 = 
        x0 = 
        fxm1 = 
        fx0 = 
        d = 
        x1 = 

        fx1 = fname(x1)
        vxk.append(x1)
        it += 1

    return x1, it, vxk


def Newton(fname, fpname, x0, tolx, tolf, nmax):
    vxk = []
    fx0 = fname(x0)
    fpx0 = fpname(x0)
    if :
        # derivata prima nulla
        return None, None, None

    d = 
    x1 = 

    fx1 = fname(x1)
    vxk.append(x1)
    it = 1

    while :  
        x0 = 
        fx0 = 
        fpx0 = 
        if : 
            # derivata prima nulla
            return None, None, None

        d = 
        x1 = 
        fx1 = fname(x1)
        it += 1

        vxk.append(x1)

    return x1, it, vxk


def NewtonMod(fname, fpname, m, x0, tolx, tolf, nmax):
    vxk = []
    fx0 = 
    fpx0 = 
    if :  
        # derivata prima nulla in x0
        return None, None, None

    d = 
    x1 = 

    fx1 = 
    vxk.append(x1)
    it = 1

    while : 
        x0 = 
        fx0 = 
        fpx0 = fpname(x0)
        if : 
            # derivata prima nulla in x0
            return None, None, None
        d = 
        """
           x1= ascissa del punto di intersezione tra  la retta che passa per il punto
           (xi,f(xi)) ed è tangente alla funzione f(x) nel punto (xi.f(xi))  e l'asse x
        """
        x1 = 
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
    if : 
        return None, None, None

    s = 

    it = 1
    x1 = 
    fx1 = fun(x1)
    vxm = [np.linalg.norm(s, 1) / np.linalg.norm(x1, 1)]

    while : 
        x0 = 
        it += 1

        matjac = jac(x0)
        if :
            return None, None, None

        s = 

        x1 = 
        fx1 = fun(x1)
        vxm.append(np.linalg.norm(s, 1) / np.linalg.norm(x1, 1))

    return x1, it, vxm


def NewtonSysCorde(fun, jac, x0, tolx, tolf, nmax):
    matjac = jac(x0)
    if :
        return None, None, None

    s = 
    it = 1
    x1 = 
    fx1 = fun(x1)

    vxm = [np.linalg.norm(s, 1) / np.linalg.norm(x1, 1)]

    while :
        x0 = 
        it += 1

        if : 
            return None, None, None

        s = 

        x1 = 
        fx1 = fun(x1)
        vxm.append(np.linalg.norm(s, 1) / np.linalg.norm(x1, 1))

    return x1, it, vxm


def NewtonSysSham(fun, jac, x0, tolx, tolf, nmax):
    matjac = jac(x0)
    if :
        return None, None, None

    s = 

    it = 1
    x1 = 
    fx1 = fun(x1)

    vxm = [np.linalg.norm(s, 1) / np.linalg.norm(x1, 1)]
    update = 10

    while :
        x0 = 
        it += 1
        if it % update == 0:
            matjac = 

            if :
                return None, None, None
        else:
            s = 

        x1 = 
        fx1 = fun(x1)
        vxm.append(np.linalg.norm(s, 1) / np.linalg.norm(x1, 1))

    return x1, it, vxm


def NewtonMin(fgrad, Hess, x0, tolx, tolf, nmax):
    matHess = 
    if :
        # la matrice non è a rango max
        return None, None, None

    grad_x0 = fgrad(x0)
    s = 

    it = 1
    x1 =
    grad_x1 = fgrad(x1)
    vxk = [np.linalg.norm(s, 1)]

    while :
        x0 = 
        it += 1
        matHess =
        grad_x0 = 

        if np.linalg.det(Hess) == 0:
            # la matrice non è a rango max
            return None, None, None

        s = 
        x1 =
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

    if :
        # la matrice non è a rango max
        return None, None, None

    s = 

    it = 1
    x1 =
    grad_x1 = np.array([fgrad[0](x1[0], x1[1]), fgrad[1](x1[0], x1[1])])
    vxk = [np.linalg.norm(s, 1)]

    while :
        x0 = x1
        grad_x0 = grad_x1
        it += 1

        matHess = np.array(
            [
                [Hess[0, 0](x0[0], x0[1]), Hess[0, 1](x0[0], x0[1])],
                [Hess[1, 0](x0[0], x0[1]), Hess[1, 1](x0[0], x0[1])],
            ]
        ) 

        if np.linalg.det(Hess) == 0:
            # la matrice non è a rango max
            return None, None, None

        s = 
        x1 = 
        grad_x1 = np.array([fgrad[0](x1[0], x1[1]), fgrad[1](x1[0], x1[1])])

        vxk.append(np.linalg.norm(s, 1))
    return x1, it, vxk


def Jacobi(A, b, x0, toll, itmax):
    err = 1000
    d = np.diag(A)
    n = A.shape[0]
    invM = np.diag(1 / d)
    E = 
    F = 
    N = 
    T = 

    autoval = np.linalg.eigvals(T)
    raggiospettrale = 
    it = 0
    verr = []

    while :
        x =
        err = np.linalg.norm(x - x0) / np.linalg.norm(x)
        verr.append(err)
        x0 = x.copy()
        it += 1
    return x, it, verr


def GaussSeidel(A, b, x0, toll, itmax):
    err = 1000
    d = np.diag(A)

    D = 
    E = 
    F = 

    M = 
    N = 
    T = 

    autoval = np.linalg.eigvals(T)
    raggiospettrale = 

    it = 0
    verr = []

    while :
        tmp = 
        x, _ = 
        err = np.linalg.norm(x - x0) / np.linalg.norm(x)
        verr.append(err)
        x0 = x.copy()
        it += 1

    return x, it, verr


def GaussSeidelSor(A, b, x0, toll, itmax, omega):
    err = 1000
    d = 
    D = 
    Dinv =
    E = 
    F = 

    Momega = D + omega * E
    Nomega = (1 - omega) * D - omega * F
    T = 

    autoval = np.linalg.eigvals(T)
    raggiospettrale = np.max(np.abs(autoval))

    M = D + E
    N = -F
    it = 0
    xold = x0.copy()
    xnew = x0.copy()
    verr = []

    while :
        tmp = 
        xtilde, flag = 
        xnew = 

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
    p =
    it = 0

    nb = np.linalg.norm(b)
    err = np.linalg.norm(r) / nb
    vsol = [x]
    verr = [err]

    while :
        it += 1
        Ap = 

        alpha = 
        x = 
        vsol.append(x)

        r =

        err = np.linalg.norm(r) / nb
        verr.append(err)
        # Direzione opposta alla direzione del gradiente
        p = 

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
        Ap = 

        alpha =
        x = 
        vsol.append(x)

        rtr = r.T @ r
        r += alpha * Ap
        gamma = 
        err = np.linalg.norm(r) / nb
        verr.append(err)

        # La nuova direzione appartiene al piano individuato da -r e p. gamma è scelto in maniera tale che la nuova direzione
        # sia coniugata rispetto alla direzione precedente( che geometricamente significa che punti verso il centro)
        p = 

    return x, verr, vsol, it


def eqnorm(A, b):
    G = 
    f = 

    L = 
    U = 

    z, flag =
    if flag == 0:  # optional
        x, flag = 
    return x


def QRLS(A, b):
    n = A.shape[1]

    Q, R = spl.qr(A)
    h = 

    x, flag = 
    residuo = np.linalg.norm(h[n:]) ** 2
    return x, residuo


def SVDLS(A, b):
    m, n = A.shape
    U, s, VT = spl.svd(A)
    V = VT.T
    thresh = np.spacing(1) * m * s[0]
    k = np.count_nonzero(s > thresh)

    d =
    d1 = 
    s1 = 

    c = 
    x = V[:, :k] @ c
    residuo = np.linalg.norm(d[k:]) ** 2

    return x, residuo


def Plagr(xnodi, j):
    xzeri = np.zeros_like(xnodi)
    n = xnodi.size

    if j == 0:
        xzeri = xnodi[1:n]
    else:
        xzeri = 

    num = 
    den = 

    p = num / den
    return p


def InterpL(x, y, xx):
    n = x.size
    m = xx.size
    L = np.zeros((m, n))

    for j in range(n):
        p = 
        L[:, j] = 

    return L @ y
