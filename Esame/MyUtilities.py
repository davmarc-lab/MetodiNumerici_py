import numpy as np
import numpy.linalg as npl


def norm1(x):
    return np.max(np.sum(np.abs(x), axis=0))


def norm2(A):
    M = np.dot(A.T, A)
    autoval = np.linalg.eigvals(M)
    return np.sqrt(np.max(np.abs(autoval)))


def norminf(x):
    return np.max(np.sum(np.abs(x), axis=1))


def cond(A):
    return norminf(A) * norminf(np.linalg.inv(A))


def errrel(x, px):
    return npl.norm(px - x, np.inf) / npl.norm(x, np.inf)

# Approssimazione ai minimi quadrati di una funzione
"""
 - x è il dominio della funzione
 - y sono tutte le `f(x)`
 - A è la matrice di vander con le `x`
 - b = y
 - n è il grado del polinomio (retta) di regressione
 - n1 sono i gradi di libertà
"""

n = 1
m = x.shape[0]
n1 = n + 1
A = np.vander(x, increasing=True)[:, :n1]

# Basta un metodo di risoluzione del sistema `Ax=b`
xe = eqnorm(A, y)
re = npl.norm(A @ xe - y.reshape(m, 1)) ** 2
xq, rq = QRLS(A, y)
xs, rs = SVDLS(A, y)

print("Residuo EQN: ", re)
print("Residuo QRLS: ", rq)
print("Residuo SVDLS: ", rs)

# Trovo le coordinate `x` del polinomio di regressione sfruttando il dominio della funzione stessa
xv = np.linspace(np.min(x), np.max(x), 100)

polEQ = np.polyval(np.flip(xe), xv)
polQRLS = np.polyval(np.flip(xq), xv)
polSVDLS = np.polyval(np.flip(xs), xv)

plt.plot(x, y, 'ro', xv, polEQ, xv, polQRLS, xv, polSVDLS)
plt.show()