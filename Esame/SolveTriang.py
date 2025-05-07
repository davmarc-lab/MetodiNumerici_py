import numpy as np

def USolve(U, b):
    """
    Risoluzione con procedura backward di Ux=b con U triangolare superiore
     Input: U matrice triangolare superiore
            b termine noto
    Output: x: soluzione del sistema lineare
            flag=  0, se sono soddisfatti i test di applicabilità
                   1, se non sono soddisfatti
    """
    # Test dimensione
    m, n = U.shape
    flag = 0
    if n != m:
        print("Errore: Matrice non quadrata")
        flag = 1
        x = []
        return x, flag

    # Test singolarità
    if np.all(np.diag(U)) != True:
        print("el. diag. nullo - matrice triangolare superiore")
        x = []
        flag = 1
        return x, flag
    # Preallocazione vettore soluzione
    x = np.zeros((n, 1))

    for i in range(n - 1, -1, -1):
        s = np.dot(U[i, i + 1 : n], x[i + 1 : n])
        x[i] = (b[i] - s) / U[i, i]

    return x, flag


def LSolve(L, b):
    """
    Risoluzione con procedura forward di Lx=b con L triangolare inferiore
     Input: L matrice triangolare inferiore
            b termine noto
    Output: x: soluzione del sistema lineare
            flag=  0, se sono soddisfatti i test di applicabilità
                   1, se non sono soddisfatti
    """
    # test dimensione
    m, n = L.shape
    flag = 0
    if n != m:
        print("Errore: Matrice non quadrata")
        flag = 1
        x = []
        return x, flag

    # Test singolaritò
    if np.all(np.diag(L)) != True:
        print("el. diag. nullo - matrice triangolare inferiore")
        x = []
        flag = 1
        return x, flag

    # Preallocazione vettore soluzione
    x = np.zeros((n, 1))

    for i in range(n):
        s = np.dot(L[i, :i], x[:i])
        x[i] = (b[i] - s) / L[i, i]

    return x, flag
