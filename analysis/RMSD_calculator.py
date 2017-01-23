import numpy as np
import pickle


def center(X):
        '''Expects n x 3 geometry'''
        return np.mean(X, axis = 0)

def focus(X, Y):
    '''Receives 2 Fragment objects and centers them around 0'''
    Xc = center(X)
    Yc = center(Y)

    X = X - Xc
    Y = Y - Yc

    result = [X, Y, Xc, Yc]
    return result

def rotate_pair(X, Y):

    # Put them in 0, 0, 0
    X, Y, Xc, Yc = focus(X, Y)

    # Received as 5 x 3
    
    #R = Y.T * X # 3 x 5 * 5 x 3 = 3 x 3
    R = np.dot(Y.T, X) # 3 x 5 * 5 x 3 = 3 x 3

    V, S, Wt = np.linalg.svd(R)


    #Z = np.diag([1, 1, -1])

    U = Wt.T * V.T

    if np.linalg.det(U) < 0:
        #print "Reflection catch"
        Wt[2] = -Wt[2]
        U = Wt.T * V.T

    # Rotate Y using U
    Yr = np.dot(Y, U.T) # + f1_center  # 5 x 3 * 3 x 3 = 5 x 3


    result = [X, Yr, U.T, S]
    return result



def compute_RMSD(X, Y, S):

    # Received as n x 3

    # Transpose to 3 x 5
    X = X.T
    Y = Y.T

    n = float(X.shape[1])

    E0 = np.sum(np.square(np.linalg.norm(X, axis = 0)) + np.square(np.linalg.norm(Y, axis = 0)))

    RMSD = np.sqrt((1 / n) * (E0 - 2 * np.sum(S)))

    return RMSD


def RMSD_func(X, Y):
    X = X.T
    Y = Y.T
    Xc = np.mean(X, 1)
    Yc = np.mean(Y, 1)
    Xt = X.T
    R = np.dot(Y, Xt)
    V, S, W = np.linalg.svd(R)
    d = [ 1, 1, -1 ]
    Z = np.diag(d)
    U = np.dot(W, V.T)
    det = np.linalg.det(U)
    if det == -1:
        U = np.dot(U, Z)
    RMSD = np.sqrt(np.sum(abs((Xc ** 2 - np.dot(U, Yc) ** 2 - 2 * Xc * Yc)), 0))
    return RMSD


#a = np.asmatrix(5 * np.random.random_sample((5, 3)) - 5)
#b = np.asmatrix(5 * np.random.random_sample((5, 3)) - 5)
#
#X, Yr, Ut, S = rotate_pair(a, b)
#
#print compute_RMSD(X, Yr, S)
#
#A = np.squeeze(np.asarray(X))
#B = np.squeeze(np.asarray(Yr))
#
#print A
#print X
#
#print B
#print Yr
#
#print RMSD_func(A, B)


X = np.load("X.npy")
Y = np.load("Y.npy")

print X
print Y

X, Yr, Ut, S = rotate_pair(X, Y)
print compute_RMSD(X, Yr, S)



















