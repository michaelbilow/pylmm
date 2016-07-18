from scipy import linalg
import numpy as np


def matrixMult(A, B):
    # If there is no fblas then we will revert to np.dot()
    try:
        linalg.fblas
    except AttributeError:
        return np.dot(A, B)

    # If the matrices are in Fortran order then the computations will be faster
    # when using dgemm.  Otherwise, the function will copy the matrix and that takes time.
    if not A.flags['F_CONTIGUOUS']:
        AA = A.T
        transA = True
    else:
        AA = A
        transA = False

    if not B.flags['F_CONTIGUOUS']:
        BB = B.T
        transB = True
    else:
        BB = B
        transB = False

    return linalg.fblas.dgemm(alpha=1., a=AA, b=BB, trans_a=transA, trans_b=transB)


def calculateKinship(W, center=False):
    """
         W is an n x m matrix encoding SNP minor alleles.

         This function takes a matrix oF SNPs, imputes missing values with the maf,
         normalizes the resulting vectors and returns the RRM matrix.
      """
    n = W.shape[0]
    m = W.shape[1]
    keep = []
    for i in range(m):
        mn = W[True - np.isnan(W[:, i]), i].mean()
        W[np.isnan(W[:, i]), i] = mn
        vr = W[:, i].var()
        if vr == 0: continue

        keep.append(i)
        W[:, i] = (W[:, i] - mn) / np.sqrt(vr)

    W = W[:, keep]
    K = matrixMult(W, W.T) * 1.0 / float(m)
    if center:
        P = np.diag(np.repeat(1, n)) - 1 / float(n) * np.ones((n, n))
        S = np.trace(matrixMult(matrixMult(P, K), P))
        K_n = (n - 1) * K / S
        return K_n
    return K


def GWAS(Y, X, K, Kva=[], Kve=[], X0=None, REML=True, refit=False):
    """

         Performs a basic GWAS scan using the LMM.  This function
         uses the LMM module to assess association at each SNP and
         does some simple cleanup, such as removing missing individuals
         per SNP and re-computing the eigen-decomp

         Y - n x 1 phenotype vector
         X - n x m SNP matrix
         K - n x n kinship matrix
         Kva,Kve = linalg.eigh(K) - or the eigen vectors and values for K
         X0 - n x q covariate matrix
         REML - use restricted maximum likelihood
         refit - refit the variance component for each SNP

      """
    n = X.shape[0]
    m = X.shape[1]

    if X0 == None:
        X0 = np.ones((n, 1))

    # Remove missing values in Y and adjust associated parameters
    v = np.isnan(Y)
    if v.sum():
        keep = True - v
        keep = keep.reshape((-1,))
        Y = Y[keep]
        X = X[keep, :]
        X0 = X0[keep, :]
        K = K[keep, :][:, keep]
        Kva = []
        Kve = []

    if len(Y) == 0:
        return np.ones(m) * np.nan, np.ones(m) * np.nan

    L = LMM(Y, K, Kva, Kve, X0)
    if not refit: L.fit()

    PS = []
    TS = []

    n = X.shape[0]
    m = X.shape[1]

    for i in range(m):
        x = X[:, i].reshape((n, 1))
        v = np.isnan(x).reshape((-1,))
        if v.sum():
            keep = True - v
            xs = x[keep, :]
            if xs.var() == 0:
                PS.append(np.nan)
                TS.append(np.nan)
                continue

            Ys = Y[keep]
            X0s = X0[keep, :]
            Ks = K[keep, :][:, keep]
            Ls = LMM(Ys, Ks, X0=X0s)
            if refit:
                Ls.fit(X=xs)
            else:
                Ls.fit()
            ts, ps = Ls.association(xs, REML=REML)
        else:
            if x.var() == 0:
                PS.append(np.nan)
                TS.append(np.nan)
                continue

            if refit:
                L.fit(X=x)
            ts, ps = L.association(x, REML=REML)

        PS.append(ps)
        TS.append(ts)

    return TS, PS