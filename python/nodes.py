import numpy as np
from numba import njit


def GetWeights(Size):
    N = Size - 1
    assert(N % 2 == 0)

    coef = np.empty([Size, Size])

    for k in range(N // 2 + 1):
        coef[2*k, N] = 1./N
        coef[2*k, 0] = 1./N

        coef[2*k, N//2] = 2./N * (2 * ((k+1) % 2) - 1)

        for n in range(1, N // 2):
            coef[2*k, n] = coef[2*k, N-n] = 2./N * np.cos(n*k*np.pi*2/N)

    wgt = np.zeros([Size])

    for i in range(Size):
        wgt[i] += coef[0, i]
        wgt[i] += coef[N, i] / (1. - N*N)
        for k in range(1, N // 2):
            w = 2./(1 - 4*k*k)
            wgt[i] += w*coef[2*k, i]

        wgt[i] *= 0.5  # for interval (0,1)

    return wgt


# Get vector with nodes by definition between 0 and 1
def GetNodes(Size):
    assert((Size - 1) % 2 == 0)

    x = np.linspace(0, np.pi, Size)
    xi = (1 - np.cos(x)) / 2

    return xi


# Evaluate cheb pols to Size at point x, x is between 0 and 1
def getPols(Size, x):
    if np.size(x) == 1:
        assert(x >= 0 and x <= 1)

    C = 2*(2*x-1)

    pol = np.empty((Size, np.size(x)))
    pol = np.squeeze(pol)

    if(Size >= 1):
        pol[0] = 1
    if(Size >= 2):
        pol[1] = C/2

    for i in range(2, Size):
        pol[i] = C*pol[i-1] - pol[i-2]
    return pol


@njit
def getPolsFast(Size, x):
    assert(x >= 0 and x <= 1)
    assert(Size >= 2)

    C = 2*(2*x-1)

    pol = np.empty(Size)

    pol[0] = 1
    pol[1] = C/2

    for i in range(2, Size):
        pol[i] = C*pol[i-1] - pol[i-2]
    return pol


# Evaluate sum of cheb pols to Size at vector x, x els are between 0 and 1
def getPolsSum(Size, x):
    assert(Size > 2)

    polSum = np.empty([Size])

    pol0 = 0*x + 1
    pol1 = 2*x - 1
    C    = 2 * pol1

    for i in range(2, Size):
        polSum[i-2] = pol0.sum()

        pol2 = C * pol1 - pol0

        pol0 = pol1
        pol1 = pol2

    polSum[Size-2] = pol0.sum()
    polSum[Size-1] = pol1.sum()

    return polSum


# Transformation matrix between cheb. nodes and cheb. coeficients
def GetCoefs(oldSize, isInverse=False):
    N = oldSize - 1
    assert(N % 2 == 0)

    coef = np.empty([oldSize, oldSize])

    mul = 1
    C = 1./N
    if isInverse:
        C = 1./2

    # isInverse = false;
    for k in range(N+1):
        if (not isInverse):
            coef[k, N] = C
            coef[k, 0] = C * (-1 if k % 2 == 1 else 1)

        else:
            mul = -1 if k % 2 == 1 else 1
            coef[N-k, N] = C * mul
            coef[N-k, 0] = C

        for n in range(1, N):
            el = np.cos(n*k*np.pi / N) * 2.*C * mul
            if(not isInverse):
                coef[k, N-n] = el
            else:
                coef[N-k, N-n] = el

    return coef


# with better normalization of the borders
def GetCoefsCheb(oldSize):
    coef = GetCoefs(oldSize)

    coef[0, :] *= 0.5
    coef[-1, :] *= 0.5

    return coef


# Evaluate Cheb. pol at point x
def evalPol(polCoef, x):
    assert(x >= 0 and x <= 1)
    pols = getPols(len(polCoef), x)
    return np.dot(pols, polCoef)

#
#//Get Interpolation vector at point x
#VectorXd interpol(const VectorXd &xi, double x)
#{
#    double Norm = (xi[xi.size()-1] - xi[0])/2;
#    VectorXd coefs(xi.size());
#    for(int i = 0; i < xi.size(); ++i) {
#        double num = 1, den = 1;
#        for(int j = 0; j < xi.size(); ++j)
#            if(j != i) {
#                num *= (x     - xi(j))/Norm;
#                den *= (xi(i) - xi(j))/Norm;
#                //cout << x  <<" "<< xi(i) <<" "<< xi(j) <<" : " << num <<" "<< den <<  endl;
#            }
#        //cout << num <<" "<< den << endl;
#        coefs(i) = num/den;
#    }
#    return coefs;
#}
#
#
#//Get interpolated function value at point x
#double interpol(VectorXd xi, VectorXd vals, double x)
#{
#    VectorXd coefs = interpol(xi, x);
#    //cout << "Radecek " << coefs << endl;
#    return coefs.dot(vals);
#}


