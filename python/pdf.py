#!/usr/bin/env python3

from numba import njit
import numpy as np

# Nodes & weights for integration
nodes = np.array([ -0.991455371120813, -0.949107912342759, -0.864864423359769, -0.741531185599394, -0.586087235467691, -0.405845151377397, -0.207784955007898, 0.000000000000000, +0.207784955007898, +0.405845151377397, +0.586087235467691, +0.741531185599394, +0.864864423359769, +0.949107912342759, +0.991455371120813])
weights = np.array([ 0.022935322010529,0.063092092629979,0.104790010322250,0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,0.063092092629979,0.022935322010529])


# Exponential
@njit
def E(p, x):
    if (x - p[0]) * p[1].real >= 0:
        r = abs(p[1].real) * np.exp(-p[1]*(x-p[0]))
    else:
        r = 0
    return r

# Gauss
@njit
def G(p, x):
    r = 1/(np.sqrt(2*np.pi) * p[1]) * np.exp( -1./2 * ((x-p[0])/p[1])**2)
    return r


# e = [pos, slope] slope>0 -> tail goes to +infty
# Exp x Exp convolution
@njit
def EE_conv(e1, e2, x):
    r = e2[1].real/(e2[1] - e1[1])*E([e1[0]+e2[0], e1[1]], x) - e1[1].real/(e2[1] - e1[1])*E([e1[0]+e2[0], e2[1]], x)
    return r

# Gaus x Exp convolution
#@njit
def GE_conv(g, e, x):
    m1     = g[0]
    sigma  = g[1]
    m2     = e[0]
    c      = e[1]

    #from scipy.special import erfc
    from math import erfc
    S = 1 if c.real > 0 else -1
    r = abs(c.real) * np.exp(c*(m1+m2) + (sigma**2 * c**2)/2 - c*x) * 0.5 * erfc(S*( (m1+m2) + c*sigma**2 - x)/(np.sqrt(2) * sigma))

    return r

# Gaus x Exp x Exp convolution
#@njit
def GEE_conv(g, e1, e2, x):
    r = e2[1].real/(e2[1] - e1[1])*GE_conv(g, [e1[0]+e2[0], e1[1]], x) - e1[1].real/(e2[1] - e1[1])*GE_conv(g, [e1[0]+e2[0], e2[1]], x)
    return r

# Resolution function (Guass with exp tails twice)
def resFun(x, pars):
    fTMains = pars["fTMains"]
    fTR     = pars["fTR"]
    sigmas  = pars["sigmas"]
    cLMs    = pars["cLMs"]
    cRMs    = pars["cRMs"]
    mus     = pars["mus"]
    
    bigSigmaFrac = pars['bigSigmaFrac']
    fTBigs = pars['fTBigs']
    bigSigmaScale = pars['bigSigmaScale']
    cLBs = pars['cLBs']
    cRBs = pars['cRBs']
    
    ret = 0
    
    A = np.array

    # small gauss
    ret += (1-bigSigmaFrac) * (1-fTMains) * G(A([mus, sigmas]), x)
    ret += (1-bigSigmaFrac) * fTMains*(1.-fTR)* GE_conv([mus,sigmas], [mus, -cLMs], x) 
    ret += (1-bigSigmaFrac) * fTMains*(fTR)* GE_conv([mus,sigmas], [mus, +cRMs], x)
    # big gauss
    ret += bigSigmaFrac * (1-fTBigs) * G(A([mus, bigSigmaScale*sigmas]), x)
    ret += bigSigmaFrac * fTBigs*(1.-fTR)* GE_conv([mus,bigSigmaScale*sigmas], [mus, -cLBs], x)
    ret += bigSigmaFrac * fTBigs*(fTR)* GE_conv([mus,bigSigmaScale*sigmas], [mus, +cRBs], x)
    
    return ret


# Resolution function convoluted with (complex) exp
def FunTemp(x, pars, c):
    fTMains = pars["fTMains"]
    fTR     = pars["fTR"]
    sigmas  = pars["sigmas"]
    cLMs    = pars["cLMs"]
    cRMs    = pars["cRMs"]
    mus     = pars["mus"]
    
    bigSigmaFrac = pars['bigSigmaFrac']
    fTBigs = pars['fTBigs']
    bigSigmaScale = pars['bigSigmaScale']
    cLBs = pars['cLBs']
    cRBs = pars['cRBs']
    
    ret = 0
    
    # small gauss
    ret += (1-bigSigmaFrac) * (1-fTMains) * GE_conv([mus, sigmas], [0,c], x)
    ret += (1-bigSigmaFrac) * fTMains*(1.-fTR)* GEE_conv([mus,sigmas], [mus, -cLMs], [0, c], x) 
    ret += (1-bigSigmaFrac) * fTMains*(fTR)* GEE_conv([mus,sigmas], [mus, +cRMs], [0, c], x)
    # big gauss
    ret += bigSigmaFrac * (1-fTBigs) * GE_conv([mus, bigSigmaScale*sigmas], [0, c], x)
    ret += bigSigmaFrac * fTBigs*(1.-fTR)* GEE_conv([mus,bigSigmaScale*sigmas], [mus, -cLBs], [0, c], x)
    ret += bigSigmaFrac * fTBigs*(fTR)* GEE_conv([mus,bigSigmaScale*sigmas], [mus, +cRBs], [0, c], x)
    
    return ret

# Decay funciton with theta-dep
def funDecayTheta(x, pars, tau=1.520, theta=np.pi/2, K=0.22359):
    #K = 0.22534  # currently hard-coded
    #K = 0.22359  # using nominal BoostVector and Y4S mass
    #tau = 1.520  # ps

    
    cp = +1./tau * 1./(1 + K*np.cos(theta))
    cm = -1./tau * 1./(1 - K*np.cos(theta))

    r = 1./(2*tau) * (1./abs(cp)*FunTemp(x, pars, cp) + 1./abs(cm)*FunTemp(x, pars, cm))
    return r


# Decay funciton (theta integrated)
def funDecay(x, pars, tau=1.520, K=0.22359):

    res = 0
    for n,w in zip(nodes, weights): #integral between cosTh = -1 and 1
        res += w*funDecayTheta(x, pars, tau=tau, theta=np.arccos(n), K=K) * 3./4 * (1 - n*n)

    return res


# normalized to 1./(1 + dm^2 tau^2) : complex
def funOscTheta(x, pars, tau=1.520, dm=0.510, theta=np.pi/2):
    K = 0.22534  # currently hard-coded
    #tau = 1.520  # ps
    #dm = 0.510   # ps-1

    cp = +1./tau * 1./(1 + K*np.cos(theta)) * (+1 - 1j*dm * tau)
    cm = +1./tau * 1./(1 - K*np.cos(theta)) * (-1 - 1j*dm * tau)

    A  = 1./(2*tau) * 1./(1 + 1j*dm*K*tau*np.cos(theta))

    r = A * (1/abs(cp.real)*FunTemp(x, pars, cp) + 1/abs(cm.real)*FunTemp(x, pars, cm))
    return r


# normalized to 1./(1 + dm^2 tau^2) : complex
def funOsc(x, pars, tau=1.520, dm=0.510):

    res = 0
    for n, w in zip(nodes, weights):
        res += w*funOscTheta(x, pars, tau=tau, dm=dm, theta=np.arccos(n)) * 3./4 * (1 - n*n)

    return res





# beta = 276.11e-3
# c = 299.792458 # um/ps
# bgc = beta/(1 - beta**2)**0.5 * c
#
# pars = {
#     'mus': 2.5998222326911318 / bgc,
#     'sigmas': 19.488840799952847 / bgc,
#     'fTMains': 0.6879121782503139,
#     'fTR': 0.5672982488332009,
#     'cLMs': 0.03847980216578396 * bgc,
#     'cRMs': 0.02771982014064649 * bgc,
#     'bigSigmaScale': 6.947680846616044,
#     'bigSigmaFrac': 0*0.16756917570490554,
#     'fTBigs': 0.6268934047143494,
#     'cLBs': 0.0003812636576085269 * bgc,
#     'cRBs': 0.006550526360797996 * bgc
# }


# a = 0.10
# def numConv(f1, f2, x):
#    s = 0.
#    L = 1000.
#    N = int(1e6)
#    for y in np.linspace(-L, L, N):
#        s += f1(y) * f2(x-y)
#    return s * (2*L / N)
#
# for x in np.linspace(-25, 25, 5):
#    #print(x, resFun(x, pars))
#    m1 = 3.
#    m2 = 1.
#    a = -0.1j
#    f1 = lambda y : E(np.array([m1, +0.2]), y)
#    f2 = lambda y : E(np.array([m2, +0.1+a]), y)
#    print(numConv(f1, f2, x), EE_conv([m1, +0.2], [m2, +0.1+a], x) )
