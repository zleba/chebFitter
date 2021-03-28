from nodes import GetWeights, GetNodes, getPols, getPolsFast, GetCoefsCheb

from numba import njit

import numpy as np

# return values of -2*log(p(x)), where p(x) is normalized to 1 over the fitted range


@njit
def calcLoopFast(Size, dataNorm):
    polSum = np.zeros(Size)

    for xx in dataNorm:
        polSum += getPolsFast(Size, xx)
    return polSum


class chebFitter:

    def getLogFunction(self, pars):
        fVals = np.empty_like(self.nodes)

        # calc function values
        for i in range(fVals.size):
            fVals[i] = self.myFun(self.nodes[i], pars)

        # normalize the function
        Int = np.dot(fVals, self.weights)

        fVals = -2 * np.log(fVals / Int)  # normalize by integral

        return fVals

    # get data transformed into the grid such that (chebFunVals dot dataGrid) == logL
    def getDataGrid(self, data):
        a = self.nodes[0]
        b = self.nodes[-1]

        # normalize between 0 and 1
        dataNorm = (data - a) / (b - a)

        polSum = calcLoopFast(self.nodes.size, dataNorm)
        print("Done")

        # transform to the basis of the cheb nodes
        gridVals = self.coefsMat @ polSum

        return gridVals

    def __init__(self, Size, xMin, xMax, data, fun):

        self.myFun = fun

        self.useCheb = True

        # loading the Cheb nodes
        self.nodes = (xMax-xMin) * GetNodes(Size) + xMin

        # loding the weights for integration
        self.weights = (xMax - xMin) * GetWeights(Size)

        # calculate the transformation matrix from pol coefs to grid points
        self.coefsMat = np.transpose(GetCoefsCheb(Size))

        print("Loading data grid")
        dataSel = data[(xMin < data) & (data < xMax)]
        self.dataGrid = self.getDataGrid(dataSel)
        self.data = dataSel
        # tie(dataGrid, dataGridCov) = getDataGridWithCov();

    # function assumed to be normalized !!!
    def getLogLikelihoodSlow(self, pars):

        L = 0
        # pragma omp parallel for  reduction(+: L)
        for d in self.data:
            v = self.myFun(d, pars)
            # cout << v << endl;
            L += -2*np.log(v)

        return L

    # evaluation using cheb pols
    def getLogLikelihoodFast(self, pars):

        funVals = self.getLogFunction(pars)
        LL = np.dot(funVals, self.dataGrid)

        return LL

    def eval(self, pars):

        if self.useCheb:
            return self.getLogLikelihoodFast(pars)
        else:
            return self.getLogLikelihoodSlow(pars)

    def evalVec(self, pars):
        p = dict(zip(self.parsNames, pars))
        return self.eval(p)

    def fitData(self, pars, limits=None, useCheb=True):

        parsVec = []
        limitsVec = []
        self.parsNames = []
        for k in pars:
            self.parsNames.append(k)
            parsVec.append(pars[k])
            if limits is not None and k in limits:
                limitsVec.append(limits[k])
            else:
                limitsVec.append((None, None))

        self.useCheb = useCheb

        from iminuit import Minuit

        m = Minuit(self.evalVec, parsVec, name=self.parsNames)
        m.limits = limitsVec
        m.migrad()  # run optimiser
        for i in range(len(m.values)):
            print(np.sqrt(m.covariance[i,i]))
     
        return dict([[self.parsNames[i], (m.values[i], m.errors[i]) ] for i in range(len(m.values))])
    
    def funFast(self, x, pars):
        if not hasattr(self, 'parsOld') or self.parsOld != pars:
            funVals = self.getLogFunction(pars)
            self.coefs = np.transpose(self.coefsMat) @ funVals
            self.parsOld = pars

        a = self.nodes[0]
        b = self.nodes[-1]

        xRel = (x - a) / (b - a)

        pols = getPols(len(self.nodes), xRel)
        res = np.dot(self.coefs, pols)

        return np.exp(-0.5*res)
