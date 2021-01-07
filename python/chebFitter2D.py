from nodes import GetWeights, GetNodes, getPols, getPolsFast, GetCoefsCheb
from numba import njit
import numpy as np

# Hardcoded outer product for the two matrices (I haven't found the numpy version)
@njit
def outer2D(a, b):
    sa = a.shape[0]
    sb = b.shape[0]

    m = np.empty((sa*sb, sa*sb))

    for i1 in range(sa):
        for j1 in range(sa):
            for i2 in range(sb):
                for j2 in range(sb):
                    i = i1*sb + i2
                    j = j1*sb + j2
                    m[i, j] = a[i1, j1] * b[i2, j2]
    return m


@njit
def calcLoopFast(SizeX, SizeY, dataNorm):
    polSum = np.zeros((SizeX, SizeY))

    for xx in dataNorm:
        rX = getPolsFast(SizeX, xx[0])
        rY = getPolsFast(SizeY, xx[1])
        # polSum += np.outer(rX, rY)
        for i in range(SizeX):
            for j in range(SizeY):
                polSum[i, j] += rX[i] * rY[j]

    return polSum.flatten()


class chebFitter:

    def getLogFunction(self, pars):
        fVals = np.empty(self.nodesX.size * self.nodesY.size)

        # calc function values
        for i in range(self.nodesX.size):
            for j in range(self.nodesY.size):
                fVals[i*self.nodesY.size+j] = self.myFun([self.nodesX[i], self.nodesY[j]], pars)

        # normalize the function
        Int = np.dot(fVals, self.weights)

        fVals = -2 * np.log(fVals / Int)  # normalize by integral

        return fVals

    # get data transformed into the grid such that (chebFunVals dot dataGrid) == logL
    def getDataGrid(self, data):
        xMin = self.nodesX[0]
        xMax = self.nodesX[-1]
        yMin = self.nodesY[0]
        yMax = self.nodesY[-1]

        # normalize between 0 and 1
        dataNorm = np.empty_like(data)
        dataNorm[:, 0] = (data[:, 0] - xMin) / (xMax - xMin)
        dataNorm[:, 1] = (data[:, 1] - yMin) / (yMax - yMin)

        polSum = calcLoopFast(self.nodesX.size, self.nodesY.size, dataNorm)
        print("Done")

        print(self.coefsMat.shape)
        print(polSum.shape)

        # transform to the basis of the cheb nodes
        gridVals = self.coefsMat @ polSum

        return gridVals

    def __init__(self, SizeX, xMin, xMax, SizeY, yMin, yMax, data, fun):

        self.myFun = fun

        self.useCheb = True

        # loading the Cheb nodes
        self.nodesX = (xMax-xMin) * GetNodes(SizeX) + xMin
        self.nodesY = (yMax-yMin) * GetNodes(SizeY) + yMin

        # loding the weights for integration
        self.weightsX = (xMax - xMin) * GetWeights(SizeX)
        self.weightsY = (yMax - yMin) * GetWeights(SizeY)
        self.weights  = np.outer(self.weightsX, self.weightsY).flatten()

        # calculate the transformation matrix from pol coefs to grid points
        self.coefsMatX = np.transpose(GetCoefsCheb(SizeX))
        self.coefsMatY = np.transpose(GetCoefsCheb(SizeY))
        self.coefsMat  = outer2D(self.coefsMatX, self.coefsMatY)

        print("Loading data grid")
        self.dataGrid = self.getDataGrid(data)
        self.data = data
        # tie(dataGrid, dataGridCov) = getDataGridWithCov();

    # function assumed to be normalized !!!
    def getLogLikelihoodSlow(self, pars):

        L = 0
        # pragma omp parallel for  reduction(+: L)
        for d in self.data:
            v = self.myFun(d, pars)
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

    def fitData(self, pars, useCheb=True):

        parsVec = []
        self.parsNames = []
        for k in pars:
            self.parsNames.append(k)
            parsVec.append(pars[k])

        self.useCheb = useCheb

        from iminuit import Minuit

        m = Minuit.from_array_func(self.evalVec, parsVec, name=self.parsNames)
        m.migrad()  # run optimiser

        return dict([[e, m.values[e]] for e in m.values])

    def funFast(self, x, pars):

        # Calculate only if the parameters changed
        if not hasattr(self, 'parsOld') or self.parsOld != pars:
            funVals = self.getLogFunction(pars)
            self.coefs = np.transpose(self.coefsMat) @ funVals
            self.parsOld = pars

        xMin = self.nodesX[0]
        xMax = self.nodesX[-1]
        yMin = self.nodesY[0]
        yMax = self.nodesY[-1]

        xRel = np.empty_like(x)
        xRel[0] = (x[0] - xMin) / (xMax - xMin)
        xRel[1] = (x[1] - yMin) / (yMax - yMin)

        polsX = getPols(len(self.nodesX), xRel[0])
        polsY = getPols(len(self.nodesY), xRel[1])
        pols  = np.outer(polsX, polsY).flatten()
        res = np.dot(self.coefs, pols)

        return np.exp(-0.5*res)

    # Get function projected to the x-axis
    def funFastProjX(self, x, pars):

        # Calculate only if the parameters changed
        if not hasattr(self, 'parsXOld') or self.parsXOld != pars:
            f = np.exp(-0.5 * self.getLogFunction(pars))
            funVals = f.reshape((self.nodesX.size, self.nodesY.size))

            funVals1D = (self.weightsY * funVals).sum(axis=1)

            self.coefsX = np.transpose(self.coefsMatX) @ funVals1D
            self.parsXOld = pars

        xMin = self.nodesX[0]
        xMax = self.nodesX[-1]

        xRel = (x - xMin) / (xMax - xMin)

        polsX = getPols(len(self.nodesX), xRel)
        res = np.dot(self.coefsX, polsX)
        return res

    # Get function projected to the y-axis
    def funFastProjY(self, y, pars):

        # Calculate only if the parameters changed
        if not hasattr(self, 'parsYOld') or self.parsYOld != pars:
            f = np.exp(-0.5 * self.getLogFunction(pars))
            funVals = f.reshape((self.nodesX.size, self.nodesY.size))

            funVals1D = (self.weightsX.reshape(-1, 1) * funVals).sum(axis=0)

            self.coefsY = np.transpose(self.coefsMatY) @ funVals1D
            self.parsYOld = pars

        yMin = self.nodesY[0]
        yMax = self.nodesY[-1]

        yRel = (y - yMin) / (yMax - yMin)

        polsY = getPols(len(self.nodesY), yRel)
        res = np.dot(self.coefsY, polsY)
        return res
