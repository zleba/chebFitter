#include "chebFitter2D.h"
#include "nodes.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>


#include "TGraph.h"
#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <Eigen/Core>
//#include <Eigen/kroneckerProduct>

using Eigen::MatrixXd;

using namespace std;


// tensor(outer) product of two vectors
VectorXd kroneckerProduct(VectorXd a, VectorXd b)
{
   VectorXd res(a.size()*b.size());

   //calc function values
   for(int i = 0; i < a.size(); ++i)
   for(int j = 0; j < b.size(); ++j) {
      int iGl = i*b.size() + j;
      res[iGl] = a[i]*b[j];
   }
   return res;
}

// tensor(outer) product of two matrices
MatrixXd kroneckerProduct(MatrixXd a, MatrixXd b)
{
   MatrixXd res(a.rows()*b.rows(), a.cols()*b.cols());

   //calc function values
   for(int i = 0; i < a.rows(); ++i)
   for(int j = 0; j < a.cols(); ++j)
   for(int k = 0; k < b.rows(); ++k)
   for(int l = 0; l < b.cols(); ++l) {
      int iGl = i*b.rows() + k;
      int jGl = j*b.cols() + l;
      res(iGl,jGl) = a(i,j) * b(k,l);
   }
   return res;
}





//return values of -2*log(p(x)), where p(x) is normalized to 1 over the fitted range
VectorXd chebFitter2D::getLogFunction(vector<double> pars) const
{
   VectorXd fVals(nodesX.size() * nodesY.size());

   //calc function values
   for(int i = 0; i < nodesX.size(); ++i)
   for(int j = 0; j < nodesY.size(); ++j) {
      int iGl = i*nodesY.size() + j;
      fVals[iGl] = myFun({nodesX[i], nodesY[j]}, pars);
   }

   //normalize the function
   double I = fVals.dot(weights);
   //cout << "Integral " << I << endl;
   fVals = -2 * log(fVals.array() / I); //normalize by integral

   return fVals; //Radek change

}

void chebFitter2D::init(int SizeX, double xMin, double xMax,  int SizeY, double yMin, double yMax)
{
   // loading the Cheb nodes
   nodesX = (xMax-xMin) * GetNodes(SizeX).array() + xMin;
   nodesY = (yMax-yMin) * GetNodes(SizeY).array() + yMin;

   // loding the weights for integration
   weightsX = (xMax - xMin) * GetWeights(SizeX);
   weightsY = (yMax - yMin) * GetWeights(SizeY);
   weights  = kroneckerProduct(weightsX, weightsY);


   // calculate the transformation matrix from pol coefs to grid points
   coefsMatX = GetCoefsCheb(SizeX).transpose();
   coefsMatY = GetCoefsCheb(SizeY).transpose();
   coefsMat  = kroneckerProduct(coefsMatX, coefsMatY);

   cout << "Loading data grid" << endl;
   dataGrid = getDataGrid();

}

//function assumed to be normalized !!!
double chebFitter2D::getLogLikelihoodSlow(vector<double> pars) const
{

   double L = 0;
   //#pragma omp parallel for  reduction(+: L)
   for(const auto &d : data) {
      double v = myFun(d, pars);
      //cout << v << endl;
      L += -2*log(v);
   }

   //cout << pars[0] <<" "<< pars[1] << " : " << L << endl;

   return L; 

}

//evaluation using cheb pols
double chebFitter2D::getLogLikelihoodFast(vector<double> pars) const
{
   VectorXd funVals = getLogFunction(pars);
   double LL = funVals.dot(dataGrid);

   for(auto p : pars) cout << p << " ";
   cout << " : "<<LL << endl;
   return LL;
}

double chebFitter2D::operator()(const double *par) const
{
   /*
   int N = 1e6;
   double a = -1, b = 1;
   double step = (b-a)/N;
   double sum = 0;
   for(int i = 0; i <=N; ++i) {
      double v = a + step*i;
      sum += myFun(v, {par[0], par[1]}) * step;
   }
   cout << "myRes " << sum << endl;
   */

   vector<double> pars(par, par + nPars);



   double LL = useCheb ? getLogLikelihoodFast(pars) : getLogLikelihoodSlow(pars);

   return LL;
   //return getLogLikelihoodSlow({par[0], par[1]});
}

// get data transformed into the grid such that (chebFunVals dot dataGrid) == logL
VectorXd chebFitter2D::getDataGrid() const
{
   double xMin = nodesX[0];
   double xMax = nodesX[nodesX.size()-1];
   double yMin = nodesY[0];
   double yMax = nodesY[nodesY.size()-1];


   VectorXd polSum = VectorXd::Zero(nodesX.size() * nodesY.size());
   for(auto xy : data) {
      double x = (xy[0] - xMin) / (xMax - xMin); //normalize between 0 and 1
      double y = (xy[1] - yMin) / (yMax - yMin); //normalize between 0 and 1
      VectorXd polX = getPols(nodesX.size(), x);
      VectorXd polY = getPols(nodesY.size(), y);

      polSum += kroneckerProduct(polX, polY);
   }
   cout << "Done " << endl;



   //transform to the basis of the cheb nodes
   //MatrixXd coefs = GetCoefsCheb(polSum.size());
   VectorXd gridVals = coefsMat * polSum;

   return gridVals;
}



//Minimize using ROOT minimizer
vector<double> chebFitter2D::fitData(vector<par> pars, bool UseCheb)
{
   nPars = pars.size();
   useCheb = UseCheb;

   ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

   // set tolerance , etc...
   minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
   minimum->SetMaxIterations(100000);  // for GSL
   minimum->SetTolerance(10.0);
   //minimum->SetPrecision(1e-5);

   minimum->SetPrintLevel(3);
   minimum->SetStrategy(2);
   minimum->SetErrorDef(1);

   // create function wrapper for minimizer
   // a IMultiGenFunction type
   ROOT::Math::Functor f(*this, pars.size());

   minimum->SetFunction(f);

   // Set the free variables to be minimized !
   for(int i = 0; i < pars.size(); ++i) {
      double step = (pars[i].vMax - pars[i].vMin) / 100;
      minimum->SetLimitedVariable (i, pars[i].name.Data(), pars[i].v, step, pars[i].vMin, pars[i].vMax);
      //minimum->SetVariable (i, pars[i].name.Data(), pars[i].v, step);
   }

   // do the minimization
   minimum->Minimize();
   //cout << minimum->X()[0] <<" "<< minimum->X()[1]  << endl;
   return vector<double>(minimum->X(), minimum->X() + pars.size());
}

/*
double chebFitter2D::getFunctionFast(vector<double> pars, double x)
{
    static VectorXd funVals = getLogFunction(pars);

    //cout << funVals << endl;

    return  exp(-0.5*interpol(nodes, funVals, x));

}
*/
