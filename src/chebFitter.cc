#include "chebFitter.h"
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

using Eigen::MatrixXd;

using namespace std;



//return values of -2*log(p(x)), where p(x) is normalized to 1 over the fitted range
VectorXd chebFitter::getLogFunction(vector<double> pars) const
{
   VectorXd fVals(nodes.size());

   //calc function values
   fVals = nodes.unaryExpr([&](double x) { return myFun(x, pars); });

   //normalize the function
   double I = fVals.dot(weights);
   //cout << "Integral " << I << endl;
   fVals = -2 * log(fVals.array() / I); //normalize by integral

   return fVals;

}

void chebFitter::init(int Size, double xMin, double xMax)
{
   // loading the Cheb nodes
   nodes = (xMax-xMin) * GetNodes(Size).array() + xMin;

   // loding the weights for integration
   weights = (xMax - xMin) * GetWeights(Size);

   cout << "Loading data grid" << endl;
   dataGrid = getDataGrid();
   //tie(dataGrid, dataGridCov) = getDataGridWithCov();

   /*
   TGraph *gr = new TGraph;
   for(int i = 0; i < dataGrid.size(); ++i) {
      gr->SetPoint(i, nodes[i], dataGrid[i]);
   }
   gr->SetName("gr");
   gr->SaveAs("test.root");
   */
}

//function assumed to be normalized !!!
double chebFitter::getLogLikelihoodSlow(vector<double> pars) const
{

   double L = 0;
   //#pragma omp parallel for  reduction(+: L)
   for(double d : data) {
      double v = myFun(d, pars);
      //cout << v << endl;
      L += -2*log(v);
   }

   cout << pars[0] <<" "<< pars[1] << " : " << L << endl;

   return L; 

}

//evaluation using cheb pols
double chebFitter::getLogLikelihoodFast(vector<double> pars) const
{
   VectorXd funVals = getLogFunction(pars);
   return funVals.dot(dataGrid);
}

double chebFitter::operator()(const double *par) const
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



   return getLogLikelihoodFast({par[0], par[1]});
   //return getLogLikelihoodSlow({par[0], par[1]});
}

// get data transformed into the grid such that (chebFunVals dot dataGrid) == logL
VectorXd chebFitter::getDataGrid() const
{
   double a = nodes[0];
   double b = nodes[nodes.size()-1];


   VectorXd polSum = VectorXd::Zero(nodes.size());
   for(double x : data) {
      double xx = (x - a) / (b - a); //normalize between 0 and 1
      polSum += getPols(nodes.size(), xx);
   }
   cout << "Done " << endl;



   /*
   //new method
   VectorXd dataNew(data.size());
   for(int i = 0; i < data.size(); ++i)
      dataNew[i] = (data[i] - a) / (b - a);

   cout << "Starting " << endl;
   VectorXd polSum = getPolsSum(nodes.size(), dataNew);
   */


   //transform to the basis of the cheb nodes
   MatrixXd coefs = GetCoefsCheb(polSum.size());
   VectorXd gridVals = coefs.transpose() * polSum;

   return gridVals;
}


// get data transformed into the grid such that (chebFunVals dot dataGrid) == logL
pair<VectorXd,MatrixXd> chebFitter::getDataGridWithCov() const
{
   double a = nodes[0];
   double b = nodes[nodes.size()-1];

   MatrixXd polSum2 = MatrixXd::Zero(nodes.size(), nodes.size());
   VectorXd polSum  = VectorXd::Zero(nodes.size());
   for(double x : data) {
      double xx = (x - a) / (b - a); //normalize between 0 and 1
      VectorXd pol = getPols(nodes.size(), xx);
      polSum  += pol;
      polSum2 += pol * pol.transpose();
   }
   cout << "Done " << endl;


   //transform to the basis of the cheb nodes
   MatrixXd coefs = GetCoefsCheb(polSum.size()).transpose();
   VectorXd gridVals = coefs * polSum;
   MatrixXd gridValsCov = coefs * polSum2 * coefs.transpose();

   return make_pair(gridVals, gridValsCov);
}


//Minimize using ROOT minimizer
vector<double> chebFitter::getMinimum(vector<par> pars)
{
    ROOT::Math::Minimizer* minimum =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
 
    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(100000);  // for GSL
    minimum->SetTolerance(0.01);
    minimum->SetPrintLevel(1);
 
    // create function wrapper for minimizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(*this, pars.size());
 
    minimum->SetFunction(f);
 
    // Set the free variables to be minimized !
    for(int i = 0; i < pars.size(); ++i) {
       double step = (pars[i].vMax - pars[i].vMin) / 10;
       minimum->SetLimitedVariable (i, pars[i].name.Data(), pars[i].v, step, pars[i].vMin, pars[i].vMax);
    }

    // do the minimization
    minimum->Minimize();
    cout << minimum->X()[0] <<" "<< minimum->X()[1]  << endl;
    return vector<double>(minimum->X(), minimum->X() + pars.size());
}
