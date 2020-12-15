#include "chebFitter.h"
#include "nodes.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>



using namespace std;

double myFun(double x, vector<double> pars)
{
    double a = -1 - pars[0];
    double b = 1 - pars[0];
    double N = sqrt(M_PI/2) * pars[1]* ( erf(b/sqrt(2)/pars[1])  - erf(a/sqrt(2)/pars[1]) );
    double f = 1./N * exp( -1./2 * pow((x-pars[0])/pars[1], 2));
    assert(f >= 0);
   return f;
}



//return values of -log(p(x)), where p(x) is normalized to 1 over the fitted range
vector<double> chebFitter::getLogFunction(vector<double> pars)
{
   vector<double> fVals(nodes.size());

   //calc function values
   for(int i = 0; i < nodes.size(); ++i)
      fVals[i] = myFun(nodes[i], pars);

   //normalize the function
   double I = 0;
   for(int i = 0; i < nodes.size(); ++i)
      I += fVals[i] * weights[i];

   //cout << "Integral " << I << endl;
   for(auto &el : fVals)
      el = -2 * log(el/I);

   return fVals;

}

void chebFitter::init(int Size, double xMin, double xMax)
{
   // loading the Cheb nodes
   nodes = GetNodes(Size);
   for(auto & el : nodes)
      el = xMin + (xMax - xMin) * el;

   // loding the weights for integration
   weights = GetWeights(Size);
   for(auto & el : weights)
      el = (xMax - xMin) * el;

   cout << "Loading data grid" << endl;
   dataGrid = getDataGrid();

}

//function assumed to be normalized !!!
double chebFitter::getLogLikelihoodSlow(vector<double> pars)
{

   double L = 0;
   //#pragma omp parallel for  reduction(+: L)
   for(double d : data) {
      double v = myFun(d, pars);
      //cout << v << endl;
      L += -2*log(v);
   }


   //cout << pars[0] <<" "<< pars[1] << " : " << L << endl;

   //exit(0);
   return L; 

}

//function assumed to be normalized !!!
double chebFitter::getLogLikelihoodFast(vector<double> pars)
{
   vector<double> funVals = getLogFunction(pars);

   double L = 0;
   for(int i = 0; i < funVals.size(); ++i) {
      L += funVals[i] * dataGrid[i];
   }
   return L; 

}

double chebFitter::operator()(const double *par)
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
vector<double> chebFitter::getDataGrid()
{
   double a = nodes.front();
   double b = nodes.back();


   vector<double> polSum(nodes.size(), 0.);
   for(double x : data) {
      double xx = (x - a) / (b - a);
      vector<double> pol =  getPols(nodes.size(), xx);

      for(int i = 0; i < pol.size(); ++i)
         polSum[i] += pol[i];

   }
   cout << "Done " << endl;

   //transform to the basis of the cheb nodes
   vector<vector<double>> coefs = GetCoefsCheb(polSum.size());

   vector<double> gridVals(polSum.size(), 0.0);
   for(int i = 0; i < gridVals.size(); ++i) {
      for(int j = 0; j < polSum.size(); ++j) {
         gridVals[i] += coefs[j][i] * polSum[j];

      }
   }

   return gridVals;
}
