#pragma once

#include <vector>
#include <utility>
#include <Eigen/Dense>
#include <functional>
#include "TString.h"

#include "chebFitter.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

/*
//to describe each parameter to minimize
struct par {
   TString name; //name of the parameter
   double v, vMin, vMax; //starting value & limits
};
*/


struct chebFitter2D {

   VectorXd getLogFunction(std::vector<double> pars) const;
   double getLogLikelihoodSlow(std::vector<double> pars) const;
   double getLogLikelihoodFast(std::vector<double> pars) const;
   void init(int SizeX, double xMin, double xMax,  int SizeY, double yMin, double yMax);
   VectorXd getDataGrid() const;

   //double getFunctionFast(std::vector<double> pars, double x);

   std::vector<double> fitData(std::vector<par> pars, bool UseCheb = true);

   double operator()(const double *par) const;


   std::vector<std::vector<double>> data;     //vector with the data points to be fitted
   VectorXd dataGrid; //vector with the data points related to the cheb nodes (dataGrid.size = nodes.size)
   MatrixXd dataGridCov;
   MatrixXd coefsMatX, coefsMatY, coefsMat; // transofrmation matrix from chebPol to gridPoints

   VectorXd nodesX, nodesY;   // vector with cheb nodes
   VectorXd weightsX, weightsY, weights; // vector with cheb weights for integration
   int nPars;
   bool useCheb; //use the Cheb pols to get maximum likelihood

   std::function<double(std::vector<double>, std::vector<double>)> myFun;

};
