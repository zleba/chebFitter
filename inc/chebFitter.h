#pragma once

#include <vector>
#include <utility>
#include <Eigen/Dense>
#include <functional>
#include "TString.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

//to describe each parameter to minimize
struct par {
   TString name; //name of the parameter
   double v, vMin, vMax; //starting value & limits
};


struct chebFitter {

   VectorXd getLogFunction(std::vector<double> pars) const;
   double getLogLikelihoodSlow(std::vector<double> pars) const;
   double getLogLikelihoodFast(std::vector<double> pars) const;
   void init(int Size, double xMin, double xMax);
   VectorXd getDataGrid() const;
   std::pair<VectorXd,MatrixXd> getDataGridWithCov() const;

   double getFunctionFast(std::vector<double> pars, double x);

   std::vector<double> fitData(std::vector<par> pars, bool UseCheb = true);

   double operator()(const double *par) const;


   std::vector<double> data;     //vector with the data points to be fitted
   VectorXd dataGrid; //vector with the data points related to the cheb nodes (dataGrid.size = nodes.size)
   MatrixXd dataGridCov;
   MatrixXd coefsMat; // transofrmation matrix from chebPol to gridPoints

   VectorXd nodes;   // vector with cheb nodes
   VectorXd weights; // vector with cheb weights for integration
   int nPars;
   bool useCheb;

   std::function<double(double, std::vector<double>)> myFun;

};
