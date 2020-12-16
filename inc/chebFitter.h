#include <vector>
#include <utility>
#include <Eigen/Dense>
#include <functional>
#include "TString.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

struct par {
   TString name;
   double v, vMin, vMax;
};


struct chebFitter {

   VectorXd getLogFunction(std::vector<double> pars) const;
   double getLogLikelihoodSlow(std::vector<double> pars) const;
   double getLogLikelihoodFast(std::vector<double> pars) const;
   void init(int Size, double xMin, double xMax);
   VectorXd getDataGrid() const;
   std::pair<VectorXd,MatrixXd> getDataGridWithCov() const;

   std::vector<double> getMinimum(std::vector<par> pars);

   double operator()(const double *par) const;


   std::vector<double> data;     //vector with the data points to be fitted
   VectorXd dataGrid; //vector with the data points related to the cheb nodes (dataGrid.size = nodes.size)
   MatrixXd dataGridCov;

   VectorXd nodes;   // vector with cheb nodes
   VectorXd weights; // vector with cheb weights for integration

   std::function<double(double, std::vector<double>)> myFun;

};
