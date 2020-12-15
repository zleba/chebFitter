#include <vector>


struct chebFitter {

   std::vector<double> getLogFunction(std::vector<double> pars);
   double getLogLikelihoodSlow(std::vector<double> pars);
   double getLogLikelihoodFast(std::vector<double> pars);
   void init(int Size, double xMin, double xMax);
   std::vector<double> getDataGrid();


   double operator()(const double *par);



   std::vector<double> data;     //vector with the data points to be fitted
   std::vector<double> dataGrid; //vector with the data points related to the cheb nodes (dataGrid.size = nodes.size)

   std::vector<double> nodes;   // vector with cheb nodes
   std::vector<double> weights; // vector with cheb weights for integration




};
