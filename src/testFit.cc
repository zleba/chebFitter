#include <iostream>
#include <random>
#include <vector>

#include "chebFitter.h"
#include "chebFitter2D.h"

#include <iostream>
#include <Eigen/Dense>
 
using namespace std;

using Eigen::MatrixXd;


using namespace std;

//note that normalization to 1 is necessary only for running using "slow" minimization
struct myGaus1D{
   double xMin = -1, xMax = 1;
   double operator()(double x, vector<double> pars)
   {
       double a = xMin - pars[0];
       double b = xMax - pars[0];
       double N = sqrt(M_PI/2) * pars[1]* ( erf(b/sqrt(2)/pars[1])  - erf(a/sqrt(2)/pars[1]) );
       double f = 1./N * exp( -1./2 * pow((x-pars[0])/pars[1], 2));
       assert(f >= 0);
      return f;
   }
};


//assuming no correlations between x and y
struct myGaus2D{

   myGaus1D gX, gY;

   double operator()(vector<double> xy, vector<double> pars)
   {
      double x = xy[0];
      double y = xy[1];
      double f = gX(x, {pars[0], pars[1]}) * gY(y, {pars[2], pars[3]});
      assert(f >= 0);
      return f;
   }
};

//test 1D minimization
void test1D()
{
    chebFitter fitter;
    myGaus1D myGaus;
    fitter.myFun = myGaus;

    //generate synthetic data
     std::default_random_engine generator;
     std::normal_distribution<double> distribution(0.1,1.1);
    for(int i = 0; i < 50000000; ++i) {
        double r = distribution(generator);
        if(abs(r) < 1)
            fitter.data.push_back(r);
    }

    //Init the fitter: number of cheb nodes (more= higher precision), xMin, xMax
    fitter.init(64+1, -1, 1);
    cout << "Init done " << endl;


    //fit by slow exact likelihood
    fitter.fitData({
          {"mean",  0.0,  -2, 2},
          {"sigma", 1,     0, 2} },
          true);

}

void test2D()
{
    chebFitter2D fitter;
    myGaus2D myGaus;
    fitter.myFun = myGaus; //which function to minimize

    //generate synthetic data and put them into the fitter
     std::default_random_engine generator;
     std::normal_distribution<double> distributionX(0.1,1.1);
     std::normal_distribution<double> distributionY(0.3,1.6);
    for(int i = 0; i < 5000000; ++i) {
        double x = distributionX(generator);
        double y = distributionY(generator);
        if(abs(x) < 1 && abs(y) < 1)
            fitter.data.push_back({x,y});
    }

    //Init the fitter: number of cheb nodes in x (more= higher precision), xMin, xMax, the say for other coordinate
    fitter.init(64+1, -1,1, 64+1, -1,1);
    cout << "Init done " << endl;


    //fit by slow exact likelihood
    fitter.fitData({
          {"meanX",  0.0,  -2, 2},
          {"sigmaX", 1,     0, 2},
          {"meanY",  0.0,  -2, 2},
          {"sigmaY", 1,     0, 2}},
          true);




}

int main()
{
   test2D();

   return 0;
}
