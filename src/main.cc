#include <iostream>
#include <random>
#include <vector>

#include "chebFitter.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"


#include <iostream>
#include <Eigen/Dense>
 
using namespace std;

using Eigen::MatrixXd;
 



using namespace std;

double myGaus(double x, vector<double> pars)
{
    double a = -1 - pars[0];
    double b = 1 - pars[0];
    double N = sqrt(M_PI/2) * pars[1]* ( erf(b/sqrt(2)/pars[1])  - erf(a/sqrt(2)/pars[1]) );
    double f = 1./N * exp( -1./2 * pow((x-pars[0])/pars[1], 2));
    assert(f >= 0);
   return f;
}



int main()
{
    chebFitter fitter;
    fitter.myFun = myGaus;

    //generate syntetic data
     std::default_random_engine generator;
     std::normal_distribution<double> distribution(0.1,1.1);
    for(int i = 0; i < 1000000; ++i) {
        double r = distribution(generator);
        if(abs(r) < 1)
            fitter.data.push_back(r);
    }

    fitter.init(256+1, -1, 1);
    cout << "Init done " << endl;


    //fit by slow exact likelyhood
   fitter.getMinimum({
                     {"mean",  0.1,  2, 2},
                     {"sigma", 1,    0, 2} });

    return 0;
}
