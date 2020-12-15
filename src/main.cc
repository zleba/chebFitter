#include <iostream>
#include <random>

#include "chebFitter.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"

struct RootFun{
double operator()(const double *xx )
{
   const Double_t x = xx[0];
   const Double_t y = xx[1];
   const Double_t tmp1 = y-x*x;
   const Double_t tmp2 = 1-x;
   return 100*tmp1*tmp1+tmp2*tmp2;
}

};

void getMinimum(chebFitter &cheb)
{
    ROOT::Math::Minimizer* minimum =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
 
    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.00001);
    minimum->SetPrintLevel(1);
 
    // create function wrapper for minimizer
    // a IMultiGenFunction type

    //RootFun rootFun;

    ROOT::Math::Functor f(cheb,2);
    double step[2] = {0.01,0.01};
    // starting point
 
    double variable[2] = { -1.,1.2};
    int randomSeed = -1;
    if (randomSeed >= 0) {
       TRandom2 r(randomSeed);
       variable[0] = r.Uniform(-20,20);
       variable[1] = r.Uniform(-20,20);
    }
 
    minimum->SetFunction(f);
 
    // Set the free variables to be minimized !
    minimum->SetVariable(0,"x",variable[0], step[0]);
    minimum->SetVariable(1,"y",variable[1], step[1]);
 
    // do the minimization
    minimum->Minimize();
}

using namespace std;

int main()
{
    chebFitter fitter;

    //generate syntetic data
     std::default_random_engine generator;
     std::normal_distribution<double> distribution(0.1,1.1);
    for(int i = 0; i < 1000000; ++i) {
        double r = distribution(generator);
        if(abs(r) < 1)
            fitter.data.push_back(r);
    }

    fitter.init(8*1024+1, -1, 1);
    cout << "Init done " << endl;
    fitter.getLogFunction({1, 3});

    //cout << "L slow " << fitter.getLogLikelihoodSlow({1, 3}) << endl;
    //cout << "L fast " << fitter.getLogLikelihoodFast({1, 3}) << endl;
    //return 0;


    //fit by slow exact likelyhood
   getMinimum(fitter);

    return 0;
}
