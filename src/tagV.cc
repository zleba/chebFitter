#include <iostream>
#include <random>
#include <vector>

#include "chebFitter2D.h"

 
using namespace std;

double Gaus(double x, double mu, double s)
{
    return 1./(sqrt(2*M_PI) * s) * exp(-1./2 * pow((x-mu)/s, 2));
}

//int a to b of exp(-c x^2 + 2*d -d*d/c)
double gausInt(double a, double b, double c, double d)
{
    double m = d / c;
    double s = 1./sqrt(c);

    double res = sqrt(M_PI)/(2*sqrt(c)) * /* exp(d*d/c) */ (TMath::Erf((b*c-d)/sqrt(c)) - TMath::Erf((a*c-d)/sqrt(c)));
    return res;
}


double convExpGaus(double sK, double tau, double a, double x)
{
    x -= a;


    double c = 1./(2*sK*sK);
    double d = 0.5*(x/sK/sK - 1/tau);

    double Const = 1./(sqrt(2*M_PI) *sK* tau) * exp(-1./2 * pow(x/sK,2) + d*d/c);

    //cout << "My " << gausInt(0, 1e15, c, d) << endl;
    //cout << "Const " << Const << endl;
    return Const * gausInt(0, 1e15, c, d);

}


double resFun(vector<double> x, vector<double> pars)
{
    double z     = x[0];
    double sigma = x[1];

    //parameters
    double fBig   = pars[0];
    double r      = pars[1];
    double fTR    = pars[2];
    double s0     = pars[3];
    double s1     = pars[4];
    double mu0    = pars[5];
    double mu1    = pars[6];
    double c1     = pars[7];
    double fTmax  = pars[8];
    double fTmu   = pars[9];
    double fTsigma= pars[10];



    double s  = s0 + s1 * sigma;
    double mu = mu0 + mu1 * sigma;
    double c = c1 / sigma;
    double fT = fTmax * TMath::Erf((sigma - fTmu)/ (fTsigma*sqrt(2)));


    double Delta = (1-fT) * ( (1 - fBig) * Gaus(z, mu, s) + fBig * Gaus(z, mu, r*s));

    double right = fBig * fTR     * convExpGaus(r*s, 1./c, 0., z);
    double left  = fBig * (1-fTR) * convExpGaus(r*s, 1./c, 0.,-z);

    return Delta + right + left;

}



void testTagV()
{

}

int main()
{
   testTagV();

   return 0;
}
