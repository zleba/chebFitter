#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>

#include <Eigen/Dense>
using  Eigen::VectorXd;
using  Eigen::MatrixXd;

using namespace std;

//For points between 0 and 1
VectorXd GetWeights(int Size)
{
    //const int Size = 33;
    const int N = Size - 1;
    assert(N%2 == 0);

    vector<vector<double>> coef(Size);
    for(auto & el : coef) el.resize(Size);


    for(int k = 0; k <= N/2; ++k) {
        double s = 0;
        coef[2*k][N] = 1./N;
        coef[2*k][0] = 1./N ;

        coef[2*k][N/2] = 2./N * (2* ((k+1)%2) -1);

        for(int n = 1; n <= N/2-1; ++n)
            coef[2*k][n] = coef[2*k][N-n] = 2./N * cos(n*k*M_PI*2/N);
    }

    //vector<double> wgt(Size, 0.);
    VectorXd wgt = VectorXd::Zero(Size);
    //arma::vec wgt(Size, arma::fill::zeros);


    for(int i = 0; i < Size; ++i) {
        wgt[i] += coef[0][i];
        wgt[i] += coef[N][i] / (1. - N*N);
        for(int k = 1; k <= N/2 - 1; ++k) {
           double w = 2./(1 - 4*k*k);
           wgt[i] += w*coef[2*k][i];
        }

        wgt[i] *= 0.5; //for interval (0,1)
    }
    return wgt;
    //for(int i = 0; i < Size; ++i)
        //cout << i <<" "<<setprecision(17)<< wgt[i] << endl;
}

//Get vector with nodes by definition between 0 and 1
VectorXd GetNodes(int Size)
{
    assert((Size - 1) % 2 == 0);
    //arma::vec xi(Size, arma::fill::zeros);
    VectorXd xi = VectorXd::Zero(Size);
    for(int i = 0; i < Size; ++i) {
      double Cos = cos(i /(Size-1.) * M_PI);
      xi[i] = (1-Cos)/2;
    }
    return xi;
}

//Evaluate cheb pols to Size at point x, x is between 0 and 1
VectorXd getPols(int Size, double x)
{
    VectorXd pol(Size);

    if(Size >= 1) pol[0] = 1;
    if(Size >= 2) pol[1] = 2*x-1;

    for(int i = 2; i < Size; ++i)
        pol[i] = 2*(2*x-1)*pol[i-1] - pol[i-2];
    return pol;
}

//Evaluate sum of cheb pols to Size at vector x, x els are between 0 and 1
VectorXd getPolsSum(int Size, VectorXd x)
{
   assert(Size > 2);

   VectorXd polSum(Size);

   VectorXd pol0 = 0*x.array()+ 1;
   VectorXd pol1 = 2*x.array()-1;
   VectorXd C    = 2 * pol1;

   VectorXd pol2(x.size());
    for(int i = 2; i < Size; ++i) {
       polSum(i-2) = pol0.sum();

       pol2 = C.array() * pol1.array() - pol0.array();

       pol0 = pol1;
       pol1 = pol2;
    }

    polSum(Size-2) = pol0.sum();
    polSum(Size-1) = pol1.sum();

    return polSum;
}





//Transformation matrix between cheb. nodes and cheb. coeficients
MatrixXd GetCoefs(int oldSize, bool isInverse = false)
{
    const int N = oldSize - 1;
    assert(N%2 == 0);

    MatrixXd  coef(oldSize, oldSize);

    //vector<vector<double>> coef(oldSize);
    //for(auto & el : coef) el.resize(oldSize);
    //arma::mat coef(oldSize,oldSize);

    double mul = 1;
    double C = 1./N;
    if(isInverse == true) {C = 1./2; }

    //isInverse = false;
    for(int k = 0; k <= N; ++k) {
        double s = 0;
        if(!isInverse) {
            coef(k,N) = C;
            coef(k,0) = C * (k % 2 == 1 ? -1 : 1);
        }
        else {
            mul = k % 2 == 1 ? -1 : 1;
            coef(N-k, N) = C * mul;
            coef(N-k, 0) = C ;
        }

        for(int n = 1; n <= N-1; ++n) {
            double el = cos(n*k*M_PI / N) * 2.*C * mul;
            if(!isInverse) coef(k,N-n) = el;
            else           coef(N-k,N-n) = el;
        }
    }
    
    return coef;
}




//with better normalization of the borders
MatrixXd GetCoefsCheb(int oldSize)
{
    auto coef = GetCoefs(oldSize);

    coef.row(0) *= 0.5;
    coef.row(coef.rows()-1) *= 0.5;

    return coef;
}




//Evaluate Cheb. pol at point x
double evalPol(const VectorXd &polCoef, double x)
{
    VectorXd pols = getPols(polCoef.size(), x);

    //double s = 0;
    //for(int i = 0; i < pols.size(); ++i)
       //s += pols[i] * polCoef[i];

    double s = pols.dot(polCoef);

    return s;
    //return arma::dot(pols,polCoef);
}

/*

//Get Interpolation vector at point x
arma::vec interpol(arma::vec xi, double x)
{
    arma::vec coefs(xi.n_rows);
    for(int i = 0; i < xi.n_rows; ++i) {
        double num = 1, den = 1;
        for(int j = 0; j < xi.n_rows; ++j)
            if(j != i) {
                num *= x     - xi(j);
                den *= xi(i) - xi(j);
            }
        coefs(i) = num/den;
    }
    return coefs;
}

//Get interpolated function value at point x
double interpol(arma::vec xi, arma::vec vals, double x)
{
    arma::vec coefs = interpol(xi, x);
    return arma::dot(coefs,vals);
}



//Transformation from chebNodes between 0 and 1 to chebNodes between a and b
arma::mat GetCoefs(int Size, double a, double b)
{
    arma::vec xi   = a + (b-a)*GetNodes(Size);

    arma::mat polsAll(Size,Size);
    for(int i = 0; i < Size; ++i) {
        polsAll.row(i) = getPols(Size,xi(i)).t();
    }

    //arma::mat coef = GetCoefs(Size); //now we have cheb pol coef
    arma::mat coef = GetCoefsCheb(Size); //now we have cheb pol coef

    return polsAll*coef;
}



//Transformation from chebNodes between 0 and 1 to chebNodes between a and b
arma::mat GetCoefsGeneric(int Size, int SizeI, double a, double b)
{
    //SizeI : the integration size
    arma::vec xi   = a + (b-a)*GetNodes(SizeI);

    arma::mat polsAll(SizeI,Size);
    for(int i = 0; i < SizeI; ++i) {
        polsAll.row(i) = getPols(Size,xi(i)).t();
    }

    //arma::mat coef = GetCoefs(Size); //now we have cheb pol coef
    arma::mat coef = GetCoefsCheb(Size); //now we have cheb pol coef

    return polsAll*coef;
}

*/
