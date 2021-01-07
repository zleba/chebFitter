#include <iostream>
#include <random>
#include <vector>

#include "chebFitter2D.h"

#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH1D.h"

#include <fstream>
 
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

double gaussianInt(double x, double mu, double s)
{
    double xMin = -100000;
    return 1./2 * ( TMath::Erf((x-mu)/(sqrt(2)*s))  -  TMath::Erf((xMin-mu)/(sqrt(2)*s)) );



}

double resFun(vector<double> x, vector<double> pars)
{
    double z     = x[0];
    double sigma = x[1];

    //parameters

    //double fBig   = pars[0];
    //double r      = pars[1];
    //double fTR    = pars[2];
    //double s0     = pars[3];
    //double s1     = pars[4];
    //double mu0    = pars[5];
    //double mu1    = pars[6];
    //double c1     = pars[7];
    //double fTmax  = pars[8];
    //double fTmu   = pars[9];
    //double fTsigma= pars[10];

    double mu0          =  pars[0];
    double mu1          =  pars[1];
    double mu2          =  0; //pars[2];
    double sigma0       =  pars[2];
    double sigma1       =  pars[3];
    double sigma2       =  0;//pars[5];
    double c0L          =  0;//pars[6];
    double c0R          =  0;//pars[7];
    double cL           =  pars[4];
    double cR           =  pars[4];
    double fTR          =  pars[5];
    double fTslope      =  pars[6];
    double fTmax        =  pars[7];
    double fTmu         =  pars[8];
    double fTsigma      =  pars[9];
    double bigSigmaFrac =  pars[10];
    double bigSigmaScale=  pars[11];



    //double sigmas  = s0 + s1 * sigma;
    //double mu = mu0 + mu1 * sigma;
    //double c = c1 / sigma;
    //double fT = fTmax * TMath::Erf((sigma - fTmu)/ (fTsigma*sqrt(2)));


    //double Delta = (1-fT) * ( (1 - fBig) * Gaus(z, mu, sigmas) + fBig * Gaus(z, mu, r*sigmas));
    //double right = fBig * fTR     * convExpGaus(r*sigmas, 1./c, 0., z);
    //double left  = fBig * (1-fTR) * convExpGaus(r*sigmas, 1./c, 0.,-z);

    double yrounded = sigma;

    double mus = mu0 + mu1*yrounded + mu2*yrounded*yrounded;

    double sigmas = sigma0 + sigma1*yrounded + sigma2*yrounded*yrounded;
    //sigmas = np.where(sigmas < 1e-8, 1e-8, sigmas)

    double cLs = 1./(c0L + cL*yrounded);
    double cRs = 1./(c0R + cR*yrounded);

    double fTs = (1.-fTslope*yrounded)*fTmax * gaussianInt(yrounded, fTmu, fTsigma);


    
    double ret = 0;
    ret += (1-fTs)*(1-bigSigmaFrac)*Gaus(z, mus, sigmas);
    ret += (1-fTs)*bigSigmaFrac    *Gaus(z, mus, bigSigmaScale*sigmas);

    ret += fTs*    fTR *(1-bigSigmaFrac)* convExpGaus(sigmas, 1./cRs, 0., (z - mus));
    ret += fTs*(1.-fTR)*(1-bigSigmaFrac)* convExpGaus(sigmas, 1./cLs, 0.,-(z - mus));

    //ret += fTs*(1.-fTR)*bigSigmaFrac* CB( x, mu= mus, sigma=bigSigmaScale*sigmas, alpha = aLOL, n=nLOL);
    //ret += fTs*     fTR*bigSigmaFrac* CB(-x, mu=-mus, sigma=bigSigmaScale*sigmas, alpha = aROL, n=nROL);

    /*
        cout << "Helenka " << ret << endl;

        cout << "Radek " << 
mu0       <<" "<<
mu1       <<" "<<
sigma0    <<" "<<
sigma1     <<" "<<
c0L        <<" "<<
c0R        <<" "<<
cL         <<" "<<
cR         <<" "<<
fTR        <<" "<<
fTslope    <<" "<<
fTmax      <<" "<<
fTmu        <<" "<<
fTsigma     <<" "<<
bigSigmaFrac<<" "<<
bigSigmaScale<<" "<< endl;

    if(ret != ret || ret <= 0) {
        exit(0);
    }

    */

    //exit(0);

    if(ret != ret) return 1e-14;
    ret = max(1e-14, ret);

    return ret;

}


vector<TString> getTokens(TString str)
{
    vector<TString> tokens;
    auto tempVec = str.Tokenize(",");
    for (int i = 0; i < tempVec->GetEntries(); ++i) {
        TString s = ((TObjString*)tempVec->At(i))->GetString();
        tokens.push_back(s.Strip());
    }
    return tokens;
}




vector<vector<double>> readData(TString fName, TString var1, TString var2)
{
    vector<vector<double>> data;

    ifstream  file(fName);

    int i = 0;
    int j1 = -1, j2 = -1;
    int jSize = -1;

    while (true) {

        string str;
        std::getline(file, str);
        if(!file.good()) break;

        auto tokens = getTokens(str);
        if(i == 0) {
            jSize = tokens.size();
            for(int j = 0; j < tokens.size(); ++j) {
                if(tokens[j] == var1) j1 = j;
                if(tokens[j] == var2) j2 = j;
            }
            assert(j1 >= 0 && j2 >= 0);
        }
        else { //not first line
            
            if(jSize != tokens.size()) continue;

            //cout << "test " << i << " "<< jSize <<" "<< tokens.size() << endl;
            assert(j1 < tokens.size());
            assert(j2 < tokens.size());

            double v1 = stod(tokens[j1].Data());
            double v2 = stod(tokens[j2].Data());
            data.push_back({v1, v2});
        }
        

        //for(auto e : tokens)
        //    cout << "H " << e << endl;
        //cout << endl;
        //break;
        ++i;
    }
    return data;

}

void testTagV()
{
    //loading data
    auto data = readData("/home/radek/belle/notebooks/thumair_jpsiks/notebooks/dataThibaud/datafortres.csv", "tagvlres", "tagvlerr");

    chebFitter2D fitter;
    fitter.myFun = resFun; //which function to minimize

    double zMin = -1000;
    double zMax = +1000;
    double errMin = 0.01;
    double errMax = 1000;

    for(auto d : data) {
        if(zMin < d[0]    && d[0] < zMax)
        if(errMin < d[1]  && d[1] < errMax)
            fitter.data.push_back(d);


    }

    //TH1D *h = new TH1D("h", "", 100, -1000, 1000);
    //for(auto x : data) {
    //    h->Fill(x[1]);
    //}
    //h->SaveAs("test.root");
    //return;

    cout << data.size() << endl;

    //fitter.data = data;

    fitter.init(128+1, zMin,zMax, 128+1, errMin,errMax);

    fitter.fitData({
        {"mu0",             3.6032098496724387,     1.0, 50.0},
        {"mu1",             -0.05308279842953489, -0.5, 0.5},
        {"sigma0",          2.7507848066603344,    0.01, 10},
        {"sigma1",          0.9817324995811869,    0.01, 10},
        {"c",               0.3582878946759941,    0.04, 300},
        {"fTR",             0.866923730831205,     0.02, 1},
        {"fTslope",         0,                       0, 0},
        {"fTmax",           0.27855665545957287,    0.03, 3},
        {"fTmu",            17.108847316080983,    0.1, 400},
        {"fTsigma",         4.464958554742986,     0.1, 100},
        {"bigSigmaFrac",    0.2513536792932563,    0.02, 1},
        {"bigSigmaScale",   1.5884729844719558,    0.1, 50}
    }, true);




}

int main()
{
   testTagV();

   return 0;
}
