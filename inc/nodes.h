#include<vector>

//get weights to calculate the integral over the nodes
std::vector<double> GetWeights(int Size);

//get Cheb nodes
std::vector<double> GetNodes(int Size);


//Evaluate cheb pols to Size at point x
std::vector<double> getPols(int Size, double x);


//Transformation matrix between cheb. nodes and cheb. coeficients
std::vector<std::vector<double>> GetCoefs(int oldSize, bool isInverse = false);


//with better normalization of the borders
std::vector<std::vector<double>> GetCoefsCheb(int oldSize);
