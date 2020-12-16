#include <vector>
#include <Eigen/Dense>
using Eigen::VectorXd;
using Eigen::MatrixXd;

//get weights to calculate the integral over the nodes
VectorXd GetWeights(int Size);

//get Cheb nodes
VectorXd GetNodes(int Size);


//Evaluate cheb pols to Size at point x
VectorXd getPols(int Size, double x);

//Evaluate sum of cheb pols to Size at vector x, x els are between 0 and 1
VectorXd getPolsSum(int Size, VectorXd x);


//Transformation matrix between cheb. nodes and cheb. coeficients
MatrixXd GetCoefs(int oldSize, bool isInverse = false);


//with better normalization of the borders
MatrixXd GetCoefsCheb(int oldSize);
