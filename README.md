# chebFitter
A repo to perform 1D or 2D maximum likelihood fit of a complicated function by caching it to the Chebyshev polynomial.

## Prerequisites
This program requires ROOT for Minuit minimizer and Eigen for the linear algebra

To install ROOT the binary can be downloaded from the following webpage
https://root.cern/install/all_releases/

The Eigen library can be found in the Ubuntu repository 
```
sudo apt install libeigen3-dev
```

Both ROOT and Eigen are also available withing BASF2 environment.

## Installation
Proceed in the standard way in the directory where you want to install the program:
```
git clone https://github.com/zleba/chebFitter.git
cd chebFitter
make -j4
```

The command `make` complies the program, `-j4` means running the compilation at 4 cores (can be omitted).

## Running C++ version

Please look at file `src/testFit.cc` where the chebFitter is used to fit 1D or 2D Gaussian distribution.
There is a possibility to use the classical way of evaluating the likelihood (= "slow" method) or the "fast" way where the Chebyshev polynomials are used.
The classical way is kept as a crosscheck, ideally both approaches should give identical result.

In the example the fitted function is interpolated using 64+1 Chebyshev nodes (=Cheb pol of order 64+1).
There is no difference in the results compared to the exact "slow" fit in the first 6 significant digits.


The optimal number of the Chebyshev nodes depends on the shape of the function, strongly oscillating functions or functions with sharp peak can require more points.


The advantage of the "fast" fitting is that the function does not need to be normalized to 1, i.e. the function is internally normalized to 1 automatically.

Run the example as:
```
./testFit
```
## Running python version

Run the follwing commands
```
cd python
jupyter-lab
```
This will open the JupyterLab in the browser, open the `testFitter.ipynb` notebook and run it.
Ensure that `numpy`, `iminuit`, `matplotlib` and `numba` python packages are installed.
