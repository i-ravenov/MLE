#include <matrix.h>
#include <mex.h>   

#include "nr3/nr3.h"
#include "nr3/gamma.h"
#include "nr3/stepper.h"
#include "nr3/stepperbs.h"
#include "nr3/odeint.h"
#include "nr3/hypgeo.h"

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

using namespace std;

const double eps = 0.00000001;


inline double gamma(const double xx)
{
	return exp(gammln(xx));
}

struct I1
{
	inline double operator()(double a, double b, double alpha, double beta)
	{
		if ((abs(a) < eps) && abs(b - 1) < eps) {
			return gamma(1-alpha) * gamma(1-beta) / gamma(2-alpha-beta);
		}
		else if (abs(b - 1) < eps) {
			return 1.0 / (1 - alpha) * ( pow(a, -beta) * hypgeo(1, beta, 2 - alpha, -1 / a).real() );
		}
		else {
			return 	1.0 / ((a + b)*(beta - 1)) * (  pow(a, 1 - beta)*pow(b, 1 - alpha)*hypgeo(1, 2 - alpha - beta, 2 - beta, a / (a + b)).real()
				- pow(1 + a, 1 - beta)* pow(b - 1, 1 - alpha) * hypgeo(1, 2 - alpha - beta, 2 - beta, (1 + a) / (a + b)).real()   );
		}
	}
};

struct I2 : I1
{
	inline double operator()(double a, double b, double alpha, double beta)
	{
		return I1::operator()(a, b, alpha, beta - 1) - a * I1::operator()(a, b, alpha, beta);
	}
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *Res_out_m;
    double *Res;
    double a, b, alpha, beta;

//associate inputs
    a = mxGetScalar(prhs[0]);
    b = mxGetScalar(prhs[1]);
    alpha = mxGetScalar(prhs[2]);
    beta = mxGetScalar(prhs[3]);

//associate outputs
    Res_out_m = plhs[0] = mxCreateDoubleMatrix(1,2,mxREAL);
    
//associate pointers
    Res = mxGetPr(Res_out_m);
    
//calculate integrals
    I1 i1;
    I2 i2;
    Res[0] = i1(a, b, alpha, beta);
    Res[1] = i2(a, b, alpha, beta);
    return;
}