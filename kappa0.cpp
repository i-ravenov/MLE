#include <matrix.h>
#include <mex.h>   

#include "nr3/nr3.h"
#include "nr3/gamma.h"
#include "nr3/stepper.h"
#include "nr3/stepperbs.h"
#include "nr3/odeint.h"
#include "nr3/hypgeo.h"
#include "nr3/qgaus.h"

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

const double H1 = 0.6;
const double H2 = 0.8;
const double alpha1 = 1 / (H2 - H1);
const double alpha2 = 1 / (2 - 2*H2);
const double alpha = H2*(2 * H2 - 1);
const double betaH2 = sqrt(alpha / (beta(H2 - 0.5, 2 - 2 * H2)));
const double cnst = pow(betaH2*(H2 - H1), 2);
const double eps = 0.00000001;
const double edge = 0.5;

inline double gamma(const double xx)
{
	return exp(gammln(xx));
}

struct H2F1_1 
{
	H2F1_1() {};
	inline double operator() (double t, double s) 
	{
		if (abs((t - s) / t - 1) < eps) return gamma(1-H1+H2)*gamma(2*H2-H1 - 0.5) / ( gamma(2*H2-2*H1+1) * gamma(H2-0.5) );
		return hypgeo(H1 - H2, 3. / 2. - H1, 1 - H1 + H2, (t - s) / t).real();  // 4th argument can not be equal 1
	};
};

struct H2F1_2 
{
	H2F1_2() {};
	inline double operator() (double t, double s) 
	{
		if (abs((t - s) / t - 1) < eps) return gamma(2 - H1 + H2)*gamma(2 * H2 - H1 - 0.5) / (gamma(2 * H2 - 2 * H1 + 1) * gamma(0.5 + H2));
		return hypgeo(H1 - H2 + 1, 3. / 2. - H1, 2. - H1 + H2, (t - s) / t).real();  // 4th argument can not be equal 1
	}
};

struct Phi : H2F1_1, H2F1_2 
{
	Phi() {};
	inline double operator() (double t, double s) 
	{
		return pow(t, H2 - H1) * beta(3. / 2. - H1, H2 - 1. / 2.) * H2F1_1::operator()(t, s) +
			pow(t - s, H2 - H1) * pow( (t - s) / t, 1 - H2 + H1 ) * beta(3./2. - H1, H2 + 1./2.) * H2F1_2::operator()(t,s);
	}
};

struct Good : Phi
{
	Good() {};
	inline double operator() (double y, double s, double u)
	{
        double arg1 = MIN(s, u) * (1 - y) / (MAX(s, u) - MIN(s, u) * y);
		double arg2 = MAX(s, u) * (1 - y) / (MAX(s, u) - MIN(s, u) * y);
		if (arg1 > 1) arg1 = 1;
		if (arg2 > 1) arg2 = 1;
		return pow(1 - (MIN(s, u)/MAX(s,u))*y, 2 * H1 - 1)*Phi::operator()(1.0, arg1) * Phi::operator()(1.0, arg2);
	}
};

struct I1 : Good
{
	double s;
	double u;
	I1(double _s, double _u) : s(_s), u(_u) {}

	inline double operator() (double t)
	{
		return	Good::operator()(pow(t, alpha1), s, u) / pow(1 - pow(t, alpha1), 2 * H2 - 1);
	}
};

struct I2 : Good
{
	double s;
	double u;
	I2(double _s, double _u) : s(_s), u(_u) {}

	inline double operator() (double t)
	{
		return	Good::operator()(1 - pow(t, alpha2), s, u) / pow(1 - pow(t, alpha2), 1 + H1 - H2);
	}
};

struct kappa0 
{
	kappa0() {}
	inline double operator() (double s, double u)
	{
        if (abs(s - u) < eps) return  cnst * pow(beta(1.5 - H1, H2 - 0.5), 2) * beta(H2 - H1, 1 + 2 * H1 - 2 * H2);
		double a = 0;
		I1 i1(s, u);
		double b1 = pow(edge, 1 / alpha1);
		double b2 = pow(edge, 1 / alpha2);
		I2 i2(s, u); 
		return cnst * (alpha1 * qgaus(i1, a, b1) + alpha2 * qgaus(i2, a, b2));
	}
};


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mexPrintf("start calculating kappa0");
//declare variables
    mxArray *S_in_m, *U_in_m, *Res_out_m;
    const mwSize *dims;
    double *S, *U, *Res;
    int dimx, dimy, numdims;
    int i,j;

//associate inputs
    S_in_m = mxDuplicateArray(prhs[0]);
    U_in_m = mxDuplicateArray(prhs[1]);

//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];

//associate outputs
    Res_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    
//associate pointers
    S = mxGetPr(S_in_m);
    U = mxGetPr(U_in_m);
    Res = mxGetPr(Res_out_m);
    kappa0 k0;
    double diagValue = k0(0.01,0.01);
    for(i=0;i<dimx;i++)
    {
        for(j=i;j<dimy;j++)
        {
            if( i == j ) 
                Res[i*dimy+j] = diagValue;
            else
                Res[i*dimy+j] = k0(S[i*dimy+j],U[i*dimy+j]); 
        }
        
        for(j=0;j<i;j++)
        {
            Res[i*dimy+j] = Res[j*dimy+i]; 
        }
    }
    mexPrintf("\nend calculating kappa0 \n");
    return;
}