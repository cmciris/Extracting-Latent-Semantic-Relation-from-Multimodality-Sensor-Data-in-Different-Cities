
#include "mex.h"

/* INPUT
     H is [M x N] matrix
     a is [D x 1] vector
     idx is [D x 1] vector
   OUTPUT
     w is [M x 1] vector
*/     

/* -------------------------------------------------------------------
 Main MEX function - interface to Matlab.
-------------------------------------------------------------------- */
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
  double *H;
  double *a;
  double *idx;
  double *w;
  int N,M,i,j,Da,Di;
 
  /*------------------------------------------------------------------- */
  /* Take input arguments                                               */
  /*------------------------------------------------------------------- */

  if( nrhs != 3) mexErrMsgTxt("Incorrect number of input arguments.");

  /* matrix H */
  H = mxGetPr(prhs[0]);   
  M = mxGetM(prhs[0]);     
  N = mxGetN(prhs[0]);

  a = mxGetPr(prhs[1]);
  Da = mxGetM(prhs[1]);
  idx = mxGetPr(prhs[2]);
  Di = mxGetM(prhs[2]);
/*
  
*/
  if ((Di!=Da)||(Da>N)) {
	mexPrintf("H[%d x %d], a [%d x ?] , idx [%d x ?]\n",M,N,Da,Di);
	mexErrMsgTxt("conflict dimension");
  }
  
  /* output "solution" vector alpha [dim x 1] */
  plhs[0] = mxCreateDoubleMatrix(M,1,mxREAL);
  w = mxGetPr(plhs[0]);

  for (j=0;j<M;j++) {
	w[j] = 0.0;
	for (i=0;i<Da;i++)
		w[j] += H[(int)(idx[i]-1)*M+j] * a[i];
  }

}
