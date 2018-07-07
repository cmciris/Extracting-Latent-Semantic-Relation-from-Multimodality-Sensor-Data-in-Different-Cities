
//  Forward backward procedure 
//  Copyright (C) 2012, by Trinh-Minh-Tri Do minhtrido@gmail.com
//  
//    This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//  
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//  
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <float.h>
#include <string.h> 
#include <math.h>
#include "time.h"
 
#define INF -INFINITY
#define SMALL_PROBABILITY 1e-6
#define FLAG_BLAS 0

void forward_backward (	double *X,
				double *lab_tok,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *F,double *Grad);

void forward_backward_sequence (	double *X,
					double *lab_tok,
					int T,
					int nbclass,int dim,
					double *w,
					double *f,double *Grad);

double expmx(double x);
double forward(int N,int T,double **WF,double *transScore,double **alpha);
void backward(int N,int T,double **WF,double *transScore,double **beta);
double sumlog(const double *v, int n);
double scalar_product(double *X,double *w,int dim);

/* gateway */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* gateway variables */
	/* all in double type for compatibility with matlab */
	double *X;
	double *lab_tok;
	double *beg_seq,*end_seq;
	int nbclass,dim,sizeX;
	double *w;
	
	/* workspace variable */
	int nbseq,nbtok;
	int dimW;

	/* output variable */
	double *F;
	double *Grad;
	
	/* input variables */
	X	= (double*)mxGetData(prhs[0]); nbtok = mxGetM(prhs[0]); 
	lab_tok	= (double*)mxGetData(prhs[1]);
	beg_seq	= (double*)mxGetData(prhs[2]); nbseq = mxGetM(prhs[2]); /* N lines, M columns */
	end_seq	= (double*)mxGetData(prhs[3]);
	nbclass	= (int)(*(double*)mxGetData(prhs[4]));
	dim	= (int)(*(double*)mxGetData(prhs[5]));
	w	= (double*)mxGetData(prhs[6]);
	
	dimW = nbclass*nbclass + dim*nbclass;

	/* output variable */
	plhs[0] = mxCreateDoubleMatrix(1 ,1 , mxREAL);       /* F   [1x1] */
	plhs[1] = mxCreateDoubleMatrix(1 ,dimW , mxREAL);    /* Grad   [1 x dimW] */
	
	F	= mxGetPr(plhs[0]);
	Grad	= mxGetPr(plhs[1]);

	memset(Grad, 0, (nbclass*nbclass+nbclass*dim)*sizeof(double));
	
	/* core */
	forward_backward	(X,
				lab_tok,
				beg_seq,end_seq,
				nbclass,dim,nbseq,
				w,
				F,Grad);
    return;
}
 
void forward_backward (	double *X,
				double *lab_tok,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *F,double *Grad) {
	int T0,T,i;
	*F = 0.0;
	for(i=0;i<nbseq;i++) {
		T0 = beg_seq[i] - 1;
		T  = end_seq[i] - beg_seq[i] + 1;
		forward_backward_sequence	(X+T0*dim,
						lab_tok+T0,
						T,nbclass,dim,w,
						F,Grad);
	}
}

void forward_backward_sequence (	double *X,
					double *lab_tok,
					int T,
					int N,int dim,
					double *w,
					double *F,double *Grad) {
static  double time_wf=0.0,time_Forward=0.0,time_Backward=0.0,time_grad=0.0,time_total=0.0,last_time=0.0;
	double mytime=clock();

	int y,yprev,d,t,L,L0,idx0,idx1,yBest,yPrevBest,i;
	double score_correct;
	double **alpha,**beta,**coef;
	double *wy,*gy;
	double **WF;
	int    dimNN;
	double Z;
	double *wf;
	/* variable for maxtrix multiplication */
	char *chnT = "T";
	char *chnN = "N";
	mwSignedIndex m=(mwSignedIndex)T,n=(mwSignedIndex)N,p=(mwSignedIndex)dim;
	double one = 1.0, zero = 0.0;
	int incx = 1, incy = 1;


	WF	= (double**)malloc(N*sizeof(double*));
	wf	= (double*)malloc(N*T*sizeof(double));
	alpha	= (double**)malloc(N*sizeof(double*));
	beta	= (double**)malloc(N*sizeof(double*));
	coef    = (double**)malloc(N*sizeof(double*));

#if FLAG_BLAS
	dgemm_ ( chnT, chnN, &m, &n, &p, &one, X , &p, w+N*N, &p, &zero, wf, &m);
#else
	for(y=0;y<N;y++) {
		wy = w + N*N + y*dim;
		for(t=0;t<T;t++)
			wf[y*T+t]= scalar_product(X+t*dim,wy,dim);
	}
#endif
	
	for(y=0;y<N;y++) {
		alpha[y] = (double*)malloc(T*sizeof(double));
		beta[y] = (double*)malloc(T*sizeof(double));
		coef[y] = (double*)malloc(T*sizeof(double));
		WF[y]   = wf + y*T; 
	}

	/* preparing local score */
	for(t=0;t<T;t++) {
		for(y=0;y<N;y++) {
			coef[y][t] = 0.0;
		}
	}

	time_wf += clock()-mytime;mytime=clock();
	Z = forward(N,T,WF,w,alpha); /* transition matrix is at the begining of w*/
	time_Forward += clock()-mytime;mytime=clock();
	backward(N,T,WF,w,beta);
	time_Backward += clock()-mytime;mytime=clock();

	*F += Z;

	for(t=0;t<T;t++) {
		L0 = t*dim;
		for (y=0;y<N;y++)
			coef[y][t] = expmx(Z-alpha[y][t]-beta[y][t]); /* marginal probability */
		y = lab_tok[t]-1;
		coef[y][t] -= 1.0;
		*F -= WF[y][t];
		for (y=0;y<N;y++) {
			if ((y==lab_tok[t]-1)||(coef[y][t]>SMALL_PROBABILITY)) {
				/* local feature */
				idx0 = N*N + y*dim;
#if FLAG_BLAS
				daxpy_(&dim,coef[y]+t,X+L0,&incx,Grad+idx0,&incy);
#else
				for (d=0;d<dim;d++)
				  Grad[idx0 + d] += coef[y][t] * X[L0+d];
#endif
				
				/* transition feature */
				if (t>=1) {
					for (yprev=0;yprev<N;yprev++)
						Grad[yprev + y*N] +=  expmx(Z-alpha[yprev][t-1]-beta[y][t]-WF[y][t]-w[yprev + y*N]);
				}
			}
		}
		if (t>=1) {
			y = lab_tok[t]-1;yprev = lab_tok[t-1]-1;
			Grad[yprev + y*N] -= 1.0;
			*F -= w[yprev + y*N];
		}
	}

	time_grad += clock()-mytime;mytime=clock();
	time_total = time_wf + time_Forward + time_Backward + time_grad;
	if (((time_total - last_time) / CLOCKS_PER_SEC > 10.0)&&(0)) {
		printf("time_wf		=%.2f\n",time_wf   / CLOCKS_PER_SEC);
		printf("time_Forward	=%.2f\n",time_Forward      / CLOCKS_PER_SEC);
		printf("time_Backward	=%.2f\n",time_Backward    / CLOCKS_PER_SEC);
		printf("time_grad	=%.2f\n",time_grad / CLOCKS_PER_SEC);
		printf("time_total  =%.2f\n====================\n",time_total / CLOCKS_PER_SEC);
		last_time = time_total;
		mexEvalString("drawnow;");
	}


	for(y=0;y<N;y++) {
		free(alpha[y]);
		free(beta[y]);
		free(coef[y]);
	}
	free(alpha);
	free(beta);
	free(coef);
	free(WF);
	free(wf);
}

/*=================================================================================
*/
double scalar_product(double *X,double *w,int dim) {
	double sum;
	int i;
	sum = 0.0;
	for(i=0;i<dim;i++)
			sum += w[i] * X[i];
	return sum;
}

/*=================================================================================
// forward proceduce 
// N states
// sequence of lenght T
// Using precalculated local score : WF (N,T) ; WF(n,t) = W_n * F_t + Delta(n,y_t) <= margin included
//=================================================================================*/

double forward(int N,int T,double **WF,double *transScore,double **alpha) {
	int s,t,y,k,sprev;
	static double *tmp=NULL;
	if(tmp==NULL)
		tmp = (double*)malloc(1000*sizeof(double));
	t=0;
	for(s=0;s<N;s++) {
		alpha[s][t] = WF[s][t];
	}

	for(t=1;t<T;t++) {
		for(s=0;s<N;s++) {
			for(k=0;k<N;k++) {
				tmp[k] = WF[s][t] + alpha[k][t-1] + transScore[k+s*N];
			}
			alpha[s][t] = sumlog(tmp,N);
		}
	}

	for(s=0;s<N;s++)
		tmp[s] = alpha[s][T-1];

/*	free(tmp);*/
	return sumlog(tmp,N);
}

/*=================================================================================
// backward proceduce (note that there is no prior in beta)
// N states
// sequence of lenght T
// Using precalculated local score : WF (N,T) ; WF(n,t) = W_n * F_t + Delta(n,y_t) <= margin included
//=================================================================================*/

void backward(int N,int T,double **WF,double *transScore,double **beta) {
	int s,t,y,k,snext;
	static double *tmp=NULL;
	if(tmp==NULL)
		tmp = (double*)malloc(1000*sizeof(double));
	t=T-1;
	for(s=0;s<N;s++)
		beta[s][t] = 0;

	for(t=T-2;t>=0;t--) {
		for(s=0;s<N;s++) {
			for(k=0;k<N;k++) {
				tmp[k] = WF[k][t+1] + beta[k][t+1] + transScore[s+k*N];
			}
			beta[s][t] = sumlog(tmp,N);
		}
	}
/*	free(tmp);*/
}

/*=================================================================================*/
double sumlog(const double *v, int n) {
	int i;
	double s = 0;
	double m;
	if (n==0)
		return INF;
	if (n==1)
		return v[0];

	m = v[0];
	for (i=1; i<n; i++)
		m = (m>v[i]?m:v[i]);

	for (i=0; i<n; i++)
		s += expmx(m-v[i]);

	return m + log(s);
}

/*=================================================================================*/
double expmx(double x)
{
  /* fast approximation of exp(-x) for x positive */
# define A0   (1.0)
# define A1   (0.125)
# define A2   (0.0078125)
# define A3   (0.00032552083)
# define A4   (1.0172526e-5) 
  if (x < 13.0) 
    {
      double y;
      y = A0+x*(A1+x*(A2+x*(A3+x*A4)));
      y *= y;
      y *= y;
      y *= y;
      y = 1/y;
      return y;
    }
  return 0;
}
