
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
#define FLAG_BLAS 0

void viterbi_inference (	double *X,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *lab_tok);

void viterbi_inference_sequence (	double *X,
					int T,
					int nbclass,int dim,
					double *w,
					double *lab_tok);

double scalar_product(double *X,double *w,int dim);

/* gateway */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* gateway variables */
	/* all in double type for compatibility with matlab */
	double *X;
	double *beg_seq,*end_seq;
	int nbclass,dim;
	double *w;
	
	/* workspace variable */
	int nbseq,nbtok;
	int dimW;

	/* output variable */
	double *lab_tok;
	
	/* input variables */
	X	= (double*)mxGetData(prhs[0]);
	beg_seq	= (double*)mxGetData(prhs[1]); nbseq = mxGetM(prhs[1]); /* N lines, M columns */
	end_seq	= (double*)mxGetData(prhs[2]);
	nbclass	= (int)(*(double*)mxGetData(prhs[3]));
	dim	= (int)(*(double*)mxGetData(prhs[4]));
	w	= (double*)mxGetData(prhs[5]);
	
	dimW    = nbclass*nbclass + dim*nbclass;
	nbtok   = end_seq[nbseq-1];

	/* output variable */
	plhs[0] = mxCreateDoubleMatrix(nbtok,1 , mxREAL); /* F   [1x1] */
	
	lab_tok   = mxGetPr(plhs[0]);
	
	/* core */
	viterbi_inference	(X,
				beg_seq,end_seq,
				nbclass,dim,nbseq,
				w,
				lab_tok);
    return;
}
 
void viterbi_inference (	double *X,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *lab_tok) {
	int T0,T,i;
	for(i=0;i<nbseq;i++) {
		T0 = beg_seq[i] - 1;
		T  = end_seq[i] - beg_seq[i] + 1;
		viterbi_inference_sequence	(X+T0*dim,
						T,
						nbclass,dim,w,
						lab_tok+T0);
	}
}

void viterbi_inference_sequence (	double *X,
					int T,
					int N,int dim,
					double *w,
					double *lab_tok) {
	int y,yprev,d,t,L,L0,idx0,idx1,yBest,yPrevBest,i;
	double **gamma;
	int    **prev;
	double *wy;
	double *wf;
	/* variable for maxtrix multiplication */
	char *chnT = "T";
	char *chnN = "N";
	mwSignedIndex m=T,n=N,p=dim;
/*	mwSignedIndex m=(mwSignedIndex)T,n=(mwSignedIndex)N,p=(mwSignedIndex)dim;*/
	double one = 1.0, zero = 0.0;


	wf	= (double*)malloc(N*T*sizeof(double));
#if FLAG_BLAS
	dgemm_ ( chnT, chnN, &m, &n, &p, &one, X , &p, w+N*N, &p, &zero, wf, &m);
#else
	for(y=0;y<N;y++) {
		wy = w + N*N + y*dim;
		for(t=0;t<T;t++)
			wf[y*T+t]= scalar_product(X+t*dim,wy,dim);
	}
#endif

	gamma = (double**)malloc(N*sizeof(double*));
	prev  = (int**)malloc(N*sizeof(int*));

	for(y=0;y<N;y++) {
		gamma[y] = (double*)malloc(T*sizeof(double));
		prev [y] = (int*)malloc(T*sizeof(int));
	}

	
	t=0;L0 = 0;
	for(y=0;y<N;y++) {
		wy = w + N*N + y*dim;
		gamma[y][t] = wf[y*T+t];/*scalar_product(X+L0,wy,dim);*/
	}


	/* t=1..T-1 */
	for(t=1;t<T;t++) {
		for (y=0;y<N;y++) {
			wy = w + N*N + y*dim;
			yPrevBest = 0;
			for (yprev=1;yprev<N;yprev++)
				if (gamma[yPrevBest][t-1] + w[y*N + yPrevBest] < gamma[yprev][t-1] + w[y*N + yprev]) {
					yPrevBest = yprev;
				}
			gamma[y][t] = gamma[yPrevBest][t-1] + w[y*N + yPrevBest] + wf[y*T+t];/*scalar_product(X+L0,wy,dim);*/
			prev[y][t] = yPrevBest;
		}
	}
	/* decoding : allow any state to be end state*/
	yBest = 0;t=T-1;
	for (y=1;y<N;y++) 
		if (gamma[yBest][t] < gamma[y][t])
			yBest = y;

	y = yBest;
	for(t=T-1;t>=0;t--) {
		lab_tok[t] = y+1;
		y = prev[y][t];
	}

	for(y=0;y<N;y++) {
		free(gamma[y]);
		free(prev [y]);
	}
	free(gamma);
	free(prev);

	free(wf);
}

double scalar_product(double *X,double *w,int dim) {
	double sum;
	int i;
	sum = 0.0;
	for(i=0;i<dim;i++)
			sum += w[i] * X[i];
	return sum;
}
