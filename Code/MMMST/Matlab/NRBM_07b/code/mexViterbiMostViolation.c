
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

void viterbi_most_violation (	double *X,
				double *lab_tok,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *F,double *Grad);

void viterbi_most_violation_sequence (	double *X,
					double *lab_tok,
					int T,
					int nbclass,int dim,
					double *w,
					double *f,double *Grad);

double scalar_product(double *X,double *w,int dim);

/* gateway */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* gateway variables */
	/* all in double type for compatibility with matlab */
	double *X;
	double *lab_tok;
	double *beg_seq,*end_seq;
	int nbclass,dim;
	double *w;
	
	/* workspace variable */
	int nbseq;
	int dimW;

	/* output variable */
	double *F;
	double *Grad;
	
	/* input variables */
	X	= (double*)mxGetData(prhs[0]);
	lab_tok	= (double*)mxGetData(prhs[1]);
	beg_seq	= (double*)mxGetData(prhs[2]); nbseq = mxGetM(prhs[2]); /* N lines, M columns */
	end_seq	= (double*)mxGetData(prhs[3]);
	nbclass	= (int)(*(double*)mxGetData(prhs[4]));
	dim	= (int)(*(double*)mxGetData(prhs[5]));
	w	= (double*)mxGetData(prhs[6]);
	
	dimW = nbclass*nbclass + dim*nbclass;

	/* output variable */
	plhs[0] = mxCreateDoubleMatrix(1 ,1 , mxREAL); /* F   [1x1] */
	plhs[1] = mxCreateDoubleMatrix(1 ,dimW , mxREAL); /* Grad   [1xdimW] */
	
	F   = mxGetPr(plhs[0]);
	Grad = mxGetPr(plhs[1]);
	
	/* core */
	viterbi_most_violation	(X,
				lab_tok,
				beg_seq,end_seq,
				nbclass,dim,nbseq,
				w,
				F,Grad);
    return;
}
 
void viterbi_most_violation (	double *X,
				double *lab_tok,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *F,double *Grad) {
	int T0,T,i;
	*F = 0.0;
	memset(Grad, 0, (nbclass*nbclass+nbclass*dim)*sizeof(double));
	for(i=0;i<nbseq;i++) {
		T0 = beg_seq[i] - 1;
		T  = end_seq[i] - beg_seq[i] + 1;
		viterbi_most_violation_sequence	(X+T0*dim,
						lab_tok+T0,
						T,nbclass,dim,w,
						F,Grad);
	}
}

void viterbi_most_violation_sequence (	double *X,
					double *lab_tok,
					int T,
					int N,int dim,
					double *w,
					double *F,double *Grad) {
	int y,yprev,d,t,L,L0,idx0,idx1,yBest,yPrevBest,i;
	double score_correct;
	double **gamma;
	int    **prev;
	double *wy;

	gamma = (double**)malloc(N*sizeof(double*));
	prev  = (int**)malloc(N*sizeof(int*));

	for(y=0;y<N;y++) {
		gamma[y] = (double*)malloc(T*sizeof(double));
		prev [y] = (int*)malloc(T*sizeof(int));
	}

	
	t=0;L0 = t*dim;
	for(y=0;y<N;y++) {
		wy = w + N*N + y*dim;
		gamma[y][t] = scalar_product(X+L0,wy,dim) + (y!=(int)(lab_tok[t]-1));
	}
	score_correct = gamma[(int)lab_tok[t]-1][t];


	/* t=1..T-1 */
	for(t=1;t<T;t++) {
		L0 = t*dim;
		for (y=0;y<N;y++) {
			wy = w + N*N + y*dim;
			yPrevBest = 0;
			for (yprev=1;yprev<N;yprev++)
				if (gamma[yPrevBest][t-1] + w[y*N + yPrevBest] < gamma[yprev][t-1] + w[y*N + yprev]) {
					yPrevBest = yprev;
				}
			gamma[y][t] = gamma[yPrevBest][t-1] + w[y*N + yPrevBest] + scalar_product(X+L0,wy,dim) + (y!=(int)(lab_tok[t]-1));
			prev[y][t] = yPrevBest;
		}
		y = lab_tok[t]-1;
		yprev = lab_tok[t-1]-1;
		wy = w + N*N + y*dim;
		score_correct += w[y*N + yprev] + scalar_product(X+L0,wy,dim);
	}

	/* decoding : allow any state to be end state*/
	yBest = 0;t=T-1;
	for (y=1;y<N;y++) 
		if (gamma[yBest][t] < gamma[y][t])
			yBest = y;

	if (gamma[yBest][t]>score_correct) {
		*F += gamma[yBest][t]-score_correct;
		y = yBest;
		for(t=T-1;t>=0;t--) {
			yprev = prev[y][t];
			idx0 = N*N + y*dim;
			idx1 = N*N + (lab_tok[t]-1)*dim;
			L0 = t*dim;

			for (i=0;i<dim;i++) {
				Grad[idx0+i] += X[L0+i];
				Grad[idx1+i] -= X[L0+i];
			}
			if (t>=1) {
				Grad[y*N + yprev] += 1.0;
				Grad[(int)(lab_tok[t]-1)*N + (int)(lab_tok[t-1]-1)] -= 1.0;
			}
			y = yprev;
		}
	}

	for(y=0;y<N;y++) {
		free(gamma[y]);
		free(prev [y]);
	}
	free(gamma);
	free(prev);
}

double scalar_product(double *X,double *w,int dim) {
	double sum;
	int i;
	sum = 0.0;
	for(i=0;i<dim;i++)
			sum += w[i] * X[i];
	return sum;
}
