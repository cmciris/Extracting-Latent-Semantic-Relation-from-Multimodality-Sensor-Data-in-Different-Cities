
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

void PosLearn_most_violation (	double *X,
				double *lab_tok,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *F,double *Grad,double *GX);


void PosLearn_most_violation_sequence (	double *X,
					double *lab_tok,
					int T,
					int nbclass,int dim,
					double *w,
					double *f,double *Grad,double *GX);

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
	double *GX;
	
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
	plhs[2] = mxCreateDoubleMatrix(dim ,nbtok , mxREAL); /* GX   [dim x ntok] */
	
	F	= mxGetPr(plhs[0]);
	Grad	= mxGetPr(plhs[1]);
	GX	= mxGetPr(plhs[2]);

	memset(Grad, 0, (nbclass*nbclass+nbclass*dim)*sizeof(double));
	memset(GX  , 0, (nbtok*dim)*sizeof(double));
	
	/* core */
	PosLearn_most_violation	(X,
				lab_tok,
				beg_seq,end_seq,
				nbclass,dim,nbseq,
				w,
				F,Grad,GX);
    return;
}
 
/*=======================================================================*/
void PosLearn_most_violation (	double *X,
				double *lab_tok,
				double *beg_seq,double *end_seq,
				int nbclass,int dim,int nbseq,
				double *w,
				double *F,double *Grad,double *GX) {
	int T0,T,i;
	*F = 0.0;
	for(i=0;i<nbseq;i++) {
		T0 = beg_seq[i] - 1;
		T  = end_seq[i] - beg_seq[i] + 1;
		PosLearn_most_violation_sequence	(X+T0*dim,
						lab_tok+T0,
						T,nbclass,dim,w,
						F,Grad,GX+T0*dim);
	}
}

/*=======================================================================*/
void PosLearn_most_violation_sequence (	double *X,
					double *lab_tok,
					int T,
					int N,int dim,
					double *w,
					double *F,double *Grad,double *GX) {

	int y,yprev,ytrue,d,t,L,L0,idx0,idx1,yBest,yPrevBest,i,nvio;
	double score_correct,score_best,score,coef;
	double **gamma_LR;
	int    **prev_LR;
	double **gamma_RL;
	int **prev_RL;
	double **tab_WF;
	double *wy;
	int *tab_vioT;
	int *tab_vioY;
	int **tab_coef_LR;
	int **tab_coef_RL;

	double *wf;
	/* variable for maxtrix multiplication */
	char *chnT = "T";
	char *chnN = "N";
	int m=T,n=N,p=dim;
	double one = 1.0, zero = 0.0;
	int incx = 1, incy = 1;

	gamma_LR = (double**)malloc(N*sizeof(double*));
	prev_LR  = (int**)malloc(N*sizeof(int*));
	gamma_RL = (double**)malloc(N*sizeof(double*));
	prev_RL  = (int**)malloc(N*sizeof(int*));
	tab_WF = (double**)malloc(N*sizeof(double*));
	tab_vioT = (int*)malloc(T*sizeof(int));
	tab_coef_LR = (int**)malloc(N*sizeof(int*));
	tab_coef_RL = (int**)malloc(N*sizeof(int*));

	/* ==== precompute local energy === */
/*        printf("precompute local energy\n");mexEvalString("drawnow");*/

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


	for(y=0;y<N;y++) {
		gamma_LR[y]	= (double*)malloc(T*sizeof(double));
		prev_LR [y]	= (int*)malloc(T*sizeof(int));
		gamma_RL[y]	= (double*)malloc(T*sizeof(double));
		prev_RL [y]	= (int*)malloc(T*sizeof(int));
		tab_WF[y]	= wf + y*T;
		tab_coef_LR[y]	= (int*)malloc(T*sizeof(int));
		tab_coef_RL[y]	= (int*)malloc(T*sizeof(int));
	}


	/* ==== Viterbi left right ==== */

	t=0;
	for(y=0;y<N;y++) {
		gamma_LR[y][t] = tab_WF[y][t];
	}
	y = lab_tok[t]-1;
	score_correct = gamma_LR[y][t];
	/* t=1..T-1 */
	for(t=1;t<T;t++) {
		for (y=0;y<N;y++) {
			yPrevBest = 0;
			for (yprev=1;yprev<N;yprev++)
				if (gamma_LR[yPrevBest][t-1] + w[y*N + yPrevBest] < gamma_LR[yprev][t-1] + w[y*N + yprev])
					yPrevBest = yprev;
			gamma_LR[y][t] = gamma_LR[yPrevBest][t-1] + w[y*N + yPrevBest] + tab_WF[y][t];
			prev_LR[y][t] = yPrevBest;
		}
		y = lab_tok[t]-1;
		yprev = lab_tok[t-1]-1;
		score_correct += w[y*N + yprev] + tab_WF[y][t];
	}

	/* ==== Viterbi right left ==== */

	t=T-1;
	for(y=0;y<N;y++) {
		gamma_RL[y][t] = tab_WF[y][t];
	}

	/* t=T-2..0 */
	for(t=T-2;t>=0;t--) {
		for (y=0;y<N;y++) {
			yPrevBest = 0; /* right to left */
			for (yprev=1;yprev<N;yprev++)
				if (gamma_RL[yPrevBest][t+1] + w[yPrevBest*N + y] < gamma_RL[yprev][t+1] + w[yprev*N + y]) 
					yPrevBest = yprev;
			gamma_RL[y][t] = gamma_RL[yPrevBest][t+1] + w[yPrevBest*N + y] + tab_WF[y][t];
			prev_RL[y][t] = yPrevBest;
		}
	}

	/* ==== compute most violation at each position tab_vioT ==== */

	nvio = 0;
	for(t=0;t<T;t++) {
		score_best = INF; /* minus infinity */
		ytrue = (int)(lab_tok[t]-1);
		for (y=0;y<N;y++)
			if (y!=ytrue)	{
				score = gamma_LR[y][t] + gamma_RL[y][t] - tab_WF[y][t];
				if (score > score_best) {
					score_best = score;
					tab_vioT[t] = y;
				}
			}
		
		if (score_correct-score_best<1) {
			*F += 1.0 + score_best - score_correct;
			nvio ++;
		} else
			tab_vioT[t] = -1;
	}
			
	/* ==== precompute gradient coeficients === */

	for (t=0;t<T;t++) {
		for (y=0;y<N;y++) {
			tab_coef_LR[y][t] = 0;
			tab_coef_RL[y][t] = 0;
		}
		if (tab_vioT[t]>=0) {
			tab_coef_LR[tab_vioT[t]][t] = 1;
			tab_coef_RL[tab_vioT[t]][t] = 1;
		}
	}

	for (t=0;t<T-1;t++)
		for (y=0;y<N;y++)
			if(tab_coef_RL[y][t]!=0) {
				yprev = prev_RL[y][t];
				tab_coef_RL[yprev][t+1] += tab_coef_RL[y][t];
			}

	for(t=T-1;t>0;t--)
		for (y=0;y<N;y++)
			if(tab_coef_LR[y][t]!=0) {
				yprev = prev_LR[y][t];
				tab_coef_LR[yprev][t-1] += tab_coef_LR[y][t];
			}

	/* ==== update gradient === */

	for(t=0;t<T;t++) {
		L0 = t*dim;
		for(y=0;y<N;y++) {
			/* ==== local  gradient ==== */
			coef = tab_coef_LR[y][t]+tab_coef_RL[y][t]; 
			if (y==tab_vioT[t]) 
				coef -= 1.0; /* note that the coef at [yvio][t] is counted 2 times */
			if ( coef != 0) {
				idx0 = N*N + y*dim;
#if FLAG_BLAS
				daxpy_(&dim,&coef,  X+L0,&incx,Grad+idx0,&incy);
				daxpy_(&dim,&coef,w+idx0,&incx,GX+L0    ,&incy);
#else
				for (d=0;d<dim;d++) {
				  Grad[idx0 + d] += coef * X[L0+d];
				  GX[L0 + d] += coef * w[idx0+d];
				}
#endif
			}

			/* ==== transition gradient ==== */
			if ((t>0)&&(tab_coef_LR[y][t]>0)) {
				yprev = prev_LR[y][t];
				Grad[y*N + yprev] += tab_coef_LR[y][t];
			}

			if ((t<T-1)&&(tab_coef_RL[y][t]>0)) {
				yprev = prev_RL[y][t];
				Grad[yprev*N + y] += tab_coef_RL[y][t];
			}
		}

		/* gradient for true labelling term */
		ytrue = (int)(lab_tok[t]-1);
		idx0 = N*N + ytrue*dim;
		coef = -nvio;
#if FLAG_BLAS
		daxpy_(&dim,&coef,X+L0,&incx,Grad+idx0,&incy);
		daxpy_(&dim,&coef,w+idx0,&incx,GX+L0    ,&incy);
#else
		for (d=0;d<dim;d++) {
		  Grad[idx0 + d] += coef * X[L0+d];
		  GX[L0 + d] += coef * w[idx0+d];
		}
#endif

		if (t>0) {
			yprev = (int)(lab_tok[t-1]-1);
			Grad[ytrue*N + yprev] -= nvio;
		}
	}

	for(y=0;y<N;y++) {
		free(gamma_LR[y]);
		free(prev_LR [y]);
		free(gamma_RL[y]);
		free(prev_RL [y]);
		free(tab_coef_LR[y]);
		free(tab_coef_RL[y]);
	}
	free(wf);
	free(gamma_LR);
	free(prev_LR );
	free(gamma_RL);
	free(prev_RL );
	free(tab_WF);
	free(tab_coef_LR);
	free(tab_coef_RL);
	free(tab_vioT);
}

/*============================================================================== ===
*/
double scalar_product(double *X,double *w,int dim) {
	double sum;
	int i;
	sum = 0.0;
	for(i=0;i<dim;i++)
			sum += w[i] * X[i];
	return sum;
}
