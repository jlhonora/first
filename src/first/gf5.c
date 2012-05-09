/*=========================================
 * gf5.c : C-version of GradientFunction5.m (Jose Luis)
 *
 * 2011 06 09 Only the cost function (no gradients)
 * 2011 06 13 Now it includes the gradients for rho
 *
 * Pablo Irarrazaval
 *===========================================*/

#include "mex.h" /* matlab includer */

#ifdef USE_MATH
  #define _USE_MATH_DEFINES /* it will define math constants (pi) */
  #include "math.h" /* mathematical operations */
#else
  #define SQ(x) ((x)*(x))
  #define M_PI (3.1415926535897931)
#endif

/* ==========================================
 * The computational function
 *
 * INPUTS
 * N : number of x positions (size: 1x1)
 * M : number of species (size: 1x1)
 * Nacq : number of echoes (size: 1x1)
 * x : vector with species, real and imaginary, field map and R2
 *     (size: 1 x 2N + 2N + ... + 2N (M times) + N + N)
 * kx : array with kx positions (size: Nacq x N)
 * kf : array with time map (size: Nacq x N)
 * Mre, Mim : arrays with real and imaginary part of raw data
 *            (size: Nacq x N each)
 * fo : frequency offsets (size: 1xM, number of species)
 * xpos : x positions (size: 1xN)
 *
 * OUTPUT
 * f : cost function value (pointer)
 * ========================================== */
void gf5(double *x, double *kx, double *kf, double *Mre, double *Mim,
		 double *fo, double *xpos, mwSize N, mwSize M, mwSize Nacq,
		 double *f, double *g)	 
{
	mwSize m; /* counter for species */
	mwSize n; /* counter for x positions */
	mwSize q; /* counter for k positions */
	mwSize e; /* counter for echoes */
	double ff; /* accumulative value for cost function */
	double MEr,MEi; /* temporal variables */
	double *pMEr,*pMEi;
	double mr,mi,ang,cs,sn;
	double re_rho_m, im_rho_m, re_gphi, im_gphi; 
	double Dqr,Dqi;
	double tmp_kx,tmp_kf,tmp_fo,tmp_fo_kf;
	double twopi;
	mwSize twoN,twoNm,twoNM;
	int i;
	
// 	printf("C: N = %d M = %d Nacq = %d\n",N,M,Nacq);
	
	/* A word about col order of matlab versus row order of C
	 * If size(a) = [M N] in matlab notation, then a(m,n) is
	 * a[m + n*M] in C                                       */
	
	twopi = 2*M_PI;
	twoN = 2*N;
	twoNM = twoN*M;
	
	/* initialize */
	ff = 0;
	
	if( g==NULL ) { /* only cost function */
		
		for (e=0; e<Nacq; e++) { /* For every echo */
			
			for (q=0; q<N; q++) { /* For every k position */
				
				tmp_kx = kx[e+q*Nacq];
				tmp_kf = kf[e+q*Nacq];
				
				MEr = 0; MEi = 0;
				
				for (m=0; m<M; m++) { /* For every species */
					
					tmp_fo = fo[m];
					twoNm = twoN*m;
					tmp_fo_kf = tmp_fo*tmp_kf;
					
					for (n=0; n<N; n++) { /* For every x positions */
						
						mr = x[twoNm+n];
						mi = x[twoNm+N+n];
						ang = -twopi*(xpos[n]*tmp_kx - tmp_fo_kf - x[twoNM+n]*tmp_kf);
						                             
						cs = cos(ang);
						sn = sin(ang);
						
						MEr = MEr + mr*cs - mi*sn;
						MEi = MEi + mi*cs + mr*sn;
						
					}
				}
				
				#ifdef USE_MATH
						ff = ff + pow(Mre[e+q*Nacq]-MEr, 2) + pow(Mim[e+q*Nacq]-MEi, 2);
				#else
						ff = ff + SQ(Mre[e+q*Nacq]-MEr) + SQ(Mim[e+q*Nacq]-MEi);
				#endif
						
			}
		}
	}
	else { /* cost function AND gradient */
		
		for (n=0; n<twoNM+twoN; n++) { /* makes sure g starts in zero */
			g[n] = 0; }
		/* allocates memory for temporal arrays */
		pMEr = (double *)malloc(twoNM*sizeof(double));
		pMEi = (double *)malloc(twoNM*sizeof(double));
		
		for (e=0; e<Nacq; e++) { /* For every echo */
			
			for (q=0; q<N; q++) { /* For every k position */
				
				tmp_kx = kx[e+q*Nacq];
				tmp_kf = kf[e+q*Nacq];
				
				MEr = 0; MEi = 0;
				for(i=0; i<twoNM; i++){ /* initialize temporal array */
					pMEr[i] = 0;
					pMEi[i] = 0; }
				
				for (m=0; m<M; m++) { /* For every species */
					
					tmp_fo = fo[m];
					twoNm = twoN*m;
					tmp_fo_kf = tmp_fo*tmp_kf;
					
					for (n=0; n<N; n++) { /* For every x positions */
						
						mr = x[twoNm+n];
						mi = x[twoNm+N+n];
						ang = -twopi*(xpos[n]*tmp_kx - tmp_fo_kf - x[twoNM+n]*tmp_kf);
						cs = cos(ang); sn = sin(ang);
						
						MEr = MEr + mr*cs - mi*sn;
						MEi = MEi + mi*cs + mr*sn;
						
					}
				}
				
				Dqr = Mre[e+q*Nacq]-MEr;
				Dqi = Mim[e+q*Nacq]-MEi;
				#ifdef USE_MATH
						ff = ff + pow(Dqr, 2) + pow(Dqi, 2);
				#else
						ff = ff + SQ(Dqr) + SQ(Dqi);
				#endif
				
				Dqr = Dqr*2; Dqi = Dqi*2; /* multiply by 2 only once */
				for (m=0; m<M; m++) { /* For every species */
					
					tmp_fo = fo[m];
					twoNm = twoN*m;
					tmp_fo_kf = tmp_fo*tmp_kf;
					
					for (n=0; n<N; n++) { /* For every x positions (real and imaginary) */
						
						re_gphi = 0;
						im_gphi = 0;

						ang = -twopi*(xpos[n]*tmp_kx - tmp_fo_kf - x[twoNM+n]*tmp_kf);
						cs = cos(ang); sn = sin(ang);
						
						g[twoNm+n] = g[twoNm+n] - Dqr*cs - Dqi*sn;
						g[twoNm+N+n] = g[twoNm+N+n] + Dqr*sn - Dqi*cs;

						re_rho_m = x[twoNm+n];
						im_rho_m = x[twoNm+N+n];

						re_gphi = re_gphi - twopi*tmp_kf*(im_rho_m*cs + re_rho_m*sn);
						im_gphi = im_gphi + twopi*tmp_kf*(re_rho_m*cs - im_rho_m*sn);
						g[twoNM+n] = g[twoNM+n] - Dqr*re_gphi - Dqi*im_gphi; // Phi
					}
				}						
			}
		}
	}
	
	*f = ff;

	/* testing the routine */
// 	printf("x =\n");
// 	for (n=0; n<(M*N*2+2*N); n++) {
// 		printf("%.4f ",x[n]); }
// 	printf("\n");
// 	
// 	printf("kx =\n");
// 	for (m=0; m<Nacq; m++) {
// 		for (n=0; n<N; n++) {
// 			printf("%.4f ", kx[n*Nacq + m]); }
// 		printf("\n");
// 	}
// 	
// 	printf("kf =\n");
// 	for (m=0; m<Nacq; m++) {
// 		for (n=0; n<N; n++) {
// 			printf("%.4f ", kf[n*Nacq + m]); }
// 		printf("\n");
// 	}
//     
//     printf("M =\n");
// 	for (m=0; m<Nacq; m++) {
// 		for (n=0; n<N; n++) {
// 			printf("%.4f+i%.4f ", Mre[n*Nacq + m],Mim[n*Nacq + m]); }
// 		printf("\n");
// 	}
}

/* ==========================================
 * The gateway function
 * ========================================== */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	
/* function f = GradientFunction5(x, kx, kf, M, frequency_offset, x_positions) */
	
	/* variables to pass into the computational function */
	double *x,*kx,*kf,*Mre,*Mim,*fo,*xpos; /* directly from Matlab script */
	mwSize N,M,Nacq; /* size descriptors */
	double *f,*g; /* the values returned by the computational function */
	
	/* Working variables */
	mwSize i,j;

    /* check for proper number of arguments */
    //if(nrhs!=6) {
    //    mexErrMsgTxt("gf5: Six inputs required.");   }
    //if( (nlhs!=1) && (nlhs!=2) ) {
    //    mexErrMsgTxt("gf5: One or two outputs required.");    }
	
	/* Process first input: x */
	//if( mxIsComplex(prhs[0]) ) { /* Warning for complez */
	//	mexWarnMsgTxt("gf5 WARNING: Imaginary part of x (input 1) ignored."); }
	x = mxGetPr(prhs[0]); /* creates pointer */
	i = mxGetM(prhs[0]); /* get size */
	//if( i != 1 ) {
	//	mexErrMsgTxt("gf5: x must be a row vector."); }
	
	/* Process second input: kx */
	//if( mxIsComplex(prhs[1]) ) { /* Warning for complez */
	//	mexWarnMsgTxt("gf5 WARNING: Imaginary part of kx (input 2) ignored."); }
	kx = mxGetPr(prhs[1]); /* creates pointer */
	Nacq = mxGetM(prhs[1]); /* get size */
	N = mxGetN(prhs[1]);
	
	/* Process third input: kf */
	//if( mxIsComplex(prhs[2]) ) { /* Warning for complez */
	//	mexWarnMsgTxt("gf5 WARNING: Imaginary part of kf (input 3) ignored."); }
	kf = mxGetPr(prhs[2]); /* creates pointer */
	i = mxGetM(prhs[2]); /* get size */
	j = mxGetN(prhs[2]);
	//if( (i != Nacq) || (j != N) ) {
	//	mexErrMsgTxt("gf5: kf (input 3) must have same dimensions as kx (input 2)."); }

	/* Process fourth input: M */
	//if( !mxIsComplex(prhs[3]) ) { /* Must be complez */
	//	mexErrMsgTxt("gf5: M (input 4) must be complex."); }
	Mre = mxGetPr(prhs[3]); /* creates pointers */
	Mim = mxGetPi(prhs[3]); 
	i = mxGetM(prhs[3]); /* get size */
	j = mxGetN(prhs[3]);
	//if( (i != Nacq) || (j != N) ) {
	//	mexErrMsgTxt("gf5: M (input 4) must have same dimensions as kx (input 2)."); }
	
	/* Process fifth input: fo */
	//if( mxIsComplex(prhs[4]) ) { /* Warning for complez */
	//	mexWarnMsgTxt("gf5 WARNING: Imaginary part of fo (input 5) ignored."); }
	fo = mxGetPr(prhs[4]); /* creates pointer */
	M = mxGetN(prhs[4]); /* get size */
	
	/* Process sixth input: xpos */
	//if( mxIsComplex(prhs[5]) ) { /* Warning for complez */
	//	mexWarnMsgTxt("gf5 WARNING: Imaginary part of xpos (input 6) ignored."); }
	xpos = mxGetPr(prhs[5]); /* creates pointer */
	i = mxGetM(prhs[5]); /* get size */
	j = mxGetN(prhs[5]);
	//if( (i != 1) || (j != N) ) {
	//	mexErrMsgTxt("gf5: xpos (input 6) must be a row vector with length N."); }

    /* create the output variables */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	if( nlhs == 2 ) { /* how do I know size? */
		plhs[1] = mxCreateDoubleMatrix(M*N*2+2*N,1,mxREAL); }

    /* get pointers to the output variables */
    f = mxGetPr(plhs[0]);
	if( nlhs == 2 ) {
		g = mxGetPr(plhs[1]); }
	else {
		g = NULL; }
	
    /* call the computational routine */
	gf5(x, kx, kf, Mre, Mim, fo, xpos, N, M, Nacq, f, g);

}

