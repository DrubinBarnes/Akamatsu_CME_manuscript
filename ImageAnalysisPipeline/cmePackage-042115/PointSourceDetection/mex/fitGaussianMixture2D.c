/* [prmVect prmStd covarianceMatrix residuals Jacobian] = fitGaussianMixture2D(image, prmVect, mode);
 *
 * (c) Francois Aguet, 2011 (last modified 06/26/2012)
 *
 * Compilation:
 * Mac/Linux: mex -I/usr/local/include -lgsl -lgslcblas fitGaussianMixture2D.c
 * Mac/Linux (static): mex -I/usr/local/include -I../../mex/include /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a fitGaussianMixture2D.c
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\..\extern\mex\include\gsl-1.15" -I"..\..\mex\include" "..\..\..\extern\mex\lib\gsl.lib" "..\..\..\extern\mex\lib\cblas.lib" -output fitGaussianMixture2D fitGaussianMixture2D.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "mex.h"
#include "matrix.h"
#include "stats.h"


#define refMode "xyasc"


typedef struct argStruct {
    double xi, yi, A, g, sigma2, sigma3;
} argStruct_t;

typedef int(*pfunc_t)(double*, int, argStruct_t*);

typedef struct dataStruct {
    int nx;               // width/height of input
    int np;               // # input parameters
    int ng;               // # gaussians in mixture
    int nparam;           // # of parameters to optimize
    int step;             // increment in parameter vector; i.e., 3 if 'xyA', 2 if 'xy' etc.

    double *pixels;       // input array
    double *buffer;       // buffer array for calculations
    double *gx, *gy;      // 1-D separated components of a Gaussian: 
			              // exp(-(x-x0)^2/(2*sigma^2))
    
    int *estIdx;          // indexes of prmVect that will be optimized, size: nparam
    int *idx;             // index of non-NaN pixels
    int nValid;           // number of non-NaN pixels
    double *x_init;       // initial values for optimization
    double *prmVect;      // parameter vector: 3*ng+2: x1, y1, A1, ... xn, yn, An, sigma, background
    
    pfunc_t *dfunc;       // function pointer for derivatives
    gsl_vector *residuals;
    gsl_matrix *J;        // Jacobian matrix
    double *Jbuffer;      // required, since jacobian is sometimes iteratively calculated (can't use gsl_matrix_set)
    double maxIter, eAbs, eRel; // optimiser settings, see GSL doc.
} dataStruct_t;


// Partial derivatives of the cost function
static int df_dx(double *J, int i, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    J[i] = A/s2*xi*g;
    return 0;
}

static int df_dy(double *J, int i, argStruct_t *argStruct) {
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    J[i] = A/s2*yi*g;
    return 0;
}

static int df_dA(double *J, int i, argStruct_t *argStruct) {
    J[i] = argStruct->g;
    return 0;
}

static int df_ds(double *J, int i, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s3 = argStruct->sigma3;
    J[i] += (xi*xi + yi*yi)*A/s3*g;
    return 0;
}

static int df_dc(double *J, int i, argStruct_t *argStruct) {
    J[i] = 1;
    return 0;
}



int gaussian_f(const gsl_vector *x, void *params, gsl_vector *f) {
    
    dataStruct_t *data = (dataStruct_t *)params;

    int nx = data->nx;
    int b = nx/2, i, k;
    int np = data->np;
    
    double *buffer = data->buffer;
    double *gx = data->gx;
    double *gy = data->gy;

    // update prmVect with new estimates
    for (i=0; i<data->nparam; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp, yp, A, xi, yi;
    double sigma = fabs(data->prmVect[np-2]);
    double d = 2.0*sigma*sigma;
    double c = data->prmVect[np-1];
    
    int gi, idx;
    div_t divRes;
    
    // loop through non-NaN pixels; background component of cost function
    for (i=0; i<data->nValid; ++i) {
        idx = data->idx[i];
        divRes = div(idx, nx);
        buffer[i] = c - data->pixels[idx];
    }
    
    // add individual gaussians to cost
    for (gi=0; gi<data->ng; ++gi) {
        xp = data->prmVect[3*gi];
        yp = data->prmVect[3*gi+1];
        A  = fabs(data->prmVect[3*gi+2]);
        
        // gaussian kernel
        for (i=0; i<nx; ++i) {
            k = i-b;
            xi = k-xp;
            yi = k-yp;
            gx[i] = exp(-xi*xi/d);
            gy[i] = exp(-yi*yi/d);
        }
        
        for (i=0; i<data->nValid; ++i) {
            idx = data->idx[i];
            divRes = div(idx, nx);
            buffer[i] += A*gx[divRes.quot]*gy[divRes.rem];
        }
    }
    
    for (i=0; i<data->nValid; ++i) {
        gsl_vector_set(f, i, buffer[i]);
    }
    return GSL_SUCCESS;
}


// partial derivatives only
int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J) {
    
    dataStruct_t *data = (dataStruct_t *)params;

    // initialize Jacobian
    int nJ = data->nparam * data->nValid;
    memset(data->Jbuffer, 0, nJ*sizeof(double));
    
    int nx = data->nx;
    int b = nx/2, i, k;
    
    double *gx = data->gx;
    double *gy = data->gy;
    
    /* update prmVect with new estimates */
    for (i=0; i<data->nparam; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp, yp, xi, yi;
    double sigma = fabs(data->prmVect[data->np-2]);
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    
    int idx;
    div_t divRes;
    
    // loop through {x,y,A} groups
    int gi;
    for (gi=0;gi<data->ng;++gi) {
        xp = data->prmVect[3*gi];
        yp = data->prmVect[3*gi+1];
        argStruct.A = fabs(data->prmVect[3*gi+2]);
        
        // x and y components of the gaussian kernel
        for (i=0; i<nx; ++i) {
            k = i-b;
            xi = k-xp;
            yi = k-yp;
            gx[i] = exp(-xi*xi/d);
            gy[i] = exp(-yi*yi/d);
        }
        
        for (i=0; i<data->nValid; ++i) {
            idx = data->idx[i];
            divRes = div(idx, nx);
            argStruct.xi = divRes.quot-b - xp;
            argStruct.yi = divRes.rem-b - yp;
            argStruct.g = gx[divRes.quot]*gy[divRes.rem];
            
            // df/dx, df/dy, df/dA
            for (k=0; k<data->step; ++k) {
                data->dfunc[gi*data->step + k](data->Jbuffer, i+(gi*data->step+k)*data->nValid, &argStruct);
            }
            
            // df/ds, df/dc
            for (k=data->ng*data->step; k<data->nparam; ++k) {
                data->dfunc[k](data->Jbuffer, i+k*data->nValid, &argStruct);
            }
        }
    }
    
    /* copy Jacobian */
    for (i=0; i<nJ; ++i) {
        divRes = div(i, data->nValid);
        gsl_matrix_set(J, divRes.rem, divRes.quot, data->Jbuffer[i]);
    }
    return GSL_SUCCESS;
}



int gaussian_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
    
    dataStruct_t *data = (dataStruct_t *)params;
    /* initialize Jacobian */
    int nJ = data->nparam * data->nValid;
    memset(data->Jbuffer, 0, nJ*sizeof(double));
    
    int nx = data->nx;
    int b = nx/2, i, k;
    
    double *buffer = data->buffer;
    double *gx = data->gx;
    double *gy = data->gy;
    
    /* update prmVect with new estimates */
    for (i=0; i<data->nparam; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp, yp, xi, yi, A;
    double sigma = fabs(data->prmVect[data->np-2]);
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    double c = data->prmVect[data->np-1];
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    
    int idx;
    div_t divRes;
    
    for (i=0; i<data->nValid; ++i) {
        idx = data->idx[i];
        divRes = div(idx, nx);
        buffer[i] = c - data->pixels[idx];
    }
    
    int gi;
    for (gi=0;gi<data->ng;++gi) {
        xp = data->prmVect[3*gi];
        yp = data->prmVect[3*gi+1];
        A  = fabs(data->prmVect[3*gi+2]);
        argStruct.A = A;
        
        for (i=0; i<nx; ++i) {
            k = i-b;
            xi = k-xp;
            yi = k-yp;
            gx[i] = exp(-xi*xi/d);
            gy[i] = exp(-yi*yi/d);
        }
        
        for (i=0; i<data->nValid; ++i) {
            idx = data->idx[i];
            divRes = div(idx, nx);
            
            buffer[i] += A*gx[divRes.quot]*gy[divRes.rem];
            
            argStruct.xi = divRes.quot-b - xp;
            argStruct.yi = divRes.rem-b - yp;
            argStruct.g = gx[divRes.quot]*gy[divRes.rem];
            
            // df/dx, df/dy, df/dA
            for (k=0; k<data->step; ++k) {
                data->dfunc[gi*data->step + k](data->Jbuffer, i+(gi*data->step+k)*data->nValid, &argStruct);
            }
            
            // df/ds, df/dc
            for (k=data->ng*data->step; k<data->nparam; ++k) {
                data->dfunc[k](data->Jbuffer, i+k*data->nValid, &argStruct);
            }
        }
    }
    
    for (i=0; i<data->nValid; ++i) {
        gsl_vector_set(f, i, buffer[i]);
    }
    
    // copy Jacobian; nJ = nValid*nparam
    for (i=0; i<nJ; ++i) {
        divRes = div(i, data->nValid);
        gsl_matrix_set(J, divRes.rem, divRes.quot, data->Jbuffer[i]);
    }
    return GSL_SUCCESS;
}



int MLalgo(struct dataStruct *data) {
    
    /* declare solvers */
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
       
    gsl_vector_view x = gsl_vector_view_array(data->x_init, data->nparam);
    
    const gsl_rng_type *type;
    
    gsl_rng_env_setup();
    
    type = gsl_rng_default;
    
    gsl_multifit_function_fdf f;
    f.f = &gaussian_f;
    f.df = &gaussian_df;
    f.fdf = &gaussian_fdf;
    f.n = data->nValid;
    f.p = data->nparam;
    f.params = data;
    
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, data->nValid, data->nparam);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);
    
    int status, status2;
    int iter = 0;
    gsl_vector *gradt = gsl_vector_alloc(data->nparam);
    
    do {
        iter++;
        
        status = gsl_multifit_fdfsolver_iterate(s);
        if (status)
            break;
        
        status = gsl_multifit_test_delta(s->dx, s->x, data->eAbs, data->eRel);
        gsl_multifit_gradient(s->J, s->f, gradt);
        status2 = gsl_multifit_test_gradient(gradt, data->eAbs);
        
        
    }
    while (status == GSL_CONTINUE && status2 == GSL_CONTINUE && iter < data->maxIter);
    gsl_vector_free(gradt);
    
    int i;
    for (i=0; i<data->nparam; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(s->x, i);
    }
    
    /* copy residuals */
    data->residuals = gsl_vector_alloc(data->nValid);
    gsl_vector_memcpy(data->residuals, s->f);
    
    /* copy Jacobian */
    data->J = gsl_matrix_alloc(data->nValid, data->nparam);
    gsl_matrix_memcpy(data->J, s->J);
    
    gsl_multifit_fdfsolver_free(s);
    return 0;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* inputs:
     * image array
     * prmVect
     * mode
     * {options}
     */
    
    dataStruct_t data;
    
    /* check inputs */
    if (nrhs < 3) mexErrMsgTxt("Inputs should be: data, prmVect, mode.");
    if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Data input must be double array.");
    int nx = (int)mxGetN(prhs[0]);
    int ny = (int)mxGetM(prhs[0]);
    if (nx != ny) mexErrMsgTxt("Input should be a square image.");
    int N = nx*ny;
    if (!mxIsDouble(prhs[1])) mexErrMsgTxt("Parameter vector must be double array.");
    int np = (int)mxGetNumberOfElements(prhs[1]);
    div_t divRes = div(np-2, 3);
    int ng = divRes.quot;
    if (divRes.rem != 0) mexErrMsgTxt("Invalid parameter vector length.");
    if (!mxIsChar(prhs[2])) mexErrMsgTxt("Mode needs to be a string.");
    if (nrhs < 4) {
        data.maxIter = 500;
        data.eAbs = 1e-8;
        data.eRel = 1e-8;
    } else {
        if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3])!=3) mexErrMsgTxt("Options must must be double array with 3 elements.");
        double *options = mxGetPr(prhs[3]);
        data.maxIter = options[0];
        data.eAbs = options[1];
        data.eRel = options[2];
    }

    
    /* read mode input */
    size_t nm = mxGetNumberOfElements(prhs[2]);
    char *mode;
    mode = (char*)malloc(sizeof(char)*(nm+1));
    mxGetString(prhs[2], mode, nm+1);
    
    int i;
    for (i=0; i<strlen(mode); ++i) {
        mode[i] = tolower(mode[i]);
    }
    
    /* detect parameters to optimize */
    int step = 0;
    for (i=0; i<3; ++i) {
        if (strchr(mode, refMode[i])!=NULL) {
            ++step;
        }
    }
    int nparam = step*ng;
    
    for (i=3; i<5; ++i) {
        if (strchr(mode, refMode[i])!=NULL) {
            ++nparam;
        }
    }
    data.nparam = nparam;
    data.step = step;
        
    /* allocate */
    data.nx = nx;
    data.np = np;
    data.ng = ng;
    data.pixels = mxGetPr(prhs[0]);
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*nx);
    data.estIdx = (int*)malloc(sizeof(int)*nparam);
    data.prmVect = mxGetPr(prhs[1]);
    data.dfunc = (pfunc_t*)malloc(sizeof(pfunc_t)*nparam);
    
    
    /* read mask/pixels */
    data.nValid = N;
    for (i=0; i<N; ++i) {
        data.nValid -= (int)mxIsNaN(data.pixels[i]);
    }
    data.buffer = (double*)malloc(sizeof(double)*data.nValid);
    data.idx = (int*)malloc(sizeof(int)*data.nValid);
    int *nanIdx = (int*)malloc(sizeof(int)*(N-data.nValid));
    int k = 0, l = 0;
    for (i=0; i<N; ++i) {
        if (!mxIsNaN(data.pixels[i])) {
            data.idx[k++] = i;
        } else {
            nanIdx[l++] = i;
        }
    }
    
    // np was the # elements in prmVect; now tmp variable
    np = 0;
    for (i=0; i<ng; ++i) {
        if (strchr(mode, 'x')!=NULL) {
            data.estIdx[np] = 3*i;
            data.dfunc[np++] = df_dx;
        }
        if (strchr(mode, 'y')!=NULL) {
            data.estIdx[np] = 1+3*i;
            data.dfunc[np++] = df_dy;
        }
        if (strchr(mode, 'a')!=NULL) {
            data.estIdx[np] = 2+3*i;
            data.dfunc[np++] = df_dA;
        }
    }
    int pos = step*ng;
    if (strchr(mode, 's')!=NULL) { data.estIdx[pos++] = 3*ng; data.dfunc[np++] = df_ds; }
    if (strchr(mode, 'c')!=NULL) { data.estIdx[pos] = 3*ng+1; data.dfunc[np++] = df_dc; }
            
    data.x_init = (double*)malloc(sizeof(double)*nparam);
    for (i=0; i<nparam; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }
    data.Jbuffer = (double*)malloc(sizeof(double)*nparam*data.nValid);
    
    MLalgo(&data);
    
    /* parameters */
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(1, data.np, mxREAL);
        memcpy(mxGetPr(plhs[0]), data.prmVect, data.np*sizeof(double));
    }
    
    /* standard dev. of parameters & covariance matrix */
    double RSS = 0.0;
    double* resValid = NULL;
    if (nlhs > 1) {
        resValid = (double*)malloc(data.nValid*sizeof(double));
        for (i=0; i<data.nValid; ++i) {
            resValid[i] = gsl_vector_get(data.residuals, i);
            RSS += resValid[i]*resValid[i];
        }
        gsl_matrix *covar = gsl_matrix_alloc(nparam, nparam);
        
        gsl_multifit_covar(data.J, 0.0, covar);
        double iRSS = RSS/(data.nValid - nparam - 1);
        plhs[1] = mxCreateDoubleMatrix(1, data.np, mxREAL); // expand
        double *prmStd = mxGetPr(plhs[1]);
        memset(prmStd, 0, data.np*sizeof(double)); // zero Var for fixed parameters
        for (i=0; i<nparam; ++i) {
            prmStd[data.estIdx[i]] = sqrt(iRSS*gsl_matrix_get(covar, i, i));
        }
        if (nlhs > 2) {
            plhs[2] = mxCreateDoubleMatrix(nparam, nparam, mxREAL);
            /* cov. matrix is symmetric, no need to transpose */
            memcpy(mxGetPr(plhs[2]), covar->data, nparam*nparam*sizeof(double));
        }
        gsl_matrix_free(covar);
    }
    
    /* residuals */
    if (nlhs > 3) {
        const char *fieldnames[] = {"data", "hAD", "mean", "std", "RSS"};
        mwSize dims[2] = {1, 1};
        plhs[3] = mxCreateStructArray(2, dims, 5, fieldnames);        
        mxArray *val = mxCreateDoubleMatrix(nx, nx, mxREAL);
        double* res = mxGetPr(val);
        
        double mean = 0.0, std = 0.0;
        for (i=0; i<data.nValid; ++i) {
            res[data.idx[i]] = resValid[i];
            mean += resValid[i];
        }
        std = sqrt((RSS-mean*mean/data.nValid)/(data.nValid-1));
        mean /= data.nValid;
       
        for (i=0; i<N-data.nValid; ++i) {
            res[nanIdx[i]] = mxGetNaN();
        }
        
        // A-D test, case 2: mean known
        unsigned char hAD = adtest(resValid, data.nValid, 2, 0.0, std, 0.05);
        mxSetFieldByNumber(plhs[3], 0, 0, val);   
        mxSetFieldByNumber(plhs[3], 0, 1, mxCreateDoubleScalar(hAD));
        mxSetFieldByNumber(plhs[3], 0, 2, mxCreateDoubleScalar(mean));
        mxSetFieldByNumber(plhs[3], 0, 3, mxCreateDoubleScalar(std));
        mxSetFieldByNumber(plhs[3], 0, 4, mxCreateDoubleScalar(RSS));
    }
    
    // Jacobian
    if (nlhs > 4) {
        // convert row-major double* data.J->data to column-major double*
        // expand Jacobian from (nValid x nparam) to (N x nparam)
        plhs[4] = mxCreateDoubleMatrix(N, nparam, mxREAL);
        double *J = mxGetPr(plhs[4]);
        int k;
        for (k=0; k<nparam; ++k) {
            for (i=0; i<data.nValid; ++i) {
                J[data.idx[i]+k*N] = gsl_matrix_get(data.J, i, k);
            }
            for (i=0; i<N-data.nValid; ++i) {
                J[nanIdx[i]+k*N] = mxGetNaN();
            }
        }
    }
    
    free(resValid);
    gsl_matrix_free(data.J);
    gsl_vector_free(data.residuals);
    free(data.Jbuffer);
    free(data.x_init);
    free(nanIdx);
    free(data.idx);
    free(data.buffer);
    free(data.dfunc);
    free(data.estIdx);
    free(data.gy);
    free(data.gx);
    free(mode);
}

// compiled with
// export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2012a.app/bin/maci64 && gcc -Wall -g -DARRAY_ACCESS_INLINING -I. -I../../mex/include/ -I/Applications/MATLAB_R2012a.app/extern/include -L/Applications/MATLAB_R2012a.app/bin/maci64 -lmx -lmex -lgsl -lgslcblas -lmat fitGaussianMixture2D.c
// tested with:
// valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./a.out 2>&1 | grep Mixture

/* int main(void) {
    
    dataStruct_t data;

    // skip inputs, define
    
    int nx = 15;
    int N = nx*nx;
    data.maxIter = 500;
    data.eAbs = 1e-8;
    data.eRel = 1e-8;
    
    double* px;
    px = (double*)malloc(sizeof(double)*N);
    
    int i;
    // fill with noise
    for (i=0; i<N; ++i) {
        px[i] = rand();
    }
     
    int nparam = 8;
    data.nparam = nparam;
    data.step = 3;
    
    // allocate
    data.nx = nx;
    data.np = 8;
    data.ng = 2;
    data.pixels = px;
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*nx);
    data.estIdx = (int*)malloc(sizeof(int)*nparam);
    
    data.prmVect = (double*)malloc(sizeof(double)*8);
    data.prmVect[0] = 0;
    data.prmVect[1] = 0;
    data.prmVect[2] = 5;
    data.prmVect[3] = 0;
    data.prmVect[4] = 0;
    data.prmVect[5] = 5;
    data.prmVect[6] = 1;
    data.prmVect[7] = 0;
    
    data.dfunc = (pfunc_t*)malloc(sizeof(pfunc_t)*nparam);
    
    
    // read mask/pixels
    data.nValid = N;
    for (i=0; i<N; ++i) {
        data.nValid -= (int)mxIsNaN(data.pixels[i]);
    }
    data.buffer = (double*)malloc(sizeof(double)*data.nValid);
    data.idx = (int*)malloc(sizeof(int)*data.nValid);
    int *nanIdx = (int*)malloc(sizeof(int)*(N-data.nValid));
    int k = 0, l = 0;
    for (i=0; i<N; ++i) {
        if (!mxIsNaN(data.pixels[i])) {
            data.idx[k++] = i;
        } else {
            nanIdx[l++] = i;
        }
    }
    
    for (i=0; i<8; ++i) {
        data.estIdx[i] = i;
    }
    data.dfunc[0] = df_dx;
    data.dfunc[1] = df_dy;
    data.dfunc[2] = df_dA;
    data.dfunc[3] = df_dx;
    data.dfunc[4] = df_dy;
    data.dfunc[5] = df_dA;
    data.dfunc[6] = df_ds;
    data.dfunc[7] = df_dc;
    
    data.x_init = (double*)malloc(sizeof(double)*nparam);
    for (i=0; i<nparam; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }
    data.Jbuffer = (double*)malloc(sizeof(double)*nparam*data.nValid);
    
    MLalgo(&data);
    
    gsl_matrix_free(data.J);
    gsl_vector_free(data.residuals);
    
    free(data.Jbuffer);
    free(data.x_init);
    free(nanIdx);
    free(data.idx);
    free(data.buffer);
    free(data.dfunc);
    free(data.estIdx);
    free(data.gy);
    free(data.gx);
    
    free(data.prmVect);
    free(px);
    return 0;
} */
