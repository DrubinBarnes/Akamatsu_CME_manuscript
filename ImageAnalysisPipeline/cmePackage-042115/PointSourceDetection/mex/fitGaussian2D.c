/* [prmVect prmStd covarianceMatrix residuals Jacobian] = fitGaussian2Dmex(prmVect, initValues, mode);
 *
 * (c) Francois Aguet & Sylvain Berlemont, 2011 (last modified Feb 23, 2011)
 *
 * Compilation:
 * Mac/Linux: mex -I/usr/local/include -I../../mex/include /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a fitGaussian2D.c
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -I"..\..\..\extern\mex\include\gsl-1.15" -I"..\..\mex\include" "..\..\..\extern\mex\lib\gsl.lib" "..\..\..\extern\mex\lib\cblas.lib" -output fitGaussian2D fitGaussian2D.c
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

#define NPARAMS 5
#define refMode "xyasc"


typedef struct aStruct {
    double xi, yi, A, g, sigma2, sigma3;
} argStruct_t;

typedef int(*pfunc_t)(gsl_matrix*, int, int, argStruct_t*);

typedef struct dataStruct {
    int nx, np;
    double *pixels;
    double *gx, *gy;
    int *estIdx;
    int *idx;
    int nValid; /* number of non-NaN pixels */
    double *x_init;
    double prmVect[NPARAMS];
    pfunc_t *dfunc;
    gsl_vector *residuals;
    gsl_matrix *J;
    double maxIter, eAbs, eRel;
} dataStruct_t;



static int df_dx(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    gsl_matrix_set(J, i, k, A/s2*xi*g);
    return 0;
}

static int df_dy(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s2 = argStruct->sigma2;
    gsl_matrix_set(J, i, k, A/s2*yi*g);
    return 0;
}

static int df_dA(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    gsl_matrix_set(J, i, k, argStruct->g);
    return 0;
}

static int df_ds(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    double xi = argStruct->xi;
    double yi = argStruct->yi;
    double g = argStruct->g;
    double A = argStruct->A;
    double s3 = argStruct->sigma3;
    gsl_matrix_set(J, i, k, (xi*xi + yi*yi)*A/s3*g);
    return 0;
}

static int df_dc(gsl_matrix *J, int i, int k, argStruct_t *argStruct) {
    gsl_matrix_set(J, i, k, 1);
    return 0;
}



static int gaussian_f(const gsl_vector *x, void *params, gsl_vector *f) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int b = nx/2, i, k;
    
    double *pixels = dataStruct->pixels;
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    
    /* update prmVect with new estimates */
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double A = dataStruct->prmVect[2];
    double sigma = fabs(dataStruct->prmVect[3]);
    double c = dataStruct->prmVect[4];
    
    double xi, yi;
    double d = 2.0*sigma*sigma;
    for (i=0; i<nx; ++i) {
        k = i-b;
        xi = k-xp;
        yi = k-yp;
        gx[i] = exp(-xi*xi/d);
        gy[i] = exp(-yi*yi/d);
    }
    
    div_t divRes;
    int idx;
    for (i=0; i<dataStruct->nValid; ++i) {
        idx = dataStruct->idx[i];
        divRes = div(idx, nx);
        gsl_vector_set(f, i, A*gx[divRes.quot]*gy[divRes.rem]+c - pixels[idx]);
    }
    return GSL_SUCCESS;
}



static int gaussian_df(const gsl_vector *x, void *params, gsl_matrix *J) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int b = nx/2, i, k;
    
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    
    /* update prmVect with new estimates */
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double A = dataStruct->prmVect[2];
    double sigma = fabs(dataStruct->prmVect[3]);
    
    double xi, yi;
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    argStruct.A = A;
    
    for (i=0; i<nx; ++i) {
        k = i-b;
        xi = k-xp;
        yi = k-yp;
        gx[i] = exp(-xi*xi/d);
        gy[i] = exp(-yi*yi/d);
    }
    
    div_t divRes;
    int idx;
    for (i=0; i<dataStruct->nValid; ++i) {
        idx = dataStruct->idx[i];
        divRes = div(idx, nx);
        argStruct.xi = divRes.quot-b - xp;
        argStruct.yi = divRes.rem-b - yp;
        argStruct.g = gx[divRes.quot]*gy[divRes.rem];
        
        for (k=0; k<dataStruct->np; ++k)
            dataStruct->dfunc[k](J, i, k, &argStruct);
    }
    return GSL_SUCCESS;
}



static int gaussian_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
    
    dataStruct_t *dataStruct = (dataStruct_t *)params;
    int nx = dataStruct->nx;
    int b = nx/2, i, k;
    
    double *pixels = dataStruct->pixels;
    double *gx = dataStruct->gx;
    double *gy = dataStruct->gy;
    
    /* update prmVect with new estimates */
    for (i=0; i<dataStruct->np; ++i) {
        dataStruct->prmVect[dataStruct->estIdx[i]] = gsl_vector_get(x, i);
    }
    
    double xp = dataStruct->prmVect[0];
    double yp = dataStruct->prmVect[1];
    double A = dataStruct->prmVect[2];
    double sigma = fabs(dataStruct->prmVect[3]);
    double c = dataStruct->prmVect[4];
    
    double xi, yi;
    double sigma2 = sigma*sigma;
    double d = 2.0*sigma2;
    
    argStruct_t argStruct;
    argStruct.sigma2 = sigma2;
    argStruct.sigma3 = sigma2*sigma;
    argStruct.A = A;
    
    for (i=0; i<nx; ++i) {
        k = i-b;
        xi = k-xp;
        yi = k-yp;
        gx[i] = exp(-xi*xi/d);
        gy[i] = exp(-yi*yi/d);
    }
    
    div_t divRes;
    int idx;
    for (i=0; i<dataStruct->nValid; ++i) {
        idx = dataStruct->idx[i];
        divRes = div(idx, nx);
        
        argStruct.xi = divRes.quot-b - xp;
        argStruct.yi = divRes.rem-b - yp;
        argStruct.g = gx[divRes.quot]*gy[divRes.rem];
        gsl_vector_set(f, i, A*argStruct.g+c - pixels[idx]);
        
        for (k=0; k<dataStruct->np; ++k)
            dataStruct->dfunc[k](J, i, k, &argStruct);
    }
    return GSL_SUCCESS;
}



static int MLalgo(struct dataStruct *data) {
    
    /* declare solvers */
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    
    gsl_vector_view x = gsl_vector_view_array(data->x_init, data->np);
    
    const gsl_rng_type *type;
    
    gsl_rng_env_setup();
    
    type = gsl_rng_default;
    
    gsl_multifit_function_fdf f;
    f.f = &gaussian_f;
    f.df = &gaussian_df;
    f.fdf = &gaussian_fdf;
    f.n = data->nValid;
    f.p = data->np;
    f.params = data;
    
    
    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc(T, data->nValid, data->np);
    gsl_multifit_fdfsolver_set(s, &f, &x.vector);
    
    int status, status2;
    int iter = 0;
    gsl_vector *gradt = gsl_vector_alloc(data->np);
    
    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate(s);
        if (status)
            break;
        
        status = gsl_multifit_test_delta(s->dx, s->x, data->eAbs, data->eRel);
        gsl_multifit_gradient(s->J, s->f, gradt);
        status2 = gsl_multifit_test_gradient(gradt, data->eAbs);
    }
    while ((status == GSL_CONTINUE || status2 == GSL_CONTINUE) && iter < data->maxIter);
    
    gsl_vector_free(gradt);
    
    int i;
    for (i=0; i<data->np; ++i) {
        data->prmVect[data->estIdx[i]] = gsl_vector_get(s->x, i);
    }
    data->prmVect[3] = fabs(data->prmVect[3]);
    
    /* copy residuals */
    data->residuals = gsl_vector_alloc(data->nValid);
    gsl_vector_memcpy(data->residuals, s->f);
    
    /* copy Jacobian */
    data->J = gsl_matrix_alloc(data->nValid, data->np);
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
    size_t nx = mxGetN(prhs[0]);
    size_t ny = mxGetM(prhs[0]);
    if (nx != ny) mexErrMsgTxt("Input image must be a square.");
    if (!(nx % 2)) mexErrMsgTxt("The side of the input image must be odd.");
    int N = nx*ny;
    if (mxGetNumberOfElements(prhs[1])!=NPARAMS || !mxIsDouble(prhs[1])) mexErrMsgTxt("Incorrect parameter vector format.");
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
    int np = (int)mxGetNumberOfElements(prhs[2]);
    char *mode;
    mode = (char*)malloc(sizeof(char)*(np+1));
    mxGetString(prhs[2], mode, np+1);
    
    int i;
    for (i=0; i<strlen(mode); ++i) {
        mode[i] = tolower(mode[i]);
    }
    
    np = 0; /* number of parameters to fit */
    for (i=0; i<NPARAMS; ++i) {
        if (strchr(mode, refMode[i])!=NULL) { np++; }
    }
    if (np==0) mexErrMsgTxt("Unknown mode.");
    
    /* allocate */
    data.nx = nx;
    data.np = np;
    data.pixels = mxGetPr(prhs[0]);
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*nx);
    data.estIdx = (int*)malloc(sizeof(int)*np);
    memcpy(data.prmVect, mxGetPr(prhs[1]), NPARAMS*sizeof(double));
    data.dfunc = (pfunc_t*) malloc(sizeof(pfunc_t) * np);
    
    /* read mask/pixels */
    data.nValid = N;
    for (i=0; i<N; ++i) {
        data.nValid -= (int)mxIsNaN(data.pixels[i]);
    }
    if (data.nValid < 5) mexErrMsgTxt("Input image must contain at least 5 data points.");

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
    
    np = 0;
    if (strchr(mode, 'x')!=NULL) {data.estIdx[np] = 0; data.dfunc[np++] = df_dx;}
    if (strchr(mode, 'y')!=NULL) {data.estIdx[np] = 1; data.dfunc[np++] = df_dy;}
    if (strchr(mode, 'a')!=NULL) {data.estIdx[np] = 2; data.dfunc[np++] = df_dA;}
    if (strchr(mode, 's')!=NULL) {data.estIdx[np] = 3; data.dfunc[np++] = df_ds;}
    if (strchr(mode, 'c')!=NULL) {data.estIdx[np] = 4; data.dfunc[np++] = df_dc;}
    
    data.x_init = (double*)malloc(sizeof(double)*np);
    for (i=0; i<np; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }
    
    MLalgo(&data);
    
    /* parameters */
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(1, NPARAMS, mxREAL);
        memcpy(mxGetPr(plhs[0]), data.prmVect, NPARAMS*sizeof(double));
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
        gsl_matrix *covar = gsl_matrix_alloc(np, np);
        gsl_multifit_covar(data.J, 0.0, covar);
        double iRSS = RSS/(data.nValid - data.np - 1);
        plhs[1] = mxCreateDoubleMatrix(1, data.np, mxREAL);
        double *prmStd = mxGetPr(plhs[1]);
        for (i=0; i<data.np; ++i) {
            prmStd[i] = sqrt(iRSS*gsl_matrix_get(covar, i, i));
        }
        if (nlhs > 2) {
            plhs[2] = mxCreateDoubleMatrix(np, np, mxREAL);
            /* cov. matrix is symmetric, no need to transpose */
            memcpy(mxGetPr(plhs[2]), covar->data, np*np*sizeof(double));
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
        mxSetFieldByNumber(plhs[3], 0, 1, mxCreateLogicalScalar(hAD));
        mxSetFieldByNumber(plhs[3], 0, 2, mxCreateDoubleScalar(mean));
        mxSetFieldByNumber(plhs[3], 0, 3, mxCreateDoubleScalar(std));
        mxSetFieldByNumber(plhs[3], 0, 4, mxCreateDoubleScalar(RSS));
    }
    
    /* Jacobian */
    if (nlhs > 4) {
        /* convert row-major double* data.J->data to column-major double* */
        plhs[4] = mxCreateDoubleMatrix(N, np, mxREAL);
        double *J = mxGetPr(plhs[4]);
        int k;
        for (k=0; k<np; ++k) {
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
    free(data.x_init);
    free(nanIdx);
    free(data.idx);
    free(data.dfunc);
    free(data.estIdx);
    free(data.gy);
    free(data.gx);
    free(mode);
}





/* compile with:
 * export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2010b.app/bin/maci64 && gcc -Wall -g -DARRAY_ACCESS_INLINING -I. -I/Applications/MATLAB_R2010b.app/extern/include -L/Applications/MATLAB_R2010b.app/bin/maci64 -lmx -lmex -lgsl -lgslcblas -lmat fitGaussian2D.c
*/

/*int main(void) {
    
    int nx = 15;
    int N = nx*nx;
    double* px;
    px = (double*)malloc(sizeof(double)*N);
    
    int i;
    // fill with noise
    for (i=0; i<N; ++i) {
        px[i] = rand();
    }
    
    int np = 5;
    
    dataStruct_t data;
    data.nx = nx;
    data.np = 5;
    data.pixels = px;
    data.gx = (double*)malloc(sizeof(double)*nx);
    data.gy = (double*)malloc(sizeof(double)*nx);
    
    data.estIdx = (int*)malloc(sizeof(int)*np);
    data.dfunc = (pfunc_t*) malloc(sizeof(pfunc_t) * np);
    
    // read mask/pixels
    data.nValid = N;
    data.idx = (int*)malloc(sizeof(int)*data.nValid);
    int k = 0;
    for (i=0; i<N; ++i) {
        if (!mxIsNaN(data.pixels[i])) {
            data.idx[k++] = i;
        }
    }
    
    data.prmVect[0] = 0;
    data.prmVect[1] = 0;
    data.prmVect[2] = 5;
    data.prmVect[3] = 1;
    data.prmVect[4] = 0;
    
    np = 0;
    data.estIdx[np] = 0; data.dfunc[np++] = df_dx;
    data.estIdx[np] = 1; data.dfunc[np++] = df_dy;
    data.estIdx[np] = 2; data.dfunc[np++] = df_dA;
    data.estIdx[np] = 3; data.dfunc[np++] = df_ds;
    data.estIdx[np] = 4; data.dfunc[np++] = df_dc;
    
    data.x_init = (double*)malloc(sizeof(double)*np);
    for (i=0; i<np; ++i) {
        data.x_init[i] = data.prmVect[data.estIdx[i]];
    }
    
    MLalgo(&data);
    
    gsl_matrix_free(data.J);
    gsl_vector_free(data.residuals);
    free(data.x_init);
    free(data.idx);
    free(data.dfunc);
    free(data.estIdx);
    free(data.gy);
    free(data.gx);
    free(px);
    
    return 0;
}*/
