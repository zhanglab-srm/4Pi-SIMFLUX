#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "memory.h"

#include <levmar.h>
// extern int dlevmar_der(
//       void (*func)(double *p, double *hx, int m, int n, void *adata),
//       void (*jacf)(double *p, double *j, int m, int n, void *adata),
//       double *p, double *x, int m, int n, int itmax, double *opts,
//       double *info, double *work, double *covar, void *adata);
// 
// extern int dlevmar_dif(
//       void (*func)(double *p, double *hx, int m, int n, void *adata),
//       double *p, double *x, int m, int n, int itmax, double *opts,
//       double *info, double *work, double *covar, void *adata);

#include <mex.h>

// 2D gaussian
//f = p[4]*exp(-(x-p[0])^2/p[2]-(y-p[1])^2/p[3])+p[5]
//adata: int * for width
void func(double *p, double *hx, int m, int n, void *adata)
{
    int imgwidth = *((int *)adata);
    int i;
    double cx, cy;
    for(i=0; i<n; i++){
        cy = i%imgwidth;
        cx = i/imgwidth;
        hx[i] = p[4]*exp(-(cx-p[0])*(cx-p[0])/(2.0*p[2]*p[2])-(cy-p[1])*(cy-p[1])/(2.0*p[3]*p[3]))+p[5];
    }
}

void jacf(double *p, double *j, int m, int n, void *adata)
{
    int imgwidth = *((int *)adata);
    int i, jcnt;
    double cx, cy, temp, temp2;
    for(i=jcnt=0; i<n; i++){
        cy = i%imgwidth;
        cx = i/imgwidth;
        temp = exp(-(cx-p[0])*(cx-p[0])/(2.0*p[2]*p[2])-(cy-p[1])*(cy-p[1])/(2.0*p[3]*p[3]));// EXP
        temp2 = p[4] * temp; //int*EXP
        
        j[jcnt++] = temp2*(cx-p[0])/(p[2]*p[2]); //p[0]
        j[jcnt++] = temp2*(cy-p[1])/(p[3]*p[3]); //p[1]
        j[jcnt++] = temp2*(cx-p[0])*(cx-p[0])/(p[2]*p[2]*p[2]); //p[2]
        j[jcnt++] = temp2*(cy-p[1])*(cy-p[1])/(p[3]*p[3]*p[3]); //p[3]
        j[jcnt++] = temp; //p[4]
        j[jcnt++] = 1.0; //p[5]
    }
}

//init parameter, cx = cy = center point,
//stdx = stdy = initstd
//height = max - min, bkg = min
void initpar(double *p, int width, double initstd, double *lp)
{
    double cmin = p[0], cmax = p[0];
    int m, plen = width * width;
    int maxidx = 0;
    for(m=1;m<plen; m++){
        cmin = (cmin>p[m]) ? p[m]: cmin; 
        cmax = (cmax<p[m]) ? p[m]: cmax; 
        maxidx = (cmax==p[m]) ? m: maxidx; 
    }
    lp[0] = maxidx % width;
    lp[1] = maxidx / width;
    lp[3] = lp[2] = initstd;
    lp[4] = cmax-cmin;
    lp[5] = cmin;
}

/*MEX-function
 *[lplist, infolist] = GausFit2D_CPU(imgstack, opts, initstd)
 *
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *p = NULL;
    double lp[6] = {0.0};
    int imgwidth = mxGetM(prhs[0]);
    double *opts;
    double info[10];
    double initstd;
    double *pimg, *plp, *pinfo;
    int m,n;
    int imglen = 1;
    int ret;
    
    //get number of images
    if(mxGetNumberOfDimensions(prhs[0]) ==3){
        imglen = mxGetDimensions(prhs[0])[2];
    }
    //printf("imglen = %d\n", imglen);
    
    //default opts: [1e-3, 1e-15, 1e-15, 1e-15, 1e-15];
//     opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
//     opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used
    opts = mxGetPr(prhs[1]);
    initstd = mxGetScalar(prhs[2]);
    //printf("opts: %E, %E, %E, %E, %E\n",opts[0], opts[1], opts[2], opts[3], opts[3]);
    
    //init result
    plhs[0] = mxCreateDoubleMatrix(imglen,6,false);
    plhs[1] = mxCreateDoubleMatrix(imglen,10,false);
    
    //fit gaussian
    pimg = mxGetPr(prhs[0]);
    plp = mxGetPr(plhs[0]);
    pinfo = mxGetPr(plhs[1]);
    
    for(m=0;m<imglen;m++){
        initpar(pimg, imgwidth, initstd, lp);
        //dlevmar_dif(func,lp,mxGetPr(prhs[0]),6,7*7,200,opts,info,NULL,NULL,(void *)(&imgwidth));
        ret = dlevmar_der(func, jacf, lp, pimg, 6, imgwidth*imgwidth, 
                (int)(opts[5]), opts, info, NULL,NULL,(void *)(&imgwidth));
        
        //copy results
        for(n=0;n<6;n++){
            plp[n*imglen] = lp[n];
        }
        for(n=0;n<10;n++){
            pinfo[n*imglen] = info[n];
        }
//         memcpy(plp, lp, 6*sizeof(double));
//         memcpy(pinfo, info, 10*sizeof(double));
        
        //update pointers
        pimg += imgwidth*imgwidth;
        plp ++;
        pinfo ++;
    }
}
