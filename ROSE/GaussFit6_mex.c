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

// adata struct
typedef struct{
    int width;
    double x0;
    double y0;
    double sx;
    double sy;
} ADATA;
// 2D gaussian
//f = p[0]*exp(-0.5*(x-x0)^2/sx^2-0.5*(y-y0)^2/sy^2)+p[1]
//lp:[int1-6 bkg1-6]
//adata: int * for width
void func(double *p, double *hx, int m, int n, void *adata)
{
    ADATA *padata = (ADATA *)adata;
    int imgwidth = padata->width;
    double x0 = padata->x0;
    double y0 = padata->y0;
    double sx = padata->sx;
    double sy = padata->sy;
    int i, imgsize = n;
    double cx, cy, cexp;
    
    //img1
    for(i=0; i<imgsize; i++){
        cy = i%imgwidth;
        cx = i/imgwidth;
        cexp = exp(-(cx-x0)*(cx-x0)/(2.0*sx*sx)-(cy-y0)*(cy-y0)/(2.0*sy*sy));
        //img1
        hx[i] = p[0]*cexp+p[1];
//         hx[i+imgsize*1] = p[1]*cexp+p[7];
//         hx[i+imgsize*2] = p[2]*cexp+p[8];
//         hx[i+imgsize*3] = p[3]*cexp+p[9];
//         hx[i+imgsize*4] = p[4]*cexp+p[10];
//         hx[i+imgsize*5] = p[5]*cexp+p[11];
    }
}

void jacf(double *p, double *j, int m, int n, void *adata)
{
    ADATA *padata = (ADATA *)adata;
    int imgwidth = padata->width;
    double x0 = padata->x0;
    double y0 = padata->y0;
    double sx = padata->sx;
    double sy = padata->sy;
    int i, jcnt = 0, imgsize = n;
    double cx, cy, cexp;
    
//     for(i=0; i<imgsize*m; i++){
//         j[i] = 0.0;
//     }
    //img1
    jcnt = 0;
    for(i=0; i<imgsize; i++){
        cy = i%imgwidth;
        cx = i/imgwidth;
        cexp = exp(-(cx-x0)*(cx-x0)/(2.0*sx*sx)-(cy-y0)*(cy-y0)/(2.0*sy*sy));
        //img1
        j[jcnt++] = cexp;
        j[jcnt++] = 1;
    }
}
//init parameter, cx = cy = center point,
//stdx = stdy = initstd
//height = max - min, bkg = min
void initpar(double *p, int width, double *lp)
{
    double cmin = p[0], cmax = p[0];
    int m, plen = width * width;
    for(m=1;m<plen; m++){
        cmin = (cmin>p[m]) ? p[m]: cmin; 
        cmax = (cmax<p[m]) ? p[m]: cmax; 
    }
    lp[0] = cmax-cmin;
    lp[1] = cmin;
}
// void initpar(double *p, int width, double *lp)
// {
//     int m;
//     for(m=0;m<6;m++){
//         initpar_lite(p+m*width*width, width, lp+m);
//     }
// }



/*MEX-function
 *[lplist, infolist] = GausFit6_mex(imgstack, opts, lplist)
 *lplist: [x0, y0, sx, sy]
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *p = NULL;
    double lp[2] = {0.0};
    int imgwidth = mxGetM(prhs[0]);
    double *opts;
    double info[10];
    double *plpin;
    double *pimg, *plp, *pinfo;
    int m,n,subcnt;
    int imglen = 1;
    int ret;
    ADATA adata;
    
    //get number of images
    if(mxGetNumberOfDimensions(prhs[0]) ==4){
        imglen = mxGetDimensions(prhs[0])[3];
    }
//     printf("imglen = %d\n", imglen);
    
    //default opts: [1e-3, 1e-15, 1e-15, 1e-15, 1e-15];
//     opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
//     opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used
    opts = mxGetPr(prhs[1]);
    plpin = mxGetPr(prhs[2]);
//     printf("opts: %E, %E, %E, %E, %E\n",opts[0], opts[1], opts[2], opts[3], opts[3]);
    
    //init result
    plhs[0] = mxCreateDoubleMatrix(imglen,12,false);
    plhs[1] = mxCreateDoubleMatrix(imglen,10,false);
    
    //fit gaussian
    pimg = mxGetPr(prhs[0]);
    plp = mxGetPr(plhs[0]);
    pinfo = mxGetPr(plhs[1]);
    
    for(m=0;m<imglen;m++){
        //init adata
        adata.width = imgwidth;
        adata.x0 = plpin[m];
        adata.y0 = plpin[m + imglen];
        adata.sx = plpin[m + imglen*2];
        adata.sy = plpin[m + imglen*3];
        
        for(subcnt=0; subcnt<6; subcnt++){
            initpar(pimg, imgwidth, lp);
            //dlevmar_dif(func,lp,mxGetPr(prhs[0]),6,7*7,200,opts,info,NULL,NULL,(void *)(&imgwidth));
            ret = dlevmar_der(func, jacf, lp, pimg, 2, imgwidth*imgwidth, 
                    (int)(opts[5]), opts, info, NULL,NULL,(void *)(&adata));

            //copy results
            plp[subcnt*imglen] = lp[0];
            plp[(subcnt+6)*imglen] = lp[1];
            
            
    //         memcpy(plp, lp, 6*sizeof(double));
    //         memcpy(pinfo, info, 10*sizeof(double));

            //update pointers
            pimg += imgwidth*imgwidth;
        }
        for(n=0;n<10;n++){
            pinfo[n*imglen] = info[n];
        }
        plp ++;
        pinfo ++;
    }
}
