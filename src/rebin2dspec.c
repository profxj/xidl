#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

/******************************************************************************/
IDL_LONG rebin2dspec
  (int         argc,
   void    *   argv[])
{
   IDL_LONG    ndim;
   IDL_LONG  *  dim;
   IDL_LONG  *  gdpix;
   double    *  wvl;
   double    *  wvh;
   double    *  bwv;
   float     *  fx;
   double    *  var;
   float     *  nwfx;
   double    *  nwvar;

   IDL_LONG    i,j;
   IDL_LONG    retval = 1;
   IDL_LONG    q, i1, i2, kk;
   int         flg_cont;
   double      frac, fmin, fmax;

   /* Allocate pointers from IDL */
   ndim = *((IDL_LONG *)argv[0]);
   dim = (IDL_LONG *)argv[1];
   gdpix = (IDL_LONG *)argv[2];
   wvl = (double *)argv[3];
   wvh = (double *)argv[4];
   bwv = (double *)argv[5];
   fx = (float *)argv[6];
   var = (double *)argv[7];
   nwfx = (float *)argv[8];
   nwvar = (double *)argv[9];

   /* Initial values */
   i1 = 0;
   i2 = 0;
   frac = 0.;

   /* Loop through all pixels in array */
   for (i=0; i < dim[0]; i++) {
     q = gdpix[i];
     flg_cont = 0;
     if(q == -1) flg_cont = 1;
     if(wvh[q] <= bwv[0] || wvl[q] >= bwv[dim[1]]) flg_cont = 1;
     /* No zeros */
     if (wvh[q] * wvl[q] == 0.) flg_cont = 1;

     /* Main stuff */
     if(flg_cont == 0) {
       /* Find pixel that bw is within */
       if(wvl[q] < bwv[0]) i1 = 0;
       else {
	 for(j=0; j<dim[1]; j++) {
            if(wvl[q] > bwv[j] && wvl[q] <= bwv[j+1]) {
	      i1 = j;
	      break;
	    }
	 }
       }
       if(wvh[q] > bwv[dim[1]]) i2 = dim[1]-1; 
       else {
	 for(j=i1; j<dim[1]; j++){
	   if(wvh[q] > bwv[j] && wvh[q] <= bwv[j+1]) {
	     i2 = j;
	     break;
	   }
	 }
       }
       /* Now Sum up */
       for(kk=i1; kk<=i2; kk++) {
	 /* Rejected pixels do not get added in */
	 if(var[q] > 0.) {
	   if(wvh[q] < bwv[kk+1]) fmin = wvh[q]; 
	   else fmin = bwv[kk+1];
	   if(wvl[q] > bwv[kk]) fmax = wvl[q]; 
	   else fmax = bwv[kk];
	   frac = ( fmin - fmax ) / (wvh[q]-wvl[q]);
	   nwfx[kk] = nwfx[kk] + frac * fx[q];
	 }
         /* Variance */
	 if(nwvar[kk] != -1.) {
	   if(var[q] <= 0.) {
	     if(dim[2] == 1) {
	       if(var[q] == -1.) nwvar[kk] = -1.;
	     }
	   }
	   else {
	     nwvar[kk] = nwvar[kk] + frac * var[q];
	   }
	 }
       }
     }
       
   }

   return retval;
}
