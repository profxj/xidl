#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

float vector_median
  (IDL_LONG   nData,
   float *    pData);

/******************************************************************************/
IDL_LONG arrmedsclipwgt
  (int         argc,
   void    *   argv[])
{
   IDL_LONG    ndim;
   IDL_LONG *  dimvec;
   float    *  array;
   float    *  var;
   IDL_LONG    dim;
   float       sigrej;
   float    *  medarr;
   float    *  finvar;
   float    *  weight;
   float    *  scale;

   IDL_LONG    i,j;
   IDL_LONG    nlo;
   IDL_LONG    nmid;
   IDL_LONG    nhi;
   IDL_LONG    ilo;
   IDL_LONG    imid;
   IDL_LONG    ihi;
   IDL_LONG    i1;
   IDL_LONG    indx;
   float       f1, f2, v1;
   float    *  tempvec;
   float    *  tempvar;
   float    *  goodvec;
   IDL_LONG    retval = 1;
   IDL_LONG *  svindx;
   /* float       sigma; */
   float       mval;
   /* float       mvar; */
   IDL_LONG    ngood;
   IDL_LONG    ngsv;
   /*   IDL_LONG    niter; */

   /* Allocate pointers from IDL */
   ndim = *((IDL_LONG *)argv[0]);
   dimvec = (IDL_LONG *)argv[1];
   array = (float *)argv[2];
   var = (float *)argv[3];
   dim = *((IDL_LONG *)argv[4]);
   sigrej = *((float *)argv[5]);
   medarr = (float *)argv[6];
   finvar = (float *)argv[7];
   weight = (float *)argv[8];
   scale = (float *)argv[9];

   nlo = 1;
   for (i=0; i < dim-1; i++) nlo *= dimvec[i];
   nhi = 1;
   for (i=dim; i < ndim; i++) nhi *= dimvec[i];
   nmid = dimvec[dim-1];

   /* Allocate memory for temporary vectors */
   tempvec = malloc(nmid * sizeof(float));
   tempvar = malloc(nmid * sizeof(float));
   goodvec = malloc(nmid * sizeof(float));
   svindx = malloc(nmid * sizeof(IDL_LONG));

   /* Loop through all pixels in array */
   for (ilo=0; ilo < nlo; ilo++) {
      for (ihi=0; ihi < nhi; ihi++) {
         i1 = ilo + ihi * nmid * nlo;
	 /* Initialize */
	 medarr[ilo + ihi * nlo] = 0.;
	 finvar[ilo + ihi * nlo] = -1.;
         /* Construct the vector of values from which to compute */
	 ngood = 0;
         for (imid=0; imid < nmid; imid++) {
            indx = i1 + imid * nlo;
	    if(var[indx] > 0.) {
	      tempvec[ngood] = array[indx];
	      tempvar[ngood] = var[indx];
	      svindx[ngood] = imid;
	      ngood++;
	    }
         }
         /* Compute the median value (avoid var < 0) */
	 if(ngood > 0) {
	   mval = vector_median(ngood, tempvec);
	   /* Initialize */
	   f1 = 0.;
	   f2 = 0.;
	   v1 = 0.;
	   ngsv = 0;
	   
	   /* Reject */
	   for (i=0; i < ngood; i++) {
	     if(abs(tempvec[i] - mval) < sigrej*sqrt(tempvar[i])*scale[svindx[i]]) {
	       j = svindx[i];
	       f1 += tempvec[i]*weight[j];
	       f2 += weight[j];
	       v1 += tempvar[i]*weight[j]*weight[j]*scale[j]*scale[j];
	       ngsv++;
	     }
	   }
	   if(ngsv > 0) {
	     medarr[ilo + ihi * nlo] = f1/f2;
	     finvar[ilo + ihi * nlo] = v1/f2/f2; 
	   } 
	 }
      }
   }

   /* Free temporary memory */
   free(tempvec);
   free(goodvec);

   return retval;
}

