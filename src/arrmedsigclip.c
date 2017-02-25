#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

float vector_median
  (IDL_LONG   nData,
   float *    pData);

/******************************************************************************/
IDL_LONG arrmedsigclip
  (int         argc,
   void    *   argv[])
{
   IDL_LONG    ndim;
   IDL_LONG *  dimvec;
   float    *  array;
   IDL_LONG    dim;
   float       sigrejlo;
   float       sigrejhi;
   IDL_LONG    maxiter;
   float    *  medarr;
   float       gain;
   float       rn;

   IDL_LONG    i;
   IDL_LONG    nlo;
   IDL_LONG    nmid;
   IDL_LONG    nhi;
   IDL_LONG    ilo;
   IDL_LONG    imid;
   IDL_LONG    ihi;
   IDL_LONG    i1;
   IDL_LONG    indx;
   float    *  tempvec;
   float    *  goodvec;
   IDL_LONG    retval = 1;
   float       sigma;
   float       mval;
   IDL_LONG    ngood;
   IDL_LONG    ngsv;
   IDL_LONG    niter;

   /* Allocate pointers from IDL */
   ndim = *((IDL_LONG *)argv[0]);
   dimvec = (IDL_LONG *)argv[1];
   array = (float *)argv[2];
   dim = *((IDL_LONG *)argv[3]);
   sigrejlo = *((float *)argv[4]);
   sigrejhi = *((float *)argv[5]);
   maxiter = *((IDL_LONG *)argv[6]);
   gain = *((float *)argv[7]);
   rn = *((float *)argv[8]);
   medarr = (float *)argv[9];

   nlo = 1;
   for (i=0; i < dim-1; i++) nlo *= dimvec[i];
   nhi = 1;
   for (i=dim; i < ndim; i++) nhi *= dimvec[i];
   nmid = dimvec[dim-1];

   /* Allocate memory for temporary vectors */
   tempvec = malloc(nmid * sizeof(float));
   goodvec = malloc(nmid * sizeof(float));

   /* Loop through all pixels in array */
   for (ilo=0; ilo < nlo; ilo++) {
      for (ihi=0; ihi < nhi; ihi++) {
         i1 = ilo + ihi * nmid * nlo;
         /* Construct the vector of values from which to compute */
         for (imid=0; imid < nmid; imid++) {
            indx = i1 + imid * nlo;
            tempvec[imid] = array[indx];
         }
         /* Compute the first median value */
         mval = vector_median(nmid, tempvec);

         /* Compute the first sigma value*/
	 sigma = sqrt(abs(mval/gain + rn*rn/gain/gain));

	 /* Iterate with rejection */
	 niter = 0;
	 ngsv = 0;
	 while(niter < maxiter) {
	   ngood = 0;
	   for (imid=0; imid < nmid; imid++) {
	     if(tempvec[imid] > mval - sigrejlo*sigma &&
		tempvec[imid] < mval + sigrejhi*sigma) 
	       goodvec[ngood++] = tempvec[imid];
	   }
	   if(ngood == 0 || ngood == ngsv) {
	     niter = maxiter;
	   } else {
	     if(ngood == 1) {    /* Catch single number array */
	       mval = goodvec[0];
	       niter = maxiter;
	     } else {
	       /* Compute the median value */
	       mval = vector_median(ngood, goodvec);
	       /* Compute the sigma value*/
	       sigma = sqrt(abs(mval/gain + rn*rn/gain/gain));
	       /* Increment niter */
	       niter++;
	       /* Save ngood */
	       ngsv = ngood;
	     }
	   }
	 }
	 /* Set the medarr value */
         medarr[ilo + ihi * nlo] = mval;
      }
   }

   /* Free temporary memory */
   free(tempvec);
   free(goodvec);

   return retval;
}

