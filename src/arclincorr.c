#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

/******************************************************************************/
IDL_LONG arclincorr
  (int         argc,
   void    *   argv[])
{
   IDL_LONG    ndim;
   float    *  all_chi;
   IDL_LONG *  all_nsub;
   IDL_LONG *  nguess;
   IDL_LONG    szguess;
   double   *  guess;
   IDL_LONG    nlin;
   double   *  center;
   IDL_LONG    nmod;
   double   *  model;
   double   *  weight;
   IDL_LONG    npix;

   int         i,j,k,l,m, nsub;
   IDL_LONG    retval = 1;
   double      zro, avdlmb, dlmb, nonl,fnpix;
   double      mxchi;
   double      t1, t2, mxwv, mnwv, chisq, mnsep, sep, wt;
   double   *  subwv;
   double   *  subwt;
   double   *  cent2;
   double   *  cent3; 
   double   *  cent4; 
   double   *  lambda;
   /* FILE     *  list; */

   /* Allocate pointers from IDL */
   ndim = *((IDL_LONG *)argv[0]);
   all_chi = (float *)argv[1];
   all_nsub = (IDL_LONG *)argv[2];
   nguess = (IDL_LONG *)argv[3];
   szguess = *((IDL_LONG *)argv[4]);
   guess = (double *)argv[5];
   nlin = *((IDL_LONG *)argv[6]);
   center = (double *)argv[7];
   nmod = *((IDL_LONG *)argv[8]);
   model = (double *)argv[9];
   weight = (double *)argv[10];
   npix = *((IDL_LONG *)argv[11]);


   /* Allocate memory for temporary vectors */
   subwv = malloc(nmod * sizeof(double));
   subwt = malloc(nmod * sizeof(double));
   cent2 = malloc(nlin * sizeof(double));
   cent3 = malloc(nlin * sizeof(double));
   cent4 = malloc(nlin * sizeof(double));
   lambda = malloc(nlin * sizeof(double));

   /* Create cent2, cent3 */
   for (i=0; i < nlin; i++) {
     cent2[i] = center[i]*center[i];
     cent3[i] = cent2[i]*center[i];
     if(szguess > 4) cent4[i] = cent3[i]*center[i];
   }
   /* Set fnpix */
   fnpix = (double)npix;

   /*
   cnt = 0;
   list = fopen("fort.13","w+");
   fprintf(list,"Entering the big loop! \n");
   for(i=0; i < nmod; i++) {
     dumf = (double)model[i];
     fprintf(list,"%f \n", dumf);
   }
   */
     
   /* Big Loop */
   for (i=0; i < nguess[0]; i++) {
     zro = guess[0] + (double)i*(guess[1]-guess[0])/(double)nguess[0];
     for (j=0; j < nguess[1]; j++) {
       dlmb = guess[2] + (double)j*(guess[3]-guess[2])/(double)nguess[1];
       for (k=0; k < nguess[2]; k++) {
	 nonl = guess[4] + (double)k*(guess[5]-guess[4])/(double)nguess[2];

	 /* Compute lambda */
	 for (l=0; l<nlin; l++) {
	   lambda[l] = zro + center[l]*dlmb + cent2[l]*nonl + cent3[l]*guess[6];
	   if(szguess > 4) lambda[l] = lambda[l] + cent4[l]*guess[8];
	 }

	 /* Find the subset of model lines */
	 t1 = zro;
	 t2 = zro + fnpix*dlmb + (fnpix*fnpix)*nonl + (fnpix*fnpix*fnpix)*guess[6];
	 if(szguess > 4) t2 = t2+(fnpix*fnpix*fnpix*fnpix)*guess[8];

	 if(t1 > t2) {
	   mxwv = t1;
	   mnwv = t2;
	 } else {
	   mxwv = t2;
	   mnwv = t1;
	 }
	 nsub = 0;
	 for(l=0; l<nmod; l++) {
	   if(model[l] > mnwv && model[l] < mxwv) {
	     subwv[nsub] = model[l];
	     subwt[nsub] = weight[l];
	     nsub++;
	   }
	 }
	 /* Average delta lambda */
	 avdlmb = (mxwv-mnwv)/fnpix;
	   
         /* Compute chisq */
	 chisq = 0.;
	 mxchi = 0.;
	 for(l=0; l<nsub; l++) {
	   mnsep = 9999.;
	   wt=0.;
	   for(m=0; m<nlin; m++) {
	     sep = fabs(lambda[m]-subwv[l]);
	     if( sep < mnsep) {
	       mnsep = sep;
	       wt = subwt[l];
	     }
	   }
	   t1 = mnsep*mnsep/(avdlmb*avdlmb)/4.;
	   if(t1 > mxchi) mxchi = t1;
	   chisq +=  t1;
	 }
	 /* Reject highest chi */
	 chisq = chisq - mxchi;

	 /* Save to the array */
	 all_chi[i+j*nguess[0]+k*(nguess[0]*nguess[1])] = chisq/(double)(nsub-1);
	 /* all_nsub[i + j*nguess[0] + k*(nguess[0]*nguess[1])] = nsub; */
       }
     }
   }
	   

   /* Free temporary memory */
   free(subwv);
   free(cent2);
   free(cent3);
   free(lambda);
   /* fclose(list); */

   return retval;
}

