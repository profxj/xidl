#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

/******************************************************************************/
IDL_LONG x_gsmooth
  (int         argc,
   void    *   argv[])
{
   IDL_LONG    npix;
   IDL_LONG    nsub;
   double    *  wave;
   double    *  subwv;
   double    *  dwv;
   float     *  subfx;
   float     *  fx;
   float      FWHM;

   int         i,j;
   IDL_LONG    retval = 1;
   double      sum1, sum2, sigma, invsig;
   double      dffwv;

   double      gaussf;
   /* FILE     *  list;  */

   /* Allocate pointers from IDL */
   npix = *((IDL_LONG *)argv[0]);
   nsub = *((IDL_LONG *)argv[1]);
   wave = (double *)argv[2];
   subwv = (double *)argv[3];
   dwv = (double *)argv[4];
   subfx = (float *)argv[5];
   fx = (float *)argv[6];
   FWHM = *((float *)argv[7]);

   /* Allocate memory for temporary vectors */

   /* list = fopen("fort.13","w+");
   fprintf(list,"Entering the big loop! \n"); */


   /* Loop through all pixels in array */
   for (i=0; i < npix; i++) {
     /* Calucate sigma */
     sigma = FWHM / 2.3548 * dwv[i];
     invsig = -0.5 / sigma / sigma;

     /* Smooth */
     sum1 = 0.;
     sum2 = 0.;
     for (j=0; j < nsub; j++) {
       dffwv = (subwv[j]-wave[i]);
       if(abs(dffwv) < 7.*sigma) {
	 gaussf = exp(invsig * dffwv * dffwv);
	 sum1 = sum1 + gaussf * subfx[j];
	 sum2 = sum2 + gaussf;
       } 
     }
     fx[i] = sum1 / sum2;
   }

   /* fclose(list); */
   return retval;
}
