#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

/******************************************************************************/
IDL_LONG rectify
  (int         argc,
   void    *   argv[])
{
   IDL_LONG    ndim;
   IDL_LONG  *  dim;
   double    *  val;
   double    *  org;
   double    *  map;

   int         i,j,jm,jp;
   IDL_LONG    retval = 1;
   int         indx;
   double      fp,fm, diff;
   /* FILE     *  list;  */

   /* Allocate pointers from IDL */
   ndim = *((IDL_LONG *)argv[0]);
   dim = (IDL_LONG *)argv[1];
   val = (double *)argv[2];
   org = (double *)argv[3];
   map = (double *)argv[4];

   /* Allocate memory for temporary vectors */

   /* list = fopen("fort.13","w+");
   fprintf(list,"Entering the big loop! \n"); */

   /* Loop through all pixels in array */
   for (i=0; i < dim[0]; i++) {
      for (j=0; j < dim[1]; j++) {
	indx = i + j*dim[0];
	diff = (double)j + map[indx];
	jm = (int)diff;
	jp = jm+1;
	fm = 1. - (diff - (double)jm);
	fp = 1. - fm;
	/* Conditionals */
	if(jm >= 0 && jm < dim[1]) val[i+jm*dim[0]] += org[indx]*fm;
	if(jp >= 0 && jp < dim[1]) val[i+jp*dim[0]] += org[indx]*fp;
      }
   }

   /* fclose(list); */
   return retval;
}
