#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"
#include "nr.h"

float vector_median
  (IDL_LONG   nData,
   float *    pData);

/******************************************************************************/
/* Find the median of the "nData" elements of a floating point array pData[].
 * The data vector is returned unchanged.
 */
float vector_median
  (IDL_LONG   nData,
   float *    pData)
{
   unsigned long nlong;
   unsigned long klong;
   float    retval;

   /* Numerical Recipes routines expect everything to be 1-indexed.
    * Pass the array pData as 1-indexed, and the element counter klong
    * as a 1-indexed element number.
    */
   nlong = nData;
   klong = (nlong+1)/2;
   retval = selip(klong, nlong, pData-1);
   if (nData % 2 == 0) {
      retval = 0.5 * (retval + selip(klong+1, nlong, pData-1));
   }

   return retval;
}

