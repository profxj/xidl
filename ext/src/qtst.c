#include <stdio.h>
#include "export.h"

/* float vector_median
  (IDL_LONG   nData,
   float *    pData); */

/******************************************************************************/
int qtst
  (int         argc,
   void    *   argv[])
{
   float *a;
   unsigned char      *dum;
   float *b;

   a = (float *) argv[0]; 
   dum = (unsigned char *) argv[1]; 
   b = (float *) argv[2]; 

   b[1] = 1.;

   return 1;
}

