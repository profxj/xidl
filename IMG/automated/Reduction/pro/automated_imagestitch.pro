






PRO automated_imagestitch,biases,stitched


   PRINT,'Getting patches...'
   
   sizeArray = SIZE(*biases[0])
   stitched = FLTARR(2*sizeArray[1],2*sizeArray[2])
   stitched[1024:2047,0:1023] = *biases[2]
   stitched[0:1023,0:1023] = *biases[3]
   stitched[1024:2047,1024:2047] = *biases[0]
   stitched[0:1023,1024:2047] = *biases[1]

END
