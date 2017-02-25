






PRO automated_FLWOmakebias,str,biasmedian,whereBias,plan=logstr


  
   IF ( N_ELEMENTS(whereBias) GT 10) THEN BEGIN

      whereBias = whereBias[0:10]
      SPLOG,'Considering only 10 bias frames!'

   ENDIF

   ;;make pointer array that will hold the location of each of the
   ;;stitched median combined extensions.

   biases = MAKE_ARRAY(N_ELEMENTS(whereBias),2048,2048,/FLOAT)


   ;;make pointer arrays that will hold the location of the raw data
   ;;extensions and the associated headers.
   chips = PTRARR(4)
   heas = PTRARR(4)

   ;;cycle over each seperate bias image...
   FOR i=0, N_ELEMENTS(whereBias) - 1 DO BEGIN

      ;;open files, tells what object looking at
      SPLOG,'Working on bias ',str[0].name[whereBias[i]]
      ;;read in the 4 extensions into chips
      FOR j=1,4 DO BEGIN
        
         chips[j-1] = PTR_NEW(MRDFITS(STRTRIM(str[0].name[whereBias[i]],2),j,hea,/SILENT,/FSCALE))
         heas[j-1] = PTR_NEW(hea)

      ENDFOR
      
      
      ;;pass chips into oscan for subtractions and resizings,pass back
      ;;biasimages (a pointer array)
      automated_FLWOoscan,chips,biasimage
      ;;
      automated_imagestitch,biasimage,stitchedBias
      
      biases[i,*,*] = stitchedBias[*,*]


   ENDFOR

  ;;combine the final bias
  splog, "Creating median bias... "  
  
  biasmedian = DJS_MEDIAN(biases,1) 
  
  ;;save median bias
  MWRFITS, biasmedian, logstr[0] + "medbias.fits", /create

  splog, "All done with median bias!"


END
