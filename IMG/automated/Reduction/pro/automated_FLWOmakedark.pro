





PRO automated_FLWOmakedark,str,whereDark,BIAS=bias,STATUS=status

   
   ;;start here
   ndark = N_ELEMENTS(whereDark)
   ;;make storage
   dark = MAKE_ARRAY(ndark,2048,2048,/FLOAT)

   valueDark  = PTRARR(4)

   FOR i=0,3 DO BEGIN

      valueDark[i] = PTR_NEW(FLTARR(ndark))

   ENDFOR

   valueTime = FLTARR(ndark)

   chips = PTRARR(4)
   heas = PTRARR(4)

   ;;now loop over all the images
   FOR i=0,ndark-1 DO BEGIN

      ;;open
      SPLOG,'Working on ',str.name[whereDark[i]]
      FOR j=1,4 DO BEGIN

         chips[j-1] = PTR_NEW(MRDFITS(STRTRIM(str[0].name[whereDark[i]],2),j,hea,/SILENT,/FSCALE))
         heas[j-1] = PTR_NEW(hea)

      ENDFOR

      STOP
      supl_FLWOoscan,chips,imageout

      IF KEYWORD_SET(bias) THEN BEGIN

         *imageout[0] = *imageout[0] - bias[1024:2047,0:1023]
         *imageout[1] = *imageout[1] - bias[0:1023,0:1023]
         *imageout[2] = *imageout[2] - bias[1024:2047,1024:2047]
         *imageout[3] = *imageout[3] - bias[0:1023,1024:2047]

      ENDIF


      FOR j=0,3 DO BEGIN

         holdArray = *valueDark[j]
         holdArray[i] = DJS_MEDIAN(*imageout[j])
         *valueDark[j] = holdArray
         
      ENDFOR

      valueTime[i] = str[0].extime[whereDark[i]]


   ENDFOR
   
   
   
   FOR i=0,3 DO BEGIN

      res = linfit(valueTime,*valueDark[i],YFIT=fit)
      status[i].dark[0] = res[1]
      status[i].dark[1] = res[0]


   ENDFOR




END
