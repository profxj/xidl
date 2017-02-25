





PRO automated_FLWOoscan,chips,dataimage



   ;;empty so far,this will be the output
   dataimage = PTRARR(4)
   FOR i=0,3 DO BEGIN

      dataImage[i] = PTR_NEW(MAKE_ARRAY(1024,1024,/FLOAT))

   ENDFOR
   
   
 
   ;;this is the image minus the osca


   ;;get bias info, this is the oscan array, reform apparently does
   ;;nothing


   ;;djs_median collapses an array along the specified axis by median
   ;;combining. So a [91,2048] array collapsed along 1 gives a [2048]
   ;; array combined across the rows.

   ;;undefined until now, transpose makes bia 2-d, [1,1024]. Rebin
   ;;makes it [1024,1024]

   


   FOR i=0,3 DO BEGIN

      ;;holder array is just a temp array because I suck at
      ;;pointers. Data is the image without the oscan. oscan is the
      ;;oscan. bia is the collapsed oscan. bias is the median combined
      ;;oscan stretched to the size of data. Dataimage is the
      ;;delicious final, subtracted, result.
      holderArray = *chips[i]
      data = REFORM(holderArray[8:1031,0:1023])
      holderArray = *chips[i]
      overscan = reform(holderArray[0:7,0:1023])
      bia = DJS_MEDIAN(overscan,1)
      bias = REBIN(TRANSPOSE(bia),1024,1024)
      *dataImage[i] = data - bias

   ENDFOR
  

END
