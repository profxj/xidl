






PRO automated_FLWOloopflats, str, nfilt, thisfilter, whereFlats, flat_chip, BIAS=bias, MINMAX=minmax, STATUS=status


   ;;make storage
   flat1 = MAKE_ARRAY(nfilt,2048,2048,/FLOAT)
   dump = 0
   
   chips = PTRARR(4)
   heas = PTRARR(4)

   ;;now loop over all the images
   FOR i=0,nfilt-1 DO BEGIN
      
      ;;open
      SPLOG, 'Working on image',str[0].name[whereFlats[thisFilter[i]]]
      FOR j=0,3 DO BEGIN
         
         chips[j] = PTR_NEW(MRDFITS(STRTRIM(str[0].name[whereFlats[thisFilter[i]]],2),j+1,hea,/SILENT,/FSCALE))
         heas[j] = PTR_NEW(hea)

      ENDFOR

      IF (DJS_MEDIAN(*chips[0]) LT minmax[1] AND DJS_MEDIAN(*chips[0]) GT minmax[0]) THEN BEGIN

         ;;go for oscan
         automated_FLWOoscan,chips,imageout
         
         ;;remove zero
         IF KEYWORD_SET(BIAS) THEN BEGIN

            *imageout[0] = *imageout[0] - bias[1024:2047,0:1023]
            *imageout[1] = *imageout[1] - bias[0:1023,0:1023]
            *imageout[2] = *imageout[2] - bias[1024:2047,1024:2047]
            *imageout[3] = *imageout[3] - bias[0:1023,1024:2047]

         ENDIF

         ;;remove dark current
         IF (status[0].dark[0] NE -9999) THEN BEGIN

            FOR k=0,3 DO BEGIN

               darval = status[k].dark[0]*str[k].extime[whereFlats[thisFilter[i]]] + status[k].dark[1]
               SPLOG,'Subtract dark current ',darkval
               *imageout[k] = *imageout[k] - darkval 

            ENDFOR

         ENDIF
         
         ;;go for stitching
         FOR q=0,3 DO BEGIN
            ;;normal action
            *imageout[q] = *imageout[q]/DJS_MEDIAN(*imageout[q])

         ENDFOR

         automated_imagestitch,imageout,stitchedFlat

         ;;store chip after normalization for flats
         flat1[i-dump,*,*] = stitchedFlat[*,*]


      ENDIF ELSE BEGIN

         ;;take out spit for saturated image
         SPLOG,'Dunmo flat ',str[0].name[whereFlats[thisFilter[i]]]
         ;;if bad flat, then keep same index 'img - dump'
         dump = dump + 1
         flat1 = REFORM(flat1[0:N_ELEMENTS(flat1[*,0,0]) - 2,*,*])

      ENDELSE

   ENDFOR
   
   
   ;; If there are no good flats then find last good ones and try this
   ;; over again
   IF dump EQ nfilt THEN BEGIN

      SPLOG,'No working flats. Extracting archive.'
      ;;Read in archive of superflat
      flat_chip = MRDFITS(GETENV('AUTOM_DIR') + '/archive/archivedflat' + str[0].filter[whereFlats[thisFilter[0]]] + '.fits')

   ENDIF ELSE BEGIN

      ;;here stack them with a simple median
      flat_chip = DJS_MEDIAN(flat1,1)
      ;;write to the archive for good flats
      MWRFITS,flat_chip, GETENV('AUTOM_DIR') + "/archive/archivedflat" + str[0].filter[whereFlats[thisFilter[0]]] + '.fits',/CREATE

   ENDELSE


END
