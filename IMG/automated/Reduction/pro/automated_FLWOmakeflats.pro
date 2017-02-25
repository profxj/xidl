

















PRO automated_FLWOmakeflats,str,whereFlats,bias=bias,whereScience,status=status,plan=logstr,minmax=minmax



   FOR k=0,99 DO BEGIN

      IF ( k EQ 0 ) THEN BEGIN 

         sciFilters = STRCOMPRESS(str[0].filter[whereScience[k]],/REMOVE_ALL)
         flaFilters = STRCOMPRESS(str[0].filter[whereFlats[k]],/REMOVE_ALL)

      ENDIF

      IF ( k GT 0 ) THEN BEGIN
         
         IF (k LT N_ELEMENTS(whereScience)) THEN BEGIN

            sciRes = WHERE(sciFilters EQ STRCOMPRESS(str[0].FILTER[whereScience[k]],/REMOVE_ALL),sciMatch)
            IF ( sciMatch EQ 0 ) THEN sciFilters = [sciFilters,STRCOMPRESS(str[0].FILTER[whereScience[k]],/REMOVE_ALL)]

         ENDIF

         IF (k LT N_ELEMENTS(whereFlats)) THEN BEGIN

            flaRes = WHERE(flaFilters EQ STRCOMPRESS(str[0].FILTER[whereFlats[k]],/REMOVE_ALL),flaMatch)
            IF ( flaMatch EQ 0 ) THEN flaFilters = [flaFilters,STRCOMPRESS(str[0].FILTER[whereFlats[k]],/REMOVE_ALL)]

         ENDIF


      ENDIF           
      
   ENDFOR

   FOR j=0,N_ELEMENTS(sciFilters)-1 DO BEGIN

      IF ( WHERE(flaFilters EQ sciFilters[j]) EQ -1 ) THEN BEGIN
         SPLOG,'There are science images filters that dont have corresponding flats! ' + sciFilters[j]
         SPLOG,'Reading in an archived flat and saving it for this night.'
         finalflat = MRDFITS(GETENV('AUTOM_DIR') + 'archive/archivedflat' + STRCOMPRESS(sciFilters[j],/REMOVE_ALL) + '.fits')
         MWRFITS, finalflat, sciFilters[j] + logstr[0] +'flat.fits',/CREATE
         indx = WHERE(status[0].filflat EQ "")
         status.filflat[indx] = sciFilters[j]

         FOR p=0,3 DO BEGIN
            
            MWRFITS,status[0],logstr[p] + 'status.str',/CREATE

         ENDFOR

      ENDIF


   ENDFOR


 






   ;;if no flats,archive
   IF ( whereFlats[0] EQ -1 ) THEN BEGIN

      ;;find group of common filters
      FOR i=0,N_ELEMENTS(str[0].filter)-1 DO BEGIN

         ;;may need to add an outer for loop and change str[0] to str[j]
         IF ( i EQ 0 ) THEN flat = STRCOMPRESS(str[0].filter[i],/REMOVE_ALL)
         IF ( i GT 0 ) THEN BEGIN
 
            ;;check if new filter
            res = WHERE(flat EQ STRCOMPRESS(str[0].FILTER[i],/REMOVE_ALL),numMatch)
            ;;if new, append
            IF ( numMatch EQ 0 ) THEN flat = [flat,STRCOMPRESS(str[0].FILTER[i],/REMOVE_ALL)]

         ENDIF

      ENDFOR

      FOR flfra = 0, N_ELEMENTS(flat)-1 DO BEGIN

         finalflat = MRDFITS(GETENV('AUTOM_DIR') + 'archive/FLWO_archivedflat' + STRCOMPRESS(flat[i],/REMOVE_ALL) + '.fits')
         ;;save final output
         MWRFITS, finalflat, flat[i] + logstr[0] +'flat.fits',/CREATE
         ;;update status structure
         status.filflat[i]=flat[i]
         FOR i=0,3 DO MWRFITS,status[0],logstr[i] + 'status.str',/CREATE

      ENDFOR

   ENDIF ELSE BEGIN

      ;;else, loop over flats to find groups of common filters
      FOR i=0,N_ELEMENTS(whereFlats)-1 DO BEGIN

         IF (i EQ 0) THEN flat=STRCOMPRESS(str[0].filter[whereFlats[i]],/REMOVE_ALL)
         IF (i GT 0) THEN BEGIN
            
            ;;check if new filter
            res = WHERE(flat EQ STRCOMPRESS(str[0].filter[whereFlats[i]],/REMOVE_ALL),numMatch)
            ;;if new,append
            IF (numMatch EQ 0) THEN flat = [flat,STRCOMPRESS(str[0].filter[whereFlats[i]],/REMOVE_ALL)]

         ENDIF

      ENDFOR

      SPLOG,'Filters found: ',flat





      ;;make flats
      FOR i=0,N_ELEMENTS(flat)-1 DO BEGIN

         ;;check if flat done
         sss = WHERE(STRCOMPRESS(status[0].filflat,/REMOVE_ALL) EQ flat[i],numbfl)

         IF (numbfl GT 0) THEN SPLOG,'Found flat ',flat[i] ELSE BEGIN

            SPLOG,'WORKING on filter ',flat[i]

            ;;find images for this filter
            thisFilter = WHERE(STRCOMPRESS(str[0].filter[whereFlats],/REMOVE_ALL) EQ flat[i],nfilt)

            ;;check if everything is ok (nfilt MUST be > 0)
            IF (nfilt EQ 0) THEN BEGIN

               SPLOG,'I cannot finf images with this filter: ',flat[i],". Stopping..."
               STOP

            ENDIF

            ;;call loop flats to make final image
            automated_FLWOloopflats, str, nfilt, thisFilter, whereFlats, finalflat, BIAS=bias, MINMAX=minmax, STATUS=status

            ;;save final output
            MWRFITS, finalflat, flat[i] + logstr[0] + 'flat.fits',/CREATE
            CD,GETENV('AUTOM_DIR') + 'archive',CURRENT=here
            MWRFITS, finalflat,'FLWO_archivedflats' + STRCOMPRESS(flat[i],/REMOVE_ALL) + '.fits',/CREATE
            CD,here

            ;;update status structure
            FOR k=0,3 DO BEGIN
               status[k].filflat[i] = flat[i]
               MWRFITS,status[k],logstr[k] + 'status.str',/CREATE
            ENDFOR

         ENDELSE

      ENDFOR


   ENDELSE


   SPLOG,'All done with flats!!'



END
