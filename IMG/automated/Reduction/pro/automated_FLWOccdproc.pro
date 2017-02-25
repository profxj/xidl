



PRO automated_FLWOccdproc, logstr, skyfit=skyfit, gzip=gzip, clobber=clobber




   ;;open the structures and store them in the array str
   str = REPLICATE(mrdfits(logstr[0],1,/silent),4)
   FOR i=1,3 DO BEGIN

      str[i] = mrdfits(logstr[i],1,/silent)

   ENDFOR


   ;;get the amount and location of the each image type.
   bias = where(str[0].type EQ 'BIA', nbias)
   dark = where(str[0].type EQ 'DAR', ndark)
   flat = where(str[0].type EQ 'FLA', nflats)
   science = where(str[0].type EQ 'IMA', nsci)
 

   ;;if clobber
   if keyword_set(clobber) then spawn, 'rm -f '+plan+'status.str'

   ;;Create the array that will store the status structure
   status = REPLICATE({BIASMED:0,DARK:[0.,0.],FLAT:0,FILFLAT:["","","","","",""],IMG:MAKE_ARRAY(N_ELEMENTS(str[0].name),/INTEGER,VALUE=0)},4)

   ;;Search if the status structure exists, if so store in status array
   FOR i=0,3 DO BEGIN

      aa = FILE_SEARCH(logstr[i] + 'status.str',COUNT=cc)
      IF ( cc EQ 0 ) THEN BEGIN

         k = i + 1
         k = STRTRIM(k,2)
         SPLOG,'Creating status structure ' + k + '...'
         status[i] = {BIASMED:0,DARK:[0.,0.],FLAT:0,FILFLAT:["","","","","",""],IMG:MAKE_ARRAY(N_ELEMENTS(str[i].name),/INTEGER,VALUE=0)}
         MWRFITS,status[i],logstr[i] + "status.str",/CREATE

      ENDIF

   ENDFOR


   ;;----------------
   ;;make bias
   ;;----------------




   ;;check if the bias already exists


   ;;Since there will never (probably) be a case where not all the
   ;;extensions have the medbias flag flagged, there is not really a
   ;;point to this for loop here. So, get rid of this for loop, have
   ;;the entire structure array passed into make_bias and do the
   ;;stitching there. That way the only thing that is popped out to
   ;;here is a completely median combined and stitched image. Then you
   ;;can set the flag on all the structures at once. Although after
   ;;this point I don's think that they are needed. Remember to use
   ;;the correct image size. Remember that the structures are the
   ;;extensions. Doing the combining the above way means that there
   ;;does not need to be any dumby variables passed into the
   ;;subroutines (I think). Get a definite answer to this and check
   ;;for consistency. Ask about whether you should use overscans from
   ;;both sides. Almost done, now get some sleep! see you later!



   ;; If the biasmedian exists for status[0], chances are that it
   ;; exists for all, so use the existing biasmedian.
   IF ( status[0].biasmed EQ 1 ) THEN BEGIN

      splog, "Using exisiting median bias."
      biasmedian = MRDFITS(logstr[0] + 'medbias.fits',0,/silent)
      
   ENDIF   

   ;;Self-explanatory
   IF ( status[0].biasmed EQ 2 ) THEN SPLOG, 'No bias frame. Use overscan only!'
      
   ;;if there is no existing median bias, make it
   IF (nbias GT 0 AND status[0].biasmed EQ 0 ) THEN BEGIN
      
      ;;passed: array of sturctures,undefined(for now,output),array of
      ;;bias location.
      automated_FLWOmakebias,str,biasmedian,bias,PLAN=logstr         
      SPLOG,'Done with median bias'
      status.BIASMED=1

   ENDIF

   IF(nbias EQ 0 AND status[0].biasmed EQ 0) THEN BEGIN
      splog, "No bias frame. Use overscan only"
      status.BIASMED=2
   ENDIF




   FOR i=0,3 DO BEGIN

      MWRFITS,status[i],logstr[i] + 'status.str',/CREATE

   ENDFOR


   ;;---------------
   ;;make dark
   ;;---------------


   IF ( status[0].dark[0] NE 0 ) THEN SPLOG,'Using existing dark...'
   IF ( status[0].dark[0] EQ -9999 ) THEN SPLOG,'No darks found...I can live with that for now!'


   IF ( ndark GT 0 AND status[0].dark[0] EQ 0 ) THEN BEGIN

      IF (status[0].biasmed EQ 1) THEN automated_FLWOmakedark,str,dark,status=status,bias=biasmedian
      IF ( status[0].biasmed EQ 2 ) THEN automated_FLWOmakedark,str,dark,status=status
      SPLOG,'All done with darks'
   ENDIF
 
   IF ( ndark EQ 0 and status[0].dark[0] EQ 0 ) THEN BEGIN

      SPLOG,'No darks found...'
      status.dark[0]=-9999

   ENDIF

   FOR i=0,3 DO BEGIN

      MWRFITS,status[i],logstr[i] + 'status.str',/CREATE

   ENDFOR


 
   ;;set min and max values to accept a flat as good
   minmax = [1000.,41000.]
 

   ;;check if flats already exists
   IF (status[0].flat EQ 1) THEN SPLOG,'Using existing flats.'
   IF (status[0].flat EQ 2) THEN BEGIN ;;these parts need to be tested more (MF June 2012)
      
      SPLOG,'No flats found... Try with archive!'
      automated_FLWOmakeflats,str,flat,science,status=status,plan=logstr,minmax=minmax
      status.flat = 1

   ENDIF

   IF (nflats GT 0 and status[0].flat EQ 0) THEN BEGIN

      IF (status[0].biasmed EQ 1) THEN automated_FLWOmakeflats, str, flat, science,STATUS=status, BIAS=biasmedian, PLAN=logstr,MINMAX=minmax
      IF (status[0].biasmed EQ 2) THEN automated_FLWOmakeflats, str, flat, science,STATUS=status,PLAN=logstr, MINMAX=minmax
      SPLOG, 'All done with flats!'
      status.flat = 1

   ENDIF


   IF (nflats EQ 0 AND status[0].flat EQ 0) THEN BEGIN


      SPLOG,'No flats found...I cannot reduce without flats!'
      SPLOG,'Exiting...Nothing will happen!'
      status.flat = 2
      RETURN

   ENDIF



   ;;--------------------
   ;;make science frame
   ;;--------------------


   IF (nsci GT 0) THEN BEGIN


      ;;loop on each science frame
      FOR i=0,nsci-1 DO BEGIN

         IF (status[0].img[i] EQ 1) THEN SPLOG,'Found previously reduced frame',str[0].name[science[i]] ELSE BEGIN

            IF (status[0].biasmed EQ 1) THEN automated_FLWOmakescience,str,science[i],STATUS=status,PLAN=logstr, BIAS=biasmedian,SKYFIT=skyfir,GZIP=gzip
            IF (status[0].biasmed EQ 2) THEN automated_FLWOmakescience,str, science[i],STATUS=status,PLAN=logstr,SKYFIT=skyfit,GZIP=gzip


            ;;update structure
            status.img[i] = 1
            FOR j=0,3 DO BEGIN
               MWRFITS,status[j],logstr[j]+'status.str',/CREATE
            ENDFOR

         ENDELSE

      ENDFOR
      

   ENDIF ELSE BEGIN

      PRINT,'No science frame!'
      RETURN

      SPLOG,'All done with data reduction!'
      
   ENDELSE
   

END
