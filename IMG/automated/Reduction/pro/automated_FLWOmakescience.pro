







PRO automated_FLWOmakescience,str,whereScience,STATUS=status,PLAN=logstr,BIAS=biasmedian,SKYFIT=skyfit,GZIP=gzip


   

   SPLOG,'Processing sci frame ', str[0].name[whereScience]


   ;;open file
   chips = PTRARR(4)
   headers = PTRARR(4)
   
   
   FOR i=0,3 DO BEGIN
      
      chips[i] = PTR_NEW(MRDFITS(STRTRIM(str[0].name[whereScience],2),i+1,header,/SILENT,/FSCALE))
      headers[i] = PTR_NEW(header)
      
   ENDFOR
   PRINT,whereScience
   ;;go for oscan subtraction
   automated_FLWOoscan,chips,imageout


   ;;remove zero if set
   IF KEYWORD_SET(BIAS) THEN BEGIN

      *imageout[0] = *imageout[0] - bias[1024:2047,0:1023]
      *imageout[1] = *imageout[1] - bias[0:1023,0:1023]
      *imageout[2] = *imageout[2] - bias[1024:2047,1024:2047]
      *imageout[3] = *imageout[3] - bias[0:1023,1024:2047]      

   ENDIF

   ;;remove dark current
   IF (status[0].dark[0] NE -9999) THEN BEGIN
      FOR i=0,3 DO BEGIN
         darkval = status[i].dark[0]*str[i].extime[whereScience] + status[i].dark[1]
         SPLOG,'Subtract dark current ',darkval
         *imageout[i] = *imageout[i] - darkval
         
      ENDFOR
   
      
   ENDIF

   ;;find corresponding flat field
   flat = WHERE(STRCOMPRESS(status[0].filflat,/REMOVE_ALL) EQ STRCOMPRESS(str[0].filter[whereScience],/REMOVE_ALL))
   nameFlat = STRCOMPRESS(status[0].filflat[flat],/REMOVE_ALL) + logstr[0] + 'flat.fits'
   flatfits = MRDFITS(nameflat[0],/SILENT)


   ;;do stuff over the chip
   ;;flat
   divideflat,*imageout[0],flatfits[1024:2047,0:1023],minval=0.001
   divideflat,*imageout[1],flatfits[0:1023,0:1023],minval=0.001
   divideflat,*imageout[2],flatfits[1024:2047,1024:2047],minval=0.001
   divideflat,*imageout[3],flatfits[0:1023,1024:2047],minval=0.001

   FOR i=0,3 DO BEGIN
      
      *imageout[i] = *imageout[i]/str[i].extime[whereScience]
      *imageout[i] = *imageout[i]*str[i].gain[whereScience]

   ENDFOR

   IF KEYWORD_SET(skyfit) THEN BEGIN

      ;;remove luminosity pattern
      SPLOG,"Skyfit does not work! :( "

   ENDIF
   
   ;;update WCS information account for the satasec (needed?)
   ;;I'm leaving this out for now since I don't know what it is.
 ;  FOR i=0,3 DO BEGIN
 ;
 ;     crpix1 = FXPAR(*headers[i],'CRPIX1')
 ;     crpix1 = crpix1 - 51
 ;     SXADDPAR,*headers[i],'CRPIX1',crpix1
 ;
 ;  ENDFOR
   
   automated_imagestitch,imageout,stitchedScience

   mainHeader = headfits(STRCOMPRESS(str[0].name[whereScience],/REMOVE_ALL))
   exHeader = HEADFITS(STRCOMPRESS(str[0].name[whereScience],/REMOVE_ALL),exten=1)
   

   ;;rewrite and modify extension headers.
   ;IF (STRCOMPRESS(FXPAR(mainHeader,'OBJECT'),/REMOVE_ALL) NE 'BL_Lac') THEN BEGIN

    ;  FXADDPAR,mainHeader,'OBJECT','BL_Lac'

   ;ENDIF
   

  ; FXADDPAR,mainHeader,'CD1_1',FXPAR(exHeader,'CD1_1')
  ; FXADDPAR,mainHeader,'CD1_2',FXPAR(exHeader,'CD1_2')
  ; FXADDPAR,mainHeader,'CD2_1',FXPAR(exHeader,'CD2_1')
  ; FXADDPAR,mainHeader,'CD2_2',FXPAR(exHeader,'CD2_2')
  ; FXADDPAR,mainHeader,'CRVAL1',FXPAR(exHeader,'CRVAL1')
  ; FXADDPAR,mainHeader,'CRVAL2',FXPAR(exHeader,'CRVAL2')
  ; FXADDPAR,mainHeader,'CRPIX1',1024.25
  ; FXADDPAR,mainHeader,'CRPIX2',1024.25

  ; FXADDPAR,mainHeader,'CD1_1',-0.00018752818859
  ; FXADDPAR,mainHeader,'CD1_2',1.4817112956E-06
  ; FXADDPAR,mainHeader,'CD2_1',1.445868427E-06
  ; FXADDPAR,mainHeader,'CD2_2',0.000187444213011
  ; FXADDPAR,mainHeader,'CRVAL1',330.649087254
  ; FXADDPAR,mainHeader,'CRVAL2',42.244975269
  ; FXADDPAR,mainHeader,'CRPIX1',1106.87126133
  ; FXADDPAR,mainHeader,'CRPIX2',961.874938629

   ;FXADDPAR,mainHeader,'CTYPE1',FXPAR(exHeader,'CTYPE1')
   ;FXADDPAR,mainHeader,'CTYPE2',FXPAR(exHeader,'CTYPE2')
   FXADDPAR,mainHeader,'BITPIX',-32
   FXADDPAR,mainHeader,'NAXIS',2,' number of bits per data pixel'
   FXADDPAR,mainHeader,'NAXIS1',2048,' length of data axis 1',AFTER='NAXIS'
   FXADDPAR,mainHeader,'NAXIS2',2048,' length of data axis 2',AFTER='NAXIS1'
   FXADDPAR,mainHeader,'BZERO',0.00000,' Set by MRD_SCALE'
   FXADDPAR,mainHeader,'BSCALE',1.00000,' Set by MRD_SCALE'
   FXADDPAR,mainHeader,'GAIN1',str[0].gain[0], ' First amp gain'
   FXADDPAR,mainHeader,'GAIN2',str[1].gain[1], ' Scendon amp gain',AFTER='GAIN1'
   FXADDPAR,mainHeader,'GAIN3',str[2].gain[2], ' Third amp gain',AFTER='GAIN2'
   FXADDPAR,mainHeader,'GAIN4',str[3].gain[3], ' Fourth amp gain',AFTER='GAIN3'
   FXADDPAR,mainHeader,'RN1',str[0].RN[0],' First amp readnoise',AFTER='GAIN4'
   FXADDPAR,mainHeader,'RN2',str[1].RN[1],' Second amp readnoise',AFTER='RN1'
   FXADDPAR,mainHeader,'RN3',str[2].RN[2],' Third amp readnoise',AFTER='RN2'
   FXADDPAR,mainHeader,'RN4',str[3].RN[3],' Fourth amp readnoise',AFTER='RN3'

   posit = STRPOS(str[0].NAME[whereScience],'.fit')
   substring = STRMID(str[0].NAME[whereScience],0,posit)

   MWRFITS,stitchedScience,substring + '_redux.fits',mainheader,/CREATE,/SILENT




   IF KEYWORD_SET(gzip) THEN BEGIN

      ;;gzip image
      COMMAND='gzip ' + substring + '_redux.fits'
      SPAWN,command
      SPLOG,'Done with ' + substring + '_redux.fits.gz'
   ENDIF ELSE SPLOG,'Done with ' + substring + '_redux.fits'

   
  

END
