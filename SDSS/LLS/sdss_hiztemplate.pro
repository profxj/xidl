;+ 
; NAME:
; sdss_hiztemplate
;    Version 1.0
;
; PURPOSE:
;    Calculates a template QSO spectrum for a set of high z SDSS QSOs
;   I am not sure how this is different than sdss_hizmedian, but I
;   suspect this is the better routine.
;
; CALLING SEQUENCE:
;  sdss_hizmedian, DR3_fits, new_template, szstart=, szend=
;
; INPUTS:
; DR3_fits -- Filename of qso file
;
; RETURNS:
;   
; OUTPUTS:
;  new_template -- Filename of median spectrum (FITS file)
;
; OPTIONAL KEYWORDS:
;  DR3PATH -- Path to the DR3 data
;  szstart  -- Starting redshift for the stack
;  szend    -- Ending redshift for the stack
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   07-May-2002 Written by JO
;-
;------------------------------------------------------------------------------
;sdss_hiztemplate,'dr3_qso.fits','foo.fits', szstart=3.5, szend=3.53
;sdss_hiztemplate,'dr6_qso.fits','hiztemplate_32_34.fits', szstart=3.2, szend=3.4

PRO sdss_hiztemplate, DR_fits, new_template, szstart=SZSTART, szend=SZEND, $
  drpath=DRPATH, WAVEMIN=wavemin, WAVEMAX=wavemax, $
  smooth=SMOOTH, MOCK=mock

;set the DR3 path
   IF NOT keyword_set(drpath) THEN drpath=getenv('SDSSPATH')+'/DR7_QSO/'
   IF NOT keyword_set (szstart) THEN szstart = 3.2
   IF NOT keyword_set(wavemin) then loglam_min = 2.9000  $ ; lambda=795 Ang
   else loglam_min = alog10(wavemin)
   
   IF NOT keyword_set(wavemax) then loglam_max = 3.2000 $ ; lambda=1585 Ang
   else loglam_max = alog10(wavemax)

   ;;open the dr3 summary fits file
   dr=xmrdfits(drpath+DR_fits,1,/silent)

   if tag_exist(dr,'ZEM') then zem = dr.zem else zem = dr.z

   ;;determine the number of qsos to search
   IF NOT keyword_set (szend) THEN szend = max(zem)
       
   range = where(zem GE szstart AND zem LT szend,nqsos)
   if nqsos EQ 0 then begin
       print, 'sdss_hiztemlate:  No quasars at your redshift of interest'
       print, 'sdss_hiztemlate:  Returning'
       return
   endif
   dr=dr[range]

;  ;; nqsos = n_elements(range)
   npix = long((loglam_max - loglam_min)/0.0001) + 1
   
   full_spec_array = fltarr(npix, nqsos)
   full_weight_array = fltarr(npix, nqsos)
   
   old_template=xmrdfits(getenv('XIDL_DIR')+'/SDSS/LLS/Full_SDSS_template.fits', $
                         /silent,hdr)
   QSOtempwav0 = 2.8
   
   ;; Loop on quasars
   FOR qq=0,nqsos-1 do begin

       if not keyword_set(MOCK) then begin
           ;;make string verisons of the plate,mjd,fiber
           splate =string(dr[qq].plate, format='(i4.4)')
           smjd  = string(dr[qq].mjd, format='(i5.5)')
           sfiber= string(dr[qq].fiberid, format='(i3.3)')
           
           ;;create the file name 
           filename = strcompress('spSpec-'+smjd+'-'+splate+'-'+sfiber+'.fit.gz', $
                                  /remove_all)
           filedir = strcompress(drpath+'spectro/1d_26/'+splate+'/1d/', $
                                 /remove_all)
           file=filedir+filename
           
           
           print,'Now doing qso '+filename+' (#'+ $
                 strcompress(string(qq),/remove_all)+' of '+ $
                 strcompress(string(nqsos-1),/remove_all)+')'
           
           ;; Detilt spectrum, and return loglam, normalized flux
           
           spspec=mrdfits(file,0,hdr,/silent)
           spec = spspec[*,0]
           specerr = spspec[*,2]
           speczem = sxpar(hdr,"Z") 
           if (abs(speczem - dr[qq].z) GT 0.01) THEN CONTINUE
           specwav0 = sxpar(hdr,"COEFF0")
       endif else begin ; MOCK
           spec=dr[qq].spec
           specerr=dr[qq].sig
           speczem = dr[qq].zem
           specwav0 = dr[qq].w0
       endelse

       logwav = findgen(n_elements(spec))*0.0001 + specwav0
       dwave=0.0001
       pixmin = long((specwav0 - QSOtempwav0 - alog10(1. + speczem))/dwave)
       tempsec = extrac(old_template,pixmin,n_elements(spec))


       good = where(logwav GT alog10((1.0 + speczem)*1215.67) $
                     AND tempsec GT 0 and specerr GT 0)
       res  =ladfit(logwav[good],spec[good]/tempsec[good])
       modelfit = poly(logwav,res)

       normflux = spec / (modelfit + (modelfit EQ 0)) 

       pixmin2 = long((loglam_min - specwav0 + alog10(1. + speczem))/dwave)
       full_spec_array[*,qq] = extrac(normflux, pixmin2, npix)
       full_weight_array[*,qq] = extrac((specerr GT 0.) * (modelfit GT 0),pixmin2,npix)
   ENDFOR

   final_template = fltarr(npix)

   for ii=0,npix-1 do begin
       good2 = where(full_weight_array[ii,*] GT 0., ngood2)
       if ngood2 LT 10 then begin
           print, 'Less than 10 good pixels at ', ii, ngood2
           continue
       endif
       
       djs_iterstat, full_spec_array[ii, good2], median=med
       
       final_template[ii] = med
   endfor
 
   smooth_template = smooth(final_template,50)
   
   
;set the template to the value of the flux at 920 angstroms for all 
;wavelengths below 920

   wav=10^(2.9 + findgen(3000)*0.0001)
   wavmin=where(wav LT 920.0,minindx)
   smooth_template[0:minindx-1]=smooth_template[minindx]
   final_template[0:minindx-1]=final_template[minindx]

   ;extend array down
   tmparr = fltarr(4000)
   tmparr2 = fltarr(4000)

   tmparr[1000:3999]=smooth_template
   tmparr[0:999]=smooth_template[0]
   tmparr2[1000:3999]=final_template
   tmparr2[0:999]=final_template[0]

   


   IF NOT keyword_set(SMOOTH) THEN BEGIN
       mwrfits, tmparr2, new_template, /create
   ENDIF ELSE BEGIN
       mwrfits, tmparr, new_template, /create
   ENDELSE

   wav2=10^(2.8 + findgen(4000)*0.0001)
   plot,wav2,tmparr2
   oplot,wav2,tmparr,color=255


;stop
   return
END

