;+ 
; NAME:
; sdss_hizmedian
;    Version 1.0
;
; PURPOSE:
;    Calculates a median QSO spectrum for a set of high z SDSS QSOs
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

PRO sdss_hizmedian, DR3_fits, new_template, szstart=SZSTART, szend=SZEND, $
                        dr3path=DR3_PATH, WAVEMIN=wavemin, WAVEMAX=wavemax
                       
;set the DR3 path
   IF NOT keyword_set(dr3path) THEN dr3path=getenv('SDSSPATH')+'/DR3_QSO/'

;open the dr3 summary fits file
   dr3=xmrdfits(dr3path+DR3_fits,1,/silent)

;determine the number of qsos to search
   IF NOT keyword_set (szstart) THEN szstart = 3.2
   IF NOT keyword_set (szend) THEN szend = max(dr3.z)
   
   range = where(dr3.z GE szstart AND dr3.z LE szend,nqsos)
   dr3=dr3[range]
;   nqsos = n_elements(range)

   IF NOT keyword_set(wavemin) then loglam_min = 2.9000  $ ; lambda=795 Ang
   else loglam_min = alog10(wavemin)

   IF NOT keyword_set(wavemax) then loglam_max = 3.2000 $  ; lambda=1585 Ang
   else loglam_max = alog10(wavemax)


   npix = long((loglam_max - loglam_min)/0.0001) + 1

   full_spec_array = fltarr(npix, nqsos)
   full_weight_array = fltarr(npix, nqsos)
   
   old_template=xmrdfits('template.fits',/silent,hdr)
   QSOtempwav0 = 2.8
   

   ;; Loop
   FOR qq=0,nqsos-1 do begin
      ;make string verisons of the plate,mjd,fiber
       splate =string(dr3[qq].plate, format='(i4.4)')
       smjd  = string(dr3[qq].mjd, format='(i5.5)')
       sfiber= string(dr3[qq].fiberid, format='(i3.3)')

      ;create the file name 
       filename = strcompress('spSpec-'+smjd+'-'+splate+'-'+sfiber+'.fit.gz', $
                              /remove_all)
       filedir = strcompress(dr3path+'spectro/1d_23/'+splate+'/1d/', $
                             /remove_all)
       file=filedir+filename
       

       print,'Now doing qso '+filename+' (#'+ $
         strcompress(string(qq),/remove_all)+' of '+ $
         strcompress(string(nqsos-1),/remove_all)+')'

       ; Detilt spectrum, and return loglam, normalized flux

       spspec=mrdfits(file,0,hdr,/silent)
       spec = spspec[*,0]
       specerr = spspec[*,2]
       speczem = sxpar(hdr,"Z") 
       if (abs(speczem - dr3[qq].z) GT 0.01) THEN CONTINUE
       specwav0 = sxpar(hdr,"COEFF0")
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
 
 mwrfits, final_template, new_template, /create

   return
END

