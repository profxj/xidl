;+ 
; NAME:
; sdss_qsotempl
;
; PURPOSE:
;    Extrapoles a continuum for a SDSS QSO in the Lya forest based on
;    a template, a fit to the QSO slope and the known average
;    absorption by the Lya forest. 
;
; CALLING SEQUENCE:
;   sdss_qsotempl
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  ORIG_ZEM =  Use the original redshifs from SDSS (in this file)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Feb-2005 Written by JO
;   14-Oct-2006 Modified by JXP
;-
;------------------------------------------------------------------------------
; help, sdss_findlls(template, 2.8, 'sdss3.fit.gz')
Function sdss_qsotempl,template,QSOtempwav0, logwav, spec, specerr, speczem, $
  TEMPFIL=tempfil, TMPWV=tmpwv, RES4=res4, $
  EPX=epx, FPX=fpx, SPECFIL=specfil

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
            'templ_spec = sdss_qsotempl(template, QSOtempwav0, [logwav,spec,speczem]) [v1.0]'
      return, -1
  endif

  if keyword_set(TEMPFIL) then template =xmrdfits(tempfil,/silent)
  
  tlogwav = findgen(n_elements(template))*0.0001 + QSOtempwav0
;;open SDSS spectrum
  if keyword_set(SPECFIL) then begin
      spspec=xmrdfits(SDSSspec,0,hdr,/silent)
      spec = spspec[*,0]
      specerr = spspec[*,2]
      speczem = sxpar(hdr,"Z")
      specwav0 = sxpar(hdr,"COEFF0")
      logwav = findgen(n_elements(spec))*0.0001 + specwav0
      speczem = sxpar(hdr,"Z") 
  endif
  dwave=0.0001


;;;redshift the template
  specwav0 = logwav[0]
  tzlogwav = alog10((1.0 + speczem) * (10.0^tlogwav))
  
;;;pick the template minimum pixel to extract from 
  pixmin = long((specwav0 - QSOtempwav0 - alog10(1. + speczem))/dwave)
  
;;;You can unweight the forest by dividing the forest by 0.7 if you include
;;;the forest in the regions which effect the slope of the template
  tempsec = extrac(template,pixmin,n_elements(spec))
;lyafrange=where(10^logwav LT (1.0 + speczem)*1215.67)
;tempsec(lyafrange)=tempsec(lyafrange) / 0.7

;;;Find the good pixels
  tempmask = (tempsec GT 0 AND specerr GT 0)
  good=where(tempmask)
;res=ladfit(logwav[good],spec[good]/tempsec[good])
;modelfit=poly(logwav,res)

;;;This part of the code slaps the slope on the SDSS template
;;;Here you can set the model matching range to either be where there are good
;;;pixels OR where there are good pixels AND wave > lya-emission

;good2=where(logwav GT llswavmax AND tempsec GT 0 and specerr GT 0)
;good2=where(logwav GT alog10((1.0 + speczem)*1250.67) and $
;            logwav LT alog10((1.0 + speczem)*1520.0) and $
;            tempsec GT 0 and specerr GT 0)
;res2=ladfit(logwav[good2],spec[good2]/tempsec[good2])
; 
;good3=where(logwav GT alog10((1.0 + speczem)*1050.0) AND logwav LT $
;                            alog10((1.0 + speczem)*1170.0) $
;            AND tempsec GT 0 and specerr GT 0)
;res3=ladfit(logwav[good3],smooth(spec[good3],15)/tempsec[good3])

good4=where(logwav GT alog10((1.0 + speczem)*950.0) and $
            logwav LT alog10((1.0 + speczem)*1800.0) and $
            tempsec GT 0 and specerr GT 0)

res4=ladfit(logwav[good4],smooth(spec[good4],15)/tempsec[good4])

modelfit2=poly(logwav,res4)

fpx = round((specwav0 - 3.58)*1e4)
IF (fpx LT 0) THEN fpx=0
lend = alog10((1.0 + speczem)*1220.0)
epx = round((lend - 3.58)*1e4)

tempspec = (modelfit2*tempsec)[0:epx-fpx]
tmpwv = 10^logwav[0:epx-fpx]


medrange=where((10^logwav[0:epx-fpx])/(1.+speczem) GT 1050.0 AND $
               (10^logwav[0:epx-fpx])/(1.+speczem) LT 1170.0)

;; new stuff, 2/28/07, finds the mean forest redshift between 1050 and 1170
;; then looks up the kirkman 05 DA, and divides that out
;; NEW new stuff, Oct 08, no longer taking out DA, so commenting out
;; corrections

;medforz = median(10^logwav[medrange]/1215.67 -1.0)
;kirkman_da = 0.0062*(1.0 + medforz)^(2.75)
;dacorr = (1.0 - kirkman_da)

;tempspec = tempspec / dacorr

  return,tempspec
END
