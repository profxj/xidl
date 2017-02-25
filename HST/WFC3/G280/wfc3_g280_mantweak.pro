;+ 
; NAME:
;  wfc3_g280_mantweak
;
; PURPOSE:
;   Manually tweak the beams to exclude certain regions of contamination
;
; CALLING SEQUENCE:
;   
;  wfc3_g280_mantweak, fil, tweak, FITSDIR=fitsdir, QADI=qadir, chk=chk
;
; INPUTS:
;  fil -- the file to read in containing the fits file to tweak and
;         the corrections to make. The file should be of the form:
;         test.fits  23-45,56-78,123-156 This will tell the file to
;         ignore columns 23 through 45, 56 though 78 and columns 123
;         though 156.
;         ALTERNATIVELY; if two inputs are set this input is
;         taken to be the fits file to reduce
;  tweak -- If two inputs are set this is a string file containing the
;           tweaks to be made (in the format described above)
;
; RETURNS:
;
; OUTPUTS:
;   Updates the FITS file with the described tweaks
;
; OPTIONAL KEYWORDS:
;  FITSDIR= -- Directory containing the FITS file
;  QADIR= -- Directory to write the QA file to
;  CHK= -- plots some checks to see if the tweak is good
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  wfc3_g280_mantweak, fil, tweak, FITSDIR=fitsdir, QADIR=qadir, CHK=chk
;
; PROCEDURES CALLED:
;  cgsetintersection
;  wfc3_g280_qa_final
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

pro wfc3_g280_mantweak, fil, tweak, FITSDIR=fitsdir, QADIR=qadir, CHK=chk

  if not keyword_set(FITSDIR) then fitsdir='Fits/'
  if not keyword_set(QADIR) then qadir='QA/'

  ;; assume that a file is input:
  if (n_params() LT 2) then $
     readcol, fil, fits,tweak, format='A,A', delim=' ', comment='#' $
  else $
     fits=fil

  for ii=0, n_elements(fits)-1 do begin
     
     ;; read in the file
     wfc3_g280_strct=mrdfits(fitsdir+fits[ii],1,h)

     ;; first clear the regions of previous bad regions
     wfc3_g280_strct.badpixa=dblarr(1000)
     wfc3_g280_strct.badpixc=dblarr(1000)

     ;; save the original structure in this variable
     owfc3_g280_strct=wfc3_g280_strct
     
     ;; create interpolated arrays
     fluxa=interpol(wfc3_g280_strct.fluxa(0:wfc3_g280_strct.cnta-1L),$
                    wfc3_g280_strct.wavea(0:wfc3_g280_strct.cnta-1L),wfc3_g280_strct.wave)
     siga=interpol(wfc3_g280_strct.flux_siga(0:wfc3_g280_strct.cnta-1L),$
                   wfc3_g280_strct.wavea(0:wfc3_g280_strct.cnta-1L),wfc3_g280_strct.wave)
     skya=interpol(wfc3_g280_strct.skya(0:wfc3_g280_strct.cnta-1L),$
                   wfc3_g280_strct.wavea(0:wfc3_g280_strct.cnta-1L),wfc3_g280_strct.wave)
     fluxc=interpol(wfc3_g280_strct.fluxc(0:wfc3_g280_strct.cntc-1L),$
                    wfc3_g280_strct.wavec(0:wfc3_g280_strct.cntc-1L),wfc3_g280_strct.wave)
     sigc=interpol(wfc3_g280_strct.flux_sigc(0:wfc3_g280_strct.cntc-1L),$
                   wfc3_g280_strct.wavec(0:wfc3_g280_strct.cntc-1L),wfc3_g280_strct.wave)
     skyc=interpol(wfc3_g280_strct.skyc(0:wfc3_g280_strct.cnta-1L),$
                   wfc3_g280_strct.wavea(0:wfc3_g280_strct.cnta-1L),wfc3_g280_strct.wave)
     
     ;; create the array of bad pixels
     badpixa=[0L]
     badpixc=[0L]
     split1=strsplit(tweak[ii],',',/extract)
     for jj=0L,n_elements(split1)-1 do begin
        split2=long(strsplit(split1(jj), '-', /extract))
        case n_elements(split2) of
           1: badreg=split2
           2: badreg=lindgen(split2[1]-split2[0]+1)+split2[0]
           else: stop
        endcase

        ;; find the intersection of this bad region and the traces
        tmpa=cgsetintersection(wfc3_g280_strct.trace_xa,badreg,pos=tbadpixa,suc=suca)
        if suca then begin
           badpixa=[badpixa,tbadpixa]
           ;; find the wavelength range of the bad region
           minwva=min(wfc3_g280_strct.wavea(tbadpixa),max=maxwva)
           tmpa=where(wfc3_g280_strct.wave ge minwva and wfc3_g280_strct.wave le maxwva, na)
           if na gt 0 then begin
              wfc3_g280_strct.flux(tmpa)=fluxc(tmpa)
              wfc3_g280_strct.sig(tmpa)=sigc(tmpa)
              wfc3_g280_strct.sky(tmpa)=skyc(tmpa)
           endif
        endif

        tmpc=cgsetintersection(wfc3_g280_strct.trace_xc,badreg,pos=tbadpixc,suc=succ)
        if succ then begin
           badpixc=[badpixc,tbadpixc]
           ;; find the wavelength range of the bad region
           minwvc=min(wfc3_g280_strct.wavec(tbadpixc),max=maxwvc)
           tmpc=where(wfc3_g280_strct.wave ge minwvc and wfc3_g280_strct.wave le maxwvc, nc)
           if nc gt 0 then begin
              wfc3_g280_strct.flux(tmpc)=fluxa(tmpc)
              wfc3_g280_strct.sig(tmpc)=siga(tmpc)
              wfc3_g280_strct.sky(tmpc)=skya(tmpc)
           endif
        endif

        if keyword_set(chk) then begin
           off=200
           if suca then xr=[minwva-off,maxwva+off]
           if succ then xr=[minwvc-off,maxwvc+off]
           if ~suca and ~succ then continue
           cgplot, wfc3_g280_strct.wave,  wfc3_g280_strct.flux, $
                   color='black', thick=2, psym=10, xthick=2, $
                   ythick=2, xr=xr, tit=wfc3_g280_strct.name
           cgplot, wfc3_g280_strct.wavea, wfc3_g280_strct.fluxa, /ov, $
                   color='red', line=2, psym=10, thick=2
           cgplot, wfc3_g280_strct.wavec, wfc3_g280_strct.fluxc, /ov, $
                   color='blue', line=2, psym=10, thick=2
           cgplot, wfc3_g280_strct.wavea, wfc3_g280_strct.flux_siga, /ov, $
                   color='red', line=1, psym=10, thick=2
           cgplot, wfc3_g280_strct.wavec, wfc3_g280_strct.flux_sigc, /ov, $
                   color='blue', line=1, psym=10, thick=2
           cgplot, wfc3_g280_strct.wave, wfc3_g280_strct.sig, /ov, $
                   color='purple', line=2, psym=10, thick=2
           stop
        endif
        
     endfor
     
     ;; save the badpixel regions to the original structure
     owfc3_g280_strct.badpixa=badpixa
     owfc3_g280_strct.badpixc=badpixc
     fitsname=fitsdir+owfc3_g280_strct.name+'_g280I.fits'
     mwrfits, owfc3_g280_strct, fitsname, /create

     ;; redo the QA plot
     wfc3_g280_strct.badpixa=badpixa
     wfc3_g280_strct.badpixc=badpixc
     wfc3_g280_qa_final, wfc3_g280_strct, 0, QADIR=qadir
     
     ;; save final image
     wfc3_g280_mkfits, wfc3_g280_strct, FITSDIR=fitsdir
     
  endfor
  
end
