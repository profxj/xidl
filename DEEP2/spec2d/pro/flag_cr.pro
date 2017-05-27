;+
; NAME:
;   flag_cr
;
; PURPOSE:
;   Flag regions that appear to be cosmic rays by setting their
;   inverse variance to zero.  To be run before fitting bsplines to
;   DEEP/DEIMOS data.
;
; CALLING SEQUENCE:
;   flag_cr, slit, newivar
;   
; INPUTS:
;   slit     - A DEEP survey spSlit structure
;   
; OUTPUTS:
;   newivar  - The new inverse variance array, with cosmic rays
;              flagged.
;   crmask   - 1 on CRs, 0 elsewhere
; 
; COMMENTS:
;   Written to keep bsplines from incorrectly fitting to cosmic rays.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Written by BFG, jul02
;-
 
pro flag_cr,  slit, newivar, crmask, medflux=medflux
  
  flux = slit.flux
  newivar =  slit.ivar
  ncols =  n_elements(flux[*, 0]) 
  nrows =  n_elements(flux[0, *])
  xidx=lindgen(ncols,nrows)
  yidx=xidx / ncols
  xidx = xidx MOD ncols

  notatedge=yidx ne 0 AND yidx ne (nrows-1) AND $
    xidx ne 0 AND xidx ne (ncols-1)


;;; DO A SIMPLE COSMIC RAY REJECTION by setting ivar=0 at cosmic rays.
;;; FIRST SMOOTH THE FLUX USING A MEDIAN FILTER.
  if n_elements(medflux) eq 0 then $
    medflux = djs_MEDIAN(flux, width=5, boundary='reflect')
  for i=0,4 do medflux[i,*]=median(flux[i,*])<medflux[i,*]
  for i=ncols-5,ncols-1 do medflux[i,*]=median(flux[i,*])<medflux[i,*]

  crmask = medflux*0b

  ivar = slit.ivar

  wh = where(ivar EQ 0., whct)
  whok=where(ivar NE 0.,okct)
  IF (whct GT 0 AND okct GT 0) THEN ivar[wh] = .3*mean(ivar[whok]) $
    else if whct GT 0 then ivar[wh]=1.E-20


;---flag outliers above 20 sigma from the local median
  whcr = WHERE( (flux-medflux GT 20./SQRT(ivar)) AND notatedge, whct)
  if whct gt 0 then begin
    ;print,  'flagged '+string(whct)+ ' CR pixels.'
    newivar[whcr] = 0.
    crmask[whcr] = 1b
  endif

  testcr=dilate(crmask,bytarr(7,7)+1)
  whcr = WHERE( (flux-medflux GT 10./SQRT(ivar)) $
                AND testcr AND notatedge, whct)

  if whct gt 0 then begin
    ;print,  'flagged '+string(whct)+ ' CR pixels.'
    newivar[whcr] = 0.
    crmask[whcr] = 1b
  endif

  testcr=dilate(crmask,bytarr(3,3)+1)
  whcr = WHERE( (flux-medflux GT 5./SQRT(ivar)) AND testcr, whct)

  if whct gt 0 then begin
    ;print,  'flagged '+string(whct)+ ' CR pixels.'
    newivar[whcr] = 0.
    crmask[whcr] = 1b
  endif
  testcr=dilate(crmask,bytarr(3,3)+1)
  whcr = WHERE( (flux-medflux GT 5./SQRT(ivar)) AND testcr, whct)

  if whct gt 0 then begin
    ;print,  'flagged '+string(whct)+ ' CR pixels.'
    newivar[whcr] = 0.
    crmask[whcr] = 1b
  endif


return

end



