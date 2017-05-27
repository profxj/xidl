;+
; NAME:
;   extract_1d
;
; PURPOSE:
;   simple extraction of 1-d spectrum from data using inverse
;   weighting and tilt of lines. This is a flux conserving extraction,
;   not an optimally weighted average
;
; CALLING SEQUENCE:
;   spec = extract_1d(slit,r1,r2, [r3,r4] [ivar = ])
; 
; INPUTS:
;   slit -- structure of spSlit or slit outputs
;   r1   -- lower row for extraction
;   r2   -- upper row for extraction
;
; OPTIONAL INPUTS:
;   r3   -- lower bound for sky subtraction region
;   r4   -- upper bound for sky subtraction region
;	
; KEYWORDS:
;
; OUTPUTS:
;   spec  -- structure of 1d extracted spectrum
;           spec.flux --flux
;           spec.lambda -- lambda
;           spec.ivar   -- ivar
;           spec.sky    -- sky spectrum
; OPTIONAL OUTPUTS:
;   ivar=ivar -- resulting inverse variance
;
; RESTRICTIONS:
;   
; EXAMPLES:
;
; COMMENTS:
;  lambda is rectified row by row with full pixel shifts to align
;  spectra before adding.  Wraps in Lambda in so doing.
;
; REVISION HISTORY:
;   MD 04Jul2002
;       08Aug2002 -- update Lambda handling and tilt, sky spectrum returned
;       10Oct2002 -- correct ivar extraction to be flux conserving
;----------------------------------------------------------------------
function extract_1d, slit, r1, r2, r3, r4

  nsky = 0
  spec = slit.flux(*, 0)*0.
  var2d = slit.ivar*0. ;initial variance array
  ivar2d = slit.ivar
  sky = spec
  var = spec
  count = spec 
  ivar = spec*0.
  ivarsky = ivar
  npix = n_elements(spec)
  nrows = (size(slit.flux, /dimens))[1]
  xx = findgen(npix)/npix*2 -1
  lambda2d = lambda_eval(slit)
  lambda = lambda2d[*, round((r1+r2)/2.)] ;lambda on the central row of object
  dldx = (lambda[npix-1]-lambda[0])/npix ;Ang/pixel
;  dxdp = slit.tiltx[0]/dldx*lambda[npix/2] ;slope of constant lambda/vert pix
  tiltx = (lambda2d(npix/2, nrows-1) - lambda2d(npix/2, 0))/nrows
  dxdp = tiltx/dldx ;slope of constant lambda/vert pix

  good = where(slit.ivar gt 0) 
  var2d[good] = 1./ivar2d[good]

  for i=r1, r2 do begin
    cshift = round((i-float(r1+r2)/2.)*dxdp)
;    spec = spec + shift(slit.flux[*, i]*slit.ivar[*, i], cshift)
; want a flux preserving integral here, not an average
    spec = spec + shift(slit.flux[*, i], cshift)
;    ivar = ivar + shift(slit.ivar[*, i], cshift)
    var = var + shift(var2d[*, i], cshift)
    count = count + shift(var2d[*, i] gt 0, cshift)

  endfor
; now give average
    good1 = where(count gt 0, ngood)
    spec[good1] = spec[good1]/count[good1]*(r2-r1+1)
    ivar[ngood] = 0.
    ivar[good1] = 1./var[good1]*(count[good1]/(r2-r1+1))^2

;  nzero = where(ivar gt 0) ;normalize regions not given 0 weight
;  spec[nzero] = spec[nzero]/ivar[nzero]
  

  if N_PARAMS() GT 3 then begin ;process sky if requested
    r3t = r3

    for i=r3, r2-1 do begin
      cshift = round((i-float(r1+r2)/2.)*dxdp)
      sky = sky + shift(slit.flux[*, i]*slit.ivar[*, i], cshift)
      ivarsky = ivarsky + shift(slit.ivar[*, i], cshift)     
      nsky = nsky+1
    endfor
  endif

  if N_PARAMS() GT 4 then begin ;process 2nd sky if requested
    r4t = r4

    for i=r2, r4-1 do begin
      cshift = round((i-float(r1+r2)/2.)*dxdp)
      sky = sky + shift(slit.flux[*, i]*slit.ivar[*, i], cshift)
      ivarsky = ivarsky + shift(slit.ivar[*, i], cshift)   
      nsky = nsky+1
    endfor
  endif


  if N_PARAMS() gt 3 then begin
    nzero = where(ivarsky gt 0)
    sky[nzero] = sky[nzero]/ivarsky[nzero] ; average sky
    spec = spec - (r2-r1)*sky ;sky subtract if requested
    ivar = ivar + ivarsky/(r2-r1)
  endif

; a sum, not an average!
;  spec = spec*(r2-r1) ;get SUMMED sky
;  ivar = ivar/(r2-r1)^2 ;correct ivar 

;stop

  ss =  {spec: spec,  lambda: lambda, ivar:ivar,  sky:sky}

return,  ss
end









