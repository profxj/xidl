;+
; NAME:
;   extract_1d_oldwave - extract_1d modified for 2D lambda array existing
;
; PURPOSE:
;   extracts 1-d spectrum from data, using inverse weighting.
;
; CALLING SEQUENCE:
;   spec = extract_1d_oldwave(slit,r1,r2, [ivar = ])
; 
; INPUTS:
;   slit -- structure of spSlit or slit outputs
;   r1   -- lower row for extraction
;   r2   -- upper row for extraction
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;   spec  -- structure of 1d extracted spectrum
;           spec.flux --flux
;           spec.lambda -- lambda
;           spec.ivar   -- ivar
; OPTIONAL OUTPUTS:
;   ivar=ivar -- resulting inverse variance
;
; RESTRICTIONS:
;   
; EXAMPLES:
;
; COMMENTS:
;  lambda is rectified row by row with full pixel shifts to align
;  spectra before adding.  Wraps in Lambda is so doing.
;
; REVISION HISTORY:
;   MD 04Jul2002
;----------------------------------------------------------------------
function extract_1d_oldwave, slit, r1, r2

  spec = slit.flux(*, 0)*0.
  ivar = spec*0.
  npix = n_elements(spec)
  xx = findgen(npix)/npix*2 -1
;  lambda2d = lambda_eval(slit)
  lambda2d = slit.lambda0#(fltarr(n_elements(slit.slitfn))+1)+slit.lambda
  lambda = lambda2d[*, round((r1+r2)/2.)] ;lambda on the central row of object
  dldx = (lambda[npix-1]-lambda[0])/npix ;Ang/pixel
;  dxdp = slit.tiltx[0]/dldx*lambda[npix/2] ;slope of constant lambda/vert pix
  dxdp =(lambda2d[npix/2., n_elements(slit.slitfn)-1]-lambda2d[npix/2., 0])/n_elements(slit.slitfn)


  for i=r1, r2 do begin
    cshift = round((i-float(r1+r2)/2.)*dxdp)
    spec = spec + shift(slit.flux[*, i]*slit.ivar[*, i], cshift)
    ivar = ivar + shift(slit.ivar[*, i], cshift)
  endfor
  spec = spec/ivar

  ss =  {spec: spec,  lambda: lambda, ivar:ivar}

return,  ss
end


