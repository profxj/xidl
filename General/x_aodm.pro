;+ 
; NAME:
; x_aodm
;   Version 1.0
;
; PURPOSE:
;    Calculates an AODM column density
;
; CALLING SEQUENCE:
;   
;   x_aodm, wave, fx, sig, rwave, colm, sig
;
; INPUTS:
;   wave - Wavelength array
;   fx - Data
;   sig - Error array
;   rwave - Rest wavelength
;
; RETURNS:
;
; OUTPUTS:
;   colm  - Column density
;   sig  - Error
;
; OPTIONAL KEYWORDS:
;   LOG - Return log answers
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   colm = x_aodm( dat[203:209], 1808.0126d)
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Dec-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_aodm, wave, fx, sig, rwave, colm, sig_colm, LOG=log, TAU=tau

  if (N_params() LT 4) then begin 
    print,'Syntax - ' + $
             'x_aodm, wave, fx, sig, rwave, colm, sig_colm, /LOG, TAU=  [V1.0]'
    return
  endif 

  ; Error catching
  if n_elements(fx) NE n_elements(sig) then $
    message, 'x_aodm: n(fx) NE n(sig) '
  npix = n_elements(fx)

  ; Get the fvalue
  getfnam, rwave, fval

  ; Constant
  cst = (10.d^14.5761)/(fval*rwave)

  ; Deal with saturated pixels
  sat = where(fx LT sig/5. OR fx LT 0.05, nsat, $
              complement=unsat, ncomplement=nunsat)
  
  ; Set nndt
  nndt = dblarr(npix)
  if nsat NE 0 then nndt[sat] = alog(5.d/sig[sat])*cst
  if nunsat NE 0 then nndt[unsat] = alog(1.d/fx[unsat])*cst

  ; Set delv
  cent = (wave[0]+wave[npix-1])/2.d
  spl = x_constants('spl')

  velo = (wave-cent)*spl/cent
  delv = dblarr(npix)
  for i=0L,npix-2 do delv[i] = velo[i+1]-velo[i]
  delv[npix-1] = delv[npix-2]

  ; Sum 

  ntot = total( nndt*delv, /double )
  tvar = total( (delv*cst*sig/fx)^2, /double )
      
  ; Log as requested
  if keyword_set( LOG ) then begin
      if ntot GT 0. then colm = alog10( ntot ) else colm = -9.99
      lgvar = ((1.d / (alog(10.0)*ntot))^2)*tvar
      sig_colm = sqrt(lgvar)
  endif else begin
      colm = ntot
      sig_colm = sqrt(tvar)
  endelse
  
  ; Free your mind
  delvarx, delv, velo, nndt

  return
end
  
     
  

