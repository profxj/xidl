;+ 
; NAME:
; x_aodm
;   Version 1.1
;
; PURPOSE:
;    Calculates an AODM column density given a wavelength, flux, and
;   error array AND the restwave of the transition.  The routine will
;   sum over the entire flux array so pass a sub-array in general.
;
; CALLING SEQUENCE:
;   
;   x_aodm, wave, fx, sig, rwave, colm, sig, /LOG, VELO=, FLG_SAT=
;
; INPUTS:
;   wave - Wavelength array
;   fx - Data
;   sig - Error array
;   rwave - Rest wavelength (Ang)
;
; RETURNS:
;
; OUTPUTS:
;   colm  - Column density
;   [sig]  - Error in the column density
;
; OPTIONAL KEYWORDS:
;   /LOG - Return log answers
;
; OPTIONAL OUTPUTS:
;   VELO=  -- Velocity array corresponding to the rest wavelength
;   FLG_SAT=  Keyword indicating whether the transtition was saturated
;   CST=   -- Constant to convert optical depth to column
;   NNDT=  -- Array of apparent optical depths
;
; COMMENTS:
;
; EXAMPLES:
;   x_aodm, wav[203:209], dat[203:209], sig[203:309], 1808.0126d, clm
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   19-Dec-2001 Written by JXP
;   20-Aug-2007 Modified flux floor, KLC
;-
;------------------------------------------------------------------------------

pro x_aodm, wave, fx, sig, rwave, colm, sig_colm, LOG=log, VELO=velo, $
            FLG_SAT=flg_sat, CST=cst, NNDT=nndt, CENT=cent

  if (N_params() LT 5) then begin 
    print,'Syntax - ' + $
      'x_aodm, wave, fx, sig, rwave, colm, [sig_colm], /LOG, FLG_SAT=,' + $
      ' VELO= [v1.1]'
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

  ; Set nndt
  nndt = dblarr(npix)

  ;; Bad pixels (unless no error array)
  bad = where(sig LE 0.,nbad)
  if nbad eq npix then clean = 1 ; no error array
  if nbad NE 0 and not keyword_set(clean) then begin
      sig[bad] = 0.
      fx[bad] = -1.
   endif
  
  ; Deal with saturated pixels
  sat = where(fx LE sig/5. OR fx LT 0.05, nsat, $
              complement=unsat, ncomplement=nunsat)

  flg_sat = 0
  if nsat NE 0 then begin
      ;; KLC 8/17/07
      ;; Floor flux to the corresponding limit
      if keyword_set(clean) then begin
         gd = lindgen(nsat) 
         ngd = nsat
      endif else gd = where(sig[sat] GT 0., ngd)
;      ;; JXP 3/28/06
;      ;; I know some old data is going to suffer from this..
;      ;; Need to fix those error arrays one of these days..
      if ngd NE 0 then begin
          sub = 0.05d > sig[sat[gd]]/5.d
          nndt[sat[gd]] = alog(1./sub)*cst
;          nndt[sat[gd]] = alog(5.d/sig[sat[gd]])*cst
          flg_sat = ngd
      endif 
;sig[sat[gd]] = median(sig)
  endif

  if nunsat NE 0 then nndt[unsat] = alog(1.d/fx[unsat])*cst


  ; Set delv
  cnst = x_constants()
  spl = cnst.c / 1e5  ; km/s

  if not keyword_set( VELO ) then begin
      if not keyword_set(CENT) then cent = (wave[0]+wave[npix-1])/2.d
      velo = (wave-cent)*spl/cent
  endif
  delv = velo - shift(velo,1)
  delv[0] = delv[1]

;  for i=0L,npix-2 do delv[i] = velo[i+1]-velo[i]
;  delv[npix-1] = delv[npix-2]

  ;; Sum 
  ntot = total( nndt*delv, /double )
  tvar = total( (delv*cst*sig/fx)^2, /double )
      
  ; Log as requested
  if keyword_set( LOG ) then begin
      if ntot GT 0. then begin
          colm = alog10( ntot ) 
          lgvar = ((1.d / (alog(10.0)*ntot))^2)*tvar
          sig_colm = sqrt(lgvar)
      endif else begin
          colm = -9.99
          sig_colm = sqrt(tvar)
      endelse
  endif else begin
      colm = ntot
      sig_colm = sqrt(tvar)
  endelse
  
  return
end
  
     
  

