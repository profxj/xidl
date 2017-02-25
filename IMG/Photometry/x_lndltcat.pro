;+ 
; NAME:
; x_lndltcat   
;   Version 1.1
;
; PURPOSE:
;    Creates a fits file of the Landolt table.  Probably needs
;  to be run only once unless a typo is corrected.
;
; CALLING SEQUENCE:
;   x_lndltcat, fil, fitsfil
;
; INPUTS:
;   fil - Name of the landolt file
;
; RETURNS:
;
; OUTPUTS:
;   fitsfil - Name of fits file to create
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  STRCT - Structre of the landolt file
;
; COMMENTS:
;
; EXAMPLES:
;   x_lndltcat, '/u/xavier/idl/xidl/IMG/Photometry/Lists/nlandolt.dat', 
;   'Landolt.fits'
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Aug-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_lndltcat, fil, fitsfil, STRCT=strct

;
  if  N_params() LT 2  then begin 
      print, 'Syntax - ' +$
        'x_lndltcat, fil, fitsfil (v1.0)'
      return
  endif 

  close, /all

;  Optional Keywords

;  Find the Standard Stars

  readcol, FORMAT='(A,A,A,F,F,F,F,F,F,I,I,F,F,F,F,F,F)', fil, $
    nam, ra, dec, V, BV, UB, VR, RI, VI, n, m, sV, sBV, sUB, sVR, sRI, sVI

;  Create Landolt structure

  nstrs = n_elements(nam)
  tmp = { lndltstr }
  strct = replicate(tmp, nstrs)

  strct.Name = strtrim(nam,2)
  strct.RA = strtrim(RA,2)
  strct.DEC = strtrim(DEC,2)
  strct.V = V
  strct.BV = BV
  strct.UB = UB
  strct.VR = VR
  strct.RI = RI
  strct.VI = VI
  strct.nobs = n
  strct.mobs = m
  strct.sig_V = sV
  strct.sig_BV = sBV
  strct.sig_UB = sUB
  strct.sig_VR = sVR
  strct.sig_RI = sRI
  strct.sig_VI = sVI

; Output

  mwrfits, strct, fitsfil, /create
  
end
