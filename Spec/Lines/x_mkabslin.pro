;+ 
; NAME:
; x_mkabslin
;
; PURPOSE:
;    Creates database of line data
;
; CALLING SEQUENCE:
;   
;   x_mkabslin
;
; INPUTS:
;   wave       - ionic transition
;
; RETURNS:
;   f          - oscillator strength
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   nam        - Name of transition
;
; COMMENTS:
;
; EXAMPLES:
;   x_mkabslin, 1215.6701, fval, name
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   16-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_mkabslin

;
  if not keyword_set( LINDAT ) then lindat = '/u/xavier/bin/junk/lin.dat'
  if not keyword_set( ATOMDAT ) then atomdat = '/u/xavier/bin/junk/atom.dat'
  if not keyword_set( OUTFIL ) then $
    outfil = getenv('XIDL_DIR')+'/xidl/Spec/Lines/all_lin.fits'

  ;; Ion
  wv = dblarr(1000)
  fv = dblarr(1000)
  gamma = fltarr(1000)
  ion = strarr(1000)

  dumc = ' '
  dumd1 = 0.0d
  dumd2 = 0.0d

  ;; Reade atom.dat
  readcol, atomdat, atm_ion, atm_wav, atm_f, atm_gamm, FORMAT='A,D,D,F'

  ;; Loop on lin.dat
  ngd = 0L
  close, 9
  openr, 9, lindat,  ERROR = err
  if(err NE 0) then stop, 'File does not exist', err
  readf, 9, nlin
  for i=1,nlin do begin
      readf, 9, format='(f9.4,1x,a10,2x,f10.5)', dumd1, dumc, dumd2
      a = where(abs(dumd1-atm_wav) LT 0.01, na)
      case na of 
          0: begin
              print, 'x_mkabslin: No Gamma value!', dumc
          end
          1: begin
              wv[ngd] = dumd1
              fv[ngd] = dumd2
              gamma[ngd] = atm_gamm[a]
              ion[ngd] = dumc
              ngd = ngd + 1
          end
          else: begin
              print, 'x_mkabslin: More than one line!', atm_wav[a]
              stop
          end
      endcase
              
  endfor

  close, 9

  ;; Structure
  tmp = { abslinstrct }
  all_lin = replicate(tmp, ngd)
  all_lin.wrest = wv[0:ngd-1]
  all_lin.ion = ion[0:ngd-1]
  all_lin.f = fv[0:ngd-1]
  all_lin.gamma = gamma[0:ngd-1]

  
  ;; Write Data file
  mwrfits, all_lin, outfil, /create, /silent
  

  return
end
