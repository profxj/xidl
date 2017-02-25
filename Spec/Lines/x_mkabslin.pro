;+ 
; NAME:
; x_mkabslin
;
; PURPOSE:
;    Creates a FITS file from the standard line lists.  Used
;    for voigt stuff.
;
; CALLING SEQUENCE:
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
  if not keyword_set( LINDAT ) then $
    lindat = getenv('XIDL_DIR')+'/Spec/Lines/Lists/grb.lst' 
  if not keyword_set( ATOMDAT ) then $
    atomdat = getenv('XIDL_DIR')+'/Spec/Lines/Lists/atom.dat'
  if not keyword_set( OUTFIL ) then $
    outfil = getenv('XIDL_DIR')+'/Spec/Lines/all_lin.fits'
  if not keyword_set( FINEDAT ) then $
    finedat = getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'

  ;; Ion
  wv = dblarr(2000)
  fv = dblarr(2000)
  gamma = fltarr(2000)
  ion = strarr(2000)
  svA = fltarr(2000)
  svj = fltarr(2000)
  svE = fltarr(2000)

  dumc = ' '
  dumd1 = 0.0d
  dumd2 = 0.0d

  ;; Read atom.dat
  readcol, atomdat, atm_ion, atm_wav, atm_f, atm_gamm, FORMAT='A,D,D,F', $
    /silent

  ;; Read fine_strct data
  readcol, finedat, fin_Z, fin_ion, fin_j, fin_wav, fin_E, fin_A, $
           FORMAT='I,I,F,D,F,F', /silent

  ;; Loop on lin.dat
  ngd = 0L
  close, 9
  openr, 9, lindat,  ERROR = err
  if(err NE 0) then stop, 'File does not exist', err
  readf, 9, nlin
  for i=1,nlin do begin
      readf, 9, format='(f9.4,1x,a11,1x,f10.5)', dumd1, dumc, dumd2
      wv[ngd] = dumd1
      fv[ngd] = dumd2
      ion[ngd] = strtrim(dumc,2)

      ;; Gamma
      a = where(abs(dumd1-atm_wav) LT 0.01, na)
      case na of 
          0: 
          1: gamma[ngd] = atm_gamm[a]
          else: begin
              print, 'x_mkabslin: More than one line!', atm_wav[a]
              stop
          end
      endcase

      ;; Fine struct
      b = where(abs(dumd1-fin_wav) LT 0.01, nb)
      case nb of 
          0: 
          1: begin
              svA[ngd] = fin_A[b]
              svj[ngd] = fin_j[b]
              svE[ngd] = fin_E[b]
          end
          else: begin
              print, 'x_mkabslin: More than one line!', atm_wav[a]
              stop
          end
      endcase
              
      ngd = ngd + 1
  endfor

  close, 9

  ;; Structure
  tmp = { nofabslinstrct }
  all_lin = replicate(tmp, ngd)
  all_lin.wrest = wv[0:ngd-1]
  all_lin.ion = strtrim(ion[0:ngd-1],2)
  all_lin.fval = fv[0:ngd-1]
  all_lin.gamma = gamma[0:ngd-1]
  all_lin.A = svA[0:ngd-1]
  all_lin.j = svj[0:ngd-1]
  all_lin.E = svE[0:ngd-1]

  ;; Write Data file
  mwrfits, all_lin, outfil, /create, /silent
  

  return
end
