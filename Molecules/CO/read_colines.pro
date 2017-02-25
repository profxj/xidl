;+ 
; NAME:
; read_colines
;  (V1.1)
;
; PURPOSE:
;    Sets up a line list of Molecular Hydrogen
;
; CALLING SEQUENCE:
;   
;   h2strct = fuse_h2lin([file]) 
;
; INPUTS:
;   [file] - H2 line list [default: $XIDL_DIR/FUSE/H2/h2sort.dat]
;
; RETURNS:
;   h2strct - Structure of molecular hydrogen
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   h2strct = fuse_h2lin()
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   09-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function read_colines, file

  if not keyword_set( FILE ) then file = $
    getenv('XIDL_DIR')+'/Molecules/CO/COsort.dat'
  
;
  readcol, file, id, iso, el, np, npp, Jp, Jpp, $
           wave, fv, lgwf, gu, gamma, label, $
           FORMAT='i,i,i,i,i,i,i,d,d,d,i,d,a'

  nlin = n_elements(el)
  ;; Create structure
  tmp = {colinstrct}
  finst = replicate(tmp, nlin)

  ;; Fill it up
  finst.c_isotope = iso
  finst.wrest = wave
  finst.f = 10.d^fv
  finst.gamma = gamma
  finst.el = el
  finst.np = np
  finst.npp = npp
  finst.Jp = Jp
  finst.Jpp = Jpp
  finst.label = label

  ;; Correct f-values
  zro = where(10.d^fv LE 1d-10, nz)
  if nz NE 0 then begin
      print, 'read_colines: Correcting f-values...'
      finst[zro].f = 10.^lgwf[zro] / wave[zro]
      stop
  endif

  ;; Name
  c12 = where(iso EQ 12)
  finst[c12].molecule = 'C12O16'
;  if nc13 NE 0 then finst[c13].molecule = 'C13O16'

  return, finst
end

