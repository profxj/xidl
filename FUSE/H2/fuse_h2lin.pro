;+ 
; NAME:
; fuse_h2lin
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

function fuse_h2lin, file

  if not keyword_set( FILE ) then file = $
    getenv('XIDL_DIR')+'/Molecules/H2/h2sort.dat'
  
;

  readcol, file, id, el, np, npp, Jp, Jpp, wave, fv, lgwf, gu, gamma, label, $
    FORMAT='i,i,i,i,i,i,d,d,f,i,d,a'

  nlin = n_elements(el)
  ;; Create structure
  tmp = {h2linstrct}
  finst = replicate(tmp, nlin)

  ;; Fill it up
  for ii=0L,nlin-1 do begin
      finst[ii].wrest = wave[ii]
      finst[ii].f = fv[ii]
      finst[ii].gamma = gamma[ii]
      finst[ii].el = el[ii]
      finst[ii].np = np[ii]
      finst[ii].npp = npp[ii]
      finst[ii].Jp = Jp[ii]
      finst[ii].Jpp = Jpp[ii]
      finst[ii].label = label[ii]
  endfor

  return, finst
end

