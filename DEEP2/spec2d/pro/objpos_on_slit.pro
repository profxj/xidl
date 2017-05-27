;+
; NAME:
;   objpos_on_slit
;
; PURPOSE:
;   find position of object(s) on a slit from design file
; 
; CALLING SEQUENCE:
;   yobj = objpos_on_slit(slitcoords, slitobjmap, slitno, nrow=nrow)
; 
; INPUTS:
;   slitcoords - slitcoords structure (bluslit) from deimos_tables
;   slitobjmap - from deimos_tables
;   slitno     - slit number (usually 0.. nslit-1)
;
; OPTIONAL INPUTS:
;   nrow       - number of rows in extracted slit
;	
; OUTPUTS:
;   yobj       - y position (in transposed extracted slit) of object(s)
;
; OPTIONAL OUPUTS:
;   nobj       - number of objects found (Can be zero)
;
; COMMENTS:
;   works fine with multiple objects, or none at all
;
; REVISION HISTORY:
;   2002-July-13  - Doug Finkbeiner, Princeton
;
;----------------------------------------------------------------------
function objpos_on_slit, slitcoords, slitobjmap, slitno, nrow=nrow, $
               nobj=nobj
  
  slitind = where(slitcoords.slitno eq slitno, ct)
  if ct eq 0 then message, 'No entry in slitcoords for this slit'
  scale_pix_per_asec = 8.52
  
  ind = where(slitobjmap.dslitid eq slitcoords[slitind].dslitid, nobj)
  if nobj eq 0 then begin 
     yobj = 0.
     return, 0
  endif 
  
  desilen = (slitobjmap[ind].topdist+slitobjmap[ind].botdist)* $
    scale_pix_per_asec

  if keyword_set(nrow) then begin 
     edgeloss = (desilen-nrow)/2
  endif else begin 
     edgeloss = 0
  endelse 
  
  yobj = slitobjmap[ind].botdist * scale_pix_per_asec - edgeloss
  
  return, yobj
end
